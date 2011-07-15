! $Id: grid_mod.f,v 1.1 2009/11/20 21:43:05 bmy Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: grid_mod.f
!
! !DESCRIPTION: Module GRID\_MOD contains variables and routines which are 
!  used to specify the parameters of a GEOS-CHEM horizontal grid.
!\\
!\\
! !INTERFACE: 
!
      MODULE GRID_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_GRID
      PUBLIC  :: COMPUTE_GRID
      PUBLIC  :: GET_AREA_M2
      PUBLIC  :: GET_AREA_CM2
      PUBLIC  :: GET_BOUNDING_BOX
      PUBLIC  :: GET_XEDGE
      PUBLIC  :: GET_XMID
      PUBLIC  :: GET_XOFFSET
      PUBLIC  :: GET_YOFFSET
      PUBLIC  :: GET_YEDGE
      PUBLIC  :: GET_YEDGE_R
      PUBLIC  :: GET_YMID
      PUBLIC  :: GET_YMID_R
      PUBLIC  :: GET_YMID_R_W
      PUBLIC  :: SET_XOFFSET
      PUBLIC  :: SET_YOFFSET
      PUBLIC  :: ITS_A_NESTED_GRID
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: INIT_GRID
!
! !REVISION HISTORY:
!  11 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Fixed typos in "grid_mod.f" (bmy, 4/28/03)
!  (2 ) Added routine GET_BOUNDING_BOX.  Now define 1x125 grid. (bmy, 12/1/04)
!  (3 ) Modified for GCAP 4x5 horizontal grid (swu, bmy, 5/24/05)
!  (4 ) Added comments re: surface area derivation (bmy, 4/20/06)
!  (5 ) Modifications for GEOS-5 nested grids (yxw, dan, bmy, 11/6/08)
!  20 Nov 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC

      !=======================================================================
      ! MODULE VARIABLES:
      !
      ! (1 ) IS_NESTED  : =T if we are using a nested-grid 
      ! (2 ) IIIPAR     : Global longitude extent [# of boxes]
      ! (3 ) JJJPAR     : Global latitude  extent [# of boxes]
      ! (4 ) I0         : Nested-grid offset in longitude (X) dimension
      ! (5 ) J0         : Nested-grid offset in latitude  (Y) dimension
      ! (6 ) XMID_G     : GLOBAL array of grid-box lon centers [degrees]
      ! (7 ) XEDGE_G    : GLOBAL array of grid-box lon edges   [degrees]
      ! (8 ) YMID_G     : GLOBAL array of grid-box lat centers [degrees]
      ! (9 ) YEDGE_G    : GLOBAL array of grid-box lat edges   [degrees]
      ! (10) YMID_R_G   : GLOBAL array of grid-box lat centers [radians]
      ! (11) YEDGE_R_G  : GLOBAL array of grid-box lat edges   [radians]
      ! (12) AREA_M2_G  : GLOBAL array of grid-box surface areas [m2]
      ! (13) AREA_CM2_G : GLOBAL array of grid-box surface areas [cm2]
      ! (14) XMID       : WINDOW array of grid-box lon centers [degrees]
      ! (15) XEDGE      : WINDOW array of grid-box lon edges   [degrees]
      ! (16) YMID       : WINDOW array of grid-box lat centers [degrees]
      ! (17) YEDGE      : WINDOW array of grid-box lat edges   [degrees]
      ! (18) YMID_R     : WINDOW array of grid-box lat centers [radians]
      ! (19) YEDGE_R    : WINDOW array of grid-box lat edges   [radians]
      ! (20) AREA_M2    : WINDOW array of grid-box surface areas [m2]
      ! (21) AREA_CM2   : WINDOW array of grid-box surface areas [cm2]
      !=======================================================================
      LOGICAL              :: IS_NESTED
      INTEGER              :: I0,           J0
      INTEGER              :: IIIPAR,       JJJPAR
      REAL*8,  ALLOCATABLE :: XMID_G(:),    XEDGE_G(:)
      REAL*8,  ALLOCATABLE :: YMID_G(:),    YEDGE_G(:)
      REAL*8,  ALLOCATABLE :: YMID_R_G(:),  YEDGE_R_G(:)
      REAL*8,  ALLOCATABLE :: AREA_M2_G(:), AREA_CM2_G(:)
      REAL*8,  ALLOCATABLE :: XMID(:),      XEDGE(:)
      REAL*8,  ALLOCATABLE :: YMID(:),      YEDGE(:)
      REAL*8,  ALLOCATABLE :: YMID_R(:),    YEDGE_R(:)
      REAL*8,  ALLOCATABLE :: YMID_R_W(:),  YEDGE_R_W(:)
      REAL*8,  ALLOCATABLE :: AREA_M2(:),   AREA_CM2(:)

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
      SUBROUTINE COMPUTE_GRID
!
! !USES:
!
#     include "F77_CMN_SIZE"  ! Size parameters
#     include "F77_CMN_GCTM"  ! Physical constants
!
! !REVISION HISTORY:
!  11 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Added fancy output (bmy, 4/26/04)
!  (2 ) Suppress some output lines (bmy, 7/20/04)
!  (3 ) Now also support 1 x 1.25 grid (bmy, 12/1/04)
!  (4 ) Now modified for GCAP 4x5 horizontal grid (swu, bmy, 5/24/05)
!  (5 ) Added comments re: surface area derivation (bmy, 4/20/06)
!  (6 ) Compute YMID, YEDGE for 0.5x0.666 nested grids (yxw, dan, bmy, 11/6/08)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: I,    J
      REAL*8        :: FMID, FEDGE

      !=================================================================
      ! COMPUTE_GRID begins here!
      !=================================================================

      ! Allocate variables on first call
      IF ( FIRST ) THEN
         CALL INIT_GRID
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Compute latitude centers & edges (algorithm from old "input.f")
      !=================================================================
      FMID  = 0.5d0 * DBLE( JJJPAR + 1 )
      FEDGE = 0.5d0 * DBLE( JJJPAR + 2 )

      DO J = 1, JJJPAR
         YMID_G(J)  = DJSIZE * ( DBLE(J) - FMID  )
         YEDGE_G(J) = DJSIZE * ( DBLE(J) - FEDGE )
      ENDDO

#if   defined( GRID4x5 ) && defined( GCAP )
      ! Overwrite YMID at poles for GCAP 4 x 5 grid (swu, bmy, 5/24/05)
      YMID_G(1)      = -88.d0
      YMID_G(JJJPAR) = +88.d0      

#elif defined( GRID4x5 )
      ! Overwrite YMID at poles for 4 x 5 grid
      YMID_G(1)      = -89.d0
      YMID_G(JJJPAR) = +89.d0

#elif defined( GRID2x25 )
      ! Overwrite YMID at poles for 2 x 2.5 grid
      YMID_G(1)      = -89.5d0
      YMID_G(JJJPAR) = +89.5d0

#elif defined ( GRID1x1 ) || defined( GRID1x125 )
      ! Overwrite YMID at poles for 1 x 1 and 1 x 1.25 grids
      YMID_G(1)      = -89.75d0
      YMID_G(JJJPAR) = +89.75d0

#elif defined ( GRID05x0666 )
      ! Overwrite YMID at poles for 0.5 x 0.666 grids
      YMID_G(1)      = -89.875d0
      YMID_G(JJJPAR) = +89.875d0

#endif

      ! Overwrite YEDGE at poles
      YEDGE_G(1)        = -90d0
      YEDGE_G(JJJPAR+1) = +90d0

      ! Compute latitude center/edges in radians
      DO J = 1, JJJPAR
         YMID_R_G(J)  = ( PI / 180d0 ) * YMID_G(J)
         YEDGE_R_G(J) = ( PI / 180d0 ) * YEDGE_G(J)
      ENDDO
         
      ! Overwrite RLATV at N. pole
      YEDGE_R_G(JJJPAR+1) = PI / 2d0

      !=================================================================
      ! Compute longitude centers & edges (algorithm from old "input.f")
      !=================================================================      
      XMID_G(1)  = -180d0
      XEDGE_G(1) = XMID_G(1) - ( DISIZE / 2d0 )

      DO I = 1, IIIPAR-1
         XMID_G(I+1)  = XMID_G(I)  + DISIZE
      ENDDO

      DO I = 1, IIIPAR
         XEDGE_G(I+1) = XEDGE_G(I) + DISIZE 
      ENDDO
      
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
      !
      ! (bmy, 4/20/06)
      !=================================================================  
      DO J = 1, JJJPAR

         ! Grid box surface areas [m2]
         AREA_M2_G(J) = 2d0 * PI * Re * Re / DBLE( IIIPAR ) *
     &                  ( SIN( YEDGE_R_G(J+1) ) - SIN( YEDGE_R_G(J) ) ) 

         ! Grid box surface areas [cm2]
         AREA_CM2_G(J) = AREA_M2_G(J) * 1d4

      ENDDO

      !=================================================================
      ! Save to local size arrays so that we can index for all grids
      !=================================================================
      
      ! XMID
      DO I = 1, IIPAR
         XMID(I) = XMID_G(I+I0)
      ENDDO

      ! XEDGE
      DO I = 1, IIPAR+1
         XEDGE(I) = XEDGE_G(I+I0)
      ENDDO

      ! YMID, YMID_R, AREA_M2, AREA_CM2
      DO J = 1, JJPAR
         YMID(J)     = YMID_G(J+J0)
         YMID_R(J)   = YMID_R_G(J+J0)
         AREA_M2(J)  = AREA_M2_G(J+J0)
         AREA_CM2(J) = AREA_CM2_G(J+J0)
      ENDDO

#if   defined( GRID05x0666 )
      ! This only needs to be done for GEOS-5 nested grid (dan, bmy, 11/18/08)
      DO J = 0, JJPAR+1
         YMID_R_W(J) = YMID_R_G(J+J0)
      ENDDO
#endif

      ! YEDGE, YEDGE_R
      DO J = 1, JJPAR+1
         YEDGE(J)   = YEDGE_G(J+J0)
         YEDGE_R(J) = YEDGE_R_G(J+J0)
      ENDDO

      !=================================================================
      ! Echo info to stdout
      !=================================================================
      WRITE( 6, '(''Nested-Grid X-offset [boxes]:'', i4 )' ) I0
      WRITE( 6, '(''Nested-Grid Y-offset [boxes]:'', i4 )' ) J0
      WRITE( 6, '(a)' )
      WRITE( 6, '(''Grid box longitude centers [degrees]: '')' )
      WRITE( 6, '(8(f8.3,1x))' ) ( XMID(I),  I=1,IIPAR )
      WRITE( 6, '(a)' )
      WRITE( 6, '(''Grid box longitude edges [degrees]: '')' )
      WRITE( 6, '(8(f8.3,1x))' ) ( XEDGE(I), I=1,IIPAR+1 )
      WRITE( 6, '(a)' )
      WRITE( 6, '(''Grid box latitude centers [degrees]: '')' )
      WRITE( 6, '(8(f8.3,1x))' ) ( YMID(J),  J=1,JJPAR )
      WRITE( 6, '(a)' )
      WRITE( 6, '(''Grid box latitude edges [degrees]: '')' )
      WRITE( 6, '(8(f8.3,1x))' ) ( YEDGE(J), J=1,JJPAR+1 )

      !=================================================================
      ! Deallocate global arrays -- we don't need these anymore
      !=================================================================
      IF ( ALLOCATED( XMID_G     ) ) DEALLOCATE( XMID_G     )
      IF ( ALLOCATED( XEDGE_G    ) ) DEALLOCATE( XEDGE_G    )
      IF ( ALLOCATED( YMID_G     ) ) DEALLOCATE( YMID_G     )
      IF ( ALLOCATED( YEDGE_G    ) ) DEALLOCATE( YEDGE_G    )
      IF ( ALLOCATED( YMID_R_G   ) ) DEALLOCATE( YMID_R_G   )
      IF ( ALLOCATED( YEDGE_R_G  ) ) DEALLOCATE( YEDGE_R_G  )
      IF ( ALLOCATED( AREA_M2_G  ) ) DEALLOCATE( AREA_M2_G  )
      IF ( ALLOCATED( AREA_CM2_G ) ) DEALLOCATE( AREA_CM2_G )

      END SUBROUTINE COMPUTE_GRID
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      SUBROUTINE SET_XOFFSET( X_OFFSET )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: X_OFFSET  ! Value to assign to I0
!
! !REVISION HISTORY:
!  11 Mar 2003 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      I0 = X_OFFSET

      END SUBROUTINE SET_XOFFSET
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
      SUBROUTINE SET_YOFFSET( Y_OFFSET )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: Y_OFFSET  ! Value to assign to J0
!
! !REVISION HISTORY:
!  11 Mar 2003 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC

      J0 = Y_OFFSET

      END SUBROUTINE SET_YOFFSET
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_xoffset
!
! !DESCRIPTION: Function GET\_XOFFSET returns the nested-grid longitude offset
!  to the calling program. (bmy, 3/11/03)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_XOFFSET( GLOBAL ) RESULT( X_OFFSET )
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
!  11 Mar 2003 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
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

      END FUNCTION GET_XOFFSET
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_xoffset
!
! !DESCRIPTION: Function GET\_XOFFSET returns the nested-grid longitude offset
!  to the calling program. (bmy, 3/11/03)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_YOFFSET( GLOBAL ) RESULT( Y_OFFSET )
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
!  11 Mar 2003 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
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

      END FUNCTION GET_YOFFSET
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_xmid
!
! !DESCRIPTION: Function GET\_XMID returns the longitude in degrees at the 
!  center of a GEOS-Chem grid box.  Works for nested-grids too.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_XMID( I ) RESULT( X )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I  ! Longitude index
!
! !RETURN VALUE:
!
      REAL*8              :: X  ! Corresponding lon value @ grid box ctr
!
! !REVISION HISTORY:
!  11 Mar 2003 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      X = XMID( I )

      END FUNCTION GET_XMID
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_xedge
!
! !DESCRIPTION: Function GET\_XEDGE returns the longitude in degrees at the 
!  western edge of a GEOS-Chem grid box.  Works for nested-grids too.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_XEDGE( I ) RESULT( X )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: I  ! Longitude index
!
! !RETURN VALUE:
!
      REAL*8              :: X  ! Corresponding lon value @ W edge of grid box
!
! !REVISION HISTORY:
!  11 Mar 2003 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      X = XEDGE( I )

      END FUNCTION GET_XEDGE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ymid
!
! !DESCRIPTION: Function GET\_YMID returns the latitude in degrees at the 
!  center of a GEOS-Chem grid box.  Works for nested-grids too.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_YMID( J ) RESULT( Y )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: J  ! Latitude index
!
! !RETURN VALUE:
!
      REAL*8              :: Y  ! Latitude value at @ grid box ctr [degrees]
!
! !REVISION HISTORY:
!  11 Mar 2003 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      Y = YMID( J )

      END FUNCTION GET_YMID
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_yedge
!
! !DESCRIPTION: Function GET\_YEDGE returns the latitude in degrees at the 
!  southern edge of a GEOS-Chem grid box.  Works for nested-grids too.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_YEDGE( J ) RESULT( Y )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: J  ! Latitude index
!
! !RETURN VALUE:
!
      REAL*8              :: Y  ! Latitude value @ S edge of grid box [degrees]
!
! !REVISION HISTORY:
!  11 Mar 2003 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      Y = YEDGE( J )

      END FUNCTION GET_YEDGE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ymid
!
! !DESCRIPTION: Function GET\_YMID\_R returns the latitude in radians at the 
!  center of a GEOS-Chem grid box.  Works for nested-grids too.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_YMID_R( J ) RESULT( Y )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: J  ! Latitude index
!
! !RETURN VALUE:
!
      REAL*8              :: Y  ! Latitude value at @ grid box ctr [radians]
!
! !REVISION HISTORY:
!  11 Mar 2003 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      Y = YMID_R( J )

      END FUNCTION GET_YMID_R
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ymid
!
! !DESCRIPTION: Function GET\_YMID\_R\_W returns the latitude in radians at 
!  the center of a GEOS-Chem grid box for the GEOS-5 nested grid.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_YMID_R_W( J ) RESULT( Y )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: J  ! Latitude index
!
! !RETURN VALUE:
!
      REAL*8              :: Y  ! Latitude value at @ grid box ctr [radians]
!
! !REVISION HISTORY:
!  06 Nov 2008 - D. Chen & R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      Y = YMID_R_W( J )

      END FUNCTION GET_YMID_R_W
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_yedge_r
!
! !DESCRIPTION: Function GET\_YEDGE\_R returns the latitude in radians at the 
!  southern edge of a GEOS-Chem grid box.  Works for nested-grids too.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_YEDGE_R( J ) RESULT( Y )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: J  ! Latitude index
!
! !RETURN VALUE:
!
      REAL*8              :: Y  ! Latitude value @ S edge of grid box [radians]
!
! !REVISION HISTORY:
!  11 Mar 2003 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC

      Y = YEDGE_R( J )

      END FUNCTION GET_YEDGE_R
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_area_m2
!
! !DESCRIPTION: Function GET\_AREA\_M2 returns the surface area [m2] of a 
!  GEOS-Chem grid box.  Works for nested grids too.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_AREA_M2( J ) RESULT( A )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: J  ! Latitude index
!
! !RETURN VALUE:
!
      REAL*8              :: A  ! Grid box surface area [m2]
!
! !REVISION HISTORY:
!  11 Mar 2003 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      A = AREA_M2( J )

      END FUNCTION GET_AREA_M2
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
      FUNCTION GET_AREA_CM2( J ) RESULT( A )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: J  ! Latitude index
!
! !RETURN VALUE:
!
      REAL*8              :: A  ! Grid box surface area [cm2]
!
! !REVISION HISTORY:
!  11 Mar 2003 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      A = AREA_CM2( J )

      END FUNCTION GET_AREA_CM2
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
      SUBROUTINE GET_BOUNDING_BOX( COORDS, INDICES )
!
! !USES:
!
      USE ERROR_MOD, ONLY : ERROR_STOP

#     include "F77_CMN_SIZE"    ! Size parameters
!
! !INPUT PARAMETERS: 
!
      REAL*8,  INTENT(IN)  :: COORDS(4)   ! (/LON_LL, LAT_LL, LON_UR, LAT_UR/)
!
! !INPUT/OUTPUT PARAMETERS: 
!
      INTEGER, INTENT(OUT) :: INDICES(4)  ! (/I_LL, J_LL, I_UR, J_UR/)
!
! !REVISION HISTORY:
!  01 Dec 2004 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER              :: I, J
      CHARACTER(LEN=255)   :: LOCATION

      !=================================================================
      ! GET_BOUNDING_BOX begins here!
      !=================================================================
      
      ! Location
      LOCATION = 'GET_BOUNDING_BOX (grid_mod.f)'

      ! Initialize
      INDICES(:) = 0

      !=================================================================
      ! Longitude search
      !=================================================================
      DO I = 1, IIPAR

         ! Locate index corresponding to the lower-left longitude
         IF ( COORDS(1) >  XEDGE(I  )   .and. 
     &        COORDS(1) <= XEDGE(I+1) ) INDICES(1) = I

         ! Locate index corresponding to upper-right longitude
         IF ( COORDS(3) >  XEDGE(I  )   .and. 
     &        COORDS(3) <= XEDGE(I+1) ) INDICES(3) = I

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
      DO J = 1, JJPAR

         ! Locate index corresponding to the lower-left latitude
         IF ( COORDS(2) >  YEDGE(J  )   .and. 
     &        COORDS(2) <= YEDGE(J+1) ) INDICES(2) = J

         ! Locate index corresponding to the upper-right latitude
         IF ( COORDS(4) >  YEDGE(J  )   .and. 
     &        COORDS(4) <= YEDGE(J+1) ) INDICES(4) = J

      ENDDO

      ! Error check lower-left longitude
      IF ( INDICES(2) == 0 ) THEN
         CALL ERROR_STOP( 'Invalid lower-left lat index!',  LOCATION )
      ENDIF

      ! Error check upper-right longitude
      IF ( INDICES(4) == 0 ) THEN
         CALL ERROR_STOP( 'Invalid upper-right lat index!', LOCATION )
      ENDIF

      END SUBROUTINE GET_BOUNDING_BOX
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
!  11 Mar 2003 - R. Yantosca - Initial version
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
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
      SUBROUTINE INIT_GRID
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "F77_CMN_SIZE"
!
! !REVISION HISTORY:
!  11 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Fixed typos that caused AREA_CM2_G and AREA_CM2 to be initialized 
!        before they were allocated. (bmy, 4/28/03)
!  (2 ) Now define IIIPAR & JJJPAR for 1 x 1.25 grid (bmy, 12/1/04)
!  (3 ) Modified for GCAP 4x5 horizontal grid (swu, bmy, 5/24/05)
!  (4 ) Modifications for 0.5 x 0.666 nested grids (dan, bmy, 11/6/08)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
      INTEGER :: AS

      !=================================================================
      ! INIT_GRID begins here!
      !=================================================================

      ! Define global sizes for grid.  We need to redefine these here
      ! since for the nested grid, we set IIPAR=IIPAR and JJPAR=JJPAR
#if   defined( GRID1x1 )
      IIIPAR = 360
      JJJPAR = 181
#elif defined( GRID05x0666 )
      IIIPAR = 540
      JJJPAR = 361
#elif defined( GRID1x125 )
      IIIPAR = 288
      JJJPAR = 181
#elif defined( GRID2x25 )
      IIIPAR = 144
      JJJPAR = 91
#elif defined( GRID4x5 ) && defined( GCAP )
      IIIPAR = 72
      JJJPAR = 45
#elif defined( GRID4x5 )
      IIIPAR = 72
      JJJPAR = 46
#endif

      !=================================================================
      ! Allocate global-sized arrays (e.g. use IIIPAR, JJJPAR)
      !=================================================================
      ALLOCATE( XMID_G( IIIPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'XMID_G' )
      XMID_G = 0

      ALLOCATE( XEDGE_G( IIIPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'XEDGE_G' )
      XEDGE_G = 0d0

      ALLOCATE( YMID_G( JJJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'YMID_G' )
      YMID_G = 0d0

      ALLOCATE( YEDGE_G( JJJPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'YEDGE_G' )
      YEDGE_G = 0d0

      ALLOCATE( YMID_R_G( JJJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'YMID_R_G' )
      YMID_R_G = 0d0

      ALLOCATE( YEDGE_R_G( JJJPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'YEDGE_R_G' )
      YEDGE_R_G = 0d0

      ALLOCATE( AREA_M2_G( JJJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AREA_M2_G' )
      AREA_M2_G = 0d0

      ALLOCATE( AREA_CM2_G( JJJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AREA_CM2_G' )
      AREA_CM2_G = 0d0      

      !=================================================================
      ! Allocate window-sized arrays (e.g. use IIPAR, JJPAR)
      !=================================================================
      ALLOCATE( XMID( IIPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'XMID' )
      XMID = 0

      ALLOCATE( XEDGE( IIPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'XEDGE' )
      XEDGE = 0d0

      ALLOCATE( YMID( JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'YMID' )
      YMID = 0d0

      ALLOCATE( YEDGE( JJPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'YEDGE' )
      YEDGE = 0d0

      ALLOCATE( YMID_R( 1:JJPAR ), STAT=AS )               
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'YMID_R' )
      YMID_R = 0d0

      ALLOCATE( YMID_R_W( 0:JJPAR+1 ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'YMID_R_W' )
      YMID_R_W = 0d0

      ALLOCATE( YEDGE_R( JJPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'YEDGE_R' )
      YEDGE_R = 0d0

      ALLOCATE( AREA_M2( JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AREA_M2' )
      AREA_M2 = 0d0

      ALLOCATE( AREA_CM2( JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AREA_CM2' )
      AREA_CM2 = 0d0

      !=================================================================
      ! Also test for 1x1 nested-grid (smaller than global size)
      !=================================================================
      IF ( IIPAR == IIIPAR .and. JJPAR == JJJPAR ) THEN
         IS_NESTED = .FALSE.
      ELSE
         IS_NESTED = .TRUE.
      ENDIF

      ! Return to calling program
      END SUBROUTINE INIT_GRID
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
      SUBROUTINE CLEANUP_GRID
!
! !REVISION HISTORY:
!  11 Mar 2003 - R. Yantosca - Initial version
!  20 Nov 2009 - D. Chen     - Now also deallocate YMID_R_W
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!EOP
!------------------------------------------------------------------------------
!BOC
      IF ( ALLOCATED( XMID       ) ) DEALLOCATE( XMID       )
      IF ( ALLOCATED( XEDGE      ) ) DEALLOCATE( XEDGE      )
      IF ( ALLOCATED( YMID       ) ) DEALLOCATE( YMID       )
      IF ( ALLOCATED( YEDGE      ) ) DEALLOCATE( YEDGE      )
      IF ( ALLOCATED( YMID_R     ) ) DEALLOCATE( YMID_R     )
      IF ( ALLOCATED( YMID_R_W   ) ) DEALLOCATE( YMID_R_W   )  
      IF ( ALLOCATED( YEDGE_R    ) ) DEALLOCATE( YEDGE_R    )
      IF ( ALLOCATED( AREA_M2    ) ) DEALLOCATE( AREA_M2    )
      IF ( ALLOCATED( AREA_CM2   ) ) DEALLOCATE( AREA_CM2   )

      END SUBROUTINE CLEANUP_GRID
!EOC
      END MODULE GRID_MOD
