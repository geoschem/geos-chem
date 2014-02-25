!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_grid_mod
!
! !DESCRIPTION: Module GLOBAL\_GRID\_MOD contains variables and routines which
!  are used to specify the parameters of a GEOS-Chem global horizontal grid.
!\\  
!\\
! !INTERFACE: 
!
MODULE Global_Grid_Mod
! 
! !USES:
!
  USE CMN_SIZE_Mod             ! Size parameters
  USE CMN_GCTM_Mod             ! Physical constants
  USE Error_Mod                ! Error-handling routines

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Cleanup_Global_Grid
  PUBLIC  :: Compute_Global_Grid
  PUBLIC  :: Get_xEdge_G
  PUBLIC  :: Get_yEdge_G
  PUBLIC  :: Get_IIIPAR
  PUBLIC  :: Get_JJJPAR

!
! !REVISION HISTORY:
!  01 May 2012 - M. Payer    - Initial version, based on grid_mod.F
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC

  !=======================================================================
  ! MODULE VARIABLES:
  !
  ! IIIPAR     : Global longitude extent [# of boxes]
  ! JJJPAR     : Global latitude  extent [# of boxes]
  ! XMID_G     : Global array of grid-box lon centers   [degrees]
  ! XEDGE_G    : Global array of grid-box lon edges     [degrees]
  ! YMID_G     : Global array of grid-box lat centers   [degrees]
  ! YEDGE_G    : Global array of grid-box lat edges     [degrees]
  ! YMID_R_G   : Global array of grid-box lat centers   [radians]
  ! YEDGE_R_G  : Global array of grid-box lat edges     [radians]
  ! AREA_M2_G  : Global array of grid-box surface areas [m2]
  ! AREA_CM2_G : Global array of grid-box surface areas [cm2]
  !=======================================================================

  ! Scalars
  INTEGER              :: IIIPAR
  INTEGER              :: JJJPAR

  ! Arrays
  REAL*8,  ALLOCATABLE :: XMID_G    (:)
  REAL*8,  ALLOCATABLE :: XEDGE_G   (:)
  REAL*8,  ALLOCATABLE :: YMID_G    (:)
  REAL*8,  ALLOCATABLE :: YEDGE_G   (:)
  REAL*8,  ALLOCATABLE :: YMID_R_G  (:)
  REAL*8,  ALLOCATABLE :: YEDGE_R_G (:)
  REAL*8,  ALLOCATABLE :: AREA_M2_G (:)
  REAL*8,  ALLOCATABLE :: AREA_CM2_G(:)
  
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_global_grid
!
! !DESCRIPTION: Subroutine COMPUTE\_GLOBAL\_GRID initializes the longitude, 
!  latitude and surface area arrays. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Global_Grid( am_I_Root )
!
! !USES:
!
    USE GRID_MOD,     ONLY : GET_XOFFSET, GET_YOFFSET
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !REVISION HISTORY:
!  01 May 2012 - M. Payer    - Initial version, based on grid_mod.F
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Logical
    LOGICAL, SAVE :: FIRST = .TRUE.

    ! Scalars
    INTEGER       :: I,    J
    REAL*8        :: FMID, FEDGE

    !=================================================================
    ! COMPUTE_GRID begins here!
    !=================================================================

    ! Allocate variables on first call
    IF ( FIRST ) THEN
       CALL INIT_GLOBAL_GRID
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

#if   defined( GRID05x0666 )
    ! Overwrite YMID at poles for 0.5 x 0.666 grids
    YMID_G(1)      = -89.875d0
    YMID_G(JJJPAR) = +89.875d0
#elif defined ( GRID025x03125 )
    ! Overwrite YMID at poles for 0.25 x 0.3125 grids (bmy, 2/10/11)
    YMID_G(1)      = -89.9375d0
    YMID_G(JJJPAR) = +89.9375d0
#endif

    ! Overwrite YEDGE at poles
    YEDGE_G(1)        = -90d0
    YEDGE_G(JJJPAR+1) = +90d0

    ! Compute latitude center/edges in radians
    DO J = 1, JJJPAR
       YMID_R_G(J)  = PI_180 * YMID_G(J)
       YEDGE_R_G(J) = PI_180 * YEDGE_G(J)
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
       AREA_M2_G(J) = 2d0 * PI * Re * Re / DBLE( IIIPAR ) * &
                    ( SIN( YEDGE_R_G(J+1) ) - SIN( YEDGE_R_G(J) ) )

       ! Grid box surface areas [cm2]
       AREA_CM2_G(J) = AREA_M2_G(J) * 1d4

    ENDDO

!    ! For debug
!    IF ( am_I_Root ) THEN
!       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!       WRITE( 6, '(a)' ) 'GLOBAL GRIDS FOR NESTED SIMULATIONS'
!       WRITE( 6, '(''Global grid box longitude centers [degrees]: '')' )
!       WRITE( 6, '(8(f8.3,1x))' ) ( XMID_G(I),  I=1,IIIPAR )
!       WRITE( 6, '(a)' )
!       WRITE( 6, '(''Global grid box longitude edges [degrees]: '')' )
!       WRITE( 6, '(8(f8.3,1x))' ) ( XEDGE_G(I), I=1,IIIPAR+1 )
!       WRITE( 6, '(a)' )
!       WRITE( 6, '(''Global grid box latitude centers [degrees]: '')' )
!       WRITE( 6, '(8(f8.3,1x))' ) ( YMID_G(J),  J=1,JJJPAR )
!       WRITE( 6, '(a)' )
!       WRITE( 6, '(''Global grid box latitude edges [degrees]: '')' )
!       WRITE( 6, '(8(f8.3,1x))' ) ( YEDGE_G(J), J=1,JJJPAR+1 )
!       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
!    ENDIF

  END SUBROUTINE Compute_Global_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_iiipar
!  
! !DESCRIPTION: Function GET\_IIIPAR returns IIIPAR.
!\\
!\\
! !INTERFACE:
!
      FUNCTION Get_IIIPAR() RESULT( I )
!
! !RETURN VALUE:
!
      INTEGER              :: I  ! IIIPAR
!
! !REVISION HISTORY:
!  01 Mar 2012 - P. Kasibhatla - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      I =  IIIPAR

      END FUNCTION Get_IIIPAR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_jjjpar
!
! !DESCRIPTION: Function GET\_JJJPAR returns JJJPAR.
!\\
!\\
! !INTERFACE:
!
      FUNCTION Get_JJJPAR() RESULT( J )
!
! !RETURN VALUE:
!
      INTEGER              :: J  ! JJJPAR
!
! !REVISION HISTORY:
!  01 Mar 2012 - P. Kasibhatla - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      J =  JJJPAR

      END FUNCTION Get_JJJPAR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_xedge_g
!
! !DESCRIPTION: Function GET\_XEDGE\_G returns the longitude in degrees at the
!  western edge of a GEOS-Chem grid box
!\\
!\\
! !INTERFACE:
!
      FUNCTION Get_xEdge_G( I ) RESULT( X )
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
!  01 Mar 2012 - P. Kasibhatla - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      X = XEDGE_G( I )

      END FUNCTION Get_xEdge_G
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_yedge_g
!
! !DESCRIPTION: Function GET\_YEDGE\_G returns the latitude in degrees at the
!  southern edge of a GEOS-Chem grid box. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION Get_Yedge_G( J ) RESULT( Y )
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
!  01 Mar 2012 - P. Kasibhatla - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      Y = YEDGE_G( J )

      END FUNCTION Get_yEdge_G
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_global_grid
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_GRID initializes variables and
!  allocates module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Global_Grid
!
! !REVISION HISTORY:
!  01 May 2012 - M. Payer    - Initial version, based on grid_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    INTEGER :: AS

    !=================================================================
    ! INIT_GLOBAL_GRID begins here!
    !=================================================================

    ! Define global sizes for grid.  We need to redefine these here
    ! since for the nested grid, we set IIPAR=IIPAR and JJPAR=JJPAR
#if   defined( GRID025x03125 )
    IIIPAR = 1152
    JJJPAR = 721
#elif defined( GRID05x0666 )
    IIIPAR = 540
    JJJPAR = 361
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

  END SUBROUTINE Init_Global_Grid
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_global_grid
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_GRID deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Global_Grid
!
! !REVISION HISTORY:
!  01 May 2012 - M. Payer    - Initial version, based on grid_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( ALLOCATED( XMID_G    ) ) DEALLOCATE( XMID_G    )
    IF ( ALLOCATED( XEDGE_G   ) ) DEALLOCATE( XEDGE_G   )
    IF ( ALLOCATED( YMID_G    ) ) DEALLOCATE( YMID_G    )
    IF ( ALLOCATED( YEDGE_G   ) ) DEALLOCATE( YEDGE_G   )
    IF ( ALLOCATED( YMID_R_G  ) ) DEALLOCATE( YMID_R_G  )
    IF ( ALLOCATED( YEDGE_R_G ) ) DEALLOCATE( YEDGE_R_G )
    IF ( ALLOCATED( AREA_M2_G ) ) DEALLOCATE( AREA_M2_G )
    
  END SUBROUTINE Cleanup_Global_Grid
!EOC
END MODULE Global_Grid_Mod
