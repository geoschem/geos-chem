! $Id: comode_mod.f,v 1.1 2003/06/30 20:26:02 bmy Exp $
      MODULE COMODE_MOD
!
!******************************************************************************
!  Module COMODE_MOD contains allocatable arrays for SMVGEAR that were
!  previously contained in common blocks in header file "comode.h".
!  (bmy, 8/31/00, 2/27/03)
!  
!  In case you were wondering, "comode" stands for:
!     "COMmon blocks: Ordinary Differential Equations"
!  
!  Module Variables:
!  ============================================================================
!  (1 ) ABSHUM   : array for absolute humidity [H2O molec/cm3]
!  (2 ) AIRDENS  : array for air density [molec/cm3]
!  (3 ) CSPEC    : array of chemical species concentration [molec/cm3]
!  (4 ) ERADIUS  : array for aerosol or dust radii [cm]
!  (5 ) IXSAVE   : array of grid box longitude indices
!  (6 ) IYSAVE   : array of grid box latitude indices
!  (7 ) IZSAVE   : array of grid box altitude indices
!  (8 ) JLOP     : array of 1-D grid box indices
!  (9 ) PRESS3   : array for grid box pressure [mb]
!  (10) REMIS    : array for emissions from GEOS-CHEM [molec/cm3] 
!  (11) T3       : array for grid box temperature [K]
!  (12) TAREA    : array for surface area of aerosol or dust [cm2/cm3]
!  (13) VOLUME   : array for grid box volume [cm3]
!
!  Module Routines:
!  ============================================================================
!  (1 ) INIT_COMODE    : allocates memory for arrays
!  (2 ) CLEANUP_COMODE : deallocates memory for arrays
!
!  GEOS-CHEM modules referenced by comode_mod.f
!  ============================================================================
!  (1 ) error_mod.f    : Module containing NaN and other error check routines
!
!  NOTES:
!  (1 ) Now zero CSPEC after allocating memory (bmy, 9/8/00)
!  (2 ) Now declare more SMVGEAR arrays allocatable (bmy, 10/19/00)
!  (3 ) Updated comments (bmy, 9/4/01)
!  (4 ) Now make ERADIUS, TAREA 2-D arrays, for het chem (bmy, 11/15/01)
!  (5 ) DARSFCA is now obsolete, remove it.  Now allocate ERADIUS and
!        TAREA arrays to be of size (ITLOOP,NDUST+NAER).  (rvm, bmy, 2/27/02)
!  (5 ) Removed obsolete code from 2/02 (bmy, 4/15/02)
!  (6 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (7 ) Now references "error_mod.f" (bmy, 10/15/02)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      REAL*8,  ALLOCATABLE :: ABSHUM(:) 
      REAL*8,  ALLOCATABLE :: AIRDENS(:) 
      REAL*8,  ALLOCATABLE :: CSPEC(:,:)
      REAL*8,  ALLOCATABLE :: ERADIUS(:,:)
      INTEGER, ALLOCATABLE :: IXSAVE(:)
      INTEGER, ALLOCATABLE :: IYSAVE(:)
      INTEGER, ALLOCATABLE :: IZSAVE(:)
      INTEGER, ALLOCATABLE :: JLOP(:,:,:)
      REAL*8,  ALLOCATABLE :: PRESS3(:)      
      REAL*8,  ALLOCATABLE :: REMIS(:,:)
      REAL*8,  ALLOCATABLE :: T3(:)      
      REAL*8,  ALLOCATABLE :: TAREA(:,:)
      REAL*8,  ALLOCATABLE :: VOLUME(:)      

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS
      
!------------------------------------------------------------------------------
      
      SUBROUTINE INIT_COMODE

      !=================================================================
      ! Subroutine INIT_COMODE allocates memory for allocatable arrays
      ! that were previously contained in common blocks in "comode.h"
      !
      ! NOTES:
      ! (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
      ! (2 ) Cosmetic chagnes (bmy, 2/27/03)
      !=================================================================

      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"
#     include "comode.h" 

      ! Local variables
      INTEGER :: AS 

      ! INIT_COMODE begins here!
      WRITE( 6, 100 )
 100  FORMAT( '     - INIT_COMODE: Allocating arrays for SMVGEAR...' )

      ALLOCATE( ABSHUM( ITLOOP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ABSHUM' )
      ABSHUM = 0d0

      ALLOCATE( AIRDENS( ITLOOP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'AIRDENS' )
      AIRDENS = 0d0      

      ALLOCATE( CSPEC( ITLOOP, IGAS ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'CSPEC' )
      CSPEC = 0d0

      ALLOCATE( ERADIUS( ITLOOP, NDUST+NAER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ERADIUS' )
      ERADIUS = 0d0      

      ALLOCATE( IXSAVE( ITLOOP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IXSAVE' )
      IXSAVE = 0

      ALLOCATE( IYSAVE( ITLOOP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IYSAVE' )
      IYSAVE = 0

      ALLOCATE( IZSAVE( ITLOOP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'IZSAVE' )
      IZSAVE = 0

      ALLOCATE( JLOP( ILONG, ILAT, IPVERT ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'JLOP' )
      JLOP = 0

      ALLOCATE( PRESS3( ITLOOP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'PRESS3' )
      PRESS3 = 0d0

      ALLOCATE( REMIS( ITLOOP, MAXGL3 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'REMIS' )
      REMIS = 0d0

      ALLOCATE( T3( ITLOOP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'T3' )
      T3 = 0d0

      ALLOCATE( TAREA( ITLOOP, NDUST+NAER ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TAREA' )
      TAREA = 0d0      

      ALLOCATE( VOLUME( ITLOOP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'VOLUME' )
      VOLUME = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_COMODE

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_COMODE

      !=================================================================
      ! Subroutine INIT_COMODE deallocates memory from allocatable arrays
      ! that were previously contained in common blocks in "comode.h"
      !================================================================= 

      ! CLEANUP_COMODE begins here!
      IF ( ALLOCATED( ABSHUM  ) ) DEALLOCATE( ABSHUM  )
      IF ( ALLOCATED( AIRDENS ) ) DEALLOCATE( AIRDENS )
      IF ( ALLOCATED( CSPEC   ) ) DEALLOCATE( CSPEC   )
      IF ( ALLOCATED( ERADIUS ) ) DEALLOCATE( ERADIUS )
      IF ( ALLOCATED( IXSAVE  ) ) DEALLOCATE( IXSAVE  )
      IF ( ALLOCATED( IYSAVE  ) ) DEALLOCATE( IYSAVE  )
      IF ( ALLOCATED( IZSAVE  ) ) DEALLOCATE( IZSAVE  )
      IF ( ALLOCATED( JLOP    ) ) DEALLOCATE( JLOP    )
      IF ( ALLOCATED( PRESS3  ) ) DEALLOCATE( PRESS3  )     
      IF ( ALLOCATED( REMIS   ) ) DEALLOCATE( REMIS   )
      IF ( ALLOCATED( T3      ) ) DEALLOCATE( T3      )     
      IF ( ALLOCATED( TAREA   ) ) DEALLOCATE( TAREA   )
      IF ( ALLOCATED( VOLUME  ) ) DEALLOCATE( VOLUME  )  

      ! Return to calling program
      END SUBROUTINE CLEANUP_COMODE

!------------------------------------------------------------------------------

      END MODULE COMODE_MOD

