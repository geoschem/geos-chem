! $Id: tropopause_mod.f,v 1.4 2006/03/24 20:22:59 bmy Exp $
      MODULE TROPOPAUSE_MOD
!
!******************************************************************************
!  Module TROPOPAUSE_MOD contains routines and variables for reading and
!  returning the value of the annual mean tropopause. (bmy, 8/15/05, 11/1/05)
! 
!  Module Variables:
!  ============================================================================
!  (1 ) LMIN   (INTEGER)     : Minimum extent of annual mean tropopause
!  (2 ) LMAX   (INTEGER)     : Maximum extent of annual mean tropopause
!  (3 ) LPAUSE (INTEGER)     : Array for annual mean tropopause
!  (4 ) IFLX   (INTEGER)     : Array for tropopause flags for ND27 (OBSOLETE)
!
!  Module Routines:
!  ============================================================================
!  (1 ) READ_TROPOPAUSE      : Reads annual mean tropopause from disk
!  (2 ) GET_MIN_TPAUSE_LEVEL : Returns min extent of ann mean tropopause
!  (3 ) GET_MAX_TPAUSE_LEVEL : Returns max extent of ann mean tropopause
!  (4 ) ITS_IN_THE_TROP      : Returns TRUE if box (I,J,L) is in troposphere
!  (5 ) ITS_IN_THE_STRAT     : Returns TRUE if box (I,J,L) is in stratosphere
!  (6 ) INIT_TROPOPAUSE      : Allocates and zeroes all module arrays
!  (7 ) CLEANUP_TROPOPAUSE   : Deallocates all module arrays
!
!  GEOS-CHEM modules referenced by tropopause_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f          : Module w/ routines for binary punch file I/O
!  (2 ) directory_mod.f      : Module w/ GEOS-CHEM met field and data dirs
!  (3 ) error_mod.f          : Module w/ NaN, other error check routines
!  (4 ) transfer_mod.f       : Module w/ routines to cast & resize arrays
!
!  NOTES:
!  (1 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (2 ) Simplify counting of tropospheric boxes (bmy, 11/1/05)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "tropopause_mod.f"
      !=================================================================

      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CLEANUP_TROPOPAUSE
      PUBLIC :: GET_MIN_TPAUSE_LEVEL
      PUBLIC :: GET_MAX_TPAUSE_LEVEL
      PUBLIC :: GET_TPAUSE_LEVEL
      PUBLIC :: ITS_IN_THE_TROP
      PUBLIC :: ITS_IN_THE_STRAT
      PUBLIC :: READ_TROPOPAUSE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      
      ! Scalars
      INTEGER              :: LMIN, LMAX

      ! Arrays
      INTEGER, ALLOCATABLE :: TROPOPAUSE(:,:)
      INTEGER, ALLOCATABLE :: IFLX(:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE READ_TROPOPAUSE
!
!******************************************************************************
!  Subroutine READ_TROPOPAUSE reads in the annual mean tropopause. 
!  (qli, bmy, 12/13/99, 11/1/05)
!
!  NOTES:
!  (1 ) Call READ_BPCH2 to read in the annual mean tropopause data
!        which is stored in binary punch file format. (bmy, 12/13/99)
!  (2 ) Now also read integer flags for ND27 diagnostic -- these determine
!        how to sum fluxes from boxes adjacent to the annual mean tropoause.
!        (qli, bmy, 1/7/00)
!  (3 ) Cosmetic changes (bmy, 3/17/00)
!  (4 ) Reference F90 module "bpch2_mod" which contains routine "read_bpch2"
!        for reading data from binary punch files (bmy, 6/28/00)
!  (5 ) Call TRANSFER_2D from "transfer_mod.f" to cast data from REAL*4 to
!        INTEGER and also to resize to (IIPAR,JJPAR).  ARRAY needs to be of 
!        size (IGLOB,JGLOB).  Also updated comments and made cosmetic changes. 
!        Removed obsolete variables.(bmy, 9/26/01)
!  (6 ) Removed obsolete code from 9/01 (bmy, 10/26/01)
!  (7 ) Now read annual mean tropopause files from the ann_mean_trop_200202/
!        subdirectory of DATA_DIR (bmy, 1/24/02)
!  (8 ) Eliminated obsolete code from 1/02 (bmy, 2/27/02)
!  (9 ) Now write file name to stdout (bmy, 4/3/02)
!  (10) Now reference GEOS_CHEM_STOP from "error_mod.f", which frees all
!        allocated memory before stopping the run. (bmy, 10/15/02)
!  (11) Now call READ_BPCH2 with QUIET=.TRUE. to suppress printing of extra
!        info to stdout.  Also updated FORMAT strings. (bmy, 3/14/03)
!  (12) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (13) Now bundled into "tropopause_mod.f' (bmy, 2/10/05)
!  (14) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (15) Simplify counting of # of tropospheric boxes (bmy, 11/1/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE ERROR_MOD,     ONLY : GEOS_CHEM_STOP
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters
#     include "CMN"           ! LPAUSE, for backwards compatibility

      ! Local Variables
      LOGICAL, SAVE          :: FIRST=.TRUE.
      INTEGER                :: I, J, COUNT
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READ_TROPOPAUSE begins here!
      !
      ! Read the annual mean tropopause from disk (binary punch file 
      ! format).  Transfer data into an array of size (IIPAR,JJPAR).
      !=================================================================
      
      ! Allocate arrays
      IF ( FIRST ) THEN
         CALL INIT_TROPOPAUSE
         FIRST = .FALSE.
      ENDIF

      ! Create filename
      FILENAME = TRIM( DATA_DIR )                      // 
     &           'ann_mean_trop_200202/ann_mean_trop.' //
     &           GET_NAME_EXT() // '.' // GET_RES_EXT()

      ! Write file name to stdout
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - READ_TROPOPAUSE: Reading ', a )

      ! Annual mean tropopause is tracer #1  
      CALL READ_BPCH2( FILENAME, 'TR-PAUSE', 1, 
     &                 0d0,       IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. )

      ! Copy from REAL*4 to INTEGER and resize to (IIPAR,JJPAR)
      CALL TRANSFER_2D( ARRAY(:,:,1), TROPOPAUSE )

      !### For backwards compatibility during transition
      LPAUSE = TROPOPAUSE

!-----------------------------------------------------------------------------
! Prior to 2/10/05:
! For now don't read in IFLUX (bmy, 1/2
!      ! Integer flags for ND27 diagnostic is tracer #4
!      CALL READ_BPCH2( FILENAME, 'TR-PAUSE', 4, 0d0, 
!     &                 IGLOB,     JGLOB,     1, ARRAY, QUIET=.TRUE. )
!
!      ! Copy from REAL*4 to INTEGER and resize to (IIPAR,JJPAR)
!      CALL TRANSFER_2D( ARRAY(:,:,1), IFLX )
!-----------------------------------------------------------------------------

      !=================================================================
      ! L <  TROPOPAUSE(I,J) are tropospheric boxes  
      ! L >= TROPOPAUSE(I,J) are stratospheric boxes
      !
      ! LMIN   = level where minimum extent of the TROPOPAUSE occurs
      ! LMAX   = level where maximum extent of the TROPOPAUSE occurs
      !
      ! LMIN-1 = level where minimum extent of the TROPOSPHERE occurs
      ! LMAX-1 = level where maximum extent of the TROPOSPHERE occurs
      !
      ! Write LMAX-1 and LMIN-1 to the standard output.
      !
      ! Also make sure that LMAX-1 does not exceed LLTROP, since LLTROP 
      ! is used to dimension the chemistry arrays in "comode.h". 
      !=================================================================
      LMIN = MINVAL( TROPOPAUSE )
      LMAX = MAXVAL( TROPOPAUSE )
      
      WRITE( 6, 120 ) LMIN-1
 120  FORMAT( '     - READ_TROPOPAUSE: Minimum tropospheric extent,',
     &        ' L=1 to L=', i3 )

      WRITE( 6, 130 ) LMAX-1
 130  FORMAT( '     - READ_TROPOPAUSE: Maximum tropospheric extent,',
     &        ' L=1 to L=', i3 )
    
      IF ( LMAX-1 > LLTROP ) THEN
         WRITE( 6, '(a)' ) 'READ_TROPOPAUSE: LLTROP is set too low!' 
         WRITE( 6, 131   ) LMAX-1, LLTROP
 131     FORMAT( 'LMAX = ', i3, '  LLTROP = ', i3 )
         WRITE( 6, '(a)' ) 'STOP in READ_TROPOPAUSE.F!!!'
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! Write the number of tropopsheric and stratospheric boxes.
      ! Recall that tropospheric boxes extend up to TROPOPAUSE - 1.
      !=================================================================
      COUNT = SUM( TROPOPAUSE - 1 )

      WRITE( 6, 140 ) COUNT
 140  FORMAT( '     - READ_TROPOPAUSE: # of tropopsheric boxes:  ', i8 )
      
      WRITE( 6, 150 ) ( IIPAR * JJPAR * LLPAR ) - COUNT
 150  FORMAT( '     - READ_TROPOPAUSE: # of stratospheric boxes: ', i8 )

      ! Return to calling program
      END SUBROUTINE READ_TROPOPAUSE

!------------------------------------------------------------------------------
      
      FUNCTION GET_MAX_TPAUSE_LEVEL() RESULT( L_MAX )
!
!******************************************************************************
!  Function GET_MAX_TPAUSE_LEVEL returns GEOS-CHEM level at the highest extent
!  of the annual mean tropopause. (bmy, 2/10/05)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: L_MAX

      !=================================================================
      ! GET_MAX_TPAUSE_LEVEL begins here!
      !=================================================================
      L_MAX = LMAX

      ! Return to calling program
      END FUNCTION GET_MAX_TPAUSE_LEVEL

!------------------------------------------------------------------------------
      
      FUNCTION GET_MIN_TPAUSE_LEVEL() RESULT( L_MIN )
!
!******************************************************************************
!  Function GET_MIN_TPAUSE_LEVEL returns GEOS-CHEM level at the lowest extent
!  of the annual mean tropopause. (bmy, 2/10/05)
!
!  NOTES:
!******************************************************************************
!
      ! Function value
      INTEGER :: L_MIN

      !=================================================================
      ! GET_MIN_TPAUSE_LEVEL begins here!
      !=================================================================
      L_MIN = LMIN

      ! Return to calling program
      END FUNCTION GET_MIN_TPAUSE_LEVEL

!------------------------------------------------------------------------------

      FUNCTION GET_TPAUSE_LEVEL( I, J ) RESULT( L_TP )
!
!******************************************************************************
!  Function GET_TPAUSE_LEVEL returns the model level L_TP which contains the 
!  GEOS_CHEM annual mean tropopause at grid box (I,J).  Note that L_TP is
!  considered to be in the stratosphere.  Levels L_TP-1 and below are 
!  considered to be purely tropospheric levels. (bmy, 8/22/05)
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I, J

      ! Local variables
      INTEGER             :: L_TP
      
      !=================================================================
      ! GET_TPAUSE_LEVEL begins here!
      !=================================================================
      L_TP = TROPOPAUSE(I,J)
      
      ! Return to calling program
      END FUNCTION GET_TPAUSE_LEVEL

!------------------------------------------------------------------------------

      FUNCTION ITS_IN_THE_TROP( I, J, L ) RESULT ( IS_TROP )
!
!******************************************************************************
!  Function ITS_IN_THE_TROP returns TRUE if grid box (I,J,L) lies within
!  the troposphere, or FALSE otherwise. (bmy, 2/10/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) : GEOS-CHEM longitude index
!  (2 ) J (INTEGER) : GEOS-CHEM latitude  index
!  (3 ) L (INTEGER) : GEOS-CHEM level     index
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Return value
      LOGICAL             :: IS_TROP

      !=================================================================
      ! ITS_IN_THE_TROP begins here
      !=================================================================
      IS_TROP = ( L < TROPOPAUSE(I,J) ) 

      ! Return to calling program
      END FUNCTION ITS_IN_THE_TROP

!------------------------------------------------------------------------------

      FUNCTION ITS_IN_THE_STRAT( I, J, L ) RESULT( IS_STRAT )
!
!******************************************************************************
!  Function ITS_IN_THE_STRAT returns TRUE if grid box (I,J,L) lies within
!  the stratosphere, or FALSE otherwise. (bmy, 2/10/05)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I (INTEGER) : GEOS-CHEM longitude index
!  (2 ) J (INTEGER) : GEOS-CHEM latitude  index
!  (3 ) L (INTEGER) : GEOS-CHEM level     index
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L

      ! Return value
      LOGICAL             :: IS_STRAT

      !=================================================================
      ! ITS_IN_THE_STRAT begins here
      !=================================================================
      IS_STRAT = ( L >= TROPOPAUSE(I,J) ) 

      ! Return to calling program
      END FUNCTION ITS_IN_THE_STRAT

!------------------------------------------------------------------------------

      SUBROUTINE INIT_TROPOPAUSE
!
!******************************************************************************
!  Subroutine INIT_TROPOPAUSE allocates & zeroes module arrays. (bmy, 2/10/05)
!  
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"

      INTEGER :: AS 

      !=================================================================
      ! INIT_TROPOPAUSE
      !=================================================================
      ALLOCATE( TROPOPAUSE( IIPAR, JJPAR ), STAT=AS ) 
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TROPOPAUSE' )
      TROPOPAUSE = 0

      ! For now don't allocate IFLX
      !ALLOCATE( IFLX( IIPAR, JJPAR ), STAT=AS ) 
      !IF ( AS /= 0 ) CALL ALLOC_ERR( 'IFLX' )
      !IFLX = 0

      ! Return to calling program
      END SUBROUTINE INIT_TROPOPAUSE

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_TROPOPAUSE
!
!******************************************************************************
!  Subroutine CLEANUP_TROPOPAUSE deallocates module arrays (bmy, 2/10/05)
! 
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_TROPOPAUSE begins here!
      !=================================================================
      IF ( ALLOCATED( TROPOPAUSE ) ) DEALLOCATE( TROPOPAUSE )
      IF ( ALLOCATED( IFLX       ) ) DEALLOCATE( IFLX       ) 

      ! Return to calling program
      END SUBROUTINE CLEANUP_TROPOPAUSE 

!------------------------------------------------------------------------------

      ! End of module
      END MODULE TROPOPAUSE_MOD
