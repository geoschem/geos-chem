! $Id: toms_mod.f,v 1.2 2004/09/21 18:04:19 bmy Exp $
      MODULE TOMS_MOD
!
!******************************************************************************
!  Module TOMS_MOD contains variables and routines for reading the EP-TOMS
!  O3 column data from disk (for use w/ the FAST-J photolysis routines).
!  (mje, bmy, 7/14/03, 7/20/04)
!
!  Module Variables:
!  ============================================================================
!  (1 ) TOMS   (REAL*8) : Monthly mean TOMS O3 data [DU]
!  (2 ) DTOMS1 (REAL*8) : Time deriv. of TOMS O3 for 1st half of month [DU/day]
!  (3 ) DTOMS2 (REAL*8) : Time deriv. of TOMS O3 for 2nd half of month [DU/day]
!
!  Module Procedures:
!  ============================================================================
!  (1 ) READ_TOMS       : Routine to read TOMS data from disk
!  (2 ) INIT_TOMS       : Routine to allocate and zero module arrays
!  (3 ) CLEANUP_TOMS    : Routine to deallocate module arrays
!
!  GEOS-CHEM modules referenced by toms_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module containing routines for binary punch file I/O
!  (2 ) directory_mod.f : Module containing GEOS-CHEM data & met field dirs
!  (3 ) error_mod.f     : Module containing I/O error and NaN check routines
!  (4 ) time_mod.f      : Module containing routines to compute date & time
!  (5 ) transfer_mod.f  : Module containing routines to cast & resize arrays
!
!  NOTES:
!  (1 ) Now references "directory_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "toms_mod.f"
      !=================================================================
      PRIVATE :: INIT_TOMS

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      REAL*8, ALLOCATABLE :: TOMS(:,:)               
      REAL*8, ALLOCATABLE :: DTOMS1(:,:)             
      REAL*8, ALLOCATABLE :: DTOMS2(:,:)             
                                                     
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE READ_TOMS( THISMONTH, THISYEAR )
!
!******************************************************************************
!  Subroutine READ_TOMS reads in TOMS O3 column data from a binary punch
!  file for the given grid, month and year. (mje, bmy 12/10/02, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Current month (1-12)
!  (2 ) THISYEAR  (INTEGER) : Current year (YYYY format)
!
!  References:
!  ============================================================================
!
!  NOTES:
!  (1 ) Bundled into "toms_mod.f" (bmy, 7/14/03)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TIME_MOD,      ONLY : EXPAND_DATE
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"  ! Size parameters
!---------------------------------------------
!#     include "CMN_SETUP" ! DATA_DIR
!---------------------------------------------

      ! Arguments
      INTEGER, INTENT(IN) :: THISMONTH
      INTEGER, INTENT(IN) :: THISYEAR

      ! Local Variables
      LOGICAL             :: FIRST = .TRUE.
      INTEGER             :: YYYYMMDD
      REAL*4              :: ARRAY(IGLOB,JGLOB,1)
      REAL*8              :: XTAU
      CHARACTER(LEN=255)  :: FILENAME

      !=================================================================
      ! READ_TOMS begins here!
      !=================================================================

      ! Allocate arrays on the first call only
      IF ( FIRST ) THEN
         CALL INIT_TOMS
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! TOMS O3 column data currently exists between 1997 and 2002.  
      ! For other years, set data arrays to -999 (missing data flag)
      ! and then return to the calling program.
      !=================================================================
      IF ( THISYEAR < 1997 .or. THISYEAR > 2002 ) THEN

         ! Set to "missing data" value
         TOMS   = -999
         DTOMS1 = -999
         DTOMS2 = -999

         ! Echo info
         WRITE( 6, 100 ) THISYEAR
 100     FORMAT( '     - READ_TOMS: TOMS data for the year ', i4, 
     &           ' does not exist!' )

         ! Return to calling program
         RETURN
      ENDIF

      !=================================================================
      ! Read TOMS data from disk
      !=================================================================

      ! Get TAU0 value for first day of the MONTH 
      XTAU = GET_TAU0( THISMONTH, 1, THISYEAR )

      ! Define filename
      FILENAME = TRIM( DATA_DIR )                    // 
     &           'TOMS_200307/TOMS_O3col_YYYY.geos.' // GET_RES_EXT()

      ! Create YYYYMMDD value
      YYYYMMDD = ( THISYEAR * 10000 ) + ( THISMONTH * 100 ) + 01
      
      ! Replace YYYY token with current year
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

      ! Echo filename
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - READ_TOMS: Reading ', a )

      !-----------------------------
      ! TOMS O3 columns
      !-----------------------------
      
      ! Read data
      CALL READ_BPCH2( FILENAME, 'TOMS-O3',  1, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. )         

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), TOMS )

      !--------------------------------
      ! d(TOMS)/dT (1st half of month)
      !--------------------------------

       ! Read data
      CALL READ_BPCH2( FILENAME, 'TOMS-O3',  2, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. )         

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), DTOMS1 )
      
      !--------------------------------
      ! d(TOMS)/dT (2nd half of month)
      !--------------------------------

       ! Read data: 
      CALL READ_BPCH2( FILENAME, 'TOMS-O3',  3, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. )         

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), DTOMS2 )

      ! Return to calling program
      END SUBROUTINE READ_TOMS

!------------------------------------------------------------------------------

      SUBROUTINE INIT_TOMS
!
!******************************************************************************
!  Subroutine INIT_TOMS allocates and zeroes all module arrays (bmy, 7/14/03)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! INIT_TOMS begins here!
      !=================================================================

      ! Allocate TOMS
      ALLOCATE( TOMS( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TOMS' )
      TOMS = 0d0

      ! Allocate DTOMS
      ALLOCATE( DTOMS1( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DTOMS1' )
      DTOMS1 = 0d0

      ! Allocate DTOMS2
      ALLOCATE( DTOMS2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DTOMS2' )
      DTOMS2 = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_TOMS

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_TOMS
!
!******************************************************************************
!  Subroutine CLEANUP_TOMS deallocates all module arrays (bmy, 7/14/03)
!******************************************************************************
!
      IF ( ALLOCATED( TOMS   ) ) DEALLOCATE( TOMS   )
      IF ( ALLOCATED( DTOMS1 ) ) DEALLOCATE( DTOMS1 )
      IF ( ALLOCATED( DTOMS2 ) ) DEALLOCATE( DTOMS2 )

      ! Return to calling program
      END SUBROUTINE CLEANUP_TOMS

!------------------------------------------------------------------------------

      END MODULE TOMS_MOD
