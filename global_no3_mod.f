! $Id: global_no3_mod.f,v 1.4 2005/02/10 19:53:26 bmy Exp $
      MODULE GLOBAL_NO3_MOD
!
!******************************************************************************
!  Module GLOBAL_NO3_MOD contains variables and routines for reading the
!  global monthly mean NO3 concentration from disk.  These are needed for the 
!  offline sulfate/aerosol simulation. (bmy, 10/15/02, 1/15/05)
!
!  Module Variables:
!  ===========================================================================
!  (1 ) NO3 (REAL*8)       : Stores global monthly mean NO3 field
!  
!  Module Routines:
!  =========================================================================== 
!  (1 ) GET_GLOBAL_NO3     : Reads global monthly mean HNO3 from disk
!  (2 ) INIT_GLOBAL_NO3    : Allocates & initializes the HNO3 array
!  (3 ) CLEANUP_GLOBAL_NO3 : Deallocates the HNO3 array
!
!  GEOS-CHEM modules referenced by global_no3_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module containing routines for binary punch file I/O
!  (2 ) directory_mod.f : Module containing GEOS-CHEM data and met field dirs
!  (3 ) error_mod.f     : Module containing NaN and other error check routines
!  (4 ) transfer_mod.f  : Module containing routines to cast & resize arrays
!
!  NOTES:
!  (1 ) Adapted from "global_oh_mod.f" (bmy, 10/3/02)
!  (2 ) Minor bug fix in FORMAT statements (bmy, 3/23/03)
!  (3 ) Cosmetic changes (bmy, 3/27/03)
!  (4 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (5 ) Now suppress output from READ_BPCH2 with QUIET=T (bmy, 1/14/05)
!******************************************************************************
!     
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "global_hno3_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE :: INIT_GLOBAL_NO3

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Array to store global monthly mean OH field
      REAL*8, ALLOCATABLE :: NO3(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE GET_GLOBAL_NO3( THISMONTH )
!
!******************************************************************************
!  Subroutine GET_GLOBAL_NO3 reads monthly mean NO3 data fields.  These 
!  are needed for simulations such as offline sulfate/aerosol. 
!  (bmy, 10/15/02, 1/13/05)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) Minor bug fix in FORMAT statements (bmy, 3/23/03)
!  (2 ) Cosmetic changes (bmy, 3/27/03)
!  (3 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (4 ) Now suppress output from READ_BPCH2 with QUIET=T (bmy, 1/14/05)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D

      IMPLICIT NONE

#     include "CMN_SIZE"   ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)  :: THISMONTH

      ! Local variables
      INTEGER              :: I, J, L
      REAL*4               :: ARRAY(IGLOB,JGLOB,LGLOB)
      REAL*8               :: XTAU
      CHARACTER(LEN=255)   :: FILENAME

      ! First time flag
      LOGICAL, SAVE        :: FIRST = .TRUE. 

      !=================================================================
      ! GET_GLOBAL_NO3 begins here!
      !=================================================================

      ! Allocate NO3 array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_NO3
         FIRST = .FALSE.
      ENDIF

      ! File name
      FILENAME = TRIM( DATA_DIR ) // 'sulfate_sim_200210/NO3.' //
     &           GET_NAME_EXT()  // '.' // GET_RES_EXT()

      ! Echo some information to the standard output
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_NO3: Reading ', a )

      ! Get the TAU0 value for the start of the given month
      ! Assume "generic" year 1985 (TAU0 = [0, 744, ... 8016])
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )
 
      ! Read NO3 data from the binary punch file (tracer #5)
      CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 5,     
     &                 XTAU,      IGLOB,     JGLOB,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable H2O2
      CALL TRANSFER_3D( ARRAY, NO3 )

      ! Return to calling program
      END SUBROUTINE GET_GLOBAL_NO3

!------------------------------------------------------------------------------

      SUBROUTINE INIT_GLOBAL_NO3
!
!******************************************************************************
!  Subroutine INIT_GLOBAL_NO3 allocates the NO3 module array (bmy, 10/15/02)
!
!  NOTES:
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER             :: AS

      !=================================================================
      ! INIT_GLOBAL_H2O2 begins here!
      !=================================================================
      ALLOCATE( NO3( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NO3' )
      NO3 = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_GLOBAL_NO3
      
!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_GLOBAL_NO3
!
!******************************************************************************
!  Subroutine CLEANUP_GLOBAL_H2O2 deallocates the H2O2 array. (bmy, 10/15/02)
!******************************************************************************
!                               
      !=================================================================
      ! CLEANUP_GLOBAL_H2O2 begins here!
      !=================================================================
      IF ( ALLOCATED( NO3 ) ) DEALLOCATE( NO3 ) 
     
      ! Return to calling program
      END SUBROUTINE CLEANUP_GLOBAL_NO3

!------------------------------------------------------------------------------

      END MODULE GLOBAL_NO3_MOD
