! $Id: global_o1d_mod.f,v 1.1 2007/11/05 16:16:19 bmy Exp $
      MODULE GLOBAL_O1D_MOD
!
!******************************************************************************
!  Module GLOBAL_O1D_MOD contains variables and routines for reading the
!  global monthly mean O1D stratospheric concentration from disk.  This is 
!  used in the H2/HD simulation.  The O1D fields were obtained from Gabriele 
!  Curci GEOS-Chem simulation in the stratosphere (v5.03).
!  (hup, phs, 9/18/07)
!
!  Module Variables:
!  ===========================================================================
!  (1 ) O1D (REAL*8)       : stores global monthly mean O1D field
!  
!  Module Routines:
!  =========================================================================== 
!  (1 ) GET_O1D            : Wrapper for GET_GLOBAL_O1D
!  (2 ) GET_GLOBAL_O1D     : Reads global monthly mean O1D from disk
!  (3 ) INIT_GLOBAL_O1D    : Allocates & initializes the O1D array
!  (4 ) CLEANUP_GLOBAL_O1D : Deallocates the OH array
!
!  GEOS-Chem modules referenced by global_o1d_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f  : Module containing routines for binary punch file I/O
!  (2 ) error_mod.f  : Module containing NaN and other error-check routines
!
!  NOTES:
!  (1 ) Adapted from GLOBAL_OH_MOD module (hup, phs, 9/18/07)
!******************************************************************************
!     
      IMPLICIT NONE

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Array to store global monthly mean O1D field
      REAL*8, ALLOCATABLE :: O1D(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE GET_GLOBAL_O1D( THISMONTH )
!
!******************************************************************************
!  Subroutine GET_GLOBAL_O1D reads global O1D from binary punch files stored
!  in the /data/ctm/GEOS_MEAN directory.  This O1D data is needed for the H2/HD
!  mechanisms in Tagged H2. (hup, phs, 9/18/07)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) THISMONTH (INTEGER) : Current month number (1-12)
!
!  NOTES:
!  (1 ) GET_GLOBAL_O1D assumes that we are reading global O1D data that 
!        occupies all CTM levels.  Contact Bob Yantosca (bmy@io.harvard.edu) 
!        for IDL regridding code which will produce the appropriate O1D files.
!  (2 ) ARRAY should now be of size (IGLOB,JGLOB,LGLOB). (bmy, 1/11/02)
!  (3 ) Now point to new O1D files in the ??? subdirectory.
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D

      IMPLICIT NONE

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: THISMONTH

      ! Local variables
      INTEGER                :: I, J, L
      REAL*4                 :: ARRAY(IGLOB,JGLOB,LGLOB)
      REAL*8                 :: XTAU
      CHARACTER(LEN=255)     :: FILENAME

      ! First time flag
      LOGICAL, SAVE          :: FIRST = .TRUE. 

      !=================================================================
      ! GET_GLOBAL_O1D begins here!
      !=================================================================

      ! Allocate O1D array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_O1D
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Read Gabriele Curci's O1D (v5.03)
      !=================================================================

      FILENAME = TRIM( DATA_DIR ) // 'hydrogen_200704/stratO1D.' //
     &           GET_NAME_EXT()  // '.' // GET_RES_EXT()

      ! Echo some information to the standard output
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_O1D: Reading O1D from: ', a )

      ! Get the TAU0 value for the start of the given month
      ! Assume "generic" year 1998
      XTAU = GET_TAU0( THISMONTH, 1, 1998 )

      ! Read O1D data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'SL-AVG-$', 2,  
     &                 XTAU,      IGLOB,     JGLOB,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )
      
      ! Assign data from ARRAY to the module variable O1D
      CALL TRANSFER_3D( ARRAY, O1D )

      ! Return to calling program
      END SUBROUTINE GET_GLOBAL_O1D

!------------------------------------------------------------------------------

      SUBROUTINE INIT_GLOBAL_O1D
!
!******************************************************************************
!  Subroutine INIT_GLOBAL_O1D allocates and zeroes the O1D array, which holds 
!  global monthly mean O1D concentrations. 
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE" 

      ! Local variables
      INTEGER :: AS

      ! Allocate O1D array
      ALLOCATE( O1D( IGLOB, JGLOB, LGLOB ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'O1D' )

      ! Zero O1D array
      O1D = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_GLOBAL_O1D
      
!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_GLOBAL_O1D
!
!******************************************************************************
!  Subroutine CLEANUP_GLOBAL_O1D deallocates the O1D array.
!******************************************************************************
!                               
      IF ( ALLOCATED( O1D ) ) DEALLOCATE( O1D ) 
     
      ! Return to calling program
      END SUBROUTINE CLEANUP_GLOBAL_O1D

!------------------------------------------------------------------------------

      ! End of module
      END MODULE GLOBAL_O1D_MOD
