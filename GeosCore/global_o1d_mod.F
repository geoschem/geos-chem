!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_o1d_mod
!
! !DESCRIPTION: Module GLOBAL\_O1D\_MOD contains variables and routines for 
!  reading the global monthly mean O1D stratospheric concentration from disk.  
!  This is used in the H2/HD simulation.  The O1D fields were obtained from 
!  Gabriele Curci GEOS-Chem simulation in the stratosphere (v5.03).
!\\
!\\
! !INTERFACE: 
!
      MODULE GLOBAL_O1D_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      ! Array to store global monthly mean O1D field
      REAL*8, PUBLIC, ALLOCATABLE :: O1D(:,:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_GLOBAL_O1D
      PUBLIC :: GET_GLOBAL_O1D
      PUBLIC :: INIT_GLOBAL_O1D
!
! !REVISION HISTORY:
!  18 Sep 2007 - H. U. Price, P. Le Sager - Initial version
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_global_O1D 
!
! !DESCRIPTION: Subroutine GET\_GLOBAL\_O1D reads global O1D from binary
!  punch files stored in the /data/ctm/GEOS\_MEAN directory.  This O1D data 
!  is needed for the H2/HD mechanisms in Tagged H2.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GLOBAL_O1D( THISMONTH )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0
      USE BPCH2_MOD,     ONLY : READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D

      USE CMN_SIZE_MOD                  ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: THISMONTH   ! Current month
! 
! !REVISION HISTORY: 
!  18 Sep 2007 - H. U. Price, P. Le Sager - Initial version
!  (1 ) GET_GLOBAL_O1D assumes that we are reading global O1D data that 
!        occupies all CTM levels.  Contact Bob Yantosca (bmy@io.harvard.edu) 
!        for IDL regridding code which will produce the appropriate O1D files.
!  (2 ) ARRAY should now be of size (IIPAR,JJPAR,LGLOB). (bmy, 1/11/02)
!  (3 ) Now point to new O1D files in the ??? subdirectory.
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, L
      REAL*4             :: ARRAY(IIPAR,JJPAR,LGLOB)
      REAL*8             :: XTAU
      CHARACTER(LEN=255) :: FILENAME

      ! First time flag
      LOGICAL, SAVE      :: FIRST = .TRUE. 

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
     &                 XTAU,      IIPAR,     JJPAR,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )
      
      ! Assign data from ARRAY to the module variable O1D
      CALL TRANSFER_3D( ARRAY, O1D )

      END SUBROUTINE GET_GLOBAL_O1D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_global_o1d
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_O1D allocates and zeroes
!  all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GLOBAL_O1D
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD 
! 
! !REVISION HISTORY: 
!  18 Sep 2007 - H. U. Price, P. Le Sager - Initial version
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      ! Allocate O1D array
      ALLOCATE( O1D( IIPAR, JJPAR, LGLOB ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'O1D' )

      ! Zero O1D array
      O1D = 0d0

      END SUBROUTINE INIT_GLOBAL_O1D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_global_O1D
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_O1D deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GLOBAL_O1D
! 
! !REVISION HISTORY: 
!  18 Sep 2007 - H. U. Price, P. Le Sager - Initial version
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC      
      IF ( ALLOCATED( O1D ) ) DEALLOCATE( O1D ) 
     
      END SUBROUTINE CLEANUP_GLOBAL_O1D
!EOC
      END MODULE GLOBAL_O1D_MOD
