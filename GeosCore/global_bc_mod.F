!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_bc_mod
!
! !DESCRIPTION:  Module GLOBAL\_BC\_MOD contains variables and routines for
!  reading the global monthly mean OC concentration from disk. Based on module
!  GLOBAL\_OH\_MOD. (clf, 1/19/2011).
!\\
!\\
! !INTERFACE:
!
      MODULE GLOBAL_BC_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: GET_GLOBAL_BC
      PUBLIC :: CLEANUP_GLOBAL_BC
!
! !PUBLIC DATA MEMBERS:
!
      PUBLIC :: BC
!
! !REVISION HISTORY:
!  19 January 2011 - C.L. Friedman - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Array to store global monthly mean BC field
      REAL*8, ALLOCATABLE :: BC(:,:,:)
!      REAL*8, ALLOCATABLE :: BCPHIL(:,:,:)
!      REAL*8, ALLOCATABLE :: BCPHOB(:,:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_global_bc
!
! !DESCRIPTION: Subroutine GET\_GLOBAL\_BC reads global BC from binary punch
!  files stored on disk. BC data is needed for partitioning of gas phase
!  organics onto BC particles (e.g., POPs). (clf, 1/19/2011)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GLOBAL_BC( THISMONTH, THISYEAR )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D

      USE CMN_SIZE_MOD        ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: THISMONTH, THISYEAR
!
! !REVISION HISTORY: 
!  19 January 2011 - C.L. Friedman - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER              :: I, J, L
      INTEGER              :: YEAR
      REAL*4               :: ARRAY(IGLOB,JGLOB,LGLOB)
      REAL*8               :: XTAU
      CHARACTER(LEN=255)   :: FILENAME

      ! First time flag
      LOGICAL, SAVE        :: FIRST = .TRUE. 

      !=================================================================
      ! GET_GLOBAL_BC begins here!
      !=================================================================

      ! Allocate BC array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_BC
         FIRST = .FALSE.
      ENDIF

      ! Filename
      FILENAME = TRIM( DATA_DIR ) // 'POPs_201209/BCPO_FC_' //
     &           GET_RES_EXT()    // '.bpch'

      ! Echo some information to the standard output
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_BC: Reading BC from: ', a )

      ! Data is only available for 2005-2009
      IF ( THISYEAR < 2005 ) THEN
         YEAR = 2005
      ELSE IF ( THIS YEAR > 2009 ) THEN
         YEAR = 2009
      ELSE
         YEAR = THISYEAR
      ENDIF

      ! Get the TAU0 value for the start of the given month
      XTAU = GET_TAU0( THISMONTH, 1, YEAR )

      ! Read in BC data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'IJ-24H-$', 13,     
     &                 XTAU,      IGLOB,     JGLOB,      
     &                 LGLOB,     ARRAY,     QUIET=.FALSE. )

      ! Assign data from ARRAY to the module variable BC
      CALL TRANSFER_3D( ARRAY, BC )

      END SUBROUTINE GET_GLOBAL_BC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  init_global_bc
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_BC allocates and zeroes the BC array, 
! which holds global monthly mean BC concentrations. (clf, 1/19/2011)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GLOBAL_BC
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD
!
! !REVISION HISTORY: 
!  19 January 2011 - C.L. Friedman - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_GLOBAL_BC begins here!
      !=================================================================

      ! Allocate BC array
      ALLOCATE( BC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BC' )

!      ALLOCATE( BCPHIL( IIPAR, JJPAR, LLPAR ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCPHIL' )

!      ALLOCATE( BCHPOB( IIPAR, JJPAR, LLPAR ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'BCPHOB' )

      ! Zero BC array
      BC = 0d0
!      BCPHIL = 0d0
!      BCPHOB = 0d0

      END SUBROUTINE INIT_GLOBAL_BC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  cleanup_global_bc
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_BC deallocates the BC array.
!  (clf, 1/19/2011)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GLOBAL_BC
!
! !REVISION HISTORY: 
!  19 January 2011 - C.L. Friedman - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! CLEANUP_GLOBAL_BC begins here!
      !=================================================================
      IF ( ALLOCATED( BC ) ) DEALLOCATE( BC ) 
!      IF ( ALLOCATED( BCPHIL ) ) DEALLOCATE( BCPHIL ) 
!      IF ( ALLOCATED( BCPHOB ) ) DEALLOCATE( BCPHOB )

      END SUBROUTINE CLEANUP_GLOBAL_BC
!EOC      
      END MODULE GLOBAL_BC_MOD
