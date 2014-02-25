!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_oc_mod
!
! !DESCRIPTION:  Module GLOBAL\_OC\_MOD contains variables and routines for
!  reading the global monthly mean OC concentration from disk. Based on module
!  GLOBAL\_OH\_MOD. (clf, 1/19/2011).
!\\
!\\
! !INTERFACE:
!
      MODULE GLOBAL_OC_MOD
! 
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: GET_GLOBAL_OC
      PUBLIC :: CLEANUP_GLOBAL_OC
!
! !PUBLIC DATA MEMBERS:
!
      PUBLIC :: OC
!
! !REVISION HISTORY:
!  19 January 2011 - C.L. Friedman - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Array to store global monthly mean OC field
      REAL*8, ALLOCATABLE :: OC(:,:,:)
!      REAL*8, ALLOCATABLE :: OCPHIL(:,:,:)
!      REAL*8, ALLOCATABLE :: OCPHOB(:,:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_global_oc
!
! !DESCRIPTION: GET\_GLOBAL\_OC reads global OC from binary punch files stored
!  on disk. OC data is needed for partitioning of gas phase organics into OC
!  particles (e.g., POPs). (clf, 1/19/2011)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GLOBAL_OC( THISMONTH, THISYEAR )
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
      ! GET_GLOBAL_OC begins here!
      !=================================================================

      ! Allocate OC array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_OC
         FIRST = .FALSE.
      ENDIF

      ! Filename
      FILENAME = TRIM( DATA_DIR ) // 'POPs_201209/OCPO_FC_' //
     &           GET_RES_EXT()    // '.bpch'

      ! Echo some information to the standard output
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_OC: Reading OC from: ', a )

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

      ! Read in OC data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'IJ-24H-$', 14,     
     &                 XTAU,      IGLOB,     JGLOB,      
     &                 LGLOB,     ARRAY,     QUIET=.FALSE. )

      ! Assign data from ARRAY to the module variable OC
      CALL TRANSFER_3D( ARRAY, OC )

      END SUBROUTINE GET_GLOBAL_OC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  init_global_oc
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_OC allocates and zeroes the OC array, 
! which holds global monthly mean OC concentrations. (clf, 1/19/2011)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GLOBAL_OC
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
      ! INIT_GLOBAL_OC begins here!
      !=================================================================

      ! Allocate OC array
      ALLOCATE( OC( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OC' )

!      ALLOCATE( OCPHIL( IIPAR, JJPAR, LLPAR ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCPHIL' )

!      ALLOCATE( OCHPOB( IIPAR, JJPAR, LLPAR ), STAT=AS )
!      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OCPHOB' )

      ! Zero OC arrays
      OC = 0d0
!      OCPHIL = 0d0
!      OCPHOB = 0d0

      END SUBROUTINE INIT_GLOBAL_OC    
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  cleanup_global_oc
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_OC deallocates the OC array.
!  (clf, 1/19/2011)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GLOBAL_OC
!
! !REVISION HISTORY: 
!  19 January 2011 - C.L. Friedman - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! CLEANUP_GLOBAL_OC begins here!
      !=================================================================
      IF ( ALLOCATED( OC ) ) DEALLOCATE( OC ) 
!      IF ( ALLOCATED( OCPHIL ) ) DEALLOCATE( OCPHIL ) 
!      IF ( ALLOCATED( OCPHOB ) ) DEALLOCATE( OCPHOB ) 
     
      END SUBROUTINE CLEANUP_GLOBAL_OC
!EOC
      END MODULE GLOBAL_OC_MOD
