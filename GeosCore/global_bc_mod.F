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
! !REMARKS:
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
      REAL*8, ALLOCATABLE, TARGET :: BC(:,:,:)

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
      SUBROUTINE GET_GLOBAL_BC( am_I_Root, Input_Opt, 
     &                          THISMONTH, THISYEAR, RC )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(IN)  :: THISMONTH
      INTEGER,        INTENT(IN)  :: THISYEAR
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?!
! !INPUT PARAMETERS: 
!
!
! !REVISION HISTORY: 
!  19 January 2011 - C.L. Friedman - Initial Version
!  20 Jun 2014 - R. Yantosca - Now pass am_I_Root, Input_Opt, RC via arg list
!  05 Sep 2014 - R. Yantosca - Add hack to read GCAP file names
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

      ! Assume success
      RC = GIGC_SUCCESS

      ! Allocate BC array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_BC
         FIRST = .FALSE.
      ENDIF

#if defined( GCAP )

      !%%% NOTE: For some reason the BC/OC files for GCAP are different
      !%%% than for the other met fields.  Put in a hack to allow unit
      !%%% testing to progress.  This will be a moot point with the
      !%%% upcoming implementation of HEMCO. (bmy, 9/5/14)

      ! Filename
      FILENAME = TRIM( Input_Opt%DATA_DIR ) // 
     &           'POPs_201209/BCPO_GCAP.bpch'

#else

      ! Filename
      FILENAME = TRIM( Input_Opt%DATA_DIR ) // 'POPs_201209/BCPO_FC_' //
     &           GET_RES_EXT()              // '.bpch'

#endif

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

      ! Zero BC array
      BC = 0d0

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

      END SUBROUTINE CLEANUP_GLOBAL_BC
!EOC      
      END MODULE GLOBAL_BC_MOD
