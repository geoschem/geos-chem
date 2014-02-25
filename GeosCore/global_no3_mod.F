!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_NO3_mod
!
! !DESCRIPTION: Module GLOBAL\_NO3\_MOD contains variables and routines for 
!  reading the global monthly mean NO3 concentration from disk.  These are 
!  needed for the offline sulfate/aerosol simulation.
!\\
!\\
! !INTERFACE: 
!
      MODULE GLOBAL_NO3_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      ! Array to store global monthly mean OH field
      REAL*8, PUBLIC, ALLOCATABLE :: NO3(:,:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: GET_GLOBAL_NO3
      PUBLIC  :: CLEANUP_GLOBAL_NO3
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE :: INIT_GLOBAL_NO3
!
! !REVISION HISTORY:
!  15 Oct 2002 - R. Yantosca - Initial version
!  (1 ) Adapted from "global_oh_mod.f" (bmy, 10/3/02)
!  (2 ) Minor bug fix in FORMAT statements (bmy, 3/23/03)
!  (3 ) Cosmetic changes (bmy, 3/27/03)
!  (4 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (5 ) Now suppress output from READ_BPCH2 with QUIET=T (bmy, 1/14/05)
!  (6 ) Now read from "sulfate_sim_200508/offline" directory (bmy, 8/1/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Bug fix: now zero ARRAY (phs, 1/22/07)
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
! !IROUTINE: get_global_NO3
!
! !DESCRIPTION: Subroutine GET\_GLOBAL\_NO3 reads monthly mean NO3 data fields.
!  These are needed for simulations such as offline sulfate/aerosol. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GLOBAL_NO3( THISMONTH )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0
      USE BPCH2_MOD,     ONLY : READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D_TROP

      USE CMN_SIZE_MOD                  ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: THISMONTH   ! Current month
! 
! !REVISION HISTORY: 
!  15 Oct 2002 - R. Yantosca - Initial version
!  (1 ) Minor bug fix in FORMAT statements (bmy, 3/23/03)
!  (2 ) Cosmetic changes (bmy, 3/27/03)
!  (3 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (4 ) Now suppress output from READ_BPCH2 with QUIET=T (bmy, 1/14/05)
!  (5 ) GEOS-3 & GEOS-4 data comes from model runs w/ 30 levels.  Also now 
!        read from "sulfate_sim_200508/offline" directory.  Also now read
!        up to LLTROP levels.  Now reference TRANSFER_3D_TROP from 
!        "transfer_mod.f". (bmy, 8/1/05)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Now zero local variable ARRAY (phs, 1/22/07)
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4             :: ARRAY(IIPAR,JJPAR,LLTROP)
      REAL*8             :: XTAU
      CHARACTER(LEN=255) :: FILENAME

      ! First time flag
      LOGICAL, SAVE      :: FIRST = .TRUE. 

      !=================================================================
      ! GET_GLOBAL_NO3 begins here!
      !=================================================================

      ! Allocate NO3 array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_NO3
         FIRST = .FALSE.
      ENDIF

      ! File name
      FILENAME = TRIM( DATA_DIR )                       // 
     &           'sulfate_sim_200508/offline/NO3.'      //
     &           GET_NAME_EXT() // '.' // GET_RES_EXT()

      ! Echo some information to the standard output
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_NO3: Reading ', a )

      ! Get the TAU0 value for the start of the given month
      ! Assume "generic" year 1985 (TAU0 = [0, 744, ... 8016])
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

      ! Zero ARRAY so that we avoid random data between 
      ! levels LLTROP_FIX and LLTROP (phs, 1/22/07)
      ARRAY = 0e0
 
      ! Read NO3 data from the binary punch file (tracer #5)
      ! NOTE: NO3 data is only defined w/in the tropopause, so set the 3rd
      ! dim of ARRAY to LLTROP_FIX (i.e, case of annual mean tropopause). 
      ! This is backward compatibility with offline data set. (phs, 1/22/07)
      CALL READ_BPCH2( 
     &         FILENAME,   'CHEM-L=$',                5,     
     &         XTAU,        IIPAR,                    JJPAR,      
     &         LLTROP_FIX,  ARRAY(:,:,1:LLTROP_FIX),  QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable H2O2
      ! Levels between LLTROP_FIX and LLROP are 0
      CALL TRANSFER_3D_TROP( ARRAY, NO3 )

      END SUBROUTINE GET_GLOBAL_NO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_global_NO3
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_NO3 allocates and zeroes
!  all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GLOBAL_NO3
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD 
! 
! !REVISION HISTORY: 
!  15 Oct 2002 - R. Yantosca - Initial version
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Now allocate NO3 array up to LLTROP levels (bmy, 8/1/05)
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_GLOBAL_NO3 begins here!
      !=================================================================
      ALLOCATE( NO3( IIPAR, JJPAR, LLTROP ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'NO3' )
      NO3 = 0d0

      END SUBROUTINE INIT_GLOBAL_NO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_global_no3
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_NO3 deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GLOBAL_NO3
! 
! !REVISION HISTORY: 
!  15 Oct 2002 - R. Yantosca - Initial version
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC      
      !=================================================================
      ! CLEANUP_GLOBAL_H2O2 begins here!
      !=================================================================
      IF ( ALLOCATED( NO3 ) ) DEALLOCATE( NO3 ) 
     
      END SUBROUTINE CLEANUP_GLOBAL_NO3
!EOC
      END MODULE GLOBAL_NO3_MOD
