!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_oh_mod
!
! !DESCRIPTION: Module GLOBAL\_OH\_MOD contains variables and routines for 
!  reading the global monthly mean OH concentration from disk. 
!\\
!\\
! !INTERFACE: 
!
      MODULE GLOBAL_OH_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      ! Array to store global monthly mean OH field [molec/cm3]
      REAL*8, PUBLIC, ALLOCATABLE :: OH(:,:,:)
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CLEANUP_GLOBAL_OH
      PUBLIC :: GET_GLOBAL_OH
      PUBLIC :: INIT_GLOBAL_OH
!
! !REVISION HISTORY:
!  28 Jul 2000 - R. Yantosca - Initial version
!  (1 ) Updated comments (bmy, 9/4/01)
!  (2 ) Now use routines from "transfer_mod.f" to regrid OH to 30 levels
!        for reduced GEOS-3 grid.  Also size OH array properly. (bmy, 1/14/02)
!  (3 ) Eliminate obsolete code from 11/01 (bmy, 2/27/02)
!  (4 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections.  Updated comments (bmy, 5/28/02)
!  (5 ) Now use updated OH fields (bmy, 10/2/02)
!  (6 ) Now references "error_mod.f" (bmy, 10/15/02)
!  (7 ) Minor bug fixes in FORMAT statements (bmy, 3/23/03)
!  (8 ) Cosmetic changes to simplify output (bmy, 3/27/03)
!  (9 ) Bug fix: OH should be (IIPAR,JJPAR,LLPAR) (bmy, 5/4/04)
!  (10) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
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
! !IROUTINE: get_global_oh
!
! !DESCRIPTION: Subroutine GET\_GLOBAL\_OH reads global OH from binary 
!  punch files stored in the /data/ctm/GEOS\_MEAN directory.  This OH data 
!  is needed as oxidant for various offline chemistry mechanisms.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GLOBAL_OH( am_I_Root, Input_Opt, 
     &                          State_Met, THISMONTH, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE BPCH2_MOD,            ONLY : GET_NAME_EXT
      USE BPCH2_MOD,            ONLY : GET_RES_EXT
      USE BPCH2_MOD,            ONLY : GET_TAU0
      USE BPCH2_MOD,            ONLY : READ_BPCH2
      USE DAO_MOD,              ONLY : AIRDEN_FULLGRID
      USE DAO_MOD,              ONLY : AIRQNT_FULLGRID
      USE DIRECTORY_MOD,        ONLY : OH_DIR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod,   ONLY : OptInput
      USE GIGC_State_Met_Mod,   ONLY : MetState
      USE TRANSFER_MOD,         ONLY : TRANSFER_3D
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Obj for input options
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Obj for met fields
      INTEGER,        INTENT(IN)    :: THISMONTH   ! Current month
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  28 Jul 2000 - R. Yantosca - Initial version
!  (1 ) GET_GLOBAL_OH assumes that we are reading global OH data that occupies
!        all CTM levels.  Contact Bob Yantosca (bmy@io.harvard.edu) for IDL
!        regridding code which will produce the appropriate OH files.
!  (2 ) Now use version of GET_TAU0 with 3 arguments.  Now call READ_BPCH2
!        with IIPAR,JJPAR,LGLOB.  Call TRANSFER_3D to cast from REAL*4 to
!        REAL*8 and to regrid to 30 levels for GEOS-3 (if necessary).
!        ARRAY should now be of size (IIPAR,JJPAR,LGLOB). (bmy, 1/11/02)
!  (3 ) Now point to new OH files in the v4-26 subdirectory.  Also eliminated
!        obsolete code from 11/01. (bmy, 2/27/02)
!  (4 ) Now point to OH files in the v4-33 subdirectory. (bmy, 10/2/02)
!  (5 ) Replace missing commas in the FORMAT statement (bmy, 3/23/03)
!  (6 ) Cosmetic changes to simplify output (bmy, 3/27/03)
!  (7 ) Add Mat's OH as an option.  Also read bpch file quietly (bmy, 5/4/04)
!  (8 ) Now use OH_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!  23 Oct 2013 - R. Yantosca - Now accept Input_Opt argument
!  29 Oct 2013 - R. Yantosca - Remove TRANSFER_3D_NOLUMP routine, we can just
!                              instead do a direct cast assignment
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
      ! GET_GLOBAL_OH begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Allocate OH array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_OH
         FIRST = .FALSE.
      ENDIF

#if defined( NESTED_NA )
      ! Filename
      FILENAME = TRIM( OH_DIR ) // 'OH_3Dglobal.' // GET_NAME_EXT() // 
     &           '.'            // GET_RES_EXT()  // '_NA'
#else
      ! Filename
      FILENAME = TRIM( OH_DIR ) // 'OH_3Dglobal.' // GET_NAME_EXT() // 
     &                              '.'           // GET_RES_EXT()
#endif

      ! Echo some information to the standard output
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_OH: Reading OH from: ', a )

      ! Get the TAU0 value for the start of the given month
      ! Assume "generic" year 1985 (TAU0 = [0, 744, ... 8016])
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

      ! Read OH data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'CHEM-L=$', 1,     
     &                 XTAU,      IIPAR,     JJPAR,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

#if defined( ESMF_ ) || defined( EXTERNAL_GRID ) || defined( EXTERNAL_FORCING )
      !-----------------------------------------------------------------
      !     %%%%% CONNECTING TO GEOS-5 GCM via ESMF INTERFACE %%%%%
      !
      ! Do not use routines AIRQNT_FULLGRID or TRANSFER_3D_NO_LUMP when
      ! we are connecting to the GEOS-5 GCM via ESMF. (mpayer, 3/14/13)
      !-----------------------------------------------------------------

      ! Assign data from ARRAY to the module variable OH
      CALL TRANSFER_3D( ARRAY, OH )

#else
      !-----------------------------------------------------------------
      !                %%%%% TRADITIONAL GEOS-Chem %%%%%
      !
      ! Include fix for vertical regridding of OH for offline simulations
      ! provided by Chris Holmes (mpayer, 3/14/13)
      !-----------------------------------------------------------------
#if defined( GRIDREDUCED )

      ! Calculate the number density of air on the full vertical grid
      CALL AIRQNT_FULLGRID( am_I_Root, Input_Opt, RC )

      ! Convert molec/cm3 -> mixing ratio (molec m3 kg-1 cm-3)
      ! because regridding requires mixing ratio quantity
      ! Note that peculiar units don't matter because we immediately
      ! convert back after regridding
      DO L=1, LGLOB
         ARRAY(:,:,L) = ARRAY(:,:,L) / AIRDEN_FULLGRID(L,:,:)
      ENDDO
         
      ! Assign data from ARRAY to the module variable OH
      CALL TRANSFER_3D( ARRAY, OH )
      
      ! convert mixing ratio (molec m3 kg-1 cm-3) -> molec/cm3
      DO L=1, LLPAR
         OH(:,:,L) = OH(:,:,L) * State_Met%AIRDEN(L,:,:)
      ENDDO

#else

      ! Assign data from ARRAY to the module variable OH
      OH = ARRAY

#endif
#endif

      END SUBROUTINE GET_GLOBAL_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_global_oh
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_OH allocates and zeroes
!  all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GLOBAL_OH
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD 
! 
! !REVISION HISTORY: 
!  28 Jul 2000 - R. Yantosca - Initial version
!  (1 ) OH array now needs to be sized (IIPAR,JJPAR,LGLOB) (bmy, 1/14/02)
!  (2 ) Also eliminated obsolete code from 11/01 (bmy, 2/27/02)
!  (3 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (4 ) OH should be (IIPAR,JJPAR,LLPAR): avoid subscript errors (bmy, 5/4/04)
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS

      !=================================================================
      ! INIT_GLOBAL_OH begins here!
      !=================================================================

      ! Allocate OH array
      ALLOCATE( OH( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'OH' )

      ! Zero OH array
      OH = 0d0

      END SUBROUTINE INIT_GLOBAL_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_global_oh
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_OH deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GLOBAL_OH
! 
! !REVISION HISTORY: 
!  28 Jul 2000 - R. Yantosca - Initial version
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_GLOBAL_OH begins here!
      !=================================================================
      IF ( ALLOCATED( OH ) ) DEALLOCATE( OH ) 
     
      END SUBROUTINE CLEANUP_GLOBAL_OH
!EOC
      END MODULE GLOBAL_OH_MOD
