!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_o3_mod
!
! !DESCRIPTION: Module GLOBAL\_O3\_MOD contains variables and routines for 
!  reading the global monthly mean O3 concentration from disk.  These are 
!  needed for the offline sulfate/aerosol simulation.
!\\
!\\
! !INTERFACE: 
!
      MODULE GLOBAL_O3_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC              :: CLEANUP_GLOBAL_O3
      PUBLIC              :: GET_GLOBAL_O3
!
! !PUBLIC DATA MEMBERS:
!            
      PUBLIC              :: O3
      REAL*8, ALLOCATABLE :: O3(:,:,:)        ! Global monthly mean OH field
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE             :: INIT_GLOBAL_O3
!
! !REVISION HISTORY:
!  (1 ) Now references "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now reads O3 data from "sulfate_sim_200508/offline" dir (bmy, 8/30/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Bug fixes in GET_GLOBAL_O3 (bmy, 12/1/05)
!  (5 ) Now reads O3 from MERGE files, which include stratospheric O3 from
!        COMBO, for GEOS-3 and GEOS-4 met fields (phs, 1/19/07)
!  (6 ) Bug fix in GET_GLOBAL_O3 (bmy, 1/14/09)
!  13 Aug 2010 - R. Yantosca - Added modifications for MERRA
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  08 Feb 2012 - R. Yantosca - Add modifications for GEOS-5.7.x
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
! !IROUTINE: get_global_o3
!
! !DESCRIPTION: Subroutine GET\_GLOBAL\_O3 reads monthly mean O3 data fields.  
!  These are needed for simulations such as offline sulfate/aerosol. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GLOBAL_O3( THISMONTH )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_3D

      USE CMN_SIZE_MOD                      ! Size parameters
      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)  :: THISMONTH       ! Current month
!
! !REVISION HISTORY: 
!  23 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Minor bug fix in FORMAT statements (bmy, 3/23/03)
!  (2 ) Cosmetic changes (bmy, 3/27/03)
!  (3 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (4 ) Now reads O3 data from "sulfate_sim_200508/offline" dir (bmy, 8/30/05)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (6 ) Tracer number for O3 is now 51.  Also need to call TRANSFER_3D_TROP
!        since the new O3 data file only goes up to LLTROP. (bmy, 11/18/05)
!  (7 ) Modified to include stratospheric O3 -- Requires access to new
!        MERGE.O3* files. (phs, 1/19/07)
!  (8 ) Renamed GRID30LEV to GRIDREDUCED (bmy, 2/7/07)
!  (9 ) Bug fix: don't call TRANSFER_3D if you use GRIDREDUCED (bmy, 1/14/09)
!  13 Aug 2010 - R. Yantosca - Rewrote logic more cleanly
!  13 Aug 2010 - R. Yantosca - Treat MERRA in same way as GEOS-5 
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!  19 Aug 2010 - R. Yantosca - Removed hardwiring of data directory
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*4               :: ARRAY(IIPAR,JJPAR,LGLOB)
      REAL*4               :: ARRAY2(IIPAR,JJPAR,LLPAR)
      REAL*8               :: XTAU
      CHARACTER(LEN=255)   :: FILENAME

      ! First time flag
      LOGICAL, SAVE        :: FIRST = .TRUE. 

      !=================================================================
      ! GET_GLOBAL_O3 begins here!
      !=================================================================

      ! Allocate O3 array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_O3
         FIRST = .FALSE.
      ENDIF

#if   defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

      !----------------------------------------------------------------
      ! GEOS-5 or MERRA
      !----------------------------------------------------------------
#if   defined( GRIDREDUCED )

      ! Filename for 47-level model
#if defined( GRID05x0666 ) && defined( NESTED_NA )
      ! Need to handle the nested-resolution file, Y. Zhang 11/1/2011
      FILENAME = TRIM( DATA_DIR )                           // 
     &           'mercury_201203/MERGE.O3.47L.'             //
     &           GET_NAME_EXT() // '.' // GET_RES_EXT()
#else
      FILENAME = TRIM( DATA_DIR )                           // 
     &           'sulfate_sim_200508/offline/MERGE.O3.47L.' //
     &           GET_NAME_EXT() // '.' // GET_RES_EXT()
#endif

#else

      ! Filename for 72-level model
      FILENAME = TRIM( DATA_DIR )                           // 
     &           'sulfate_sim_200508/offline/MERGE.O3.'     //
     &           GET_NAME_EXT() // '.' // GET_RES_EXT()

#endif


#else

      !----------------------------------------------------------------
      ! GEOS-4 or GEOS-3
      !----------------------------------------------------------------
#if   defined( GRIDREDUCED )

      ! Filename for 30-level model
      FILENAME = TRIM( DATA_DIR )                           // 
     &           'sulfate_sim_200508/offline/MERGE.O3.30L.' // 
     &           GET_NAME_EXT() // '.' // GET_RES_EXT()


#else

      ! Filename for full vertical grid
      FILENAME = TRIM( DATA_DIR )                           // 
     &           'sulfate_sim_200508/offline/MERGE.O3.'     // 
     &           GET_NAME_EXT() // '.' // GET_RES_EXT()

#endif

#endif


      ! Echo some information to the standard output
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_O3: Reading ', a )

      ! Get the TAU0 value for the start of the given month
      ! Assume "generic" year 1985 (TAU0 = [0, 744, ... 8016])
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )
 
#if   defined( GRIDREDUCED )
 
      ! Read O3 data (v/v) from the binary punch file (tracer #51)
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 51,     
     &                 XTAU,      IIPAR,     JJPAR,      
     &                 LLPAR,     ARRAY2,    QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable O3
      ! (don't have to fold layers in the stratosphere)
      O3 = ARRAY2

#else

      ! Read O3 data (v/v) from the binary punch file (tracer #51)
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 51,     
     &                 XTAU,      IIPAR,     JJPAR,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable O3
      ! (folding layers in the stratosphere)
      CALL TRANSFER_3D( ARRAY, O3 )

#endif

      ! Return to calling program
      END SUBROUTINE GET_GLOBAL_O3
!EOC      
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_global_o3
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_O3 allocates the O3 module array.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GLOBAL_O3
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

      USE CMN_SIZE_MOD  ! Size parameters
! 
! !REVISION HISTORY: 
!  13 Jul 2004 - R. Yantosca - Initial version
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 7/13/04)
!  (2 ) Now dimension O3 with LLTROP (bmy, 12/1/05)
!  (3 ) Now dimension O3 with LLPAR (phs, 1/19/07)
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local variables
      INTEGER             :: AS

      !=================================================================
      ! INIT_GLOBAL_O3 begins here!
      !=================================================================
      ALLOCATE( O3( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'O3' )
      O3 = 0d0

      END SUBROUTINE INIT_GLOBAL_O3
!EOC      
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_global_o3
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_O3 deallocates the O3 array.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GLOBAL_O3
!
! !REVISION HISTORY: 
!  13 Jul 2004 - R. Yantosca - Initial version
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_GLOBAL_O3 begins here!
      !=================================================================
      IF ( ALLOCATED( O3 ) ) DEALLOCATE( O3 ) 
     
      END SUBROUTINE CLEANUP_GLOBAL_O3
!EOC
      END MODULE GLOBAL_O3_MOD
