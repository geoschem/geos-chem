!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_hno3_mod
!
! !DESCRIPTION: Module GLOBAL\_HNO3\_MOD contains variables and routines for
! reading the global monthly mean HNO3 fields from disk. (bmy, 10/15/02, 2/7/07)
!\\
!\\
! !INTERFACE: 
!
      MODULE GLOBAL_HNO3_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS: 
!
      PUBLIC :: GET_HNO3_UGM3
      PUBLIC :: GET_HNO3_VV
      PUBLIC :: GET_GLOBAL_HNO3
      PUBLIC :: CLEANUP_GLOBAL_HNO3
! 
! !REVISION HISTORY: 
!  (1 ) Minor bug fix in FORMAT statement (bmy, 3/23/03)
!  (2 ) Cosmetic changes (bmy, 3/27/03)
!  (3 ) Now references "directory_mod.f" & "tracer_mod.f" (bmy, 7/20/04)
!  (4 ) Now suppress output from READ_BPCH2 with QUIET=T (bmy, 1/14/05)
!  (5 ) Now read total gas + aerosol HNO3 data (bec, bmy, 4/13/05)
!  (6 ) Now read files from "sulfate_sim_200508/offline" dir (bmy, 8/1/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!   8 Feb 2012 - R. Yantosca - Add modifications for GEOS_5.7.x
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  05 Sep 2013 - M. Sulprizio- Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
!
! !PRIVATE TYPES:
!
      ! Array to store global monthly mean HNO3 field
      REAL*8, ALLOCATABLE :: HNO3(:,:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_hno3_ugm3
!
! !DESCRIPTION: Function GET\_HNO3\_UGM3 converts monthly mean HNO3 mixing
! ratio from [v/v] to [ug/m3].  This is necessary for the RPMARES code.  We
! allow HNO3 concentrations to evolve but relax back to the monthly mean value
! every 3 hours. (bmy, 10/15/02, 7/20/04)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_HNO3_UGM3( I, J, L, State_Met ) RESULT( HNO3_UGM3 )
!
! !USES:
!
      USE CMN_SIZE_MOD             ! Size parameters
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TRACER_MOD,         ONLY : TCVV
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)  :: I, J, L     ! Lon, lat, alt indices
      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !RETURN VALUE:
!
      REAL*8                      :: HNO3_UGM3
!
! !REVISION HISTORY:
!  (1 ) Now references TCVV from "tracer_mod.f" (bmy, 7/20/04)
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  05 Sep 2013 - M. Sulprizio- Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! GET_HNO3_UGM3 begins here!
      !=================================================================

      ! First convert HNO3 from [v/v] to [kg]
      HNO3_UGM3 = HNO3(I,J,L) * State_Met%AD(I,J,L) / ( 28.97d0 / 63d0 )

      ! Then convert HNO3 from [kg] to [ug/m3]
      HNO3_UGM3 = HNO3_UGM3 * 1.d9 / State_Met%AIRVOL(I,J,L)

      END FUNCTION GET_HNO3_UGM3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_hno3_vv
!
! !DESCRIPTION: Function GET\_HNO3\_VV returns HNO3 concentrations in units
! of v/v (bec, 12/23/11)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_HNO3_VV( I, J, L ) RESULT( HNO3_VV )
!
! !USES:
!
      USE CMN_SIZE_MOD      ! Size parameters
      USE TRACER_MOD,  ONLY : TCVV
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I, J, L     ! Lon, lat, alt indices
!
! !RETURN VALUE:
!
      REAL*8              :: HNO3_VV
!
! !REVISION HISTORY:
!  (1 ) Now references TCVV from "tracer_mod.f" (bmy, 7/20/04)
!  05 Sep 2013 - M. Sulprizio- Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
      !=================================================================
      ! GET_HNO3_vv begins here!
      !=================================================================

      ! Convert HNO3 from [ppbv] to [v/v] 
      HNO3_VV = HNO3(I,J,L) * 1.d-9 !HNO3 is in ppbv, not v/v (bec, 9/28/09)

      END FUNCTION GET_HNO3_VV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_global_hno3
!
! !DESCRIPTION: Subroutine GET\_GLOBAL\_HNO3 reads global OH from binary punch
!  files stored in the data directory.  This is needed for the offline sulfate
!  simulation. (bmy, 10/3/02, 2/7/07)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_GLOBAL_HNO3( THISMONTH )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,     READ_BPCH2
      USE ERROR_MOD,     ONLY : ERROR_STOP
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D,  TRANSFER_3D
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: THISMONTH   ! Current month number (1-12)
!
! !REVISION HISTORY:
!  (1 ) Bug fix in FORMAT statement: Replace missing commas (bmy, 3/23/03)
!  (2 ) Cosmetic changes (bmy, 3/27/03)
!  (3 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (4 ) Now suppress output from READ_BPCH2 with QUIET=T (bmy, 1/14/05)
!  (5 ) Now read total gas + aerosol HNO3 data (bec, bmy, 4/13/05)
!  (6 ) GEOS-3 and GEOS-4 data comes from model runs w/ 30 layers.  Also now
!        read from "sulfate_sim_200508/offline" directory (bmy, 8/1/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (8 ) Renamed GRID30LEV to GRIDREDUCED (bmy, 2/7/07)
!  (9 ) Treat MERRA in the same way as for GEOS-5 (bmy, 8/13/10)
!  08 Feb 2012 - R. Yantosca - Treat GEOS-5.7.x in the same way as MERRA
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  05 Sep 2013 - M. Sulprizio- Added ProTeX headers
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J, L
      REAL*4              :: ARRAY(IIPAR,JJPAR,LGLOB)
      REAL*8              :: XTAU
      CHARACTER(LEN=255)  :: FILENAME

      ! First time flag
      LOGICAL, SAVE       :: FIRST = .TRUE. 

      !=================================================================
      ! GET_GLOBAL_HNO3 begins here!
      !=================================================================

      ! Allocate OH array, if this is the first call
      IF ( FIRST ) THEN
         CALL INIT_GLOBAL_HNO3
         FIRST = .FALSE.
      ENDIF

      ! File name for modified HNO3 (total gas + aerosol nitrate)
      ! after sea-salt chemistry (bec, bmy, 4/13/05, 8/1/05)
      FILENAME = TRIM( DATA_DIR )                      // 
     &           'sulfate_sim_200508/offline/THNO3.'   // 
     &           GET_NAME_EXT() // '.' // GET_RES_EXT()!

      ! Echo some information to the standard output
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - GET_GLOBAL_HNO3: Reading ', a )

      ! Get the TAU0 value for the start of the given month
      ! Assume "generic" year 1985 (TAU0 = [0, 744, ... 8016])
      XTAU = GET_TAU0( THISMONTH, 1, 1985 )

#if   defined( GEOS_4 ) || defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )

      !-------------------------------------------------------
      ! GEOS-3 / GEOS-4 data come from the 30 level model run
      !-------------------------------------------------------

#if   defined( GRIDREDUCED )

      ! Read HNO3 data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$',          7,     
     &                 XTAU,      IIPAR,              JJPAR,      
     &                 LLPAR,     ARRAY(:,:,1:LLPAR), QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable HNO3
      DO L = 1, LLPAR
         CALL TRANSFER_2D( ARRAY(:,:,L), HNO3(:,:,L) )
      ENDDO

#else
 
      ! If LLPAR is not 30 levels then stop with error
      CALL ERROR_STOP( 'Must use 30 levels for offline aerosol sim!',
     &                 'GET_GLOBAL_HNO3 ("global_hno3_mod.f!")' )

#endif

#else

      !-------------------------------------------------------
      ! Data for other GEOS grids have LGLOB levels 
      !-------------------------------------------------------

      ! Read HNO3 data from the binary punch file
      CALL READ_BPCH2( FILENAME, 'IJ-AVG-$', 7,     
     &                 XTAU,      IIPAR,     JJPAR,      
     &                 LGLOB,     ARRAY,     QUIET=.TRUE. )

      ! Assign data from ARRAY to the module variable HNO3
      CALL TRANSFER_3D( ARRAY, HNO3 )

#endif

      END SUBROUTINE GET_GLOBAL_HNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_global_hno3
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_HNO3 allocates and zeroes the HNO3
!  array (bmy, 10/15/02)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GLOBAL_HNO3
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,    ONLY : ALLOC_ERR
!
! !REVISION HISTORY:
!  (1 ) Now references ALLOC_ERR from "error_mod.f" (bmy, 10/15/02)
!  (2 ) Now dimension HNO3 as (IIPAR,JJPAR,LLPAR) (bmy, 8/1/05)
!  05 Sep 2013 - M. Sulprizio- Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: AS

      !=================================================================
      ! INIT_GLOBAL_HNO3 begins here!
      !=================================================================
      ALLOCATE( HNO3( IIPAR, JJPAR, LLPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'HNO3' )
      HNO3 = 0d0

      END SUBROUTINE INIT_GLOBAL_HNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_global_hno3
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_HNO3 deallocates the HNO3 array.
! (bmy, 10/15/02)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_GLOBAL_HNO3
!
! !REVISION HISTORY:
!  05 Sep 2013 - M. Sulprizio- Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC                            
      !=================================================================
      ! CLEANUP_GLOBAL_HNO3 begins here!
      !=================================================================
      IF ( ALLOCATED( HNO3 ) ) DEALLOCATE( HNO3 ) 
     
      END SUBROUTINE CLEANUP_GLOBAL_HNO3
!EOC
      END MODULE GLOBAL_HNO3_MOD
