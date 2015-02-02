! $Id: uvalbedo_mod.f,v 1.1 2009/09/16 14:05:58 bmy Exp $
      MODULE UVALBEDO_MOD
!
!******************************************************************************
!  Module UVALBEDO_MOD contains variables and routines for reading the UV
!  Albedo data from disk (for use w/ the FAST-J photolysis routines).
!  (bmy, 4/19/02, 10/3/05)
!
!  Module Variables:
!  ============================================================================
!  (1 ) UVALBEDO (REAL(fp)) : Array to hold UV Albedo data from disk
!
!  Module Routines:
!  ============================================================================
!  (1 ) READ_UVALBEDO     : Routine to allocate UVALBEDO array and read data
!  (2 ) CLEANUP_UVALBEDO  : Routine to deallocate UVALBEDO array
!
!  GEOS-CHEM modules referenced by biomass_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module containing routines for binary punch file I/O
!  (2 ) directory_mod.f : Module containing GEOS-CHEM data & met field dirs
!  (3 ) error_mod.f     : Module containing NaN and other error check routines
!  (4 ) transfer_mod.f  : Module containing routines to cast & resize arrays
!
!  NOTES:
!  (1 ) Now read uvalbedo file directly from DATA_DIR/uvalbedo_200111
!        subdirectory.  (bmy, 4/2/02)
!  (2 ) Now divide module header into MODULE PRIVATE, MODULE VARIABLES, and
!        MODULE ROUTINES sections. (bmy, 5/28/02)
!  (3 ) Now references "error_mod.f" (bmy, 10/15/02)
!  (4 ) Minor modification in READ_UVALBEDO (bmy, 3/14/03)
!  (5 ) Now references "directory_mod.f" (bmy, 7/20/04)
!  (6 ) Bug fix for GCAP grid in READ_UVALBEDO (bmy, 8/16/05)
!  (7 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  24 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  12 Jan 2015 - R. Yantosca - Remove CLEANUP_UVALBEDO routine
!******************************************************************************
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement 
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE READ_UVALBEDO( am_I_Root, Input_Opt, 
     &                          MONTH,     State_Met, RC )
!
!******************************************************************************
!  Subroutine READ_UVALBEDO reads in UV albedo data from a binary punch
!  file for the given grid, model, and month. (bmy, 2/2/00, 10/3/05)  
!
!  Arguments as Input:
!  ==========================================================================
!  (1 ) MONTH    (INTEGER) : Current month (1-12)
!  (2 ) UVALBEDO (REAL(fp) ) : Array with UV albedo data 
!
!  Reference:
!  ==========================================================================
!  Herman, J.R and Celarier, E.A., "Earth surface reflectivity climatology
!     at 340-380 nm from TOMS data", JGR, Vol. 102, D23, pp. 28003-28011, 
!     Dec 20, 1997.
!
!  NOTES:
!  (1 ) Call READ_BPCH2 to read in the UV albedo data from the binary punch 
!        file. (bmy, 2/2/00)
!  (2 ) Cosmetic changes (bmy, 3/17/00)
!  (3 ) Reference F90 module "bpch2_mod" which contains routine "read_bpch2"
!        for reading data from binary punch files (bmy, 6/28/00)
!  (4 ) Remove IOS variable -- it wasn't used (bmy, 9/13/00)
!  (5 ) Now use GET_TAU0 to return the TAU0 values for 1985.  Also use 
!        TRANSFER_2D from "transfer_mod.f" to copy data from an array of 
!        size (IIPAR,JJPAR) to an array of size (IIPAR,JJPAR).  ARRAY needs 
!        to be of size (IIPAR,JJPAR).  Also updated comments and made 
!        cosmetic changes. (bmy, 9/26/01)
!  (6 ) Removed obsolete code from 9/01 (bmy, 10/24/01)
!  (7 ) Now echo FILENAME to the std output (bmy, 11/15/01)
!  (8 ) Bundled into "uvalbedo_mod.f" (bmy, 1/15/02)
!  (9 ) Now read uvalbedo file directly from DATA_DIR/uvalbedo_200111
!        subdirectory.  (bmy, 4/2/02)
!  (10) Now references ALLOC_ERR from "error_mod.f".  Also eliminated obsolete
!        code from 4/02.  Updated comments, cosmetic changes. (bmy, 10/15/02)
!  (11) Now call READ_BPCH2 with QUIET=.TRUE. to suppress printing of extra 
!        info to stdout.  Also made cosmetic changes. (bmy, 3/14/03)
!  (12) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (13) Read proper filename for GCAP or GEOS grids (swu, bmy, 8/15/05) 
!  (14) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,          ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,          ONLY : GET_TAU0,        READ_BPCH2
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      INTEGER,        INTENT(IN)    :: MONTH       ! Current month
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?

!
! !LOCAL VARIABLES:
!      
      LOGICAL              :: FIRST = .TRUE.
      INTEGER              :: AS
      REAL(f4)             :: ARRAY(IIPAR,JJPAR,1)
      REAL(f8)             :: XTAU
      CHARACTER(LEN=255)   :: FILENAME

      !=================================================================
      ! READ_UVALBEDO begins here!
      !=================================================================

      ! Assume success
      RC  =  GIGC_SUCCESS

      !=================================================================
      ! Read UVALBEDO data from disk
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Create filename
      FILENAME = TRIM( Input_Opt%DATA_DIR )  // 
     &           'uvalbedo_200111/uvalbedo.' // GET_NAME_EXT_2D() //
     &           '.'                         // GET_RES_EXT()

      ! Echo filename
      IF ( am_I_Root ) THEN
         WRITE( 6, 110 ) TRIM( FILENAME )
 110     FORMAT( '     - READ_UVALBEDO: Reading ', a )
      ENDIF

      ! Get TAU0 value for first day of the MONTH -- use generic year 1985
      XTAU = GET_TAU0( MONTH, 1, 1985 )

      ! Read data: UV albedos are tracer #1, category name "UVALBEDO"
      CALL READ_BPCH2( FILENAME, 'UVALBEDO', 1, 
     &                 XTAU,      IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )         

      ! Transfer data from REAL*4 to REAL(fp) and to size (IIPAR,JJPAR)
      State_Met%UVALBEDO = ARRAY(:,:,1)

      ! Return to calling program
      END SUBROUTINE READ_UVALBEDO

!------------------------------------------------------------------------------

      END MODULE UVALBEDO_MOD
