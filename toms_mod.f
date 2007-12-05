! $Id: toms_mod.f,v 1.8 2007/12/05 16:10:11 bmy Exp $
      MODULE TOMS_MOD
!
!******************************************************************************
!  Module TOMS_MOD contains variables and routines for reading the TOMS/SBUV
!  O3 column data from disk (for use w/ the FAST-J photolysis routines).
!  (mje, bmy, 7/14/03, 12/5/07)
!
!  Module Variables:
!  ============================================================================
!  (1 ) TOMS   (REAL*8) : Monthly mean TOMS O3 data [DU]
!  (2 ) DTOMS1 (REAL*8) : Time deriv. of TOMS O3 for 1st half of month [DU/day]
!  (3 ) DTOMS2 (REAL*8) : Time deriv. of TOMS O3 for 2nd half of month [DU/day]
!
!  Module Procedures:
!  ============================================================================
!  (1 ) READ_TOMS       : Routine to read TOMS data from disk
!  (2 ) INIT_TOMS       : Routine to allocate and zero module arrays
!  (3 ) CLEANUP_TOMS    : Routine to deallocate module arrays
!
!  GEOS-CHEM modules referenced by toms_mod.f
!  ============================================================================
!  (1 ) bpch2_mod.f     : Module containing routines for binary punch file I/O
!  (2 ) directory_mod.f : Module containing GEOS-CHEM data & met field dirs
!  (3 ) error_mod.f     : Module containing I/O error and NaN check routines
!  (4 ) time_mod.f      : Module containing routines to compute date & time
!  (5 ) transfer_mod.f  : Module containing routines to cast & resize arrays
!
!  References:
!  ============================================================================
!  TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 3.
!  Resolution:  5 x 10 deg.
!
!  Source: http://code916.gsfc.nasa.gov/Data_services/merged/index.html
!
!  Contact person for the merged data product:
!  Stacey Hollandsworth Frith (smh@hyperion.gsfc.nasa.gov)
!
!  NOTES:
!  (1 ) Now references "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now can read files for GEOS or GCAP grids (bmy, 8/16/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now always use 2002 TOMS O3 data for GCAP (swu, bmy, 10/3/06)
!  (5 ) Now reads from TOMS_200701 directory, w/ updated data (bmy, 2/1/07)
!  (6 ) Now don't replace any tokens in the DATA_DIR variable (bmy, 12/5/07)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "toms_mod.f"
      !=================================================================
      
      ! Make everything PRIVATE ...
      PRIVATE

      ! ... except these variables ...
      PUBLIC :: TOMS
      PUBLIC :: DTOMS1
      PUBLIC :: DTOMS2

      ! ... and these routines
      PUBLIC :: CLEANUP_TOMS
      PUBLIC :: READ_TOMS

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Arrays
      REAL*8, ALLOCATABLE :: TOMS(:,:)               
      REAL*8, ALLOCATABLE :: DTOMS1(:,:)             
      REAL*8, ALLOCATABLE :: DTOMS2(:,:)             
       
      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------

      SUBROUTINE READ_TOMS( THISMONTH, THISYEAR )
!
!******************************************************************************
!  Subroutine READ_TOMS reads in TOMS O3 column data from a binary punch
!  file for the given grid, month and year. (mje, bmy 12/10/02, 2/12/07)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) THISMONTH (INTEGER) : Current month (1-12)
!  (2 ) THISYEAR  (INTEGER) : Current year (YYYY format)
!
!  NOTES:
!  (1 ) Bundled into "toms_mod.f" (bmy, 7/14/03)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now can read files for GEOS or GCAP grids (bmy, 8/16/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Now always use 2002 TOMS O3 data for GCAP (swu, bmy, 10/3/06)
!  (6 ) Now reads from TOMS_200701 directory, w/ updated data.  Also always
!        use 1979 data prior to 1979 or 2005 data after 2005. (bmy, 2/12/07)
!  (7 ) Bug fix: don't include DATA_DIR in filename, just in case someone's 
!        file path has replaceable tokens (e.g. hh, mm, MM etc.) (bmy, 12/5/07)
!******************************************************************************
!
      ! References to F90 modules
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D, GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0,        READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TIME_MOD,      ONLY : EXPAND_DATE
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"      ! Size parameters

      ! Arguments
      INTEGER, INTENT(IN)    :: THISMONTH
      INTEGER, INTENT(IN)    :: THISYEAR

      ! Local Variables
      LOGICAL                :: FIRST = .TRUE.
      INTEGER                :: YYYYMMDD, YEAR
      REAL*4                 :: ARRAY(IGLOB,JGLOB,1)
      REAL*8                 :: XTAU
      CHARACTER(LEN=255)     :: FILENAME

      !=================================================================
      ! READ_TOMS begins here!
      !
      ! TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 3.
      ! Resolution:  5 x 10 deg.
      !
      ! Methodology (bmy, 2/12/07)
      ! ----------------------------------------------------------------
      ! FAST-J comes with its own default O3 column climatology (from 
      ! McPeters 1992 & Nagatani 1991), which is stored in the input 
      ! file "jv_atms.dat".  These "FAST-J default" O3 columns are used 
      ! in the computation of the actinic flux and other optical 
      ! quantities for the FAST-J photolysis.  
      !
      ! The TOMS/SBUV O3 columns and 1/2-monthly O3 trends (contained 
      ! in the TOMS_200701 directory) are read into GEOS-Chem by routine 
      ! READ_TOMS in "toms_mod.f".  Missing values (i.e. locations where 
      ! there are no data) in the TOMS/SBUV O3 columns are defined by 
      ! the flag -999.  
      ! 
      ! After being read from disk in routine READ_TOMS, the TOMS/SBUV 
      ! O3 data are then passed to the FAST-J routine "set_prof.f".  In 
      ! "set_prof.f", a test is done to make sure that the TOMS/SBUV O3 
      ! columns and 1/2-monthly trends do not have any missing values 
      ! for (lat,lon) location for the given month.  If so, then the 
      ! TOMS/SBUV O3 column data is interpolated to the current day and 
      ! is used to weight the "FAST-J default" O3 column.  This 
      ! essentially "forces" the "FAST-J default" O3 column values to 
      ! better match the observations, as defined by TOMS/SBUV.
      !
      ! If there are no TOMS/SBUV O3 columns (and 1/2-monthly trends) 
      ! at a (lat,lon) location for given month, then FAST-J will revert 
      ! to its own "default" climatology for that location and month.  
      ! Therefore, the TOMS O3 can be thought of as an  "overlay" data 
      ! -- it is only used if it exists.
      !
      ! Note that there are no TOMS/SBUV O3 columns at the higher 
      ! latitudes.  At these latitudes, the code will revert to using 
      ! the "FAST-J default" O3 columns.
      !
      ! As of February 2007, we have TOMS/SBUV data for 1979 thru 2005.  
      ! 2006 TOMS/SBUV data is incomplete as of this writing.  For years
      ! 2006 and onward, we use 2005 TOMS O3 columns.
      !
      ! This methodology was originally adopted by Mat Evans.  Symeon 
      ! Koumoutsaris was responsible for creating the downloading and 
      ! processing the TOMS O3 data files from 1979 thru 2005 in the 
      ! TOMS_200701 directory.
      !=================================================================

      ! Allocate arrays on the first call only
      IF ( FIRST ) THEN
         CALL INIT_TOMS
         FIRST = .FALSE.
      ENDIF

      ! Always use 2002 data for GCAP
#if   defined ( GCAP )
      YEAR = 2002
#else 
      YEAR = THISYEAR
#endif

      ! Use 1979 data prior to 1979
      IF ( YEAR < 1979 ) THEN
         WRITE( 6, 100 ) YEAR
         YEAR = 1979
      ENDIF

      ! Use 2005 data after 2005
      IF ( YEAR > 2005 ) THEN
         WRITE( 6, 105 ) YEAR
         YEAR = 2005
      ENDIF
      
      ! FORMAT statemetns
 100  FORMAT( '     - READ_TOMS: No data for ',i4,', using 1979!' )
 105  FORMAT( '     - READ_TOMS: No data for ',i4,', using 2005!' )

      !=================================================================
      ! Read TOMS data from disk
      !=================================================================

      ! Get TAU0 value for first day of the MONTH 
      XTAU     = GET_TAU0( THISMONTH, 1, YEAR )

      ! Create YYYYMMDD value
      YYYYMMDD = ( YEAR * 10000 ) + ( THISMONTH * 100 ) + 01
     
!-----------------------------------------------------------------------------
! Prior to 12/5/07:
! Bug fix: don't include DATA_DIR in filename, just in case someone's file
! path contains replaceable tokens (e.g. hh, mm, MM etc.) (bmy, 12/5/07)
!      ! Define filename
!      FILENAME = TRIM( DATA_DIR )               // 
!     &           'TOMS_200701/TOMS_O3col_YYYY.' // GET_NAME_EXT_2D() //
!     &           '.'                            // GET_RES_EXT()
!-----------------------------------------------------------------------------

      ! Define filename (with replaceable tokens)
      FILENAME = 'TOMS_200701/TOMS_O3col_YYYY.' // GET_NAME_EXT_2D() //
     &           '.'                            // GET_RES_EXT()

      ! Replace YYYY token with current year
      CALL EXPAND_DATE( FILENAME, YYYYMMDD, 000000 )

      ! Now prefix the data directory
      FILENAME = TRIM( DATA_DIR ) // TRIM( FILENAME )

      ! Echo filename
      WRITE( 6, 110 ) TRIM( FILENAME )
 110  FORMAT( '     - READ_TOMS: Reading ', a )

      !-----------------------------
      ! TOMS O3 columns
      !-----------------------------
      
      ! Read data
      CALL READ_BPCH2( FILENAME, 'TOMS-O3',  1, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. )         

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), TOMS )

      !--------------------------------
      ! d(TOMS)/dT (1st half of month)
      !--------------------------------

       ! Read data
      CALL READ_BPCH2( FILENAME, 'TOMS-O3',  2, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. )         

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), DTOMS1 )
      
      !--------------------------------
      ! d(TOMS)/dT (2nd half of month)
      !--------------------------------

       ! Read data: 
      CALL READ_BPCH2( FILENAME, 'TOMS-O3',  3, 
     &                 XTAU,      IGLOB,     JGLOB,     
     &                 1,         ARRAY,     QUIET=.TRUE. )         

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), DTOMS2 )

      ! Return to calling program
      END SUBROUTINE READ_TOMS

!------------------------------------------------------------------------------

      SUBROUTINE INIT_TOMS
!
!******************************************************************************
!  Subroutine INIT_TOMS allocates and zeroes all module arrays (bmy, 7/14/03)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"  ! Size parameters

      ! Local variables
      INTEGER :: AS

      !=================================================================
      ! INIT_TOMS begins here!
      !=================================================================

      ! Allocate TOMS
      ALLOCATE( TOMS( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TOMS' )
      TOMS = 0d0

      ! Allocate DTOMS
      ALLOCATE( DTOMS1( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DTOMS1' )
      DTOMS1 = 0d0

      ! Allocate DTOMS2
      ALLOCATE( DTOMS2( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'DTOMS2' )
      DTOMS2 = 0d0

      ! Return to calling program
      END SUBROUTINE INIT_TOMS

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_TOMS
!
!******************************************************************************
!  Subroutine CLEANUP_TOMS deallocates all module arrays (bmy, 7/14/03)
!
!  NOTES:
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_TOMS begins here!
      !=================================================================
      IF ( ALLOCATED( TOMS   ) ) DEALLOCATE( TOMS   )
      IF ( ALLOCATED( DTOMS1 ) ) DEALLOCATE( DTOMS1 )
      IF ( ALLOCATED( DTOMS2 ) ) DEALLOCATE( DTOMS2 )

      ! Return to calling program
      END SUBROUTINE CLEANUP_TOMS

!------------------------------------------------------------------------------

      ! End of module
      END MODULE TOMS_MOD
