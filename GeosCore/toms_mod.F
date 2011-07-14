!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: toms_mod
!
! !DESCRIPTION: Module TOMS\_MOD contains variables and routines for reading 
!  the TOMS/SBUV O3 column data from disk (for use w/ the FAST-J photolysis 
!  routines).
!\\
!\\
! !INTERFACE: 
!
      MODULE TOMS_MOD
!
! !USES:
!
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC DATA MEMBERS:
!
      REAL*8, PUBLIC, ALLOCATABLE :: TOMS(:,:)               
      REAL*8, PUBLIC, ALLOCATABLE :: DTOMS1(:,:)             
      REAL*8, PUBLIC, ALLOCATABLE :: DTOMS2(:,:)             
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC                      :: CLEANUP_TOMS
      PUBLIC                      :: READ_TOMS
!
! !PRIVATE MEMBER FUNCTIONS:
! 
      PRIVATE                     :: INIT_TOMS
!
! !REMARKS:
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
! !REVISION HISTORY:
!  14 Jul 2003 - R. Yantosca - Initial version
!  (1 ) Now references "directory_mod.f" (bmy, 7/20/04)
!  (2 ) Now can read files for GEOS or GCAP grids (bmy, 8/16/05)
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (4 ) Now always use 2002 TOMS O3 data for GCAP (swu, bmy, 10/3/06)
!  (5 ) Now reads from TOMS_200701 directory, w/ updated data (bmy, 2/1/07)
!  (6 ) Now don't replace any tokens in the DATA_DIR variable (bmy, 12/5/07)
!  (7 ) Latest year of TOMS data is now 2007 (bmy, 1/14/09)
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_toms
!
! !DESCRIPTION: Subroutine READ\_TOMS reads in TOMS O3 column data from a 
!  binary punch file for the given grid, month and year. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READ_TOMS( THISMONTH, THISYEAR )
!
! !USES:
!
      USE BPCH2_MOD,     ONLY : GET_NAME_EXT_2D
      USE BPCH2_MOD,     ONLY : GET_RES_EXT
      USE BPCH2_MOD,     ONLY : GET_TAU0
      USE BPCH2_MOD,     ONLY : READ_BPCH2
      USE DIRECTORY_MOD, ONLY : DATA_DIR
      USE TIME_MOD,      ONLY : EXPAND_DATE
      USE TRANSFER_MOD,  ONLY : TRANSFER_2D

#     include "CMN_SIZE"                    ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: THISMONTH   ! Current month
      INTEGER, INTENT(IN)    :: THISYEAR    ! Current year
!
! !REMARKS:
!  TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 3.
!  Resolution:  5 x 10 deg.
!                                                                             .
!  Methodology (bmy, 2/12/07)
!  ----------------------------------------------------------------
!  FAST-J comes with its own default O3 column climatology (from 
!  McPeters 1992 & Nagatani 1991), which is stored in the input 
!  file "jv_atms.dat".  These "FAST-J default" O3 columns are used 
!  in the computation of the actinic flux and other optical 
!  quantities for the FAST-J photolysis.  
!                                                                             .
!  The TOMS/SBUV O3 columns and 1/2-monthly O3 trends (contained 
!  in the TOMS_200701 directory) are read into GEOS-Chem by routine 
!  READ_TOMS in "toms_mod.f".  Missing values (i.e. locations where 
!  there are no data) in the TOMS/SBUV O3 columns are defined by 
!  the flag -999.  
!                                                                             .
!  After being read from disk in routine READ_TOMS, the TOMS/SBUV 
!  O3 data are then passed to the FAST-J routine "set_prof.f".  In 
!  "set_prof.f", a test is done to make sure that the TOMS/SBUV O3 
!  columns and 1/2-monthly trends do not have any missing values 
!  for (lat,lon) location for the given month.  If so, then the 
!  TOMS/SBUV O3 column data is interpolated to the current day and 
!  is used to weight the "FAST-J default" O3 column.  This 
!  essentially "forces" the "FAST-J default" O3 column values to 
!  better match the observations, as defined by TOMS/SBUV.
!                                                                             .
!  If there are no TOMS/SBUV O3 columns (and 1/2-monthly trends) 
!  at a (lat,lon) location for given month, then FAST-J will revert 
!  to its own "default" climatology for that location and month.  
!  Therefore, the TOMS O3 can be thought of as an  "overlay" data 
!  -- it is only used if it exists.
!                                                                             .
!  Note that there are no TOMS/SBUV O3 columns at the higher 
!  latitudes.  At these latitudes, the code will revert to using 
!  the "FAST-J default" O3 columns.
!                                                                             .
!  As of February 2007, we have TOMS/SBUV data for 1979 thru 2005.  
!  2006 TOMS/SBUV data is incomplete as of this writing.  For years
!  2006 and onward, we use 2005 TOMS O3 columns.
!                                                                             .
!  This methodology was originally adopted by Mat Evans.  Symeon 
!  Koumoutsaris was responsible for creating the downloading and 
!  processing the TOMS O3 data files from 1979 thru 2005 in the 
!  TOMS_200701 directory.
!  
! !REVISION HISTORY: 
!  10 Dec 2002 - M. Evans - Initial version
!  (1 ) Bundled into "toms_mod.f" (bmy, 7/14/03)
!  (2 ) Now references DATA_DIR from "directory_mod.f" (bmy, 7/20/04)
!  (3 ) Now can read files for GEOS or GCAP grids (bmy, 8/16/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Now always use 2002 TOMS O3 data for GCAP (swu, bmy, 10/3/06)
!  (6 ) Now reads from TOMS_200701 directory, w/ updated data.  Also always
!        use 1979 data prior to 1979 or 2005 data after 2005. (bmy, 2/12/07)
!  (7 ) Bug fix: don't include DATA_DIR in filename, just in case someone's 
!        file path has replaceable tokens (e.g. hh, mm, MM etc.) (bmy, 12/5/07)
!  (8 ) Latest year of TOMS data is now 2007 (bmy, 1/14/09)
!  (9 ) Updated TOMS data in TOMS_200906. Latest year is 2008. (ccc, 6/15/09)
!  08 Dec 2009 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL            :: FIRST = .TRUE.
      INTEGER            :: YYYYMMDD, YEAR
      REAL*4             :: ARRAY(IIPAR,JJPAR,1)
      REAL*8             :: XTAU
      CHARACTER(LEN=255) :: FILENAME

      !=================================================================
      ! Initialization
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

      ! Use 2008 data after 2008
      IF ( YEAR > 2008 ) THEN
         WRITE( 6, 105 ) YEAR
         YEAR = 2008
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
     
      ! Define filename (with replaceable tokens)
      FILENAME = 'TOMS_200906/TOMS_O3col_YYYY.' // GET_NAME_EXT_2D() //
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
     &                 XTAU,      IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )         

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), TOMS )

      !--------------------------------
      ! d(TOMS)/dT (1st half of month)
      !--------------------------------

       ! Read data
      CALL READ_BPCH2( FILENAME, 'TOMS-O3',  2, 
     &                 XTAU,      IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )         

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), DTOMS1 )
      
      !--------------------------------
      ! d(TOMS)/dT (2nd half of month)
      !--------------------------------

       ! Read data: 
      CALL READ_BPCH2( FILENAME, 'TOMS-O3',  3, 
     &                 XTAU,      IIPAR,     JJPAR,     
     &                 1,         ARRAY,     QUIET=.TRUE. )         

      ! Cast to REAL*8 and resize if necessary
      CALL TRANSFER_2D( ARRAY(:,:,1), DTOMS2 )

      END SUBROUTINE READ_TOMS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_toms
!
! !DESCRIPTION: Subroutine INIT\_TOMS allocates and zeroes all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_TOMS
!
! !USES:
!
      USE ERROR_MOD, ONLY : ALLOC_ERR

#     include "CMN_SIZE"  ! Size parameters
! 
! !REVISION HISTORY: 
!  14 Jul 2003 - R. Yantosca - Initial version
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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

      END SUBROUTINE INIT_TOMS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_toms
!
! !DESCRIPTION: Subroutine CLEANUP\_TOMS deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_TOMS
! 
! !REVISION HISTORY: 
!  14 Jul 2003 - R. Yantosca - Initial version
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_TOMS begins here!
      !=================================================================
      IF ( ALLOCATED( TOMS   ) ) DEALLOCATE( TOMS   )
      IF ( ALLOCATED( DTOMS1 ) ) DEALLOCATE( DTOMS1 )
      IF ( ALLOCATED( DTOMS2 ) ) DEALLOCATE( DTOMS2 )

      END SUBROUTINE CLEANUP_TOMS
!EOC
      END MODULE TOMS_MOD
