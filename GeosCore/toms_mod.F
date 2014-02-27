!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      USE CMN_SIZE_MOD                    ! Size parameters

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
! 
      PUBLIC :: INIT_TOMS
      PUBLIC :: READ_TOMS
      PUBLIC :: COMPUTE_OVERHEAD_O3
      PUBLIC :: GET_OVERHEAD_O3
      PUBLIC :: CLEANUP_TOMS
!
! !PUBLIC DATA MEMBERS:
!
      ! First & last years for which TOMS/SBUV data is is available
      ! (update these as new data is added to the archive)
      INTEGER, PUBLIC, PARAMETER :: FIRST_TOMS_YEAR = 1979
#if defined( GRID4x5 ) || defined( GRID2x25 ) || defined( GRID025x03125 )
      INTEGER, PUBLIC, PARAMETER :: LAST_TOMS_YEAR  = 2010
#else
      INTEGER, PUBLIC, PARAMETER :: LAST_TOMS_YEAR  = 2008
#endif
!
! !REMARKS:
!  References:
!  ============================================================================
!  Version 8 Merged Ozone Data Sets
!  Total Ozone Revision 05 
!  DATA THROUGH: MAR 2009
!  LAST MODIFIED: 01 MAY 2009
!                                                                             .
!  http://acdb-ext.gsfc.nasa.gov/Data_services/merged/index.html
!                                                                             .
!  TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 5.
!  Resolution:  5 x 10 deg.
!                                                                             .
!  * Includes reprocessed N16 and N17 SBUV/2 data using latest calibration. 
!  * OMI data updated from Collection 2 to Collection 3. 
!  * New offsets derived based on revised data sets. 
!  * 1970-1972 N4 BUV data added with no adjustments. User may wish to apply 
!    offset based on Comparisons between BUV and Dobson Measurements. 
!                                                                             .
!  Responsible NASA official:
!  Dr. Richard Stolarski (Richard.S.Stolarski@nasa.gov)
!  Stacey Frith          (Stacey.M.Frith@nasa.gov     )
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
!  06 Mar 2012 - R. Yantosca - Added function GET_TOTAL_O3
!  06 Mar 2012 - R. Yantosca - Added parameters FIRST_TOMS_YEAR, LAST_TOMS_YEAR
!  06 Mar 2012 - R. Yantosca - Updated comments
!  06 Mar 2012 - R. Yantosca - Now make TOMS, DTOMS1, DTOMS2 arrays PRIVATE
!  06 Mar 2012 - R. Yantosca - Add TO3_DAILY as a PRIVATE module array
!  25 Jun 2012 - S. Kim      - Now reads from TOMS_201203 directory, w/ 
!                              updated data
!  03 Jul 2012 - R. Yantosca - Restrict reading from TOMS_201203 directory
!                              to GEOS-5.7.2 met data for the time being.
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      REAL*8, PRIVATE, ALLOCATABLE :: TOMS(:,:)               
      REAL*8, PRIVATE, ALLOCATABLE :: DTOMS1(:,:)             
      REAL*8, PRIVATE, ALLOCATABLE :: DTOMS2(:,:)     
#if   defined( DEVEL )
      REAL*8, PUBLIC,  ALLOCATABLE :: TO3_DAILY(:,:)
#else
      REAL*8, PRIVATE, ALLOCATABLE :: TO3_DAILY(:,:)
#endif

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      SUBROUTINE READ_TOMS( THISMONTH, THISYEAR, USE_O3_FROM_MET )
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
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: THISMONTH       ! Current month
      INTEGER, INTENT(IN)    :: THISYEAR        ! Current year
      LOGICAL, INTENT(IN)    :: USE_O3_FROM_MET ! Use TO3 from met fields?
!
! !REMARKS:
!  TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 5.
!  Resolution:  5 x 10 deg.
!                                                                             .
!  Methodology
!  ------------------------------------------------------------------------
!  FAST-J comes with its own default O3 column climatology (from McPeters 
!  1992 & Nagatani 1991), which is stored in the input file "jv_atms.dat".  
!  These "FAST-J default" O3 columns are used in the computation of the 
!  actinic flux and other optical quantities for the FAST-J photolysis.  
!                                                                             .
!  The TOMS/SBUV O3 columns and 1/2-monthly O3 trends (contained in the 
!  TOMS_200906 directory) are read into GEOS-Chem by routine READ_TOMS in 
!  "toms_mod.f".  Missing values (i.e. locations where there are no data) 
!  in the TOMS/SBUV O3 columns are defined by the flag -999.  
!                                                                             .
!  After being read from disk in routine READ_TOMS, the TOMS/SBUV O3 data 
!  are then passed to the FAST-J routine "set_prof.F".  In "set_prof.F", a 
!  test is done to make sure that the TOMS/SBUV O3 columns and 1/2-monthly 
!  trends do not have any missing values for (lat,lon) location for the given
!  month.  If so, then the TOMS/SBUV O3 column data is interpolated to the 
!  current day and is used to weight the "FAST-J default" O3 column.  This 
!  essentially "forces" the "FAST-J default" O3 column values to better match
!  the observations, as defined by TOMS/SBUV.
!                                                                             .
!  If there are no TOMS/SBUV O3 columns (and 1/2-monthly trends) at a (lat,
!  lon) location for given month, then FAST-J will revert to its own "default"
!  climatology for that location and month.  Therefore, the TOMS O3 can be 
!  thought of as an  "overlay" data -- it is only used if it exists.
!                                                                             .
!  Note that there are no TOMS/SBUV O3 columns at the higher latitudes.  
!  At these latitudes, the code will revert to using the "FAST-J default" 
!  O3 columns.
!                                                                             .
!  As of March 2012, we have TOMS/SBUV data for 1979 thru 2008.  We will
!  update to the latest TOMS/SBUV data set shortly.
!                                                                             .
!  This methodology was originally adopted by Mat Evans.
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
!  03 Jul 2012 - R. Yantosca - Restrict reading from TOMS_201203 directory
!                              to GEOS-5.7.2 met data for the time being.
!  13 Dec 2013 - M. Sulprizio- Now pass USE_O3_FROM_MET logical flag so that
!                              we bypass reading TOMS O3 data when using TO3
!                              from the met fields.
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

      ! Only read TOMS if not using TO3 from the met fields (mpayer, 12/13/13)
      IF ( USE_O3_FROM_MET ) RETURN

      ! Always use 2002 data for GCAP
#if   defined ( GCAP )
      YEAR = 2002
#else 
      YEAR = THISYEAR
#endif

      ! Use 1979 data prior to 1979
      IF ( YEAR < FIRST_TOMS_YEAR ) THEN
         WRITE( 6, 100 ) YEAR
         YEAR = FIRST_TOMS_YEAR
      ENDIF

      ! Use 2010 data after 2010
      IF ( YEAR > LAST_TOMS_YEAR ) THEN
         WRITE( 6, 105 ) YEAR
         YEAR = LAST_TOMS_YEAR
      ENDIF
      

      ! FORMAT statemetns
 100  FORMAT( '     - READ_TOMS: No data for ',i4,', using 1979!' )
 105  FORMAT( '     - READ_TOMS: No data for ',i4,', using 2010!' )

      !=================================================================
      ! Read TOMS data from disk
      !=================================================================

      ! Get TAU0 value for first day of the MONTH 
      XTAU     = GET_TAU0( THISMONTH, 1, YEAR )

      ! Create YYYYMMDD value
      YYYYMMDD = ( YEAR * 10000 ) + ( THISMONTH * 100 ) + 01
     
      ! Define filename (with replaceable tokens)
#if   !defined( GCAP )

#if defined( GRID4x5 ) || defined( GRID2x25 ) || defined( GRID025x03125 )
      FILENAME = 'TOMS_201203/TOMS_O3col_YYYY.' // GET_NAME_EXT_2D() //
     &           '.'                            // GET_RES_EXT()
#else
      FILENAME = 'TOMS_200906/TOMS_O3col_YYYY.' // GET_NAME_EXT_2D() //
     &           '.'                            // GET_RES_EXT()
#endif

#else
      FILENAME = 'TOMS_200701/TOMS_O3col_YYYY.' // GET_NAME_EXT_2D() //
     &           '.'                            // GET_RES_EXT()
#endif

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
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: compute_overhead_o3
!
! !DESCRIPTION: Subroutine COMPUTE\_OVERHEAD\_O3 returns the resultant total 
!  overhead O3 column for the FAST-J photolysis.  This will be one of two
!  options:

!  \begin{enumerate}
!  \item Default: TOMS/SBUV overhead O3 columns.  These will be used be
!        the FAST-J routine set\_prof.F to overwrite the existing FAST-J
!        climatology (cf McPeters \& Nagatani 1992).  Missing data (i.e.
!        for months \& locations where TOMS/SBUV data does not exist) 
!        is denoted by the value -999; FAST-J will skip over these points.
!  \item Overhead O3 columns taken directly from the met fields.  These 
!        will be returned if the flag  USE\_O3\_FROM\_MET is set to TRUE.
!  \end{enumerate}
!
! !INTERFACE:
!
      SUBROUTINE COMPUTE_OVERHEAD_O3( DAY, USE_O3_FROM_MET, TO3 )
!
! !INPUT PARAMETERS: 
!
      ! Day of month 
      INTEGER, INTENT(IN)  :: DAY      

      ! Switch to denote if we should use the default TO3
      ! directly from the met fields
      LOGICAL, INTENT(IN)  :: USE_O3_FROM_MET
      
      ! TO3 from the met fields [Dobsons]
      REAL*8,  INTENT(IN)  :: TO3(IIPAR,JJPAR)
!
! !REMARKS:
! Reference for the TOMS/SBUV merged O3 columns:
!                                                                             .
! 1985 - 2005 are taken from:
!                                                                             .
!   http://code916.gsfc.nasa.gov/Data_services/merged/index.html
!                                                                             .
!   TOMS/SBUV MERGED TOTAL OZONE DATA, Version 8, Revision 3.
!   Resolution:  5 x 10 deg.
!                                                                             .
!   Contact person for the merged data product:
!   Stacey Hollandsworth Frith (smh@hyperion.gsfc.nasa.gov)
!                                                                             .
! 2006 and 2007 are taken from:
!                                                                             .
!    http://code916.gsfc.nasa.gov/Data_services/merged/index.html
!                                                                             .
!    Version 8 Merged Ozone Data Sets
!    Revision 04 
!    DATA THROUGH: SEP 2008
!    LAST MODIFIED: 20 OCT 2008
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
!  06 Mar 2012 - R. Yantosca - Initial version, pulled code out from 
!                              the FAST-J routine SET_PROF; based on the
!                              GEOS-Chem column code routine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: I, J

      ! Initialize
      TO3_DAILY = 0d0

      !=================================================================
      ! Now weight the O3 column by the observed monthly mean TOMS.
      ! Missing data is denoted by the flag -999. (mje, bmy, 7/15/03)
      !=================================================================
      IF ( USE_O3_FROM_MET ) THEN

         !---------------------------------------------------------------
         ! Here we are using the overhead O3 from the meteorology;
         ! we won't overwrite this with TOMS/SBUV O3 columns
         !---------------------------------------------------------------
         IF ( FIRST ) THEN
            WRITE( 6, '(a)' ) REPEAT( '%', 45 )
            WRITE( 6, 100   )
 100        FORMAT( '%%% USING O3 COLUMNS FROM THE MET FIELDS! %%% ' )
            WRITE( 6, '(a)' ) REPEAT( '%', 45 )
            FIRST = .FALSE.
         ENDIF

         ! Get the overhead O3 column directly from the met field O3
         TO3_DAILY = TO3

      ELSE

         !---------------------------------------------------------------
         ! Here we are returning the default FAST-J overhead O3
         ! climatology with the TOMS/SBUV O3 columns (where data exists)
         !---------------------------------------------------------------
         IF ( DAY <= 15 ) THEN 

            ! Interpolate O3 to current day (w/in 1st half of month)
!$OMP PARALLEL DO
!$OMP+PRIVATE( I, J )
!$OMP+DEFAULT( SHARED )
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               IF ( TOMS(I,J)   > -999d0  .AND. 
     &              DTOMS1(I,J) > -999d0 ) THEN  
                  TO3_DAILY(I,J) = TOMS(I,J) 
     &                           + DTOMS1(I,J) * ( DAY - 15 )
               ENDIF            
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ELSE
  
            ! Interpolate O3 to current day (w/in 2nd half of month)
!$OMP PARALLEL DO
!$OMP+PRIVATE( I, J )
!$OMP+DEFAULT( SHARED )
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               IF ( TOMS(I,J)   > -999d0  .AND. 
     &              DTOMS2(I,J) > -999d0 ) THEN  
                  TO3_DAILY(I,J) = TOMS(I,J) 
     &                           + DTOMS2(I,J) * ( DAY - 15 )
               ENDIF
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

         ENDIF

      ENDIF

      END SUBROUTINE COMPUTE_OVERHEAD_O3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_overhead_O3
!
! !DESCRIPTION: Function GET\_OVERHEAD\_O3 returns the total overhead O3
!  column [DU] (which is taken either from TOMS/SBUV or directly from the
!  met fields) at a given surface grid box location (I,J).
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_OVERHEAD_O3( I, J ) RESULT( OVERHEAD_O3 )
!
! !INPUT PARAMETERS: 
!
      INTEGER :: I             ! Grid box longitude index
      INTEGER :: J             ! Grid box latitude index
!
! !RETURN VALUE:
!
      REAL*8  :: OVERHEAD_O3   ! Total overhead O3 column [DU]
! 
! !REVISION HISTORY: 
!  06 Mar 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Return overhead O3 to calling routine
      OVERHEAD_O3 = TO3_DAILY(I,J)

      END FUNCTION GET_OVERHEAD_O3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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

      USE CMN_SIZE_MOD  ! Size parameters
! 
! !REVISION HISTORY: 
!  14 Jul 2003 - R. Yantosca - Initial version
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!  06 Mar 2012 - R. Yantosca - Now allocate TO3_DAILY
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

      ALLOCATE( TO3_DAILY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TO3_DAILY' )
      TO3_DAILY = 0d0

      END SUBROUTINE INIT_TOMS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!  06 Mar 2012 - R. Yantosca - Now deallocate TO3_DAILY
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_TOMS begins here!
      !=================================================================
      IF ( ALLOCATED( TOMS      ) ) DEALLOCATE( TOMS      )
      IF ( ALLOCATED( DTOMS1    ) ) DEALLOCATE( DTOMS1    )
      IF ( ALLOCATED( DTOMS2    ) ) DEALLOCATE( DTOMS2    )
      IF ( ALLOCATED( TO3_DAILY ) ) DEALLOCATE( TO3_DAILY )

      END SUBROUTINE CLEANUP_TOMS
!EOC
      END MODULE TOMS_MOD
