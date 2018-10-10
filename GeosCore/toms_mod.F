!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: toms_mod.F
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
      USE CMN_SIZE_MOD     ! Size parameters
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

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
      INTEGER, PUBLIC, PARAMETER :: LAST_TOMS_YEAR  = 2010
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
!  21 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  16 Mar 2015 - R. Yantosca - Retire bpch input; now read TOMS via HEMCO
!  03 Nov 2016 - B. Henderson- Added fix for weird O3 in GEOS-5 simulations
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
      ! Arrays
      REAL(fp), PRIVATE, ALLOCATABLE :: TO3_DAILY(:,:)
      REAL(f4), PRIVATE, ALLOCATABLE :: STOMS(:,:)               

      ! Pointers to fields in the HEMCO data structure
      REAL(f4), PRIVATE, POINTER     :: TOMS(:,:)               
      REAL(f4), PRIVATE, POINTER     :: TOMS1(:,:)               
      REAL(f4), PRIVATE, POINTER     :: TOMS2(:,:)               
      REAL(f4), PRIVATE, POINTER     :: DTOMS1(:,:)             
      REAL(f4), PRIVATE, POINTER     :: DTOMS2(:,:)     

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
      SUBROUTINE READ_TOMS( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE HCO_Interface_Mod,  ONLY : HcoState
      USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
      USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
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
!  ##########################################################################
!  #####    NOTE: BINARY PUNCH INPUT IS BEING PHASED OUT.  THIS DATA    #####
!  #####    WILL EVENTUALLY BE READ IN FROM netCDF FILES VIA HEMCO!     #####
!  #####       -- Bob Yantosca (05 Mar 2015)                            #####
!  ##########################################################################
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
!  05 Mar 2014 - M. Sulprizio- Now regrid TOMS O3 from the files on GEOS 1x1
!                              grid (J. Lin)
!  20 Jun 2014 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  23 Jun 2014 - R. Yantosca - Now use Input_Opt%DATA_DIR_1x1
!  06 Nov 2014 - R. Yantosca - Replace TRANSFER_2D with direct casts
!  17 Dec 2014 - R. Yantosca - Leave time/date variables as 8-byte
!  05 Mar 2015 - R. Yantosca - Now read data w/r/t ExtData/CHEM_INPUTS
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Strings     
      CHARACTER(LEN=255) :: ThisLoc 
      CHARACTER(LEN=255) :: ErrMsg

      !=================================================================
      ! READ_TOMS begins here
      !=================================================================
      
      ! Assume success
      RC = GC_SUCCESS
      
      ! Exit if we are not using TOMS overhead O3 columns
      IF ( .not. Input_Opt%USE_TOMS_O3 ) RETURN

      ! Initialize
      ErrMsg  = ''
      ThisLoc = ' -> at READ_TOMS (in module GeosCore/toms_mod.F)'

      !-----------------------------------------------------------------
      ! Read TOMS O3 columns [dobsons]
      !-----------------------------------------------------------------
      CALL HCO_GetPtr( am_I_Root, HcoState, 'TOMS_O3_COL', TOMS, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to field: TOMS_O3_COL'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !-----------------------------------------------------------------
      ! Read TOMS O3 columns first day [dobsons]
      !-----------------------------------------------------------------
      CALL HCO_GetPtr( am_I_Root, HcoState, 'TOMS1_O3_COL', TOMS1, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to field: TOMS1_O3_COL!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !-----------------------------------------------------------------
      ! Read TOMS O3 columns last day [dobsons]
      !-----------------------------------------------------------------
      CALL HCO_GetPtr( am_I_Root, HcoState, 'TOMS2_O3_COL', TOMS2, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to field: TOMS2_O3_COL!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !-----------------------------------------------------------------
      ! Read d(TOMS)/dt, 1st half of the month [dobsons/day]
      !-----------------------------------------------------------------
      CALL HCO_GetPtr( am_I_Root, HcoState, 'DTOMS1_O3_COL', DTOMS1, RC)
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to field: DTOMS1_O3_COL'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      !-----------------------------------------------------------------
      ! Read d(TOMS)/dt, 2nd half of the month [dobsons/day]
      !-----------------------------------------------------------------
      CALL HCO_GetPtr( am_I_Root, HcoState,'DTOMS2_O3_COL', DTOMS2, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Cannot get pointer to field: DTOMS2_O3_COL!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

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
      SUBROUTINE COMPUTE_OVERHEAD_O3( am_I_Root,       DAY,
     &                                USE_O3_FROM_MET, TO3 )
!
! !INPUT PARAMETERS: 
!
      ! Root CPU?
      LOGICAL, INTENT(IN)  :: am_I_Root

      ! Day of month 
      INTEGER, INTENT(IN)  :: DAY      

      ! Switch to denote if we should use the default TO3
      ! directly from the met fields
      LOGICAL, INTENT(IN)  :: USE_O3_FROM_MET
      
      ! TO3 from the met fields [Dobsons]
      REAL(fp),  INTENT(IN)  :: TO3(IIPAR,JJPAR)
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
!  22 Oct 2014 - C. Keller   - Added am_I_Root parameter.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: I, J

      ! Initialize
      TO3_DAILY = 0e+0_fp

      !=================================================================
      ! Now weight the O3 column by the observed monthly mean TOMS.
      ! Missing data is denoted by the flag -999. (mje, bmy, 7/15/03)
      !=================================================================
      IF ( USE_O3_FROM_MET ) THEN

         !---------------------------------------------------------------
         ! Here we are using the overhead O3 from the meteorology;
         ! we won't overwrite this with TOMS/SBUV O3 columns
         !---------------------------------------------------------------
         IF ( am_I_Root .AND. FIRST ) THEN
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
            ! Calc difference
!$OMP PARALLEL DO
!$OMP+PRIVATE( I, J )
!$OMP+DEFAULT( SHARED )
            DO J = 1, JJPAR
            DO I = 1, IIPAR
                 STOMS(I,J) = (TOMS2(I,J)-TOMS1(I,J))/30.0_fp
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
            ! Interpolate O3 to current day (w/in 2nd half of month)
!$OMP PARALLEL DO
!$OMP+PRIVATE( I, J )
!$OMP+DEFAULT( SHARED )
            DO J = 1, JJPAR
            DO I = 1, IIPAR
                TO3_DAILY(I,J) = TOMS1(I,J) + (DAY - 1) * STOMS(I,J)
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

      ENDIF

      END SUBROUTINE COMPUTE_OVERHEAD_O3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_overhead_o3
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
      REAL(fp)  :: OVERHEAD_O3   ! Total overhead O3 column [DU]
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
      SUBROUTINE INIT_TOMS( am_I_Root, Input_Opt, 
     &                      State_Chm, State_Diag, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD  
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Diag_Mod,     ONLY : DgnState
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)           :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC         ! Failure or success
! 
! !REVISION HISTORY: 
!  14 Jul 2003 - R. Yantosca - Initial version
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!  06 Mar 2012 - R. Yantosca - Now allocate TO3_DAILY
!  16 Mar 2015 - R. Yantosca - TOMS, DTOMS1, DTOMS2 are now pointers
!  16 Mar 2015 - R. Yantosca - Add Input_Opt, RC arguments for std interface
!  03 Nov 2016 - B. Henderson- Add extra fields to fix OH issue for GEOS-5
!  07 Aug 2018 - H.P. Lin    - Now accepts State_Chm, State_Diag to unify input
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !=================================================================
      ! INIT_TOMS begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS
      
      ! Allocate arrays
      IF ( .not. ALLOCATED( TO3_DAILY ) ) THEN
         ALLOCATE( TO3_DAILY( IIPAR, JJPAR ), STAT=RC )        
         CALL GC_CheckVar( 'toms_mod.F:TO3_DAILY', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         TO3_DAILY = 0.0_fp
      ENDIF

      IF ( .not. ALLOCATED( STOMS ) ) THEN
         ALLOCATE( STOMS( IIPAR, JJPAR ), STAT=RC )
         CALL GC_CheckVar( 'toms_mod.F:STOMS', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         STOMS = 0e+0_fp
      ENDIF
         
      ! Initialize pointers
      TOMS   => NULL()
      TOMS1  => NULL()
      TOMS2  => NULL()
      DTOMS1 => NULL()
      DTOMS2 => NULL()

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
      SUBROUTINE CLEANUP_TOMS( am_I_Root, RC )
!
! !USES:
!
      USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  14 Jul 2003 - R. Yantosca - Initial version
!  01 Dec 2010 - R. Yantosca - Added ProTeX headers
!  06 Mar 2012 - R. Yantosca - Now deallocate TO3_DAILY
!  03 Nov 2016 - R. Yantosca - Deallocate and nullify extra fields
!  10 Oct 2018 - R. Yantosca - Add am_I_Root, RC; Improve error trapping
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP_TOMS begins here!
      !=================================================================

      ! Assume success
      RC = GC_SUCCESS

      ! Deallocate variables
      IF ( ALLOCATED( TO3_DAILY ) ) THEN
         DEALLOCATE( TO3_DAILY, STAT=RC )     
         CALL GC_CheckVar( 'toms_mod.F:TO3_DAILY', 2, RC )
         RETURN
      ENDIF

      IF ( ALLOCATED( STOMS ) ) THEN
         DEALLOCATE( STOMS, STAT=RC )     
         CALL GC_CheckVar( 'toms_mod.F:STOMS', 2, RC )
         RETURN
      ENDIF

      ! Free pointers
      TOMS   => NULL()
      TOMS1  => NULL()
      TOMS2  => NULL()
      DTOMS1 => NULL()
      DTOMS2 => NULL()

      END SUBROUTINE CLEANUP_TOMS
!EOC
      END MODULE TOMS_MOD
