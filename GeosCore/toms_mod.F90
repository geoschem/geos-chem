!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: toms_mod.F90
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
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: READ_TOMS
  PUBLIC :: COMPUTE_OVERHEAD_O3
  PUBLIC :: GET_OVERHEAD_O3
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
!  See https://github.com/geoschem/geos-chem for complete history
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
! !IROUTINE: read_toms
!
! !DESCRIPTION: Subroutine READ\_TOMS reads in TOMS O3 column data from a
!  binary punch file for the given grid, month and year.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_TOMS( Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Interface_Mod,  ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
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
!  "toms_mod.F90".  Missing values (i.e. locations where there are no data)
!  in the TOMS/SBUV O3 columns are defined by the flag -999.
!                                                                             .
!  After being read from disk in routine READ_TOMS, the TOMS/SBUV O3 data
!  are then passed to the FAST-J routine "set_prof.F90".  In "set_prof.F90", a
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
!  See https://github.com/geoschem/geos-chem for complete history
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
    ThisLoc = ' -> at READ_TOMS (in module GeosCore/toms_mod.F90)'

    !-----------------------------------------------------------------
    ! Read TOMS O3 columns [dobsons]
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( HcoState, 'TOMS_O3_COL', State_Chm%TOMS, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to field: TOMS_O3_COL'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Read TOMS O3 columns first day [dobsons]
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( HcoState, 'TOMS1_O3_COL', State_Chm%TOMS1, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to field: TOMS1_O3_COL!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Read TOMS O3 columns last day [dobsons]
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( HcoState, 'TOMS2_O3_COL', State_Chm%TOMS2, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to field: TOMS2_O3_COL!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Read d(TOMS)/dt, 1st half of the month [dobsons/day]
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( HcoState, 'DTOMS1_O3_COL', State_Chm%DTOMS1, RC)
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to field: DTOMS1_O3_COL'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Read d(TOMS)/dt, 2nd half of the month [dobsons/day]
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( HcoState,'DTOMS2_O3_COL', State_Chm%DTOMS2, RC )
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
!        the FAST-J routine set\_prof.F90 to overwrite the existing FAST-J
!        climatology (cf McPeters \& Nagatani 1992).  Missing data (i.e.
!        for months \& locations where TOMS/SBUV data does not exist)
!        is denoted by the value -999; FAST-J will skip over these points.
!  \item Overhead O3 columns taken directly from the met fields.  These
!        will be returned if the flag  USE\_O3\_FROM\_MET is set to TRUE.
!  \end{enumerate}
!
! !INTERFACE:
!
  SUBROUTINE COMPUTE_OVERHEAD_O3( Input_Opt, State_Grid, State_Chm, DAY, &
                                  USE_O3_FROM_MET, TO3 )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Chm_Mod,  ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt       ! Input Options object
    TYPE(GrdState), INTENT(IN) :: State_Grid      ! Grid State object
    TYPE(ChmState), INTENT(IN) :: State_Chm       ! Chemistry State object
    INTEGER,        INTENT(IN) :: DAY             ! Day of month
    LOGICAL,        INTENT(IN) :: USE_O3_FROM_MET ! Use TO3 directly from met?
    REAL(fp),       INTENT(IN) :: TO3(State_Grid%NX,State_Grid%NY) ! Met TO3
                                                                   ! [Dobsons]
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
!  READ_TOMS in "toms_mod.F90".  Missing values (i.e. locations where
!  there are no data) in the TOMS/SBUV O3 columns are defined by
!  the flag -999.
!                                                                             .
!  After being read from disk in routine READ_TOMS, the TOMS/SBUV
!  O3 data are then passed to the FAST-J routine "set_prof.F90".  In
!  "set_prof.F90", a test is done to make sure that the TOMS/SBUV O3
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE :: FIRST = .TRUE.
    INTEGER       :: I, J

    ! Initialize
    State_Chm%TO3_DAILY = 0e+0_fp

    !=================================================================
    ! Now weight the O3 column by the observed monthly mean TOMS.
    ! Missing data is denoted by the flag -999. (mje, bmy, 7/15/03)
    !=================================================================
    IF ( USE_O3_FROM_MET ) THEN

       !---------------------------------------------------------------
       ! Here we are using the overhead O3 from the meteorology;
       ! we won't overwrite this with TOMS/SBUV O3 columns
       !---------------------------------------------------------------
       IF ( FIRST .and. Input_Opt%amIRoot ) THEN
          WRITE( 6, '(a)' ) REPEAT( '%', 45 )
          WRITE( 6, 100   )
100       FORMAT( '%%% USING O3 COLUMNS FROM THE MET FIELDS! %%% ' )
          WRITE( 6, '(a)' ) REPEAT( '%', 45 )
          FIRST = .FALSE.
       ENDIF

       ! Get the overhead O3 column directly from the met field O3
       State_Chm%TO3_DAILY = TO3

    ELSE

       !---------------------------------------------------------------
       ! Here we are returning the default FAST-J overhead O3
       ! climatology with the TOMS/SBUV O3 columns (where data exists)
       !---------------------------------------------------------------
       ! Calc difference
       !$OMP PARALLEL DO     &
       !$OMP PRIVATE( I, J ) &
       !$OMP DEFAULT( SHARED )
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          State_Chm%STOMS(I,J) = (State_Chm%TOMS2(I,J)-State_Chm%TOMS1(I,J))/30.0_fp
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       ! Interpolate O3 to current day (w/in 2nd half of month)
       !$OMP PARALLEL DO     &
       !$OMP PRIVATE( I, J ) &
       !$OMP DEFAULT( SHARED )
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          State_Chm%TO3_DAILY(I,J) = State_Chm%TOMS1(I,J) + (DAY - 1) * State_Chm%STOMS(I,J)
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
  FUNCTION GET_OVERHEAD_O3( State_Chm, I, J ) RESULT( OVERHEAD_O3 )
!
! !USES:
!
    USE State_Chm_Mod,  ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
    INTEGER :: I             ! Grid box longitude index
    INTEGER :: J             ! Grid box latitude index
!
! !RETURN VALUE:
!
    REAL(fp)  :: OVERHEAD_O3   ! Total overhead O3 column [DU]
!
! !REVISION HISTORY:
!  06 Mar 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    OVERHEAD_O3 = State_Chm%TO3_DAILY(I,J)

  END FUNCTION GET_OVERHEAD_O3
END MODULE TOMS_MOD
