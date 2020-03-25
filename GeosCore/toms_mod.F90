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
  USE ErrCode_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: INIT_TOMS
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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Arrays
  REAL(fp), PRIVATE, ALLOCATABLE :: TO3_DAILY(:,:)

CONTAINS
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
  SUBROUTINE COMPUTE_OVERHEAD_O3( Input_Opt, State_Grid, DAY, &
                                  USE_O3_FROM_MET, TO3, RC )
!
! !USES:
!
    USE Input_Opt_Mod,     ONLY : OptInput
    USE State_Grid_Mod,    ONLY : GrdState
    USE HCO_Calc_Mod,      ONLY : Hco_EvalFld
    USE HCO_Interface_Mod, ONLY : HcoState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt       ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid      ! Grid State object
    INTEGER,        INTENT(IN)  :: DAY             ! Day of month
    LOGICAL,        INTENT(IN)  :: USE_O3_FROM_MET ! Use TO3 directly from met?
    REAL(fp),       INTENT(IN)  :: TO3(State_Grid%NX,State_Grid%NY) ! Met TO3
                                                                    ! [Dobsons]
!
! !OUTPUT PARAMETERS:
!
      INTEGER,      INTENT(OUT) :: RC              ! Success or failure?!
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
    LOGICAL, SAVE         :: FIRST = .TRUE.
    INTEGER               :: I, J
    CHARACTER(LEN=255)    :: ErrMsg
    REAL(fp), ALLOCATABLE :: TOMS1(:,:)
    REAL(fp), ALLOCATABLE :: TOMS2(:,:)

    ! Assume success
    RC = GC_SUCCESS

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
       IF ( FIRST .and. Input_Opt%amIRoot ) THEN
          WRITE( 6, '(a)' ) REPEAT( '%', 45 )
          WRITE( 6, 100   )
100       FORMAT( '%%% USING O3 COLUMNS FROM THE MET FIELDS! %%% ' )
          WRITE( 6, '(a)' ) REPEAT( '%', 45 )
          FIRST = .FALSE.
       ENDIF

       ! Get the overhead O3 column directly from the met field O3
       TO3_DAILY = TO3

    ELSE

       ! Evalulate the first day TOMS O3 columns from HEMCO
       ALLOCATE( TOMS1( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'toms_mod.F: TOMS1', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       TOMS1 = 0.0_fp
       CALL HCO_EvalFld( Input_Opt%amIRoot, HcoState, 'TOMS1_O3_COL', &
                         TOMS1, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find TOMS1_O3_COL in HEMCO data list!'
          CALL GC_Error( ErrMsg, RC, 'toms_mod.F' )
          RETURN
       ENDIF
       
       ! Evalulate the last day TOMS O3 columns from HEMCO
       ALLOCATE( TOMS2( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'toms_mod.F: TOMS2', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       TOMS2 = 0.0_fp
       CALL HCO_EvalFld( Input_Opt%amIRoot, HcoState, 'TOMS2_O3_COL', &
                         TOMS2, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Could not find TOMS2_O3_COL in HEMCO data list!'
          CALL GC_Error( ErrMsg, RC, 'toms_mod.F' )
          RETURN
       ENDIF

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
          STOMS(I,J) = (TOMS2(I,J)-TOMS1(I,J))/30.0_fp
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       ! Interpolate O3 to current day (w/in 2nd half of month)
       !$OMP PARALLEL DO     &
       !$OMP PRIVATE( I, J ) &
       !$OMP DEFAULT( SHARED )
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          TO3_DAILY(I,J) = TOMS1(I,J) + (DAY - 1) &
                           * ( (TOMS2(I,J)-TOMS1(I,J))/30.0_fp )
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       ! Clean up
       IF ( ALLOCATED( TOMS1 ) ) DEALLOCATE ( TOMS1 )
       IF ( ALLOCATED( TOMS2 ) ) DEALLOCATE ( TOMS2 )

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
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

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
  SUBROUTINE INIT_TOMS( Input_Opt,  State_Chm, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Failure or success
!
! !REVISION HISTORY:
!  14 Jul 2003 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
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
       ALLOCATE( TO3_DAILY( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'toms_mod.F90:TO3_DAILY', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       TO3_DAILY = 0.0_fp
    ENDIF

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
  SUBROUTINE CLEANUP_TOMS( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  14 Jul 2003 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
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
       CALL GC_CheckVar( 'toms_mod.F90:TO3_DAILY', 2, RC )
       RETURN
    ENDIF

  END SUBROUTINE CLEANUP_TOMS
!EOC
END MODULE TOMS_MOD
