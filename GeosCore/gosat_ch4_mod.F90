!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gosat_ch4_mod.F90
!
! !DESCRIPTION: Module GOSAT\_CH4\_MOD
!\\
!\\
! !INTERFACE:
!
MODULE GOSAT_CH4_MOD
!
! !USES:
!
  USE m_netcdf_io_open       ! netCDF open
  USE m_netcdf_io_get_dimlen ! netCDF dimension queries
  USE m_netcdf_io_read       ! netCDF data reads
  USE m_netcdf_io_close      ! netCDF close
  USE PRECISION_MOD          ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CALC_GOSAT_CH4_FORCE
!
! !REVISION HISTORY:
!  16 Jun 2017 - M. Sulprizio- Initial version based on GOSAT CH4 observation
!                              operator from GC Adjoint v35j with updates from
!                              M. Sulprizio, J.D. Maasakkers, A. Turner, and
!                              K. Wecht
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  INTEGER,  PARAMETER  :: MAXLEV = 20
  INTEGER,  PARAMETER  :: MAXGOS = 100000
  INTEGER,  PARAMETER  :: NT_FD = 335 !1055 ! Index for the FD test
  REAL(fp), PARAMETER  :: GC_XCH4_ERROR = 0e+0_fp !ppb
  REAL(fp), PARAMETER  :: PMAG = 1e+0_fp ! Perturbation magnitude (%)
  LOGICAL,  PARAMETER  :: LDCH4SAT    = .TRUE.
  LOGICAL,  PARAMETER  :: EXT_OBS_MAT = .FALSE. ! Read external?
!
! !MODULE VARIABLES
!
  ! Record to store data from each GOS obs
  TYPE GOS_CH4_OBS
     INTEGER           :: LGOS(1)
     REAL(fp)          :: LAT(1)
     REAL(fp)          :: LON(1)
     INTEGER           :: YEAR(1)
     INTEGER           :: MONTH(1)
     INTEGER           :: DAY(1)
     INTEGER           :: HOUR(1)
     INTEGER           :: MINUTE(1)
     INTEGER           :: SEC(1)
     REAL(fp)          :: TIME(1)
     REAL(fp)          :: CH4(1)
     REAL(fp)          :: CH4_ERROR(1)
     REAL(fp)          :: PRES(MAXLEV)
     REAL(fp)          :: PRIOR(MAXLEV)
     REAL(fp)          :: AVG_KERNEL(MAXLEV)
     REAL(fp)          :: P_WEIGHT(MAXLEV)
     INTEGER           :: QFLAG(1)
     INTEGER           :: GLINT(1)
     INTEGER           :: GAIN(1)
     CHARACTER(LEN=22) :: EXP_ID(1)
  ENDTYPE GOS_CH4_OBS

  TYPE(GOS_CH4_OBS)    :: GOS(MAXGOS)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_gos_ch4_obs
!
! !DESCRIPTION: Subroutine READ\_GOS\_CH4\_OBS reads the file and passes back
!  info contained therein. (dkh, 10/12/10)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_GOS_CH4_OBS( YYYYMMDD, NGOS )
!
! !USES:
!
    USE TIME_MOD,        ONLY : EXPAND_DATE
    USE ERROR_MOD,       ONLY : ALLOC_ERR
    USE ERROR_MOD,       ONLY : GEOS_CHEM_STOP
    USE TIME_MOD,        ONLY : GET_YEAR
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)      :: YYYYMMDD  ! Current date
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT)     :: NGOS      ! Number of GOS retrievals
!
! !REVISION HISTORY:
!  (1 )Based on READ_TES_NH3 OBS (dkh, 04/26/10)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: LGOS
    INTEGER                        :: N, J
    INTEGER                        :: YEAR
    CHARACTER(LEN=4)               :: CYEAR
    LOGICAL                        :: EXIST_VAR
    REAL*4, ALLOCATABLE            :: pres(:,:)
    REAL*4, ALLOCATABLE            :: prior(:,:)
    REAL*4, ALLOCATABLE            :: ak(:,:)
    REAL*4, ALLOCATABLE            :: pres_w(:,:)
    REAL*4, ALLOCATABLE            :: lat(:)
    REAL*4, ALLOCATABLE            :: lon(:)
    REAL*4, ALLOCATABLE            :: time(:)
    REAL*4, ALLOCATABLE            :: xch4(:)
    REAL*4, ALLOCATABLE            :: xch4_error(:)
    INTEGER,  ALLOCATABLE          :: qflag(:)
    INTEGER,  ALLOCATABLE          :: glint(:)
    INTEGER,  ALLOCATABLE          :: gain(:)
    CHARACTER(LEN=22),ALLOCATABLE  :: exp_id(:)

    INTEGER                        :: NLEV, StrLen
    INTEGER                        :: I, L, AS

    ! For reading netCDF file
    INTEGER            :: fId              ! netCDF file ID
    INTEGER            :: Status           ! 0 means variable in file
    INTEGER            :: X, Y, Z, T       ! netCDF file dimensions
    INTEGER            :: time_index       ! Read this slice of data
    INTEGER            :: st1d(1), ct1d(1) ! Start + count for 1D arrays
    INTEGER            :: st2d(2), ct2d(2) ! Start + count for 2D arrays
    CHARACTER(LEN=16)  :: stamp            ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file          ! netCDF file name
    CHARACTER(LEN=255) :: v_name           ! netCDF variable name
    CHARACTER(LEN=255) :: dir              ! Data directory path
    CHARACTER(LEN=255) :: errMsg           ! Error message
    CHARACTER(LEN=255) :: caller           ! Name of this routine

    !=================================================================
    ! READ_GOS_CH4_OBS begins here!
    !=================================================================

    caller = 'READ_GOS_CH4_OBS in gosat_ch4_mod.F90'

    ! Get current year
    YEAR = GET_YEAR()
    WRITE( CYEAR, '(i4)' ) YEAR

    ! Filename
    nc_file = 'ESACCI-GHG-L2-CH4-GOSAT-OCPR-YYYYMMDD-fv7.0.nc'
    CALL Expand_Date( nc_file, YYYYMMDD, 9999 )

    ! Construct complete file path
    dir = '/n/regal/jacob_lab/msulprizio/ExtData/Obs/ULGOSAT_v7/' // &
           CYEAR // '/'
    nc_file = TRIM( dir ) // TRIM( nc_file )
    WRITE( 6, 10 ) TRIM( nc_file )
10  FORMAT( '     - Reading ', a)

    ! Make sure the file exists (ajt, 03/31/2013)
    INQUIRE( FILE=TRIM( nc_file ), EXIST=EXIST_VAR )
    IF ( .not. EXIST_VAR ) THEN
       NGOS = -1
       RETURN
    ENDIF

    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'm',   NLEV )
    CALL NcGet_DimLen( fId, 'n',   NGOS )
    CALL NcGet_DimLen( fId, 'string_length', StrLen )

    IF ( NLEV > MAXLEV ) THEN
       print*,' # Levels this day = ', NLEV
       print*, 'WARNING: NLEV > MAXLEV. Need to increase'
       print*, ' MAXLEV in gosat_ch4_mod.F90.'
       CALL GEOS_CHEM_STOP
    ENDIF

    print*,' # GOSAT Observations this day = ', NGOS
    print*, 'levels', NLEV

    !--------------------------------
    ! Read 1-D Data
    !--------------------------------

    ALLOCATE( lat(         NGOS ), STAT=AS )
    ALLOCATE( lon(         NGOS ), STAT=AS )
    ALLOCATE( time(        NGOS ), STAT=AS )
    ALLOCATE( qflag(       NGOS ), STAT=AS )
    ALLOCATE( xch4(        NGOS ), STAT=AS )
    ALLOCATE( xch4_error(  NGOS ), STAT=AS )
    ALLOCATE( glint(       NGOS ), STAT=AS )
    ALLOCATE( gain(        NGOS ), STAT=AS )
    ALLOCATE( exp_id(      NGOS ), STAT=AS )
    ALLOCATE( prior( NLEV, NGOS ), STAT=AS )
    ALLOCATE( pres(  NLEV, NGOS ), STAT=AS )
    ALLOCATE( ak(    NLEV, NGOS ), STAT=AS )
    ALLOCATE( pres_w(NLEV, NGOS ), STAT=AS )

    ! Latitude
    v_name = 'latitude'
    st1d   = (/ 1    /)
    ct1d   = (/ NGOS /)
    CALL NcRd( lat, fId, TRIM(v_name), st1d, ct1d )

    ! Longitude
    v_name = 'longitude'
    st1d   = (/ 1    /)
    ct1d   = (/ NGOS /)
    CALL NcRd( lon, fId, TRIM(v_name), st1d, ct1d )

    ! Time (seconds since 1970-01-01 00:00)
    v_name = 'time'
    st1d   = (/ 1    /)
    ct1d   = (/ NGOS /)
    CALL NcRd( time, fId, TRIM(v_name), st1d, ct1d )

    ! Qflag (0=good, 1=bad)
    v_name = 'xch4_quality_flag'
    st1d   = (/ 1    /)
    ct1d   = (/ NGOS /)
    CALL NcRd( qflag, fId, TRIM(v_name), st1d, ct1d )

    ! Xch4 (ppb)
    v_name = 'xch4'
    st1d   = (/ 1    /)
    ct1d   = (/ NGOS /)
    CALL NcRd( xch4, fId, TRIM(v_name), st1d, ct1d )

    ! Xch4_Error (ppb)
    v_name = 'xch4_uncertainty'
    st1d   = (/ 1    /)
    ct1d   = (/ NGOS /)
    CALL NcRd( xch4_error, fId, TRIM(v_name), st1d, ct1d )

    ! Retrieval type flag (0 = land, 1 = glint)
    v_name = 'retr_flag'
    st1d   = (/ 1    /)
    ct1d   = (/ NGOS /)
    CALL NcRd( glint, fId, TRIM(v_name), st1d, ct1d )

    ! Gain mode (1 = high gain, 0 = medium gain)
    v_name = 'gain'
    st1d   = (/ 1    /)
    ct1d   = (/ NGOS /)
    CALL NcRd( gain, fId, TRIM(v_name), st1d, ct1d )

    ! Exposure ID
    v_name = 'exposure_id'
    st2d   = (/ 1, 1    /)
    ct2d   = (/ StrLen, NGOS /)
    CALL NcRd( exp_id, fId, TRIM(v_name), st2d, ct2d )

    !--------------------------------
    ! Read 2D Data
    !--------------------------------

    ! APRIORI (ppb)
    v_name = 'ch4_profile_apriori'
    st2d   = (/ 1,    1    /)
    ct2d   = (/ NLEV, NGOS /)
    CALL NcRd( prior, fId, TRIM(v_name), st2d, ct2d )

    ! Pressure (hPa)
    v_name = 'pressure_levels'
    st2d   = (/ 1,    1    /)
    ct2d   = (/ NLEV, NGOS /)
    CALL NcRd( pres, fId, TRIM(v_name), st2d, ct2d )

    ! Averaging Kernel (unitless)
    v_name = 'xch4_averaging_kernel'
    st2d   = (/ 1,    1    /)
    ct2d   = (/ NLEV, NGOS /)
    CALL NcRd( ak, fId, TRIM(v_name), st2d, ct2d )

    ! Pressure weights (unitless)
    v_name = 'pressure_weight'
    st2d   = (/ 1,    1    /)
    ct2d   = (/ NLEV, NGOS /)
    CALL NcRd( pres_w, fId, TRIM(v_name), st2d, ct2d )

    !--------------------------------
    ! Place data into GOS structure
    !--------------------------------
    DO N = 1, NGOS

       ! 0-D data
       GOS(N)%LGOS(1)      = NLEV
       GOS(N)%LAT(1)       = lat(N)
       GOS(N)%LON(1)       = lon(N)
       GOS(N)%TIME(1)      = time(N)
       GOS(N)%CH4(1)       = xch4(N)
       GOS(N)%CH4_ERROR(1) = xch4_error(N)
       GOS(N)%QFLAG(1)     = qflag(N)
       GOS(N)%GLINT(1)     = glint(N)
       GOS(N)%GAIN(1)      = gain(N)
       GOS(N)%EXP_ID(1)    = exp_id(N)

       ! 1-D data
       LGOS = NLEV
       GOS(N)%PRES(1:LGOS)       = pres(1:LGOS,N)
       GOS(N)%PRIOR(1:LGOS)      = prior(1:LGOS,N)
       GOS(N)%AVG_KERNEL(1:LGOS) = ak(1:LGOS,N)
       GOS(N)%P_WEIGHT(1:LGOS)  = pres_w(1:LGOS,N)

    ENDDO

    !--------------------------------
    ! Close netCDF file
    !--------------------------------

    ! Echo info
    WRITE( 6, 20 ) YYYYMMDD
20  FORMAT( '     - Finished reading GOSAT CH4 observations for ', i8)

    ! Close netCDF file
    CALL NcCl( fId )
    
  END SUBROUTINE READ_GOS_CH4_OBS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_gosat_ch4_force
!
! !DESCRIPTION: Subroutine CALC\_GOS\_CH4\_FORCE calculates the adjoint forcing
!  from the GOSAT CH4 observations and updates the cost function.
!  (dkh, 10/12/10)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_GOSAT_CH4_FORCE( Input_Opt, State_Chm, State_Grid, &
                                   State_Met )
!
! !USES:
!
    USE ERROR_MOD,          ONLY : IT_IS_NAN
    USE ERROR_MOD,          ONLY : IT_IS_FINITE
    USE GC_GRID_MOD,        ONLY : GET_IJ
    USE Input_Opt_Mod,      ONLY : OptInput
    USE JULDAY_MOD,         ONLY : CALDATE,     JULDAY
    USE TIME_MOD
    USE PhysConstants,      ONLY : XNUMOLAIR, AIRMW
    USE State_Chm_Mod,      ONLY : ChmState, Ind_
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
    TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !REVISION HISTORY:
!  16 Jun 2017 - M. Sulprizio- Initial version based on GOSAT CH4 observation
!                              operator from GC Adjoint v35j with updates from
!                              M. Sulprizio, J.D. Maasakkers, A. Turner, and
!                              K. Wecht
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp), SAVE     :: COST_FUNC   ! Cost function [unitless]
    INTEGER            :: id_CH4
    INTEGER            :: NTSTART, NTSTOP, NT
    INTEGER            :: YYYYMM,  YYYYMMDD, HHMMSS
    INTEGER            :: YEAR,    MON,      DAY
    INTEGER            :: HOUR,    MIN,      SEC
    INTEGER            :: IIJJ(2), I,      J,     N
    INTEGER            :: L,       LL,     LGOS
    INTEGER            :: JLOOP,   NOBS,   IND
    INTEGER            :: INDS(MAXGOS)
    REAL(fp)           :: REF_DATE, TIME
    REAL(fp)           :: GC_PRES(State_Grid%NZ)
    REAL(fp)           :: GC_PEDGE(State_Grid%NZ+1)
    REAL(fp)           :: GC_CH4_NATIVE(State_Grid%NZ)
    REAL(fp)           :: GC_CH4(MAXLEV)
    REAL(fp)           :: GC_CH4_cm2(MAXLEV)
    REAL(fp)           :: GC_PSURF
    REAL(fp)           :: MAP(State_Grid%NZ,MAXLEV)
    REAL(fp)           :: CH4_HAT(MAXLEV)
    REAL(fp)           :: NEW_COST(MAXGOS)
    REAL(fp)           :: OLD_COST
    REAL(fp), SAVE     :: TIME_FRAC(MAXGOS)
    INTEGER,  SAVE     :: NGOS
    REAL(fp)           :: frac, RHS, LHS
    REAL(fp)           :: CH4_PRIOR(MAXLEV)
    REAL(fp)           :: CH4_PRIOR_cm2(MAXLEV)
    REAL(fp)           :: molecongos(MAXLEV)
    REAL(fp)           :: GOS_XCH4, GC_XCH4, GOS_XCH4_ERROR
    REAL(fp)           :: p(MAXLEV), h(MAXLEV)
    REAL(fp)           :: ak(MAXLEV), prior(MAXLEV)
    REAL(fp)           :: pres_w(MAXLEV)
    REAL(fp)           :: XCH4m, XCH4a
    REAL(fp)           :: SATELLITE_BIAS(3) ! Hardcode for now (mps,6/16/17)
    REAL(fp)           :: MEAN_MODEL_BIAS   ! Hardcode for now (mps,6/16/17)
    REAL(fp)           :: FORCE
    REAL(fp)           :: DIFF
    REAL(fp)           :: S_OBS

    ! --- zyz --- Sept 19, 2018
    ! fix the problem that some GOSAT record has negative pressure
    ! values if the surface has high altidue
    INTEGER            :: L0

    ! For miscellaneous
    LOGICAL, SAVE      :: FIRST = .TRUE.
    INTEGER            :: IOS
    INTEGER, SAVE      :: TotalObs = 0
    CHARACTER(LEN=255) :: FILENAME

    !=================================================================
    ! CALC_GOS_CH4_FORCE begins here!
    !=================================================================

    print*, '     - CALC_GOS_CH4_FORCE '

    ! Initialize species ID flag
    id_CH4     = Ind_('CH4'       )

    ! Reset
    NEW_COST(:) = 0.0_fp

    ! Hardcode bias values for now (mps, 6/16/17)
    SATELLITE_BIAS(1) = 0.0e+0_fp
    SATELLITE_BIAS(2) = -0.075272936e+0_fp
    SATELLITE_BIAS(3) = 0.003712500e+0_fp
    MEAN_MODEL_BIAS   = 0.0e+0_fp

    ! Open files for diagnostic output
    IF ( FIRST ) THEN

       ! For recording model and observation values
       FILENAME = 'sat_obs.gosat.00.m'
       FILENAME = TRIM( Input_Opt%RUN_DIR ) //  TRIM( FILENAME )
       OPEN( 117,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN', &
             IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       ! Write header of sat_obs.NN.m
       WRITE( 117, 281 ) '       NNN', &
                         '   I', '   J', '     LON','     LAT', &
                         'YYYY','MM', 'DD', 'hh', 'mm', 'ss',   &
                         '         TAU', '       GOSAT', &
                         '       model', '       S_OBS', &
                         '    COST_FUN', 'GLINT', 'GAIN', &
                         'EXPOSURE_ID'
281    FORMAT( A10,2x,A4,2x,A4,2x,A8,2x,A8,2x,A4,2x,A2,2x,A2,2x,A2,2x, &
               A2,2x,A2,2x,A12,2x,A12,2x,A12,2x,A12,2x,A12,2x, &
               A5,2x,A4,2x,A22)

       ! Set Total Observations = 0
       TotalObs = 0

    ENDIF

    ! Save a value of the cost function first
    OLD_COST = COST_FUNC

    ! Read Observations at first call and at end of the day
    IF ( FIRST .OR. ITS_A_NEW_DAY() ) THEN

       ! Read the GOS CH4 file for this day
       YYYYMMDD = 1d4*GET_YEAR() + 1d2*GET_MONTH() + GET_DAY()
       CALL READ_GOS_CH4_OBS( YYYYMMDD, NGOS )

       IF ( FIRST ) FIRST = .FALSE.

       ! Make sure there are observations on this day
       IF ( NGOS .EQ. -1 ) RETURN

       ! TIME is seconds since 1970-01-01 00:00:00
       ! Convert to calendar date by expressing as days since REF_DATE
       ! and adding Julian day value for REF_DATE (mpayer 9/12/12)
       REF_DATE = JULDAY( 1970, 1, 1.d0 )
       DO N = 1, NGOS
          TIME = ( GOS(N)%TIME(1) / 86400.d0 ) + REF_DATE
          CALL CALDATE ( TIME, YYYYMMDD, HHMMSS )
          CALL YMD_EXTRACT( YYYYMMDD, YEAR, MON, DAY )
          CALL YMD_EXTRACT( HHMMSS,   HOUR, MIN, SEC )
          GOS(N)%YEAR(1)   = YEAR
          GOS(N)%MONTH(1)  = MON
          GOS(N)%DAY(1)    = DAY
          GOS(N)%HOUR(1)   = HOUR
          GOS(N)%MINUTE(1) = MIN
          GOS(N)%SEC(1)    = SEC
       ENDDO

       IF ( FIRST ) FIRST = .FALSE.

    ENDIF

    ! Get indices of GOS observations in the current hour
    !   At the start of each hour, assimilate observations that
    !   were made in the previous 60 minutes.
    !   For example, at time 18:00, assimilate observations
    !   made from 18:00 - 18:59
    INDS(:) = 0
    NOBS    = 0
    !print*,'Looking for observations at MONTH, DAY, HOUR = ', &
    !          GET_MONTH(), GET_DAY(), GET_HOUR()
    DO NT = 1, NGOS
       IF ( GOS(NT)%MONTH(1) .EQ. GET_MONTH() .AND. &
            GOS(NT)%DAY(1)   .EQ. GET_DAY()   .AND. &
            GOS(NT)%HOUR(1)  .EQ. GET_HOUR()  ) THEN
          NOBS = NOBS + 1
          INDS(NOBS) = NT
          !print*,'Found a good observation! NT = ', NT
       ENDIF
    ENDDO

    IF ( NOBS == 0 ) THEN
       print*, ' No matching GOSAT CH4 obs for this hour'
       RETURN
    ENDIF

    print*, ' for hour range: ', GET_HOUR(), GET_HOUR()+1
    print*, ' found # GOSAT observations: ', NOBS


    !! need to update this in order to do i/o with this loop parallel
    !!      ! Now do a parallel loop for analyzing data
    !!!$OMP PARALLEL DO &
    !!!$OMP DEFAULT( PRIVATE ) &
    !!!$OMP PRIVATE( IND, NT, MAP, LGOS, IIJJ,  I, J,  L,   LL, JLOOP ) &
    !!!$OMP PRIVATE( GC_CH4, FORCE, CH4_PRIOR, GC_PRES, FILENAME      ) &
    !!!$OMP PRIVATE( GC_PEDGE, GC_PSURF, GC_CH4_NATIVE, GOS_XCH4      ) &
    !!!$OMP PRIVATE( GOS_XCH4_ERROR, S_OBS, h, p, XCH4a, XCH4m        ) &
    !!!$OMP PRIVATE( GC_XCH4, DIFF, DIFF_ADJ, GC_XCH4_ADJ             ) &
    !!!$OMP PRIVATE( GC_CH4_NATIVE_ADJ, GC_CH4_ADJ, TotalObs          )
    DO IND = 1, NOBS

       NT = INDS(IND)
       !print*, '     - CALC_GOS_CH4_FORCE: analyzing record ', NT

       ! Quality screening (0=good, 1=bad)
       IF ( GOS(NT)%QFLAG(1) .NE. 0 ) THEN
          print*, ' Bad QFLAG, skipping record         ', NT
          CYCLE
       ENDIF

       !! Skip glint mode observations (0 = land, 1 = glint)
       !IF ( GOS(NT)%GLINT(1) .EQ. 1 ) THEN
       !   print*, ' GLINT obs, skipping record         ', NT
       !   CYCLE
       !ENDIF

       ! Check for NaN in data
       IF ( IT_IS_NAN( GOS(NT)%CH4(1) ) ) THEN
          print*, ' XCH4 is NaN, skipping record       ', NT
          CYCLE
       ENDIF
       IF ( IT_IS_NAN( GOS(NT)%PRIOR(1) ) ) THEN
          print*, ' PRIOR is NaN, skipping record      ', NT
          CYCLE
       ENDIF
       IF ( IT_IS_NAN( GOS(NT)%AVG_KERNEL(1) ) ) THEN
          print*, ' AVG_KERNEL is NaN, skipping record ', NT
          CYCLE
       ENDIF

       ! Check for infinity in data
       IF ( .not. IT_IS_FINITE( GOS(NT)%CH4(1) ) ) THEN
          print*, ' XCH4 is infinity, skipping record  ', NT
          CYCLE
       ENDIF

       ! Check for negative/zero data
       IF ( GOS(NT)%CH4(1) <= 0d0 ) THEN
          print*, ' XCH4 is <= 0, skipping record      ', NT
          CYCLE
       ENDIF

       ! Check for negative/zero error values
       ! Skip zero error values for now to avoid returning S_OBS=0
       ! which will cause COST_FUN=infinity (mps, 6/6/16)
       IF ( GOS(NT)%CH4_ERROR(1) <= 0d0 ) THEN
          print*, ' XCH4 ERROR is <= 0, skipping record ', NT
          CYCLE
       ENDIF

       ! zyz, Oct 4 2018
       ! GOSAT record has negative pressure value when the surface
       ! has high altitude such as Tibet
       ! Fix it by finding the lowest valid layer L0
       L0 = GOS(NT)%LGOS(1)
       DO L = 1, GOS(NT)%LGOS(1)
          IF (GOS(NT)%PRES(L) .gt. 0.0_fp) THEN
             L0 = L
             EXIT
          ENDIF
       ENDDO
       ! Only find negative pressure at 1st layer so far.
       ! But use this condition to prevent really invalid cases
       ! (zyz)
       IF (L0 .gt. 3 ) CYCLE

       ! Skip observations outside the domain
       IF ( GOS(NT)%LAT(1) < State_Grid%YMin .OR. &
            GOS(NT)%LAT(1) > State_Grid%YMax .OR. &
            GOS(NT)%LON(1) < State_Grid%XMin .OR. &
            GOS(NT)%LON(1) > State_Grid%XMax ) THEN
          print*, ' Outside nested domain, skipping record ', NT
          CYCLE
       ENDIF

       ! Skip data at high lat (ajt, 7/23/13)
       !JDMIF ( GOS(NT)%LAT(1) >   60d0 ) THEN
       !JDMprint*, ' Skip all data above 60N ', NT
       !JDMCYCLE
       !JDMENDIF
       !Turn of this for now to reproduce plots made by AJT
       !Should probably turn this back on with SH once we start
       !To do the actual adjoint run

       ! Get grid box of current record
       IIJJ = GET_IJ( REAL(GOS(NT)%LON(1),4), REAL(GOS(NT)%LAT(1),4), &
                      State_Grid )
       I    = IIJJ(1)
       J    = IIJJ(2)

       ! skip observations where the GOSAT surface pressure is
       ! less than the model
       ! Use L0 for lowest valid layer in GOSAT (zyz)
       IF ( (GOS(NT)%PRES(L0) - State_Met%PEDGE(I,J,1)) > 50e0 ) THEN
          print*, ' Psurf threshold not met, skipping record ', NT
          CYCLE
       ENDIF

       !------------------------------
       ! Begin good observations
       !------------------------------
       print*,' Begin assimilating good observation. NT = ', NT

       ! For safety, initialize these up to LGOS
       LGOS            = 0
       GC_CH4(:)       = 0.0_fp
       MAP(:,:)        = 0.0_fp
       FORCE           = 0.0_fp
       CH4_PRIOR(:)    = 0.0_fp

       ! Copy variable names to make coding a bit cleaner
       ! Use L0 for lowest layer in GOSAT record (zyz)
       LGOS = GOS(NT)%LGOS(1)
       DO L = L0, LGOS
          CH4_PRIOR(L) = GOS(NT)%PRIOR(L) * 1d-9 ! [ppb] --> [v/v]
       ENDDO

       ! Get GC pressure levels (mbar)
       DO L = 1, State_Grid%NZ
          GC_PRES(L) = State_Met%PMID(I,J,L)
       ENDDO

       ! Get GC pressure edges (mbar)
       DO L = 1, State_Grid%NZ+1
          GC_PEDGE(L) = State_Met%PEDGE(I,J,L)
       ENDDO

       ! Get GC surface pressure (mbar)
       GC_PSURF = State_Met%PEDGE(I,J,1)

       ! Calculate the interpolation weight matrix
       MAP(:,:) = 0.0_fp
       CALL GET_INTMAP( State_Grid, GC_PRES, GC_PSURF, GOS(NT)%PRES, &
                        L0, LGOS, MAP )

       ! Get CH4 values at native model resolution
       GC_CH4_NATIVE(:) = 0.0_fp

       ! Get species concentrations
       ! Convert from [kg/box] --> [v/v]
       GC_CH4_NATIVE(:) = State_Chm%Species(I,J,:,id_CH4) * ( AIRMW / &
                          State_Chm%SpcData(id_CH4)%Info%emMW_g )

       ! Interpolate GC CH4 column to GOSAT grid
       ! Use L0 for lowest valid layer for GOSAT (zyz)
       DO LL = L0, LGOS
          GC_CH4(LL) = 0.0_fp
          DO L = 1, State_Grid%NZ
             GC_CH4(LL) = GC_CH4(LL) + MAP(L,LL) * GC_CH4_NATIVE(L)
          ENDDO
       ENDDO

       ! GOSAT Proxy XCH4 [ppb] --> [v/v]
       GOS_XCH4       = GOS(NT)%CH4(1)       * 1e-9_fp
       GOS_XCH4_ERROR = GOS(NT)%CH4_ERROR(1) * 1e-9_fp

       ! Remove any GOSAT bias
       GOS_XCH4 = GOS_XCH4 + 1e-9 * ( SATELLITE_BIAS(1) &
                           + SATELLITE_BIAS(2)*(GOS(NT)%LAT(1)) &
                           + SATELLITE_BIAS(3)*(GOS(NT)%LAT(1))**2 )

       ! Get the S_obs, assume stddev adds in quadrature, variance
       ! adds linearly.  (ajt, 03/27/2013)
       S_OBS = GOS_XCH4_ERROR**2 + (GC_XCH4_ERROR * 1e-9_fp)**2

       !--------------------------------------------------------------
       ! Apply GOSAT observation operator
       !
       !   Xch4_m = Xch4_a + SUM_j( h_j * a_j * (x_m - x_a) )
       !
       !   Xch4_m  - model XCH4
       !   Xch4_a  - apriori XCH4 = h^T * x_a
       !   h       - pressure weighting function
       !   a       - column averaging kernel
       !   x_m     - model CH4 [v/v]
       !   x_a     - apriori CH4 [v/v]
       !
       !   The pressure weighting function is defined in Connor et al. 2008
       !     and the OCO-2 ATBD
       !--------------------------------------------------------------

       ! Pressure weighting array
       ! Use L0 for lowest valid layer in GOSAT
       h(:)          = 0.0_fp
       p(L0:LGOS)     = GOS(NT)%PRES(L0:LGOS)
       ak(L0:LGOS)    = GOS(NT)%AVG_KERNEL(L0:LGOS)
       prior(L0:LGOS) = GOS(NT)%PRIOR(L0:LGOS) * 1e-9_fp  ! [ppb] --> [v/v]

       ! Need to integrate from the toa to surface (ajt, 05/21/13)
       IF (LGOS .GT. 1) THEN
          IF(GOS(NT)%PRES(L0+1) .LT. GOS(NT)%PRES(L0)) THEN
             p(1:LGOS) = p(LGOS:1:-1)
          ENDIF
       ENDIF

       L = 1
       h(L) = 1./GOS(NT)%PRES(L0) * ABS( &
              ( -1e0*p(L) + ( (p(L+1)-p(L))/(LOG(p(L+1)/p(L))) ) ) )
       L = LGOS - L0 + 1
       h(L) = 1./GOS(NT)%PRES(L0) * ABS( &
              (  p(L) - ( (p(L)-p(L-1))/(LOG(p(L)/p(L-1))) ) ) )
       DO L=2,LGOS-L0
          h(L) = 1./GOS(NT)%PRES(L0) * ABS( &
                 ( -1e0*p(L) + ( (p(L+1)-p(L))/(LOG(p(L+1)/p(L))) ) ) + &
                 (      p(L) - ( (p(L)-p(L-1))/(LOG(p(L)/p(L-1))) ) )   )
       ENDDO

       ! Now return to the orientation of the other variables
       IF (LGOS .GT. 1) THEN
          IF(GOS(NT)%PRES(L0+1) .LT. GOS(NT)%PRES(L0)) THEN
             h(1:LGOS) = h(LGOS:1:-1)
             p(1:LGOS) = p(LGOS:1:-1)
          ENDIF
       ENDIF

       ! Calculate Xch4_a
       XCH4a = 0.0_fp
       DO L = L0,LGOS
          XCH4a = XCH4a + h(L) * prior(L)
       ENDDO

       ! Calculate Xch4_m
       XCH4m = 0.0_fp
       XCH4m = XCH4a
       DO L = L0, LGOS
          XCH4m = XCH4m + ( h(L) * ak(L) * ( GC_CH4(L) - prior(L) ) )
       ENDDO
       GC_XCH4 = 0.0_fp
       GC_XCH4 = XCH4m

       ! Remove mean bias from GEOS-Chem XCH4
       GC_XCH4 = GC_XCH4 - ( 1e-9_fp * MEAN_MODEL_BIAS )

       !--------------------------------------------------------------
       ! Calculate cost function, given S is error in vmr
       ! J = 1/2 [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
       !--------------------------------------------------------------

       ! Calculate difference between modeled and observed profile
       DIFF = GC_XCH4 - GOS_XCH4

       ! Calculate 1/2 * DIFF^T * S_{obs}^{-1} * DIFF
       ! Need to account for the model error (ajt, 03/27/2013)
       FORCE        = ( 1.0_fp / (S_OBS) ) * DIFF
       NEW_COST(NT) = NEW_COST(NT) + 0.5e0 * DIFF * FORCE

       TotalObs = TotalObs + 1

       ! Record information for satellite diagnostics
       IF ( LDCH4SAT ) THEN
          WRITE( 117, 282 ) TotalObs, I, J, GOS(NT)%LON(1), &
                GOS(NT)%LAT(1),GOS(NT)%YEAR(1), &
                GOS(NT)%MONTH(1),GOS(NT)%DAY(1), GOS(NT)%HOUR(1), &
                GOS(NT)%MINUTE(1), GOS(NT)%SEC(1), &
                GET_TAU(), GOS_XCH4, GC_XCH4, S_OBS, 0.5d0*FORCE*DIFF, &
                GOS(NT)%GLINT(1),GOS(NT)%GAIN(1),GOS(NT)%EXP_ID(1)
       ENDIF

    ENDDO  ! NT
    !!$OMP END PARALLEL DO

    ! Number of observations processed in total
    !TotalObs = TotalObs + NOBS

    ! Update cost function
    COST_FUNC = COST_FUNC + SUM(NEW_COST(:))

282 FORMAT( I10,2x,I4,2x,I4,2x,F8.3,2x,F8.4,2x,I4,2x,I2,2x,I2,2x,I2, &
            2x,I2,2x,I2,2x,F12.3,2x,E12.6,2x,E12.6,2x,E12.6,2x,E12.6, &
            2x,I5,2x,I5,2x,A22)

    print*, ' Updated value of COST_FUNC = ', COST_FUNC
    print*, ' GOS contribution           = ', COST_FUNC - OLD_COST
    print*, ' Number of observations this hour = ', NOBS
    print*, ' Number of observations total     = ', TotalObs

  END SUBROUTINE CALC_GOSAT_CH4_FORCE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_intmap
!
! !DESCRIPTION: Function GET\_INTMAP linearly interpolates column quatities
!   based upon the centered (average) pressue levels.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_INTMAP( State_Grid, GCPCEN, GCPSURF, GOSPEDGE, &
                         L0, L1, INTMAP )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid
    REAL(fp),       INTENT(IN)  :: GCPCEN(State_Grid%NZ)
    REAL(fp),       INTENT(IN)  :: GCPSURF
    REAL(fp),       INTENT(IN)  :: GOSPEDGE(MAXLEV)
    INTEGER,        INTENT(IN)  :: L0, L1
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: INTMAP(State_Grid%NZ,MAXLEV)
!
! !REVISION HISTORY:
!  16 Jun 2017 - M. Sulprizio- Initial version based on GOSAT CH4 observation
!                              operator from GC Adjoint v35j
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: LGC, LTM
    REAL(fp) :: DIFF, DELTA_SURFP
    REAL(fp) :: LOW, HI

    !=================================================================
    ! GET_INTMAP begins here!
    !=================================================================

    ! Initialize
    INTMAP(:,:) = 0e+0_fp

    ! Loop over each pressure level of GOS grid
    DO LTM = L0, L1

       ! Find the levels from GC that bracket level LTM
       DO LGC = 1, State_Grid%NZ-1

          LOW = GCPCEN(LGC+1)
          HI  = GCPCEN(LGC)

          ! Match GEOS-Chem level to GOS level
          IF ( GOSPEDGE(LTM) <= HI .and. GOSPEDGE(LTM)  > LOW) THEN

             DIFF             = HI - LOW
             INTMAP(LGC+1,LTM) = ( HI - GOSPEDGE(LTM)  ) / DIFF
             INTMAP(LGC  ,LTM) = ( GOSPEDGE(LTM) - LOW ) / DIFF

          ENDIF

       ENDDO

    ENDDO

    ! Correct for case where GOSAT pressure is higher than the
    ! highest GC pressure center.  In this case, just 1:1 map.
    DO LTM = L0, L1
       IF ( GOSPEDGE(LTM) > GCPCEN(1) ) THEN
          INTMAP(:,LTM) = 0e+0_fp
          INTMAP(1,LTM) = 1e+0_fp
       ENDIF
    ENDDO

  END SUBROUTINE GET_INTMAP
!EOC
END MODULE GOSAT_CH4_MOD
