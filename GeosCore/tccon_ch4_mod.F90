!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tccon_ch4_mod.F90
!
! !DESCRIPTION: Module TCCON\_CH4\_MOD
!\\
!\\
! !INTERFACE:
!
MODULE TCCON_CH4_MOD
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
  PUBLIC  :: CALC_TCCON_CH4_FORCE
!
! !REVISION HISTORY:
!  17 Aug 2017 - M. Sulprizio- Initial version based on TCCON CH4 observation
!                              operator from GC Adjoint v35j with updates from
!                              M. Sulprizio, J.D. Maasakkers, and A. Turner
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  INTEGER,  PARAMETER  :: MAXLEV = 71
  INTEGER,  PARAMETER  :: MAXTCC = 10000
  INTEGER,  PARAMETER  :: NT_FD = 237 !142 ! Index for the FD test
  REAL(fp), PARAMETER  :: GC_XCH4_ERROR = 0e+0_fp !ppb
  REAL(fp), PARAMETER  :: PMAG = 1e0_fp ! Perturbation magnitude (%)
  LOGICAL,  PARAMETER  :: LDCH4SAT    = .TRUE.
  LOGICAL,  PARAMETER  :: FORCE_TCCON = .FALSE. ! Assimilate?
!
! !MODULE VARIABLES
!
  ! Record to store data from each TCCON obs
  TYPE TCCON_CH4_OBS
     INTEGER           :: LTCC(1)
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
     REAL(fp)          :: XH2O(1)
     REAL(fp)          :: XH2OE(1)
     REAL(fp)          :: PSURF(1)
     REAL(fp)          :: MH2O(1)
     REAL(fp)          :: MAIR(1)
     REAL(fp)          :: AVNUM(1)
     REAL(fp)          :: PRES(MAXLEV)
     REAL(fp)          :: PRIOR(MAXLEV)
     REAL(fp)          :: AVG_KERNEL(MAXLEV)
     REAL(fp)          :: GRAVITY(MAXLEV)
     REAL(fp)          :: H2O(MAXLEV)
     CHARACTER(LEN=2)  :: SITE(MAXLEV)
     INTEGER           :: QFLAG(1)
  ENDTYPE TCCON_CH4_OBS

  TYPE(TCCON_CH4_OBS)             :: TCC(MAXTCC)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_tcc_ch4_obs
!
! !DESCRIPTION: Subroutine READ\_TCC\_CH4\_OBS reads the file and passes back
!  info contained therein. (ajt, 01/13/14)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_TCC_CH4_OBS( YYYYMMDD, NTCC )
!
! ! USES:
!
    USE TIME_MOD,               ONLY : EXPAND_DATE
    USE ERROR_MOD,              ONLY : ALLOC_ERR
    USE ERROR_MOD,              ONLY : GEOS_CHEM_STOP
    USE TIME_MOD,               ONLY : GET_YEAR
!
! !INPUT PARAMETERS:
!
    INTEGER,            INTENT(IN)  :: YYYYMMDD ! Current date
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(OUT) :: NTCC     ! Number of TCC retrievals
!
! !REVISION HISTORY:
!  12 Jan 2014 - A. Turner - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                       :: LTCC
    INTEGER                       :: N, J
    INTEGER                       :: MYEAR
    CHARACTER(LEN=4)              :: CYEAR
    LOGICAL                       :: EXIST_VAR
    REAL*4                        :: temppres(MAXLEV)
    REAL*4                        :: tempprior(MAXLEV)
    REAL*4                        :: tempak(MAXLEV)
    REAL*4, ALLOCATABLE           :: pres(:,:)
    REAL*4, ALLOCATABLE           :: prior(:,:)
    REAL*4, ALLOCATABLE           :: ak(:,:)
    REAL*4, ALLOCATABLE           :: pres_w(:,:)
    REAL*4, ALLOCATABLE           :: lat(:)
    REAL*4, ALLOCATABLE           :: lon(:)
    INTEGER, ALLOCATABLE          :: year(:)
    INTEGER, ALLOCATABLE          :: month(:)
    INTEGER, ALLOCATABLE          :: day(:)
    INTEGER, ALLOCATABLE          :: hour(:)
    INTEGER, ALLOCATABLE          :: minute(:)
    INTEGER, ALLOCATABLE          :: second(:)
    REAL*4, ALLOCATABLE           :: xch4(:)
    REAL*4, ALLOCATABLE           :: xch4_error(:)
    REAL*4, ALLOCATABLE           :: xh2o(:)
    REAL*4, ALLOCATABLE           :: xh2o_error(:)
    REAL*4, ALLOCATABLE           :: mass_h2o(:)
    REAL*4, ALLOCATABLE           :: mass_dry_air(:)
    REAL*4, ALLOCATABLE           :: AvNum(:)
    REAL*4, ALLOCATABLE           :: grav_apri(:,:)
    REAL*4, ALLOCATABLE           :: tccon_h2o(:,:)
    INTEGER, ALLOCATABLE          :: qflag(:)
    CHARACTER(LEN=2), ALLOCATABLE :: TCCON_site(:)

    INTEGER                       :: NLEV, StrLen
    INTEGER                       :: I, L, AS

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
    ! READ_TCC_CH4_OBS begins here!
    !=================================================================

    caller = 'READ_TCC_CH4_OBS in tccon_ch4_mod.F90'

    ! Get current year
    MYEAR = GET_YEAR()
    WRITE( CYEAR, '(i4)' ) MYEAR

    ! Filename
    nc_file = 'TCCON_dat_YYYYMMDD.ncdf'
    CALL Expand_Date( nc_file, YYYYMMDD, 9999 )

    ! Construct complete filename
    dir = '/n/regal/jacob_lab/msulprizio/ExtData/Obs/TCCON/' // CYEAR // '/'
    nc_file = TRIM( dir ) // TRIM( nc_file )
    WRITE( 6, 10 ) TRIM( nc_file )
10  FORMAT( '     - Reading ', a)

    ! Make sure the file exists (ajt, 03/31/2013)
    INQUIRE( FILE=nc_file, EXIST=EXIST_VAR )
    IF ( .not. EXIST_VAR ) THEN
       NTCC = -1d0
       RETURN
    ENDIF

    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'numexp',   NTCC )
    CALL NcGet_DimLen( fId, 'levels',   NLEV )
    CALL NcGet_DimLen( fId, 'string_length', StrLen )

    IF ( NLEV > MAXLEV ) THEN
       print*,' # Levels this day = ', NLEV
       print*, 'WARNING: NLEV > MAXLEV. Need to increase'
       print*, ' MAXLEV in tccon_ch4_mod.F90.'
       CALL GEOS_CHEM_STOP
    ENDIF

    print*,' # TCCON Observations this day = ', NTCC
    print*, 'levels', NLEV

    !--------------------------------
    ! Read 1-D Data
    !--------------------------------

    ALLOCATE( lat(             NTCC ), STAT=AS )
    ALLOCATE( lon(             NTCC ), STAT=AS )
    ALLOCATE( year(            NTCC ), STAT=AS )
    ALLOCATE( month(           NTCC ), STAT=AS )
    ALLOCATE( day(             NTCC ), STAT=AS )
    ALLOCATE( hour(            NTCC ), STAT=AS )
    ALLOCATE( minute(          NTCC ), STAT=AS )
    ALLOCATE( second(          NTCC ), STAT=AS )
    ALLOCATE( qflag(           NTCC ), STAT=AS )
    ALLOCATE( xch4(            NTCC ), STAT=AS )
    ALLOCATE( xch4_error(      NTCC ), STAT=AS )
    ALLOCATE( xh2o(            NTCC ), STAT=AS )
    ALLOCATE( xh2o_error(      NTCC ), STAT=AS )
    ALLOCATE( mass_h2o(        NTCC ), STAT=AS )
    ALLOCATE( mass_dry_air(    NTCC ), STAT=AS )
    ALLOCATE( AvNum(           NTCC ), STAT=AS )
    ALLOCATE( TCCON_site(      NTCC ), STAT=AS )
    ALLOCATE( prior(     NLEV, NTCC ), STAT=AS )
    ALLOCATE( pres(      NLEV, NTCC ), STAT=AS )
    ALLOCATE( ak(        NLEV, NTCC ), STAT=AS )
    ALLOCATE( pres_w(    NLEV, NTCC ), STAT=AS )
    ALLOCATE( grav_apri( NLEV, NTCC ), STAT=AS )
    ALLOCATE( tccon_h2o( NLEV, NTCC ), STAT=AS )

    ! Latitude
    v_name = 'lat'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( lat, fId, TRIM(v_name), st1d, ct1d )

    ! Longitude
    v_name = 'lon'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( lon, fId, TRIM(v_name), st1d, ct1d )

    ! Year
    v_name = 'year'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( year, fId, TRIM(v_name), st1d, ct1d )

    ! Month
    v_name = 'month'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( month, fId, TRIM(v_name), st1d, ct1d )

    ! Day
    v_name = 'day'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( day, fId, TRIM(v_name), st1d, ct1d )

    ! Hour
    v_name = 'hour'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( hour, fId, TRIM(v_name), st1d, ct1d )

    ! Minute
    v_name = 'minute'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( minute, fId, TRIM(v_name), st1d, ct1d )

    ! Second
    v_name = 'second'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( second, fId, TRIM(v_name), st1d, ct1d )

    ! Qflag (0=good)
    v_name = 'quality_flag'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( qflag, fId, TRIM(v_name), st1d, ct1d )

    ! XCH4 (ppm)
    v_name = 'xch4_ppm'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( xch4, fId, TRIM(v_name), st1d, ct1d )

    ! XCH4_Error (ppm)
    v_name = 'xch4_ppm_error'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( xch4_error, fId, TRIM(v_name), st1d, ct1d )

    ! XH2O (ppm)
    v_name = 'xh2o_ppm'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( xh2o, fId, TRIM(v_name), st1d, ct1d )

    ! XH2O_Error (ppm)
    v_name = 'xh2o_ppm_error'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( xh2o_error, fId, TRIM(v_name), st1d, ct1d )

    ! Mass of water (kg/mol)
    v_name = 'mass_h2o'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( mass_h2o, fId, TRIM(v_name), st1d, ct1d )

    ! Mass of dry air (kg/mol)
    v_name = 'mass_dry_air'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( mass_dry_air, fId, TRIM(v_name), st1d, ct1d )

    ! Avogadro's number (molecules/mol)
    v_name = 'avogadros_number'
    st1d   = (/ 1    /)
    ct1d   = (/ NTCC /)
    CALL NcRd( AvNum, fId, TRIM(v_name), st1d, ct1d )

    ! Site
    v_name = 'site'
    st2d   = (/ 1, 1    /)
    ct2d   = (/ StrLen, NTCC /)
    CALL NcRd( TCCON_site, fId, TRIM(v_name), st2d, ct2d )

    !--------------------------------
    ! Read 2-D Data
    !--------------------------------

    ! APRIORI (ppb)
    v_name = 'ch4_apriori'
    st2d   = (/ 1,    1    /)
    ct2d   = (/ NLEV, NTCC /)
    CALL NcRd( prior, fId, TRIM(v_name), st2d, ct2d )

    ! Pressure (hPa)
    v_name = 'pressure_apriori'
    st2d   = (/ 1,    1    /)
    ct2d   = (/ NLEV, NTCC /)
    CALL NcRd( pres, fId, TRIM(v_name), st2d, ct2d )

    ! Averaging Kernel (unitless)
    v_name = 'ch4_ak'
    st2d   = (/ 1,    1    /)
    ct2d   = (/ NLEV, NTCC /)
    CALL NcRd( ak, fId, TRIM(v_name), st2d, ct2d )

    ! Pressure weights (unitless)
    v_name = 'pressure_ak'
    st2d   = (/ 1,    1    /)
    ct2d   = (/ NLEV, NTCC /)
    CALL NcRd( pres_w, fId, TRIM(v_name), st2d, ct2d )

    ! Gravity from the apriori (m/s2)
    v_name = 'grav_apriori'
    st2d   = (/ 1,    1    /)
    ct2d   = (/ NLEV, NTCC /)
    CALL NcRd( grav_apri, fId, TRIM(v_name), st2d, ct2d )

    ! H2O from the apriori (m/s2)
    v_name = 'h2o_apriori'
    st2d   = (/ 1,    1    /)
    ct2d   = (/ NLEV, NTCC /)
    CALL NcRd( tccon_h2o, fId, TRIM(v_name), st2d, ct2d )

    !--------------------------------
    ! Place data into TCC structure
    !--------------------------------
    DO N = 1, NTCC

       ! 0-D data
       TCC(N)%LTCC(1)      = NLEV
       TCC(N)%LAT(1)       = lat(N)
       TCC(N)%LON(1)       = lon(N)
       TCC(N)%YEAR(1)      = year(N)
       TCC(N)%MONTH(1)     = month(N)
       TCC(N)%DAY(1)       = day(N)
       TCC(N)%HOUR(1)      = hour(N)
       TCC(N)%MINUTE(1)    = minute(N)
       TCC(N)%SEC(1)       = second(N)
       TCC(N)%CH4(1)       = xch4(N) * 1d3        ! ppm --> ppb
       TCC(N)%CH4_ERROR(1) = xch4_error(N) * 1d3  ! ppm --> ppb
       TCC(N)%XH2O(1)      = xh2o(N) * 1d-6       ! ppm --> v/v
       TCC(N)%XH2OE(1)     = xh2o_error(N) * 1d-6 ! ppm --> v/v
       TCC(N)%PSURF(1)     = pres(1,N)
       TCC(N)%QFLAG(1)     = qflag(N)
       TCC(N)%AVNUM(1)     = AvNum(N)                   ! molecules/mole
       TCC(N)%MH2O(1)      = mass_h2o(N) / AvNum(N)     ! kg/molecule
       TCC(N)%MAIR(1)      = mass_dry_air(N) / AvNum(N) ! kg/molecule
       TCC(N)%SITE(1)      = TCCON_site(N)

       ! Time fraction of day
       TCC(N)%TIME(1) = ( REAL(HOUR(N))*3600. + REAL(MINUTE(N))*60. &
                        + REAL(SECOND(N)) ) / 86400d0

       ! 1-D data
       LTCC = NLEV
       TCC(N)%PRES(1:LTCC)       = pres(1:LTCC,N)
       TCC(N)%PRIOR(1:LTCC)      = prior(1:LTCC,N) * 1d-9 ! ppb --> v/v
       TCC(N)%AVG_KERNEL(1:LTCC) = ak(1:LTCC,N)
       TCC(N)%GRAVITY(1:LTCC)    = grav_apri(1:LTCC,N)
       TCC(N)%H2O(1:LTCC)        = tccon_h2o(1:LTCC,N)

       !!  Reverse indices so that L=1 is surface
       !LTCC = NLEV
       !temppres(:)  = 0d0
       !tempprior(:) = 0d0
       !tempak(:)    = 0d0
       !DO L = 1, LTCC
       !   tempak(L)    = ak(LTCC+1-L,N)
       !   tempprior(L) = prior(LTCC+1-L,N)
       !   temppres(L)  = pres(LTCC+1-L,N) * 1d-2  ! [Pa] --> [hPa]
       !ENDDO
       !
       !TCC(N)%PRES(1:LTCC)       = temppres(1:LTCC)
       !TCC(N)%PRIOR(1:LTCC)      = tempprior(1:LTCC)
       !TCC(N)%AVG_KERNEL(1:LTCC) = tempak(1:LTCC)

       ! Store the surface pressure
       TCC(N)%PSURF(1) = TCC(N)%PRES(1)

    ENDDO

    !--------------------------------
    ! Close netCDF file
    !--------------------------------

    ! Echo info
    WRITE( 6, 20 ) YYYYMMDD
20  FORMAT( '     - Finished reading TCCON CH4 observations for ', i8)

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE READ_TCC_CH4_OBS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_tccon_ch4_force
!
! !DESCRIPTION: Subroutine CALC\_TCCON\_CH4\_FORCE calculates the adjoint
!  forcing from the TCCON CH4 observations and updates the cost function.
!  (dkh, 10/12/10)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_TCCON_CH4_FORCE( Input_Opt, State_Chm, State_Grid, &
                                   State_Met )
!
! !USES:
!
    USE ERROR_MOD,          ONLY : IT_IS_NAN
    USE ERROR_MOD,          ONLY : IT_IS_FINITE
    USE GC_GRID_MOD,        ONLY : GET_IJ
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
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
!  17 Aug 2017 - M. Sulprizio- Initial version based on TCCON CH4 observation
!                              operator from GC Adjoint v35j with updates from
!                              M. Sulprizio, J.D. Maasakkers, and A. Turner
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp), SAVE     :: COST_FUNC   ! Cost function [unitless]
    INTEGER            :: id_CH4
    INTEGER            :: NTSTART, NTSTOP, NT, YYYYMMDD
    INTEGER            :: IIJJ(2), I,      J
    INTEGER            :: L,       LL,     LTCC
    INTEGER            :: JLOOP,   NOBS,   IND
    INTEGER            :: INDS(MAXTCC)
    REAL(fp)           :: GC_PRES(State_Grid%NZ)
    REAL(fp)           :: GC_PEDGE(State_Grid%NZ+1)
    REAL(fp)           :: GC_CH4_NATIVE(State_Grid%NZ)
    REAL(fp)           :: GC_H2O_NATIVE(State_Grid%NZ)
    REAL(fp)           :: GC_CH4(MAXLEV)
    REAL(fp)           :: GC_H2O(MAXLEV)
    REAL(fp)           :: GC_CH4_cm2(MAXLEV)
    REAL(fp)           :: GC_PSURF
    REAL(fp)           :: MAP(State_Grid%NZ,MAXLEV)
    REAL(fp)           :: CH4_HAT(MAXLEV)
    REAL(fp)           :: NEW_COST(MAXTCC)
    REAL(fp)           :: OLD_COST
    REAL(fp), SAVE     :: TIME_FRAC(MAXTCC)
    INTEGER,SAVE       :: NTCC
    REAL(fp)           :: frac, RHS, LHS
    REAL(fp)           :: CH4_PRIOR(MAXLEV)
    REAL(fp)           :: CH4_PRIOR_cm2(MAXLEV)
    REAL(fp)           :: molecongos(MAXLEV)
    REAL(fp)           :: TCC_XCH4, GC_XCH4, TCC_XCH4_ERROR
    REAL(fp)           :: TCC_P_EDGE(MAXLEV)
    REAL(fp)           :: p(MAXLEV), h(MAXLEV), psurf
    REAL(fp)           :: ak(MAXLEV), prior(MAXLEV)
    REAL(fp)           :: XCH4m, XCH4a
    REAL(fp)           :: SATELLITE_BIAS(3) ! Hardcode for now (mps,6/16/17)
    REAL(fp)           :: MEAN_MODEL_BIAS   ! Hardcode for now (mps,6/16/17)
    REAL(fp)           :: FORCE
    REAL(fp)           :: DIFF
    REAL(fp)           :: S_OBS

    ! For miscellaneous
    LOGICAL, SAVE      :: FIRST = .TRUE.
    LOGICAL            :: WRITE_FILE = .TRUE.
    INTEGER            :: IOS
    INTEGER, SAVE      :: TotalObs = 0
    CHARACTER(LEN=255) :: FILENAME

    !=================================================================
    ! CALC_TCCON_CH4_FORCE begins here!
    !=================================================================

    print*, '     - CALC_TCCON_CH4_FORCE '

    ! Initialize species ID flag
    id_CH4     = Ind_('CH4'       )

    ! Reset
    NEW_COST(:) = 0.0_fp

    ! Hardcode bias values for now (mps, 6/16/17)
    SATELLITE_BIAS(1) = 0.0e+0_fp
    SATELLITE_BIAS(2) = 0.0e+0_fp
    SATELLITE_BIAS(3) = 0.0e+0_fp
    MEAN_MODEL_BIAS   = 0.0e+0_fp

    ! Open files for diagnostic output
    IF ( FIRST ) THEN

       ! For recording model and observation values
       FILENAME = 'tccon_obs.00.m'
       FILENAME = TRIM( Input_Opt%RUN_DIR ) //  TRIM( FILENAME )

       ! Check if file exists
       INQUIRE( FILE=TRIM(FILENAME), EXIST=WRITE_FILE )

       IF ( WRITE_FILE ) THEN
          OPEN( 119,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN', &
                IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
       ELSE
          OPEN( 119,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN', &
                IOSTAT=IOS, FORM='FORMATTED',    ACCESS='APPEND' )

          ! Write header of tccon_obs.NN.m
          WRITE( 119, 281 ) '       NNN', &
                            '   I', '   J', '     LON','     LAT','YYYY', &
                            'MM', 'DD', 'hh', 'mm', 'ss', 'site', &
                            '         TAU', '       TCCON', &
                            '       model', '       S_OBS', &
                            '    COST_FUN', 'XH2O', 'XH2O_ERROR'
       ENDIF

281    FORMAT( A10,2x,A4,2x,A4,2x,A8,2x,A8,2x,A4,2x,A2,2x,A2,2x,A2,2x, &
               A2,2x,A2,2x,A4,2x,A12,2x,A36,2x,A36,2x,A36,2x,A36,2x, &
               A36,2x,A36 )

       ! Set Total Observations = 0
       TotalObs = 0

    ENDIF

    ! Save a value of the cost function first
    OLD_COST = COST_FUNC

    ! Read Observations at first call and at end of the day
    IF ( FIRST .OR. ITS_A_NEW_DAY() ) THEN

       ! Read the TCC CH4 file for this month
       YYYYMMDD = 1d4*GET_YEAR() + 1d2*GET_MONTH() + GET_DAY()
       CALL READ_TCC_CH4_OBS( YYYYMMDD, NTCC )

       IF ( FIRST ) FIRST = .FALSE.

       ! Make sure there are observations on this day
       IF ( NTCC .EQ. -1d0 ) RETURN

    ENDIF

    ! Get indices of TCC observations in the current hour
    !   At the start of each hour, assimilate observations that
    !   were made in the previous 60 minutes.
    !   For example, at time 18:00, assimilate observations
    !   made from 18:00 - 18:59
    INDS(:) = 0d0
    NOBS    = 0
    DO NT=1,NTCC
       IF ( TCC(NT)%MONTH(1) .EQ. GET_MONTH() .AND. &
            TCC(NT)%DAY(1)   .EQ. GET_DAY()   .AND. &
            TCC(NT)%HOUR(1)  .EQ. GET_HOUR()  ) THEN
          NOBS = NOBS + 1
          INDS(NOBS) = NT
          !print*,'Found a good observation! NT = ', NT
       ENDIF
    ENDDO

    IF ( NOBS == 0 ) THEN
       print*, ' No matching TCCON CH4 obs for this hour'
       RETURN
    ENDIF

    print*, ' for hour range: ', GET_HOUR(), GET_HOUR()+1
    print*, ' found # observations: ', NOBS

    !! need to update this in order to do i/o with this loop parallel
    !!      ! Now do a parallel loop for analyzing data
    !!!$OMP PARALLEL DO        &
    !!!$OMP DEFAULT( PRIVATE ) &
    !!!$OMP PRIVATE( IND, NT, MAP, LTCC, IIJJ,  I, J,  L,   LL, JLOOP ) &
    !!!$OMP PRIVATE( GC_CH4, FORCE, CH4_PRIOR, GC_PRES, FILENAME      ) &
    !!!$OMP PRIVATE( GC_PEDGE, GC_PSURF, GC_CH4_NATIVE, TCC_XCH4      ) &
    !!!$OMP PRIVATE( GC_H2O, GC_H2O_NATIVE,                           ) &
    !!!$OMP PRIVATE( TCC_XCH4_ERROR, S_OBS, h, psurf, p, XCH4a, XCH4m ) &
    !!!$OMP PRIVATE( GC_XCH4, DIFF, DIFF_ADJ, GC_XCH4_ADJ             ) &
    !!!$OMP PRIVATE( GC_CH4_NATIVE_ADJ, GC_CH4_ADJ, TotalObs          )
    DO IND = 1, NOBS

       NT = INDS(IND)
       !print*, '     - CALC_TCCON_CH4_FORCE: analyzing record ', NT

       ! quality screening (0 is good!)
       IF ( TCC(NT)%QFLAG(1) .NE. 0 ) THEN
          print*, ' BAD QFLAG, skipping record        ', NT
          CYCLE
       ENDIF

       ! Check for NaN in data
       IF ( IT_IS_NAN( TCC(NT)%CH4(1) ) ) THEN
          print*, ' XCH4 is NaN, skipping record      ', NT
          CYCLE
       ENDIF

       ! Check for infinity in data
       IF ( .not. IT_IS_FINITE( TCC(NT)%CH4(1) ) ) THEN
          print*, ' XCH4 is infinity, skipping record ', NT
          CYCLE
       ENDIF

       ! Check for negative/zero data
       IF ( TCC(NT)%CH4(1) <= 0d0 ) THEN
          print*, ' XCH4 is <= 0, skipping record     ', NT
          CYCLE
       ENDIF

       ! Skip observations outside the nested domain
       IF ( TCC(NT)%LAT(1) < State_Grid%YMin .OR. &
            TCC(NT)%LAT(1) > State_Grid%YMax .OR. &
            TCC(NT)%LON(1) < State_Grid%XMin .OR. &
            TCC(NT)%LON(1) > State_Grid%XMax ) THEN
          print*, ' Outside nested domain, skipping record ', NT
          CYCLE
       ENDIF

       ! Get grid box of current record
       IIJJ  = GET_IJ( REAL(TCC(NT)%LON(1),4), REAL(TCC(NT)%LAT(1),4), &
                       State_Grid )
       I     = IIJJ(1)
       J     = IIJJ(2)

       !! skip observations where the TCCON surface pressure is
       !! less than the model.  Use w/ updated interp (ajt, 11/6/13)
       !IF ( (TCC(NT)%PSURF(1) - GET_PEDGE(I,J,1)) .GT. 50d0 ) THEN
       !   print*, ' Psurf threshold not met, skipping record ', NT
       !   CYCLE
       !ENDIF

       ! begin good observations
       print*,'Begin assimilating good observation. NT = ', NT

       ! For safety, initialize these up to LTCC
       LTCC            = 0
       GC_CH4(:)       = 0.0_fp
       MAP(:,:)        = 0.0_fp
       FORCE           = 0.0_fp
       CH4_PRIOR(:)    = 0.0_fp

       ! Copy variable names to make coding a bit cleaner
       LTCC = TCC(NT)%LTCC(1)
       DO L = 1, LTCC
          CH4_PRIOR(L) = TCC(NT)%PRIOR(L)
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
       CALL GET_INTMAP( State_Grid, GC_PRES, TCC(NT)%PRES, MAP )

       ! Get CH4 values at native model resolution
       GC_CH4_NATIVE(:) = 0.0_fp
       GC_H2O_NATIVE(:) = 0.0_fp

       ! Get species concentrations
       ! Convert from [kg/box] --> [v/v]
       GC_CH4_NATIVE(:) = State_Chm%Species(I,J,:,id_CH4) * ( AIRMW / &
                          State_Chm%SpcData(id_CH4)%Info%emMW_g )

       GC_H2O_NATIVE(:) = State_Met%AVGW(I,J,:)

       !! Interpolate GC CH4 column to TCCON grid
       !! use the method of Parker et al., (2011)
       !GC_CH4 = INTERP_GC_TCC( LTCC, State_Grid%NZ, GC_PRES, GC_PEDGE, &
       !            TCC(NT)%PRES(1:LTCC), GC_CH4_NATIVE )

       ! Old interpolation method from kjw
       DO LL = 1, LTCC
          GC_CH4(LL) = 0.0_fp
          GC_H2O(LL) = 0.0_fp
          DO L = 1, State_Grid%NZ
             GC_CH4(LL) = GC_CH4(LL) + MAP(L,LL) * GC_CH4_NATIVE(L)
             GC_H2O(LL) = GC_H2O(LL) + MAP(L,LL) * GC_H2O_NATIVE(L)
             !GC_H2O(LL) = TCC(NT)%H2O(LL)
          ENDDO
       ENDDO

       ! Store the TCCON XCH4 proxy in [v/v]
       TCC_XCH4       = TCC(NT)%CH4(1)       * 1e-9_fp
       TCC_XCH4_ERROR = TCC(NT)%CH4_ERROR(1) * 1e-9_fp

       ! Remove any TCCON bias (probably not...)
       !TCC_XCH4 = TCC_XCH4 + 1d-9 * ( SATELLITE_BIAS(1) &
       !                    + SATELLITE_BIAS(2)*(TCC(NT)%LAT(1)) &
       !                    + SATELLITE_BIAS(3)*(TCC(NT)%LAT(1))**2 )

       ! Get the S_obs, assume stddev adds in quadrature, variance
       ! adds linearly.  (ajt, 03/27/2013)
       S_OBS = TCC_XCH4_ERROR**2 + (GC_XCH4_ERROR * 1e-9_fp)**2

       !--------------------------------------------------------------
       ! Apply TCCON observation operator
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
       h(:)          = 0.0_fp
       psurf         = TCC(NT)%PSURF(1)
       p(1:LTCC)     = TCC(NT)%PRES(1:LTCC)
       ak(1:LTCC)    = TCC(NT)%AVG_KERNEL(1:LTCC)
       prior(1:LTCC) = TCC(NT)%PRIOR(1:LTCC)

       ! Need to integrate from the toa to surface (ajt, 05/21/13)
       IF (LTCC .GT. 1) THEN
          IF(TCC(NT)%PRES(2) .LT. TCC(NT)%PRES(1)) THEN
             p(1:LTCC) = p(LTCC:1:-1)
          ENDIF
       ENDIF

       L = 1
       h(L) = 1./TCC(NT)%PSURF(1) * ABS( &
              ( -1e0*p(L) + ( (p(L+1)-p(L))/(LOG(p(L+1)/p(L))) ) ) )
       L = LTCC
       h(L) = 1./TCC(NT)%PSURF(1) * ABS( &
              (  p(L) - ( (p(L)-p(L-1))/(LOG(p(L)/p(L-1))) ) ) )
       DO L=2,LTCC-1
          h(L) = 1./TCC(NT)%PSURF(1) * ABS( &
                 ( -1e0*p(L) + ( (p(L+1)-p(L))/(LOG(p(L+1)/p(L))) ) ) + &
                 (      p(L) - ( (p(L)-p(L-1))/(LOG(p(L)/p(L-1))) ) )   )
       ENDDO

       ! Now return to the orientation of the other variables
       IF (LTCC .GT. 1) THEN
          IF(TCC(NT)%PRES(2) .LT. TCC(NT)%PRES(1)) THEN
             h(1:LTCC) = h(LTCC:1:-1)
             p(1:LTCC) = p(LTCC:1:-1)
          ENDIF
       ENDIF

       ! Get the TCCON pressure edges
       TCC_P_EDGE = TCC_PEDGE( LTCC, TCC(NT)%PRES(1:LTCC) )

       ! Compute h and convert the prior to a dry-air mole fraction
       DO L = 1, LTCC
          !h(L) = ( ( TCC_P_EDGE(L) - TCC_P_EDGE(L+1) )      &
          !     / ( TCC(NT)%GRAV(L) * TCC(NT)%MAIR(1)        &
          !     * (1e0 + GC_H2O(L) / ( 1e0 - GC_H2O(L) )     &
          !     * (TCC(NT)%MH2O(1) / TCC(NT%MAIR(1)) ) ) ) ) &
          !     / ( psurf                                    &
          !     / ( TCC(NT)%GRAV(L) * TCC(NT)%MAIR(1)        &
          !     * (1e0 + GC_H2O(L) / ( 1e0 - GC_H2O(L) )     &
          !     * (TCC(NT)%MH2O(1) / TCC(NT%MAIR(1)) ) ) ) )
          !h(L) = ( TCC_P_EDGE(L) - TCC_P_EDGE(L+1) )        &
          !     / ( TCC(NT)%GRAV(L) * TCC(NT)%MAIR(1)        &
          !     * (1e0 + GC_H2O(L) / ( 1e0 - GC_H2O(L) )     &
          !     * (TCC(NT)%MH2O(1) / TCC(NT%MAIR(1)) ) ) )

          prior(L) = prior(L) / ( 1.0_fp - GC_H2O(L) )
       ENDDO

       ! Calculate Xch4_a
       XCH4a  = 0.0_fp
       DO L = 1,LTCC
          XCH4a  = XCH4a + h(L) * prior(L)
       ENDDO

       ! Calculate Xch4_m
       XCH4m = 0.0_fp
       XCH4m = XCH4a
       DO L = 1, LTCC
          XCH4m = XCH4m + ( h(L) * ak(L) * ( GC_CH4(L) - prior(L) ) )
       ENDDO
       GC_XCH4 = 0.0_fp
       GC_XCH4 = XCH4m

       ! Remove any GEOS-Chem bias
       GC_XCH4 = GC_XCH4 + ( 1e-9_fp * MEAN_MODEL_BIAS )

       !--------------------------------------------------------------
       ! Calculate cost function, given S is error in vmr
       ! J = 1/2 [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
       !--------------------------------------------------------------

       ! Calculate difference between modeled and observed profile
       DIFF = GC_XCH4 - TCC_XCH4

       ! Calculate 1/2 * DIFF^T * S_{obs}^{-1} * DIFF
       ! Need to account for the model error (ajt, 03/27/2013)
       FORCE        = ( 1.0_fp / (S_OBS) ) * DIFF
       NEW_COST(NT) = NEW_COST(NT) + 0.5e0 * DIFF * FORCE

       TotalObs = TotalObs + 1

       ! Record information for satellite diagnostics
       IF ( LDCH4SAT ) THEN
          WRITE( 119, 282 ) TotalObs, I, J, TCC(NT)%LON(1),           &
               TCC(NT)%LAT(1),TCC(NT)%YEAR(1),                        &
               TCC(NT)%MONTH(1),TCC(NT)%DAY(1), TCC(NT)%HOUR(1),      & 
               TCC(NT)%MINUTE(1), TCC(NT)%SEC(1), TCC(NT)%SITE(1),    &
               GET_TAU(), TCC_XCH4, GC_XCH4, S_OBS, 0.5d0*FORCE*DIFF, &
               TCC(NT)%XH2O(1), TCC(NT)%XH2OE(1)
       ENDIF

    ENDDO  ! NT
    !!!$OMP END PARALLEL DO

    ! Number of observations processed in total
    !TotalObs = TotalObs + NOBS

    ! Update cost function
    IF ( FORCE_TCCON ) THEN
       COST_FUNC = COST_FUNC + SUM(NEW_COST(:))
    ELSE
       COST_FUNC = OLD_COST
    ENDIF

282 FORMAT( I10,2x,I4,2x,I4,2x,F8.3,2x,F8.4,2x,I4,2x,I2,2x,I2,2x,I2,  &
            2x,I2,2x,I2,2x,A4,2x,F12.3,2x,E36.30,2x,E36.30,2x,E36.30, &
            2x,E36.30,2x,E36.30,2x,E36.30)

    print*, ' Updated value of COST_FUNC = ', COST_FUNC
    print*, ' TCC contribution           = ', COST_FUNC - OLD_COST
    print*, ' Number of observations this hour = ', NOBS
    print*, ' Number of observations total     = ', TotalObs

  END SUBROUTINE CALC_TCCON_CH4_FORCE
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
  SUBROUTINE GET_INTMAP( State_Grid, GCPCEN, TCCPCEN, HINTERPZ )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid
    REAL(fp),       INTENT(IN)  :: GCPCEN(State_Grid%NZ)
    REAL(fp),       INTENT(IN)  :: TCCPCEN(MAXLEV)
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT) :: HINTERPZ(State_Grid%NZ,MAXLEV)
!
! !REVISION HISTORY:
!  17 Aug 2017 - M. Sulprizio- Initial version based on TCCON CH4 observation
!                              operator from GC Adjoint v35j
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: LGC, LTM, NGC, NSAT
    REAL(fp) :: DIFF, HI, MID, LOW
    LOGICAL  :: WHICH_MAP = .TRUE.

    !=================================================================
    ! GET_HINTERPZ_2 begins here!
    !=================================================================

    HINTERPZ(:,:) = 0e+0_fp

    !! Determine what to map onto
    !IF ( MAXLEV > State_Grid%NZ ) THEN
    !   WHICH_MAP = .FALSE.
    !ENDIF

    ! Determine the looping
    IF ( WHICH_MAP ) THEN
       NGC  = State_Grid%NZ-1
       NSAT = MAXLEV
    ELSE
       NSAT = MAXLEV-1
       NGC  = State_Grid%NZ
    ENDIF

    ! Loop over each pressure level of TCCON grid
    DO LTM = 1, NSAT

       ! Find the levels from TCCON that bracket GC
       DO LGC = 1, NGC

          ! Find the bounding values
          IF (MAXLEV > State_Grid%NZ) THEN
             LOW = GCPCEN(LGC+1)
             HI  = GCPCEN(LGC)
             MID = TCCPCEN(LTM)
          ELSE
             LOW = TCCPCEN(LTM+1)
             HI  = TCCPCEN(LTM)
             MID = GCPCEN(LGC)
          ENDIF

          ! Match TCCON level to the GEOS-Chem level
          IF ( ( MID <= HI ) .AND. (  MID  > LOW ) ) THEN
             DIFF = HI - LOW
             IF (MAXLEV > State_Grid%NZ) THEN
                HINTERPZ(LGC+1,LTM) = ( HI - MID  ) / DIFF
                HINTERPZ(LGC  ,LTM) = ( MID - LOW ) / DIFF
             ELSE
                HINTERPZ(LGC,LTM+1) = ( HI - MID  ) / DIFF
                HINTERPZ(LGC,LTM  ) = ( MID - LOW ) / DIFF
             ENDIF
          ENDIF

       ENDDO
    ENDDO

    ! Correct for case where TCCON pressure is higher than the
    ! highest GC pressure center.  In this case, just 1:1 map.
    DO LTM = 1, MAXLEV
       IF ( TCCPCEN(LTM) > GCPCEN(1) ) THEN
          HINTERPZ(:,LTM) = 0e+0_fp
          HINTERPZ(1,LTM) = 1e+0_fp
       ENDIF
    ENDDO

  END SUBROUTINE GET_INTMAP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tcc_pedge
!
! !DESCRIPTION: Function TCC\_PEDGE get the pressure edges for the TCCON
!  pressure grid
!\\
!\\
! !INTERFACE:
!
  FUNCTION TCC_PEDGE( LTCC, TCC_P_I )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: LTCC              ! # TCCON presure levels
    REAL(fp), INTENT(IN) :: TCC_P_I(LTCC)     ! TCCON pressure levels [mbar]
!
! !RETURN VALUE:
!
    REAL(fp)             :: TCC_PEDGE(LTCC+1) ! TCCON pressure edges
!
! !REVISION HISTORY:
!  17 Aug 2017 - M. Sulprizio- Initial version based on TCCON CH4 observation
!                              operator from GC Adjoint v35j
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: LL
    REAL(fp)             :: p(LTCC)
    REAL(fp)             :: PEDGE(LTCC+1)
    REAL(fp)             :: LAYER_BOT(LTCC)
    REAL(fp)             :: LAYER_TOP(LTCC)
    REAL(fp)             :: TCC_P(LTCC)
    REAL(fp)             :: test1, test2

    !=================================================================
    ! TCC_PEDGE begins here!
    !=================================================================

    ! Reverse the array orders so toa is first
    TCC_P = TCC_P_I
    IF (LTCC .GT. 1) THEN
       IF (TCC_P_I(2) .LT. TCC_P_I(1)) TCC_P = TCC_P(LTCC:1:-1)
    ENDIF

    ! Determine how much a level influences the levels above and below it
    DO LL = 1, LTCC

       ! Make coding cleaner
       p = TCC_P

       ! Different cases: (a) toa, (b) surface, (c) else
       IF (LL .EQ. 1) THEN
          test1 = ABS( -1e0*p(LL) + ( ( p(LL+1) - p(LL) ) &
                  / LOG( p(LL+1) / p(LL) ) ) )
          test2 = 0e+0_fp
       ELSE IF (LL .EQ. LTCC) THEN
          test1 = 0e+0_fp
          test2 = ABS( p(LL) - ( ( p(LL) - p(LL-1) ) &
                  / LOG( p(LL) / p(LL-1) ) ) )
       ELSE
          test1 = ABS( -1e0*p(LL) + ( ( p(LL+1) - p(LL) ) &
                  / LOG( p(LL+1) / p(LL) ) ) )
          test2 = ABS( p(LL) - ( ( p(LL) - p(LL-1) ) &
                  / LOG( p(LL) / p(LL-1) ) ) )
       ENDIF

       ! Get the three layers (bottom, top, and center)
       LAYER_BOT(LL) = ( TCC_P(LL) + test1 )
       LAYER_TOP(LL) = ( TCC_P(LL) - test2 )
    ENDDO

    ! Get the edges
    DO LL = 1, LTCC+1

       IF (LL .EQ. LTCC+1) THEN
          PEDGE(LL) = LAYER_BOT(LL-1)
       ELSE
          PEDGE(LL) = LAYER_TOP(LL)
       ENDIF

    ENDDO

    ! Return to the original grid formatting
    IF (LTCC .GT. 1) THEN
       IF (TCC_P_I(2) .LT. TCC_P_I(1)) PEDGE = PEDGE(LTCC+1:1:-1)
    ENDIF

    ! Store the output
    TCC_PEDGE = PEDGE

  END FUNCTION TCC_PEDGE
!EOC
END MODULE TCCON_CH4_MOD
