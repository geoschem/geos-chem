!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gosat_ch4_mod.F
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Parameters
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
!  (2 ) Add calculation of S_OER_INV, though eventually we probably want to 
!        do this offline. (dkh, 05/04/10) 
!  02 May 2016 - M. Sulprizio- Update code to read data from SRON
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

      caller = 'READ_GOS_CH4_OBS in gosat_ch4_mod.F'

      ! Get current year
      YEAR = GET_YEAR()
      WRITE( CYEAR, '(i4)' ) YEAR

      ! Filename
      nc_file = 'ESACCI-GHG-L2-CH4-GOSAT-OCPR-YYYYMMDD-fv7.0.nc'
      CALL Expand_Date( nc_file, YYYYMMDD, 9999 ) 

      ! Construct complete file path
      dir = '/n/regal/jacob_lab/msulprizio/ExtData/Obs/ULGOSAT_v7/' //
     &       CYEAR // '/' 
      nc_file = TRIM( dir ) // TRIM( nc_file )
      WRITE( 6, 10 ) TRIM( nc_file )
 10   FORMAT( '     - Reading ', a)

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
         print*, ' MAXLEV in gosat_ch4_mod.f.'
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
 20   FORMAT( '     - Finished reading GOSAT CH4 observations for ', i8)

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
      SUBROUTINE CALC_GOSAT_CH4_FORCE( Input_Opt, State_Chm, State_Met,
     &                                 COST_FUNC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE ERROR_MOD,          ONLY : IT_IS_FINITE
      USE GC_GRID_MOD,        ONLY : GET_IJ
      USE Input_Opt_Mod,      ONLY : OptInput
      USE JULDAY_MOD,         ONLY : CALDATE,     JULDAY
      USE TIME_MOD
      USE PhysConstants,      ONLY : XNUMOLAIR, AIRMW
      USE State_Chm_Mod,      ONLY : ChmState, Ind_
      USE State_Met_Mod,      ONLY : MetState
!------------------------------------------------------------------------------
! For adjoint model only:
!      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
!      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
!      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
!      USE CHECKPT_MOD,        ONLY : CHK_STT
!      USE COMODE_MOD,         ONLY : CSPEC, JLOP
!      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
!      USE LOGICAL_ADJ_MOD,    ONLY : LANAL_INV
!      USE ADJ_ARRAYS_MOD,     ONLY : EMS_PERT_NUM
!      USE ADJ_ARRAYS_MOD,     ONLY : PERT_EMS
!      USE ADJ_ARRAYS_MOD,     ONLY : NCLUSTEM
!      USE ADJ_ARRAYS_MOD,     ONLY : OBS_STT
!      USE ADJ_ARRAYS_MOD,     ONLY : SATELLITE_BIAS, MEAN_MODEL_BIAS
!------------------------------------------------------------------------------
!
! !INPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
      TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
      TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      REAL(fp), INTENT(INOUT)    :: COST_FUNC   ! Cost function [unitless]
!
! !REVISION HISTORY:
!  16 Jun 2017 - M. Sulprizio- Initial version based on GOSAT CH4 observation
!                              operator from GC Adjoint v35j with updates from
!                              M. Sulprizio, J.D. Maasakkers, A. Turner, and
!                              K. Wecht
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
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
      REAL(fp)           :: GC_PRES(LLPAR)
      REAL(fp)           :: GC_PEDGE(LLPAR+1)
      REAL(fp)           :: GC_CH4_NATIVE(LLPAR)
      REAL(fp)           :: GC_CH4(MAXLEV)
      REAL(fp)           :: GC_CH4_cm2(MAXLEV)
      REAL(fp)           :: GC_PSURF
      REAL(fp)           :: MAP(LLPAR,MAXLEV)
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

!------------------------------------------------------------------------------
! For adjoint model only:
!      ! For adjoint of observation operator
!      REAL(fp)           :: DIFF_ADJ
!      REAL(fp)           :: GC_XCH4_ADJ
!      REAL(fp)           :: GC_CH4_ADJ(MAXLEV)
!      REAL(fp)           :: GC_CH4_cm2_ADJ(MAXLEV)
!      REAL(fp)           :: GC_CH4_NATIVE_ADJ(LLPAR)
!
!      ! For the observational error matrix
!      REAL*4             :: S_OBS_MAT(IIPAR,JJPAR)
!      INTEGER            :: FID, VID
!
!      ! For Finite Difference Testing
!      LOGICAL, SAVE      :: DO_FDTEST = .FALSE.
!      LOGICAL, SAVE      :: FDTEST_MATCH = .FALSE. ! Do not change
!      REAL(fp)           :: PERT(LLPAR)
!      REAL(fp)           :: COST_FUNC_0
!      REAL(fp)           :: COST_FUNC_pos, COST_FUNC_neg
!      REAL(fp)           :: ADJ(LLPAR), ADJ_SAVE(LLPAR)
!      REAL(fp)           :: FD_CEN(LLPAR), FD_POS(LLPAR)
!      REAL(fp)           :: FD_NEG(LLPAR)
!------------------------------------------------------------------------------

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

!------------------------------------------------------------------------------
! For adjoint model only:
!         !kjw for testing adjoint of obs operator
!         FILENAME = 'test_adjoint_obs.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
!         FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         OPEN( 116,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
!     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )
!------------------------------------------------------------------------------

         ! For recording model and observation values
!         FILENAME = 'sat_obs.gosat.NN.m'
!         CALL EXPAND_NAME( FILENAME, N_CALC )
         FILENAME = 'sat_obs.gosat.00.m'
         FILENAME = TRIM( Input_Opt%RUN_DIR ) //  TRIM( FILENAME )
         OPEN( 117,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN',
     &       IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

         ! Write header of sat_obs.NN.m
         WRITE( 117, 281 ) '       NNN',
     &                     '   I', '   J', '     LON','     LAT','YYYY',
     &                     'MM', 'DD', 'hh', 'mm', 'ss',
     &                     '         TAU', '       GOSAT', 
     &                     '       model', '       S_OBS', 
     &                     '    COST_FUN', 'GLINT', 'GAIN',
     &                     'EXPOSURE_ID'
 281     FORMAT( A10,2x,A4,2x,A4,2x,A8,2x,A8,2x,A4,2x,A2,2x,A2,2x,A2,2x,
     &           A2,2x,A2,2x,A12,2x,A12,2x,A12,2x,A12,2x,A12,2x,
     &           A5,2x,A4,2x,A22)

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

!------------------------------------------------------------------------------
! For adjoint model only:
!         ! Initialize the observational Error matrix
!         S_OBS_MAT(:,:) = -1d0
!------------------------------------------------------------------------------

         IF ( FIRST ) FIRST = .FALSE. 

      ENDIF 

      ! Get indices of GOS observations in the current hour
      !   At the start of each hour, assimilate observations that 
      !   were made in the previous 60 minutes.
      !   For example, at time 18:00, assimilate observations 
      !   made from 18:00 - 18:59
      INDS(:) = 0
      NOBS    = 0
!      print*,'Looking for observations at MONTH, DAY, HOUR = ',
!     &          GET_MONTH(), GET_DAY(), GET_HOUR()
      DO NT = 1, NGOS
         IF ( GOS(NT)%MONTH(1) .EQ. GET_MONTH() .AND.
     &        GOS(NT)%DAY(1)   .EQ. GET_DAY()   .AND.
     &        GOS(NT)%HOUR(1)  .EQ. GET_HOUR()  ) THEN
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
!!$OMP PARALLEL DO
!!$OMP+DEFAULT( PRIVATE )
!!!$OMP+PRIVATE( IND, NT, MAP, LGOS, IIJJ,  I, J,  L,   LL, JLOOP )
!!!$OMP+PRIVATE( GC_CH4, FORCE, CH4_PRIOR, GC_PRES, FILENAME      )
!!!$OMP+PRIVATE( GC_PEDGE, GC_PSURF, GC_CH4_NATIVE, GOS_XCH4      )
!!!$OMP+PRIVATE( GOS_XCH4_ERROR, S_OBS, h, p, XCH4a, XCH4m        )
!!!$OMP+PRIVATE( GC_XCH4, DIFF, DIFF_ADJ, GC_XCH4_ADJ             )
!!!$OMP+PRIVATE( GC_CH4_NATIVE_ADJ, GC_CH4_ADJ, TotalObs          )
      DO IND = 1, NOBS

         NT = INDS(IND)
         !print*, '     - CALC_GOS_CH4_FORCE: analyzing record ', NT 

         ! Quality screening (0=good, 1=bad)
         IF ( GOS(NT)%QFLAG(1) .NE. 0 ) THEN 
            print*, ' Bad QFLAG, skipping record         ', NT
            CYCLE
         ENDIF

!         ! Skip glint mode observations (0 = land, 1 = glint)
!         IF ( GOS(NT)%GLINT(1) .EQ. 1 ) THEN 
!            print*, ' GLINT obs, skipping record         ', NT
!            CYCLE
!         ENDIF

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

         ! skip Observations outside the nested domain
         ! Only skip if we're running nested (ajt, 1/26/13)
#if   defined( NESTED_NA )
         IF ( GOS(NT)%LAT(1) <   10e0 .OR. GOS(NT)%LAT(1) >  70e0  .OR. 
     &        GOS(NT)%LON(1) < -140e0 .OR. GOS(NT)%LON(1) > -40e0 ) THEN
            print*, ' Outside nested domain, skipping record ', NT
            CYCLE
         ENDIF
#endif

         ! Skip data at high lat (ajt, 7/23/13)
         !JDMIF ( GOS(NT)%LAT(1) >   60d0 ) THEN
         !JDMprint*, ' Skip all data above 60N ', NT
         !JDMCYCLE
         !JDMENDIF
	 !Turn of this for now to reproduce plots made by AJT
	 !Should probably turn this back on with SH once we start
         !To do the actual adjoint run

         ! Get grid box of current record
         IIJJ = GET_IJ( REAL(GOS(NT)%LON(1),4), REAL(GOS(NT)%LAT(1),4))
         I    = IIJJ(1)
         J    = IIJJ(2)

         ! skip observations where the GOSAT surface pressure is 
         ! less than the model
         IF ( (GOS(NT)%PRES(1) - State_Met%PEDGE(I,J,1)) > 50e0 ) THEN
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
         LGOS = GOS(NT)%LGOS(1)
         DO L = 1, LGOS
            CH4_PRIOR(L) = GOS(NT)%PRIOR(L) * 1d-9 ! [ppb] --> [v/v]
         ENDDO

         ! Get GC pressure levels (mbar) 
         DO L = 1, LLPAR
            GC_PRES(L) = State_Met%PMID(I,J,L)
         ENDDO

         ! Get GC pressure edges (mbar) 
         DO L = 1, LLPAR+1
            GC_PEDGE(L) = State_Met%PEDGE(I,J,L)
         ENDDO

         ! Get GC surface pressure (mbar) 
         GC_PSURF = State_Met%PEDGE(I,J,1) 

         ! Calculate the interpolation weight matrix 
         MAP(:,:) = 0.0_fp
         CALL GET_INTMAP( GC_PRES, GC_PSURF, GOS(NT)%PRES, LGOS, MAP )

!------------------------------------------------------------------------------
! For adjoint model only:
!         IF ( DO_FDTEST .AND. (NT .EQ. NT_FD) ) THEN
!         print*,'---------------------------------------'
!         print*,'LGOS = ',LGOS
!         FDTEST_MATCH = .TRUE.
!         WRITE(6,'(20x,20F16.8)') GOS(NT)%PRES
!         DO LL=1,LLPAR
!            WRITE(6,'(I4, 21F16.8)') LL, GET_PCENTER(I,J,LL), MAP(LL,:)
!         ENDDO
!         print*,'---------------------------------------'
!         ENDIF
!------------------------------------------------------------------------------

         ! Get CH4 values at native model resolution
         GC_CH4_NATIVE(:) = 0.0_fp

!------------------------------------------------------------------------------
! For adjoint model only:
!         ! Use actual STT in forward mode
!         IF (LANAL_INV) THEN
!            GC_CH4_NATIVE(:) = OBS_STT(I,J,:,IDTCH4)
!         ELSE
!            GC_CH4_NATIVE(:) = CHK_STT(I,J,:,IDTCH4)
!         ENDIF
!------------------------------------------------------------------------------
         ! Get species concentrations
         ! Convert from [kg/box] --> [v/v]
         GC_CH4_NATIVE(:) = State_Chm%Species(I,J,:,id_CH4) * ( AIRMW / 
     &                      State_Chm%SpcData(id_CH4)%Info%emMW_g )

 
         ! Interpolate GC CH4 column to GOSAT grid
         DO LL = 1, LGOS
            GC_CH4(LL) = 0.0_fp
            DO L = 1, LLPAR 
               GC_CH4(LL) = GC_CH4(LL) 
     &                    + MAP(L,LL) * GC_CH4_NATIVE(L) 
            ENDDO
         ENDDO

         ! GOSAT Proxy XCH4 [ppb] --> [v/v]
         GOS_XCH4       = GOS(NT)%CH4(1)       * 1e-9_fp
         GOS_XCH4_ERROR = GOS(NT)%CH4_ERROR(1) * 1e-9_fp

         ! Remove any GOSAT bias
         GOS_XCH4 = GOS_XCH4 + 1e-9 * ( SATELLITE_BIAS(1) 
     &                       + SATELLITE_BIAS(2)*(GOS(NT)%LAT(1))
     &                       + SATELLITE_BIAS(3)*(GOS(NT)%LAT(1))**2 )

         ! Get the S_obs, assume stddev adds in quadrature, variance
         ! adds linearly.  (ajt, 03/27/2013) 
         S_OBS = GOS_XCH4_ERROR**2 + (GC_XCH4_ERROR * 1e-9_fp)**2

!------------------------------------------------------------------------------
! EXT_OBS_MAT is hardcoded to FALSE at beginning of routine (mps, 6/16/17)
!         ! Specify the observational error externally
!         ! (ajt, 11/5/13)
!         IF (EXT_OBS_MAT) THEN 
!     
!            ! Read in the obs error matrix
!            FILENAME = '../ext_sat_obs.gosat.nc'
!            FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!            CALL NCDF_OPEN_FOR_READ( FID, TRIM( FILENAME ) )
!            VID = ncdf_get_varid( FID, 'obs_error' )
!            CALL NCDF_GET_VAR( FID, VID, S_OBS_MAT )
!            CALL NCDF_CLOSE( FID )
!
!            ! Get our observational error
!            IF ( (S_OBS_MAT(I,J) .GT.                  0d0)   .AND.
!     &           (S_OBS_MAT(I,J) .GT. GOS(NT)%CH4_ERROR(1)) ) THEN 
!               S_OBS = (S_OBS_MAT(I,J) * 1d-9)**2
!            ENDIF
!     
!         ENDIF
!------------------------------------------------------------------------------

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
         h(:)          = 0.0_fp
         p(1:LGOS)     = GOS(NT)%PRES(1:LGOS)
         ak(1:LGOS)    = GOS(NT)%AVG_KERNEL(1:LGOS)
         prior(1:LGOS) = GOS(NT)%PRIOR(1:LGOS) * 1e-9_fp  ! [ppb] --> [v/v]

         ! Need to integrate from the toa to surface (ajt, 05/21/13)
         IF (LGOS .GT. 1) THEN
            IF(GOS(NT)%PRES(2) .LT. GOS(NT)%PRES(1)) THEN
               p(1:LGOS) = p(LGOS:1:-1)
            ENDIF
         ENDIF

         L = 1
         h(L) = 1./GOS(NT)%PRES(1) * ABS( 
     &         ( -1e0*p(L) + ( (p(L+1)-p(L))/(LOG(p(L+1)/p(L))) ) ) )
         L = LGOS
         h(L) = 1./GOS(NT)%PRES(1) * ABS( 
     &         (  p(L) - ( (p(L)-p(L-1))/(LOG(p(L)/p(L-1))) ) ) )
         DO L=2,LGOS-1
            h(L) = 1./GOS(NT)%PRES(1) * ABS( 
     &         ( -1e0*p(L) + ( (p(L+1)-p(L))/(LOG(p(L+1)/p(L))) ) ) +  
     &         (      p(L) - ( (p(L)-p(L-1))/(LOG(p(L)/p(L-1))) ) )   )
         ENDDO

         ! Now return to the orientation of the other variables
         IF (LGOS .GT. 1) THEN 
            IF(GOS(NT)%PRES(2) .LT. GOS(NT)%PRES(1)) THEN 
               h(1:LGOS) = h(LGOS:1:-1)
               p(1:LGOS) = p(LGOS:1:-1)
            ENDIF
         ENDIF

         ! Calculate Xch4_a
         XCH4a = 0.0_fp
         DO L = 1,LGOS
            XCH4a = XCH4a + h(L) * prior(L)
         ENDDO

         ! Calculate Xch4_m
         XCH4m = 0.0_fp
         XCH4m = XCH4a
         DO L = 1, LGOS
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

!------------------------------------------------------------------------------
! For adjoint model only:
!         !--------------------------------------------------------------
!         ! Begin adjoint calculations 
!         !--------------------------------------------------------------
!
!         ! The adjoint forcing  is S_{obs}^{-1} * DIFF = FORCE
!         DIFF_ADJ = FORCE
!
!         ! Adjoint of difference
!         GC_XCH4_ADJ = DIFF_ADJ
!
!         ! Adjoint of GOSAT observation operator
!         GC_CH4_ADJ(:) = 0d0
!         DO LL = 1, LGOS
!            !GC_CH4_ADJ(LL) = h(LL) * ak(LL)*GC_XCH4_ADJ
!            GC_CH4_ADJ(LL) = pres_w(LL) * ak(LL)*GC_XCH4_ADJ
!         ENDDO
! 
!! Old interpolation from kjw
!         ! adjoint of interpolation 
!         DO L  = 1, LLPAR
!            GC_CH4_NATIVE_ADJ(L) = 0d0 
!            DO LL = 1, LGOS
!               GC_CH4_NATIVE_ADJ(L) = GC_CH4_NATIVE_ADJ(L)
!     &                              + MAP(L,LL) * GC_CH4_ADJ(LL)
!            ENDDO
!         ENDDO
!         
!!         ! Updated adjoint of interpolation
!!         GC_CH4_NATIVE_ADJ = INTERP_GC_GOS( LLPAR, LGOS,
!!     &                       GOS(NT)%PRES(1:LGOS),
!!     &                       GOS_PEDGE( LGOS, GOS(NT)%PRES(1:LGOS) ),
!!     &                       GC_PRES, GC_CH4_ADJ )
!
!         ! Adjoint of unit conversion 
!         GC_CH4_NATIVE_ADJ(:) = GC_CH4_NATIVE_ADJ(:) * TCVV_CH4
!     &                        / AD(I,J,:) 
!
!        ! Pass adjoint back to adjoint tracer array
!         STT_ADJ(I,J,:,IDTCH4)  = STT_ADJ(I,J,:,IDTCH4)
!     &                          + GC_CH4_NATIVE_ADJ(:)
!------------------------------------------------------------------------------

         TotalObs = TotalObs + 1

         ! Record information for satellite diagnostics
         IF ( LDCH4SAT ) THEN 
            WRITE( 117, 282 ) TotalObs, I, J, GOS(NT)%LON(1),
     &           GOS(NT)%LAT(1),GOS(NT)%YEAR(1), 
     &           GOS(NT)%MONTH(1),GOS(NT)%DAY(1), GOS(NT)%HOUR(1), 
     &           GOS(NT)%MINUTE(1), GOS(NT)%SEC(1),
     &           GET_TAU(), GOS_XCH4, GC_XCH4, S_OBS, 0.5d0*FORCE*DIFF,
     &           GOS(NT)%GLINT(1),GOS(NT)%GAIN(1),GOS(NT)%EXP_ID(1)
         ENDIF



      ENDDO  ! NT
!!$OMP END PARALLEL DO

      ! Number of observations processed in total
      !TotalObs = TotalObs + NOBS

      ! Update cost function 
      COST_FUNC = COST_FUNC + SUM(NEW_COST(:))

!==============================================================================
!  For adjoint model only:
!! -----------------------------------------------------------------------
!!    Use this section to test the adjoint of the GOSAT_CH4 operator by
!!          slightly perturbing model [CH4] and recording resultant change
!!          in calculated contribution to the cost function.
!!
!!    This routine will write the following information for each observation
!!          to rundir/diagadj/test_adjoint_obs.NN.m
!!
!!    The adjoint of the observation operator has been tested and validated
!!          as of 11/06/13, ajt.
!!
!!      !IF (( DO_FDTEST ) .AND. ( nboxes .gt. 1 )) THEN
!      IF ( DO_FDTEST  .AND. FDTEST_MATCH ) THEN
!
!      print*,'About to start FD testing ...'
!      print*,' -> DO_FDTEST    = ', DO_FDTEST
!      print*,' -> FDTEST_MATCH = ', FDTEST_MATCH
!
!      PERT(:) = 0D0
!      NT = NT_FD ! Just in case I forgot to change all the NT's
!
!      ! Unperturbed Cost function
!      COST_FUNC_0 = 0d0
!      CALL CALC_GOSAT_CH4_FORCE_FD( COST_FUNC_0, PERT, ADJ, NT )
!      ADJ_SAVE(:) = ADJ(:)
!
!      ! Get grid box of current record
!      IIJJ  = GET_IJ( REAL(GOS(NT)%LON(1),4), REAL(GOS(NT)%LAT(1),4))
!      I     = IIJJ(1)
!      J     = IIJJ(2)
!
!      ! Get GC pressure levels (mbar) 
!      DO L = 1, LLPAR
!         GC_PRES(L) = GET_PCENTER(I,J,L)
!      ENDDO
!
!      ! Get GC surface pressure (mbar) 
!      GC_PSURF = GET_PEDGE(I,J,1) 
!
!      ! Write identifying information to top of satellite diagnostic file
!      WRITE(116,'(a)') '----------------------------------------------'
!      WRITE(116,212) 'GC_PSURF    ', GC_PSURF
!      WRITE(116,212) 'GOSAT PSURF  ', GOS(NT_FD)%PRES(1)
!      WRITE(116,212) 'NEW_COST:   ', NEW_COST(NT_FD)
!      WRITE(116,212) 'COST_FUNC_0:', COST_FUNC_0
!      WRITE(116,'(a)') '----------------------------------------------'
!      WRITE(116,'(a)') ' '
!      WRITE(116,210) '  LG'        ,! ' box',       '  TROP',
!     &               '     GC_PRES',
!     &               '      FD_POS', '      FD_NEG', '      FD_CEN',
!     &               '         ADJ', '    COST_POS', '    COST_NEG',
!     &               '  FD_POS/ADJ', '  FD_NEG/ADJ', '  FD_CEN/ADJ'
!
!      ! Perform finite difference testing at each vertical level
!      ! and for each horizontal grid box in this observation
!      DO LL = 1, LLPAR
!
!         ! Positive perturbation to GEOS-Chem CH4 columns
!         PERT(:) = 0.0
!         PERT(LL) = 1D-2*PMAG
!         COST_FUNC_pos = 0D0
!         CALL CALC_GOSAT_CH4_FORCE_FD( COST_FUNC_pos, PERT, ADJ, NT )
!
!         ! Negative perturbation to GEOS-Chem CH4 columns
!         PERT(:) = 0.0
!         PERT(LL) = -1D-2*PMAG
!         COST_FUNC_neg = 0D0
!         CALL CALC_GOSAT_CH4_FORCE_FD( COST_FUNC_neg, PERT, ADJ, NT )
!
!         ! Calculate dJ/dCH4 from perturbations
!         FD_CEN(LL)   = ( COST_FUNC_pos - COST_FUNC_neg ) / (2D-2*PMAG)
!         FD_POS(LL)   = ( COST_FUNC_pos - COST_FUNC_0 )   / (1D-2*PMAG)
!         FD_NEG(LL)   = ( COST_FUNC_0 - COST_FUNC_neg )   / (1D-2*PMAG)
!
!         ! Write information to satellite diagnostic file
!         WRITE(116, 211)  LL,     GC_PRES(LL),
!     &                    FD_POS(LL),   FD_NEG(LL),
!     &                    FD_CEN(LL),   ADJ_SAVE(LL), 
!     &                    COST_FUNC_pos, COST_FUNC_neg, 
!     &                    FD_POS(LL)/ADJ_SAVE(LL),
!     &                    FD_NEG(LL)/ADJ_SAVE(LL),
!     &                    FD_CEN(LL)/ADJ_SAVE(LL)
!
!      ENDDO
!
!      DO_FDTEST    = .FALSE.
!      FDTEST_MATCH = .FALSE.
!
!      ENDIF   ! End if DO_FDTEST
!
! 210  FORMAT(A4,2x,A12,2x,A12,2x,A12,2x,A12,2x,A12,2x,A12,2x,A12,2x,A12,
!     &       2x,A12,2x,A12)
! 211  FORMAT(I4,2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6,
!     &        2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6,2x,D12.6)
! 212  FORMAT(A12,F22.6)
!==============================================================================

 282  FORMAT( I10,2x,I4,2x,I4,2x,F8.3,2x,F8.4,2x,I4,2x,I2,2x,I2,2x,I2,
     &        2x,I2,2x,I2,2x,F12.3,2x,E12.6,2x,E12.6,2x,E12.6,2x,E12.6,
     &        2x,I5,2x,I5,2x,A22)

      print*, ' Updated value of COST_FUNC = ', COST_FUNC 
      print*, ' GOS contribution           = ', COST_FUNC - OLD_COST  
      print*, ' Number of observations this hour = ', NOBS
      print*, ' Number of observations total     = ', TotalObs

      END SUBROUTINE CALC_GOSAT_CH4_FORCE
!EOC
!==============================================================================
! For adjoint model only:
!!------------------------------------------------------------------------------
!
!      SUBROUTINE CALC_GOSAT_CH4_FORCE_FD( COST_FUNC_A, PERT, ADJ, NT )
!!
!!******************************************************************************
!!  Subroutine CALC_GOS_CH4_FORCE calculates the adjoint forcing from the GOSAT
!!  CH4 observations and updates the cost function. (dkh, 10/12/10) 
!! 
!!
!!  Arguments as Input/Output:
!!  ============================================================================
!!  (1 ) COST_FUNC (REAL(fp)) : Cost funciton                        [unitless]
!!     
!!     
!!  NOTES:
!!  
!!******************************************************************************
!!
!      ! Reference to f90 modules
!      USE ADJ_ARRAYS_MOD,     ONLY : STT_ADJ
!      USE ADJ_ARRAYS_MOD,     ONLY : N_CALC
!      USE ADJ_ARRAYS_MOD,     ONLY : EXPAND_NAME
!      USE CHECKPT_MOD,        ONLY : CHK_STT
!      USE COMODE_MOD,         ONLY : CSPEC, JLOP
!      USE DAO_MOD,            ONLY : AD, TROPP
!      USE DAO_MOD,            ONLY : AIRDEN
!      USE DAO_MOD,            ONLY : BXHEIGHT
!      USE DIRECTORY_ADJ_MOD,  ONLY : DIAGADJ_DIR
!      USE GRID_MOD,           ONLY : GET_IJ
!      USE PRESSURE_MOD,       ONLY : GET_PCENTER, GET_PEDGE
!      USE TIME_MOD,           ONLY : GET_NYMD,    GET_NHMS
!      USE TIME_MOD,           ONLY : GET_YEAR,    GET_MONTH
!      USE TIME_MOD,           ONLY : GET_TS_CHEM, ITS_A_NEW_MONTH
!      USE TIME_MOD,           ONLY : ITS_A_NEW_DAY
!      USE TRACER_MOD,         ONLY : XNUMOLAIR
!      USE TROPOPAUSE_MOD,     ONLY : ITS_IN_THE_TROP
!
!      !ajt. Analytical inversion
!      USE LOGICAL_ADJ_MOD,    ONLY : LANAL_INV
!      
!      ! ajt. Forward mode
!      USE ADJ_ARRAYS_MOD,     ONLY : OBS_STT
!
!      ! Get the mean bias
!      USE ADJ_ARRAYS_MOD,     ONLY : SATELLITE_BIAS, MEAN_MODEL_BIAS
!
!      !ajt for netcdf
!      USE NETCDF_UTIL_MOD,    ONLY : NCDF_OPEN_FOR_READ
!      USE NETCDF_UTIL_MOD,    ONLY : NCDF_GET_VARID
!      USE NETCDF_UTIL_MOD,    ONLY : NCDF_GET_VAR
!      USE NETCDF_UTIL_MOD,    ONLY : NCDF_CLOSE
!
!#     include      "CMN_SIZE"      ! Size params
!
!      ! Arguments
!      REAL(fp), INTENT(INOUT)       :: COST_FUNC_A
!      REAL(fp), INTENT(OUT)         :: ADJ(LLPAR)
!      REAL(fp), INTENT(IN)          :: PERT(LLPAR)
!      INTEGER, INTENT(IN)         :: NT
!
!      ! Local variables for observation operator
!      INTEGER                     :: NTSTART, NTSTOP
!      INTEGER                     :: IIJJ(2), I,      J
!      INTEGER                     :: L,       LL,     LGOS
!      INTEGER                     :: JLOOP,   NOBS,   IND
!      INTEGER                     :: INDS(MAXGOS)
!      REAL(fp)                      :: GC_PRES(LLPAR)
!      REAL(fp)                      :: GC_PEDGE(LLPAR+1)
!      REAL(fp)                      :: GC_CH4_NATIVE(LLPAR)
!      REAL(fp)                      :: GC_CH4(MAXLEV)
!      REAL(fp)                      :: GC_CH4_cm2(MAXLEV)
!      REAL(fp)                      :: GC_PSURF
!      REAL(fp)                      :: MAP(LLPAR,MAXLEV)
!      REAL(fp)                      :: CH4_HAT(MAXLEV)
!      REAL(fp)                      :: NEW_COST(MAXGOS)
!      REAL(fp)                      :: OLD_COST
!      REAL(fp), SAVE                :: TIME_FRAC(MAXGOS)
!      INTEGER,SAVE                :: NGOS
!      REAL(fp)                      :: frac, RHS, LHS
!      REAL(fp)                      :: CH4_PRIOR(MAXLEV)
!      REAL(fp)                      :: CH4_PRIOR_cm2(MAXLEV)
!      REAL(fp)                      :: molecongos(MAXLEV)
!      REAL(fp)                      :: GOS_XCH4, GC_XCH4, GOS_XCH4_ERROR
!      REAL(fp)                      :: p(MAXLEV), h(MAXLEV)
!      REAL(fp)                      :: ak(MAXLEV), prior(MAXLEV)
!      REAL(fp)                      :: pres_w(MAXLEV)
!      REAL(fp)                      :: XCH4m, XCH4a
!
!      ! For adjoint of observation operator
!      REAL(fp)                      :: FORCE
!      REAL(fp)                      :: DIFF
!      REAL(fp)                      :: DIFF_ADJ
!      REAL(fp)                      :: GC_XCH4_ADJ
!      REAL(fp)                      :: GC_CH4_ADJ(MAXLEV)
!      REAL(fp)                      :: GC_CH4_cm2_ADJ(MAXLEV)
!      REAL(fp)                      :: GC_CH4_NATIVE_ADJ(LLPAR)
!      REAL(fp)                      :: GC_CH4_ori(MAXLEV)
!      REAL(fp)                      :: S_OBS
!
!      ! For Observational error
!      REAL(fp)                      :: S_OBS_MAT(IIPAR,JJPAR)
!      INTEGER                     :: FID, VID
!
!      ! For miscellaneous
!      INTEGER                     :: IOS, TotalObs
!      CHARACTER(LEN=255)          :: FILENAME
!
!
!
!      !=================================================================
!      ! CALC_GOS_CH4_FORCE_FD begins here!
!      !=================================================================
!
!      print*, '     - CALC_GOS_CH4_FORCE_FD '
!    
!      ! Reset
!      COST_FUNC_A = 0d0
!      ADJ(:)      = 0d0
!
!      ! For safety, initialize these up to LLGOS
!      LGOS            = 0
!      GC_CH4(:)       = 0d0 
!      MAP(:,:)        = 0d0 
!      FORCE           = 0d0 
!      CH4_PRIOR(:)    = 0d0
!
!      ! Initialize the observational error matrix (-1 if not read
!      ! before)
!      S_OBS_MAT(:,:) = -1d0
!
!      ! Copy LGOS to make coding a bit cleaner
!      LGOS = GOS(NT)%LGOS(1)
!      DO L = 1, LGOS
!         CH4_PRIOR(L) = GOS(NT)%PRIOR(L) * 1d-9 ! [ppb] --> [v/v]
!      ENDDO
!
!      ! Get grid box of current record
!      IIJJ  = GET_IJ( REAL(GOS(NT)%LON(1),4), REAL(GOS(NT)%LAT(1),4))
!      I     = IIJJ(1)
!      J     = IIJJ(2)
!
!      ! Get GC pressure levels (mbar) 
!      DO L = 1, LLPAR
!         GC_PRES(L) = GET_PCENTER(I,J,L)
!      ENDDO
!
!      ! Get GC pressure edges (mbar) 
!      DO L = 1, LLPAR+1           
!         GC_PEDGE(L) = GET_PEDGE(I,J,L)
!      ENDDO
!
!      ! Get GC surface pressure (mbar) 
!      GC_PSURF = GET_PEDGE(I,J,1) 
!
!
!      ! Calculate the interpolation weight matrix 
!      MAP(:,:) = 0d0
!      CALL GET_INTMAP( GC_PRES, GC_PSURF, GOS(NT)%PRES, LGOS, MAP )
!
!      ! Get CH4 values at native model resolution
!      GC_CH4_NATIVE(:) = 0d0
!      IF (LANAL_INV) THEN
!         GC_CH4_NATIVE(:) = OBS_STT(I,J,:,IDTCH4)
!      ELSE
!         GC_CH4_NATIVE(:) = CHK_STT(I,J,:,IDTCH4)
!      ENDIF
!
!      ! GOSAT v4 retrieval crashes without a print statement here. This
!      ! wasn't an issue with v3 retrievals (saved in monthly blocks)...
!      ! (ajt, 03/25/2013)
!      !print*, 'code too fast here.'
!
!      ! Perturb level(s) of GC_CH4_NATIVE, given by input argument PERT
!      DO L = 1, LLPAR
!         GC_CH4_NATIVE(L) = GC_CH4_NATIVE(L) * ( 1 + PERT(L) )
!      ENDDO
!
!      ! Convert from kg/box to v/v
!      GC_CH4_NATIVE(:) = GC_CH4_NATIVE(:) * TCVV_CH4
!     &                    / AD(I,J,:) 
! 
!!      ! Interpolate GC CH4 column to GOSAT grid
!!      ! use the method of Parker et al., (2011)
!!      GC_CH4 = INTERP_GC_GOS( LGOS, LLPAR, GC_PRES, GC_PEDGE,
!!     &         GOS(NT)%PRES(1:LGOS), GC_CH4_NATIVE )
!
!! Old interpolation method from kjw
!      ! Interpolate GC CH4 column to GOSAT grid 
!      DO LL = 1, LGOS
!         GC_CH4(LL) = 0d0 
!         DO L = 1, LLPAR 
!            GC_CH4(LL) = GC_CH4(LL) 
!     &           + MAP(L,LL) * GC_CH4_NATIVE(L) 
!         ENDDO
!      ENDDO
!
!      ! GOSAT Proxy XCH4 [ppb] --> [v/v]
!      GOS_XCH4       = GOS(NT)%CH4(1) * 1d-9
!      GOS_XCH4_ERROR = GOS(NT)%CH4_ERROR(1) * 1d-9
!
!      ! Remove any GOSAT bias
!      GOS_XCH4 = GOS_XCH4 + 1d-9 * ( SATELLITE_BIAS(1) 
!     &                    + SATELLITE_BIAS(2)*(GOS(NT)%LAT(1))
!     &                    + SATELLITE_BIAS(3)*(GOS(NT)%LAT(1))**2 )
!
!      ! Get the S_obs, assume stddev adds in quadrature, variance
!      ! adds linearly.  (ajt, 03/27/2013) 
!      S_OBS = GOS_XCH4_ERROR**2 + (GC_XCH4_ERROR * 1d-9)**2
!
!      ! Specify the observational error externally
!      ! (ajt, 11/5/13)
!      IF (EXT_OBS_MAT) THEN 
!              
!         ! Read in the obs error matrix
!         FILENAME = '../ext_sat_obs.gosat.nc'
!         FILENAME = TRIM( DIAGADJ_DIR ) //  TRIM( FILENAME )
!         CALL NCDF_OPEN_FOR_READ( FID, TRIM( FILENAME ) )
!         VID = ncdf_get_varid( FID, 'obs_error' )
!         CALL NCDF_GET_VAR( FID, VID, S_OBS_MAT )
!         CALL NCDF_CLOSE( FID )
!
!         ! Get our observational error
!         IF ( (S_OBS_MAT(I,J) .GT.                  0d0)   .AND.
!     &        (S_OBS_MAT(I,J) .GT. GOS(NT)%CH4_ERROR(1)) ) THEN 
!            S_OBS = (S_OBS_MAT(I,J) * 1d-9)**2
!         ENDIF
!  
!      ENDIF
!
!      !--------------------------------------------------------------
!      ! Apply GOSAT observation operator
!      !
!      !   Xch4_m = Xch4_a + SUM_j( h_j * a_j * (x_m - x_a) )
!      !  
!      !   Xch4_m  - model XCH4
!      !   Xch4_a  - apriori XCH4 = h^T * x_a
!      !   h       - pressure weighting function
!      !   a       - column averaging kernel
!      !   x_m     - model CH4 [v/v]
!      !   x_a     - apriori CH4 [v/v]
!      !
!      !   The pressure weighting function is defined in Connor et al. 2008 
!      !     and the OCO-2 ATBD
!      !--------------------------------------------------------------
!
!      ! Pressure weighting array
!      !h(:)          = 0d0
!      p(1:LGOS)     = GOS(NT)%PRES(1:LGOS)
!      ak(1:LGOS)    = GOS(NT)%AVG_KERNEL(1:LGOS)
!      pres_w(1:LGOS)= GOS(NT)%P_WEIGHT(1:LGOS)
!      prior(1:LGOS) = GOS(NT)%PRIOR(1:LGOS) * 1d-9  ! [ppb] --> [v/v]
!
!!      ! Need to integrate from the toa to surface (ajt, 05/21/13)
!!      IF (LGOS .GT. 1) THEN
!!         IF(GOS(NT)%PRES(2) .LT. GOS(NT)%PRES(1)) THEN
!!            p(1:LGOS) = p(LGOS:1:-1)
!!         ENDIF
!!      ENDIF
!!
!!      L = 1
!!      h(L) = 1./GOS(NT)%PRES(1) * ABS(
!!     &      ( -1D0*p(L) + ( (p(L+1)-p(L))/(LOG(p(L+1)/p(L))) ) ) )
!!      L = LGOS
!!      h(L) = 1./GOS(NT)%PRES(1) * ABS(
!!     &      (  p(L) - ( (p(L)-p(L-1))/(LOG(p(L)/p(L-1))) ) ) )
!!      DO L=2,LGOS-1
!!         h(L) = 1./GOS(NT)%PRES(1) * ABS(
!!     &      ( -1D0*p(L) + ( (p(L+1)-p(L))/(LOG(p(L+1)/p(L))) ) ) +
!!     &      (      p(L) - ( (p(L)-p(L-1))/(LOG(p(L)/p(L-1))) ) )   )
!!      ENDDO
!!
!!      ! Now return to the orientation of the other variables
!!      IF (LGOS .GT. 1) THEN
!!         IF(GOS(NT)%PRES(2) .LT. GOS(NT)%PRES(1)) THEN
!!            h(1:LGOS) = h(LGOS:1:-1)
!!            p(1:LGOS) = p(LGOS:1:-1)
!!         ENDIF
!!      ENDIF
!
!      ! Calculate Xch4_a
!      XCH4a = 0d0
!      DO L = 1,LGOS
!         !XCH4a = XCH4a + h(L) * prior(L)
!         XCH4a = XCH4a + pres_w(L) * prior(L)
!      ENDDO
!
!      ! Calculate Xch4_m
!      XCH4m = 0d0
!      XCH4m = XCH4a
!      DO L = 1, LGOS
!         !XCH4m = XCH4m + ( h(L) * ak(L) * ( GC_CH4(L) - prior(L) ) )
!         XCH4m = XCH4m + ( pres_w(L) * ak(L) *
!     &           ( GC_CH4(L) - prior(L) ) )
!      ENDDO
!      GC_XCH4 = 0d0
!      GC_XCH4 = XCH4m
!
!      ! Remove any GEOS-Chem bias
!      GC_XCH4 = GC_XCH4 + ( 1d-9 * MEAN_MODEL_BIAS )
!
!
!      !--------------------------------------------------------------
!      ! Calculate cost function, given S is error in vmr
!      ! J = 1/2 [ model - obs ]^T S_{obs}^{-1} [ model - obs ]
!      !--------------------------------------------------------------
!
!      ! Calculate difference between modeled and observed profile
!      DIFF = GC_XCH4 - GOS_XCH4
!
!      ! Calculate 1/2 * DIFF^T * S_{obs}^{-1} * DIFF
!      ! Need to account for the model error (ajt, 03/27/2013)
!      FORCE        = ( 1d0 / (S_OBS) ) * DIFF
!      COST_FUNC_A  = 0.5d0 * DIFF * FORCE
!
!
!      !--------------------------------------------------------------
!      ! Begin adjoint calculations 
!      !--------------------------------------------------------------
!
!      ! The adjoint forcing  is S_{obs}^{-1} * DIFF = FORCE
!      DIFF_ADJ = FORCE
! 
!      ! Adjoint of difference
!      GC_XCH4_ADJ = DIFF_ADJ
!
!
!      ! Adjoint of GOSAT observation operator
!      GC_CH4_ADJ(:) = 0d0
!      DO LL = 1, LGOS
!         !GC_CH4_ADJ(LL) = h(LL) * ak(LL) * GC_XCH4_ADJ
!         GC_CH4_ADJ(LL) = pres_w(LL) * ak(LL) * GC_XCH4_ADJ
!      ENDDO
!
!! Old interpolation from kjw
!      ! adjoint of interpolation 
!      DO L  = 1, LLPAR
!         GC_CH4_NATIVE_ADJ(L) = 0d0 
!         DO LL = 1, LGOS
!            GC_CH4_NATIVE_ADJ(L) = GC_CH4_NATIVE_ADJ(L)
!     &                           + MAP(L,LL) * GC_CH4_ADJ(LL)
!         ENDDO
!      ENDDO
!
!!      ! Updated adjoint of interpolation
!!      GC_CH4_NATIVE_ADJ = INTERP_GC_GOS_ADJ( LLPAR, LGOS,
!!     &                    GOS(NT)%PRES(1:LGOS),
!!     &                    GOS_PEDGE( LGOS, GOS(NT)%PRES(1:LGOS) ),
!!     &                    GC_PRES, GC_CH4_ADJ, MAP )
!
!      ! Adjoint of unit conversion 
!      GC_CH4_NATIVE_ADJ(:) = GC_CH4_NATIVE_ADJ(:) * TCVV_CH4
!     &                     / AD(I,J,:) 
!
!      ! Calculate adjoint sensitivity to values in CHK_STT
!      ADJ(:) = 0d0
!      DO LL=1,LLPAR
!         IF (LANAL_INV) THEN
!            ADJ(LL) = GC_CH4_NATIVE_ADJ(LL) * OBS_STT(I,J,LL,IDTCH4)
!         ELSE
!            ADJ(LL) = GC_CH4_NATIVE_ADJ(LL) * CHK_STT(I,J,LL,IDTCH4)
!         ENDIF
!      ENDDO
!
!
!
!      ! Return to calling program
!      END SUBROUTINE CALC_GOSAT_CH4_FORCE_FD
!
!!------------------------------------------------------------------------------
!
!      FUNCTION INTERP_GC_GOS( LGOS, LGC, GC_P_I, PEDGE_I, GOS_P_I, 
!     &                        GC_CH4_NATIVE_I )
!!
!!******************************************************************************
!!  Function INTERP_GC_GOS interpolated GEOS-Chem CH4 to the GOSAT pressure grid
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) LGOS            : Number of GOSAT pressure levels     [unit]
!!  (2 ) LGC             : Number of GEOS-Chem pressure levels [unit]
!!  (3 ) GC_P_I          : GEOS-Chem pressure centers          [mbar]
!!  (4 ) PEDGE_I         : GEOS-Chem pressure edges            [mbar]
!!  (5 ) GOS_P_I         : GOSAT pressure levels               [mbar]
!!  (6 ) GC_CH4_NATIVE_I : GEOS-Chem CH4 on the native grid    [kg/box]
!!     
!!  Arguments as Output:
!!  ============================================================================
!!  (1 ) INTERP_GC_GOS   : GEOS-Chem CH4 interpolated to GOSAT [v/v]
!!     
!!  NOTES:
!!  (1 ) Converts to Layer Column Density (LCD) before performing interpolation
!!       because LCD is a monotonically increasing function of pressure.
!!  (2 ) Recommended by Rob Parker from a lecture slide by Chris Barnet.
!!
!!******************************************************************************
!!
!      ! Arguments
!      INTEGER,  INTENT(IN)           :: LGOS, LGC
!      REAL(fp),  INTENT(IN)            :: GC_CH4_NATIVE_I(LGC)
!      REAL(fp), INTENT(IN)             :: GC_P_I(LGC)
!      REAL(fp), INTENT(IN)             :: PEDGE_I(LGC+1)
!      REAL(fp), INTENT(IN)             :: GOS_P_I(LGOS)
!
!      ! Parameters
!      REAL(fp), PARAMETER              :: AV_NUM = 6.0221418d23 
!      REAL(fp), PARAMETER              :: MOLWT  = 16.0426*1d-3 
!      REAL(fp), PARAMETER              :: GRAV   = 9.81d0
!
!      ! Return value 
!      REAL(fp)                         :: INTERP_GC_GOS(LGOS)
!
!      ! Local variables 
!      INTEGER                        :: L, LL
!      REAL(fp)                         :: SCD(LGC+1)
!      REAL(fp)                         :: GC_CH4_NATIVE(LGC)
!      REAL(fp)                         :: p(LGOS)
!      REAL(fp)                         :: GC_P(LGC)
!      REAL(fp)                         :: PEDGE(LGC+1)
!      REAL(fp)                         :: GOS_P(LGOS)
!      REAL(fp)                         :: DPRES
!      REAL(fp)                         :: TEMP_SCD, TOTAL_SCD
!      REAL(fp)                         :: test1, test2
!      REAL(fp)                         :: temp1(1), temp2(1)
!      REAL(fp)                         :: LAYER_BOT(LGOS)
!      REAL(fp)                         :: LAYER_CEN(LGOS)
!      REAL(fp)                         :: LAYER_TOP(LGOS)
!      REAL(fp)                         :: PAR_COL(LGOS)
!      REAL(fp)                         :: GC_CH4(LGOS)
!
!      !=================================================================
!      ! INTER_GC_GOS begins here!
!      !=================================================================
!
!      ! Convert hPa to Pascals
!      PEDGE = PEDGE_I * 1D2
!      GOS_P = GOS_P_I * 1D2
!      GC_P  = GC_P_I  * 1D2
!
!      ! Reverse the array orders so toa is first
!      IF (LGOS .GT. 1) THEN
!         IF (GOS_P_I(2) .LT. GOS_P_I(1)) THEN
!            GOS_P         = GOS_P(LGOS:1:-1)
!         ENDIF
!      ENDIF
!      IF (LGC .GT. 1) THEN
!         IF (GC_P_I(2) .LT. GC_P_I(1)) THEN
!            PEDGE         = PEDGE(LGC+1:1:-1)
!            GC_P          = GC_P(LGC:1:-1)
!            GC_CH4_NATIVE = GC_CH4_NATIVE_I(LGC:1:-1)
!         ENDIF
!      ENDIF
!
!      ! Interpolate GC CH4 column to GOSAT grid 
!      ! use the method of Parker et al., (2011)
!      TOTAL_SCD  = 0D0
!      SCD(:)     = 0D0 ! SCD at toa is zero
!      DO L = 1, LGC
!         DPRES     = PEDGE(L+1) - PEDGE(L)
!         TEMP_SCD  = ( ( AV_NUM * GC_CH4_NATIVE(L) )
!     &               / ( MOLWT * GRAV ) ) * DPRES
!         SCD(L+1)  = TOTAL_SCD + TEMP_SCD
!         TOTAL_SCD = TOTAL_SCD + TEMP_SCD
!       ENDDO
!
!      ! Determine how much a level influences the levels above and below it
!      DO LL = 1, LGOS
!         
!         ! Make coding cleaner
!         p = GOS_P
!
!         ! Different cases: (a) toa, (b) surface, (c) else
!         IF (LL .EQ. 1) THEN
!            test1 = ABS( -1D0*p(LL) + ( ( p(LL+1) - p(LL) )
!     &              / LOG( p(LL+1) / p(LL) ) ) )
!            test2 = 0D0
!         ELSE IF (LL .EQ. LGOS) THEN
!            test1 = 0D0
!            test2 = ABS( p(LL) - ( ( p(LL) - p(LL-1) )
!     &              / LOG( p(LL) / p(LL-1) ) ) )
!         ELSE
!            test1 = ABS( -1D0*p(LL) + ( ( p(LL+1) - p(LL) )
!     &              / LOG( p(LL+1) / p(LL) ) ) )
!            test2 = ABS( p(LL) - ( ( p(LL) - p(LL-1) )
!     &              / LOG( p(LL) / p(LL-1) ) ) )
!         ENDIF
!
!         ! Get the three layers (bottom, top, and center)
!         LAYER_BOT(LL) = ( GOS_P(LL) + test1 )
!         LAYER_CEN(LL) = ( GOS_P(LL)         )
!         LAYER_TOP(LL) = ( GOS_P(LL) - test2 )
!
!         ! Interpolate to the layers
!         temp1 = INTERPOL( SCD, PEDGE, LAYER_BOT(LL), LGC+1, 1 )
!         temp2 = INTERPOL( SCD, PEDGE, LAYER_TOP(LL), LGC+1, 1 )
!
!         ! Get the interpolated VMR
!         PAR_COL(LL) = (temp1(1) - temp2(1))
!         DPRES       = ( LAYER_BOT(LL) - LAYER_TOP(LL))
!         GC_CH4(LL)  = ( (PAR_COL(LL) / DPRES) * (MOLWT * GRAV)
!     &               / AV_NUM )
!
!      ENDDO
!
!      ! Keep the original orientation
!      IF (LGOS .GT. 1) THEN
!         IF (GOS_P_I(2) .LT. GOS_P_I(1)) THEN
!            GC_CH4 = GC_CH4(LGOS:1:-1)
!         ENDIF
!      ENDIF
!
!      ! Store the output
!      INTERP_GC_GOS = GC_CH4
!
!      ! Return
!      RETURN
!
!      ! Return to calling program
!      END FUNCTION INTERP_GC_GOS
!
!!------------------------------------------------------------------------------
!
!      FUNCTION GOS_PEDGE( LGOS, GOS_P_I )
!!
!!******************************************************************************
!!  Function GOS_PEDGE get the pressure edges for the GOSAT pressure grid
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) LGOS      : Number of GOSAT pressure levels           [unit]
!!  (2 ) GOS_P_I   : GOSAT pressure levels                     [mbar]
!!     
!!  Arguments as Output:
!!  ============================================================================
!!  (1 ) GOS_PEDGE : GOSAT pressure edges                      [v/v]
!!     
!!  NOTES:
!!  (1 ) N/A
!!
!!******************************************************************************
!!
!      ! Arguments
!      INTEGER,  INTENT(IN)           :: LGOS
!      REAL(fp), INTENT(IN)             :: GOS_P_I(LGOS)
!
!      ! Return value 
!      REAL(fp)                         :: GOS_PEDGE(LGOS+1)
!
!      ! Local variables 
!      INTEGER                        :: LL
!      REAL(fp)                         :: p(LGOS)
!      REAL(fp)                         :: PEDGE(LGOS+1)
!      REAL(fp)                         :: LAYER_BOT(LGOS)
!      REAL(fp)                         :: LAYER_TOP(LGOS)
!      REAL(fp)                         :: GOS_P(LGOS)
!      REAL(fp)                         :: test1, test2
!
!      !=================================================================
!      ! INTER_GC_GOS begins here!
!      !=================================================================
!
!      ! Reverse the array orders so toa is first
!      GOS_P = GOS_P_I
!      IF (LGOS .GT. 1) THEN
!         IF (GOS_P_I(2) .LT. GOS_P_I(1)) GOS_P = GOS_P(LGOS:1:-1)
!      ENDIF
!
!      ! Determine how much a level influences the levels above and below it
!      DO LL = 1, LGOS
!         
!         ! Make coding cleaner
!         p = GOS_P
!
!         ! Different cases: (a) toa, (b) surface, (c) else
!         IF (LL .EQ. 1) THEN
!            test1 = ABS( -1D0*p(LL) + ( ( p(LL+1) - p(LL) )
!     &              / LOG( p(LL+1) / p(LL) ) ) )
!            test2 = 0D0
!         ELSE IF (LL .EQ. LGOS) THEN
!            test1 = 0D0
!            test2 = ABS( p(LL) - ( ( p(LL) - p(LL-1) )
!     &              / LOG( p(LL) / p(LL-1) ) ) )
!         ELSE
!            test1 = ABS( -1D0*p(LL) + ( ( p(LL+1) - p(LL) )
!     &              / LOG( p(LL+1) / p(LL) ) ) )
!            test2 = ABS( p(LL) - ( ( p(LL) - p(LL-1) )
!     &              / LOG( p(LL) / p(LL-1) ) ) )
!         ENDIF
!
!         ! Get the three layers (bottom, top, and center)
!         LAYER_BOT(LL) = ( GOS_P(LL) + test1 )
!         LAYER_TOP(LL) = ( GOS_P(LL) - test2 )
!      ENDDO
!      
!      ! Get the edges
!      DO LL = 1, LGOS+1
!
!         IF (LL .EQ. LGOS+1) THEN
!            PEDGE(LL) = LAYER_BOT(LL-1)
!         ELSE
!            PEDGE(LL) = LAYER_TOP(LL)
!         ENDIF
!
!      ENDDO
!
!      ! Return to the original grid formatting
!      IF (LGOS .GT. 1) THEN
!         IF (GOS_P_I(2) .LT. GOS_P_I(1)) PEDGE = PEDGE(LGOS+1:1:-1)
!      ENDIF
!
!      ! Store the output
!      GOS_PEDGE = PEDGE
!
!      ! Return
!      RETURN
!
!
!      ! Return to calling program
!      END FUNCTION GOS_PEDGE
!
!!------------------------------------------------------------------------------
!
!      FUNCTION INTERPOL( Y, X, U, Na, Nb )
!!
!!******************************************************************************
!!  Function INTERPOL does simple linear interpolation
!!
!!  Arguments as Input:
!!  ============================================================================
!!  (1 ) Y  : Input vector                                     [unit]
!!  (2 ) X  : Abscissa values for "Y"                          [unit]
!!  (3 ) U  : Abscissa values to interpolate onto              [unit]
!!  (4 ) Na : Length of "Y" and "X"                            [unit]
!!  (5 ) Nb : Length of "U"                                    [unit]
!!     
!!  Arguments as Output:
!!  ============================================================================
!!  (1 ) INTERPOL : Interpolated data                          [unit]
!!     
!!  NOTES:
!!  (1 ) N/A
!!
!!******************************************************************************
!!
!      ! Arguments
!      INTEGER,  INTENT(IN)           :: Na
!      INTEGER,  INTENT(IN)           :: Nb
!      REAL(fp),  INTENT(IN)            :: Y(Na)
!      REAL(fp),  INTENT(IN)            :: X(Na)
!      REAL(fp),  INTENT(IN)            :: U(Nb)
! 
!      ! Return value 
!      REAL(fp)                         :: INTERPOL(Nb)
!
!      ! Local variables 
!      REAL(fp)   :: X_S, X_L, X_H
!      REAL(fp)   :: Y_S, Y_L, Y_H
!      INTEGER  :: I, J
!
!      !=================================================================
!      ! INTERPOL begins here!
!      !=================================================================
!
!      INTERPOL(:) = 0D0 
!
!      ! Loop over the variable we're interpolating to
!      DO I = 1, Nb
!
!         ! The x position we want
!         X_S = U(I)
!
!         ! Is the data in forward or reverse order?
!         IF ( X(2) .GE. X(1) ) THEN
!            
!            ! Checks to see if we're outside of the "X" data
!            IF ( X_S .LT. X(1) ) THEN
!               X_L = X(1)
!               X_H = X(2)
!               Y_L = Y(1)
!               Y_H = Y(2)
!            ELSE IF ( X_S .GT. X(Na) ) THEN
!               X_L = X(Na-1)
!               X_H = X(Na)
!               Y_L = Y(Na-1)
!               Y_H = Y(Na)
!            ELSE
!               ! Loop over the original data to find brackets
!               DO J = 1, Na-1
!                  IF ( (X(J) .LE. X_S) .AND. (X(J+1) .GT. X_S) ) THEN
!                     X_L = X(J)
!                     X_H = X(J+1)
!                     Y_L = Y(J)
!                     Y_H = Y(J+1)
!                  ENDIF
!               ENDDO
!            ENDIF
!         
!         ELSE ! Data is in reverse order
!
!            ! Checks to see if we're outside of the "X" data
!            IF ( X_S .LT. X(Na) ) THEN
!               X_L = X(Na)
!               X_H = X(Na-1)
!               Y_L = Y(Na)
!               Y_H = Y(Na-1)
!            ELSE IF ( X_S .GT. X(1) ) THEN
!               X_L = X(2)
!               X_H = X(1)
!               Y_L = Y(2)
!               Y_H = Y(1)
!            ELSE
!               ! Loop over the original data to find brackets
!               DO J = 1, Na-1
!                  IF ( (X(J) .GT. X_S) .AND. (X(J+1) .LE. X_S) ) THEN
!                     X_L = X(J+1)
!                     X_H = X(J)
!                     Y_L = Y(J+1)
!                     Y_H = Y(J)
!                  ENDIF
!               ENDDO
!            ENDIF
!
!         ENDIF
!         
!         ! Interpolate
!         Y_S = ( (Y_H - Y_L) / (X_H - X_L) ) * (X_S - X_L) + Y_L
!
!         ! Save the output
!         INTERPOL(I) = Y_S
!
!      ENDDO
!
!      ! Return the data
!      RETURN
!
!
!      ! Return to calling program
!      END FUNCTION INTERPOL
!
!==============================================================================
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
      SUBROUTINE GET_INTMAP( GCPCEN, GCPSURF, GOSPEDGE, nlev, INTMAP )
!
! !USES:
!
      USE CMN_SIZE_MOD
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN)  :: GCPCEN(LLPAR)
      REAL(fp), INTENT(IN)  :: GCPSURF
      REAL(fp), INTENT(IN)  :: GOSPEDGE(MAXLEV) 
      INTEGER,  INTENT(IN)  :: nlev
!
! !OUTPUT PARAMETERS:
!
      REAL(fp), INTENT(OUT) :: INTMAP(LLPAR,MAXLEV)
!
! !REVISION HISTORY:
!  16 Jun 2017 - M. Sulprizio- Initial version based on GOSAT CH4 observation
!                              operator from GC Adjoint v35j
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
      DO LTM = 1, nlev

         ! Find the levels from GC that bracket level LTM
         DO LGC = 1, LLPAR-1

            LOW = GCPCEN(LGC+1)
            HI  = GCPCEN(LGC)

            ! Match GEOS-Chem level to GOS level
            IF ( GOSPEDGE(LTM) <= HI .and. 
     &           GOSPEDGE(LTM)  > LOW) THEN 

               DIFF             = HI - LOW  
               INTMAP(LGC+1,LTM) = ( HI - GOSPEDGE(LTM)  ) / DIFF
               INTMAP(LGC  ,LTM) = ( GOSPEDGE(LTM) - LOW ) / DIFF

            ENDIF

          ENDDO

       ENDDO

       ! Correct for case where GOSAT pressure is higher than the
       ! highest GC pressure center.  In this case, just 1:1 map. 
       DO LTM = 1, nlev
          IF ( GOSPEDGE(LTM) > GCPCEN(1) ) THEN
             INTMAP(:,LTM) = 0e+0_fp
             INTMAP(1,LTM) = 1e+0_fp
          ENDIF
       ENDDO

      END SUBROUTINE GET_INTMAP
!EOC
      END MODULE GOSAT_CH4_MOD
