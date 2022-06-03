!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: airs_ch4_mod.F90
!
! !DESCRIPTION: Module AIRS\_CH4\_MOD
!\\
!\\
! !INTERFACE: 
!
MODULE AIRS_CH4_MOD
!
! !USES:
!
#if !defined( MODEL_CESM )
  USE m_netcdf_io_open       ! netCDF open
  USE m_netcdf_io_get_dimlen ! netCDF dimension queries
  USE m_netcdf_io_read       ! netCDF data reads
  USE m_netcdf_io_close      ! netCDF close
#endif
  USE PRECISION_MOD          ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE 
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CALC_AIRS_CH4_FORCE
!
! !REVISION HISTORY:
!  20 Sept 2018 - Yuzhong Zhang - Initial version based on GOSAT_CH4_MOD.F
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  INTEGER,  PARAMETER  :: MAXLEV = 26
  INTEGER,  PARAMETER  :: MAXOBS = 200000
  LOGICAL,  PARAMETER  :: LDCH4SAT    = .TRUE.
!
! !MODULE VARIABLES
!
  ! Record to store data from each AIRS obs
  TYPE AIRS_CH4_OBS 
     INTEGER           :: LAIRS(1)
     REAL(fp)          :: LAT(1)
     REAL(fp)          :: LON(1)
     INTEGER           :: YEAR(1)
     INTEGER           :: MONTH(1)
     INTEGER           :: DAY(1)
     INTEGER           :: HOUR(1)
     INTEGER           :: MINUTE(1)
     INTEGER           :: SECOND(1)
     REAL(fp)          :: CH4(1)   ! column v/v
     REAL(fp)          :: CH4_ERROR(1) ! column v/v
     REAL(fp)          :: PRES(MAXLEV)  ! hPa
     REAL(fp)          :: PRIOR(MAXLEV) ! v/v
     REAL(fp)          :: AVG_KERNEL(MAXLEV) ! column vmr/vmr 
     REAL(fp)          :: P_WEIGHT(MAXLEV) ! pressure weight 
     REAL(fp)          :: DOF(1) !trace of avg_kernel
  ENDTYPE AIRS_CH4_OBS

  TYPE(AIRS_CH4_OBS)    :: AIRS(MAXOBS)

CONTAINS
!EOC
#if !defined( MODEL_CESM )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: read_airs_ch4_obs
!
! !DESCRIPTION: Subroutine READ\_AIRS\_CH4\_OBS reads the file and passes back
!  info contained therein. (dkh, 10/12/10) 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READ_AIRS_CH4_OBS( YYYYMMDD, NAIRS )
!
! !USES:
!
!
    USE TIME_MOD,             ONLY : EXPAND_DATE
    USE ERROR_MOD,            ONLY : ALLOC_ERR
    USE ERROR_MOD,            ONLY : GEOS_CHEM_STOP
    USE TIME_MOD,             ONLY : GET_YEAR, YMD_EXTRACT
 
#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)      :: YYYYMMDD   ! Date
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT)     :: NAIRS      ! Number of AIRS retrievals
!     
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                        :: NLEV
    INTEGER                        :: NTIME
    ! For reading netCDF file
    INTEGER            :: fId                ! netCDF file ID
    INTEGER            :: Status             ! 0 means variable in file
    INTEGER            :: X, Y, Z, T         ! netCDF file dimensions
    INTEGER            :: time_index         ! Read this slice of data
    INTEGER            :: st1d(1), ct1d(1)   ! Start + count for 1D arrays
    INTEGER            :: st2d(2), ct2d(2)   ! Start + count for 2D arrays
    INTEGER            :: RC                 ! 0 means dimension exists
    CHARACTER(LEN=16)  :: stamp              ! Time and date stamp
    CHARACTER(LEN=255) :: dir                ! Data directory path
    CHARACTER(LEN=255) :: nc_file            ! netCDF file name
    CHARACTER(LEN=255) :: v_name             ! netCDF variable name
    CHARACTER(LEN=255) :: errMsg             ! Error message
    CHARACTER(LEN=255) :: caller             ! Name of this routine
    CHARACTER(LEN=4)   :: CYEAR              ! Year in character format
    CHARACTER(LEN=255) :: DimName1, DimName2 ! netCDF dimension names
    CHARACTER(LEN=255) :: DimName3           ! netCDF dimension names

    INTEGER                :: LAIRS
    LOGICAL                :: EXIST_VAR
    REAL(fp), ALLOCATABLE  :: lat(:)
    REAL(fp), ALLOCATABLE  :: lon(:)
    REAL(fp), ALLOCATABLE  :: ch4(:)
    REAL(fp), ALLOCATABLE  :: ch4_error(:)
    REAL(fp), ALLOCATABLE  :: pres(:,:)
    REAL(fp), ALLOCATABLE  :: prior(:,:)
    REAL(fp), ALLOCATABLE  :: time(:,:)
    REAL(fp), ALLOCATABLE  :: ak(:,:)
    REAL(fp), ALLOCATABLE  :: pres_w(:,:)

    INTEGER                :: I, L, AS, N
    REAL(fp)               :: REF_DATE, TIME_JD

    !=================================================================
    ! READ_AIRS_CH4_OBS begins here!
    !=================================================================

    caller = 'READ_AIRS_CH4_OBS in airs_ch4_mod.F90'

    ! Get current year
    WRITE( CYEAR, '(i4)' ) GET_YEAR()

    ! Filename
    nc_file = 'TROPESS_AIRS_v1.6_CH4_YYYYMM.nc'
    CALL EXPAND_DATE( nc_file, YYYYMMDD, 9999 ) 

    ! Construct complete file path
    dir = '/n/seasasfs02/epenn/AIRS/'
    nc_file = TRIM( dir ) // TRIM( nc_file )
    WRITE( 6, 10 ) TRIM( nc_file )
10  FORMAT( '     - Reading ', a)

    ! Make sure the file exists (ajt, 03/31/2013)
    INQUIRE( FILE=TRIM( nc_file ), EXIST=EXIST_VAR )
    IF ( .NOT. EXIST_VAR ) THEN
       NAIRS = -1
       RETURN
    ENDIF

    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Get name of dimension 1 (Observation)
    RC = NF_Inq_DimName( fId, 1, DimName1 )
    IF ( RC /= 0 ) THEN
       PRINT*, 'Could not find name for dimension 1...'
       CALL GEOS_CHEM_STOP
    ENDIF

    ! Get name of dimension 2 (Level)
    RC = NF_Inq_DimName( fId, 2, DimName2 )
    IF ( RC /= 0 ) THEN
       PRINT*, 'Could not find name for dimension 2...'
       CALL GEOS_CHEM_STOP
    ENDIF

    ! Get name of dimension 3 (UTC\ Date)
    RC = NF_Inq_DimName( fId, 3, DimName3 )
    IF ( RC /= 0 ) THEN
       PRINT*, 'Could not find name for dimension 3...'
       CALL GEOS_CHEM_STOP
    ENDIF

    ! Get length of dimensions
    CALL NcGet_DimLen( fId, DimName1,   NLEV )
    CALL NcGet_DimLen( fId, DimName2,   NAIRS )
    CALL NcGet_DimLen( fId, DimName3,   NTIME )

    print*,'DimName3: ',DimName3
    print*,'length: ',NTIME

    IF ( NLEV > MAXLEV ) THEN
       print*,' # Levels this day = ', NLEV
       print*, 'WARNING: NLEV > MAXLEV. Need to increase'
       print*, ' MAXLEV in airs_ch4_mod.f.'
       CALL GEOS_CHEM_STOP
    ENDIF

    IF ( NAIRS > MAXOBS ) THEN
       print*,' # of observation in this month = ', NAIRS
       print*, 'WARNING: NAIRS > MAXOBS. Need to increase'
       print*, ' MAXOBS in airs_ch4_mod.f.'
       CALL GEOS_CHEM_STOP
    ENDIF

    print*,' # AIRS Observations this month = ', NAIRS
    print*, 'levels', NLEV

    !--------------------------------
    ! Allocate arrays for data to be read in
    !--------------------------------

    ALLOCATE( lat(         NAIRS ), STAT=AS )
    ALLOCATE( lon(         NAIRS ), STAT=AS )
    ALLOCATE( ch4(         NAIRS),  STAT=AS )
    ALLOCATE( ch4_error(   NAIRS ), STAT=AS )

    ALLOCATE( prior(     NLEV,  NAIRS ), STAT=AS )
    ALLOCATE( pres(      NLEV,  NAIRS ), STAT=AS )
    ALLOCATE( ak(        NLEV,  NAIRS ), STAT=AS )
    ALLOCATE( pres_w(    NLEV,  NAIRS ), STAT=AS )
    ALLOCATE( time(      NTIME, NAIRS ), STAT=AS )

    !--------------------------------
    ! Read 1-D Data
    !--------------------------------

    ! Latitude
    v_name = 'latitude'
    st1d   = (/ 1    /)
    ct1d   = (/ NAIRS /)
    CALL NcRd( lat, fId, TRIM(v_name), st1d, ct1d )

    ! Longitude
    v_name = 'longitude'
    st1d   = (/ 1    /)
    ct1d   = (/ NAIRS /)
    CALL NcRd( lon, fId, TRIM(v_name), st1d, ct1d )

    ! CH4 (ppb)
    v_name = 'xcol_ft'
    st1d   = (/ 1   /)
    ct1d   = (/ NAIRS /)
    CALL NcRd( ch4, fId, TRIM(v_name), st1d, ct1d )

    ! ch4 column error
    v_name = 'xcol_err_ft'
    st1d   = (/ 1   /)
    ct1d   = (/ NAIRS /)
    CALL NcRd( ch4_error, fId, TRIM(v_name), st1d, ct1d )

    !-------------------------------- 
    ! Read 2D Data
    !-------------------------------- 

    ! APRIORI (ppb)
    v_name = 'xa_ch4'
    st2d   = (/ 1,    1    /)
    ct2d   = (/ NLEV, NAIRS /)
    CALL NcRd( prior, fId, TRIM(v_name), st2d, ct2d )

    ! Pressure (hPa)
    v_name = 'pressure'
    st2d   = (/ 1,    1    /)
    ct2d   = (/ NLEV, NAIRS /)
    CALL NcRd( pres, fId, TRIM(v_name), st2d, ct2d )

    ! Averaging Kernel (linearized & multiplied with pres_w)
    v_name = 'ak_col_ft'
    st2d   = (/ 1,    1    /)
    ct2d   = (/ NLEV, NAIRS /)
    CALL NcRd( ak, fId, TRIM(v_name), st2d, ct2d )

    ! Pressure Weight 
    v_name = 'col_ft'
    st2d   = (/ 1,    1    /)
    ct2d   = (/ NLEV, NAIRS /)
    CALL NcRd( pres_w, fId, TRIM(v_name), st2d, ct2d )

    ! UTC date and time 
    v_name = 'UTC'
    st2d   = (/ 1,   1    /)
    ct2d   = (/ NTIME, NAIRS /)
    CALL NcRd( time, fId, TRIM(v_name), st2d, ct2d )

    !-------------------------------- 
    ! Place data into AIRS structure
    !-------------------------------- 
    DO N = 1, NAIRS

       ! 1-D data
       AIRS(N)%LAIRS(1)        = NLEV
       AIRS(N)%LAT(1)          = lat(N)
       AIRS(N)%LON(1)          = lon(N)
       AIRS(N)%CH4(1)          = ch4(N)
       AIRS(N)%CH4_ERROR(1)    = ch4_error(N)
         
       ! 2-D data
       LAIRS = NLEV
       AIRS(N)%PRIOR(1:LAIRS)      = prior(1:LAIRS, N)
       AIRS(N)%PRES(1:LAIRS)       = pres(1:LAIRS, N)
       AIRS(N)%AVG_KERNEL(1:LAIRS) = ak(1:LAIRS, N)
       AIRS(N)%P_WEIGHT(1:LAIRS)   = pres_w(1:LAIRS, N)
         
       ! derived date var
       ! time array is length 7 with format (mm, dd, YYYY, HH, MM, SS, 0)
       ! split these into separate variables and convert 
       ! from floats to integers (erp, 07/02/2020)
       AIRS(N)%MONTH(1)  = NINT(time(1,N))
       AIRS(N)%DAY(1)    = NINT(time(2,N))
       AIRS(N)%YEAR(1)   = NINT(time(3,N))
       AIRS(N)%HOUR(1)   = NINT(time(4,N))
       AIRS(N)%MINUTE(1) = NINT(time(5,N))
       AIRS(N)%SECOND(1) = NINT(time(6,N))
 
    ENDDO

    !-------------------------------- 
    ! Close netCDF file
    !-------------------------------- 

    ! Echo info
    WRITE( 6, 20 ) YYYYMMDD
20  FORMAT( '     - Finished reading AIRS CH4 observations for ', i8)

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE READ_AIRS_CH4_OBS
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_airs_ch4_force
!
! !DESCRIPTION: Subroutine CALC\_AIRS\_CH4\_FORCE calculates the adjoint forcing
!  from the GOSAT CH4 observations and updates the cost function.
!  (dkh, 10/12/10)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_AIRS_CH4_FORCE( Input_Opt, State_Chm, State_Grid, &
                                  State_Met )
!
! !USES:
!
    USE ErrCode_Mod
    USE GC_GRID_MOD,        ONLY : GET_IJ
    USE Input_Opt_Mod,      ONLY : OptInput
    USE TIME_MOD
    USE PhysConstants,      ONLY : XNUMOLAIR, AIRMW
    USE State_Chm_Mod,      ONLY : ChmState, Ind_
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input options
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
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
    INTEGER            :: id_CH4
    INTEGER,  SAVE     :: NAIRS
    INTEGER            :: IIJJ(2), I,      J,  NT
    INTEGER            :: L,       LL
    INTEGER            :: NOBS,   IND
    INTEGER            :: L0, L1
    LOGICAL            :: IFVALID
    INTEGER            :: INDS(MAXOBS)
    INTEGER            :: YYYYMMDD
    REAL(fp)           :: GC_PRES(State_Grid%NZ)
    REAL(fp)           :: GC_PEDGE(State_Grid%NZ+1)
    REAL(fp)           :: DRY_AIR(State_Grid%NZ)
    REAL(fp)           :: WATER_VAPOR(State_Grid%NZ)
    REAL(fp)           :: GC_CH4_NATIVE(State_Grid%NZ)
    REAL(fp)           :: GC_CH4(MAXLEV)
    REAL(fp)           :: GC_CH4_ORIG(MAXLEV)
    REAL(fp)           :: GC_PSURF
    REAL(fp)           :: MAP(State_Grid%NZ,MAXLEV)
    REAL(fp)           :: p(MAXLEV)
    REAL(fp)           :: h(MAXLEV)
    REAL(fp)           :: ch4
    REAL(fp)           :: ak1(MAXLEV), ak2(MAXLEV)
    REAL(fp)           :: pres_w(MAXLEV)
    REAL(fp)           :: prior(MAXLEV)
    REAL(fp)           :: TropP !tropopause pressure (hPa)
    INTEGER            :: LTrop !tropopause layer in GOSAT levels
    REAL(fp)           :: WT_LTrop
    REAL(fp)           :: GC_XCH4, GC_XCH4_ORIG, GC_XCH4_ORIG_oldpw
    REAL(fp)           :: GC_XCH4_prior, GC_XCH4_trop, GC_XCH4_strat

    ! For miscellaneous
    LOGICAL, SAVE      :: FIRST = .TRUE. 
    INTEGER            :: IOS
    INTEGER, SAVE      :: TotalObs = 0
    CHARACTER(LEN=255) :: FILENAME
    CHARACTER(LEN=63)  :: OrigUnit
    CHARACTER(LEN=255) :: ThisLoc
    CHARACTER(LEN=512) :: ErrMsg
    INTEGER            :: RC
    REAL(fp)           :: foo ! throwaway output

    !=================================================================
    ! CALC_AIRS_CH4_FORCE begins here!
    !=================================================================

    print*, '     - CALC_AIRS_CH4_FORCE '

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at CALC_AIRS_CH4_FORCE (in gosat_ch4_mod.F)'

    ! Initialize species ID flag
    id_CH4     = Ind_('CH4'       )

    ! Open files for diagnostic output
    IF ( FIRST ) THEN
       IF ( FIRST ) FIRST = .FALSE.

       FILENAME = 'sat_obs.airs.00.m'
       FILENAME = TRIM( Input_Opt%RUN_DIR ) //  TRIM( FILENAME )
       OPEN( 118,      FILE=TRIM( FILENAME   ), STATUS='UNKNOWN', &
            IOSTAT=IOS, FORM='FORMATTED',    ACCESS='SEQUENTIAL' )

       ! Write header of sat_obs.NN.m
       WRITE( 118, 281 ) '       NNN', &
                         '   I', '   J', '     LON','     LAT','YYYY', &
                         'MM', 'DD', 'hh', 'mm', 'ss',                 &
                         '         TAU', '        AIRS',               &
                         '       S_OBS',  '       model',              &
                         '  model_orig',  'model_oldpw',               & !erp
                         '       TropP', ' model_prior',               & !zyz
                         '  model_trop', ' model_strat'                  !zyz
281    FORMAT( A10,2x,A4,2x,A4,2x,A8,2x,A8,2x,A4,2x,A2,2x,A2,2x,A2,2x, &
               A2,2x,A2,2x,A12,2x,A12,2x,A12,2x,A12,2x, A12,2x,        &
               A12,2x, A12,2x, A12, 2x, A12, 2x, A12, 2x)

       ! Set Total Observations = 0
       TotalObs = 0

    ENDIF ! FIRST

    ! Read Observations at first call and at end of the day
    IF ( FIRST .OR. ITS_A_NEW_DAY() ) THEN
 
       ! Read the AIRS CH4 file for this month
       YYYYMMDD = 1d4*GET_YEAR() + 1d2*GET_MONTH() + GET_DAY()
       CALL READ_AIRS_CH4_OBS( YYYYMMDD, NAIRS )

       ! Make sure there are observations on this day
       IF ( NAIRS .EQ. -1 ) RETURN

    ENDIF

    ! Get indices of AIRS observations in the current hour
    !   At the start of each hour, assimilate observations that 
    !   were made in the previous 60 minutes.
    !   For example, at time 18:00, assimilate observations 
    !   made from 18:00 - 18:59
    INDS(:) = 0
    NOBS    = 0
    !print*,'Looking for observations at MONTH, DAY, HOUR = ', &
    !        GET_MONTH(), GET_DAY(), GET_HOUR()

    DO NT = 1, NAIRS
       IF ( AIRS(NT)%MONTH(1) .EQ. GET_MONTH() .AND. &
            AIRS(NT)%DAY(1)   .EQ. GET_DAY()   .AND. &
            AIRS(NT)%HOUR(1)  .EQ. GET_HOUR()  ) THEN
          NOBS = NOBS + 1
          INDS(NOBS) = NT
          !print*,'Found a good observation! NT = ', NT
       ENDIF
    ENDDO

    IF ( NOBS == 0 ) THEN
       print*, ' No matching AIRS CH4 obs for this hour'
       RETURN
    ENDIF
    print*, ' for day ',GET_YEAR(), GET_MONTH(), GET_DAY()
    print*, ' for hour range: ', GET_HOUR(), GET_HOUR()+1
    print*, ' found # AIRS observations: ', NOBS

    ! Convert species units to [v/v] (mps, 6/12/2020)
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'v/v dry', RC,        OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error (kg/kg dry -> v/v dry)'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

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

       ! Check for layers with valid retrieved value in data
       CALL GET_VALID_LAYERS (NT, L0, L1, IFVALID)
       IF ( .NOT. IFVALID ) THEN
          print*, 'This retrieval appears to be invalid', NT
          CYCLE
       ENDIF

       ! Skip Observations outside the domain
       IF ( AIRS(NT)%LAT(1) < State_Grid%YMin .OR. &
            AIRS(NT)%LAT(1) > State_Grid%YMax .OR. &
            AIRS(NT)%LON(1) < State_Grid%XMin .OR. &
            AIRS(NT)%LON(1) > State_Grid%XMax ) THEN
          print*, ' Outside nested domain, skipping record ', NT
          CYCLE
       ENDIF

       ! Get grid box of current record
       IIJJ = GET_IJ( REAL(AIRS(NT)%LON(1),4), &
                      REAL(AIRS(NT)%LAT(1),4), &
                      State_Grid )
       I    = IIJJ(1)
       J    = IIJJ(2)

       ! skip observations where the AIRS surface  is much
       ! lower than the model
       IF ( (AIRS(NT)%PRES(L0) - State_Met%PEDGE(I,J,1)) > 50e0 ) THEN
          print*, ' Psurf threshold not met, skipping record ', NT
          CYCLE
       ENDIF

       !------------------------------
       ! Begin good observations
       !------------------------------
       print*,' Begin assimilating good observation. NT = ', NT

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
         
       ! Get CH4 values at native model resolution [v/v]
       GC_CH4_NATIVE(:) = 0.0_fp 
       GC_CH4_NATIVE(:) = State_Chm%Species(id_CH4)%Conc(I,J,:)

       ! Get # mols of dry air in each gridbox
       DRY_AIR(:) = 0.0_fp
       DRY_AIR(:) =  State_Met%AIRDEN(I,J,:) * XNUMOLAIR * &
                     1e-6_fp * State_Met%BXHEIGHT(I,J,:)*1e2_fp

       WATER_VAPOR(:) = 0.0_fp
       WATER_VAPOR(:) = State_Met%AVGW(I,J,:) 

       ! Use short names for clarity
       p(:)      = AIRS(NT)%PRES(:)
       ak1(:)    = AIRS(NT)%AVG_KERNEL(:)
       prior(:)  = AIRS(NT)%PRIOR(:)
       pres_w(:) = AIRS(NT)%P_WEIGHT(:)

       ! Compute the AIRS layer that encloses tropopause
       TropP = State_Met%TROPP(I,J)
       CALL GET_TROP_LAYER(TropP, p, L0, L1, LTrop, WT_LTrop)

       ! Convert from AIRS AK format to GOSAT AK format
       CALL CONVERT_COL_AK(ak1, pres_w, L0, L1, ak2)

       ! *****************************
       !  SIMPLE LINEAR INTERPOLATION
       ! *****************************
       !   simple interpolation does not consider mass or dry air
       !   it may result in loss of information if satellite grid 
       !   is lower resolution than model (see Rodgers 2000, and 
       !   Keppens et al. 2019)
       !   this is what JDM, ZYZ, and AJT use 
       ! Calculate the interpolation weight matrix 
       MAP(:,:) = 0.0_fp
       CALL GET_INTMAP_OLD( State_Grid, GC_PRES, GC_PSURF, &
                            AIRS(NT)%PRES, L0, L1, MAP )

       ! Interpolate GC CH4 column to AIRS grid
       GC_CH4_ORIG(:) = 0.0_fp
       DO LL = L0, L1
          GC_CH4_ORIG(LL) = 0.0_fp
          DO L = 1, State_Grid%NZ 
             GC_CH4_ORIG(LL) = GC_CH4_ORIG(LL) &
                             + MAP(L,LL) * GC_CH4_NATIVE(L) 
          ENDDO
       ENDDO

       ! Compute the GEOS-Chem XCH4 corresponding to the observation
       ! With AIRS pressure weighting
       CALL CALC_GC_XCH4 (GC_CH4_ORIG, ak2, pres_w, prior, L0, L1, &
                          LTROP, WT_LTROP, &
                          GC_XCH4_ORIG, foo, & 
                          foo, foo)
         
       ! Calculate AJT and JDM pressure weighting
       h = 0.0_fp
       CALL CALC_PRES_WEIGHT_AJT (p, L0, L1, h)

       ! Compute the GEOS-Chem XCH4 corresponding to the observation
       ! Using AJT & JDM's pressure weighting
       CALL CALC_GC_XCH4 (GC_CH4_ORIG, ak2, h, prior, L0, L1, &
                          LTROP, WT_LTROP, &
                          GC_XCH4_ORIG_oldpw, foo, & 
                          foo, foo)

       ! **************************
       !  MASS-BASED INTERPOLATION
       ! **************************
       !   Interpolate using new (mass-based) interp 
       !   that gets concentration on pressure **edges** 
       !   Function based on Keppens et al. 2019
       !   erp, Oct 6, 2020

       CALL MASS_INTERP( GC_PEDGE, AIRS(NT)%PRES, GC_CH4_NATIVE, &
                         L0, State_Grid%NZ, AIRS(NT)%LAIRS(1),   &
                         GC_CH4 )

       ! Separate the impact prior, troposphere, and stratosphere
       CALL CALC_GC_XCH4 (GC_CH4, ak2, pres_w, prior, L0, L1, &
                          LTROP, WT_LTROP, &
                          GC_XCH4, GC_XCH4_prior, &
                          GC_XCH4_trop, GC_XCH4_strat)

       TotalObs = TotalObs + 1
       ! Record information for satellite diagnostics
       IF ( LDCH4SAT ) THEN 
          WRITE( 118, 283 ) TotalObs, I, J, AIRS(NT)%LON(1), &
             AIRS(NT)%LAT(1),AIRS(NT)%YEAR(1), &
             AIRS(NT)%MONTH(1), AIRS(NT)%DAY(1), AIRS(NT)%HOUR(1), &
             AIRS(NT)%MINUTE(1), AIRS(NT)%SECOND(1), GET_TAU(), &
             AIRS(NT)%CH4(1), AIRS(NT)%CH4_ERROR(1), GC_XCH4, &
             GC_XCH4_ORIG, GC_XCH4_ORIG_oldpw, TROPP, &
             GC_XCH4_prior, GC_XCH4_trop, GC_XCH4_strat
       ENDIF

    ENDDO  ! NT
    !!$OMP END PARALLEL DO

    ! Convert species units back to original unit (mps, 6/12/2020)
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit,  RC ) 
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

283 FORMAT( I10,2x,I4,2x,I4,2x,F8.3,2x,F8.4,2x,I4,2x,I2,2x,I2,2x,I2, &
            2x,I2,2x,I2,2x,F12.3,2x,E12.6,2x,E12.6,2x,E12.6, &
            2x, E12.6,2x, E12.6,2x,  F12.3, 2x,E12.6, 2x, E12.6, 2x, &
            E12.6,2x ) 
    print*, ' Number of observations this hour = ', NOBS
    print*, ' Number of observations total     = ', TotalObs 


  END SUBROUTINE CALC_AIRS_CH4_FORCE
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
  SUBROUTINE GET_INTMAP_OLD( State_Grid, GCPCEN, GCPSURF, AIRSPEDGE, &
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
    REAL(fp),       INTENT(IN)  :: AIRSPEDGE(MAXLEV) 
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

    ! Loop over each pressure level of AIRS grid
    DO LTM = L0, L1

       ! Find the levels from GC that bracket level LTM
       DO LGC = 1, State_Grid%NZ-1

          LOW = GCPCEN(LGC+1)
          HI  = GCPCEN(LGC)

          ! Match GEOS-Chem level to AIRS level
          IF ( AIRSPEDGE(LTM) <= HI .and. &
               AIRSPEDGE(LTM)  > LOW) THEN 

             DIFF             = HI - LOW  
             INTMAP(LGC+1,LTM) = ( HI - AIRSPEDGE(LTM)  ) / DIFF
             INTMAP(LGC  ,LTM) = ( AIRSPEDGE(LTM) - LOW ) / DIFF

          ENDIF

       ENDDO

    ENDDO

    ! zyz- Need to check what this means?
    ! Correct for case where AIRS pressure is higher than the
    ! highest GC pressure center.  In this case, just 1:1 map. 
    DO LTM = L0, L1
       IF ( AIRSPEDGE(LTM) > GCPCEN(1) ) THEN
          INTMAP(:,LTM) = 0e+0_fp
          INTMAP(1,LTM) = 1e+0_fp
       ENDIF
    ENDDO

  END SUBROUTINE GET_INTMAP_OLD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: get_trop_layer
!
! !DESCRIPTION: Find the AIRS layer that encloses tropopause and compute
!  the fraction of this layer in the tropospher
!\\
!\\
! !INTERFACE: 
!
  SUBROUTINE GET_TROP_LAYER(TROPP, AIRSP, L0, L1, LTROP, WT_LTROP)
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  :: TROPP
    REAL(fp), INTENT(IN)  :: AIRSP(MAXLEV) 
    INTEGER,  INTENT(IN)  :: L0, L1
!
! !OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(OUT) :: LTROP     !Layer encloses tropopause
    REAL(fp), INTENT(OUT) :: WT_LTROP  !Fraction in troposphere
!     
! !REVISION HISTORY:
!  21 Sept 2018 - Yuzhong Zhang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: L
    REAL(fp)              :: AIRSP_HF(MAXLEV+1)

    !Estimate the edge of AIRS layers
    AIRSP_HF(L0) = AIRSP(L0) + (AIRSP(L0)-AIRSP(L0+1))/2.0_fp
    DO L = L0, L1-1
       AIRSP_HF(L+1) = (AIRSP(L)+AIRSP(L+1))/2e+0_fp
    ENDDO
    AIRSP_HF(L1+1) = AIRSP(L1) - (AIRSP(L1-1)-AIRSP(L1))/2.0_fp

    !Find the layer that enclose the tropopause
    !Estimate the WT_LTROP as the fraction of pressure
    !difference in that layer
    IF (TROPP.GE.AIRSP_HF(L0)) THEN
       LTROP = L0 - 1
       WT_LTROP = 1.0_fp
    ELSE IF (TROPP.LT.AIRSP(L1+1)) THEN
       LTROP = L1 + 1
       WT_LTROP = 0.0_fp
    ELSE
       LTROP = L0
       DO L = L0, L1
          IF (TROPP < AIRSP_HF(L) .and. TROPP>= AIRSP_HF(L+1)) THEN
             LTROP=L
             WT_LTROP= (AIRSP_HF(L)-TROPP)/(AIRSP_HF(L)-AIRSP_HF(L+1))
             EXIT
          ENDIF
       ENDDO
    ENDIF

  END SUBROUTINE GET_TROP_LAYER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: get_valid_layers
!
! !DESCRIPTION: First check if the AIRS profile is valid; Second find the
!  lowest and highest layer with valid number
!\\
!\\
! !INTERFACE: 
! 
  SUBROUTINE GET_VALID_LAYERS (NT, L0, L1, IFVALID)
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: NT
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: L0
    INTEGER, INTENT(OUT) :: L1
    LOGICAL, INTENT(OUT) :: IFVALID
!     
! !REVISION HISTORY:
!  21 Sept 2018 - Yuzhong Zhang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: L

    L0 = 0
    L1 = 0
    IFVALID = .TRUE.

    ! check if CH4 column is missing
    IF (AIRS(NT)%CH4(1) .LT. 0d0) THEN
       IFVALID=.FALSE.
    ENDIF

    ! find the first nonmissing level
    DO L = 1, AIRS(NT)%LAIRS(1) 
       IF( AIRS(NT)%PRIOR(L) .GT. 0d0 .AND. &
           AIRS(NT)%PRES(L)  .GT. 0d0 .AND. &
           AIRS(NT)%P_WEIGHT(L) .GT. 0d0 ) THEN
          L0 = L
          EXIT
       ENDIF
    ENDDO

    ! if there are no levels, obs. is invalid
    IF (L0 .EQ. 0) THEN
       IFVALID=.FALSE.
       RETURN
    ENDIF
     
    ! find the last nonmissing level
    DO L = AIRS(NT)%LAIRS(1), 1, -1
       IF( AIRS(NT)%PRIOR(L) .GT. 0d0 .AND. &
           AIRS(NT)%PRES(L) .GT. 0d0 .AND. &
           AIRS(NT)%P_WEIGHT(L) .GT. 0d0 ) THEN
          L1 = L
          EXIT
       ENDIF
    ENDDO

    ! if there are 3 or less levels, obs. is invalid
    IF ((L1-L0).LE.3) THEN
       IFVALID = .FALSE.
       RETURN
    ENDIF

    ! if there are missing levels in the middle, 
    ! (not adjacent to the top/bottom), then obs. is invalid
    DO L = L0, L1
       IF( AIRS(NT)%PRIOR(L) .LE. 0d0 .AND. &
           AIRS(NT)%PRES(L) .LE. 0d0 .AND. &
           AIRS(NT)%AVG_KERNEL(L) .LT. -998d0 ) THEN
          IFVALID=.FALSE.
          EXIT
       ENDIF
    ENDDO

  END SUBROUTINE GET_VALID_LAYERS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ak_log_to_linear
!
! !DESCRIPTION: Convert averaging kernel to linear
!  See Zhang, L., Intercomparison methods for satellite measurements of
!  atmospheric composition: application to tropospheric ozone from
!  AIRS and OMI, Atmos. Chem. Phys. 10, 4,725-4,739, 2010
!  AK0 in d ln(vmr)/d ln(vmr)
!  AK1 in d vmr/ d vmr
!\\
!\\
! !INTERFACE: 
!
  SUBROUTINE AK_LOG_TO_LINEAR (AK0, PRIOR, L0, L1, AK1)
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN) :: AK0(MAXLEV, MAXLEV)
    REAL(fp), INTENT(IN) :: PRIOR(MAXLEV)
    INTEGER,  INTENT(IN) :: L0
    INTEGER,  INTENT(IN) :: L1
!
! !OUTPUT PARAMETERS:
!  
    REAL(fp), INTENT(OUT) :: AK1(MAXLEV, MAXLEV)
!     
! !REVISION HISTORY:
!  21 Sept 2018 - Yuzhong Zhang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: IL, JL

    AK1(:,:) = -999d0
    DO IL = L0, L1
       DO JL = L0, L1
          AK1(IL,JL) = AK0(IL,JL)*PRIOR(IL)/PRIOR(JL)
       ENDDO
    ENDDO

  END SUBROUTINE AK_LOG_TO_LINEAR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: calc_pres_weight_zyz
!
! !DESCRIPTION: Caluate the prssure weight for each AIRS layer
!\\
!\\
! !INTERFACE: 
! 
  SUBROUTINE CALC_PRES_WEIGHT_ZYZ (PRES, L0, L1, PRES_WT)
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  :: PRES(MAXLEV)
    INTEGER,  INTENT(IN)  :: L0, L1
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: PRES_WT(MAXLEV)
!     
! !REVISION HISTORY:
!  21 Sept 2018 - Yuzhong Zhang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: LL
    REAL(fp) :: TOTWT

    PRES_WT(:) = 0.0_fp
    PRES_WT(L0) = PRES(L0)-PRES(L0+1)

    DO LL = L0+1, L1-1
       PRES_WT(LL) = (PRES(LL-1)-PRES(LL+1))/2.0_fp
    ENDDO

    PRES_WT(L1) = PRES(L1-1)-PRES(L1)
    TOTWT = SUM(PRES_WT)
    PRES_WT(:) = PRES_WT(:) / TOTWT

  END SUBROUTINE CALC_PRES_WEIGHT_ZYZ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: calc_pres_wright_ajt
!
! !DESCRIPTION: Compute GC_XCH4 with the h operator used in AJT and JDM's work
!  This method seems to have two problems:
!  1. It ignores the impact of  vertical variation in specific humidity
!     on the weighting function
!  2. The layer boundary here is at the altidue midpoint
!     The paper that the UL-GOSAT paper cites uses pressure midpoint
!\\
!\\
! !INTERFACE: 
!
  SUBROUTINE CALC_PRES_WEIGHT_AJT (PRES, L0, L1, PRES_WT)
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  :: PRES(MAXLEV)
    INTEGER,  INTENT(IN)  :: L0, L1
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: PRES_WT(MAXLEV)
!     
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: p(MAXLEV), h(MAXLEV)
    INTEGER  :: L

    PRES_WT(:) = 0.0_fp
    h = 0.0_fp
    p = 0.0_fp
    ! Need to integrate from the toa to surface 
    !    so flip order of layers (ajt, 05/21/13)
    p(L0:L1) = PRES(L0:L1)
    IF (L1 .GT. 1) THEN
       IF(PRES(L0+1) .LT. PRES(L0)) THEN
          p(1:L1) = p(L1:1:-1)
       ENDIF
    ENDIF

    ! assign weight to TOA
    L = 1
    h(L) = 1./PRES(L0) * &
           ABS(( -1e0*p(L) + ( (p(L+1)-p(L))/(LOG(p(L+1)/p(L))) ) ) )
      
    ! assign weight to surface
    L = L1 - L0 + 1
    h(L) = 1./PRES(L0) * &
           ABS((  p(L) - ( (p(L)-p(L-1))/(LOG(p(L)/p(L-1))) ) ) )
      
    ! assign weights to middle layers
    DO L=2,L1-L0
       h(L) = 1./PRES(L0) * ABS( &
              ( -1e0*p(L) + ( (p(L+1)-p(L))/(LOG(p(L+1)/p(L))) ) ) + &
              (      p(L) - ( (p(L)-p(L-1))/(LOG(p(L)/p(L-1))) ) )   )
    ENDDO

    ! Now return to the orientation of the other variables
    IF (L1 .GT. 1) THEN
       ! erp includes bugfix 09/23/2019
       IF(PRES(L0+1) .LT. PRES(L0)) THEN
          h(L0:L1) = h(L1-L0+1:1:-1)
          p(L0:L1) = p(L1-L0+1:1:-1)
       ENDIF
    ENDIF

    PRES_WT(L0:L1) = h(L0:L1)

  ENDSUBROUTINE CALC_PRES_WEIGHT_AJT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: calc_col_ak
!
! !DESCRIPTION: Calcuate the column averaging kernel
!  Aj = (h'A)j / hj
!  See Worden et al, 2015 AMT, Quantifying lower tropospheric 
!  methane concentrations using GOSAT near-IR and AIRS thermal
!  TIR measurements
!\\
!\\
! !INTERFACE: 
!
  SUBROUTINE CALC_COL_AK (AK1, PRES_WT, L0, L1, AK2)
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  :: AK1(MAXLEV, MAXLEV)
    REAL(fp), INTENT(IN)  :: PRES_WT(MAXLEV)
    INTEGER,  INTENT(IN)  :: L0, L1
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: AK2(MAXLEV)
!     
! !REVISION HISTORY:
!  21 Sept 2018 - Yuzhong Zhang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER   :: I, J

    AK2(:) = 0.0_fp
    DO J = L0, L1
       AK2(J) = SUM(AK1(L0:L1,J) * PRES_WT(L0:L1))
    ENDDO
    AK2(L0:L1) = AK2(L0:L1) / PRES_WT(L0:L1)

  ENDSUBROUTINE CALC_COL_AK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: convert_col_ak
!
! !DESCRIPTION: Convert column averaging kernel - erp, July 15, 2020
!  The GOSAT column averaging kenel is defined as:
!  Aj = (h'A)j / hj
!  But AIRS column averaging kernel is defined as:
!  Aj = (h'A)j
!  So we need to convert to the GOSAT format to use the default
!  functions. 
!  See Worden et al, 2015 AMT, Quantifying lower tropospheric 
!  methane concentrations using GOSAT near-IR and AIRS thermal
!  TIR measurements
!\\
!\\
! !INTERFACE: 
!
  SUBROUTINE CONVERT_COL_AK (AK1, PRES_WT, L0, L1, AK2)
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  :: AK1(MAXLEV)
    REAL(fp), INTENT(IN)  :: PRES_WT(MAXLEV)
    INTEGER,  INTENT(IN)  :: L0, L1
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: AK2(MAXLEV)
!     
! !REVISION HISTORY:
!  13 Sept 2019 - E. Penn - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: L
      
    AK2(:) = 0.0_fp
    DO L = L0, L1
       ! prevent dividing by zero 
       IF (PRES_WT(L) .LE. 0d0) THEN 
          AK2(L) = 0.0_fp
       ELSE
          AK2(L) = AK1(L) / PRES_WT(L)
       ENDIF
    ENDDO

  ENDSUBROUTINE CONVERT_COL_AK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: calc_airs_bias_correction
!
! !DESCRIPTION: Correct bias in the AIRS ORIGINAL_SPECIES variable. 
!  From email from J Worden, September 2018. 
!\\
!\\
! !INTERFACE: 
!
  SUBROUTINE CALC_AIRS_BIAS_CORRECTION(CH4, AK, PRES, L0, L1, CH4_CORRECTED)
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  :: CH4(MAXLEV)
    REAL(fp), INTENT(IN)  :: AK(MAXLEV, MAXLEV)
    REAL(fp), INTENT(IN)  :: PRES(MAXLEV)
    INTEGER,  INTENT(IN)  :: L0, L1
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: CH4_CORRECTED(MAXLEV)
!     
! !REVISION HISTORY:
!  13 Sept 2019 - E. Penn - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I, NTROPINDS
    INTEGER  :: TROPINDS(MAXLEV)
    REAL(fp) :: BIAS(MAXLEV), DELTA(MAXLEV)
      
    ! Find indicies of pres > 150 (~trop)
    TROPINDS(:) = 0
    NTROPINDS = 1
    DO I = L0,L1
       IF ( PRES(I) > 150.0_fp ) THEN
          TROPINDS(NTROPINDS) = I
          NTROPINDS = NTROPINDS+1
       ENDIF
    ENDDO
     
    BIAS(:) = 0.0_fp
    BIAS(TROPINDS(1:NTROPINDS)) = -0.038_fp - 0.006_fp
    ! perform matrix multiplication AK * BIAS
    DO I = L0,L1
       DELTA(I) = SUM(AK(I,L0:L1) * BIAS(L0:L1))
    ENDDO
    CH4_CORRECTED(L0:L1) = CH4(L0:L1) + DELTA(L0:L1) * CH4(L0:L1)
      
  ENDSUBROUTINE CALC_AIRS_BIAS_CORRECTION
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: calc_airs_xch4
!
! !DESCRIPTION: Calcuate the XCH4 as observed by AIRS
!  Average CH4 mixing ratio weighted by pressure layer thickness
!\\
!\\
! !INTERFACE: 
!
  SUBROUTINE CALC_AIRS_XCH4 (CH4, PRES_WT, L0, L1, XCH4)
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN) :: CH4(MAXLEV)
    REAL(fp), INTENT(IN) :: PRES_WT(MAXLEV)
    INTEGER,  INTENT(IN) :: L0, L1
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: XCH4
!     
! !REVISION HISTORY:
!  21 Sept 2018 - Yuzhong Zhang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    XCH4 = SUM( CH4(L0:L1) * PRES_WT(L0:L1))
      
  ENDSUBROUTINE CALC_AIRS_XCH4
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: calc_gc_xch4
!
! !DESCRIPTION: Add computation to record separate tropospheric contribution
!  to XCH4, zyz, Sept 19, 2018
!  We can record XCH4_prior, XCH4_trop, XCH4_strat, 
!  so startospheric bias correction can be done offline
!  XCH4m = XCH4a + XCH4c 
!        = SUM(l) (wl*pl) + SUM(l) (wl*al*(ml-pl))
!        = SUM(l) ((1-al)*wl*pl) + SUM(l) (al*wl*ml)
!        = SUM(l) ((1-al)*wl*pl) +       ===> XCH4m_prior
!          SUM(l<=LTROP) (al*wl*ml) +    ===> XCH4m_trop
!          SUM(l>LTROP) (al*wl*ml)       ===> XCH4m_strat
!  .
!  wl: weight for layer l
!  pl: prior mixing ratio at layer l
!  ml: model mixing ratio at layer l
!  al: column averaging kernel at layer l
!  LTROP: layer of tropopause
!\\
!\\
! !INTERFACE: 
!
  SUBROUTINE CALC_GC_XCH4( GC_CH4, AK2, PRES_WT, PRIOR, L0, L1, &
                           LTROP, WT_LTROP, &
                           XCH4, XCH4_prior, XCH4_trop, XCH4_strat)
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  :: GC_CH4(MAXLEV), AK2(MAXLEV) 
    REAL(fp), INTENT(IN)  :: PRIOR(MAXLEV), PRES_WT(MAXLEV)
    INTEGER,  INTENT(IN)  :: L0, L1, LTROP
    REAL(fp), INTENT(IN)  :: WT_LTROP
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: XCH4, XCH4_prior, XCH4_trop, XCH4_strat
!     
! !REVISION HISTORY:
!  19 Sept 2018 - Yuzhong Zhang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: L

    ! Calculate XCH4_prior
    XCH4_PRIOR = 0.0_fp
    DO L = L0, L1
       XCH4_PRIOR = XCH4_PRIOR + (1.0_fp - AK2(L)) * PRES_WT(L) * PRIOR(L)
    ENDDO

    ! Calculate XCH4_trop
    XCH4_trop = 0.0_fp
    DO L = L0, LTROP
       XCH4_trop = XCH4_trop + AK2(L) * PRES_WT(L) * GC_CH4(L)
    ENDDO
    IF (LTROP .ge. L0) THEN
       XCH4_trop = XCH4_trop - (1.0_fp - WT_LTROP) * &
                   AK2(LTROP) * PRES_WT(LTROP) * GC_CH4(LTROP)
    ENDIF

    ! Calculate XCH4_strat
    XCH4_strat = 0.0_fp
    DO L = LTROP, L1
       XCH4_strat = XCH4_strat + AK2(L) * PRES_WT(L) * GC_CH4(L)
    ENDDO
    IF (LTROP .le. L1) THEN
       XCH4_strat = XCH4_strat - WT_LTROP * &
                    AK2(LTROP) * PRES_WT(LTROP) * GC_CH4(LTROP)
    ENDIF
 
    ! Calculate total column
    XCH4 = XCH4_PRIOR + XCH4_trop + XCH4_strat

  END SUBROUTINE CALC_GC_XCH4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: h_to_hprime
!
! !DESCRIPTION: Calculate HPRIME pressure edges from equation 11 of Keppens 2019
!\\
!\\
! !INTERFACE: 
!
  SUBROUTINE H_TO_HPRIME( OBS_PEDGE, L0, nlev_obs, HPRIME )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: nlev_obs, L0
    REAL(fp), INTENT(IN)  :: OBS_PEDGE(nlev_obs)
!         
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: HPRIME(nlev_obs+1)
      !     
! !REVISION HISTORY:
!  06 Oct 2020 - E. Penn - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: L

    ! Initialize
    HPRIME(:) = -9999.0_fp  ! this is your missing value
         
    HPRIME(L0) = OBS_PEDGE(L0)
    HPRIME(nlev_obs+1) = OBS_PEDGE(nlev_obs)
    ! Loop over each pressure level of observation grid
    DO L = L0+1, nlev_obs
       HPRIME(L) = 0.5_fp*OBS_PEDGE(L) + 0.5_fp*OBS_PEDGE(L-1)
    ENDDO
       
  END SUBROUTINE H_TO_HPRIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: extend_gc
!
! !DESCRIPTION: Extend GEOS-Chem pressure leveys so they cover the full
! vertical range of the observations
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EXTEND_GC( GC_PEDGE, OBS_PEDGE, L0, &
                        nlev_gc, nlev_obs, GC_PEDGE_EXT )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: nlev_obs, nlev_gc, L0
    REAL(fp), INTENT(IN) :: GC_PEDGE(nlev_gc+1)
    REAL(fp), INTENT(IN) :: OBS_PEDGE(nlev_obs)
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: GC_PEDGE_EXT(nlev_gc+1)
!     
! !REVISION HISTORY:
!  13 Sept 2019 - E. Penn - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
    ! Initialize
    GC_PEDGE_EXT(:) = 0.0_fp
    GC_PEDGE_EXT(:) = GC_PEDGE(:)
         
    ! if observation surf pres. is higher than model, then
    ! extend the surface model layer down to the obs. surface
    IF ( OBS_PEDGE(L0) > GC_PEDGE(1) ) THEN
       GC_PEDGE_EXT(1) = OBS_PEDGE(L0)
    ENDIF
         
    ! if observation TOA pres. is lower than model, then 
    ! extend the top model layer up to the obs. TOA
    IF ( OBS_PEDGE(nlev_obs) < GC_PEDGE(nlev_gc+1) ) THEN
       GC_PEDGE_EXT(nlev_gc+1) = OBS_PEDGE(nlev_obs)
    ENDIF
       
  END SUBROUTINE EXTEND_GC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: get_overlap_map
!
! !DESCRIPTION: OVERLAP_MAP is W in eq 13 of Keppens 2019
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_OVERLAP_MAP( GC_PEDGE, OBS_PEDGE, &
                              L0, nlev_gc, nlev_obs, OVERLAP_MAP)
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: nlev_gc
    INTEGER,  INTENT(IN)  :: nlev_obs
    REAL(fp), INTENT(IN)  :: GC_PEDGE(nlev_gc)
    REAL(fp), INTENT(IN)  :: OBS_PEDGE(nlev_obs)
    INTEGER,  INTENT(IN)  :: L0 ! lowest valid observation level 
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: OVERLAP_MAP(nlev_obs-1,nlev_gc-1)
!
! !REVISION HISTORY:
!  23 Sep 2020 - Elise Penn - get map of layer overlaps based on interpolation
!                             in Langerock et a. 2015 and used in equation 13 of
!                             Keppens et al. 2019. 
!                             See description of "mass-conserved regridding" in
!                             Keppens et al. 2019: 
!                             https://doi.org/10.5194/amt-12-4379-2019
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: LGC, LTM
    REAL(fp) :: DIFF
    REAL(fp) :: LOW_GC, HI_GC
    REAL(fp) :: LOW_OBS, HI_OBS

    !=================================================================
    ! GET_OVERLAP_MAP begins here!
    !=================================================================

    ! Initialize
    OVERLAP_MAP(:,:) = 0e+0_fp
         
    ! Loop over each pressure level of observation retrieval grid
    DO LTM = L0, nlev_obs-1

       LOW_OBS = OBS_PEDGE(LTM+1)
       HI_OBS  = OBS_PEDGE(LTM)
                 
       ! Find the levels from GC that bracket level LTM
       DO LGC = 1, nlev_gc-1

          LOW_GC = GC_PEDGE(LGC+1)
          HI_GC  = GC_PEDGE(LGC)

          ! Match GEOS-Chem level to observation level
          IF ( ( HI_OBS  <= HI_GC .and.    &
                 HI_OBS  >  LOW_GC ) .or.  &
               ( LOW_OBS <= HI_GC .and.    &
                 LOW_OBS >  LOW_GC ) .or.  &
               ( HI_GC   <= HI_OBS .and.   &
                 HI_GC   >  LOW_OBS ) .or. &
               ( LOW_GC  <= HI_OBS .and.   &
                 LOW_GC  >  LOW_OBS ) ) THEN

             DIFF             = HI_GC- LOW_GC
             OVERLAP_MAP(LTM,LGC) = ( MIN(HI_OBS,HI_GC) - &
                                      MAX(LOW_OBS,LOW_GC) ) / DIFF
             
          ENDIF

       ENDDO

    ENDDO
         
  END SUBROUTINE GET_OVERLAP_MAP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: mass_interp
!
! !DESCRIPTION: Perform interpolation from model levels to the levels of your
!  observation, erp, Oct 6, 2020
!  The interpolation redistributes mass between model and 
!  observation layers, then returns it to the edges of 
!  the layers for application of the AK and pressure weights.   
!  Based on equation 13 of Keppens 2019
!\\
!\\
! !INTERFACE: 
!
  SUBROUTINE MASS_INTERP( GC_PEDGE, OBS_PEDGE, GC_CH4_NATIVE, L0, &
                          nlev_gc, nlev_obs, CH4_INTERP_EDGES )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN) :: nlev_obs, nlev_gc, L0
    REAL(fp), INTENT(IN) :: GC_PEDGE(nlev_gc+1)
    REAL(fp), INTENT(IN) :: GC_CH4_NATIVE(nlev_gc)
    REAL(fp), INTENT(IN) :: OBS_PEDGE(nlev_obs)
!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: CH4_INTERP_EDGES(nlev_obs)
!     
! !REVISION HISTORY:
!  06 Oct 2020 - E. Penn - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: GC_PEDGE_EXT(nlev_gc+1), OBS_HPRIME(nlev_obs+1)
    REAL(fp) :: OVERLAP_MAP(nlev_obs,nlev_gc)
    REAL(fp) :: M_obs(nlev_obs), M_gc(nlev_gc)
    REAL(fp) :: CH4_INTERP_INTEGRATED(nlev_obs)
    REAL(fp) :: CH4_MODEL_INTEGRATED(nlev_gc)

    INTEGER  :: LGC, LOBS
         
    ! Initialize
    CH4_INTERP_EDGES(:) = 0.0_fp ! note your missing value is 0
         
    ! extend GEOS-Chem so it covers the full vertical range 
    !   of the observations
    CALL EXTEND_GC( GC_PEDGE, OBS_PEDGE, L0, nlev_gc+1, nlev_obs, &
                    GC_PEDGE_EXT )

    ! calculate HPRIME pressure edges from equation 11
    CALL H_TO_HPRIME( OBS_PEDGE, L0, nlev_obs, OBS_HPRIME )

    ! OVERLAP_MAP is W in eq 13
    CALL GET_OVERLAP_MAP( GC_PEDGE_EXT, OBS_HPRIME, L0, nlev_gc+1, &
                          nlev_obs+1, OVERLAP_MAP )
         
    ! M_gc and M_obs are M_in and M_out from eq 14
    ! They are diagonal matrices, so we can use a vector
    M_gc(:) = 0.0_fp
    DO LGC = 1, nlev_gc
       M_gc(LGC) = GC_PEDGE_EXT(LGC) - GC_PEDGE_EXT(LGC+1)
    ENDDO
    M_obs(:) = 0.0_fp
    DO LOBS = L0, nlev_obs
       M_obs(LOBS) = OBS_HPRIME(LOBS) - OBS_HPRIME(LOBS+1)
    ENDDO
         
    ! Intermediate steps for eq 14:
    ! 1) M_in * x
    CH4_MODEL_INTEGRATED(:) = 0.0_fp
    DO LGC = 1, nlev_gc
       CH4_MODEL_INTEGRATED(LGC) = M_gc(LGC) * GC_CH4_NATIVE(LGC)
    ENDDO

    ! 2) W * M_in * x (matrix multiply W and M_in*x)
    CH4_INTERP_INTEGRATED(:) = 0.0_fp
    DO LGC = 1, nlev_gc
       DO LOBS = L0, nlev_obs
          CH4_INTERP_INTEGRATED(LOBS) = CH4_INTERP_INTEGRATED(LOBS) + &
               OVERLAP_MAP(LOBS,LGC) * CH4_MODEL_INTEGRATED(LGC)
       ENDDO
    ENDDO

    ! 3) inv(M_out) * W * M_in * x
    DO LOBS = L0, nlev_obs
       ! inv(M_out) = 1/M_out because it is diagonal
       CH4_INTERP_EDGES(LOBS) = 1.0_fp/M_obs(LOBS) * CH4_INTERP_INTEGRATED(LOBS)
    ENDDO
         
  END SUBROUTINE MASS_INTERP
!EOC
END MODULE AIRS_CH4_MOD
