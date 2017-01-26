!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: merra2_read_mod.F90
!
! !DESCRIPTION: Module MERRA2\_READ\_MOD contains subroutines for reading the 
!  MERRA2 data from disk (in netCDF format).
!\\
!\\
! !INTERFACE: 
!
MODULE Merra2_Read_Mod
!
! !USES:
!
  ! NcdfUtil modules for netCDF I/O
  USE m_netcdf_io_open                    ! netCDF open
  USE m_netcdf_io_get_dimlen              ! netCDF dimension queries
  USE m_netcdf_io_read                    ! netCDF data reads
  USE m_netcdf_io_close                   ! netCDF close

  ! GEOS-Chem modules
  USE CMN_SIZE_MOD                        ! Size parameters
#if defined( BPCH_DIAG )
  USE CMN_DIAG_MOD                        ! Diagnostic arrays & counters
  USE DIAG_MOD,      ONLY : AD21          ! Array for ND21 diagnostic  
  USE DIAG_MOD,      ONLY : AD66          ! Array for ND66 diagnostic  
  USE DIAG_MOD,      ONLY : AD67          ! Array for ND67 diagnostic
#endif
  USE ERROR_MOD,     ONLY : ERROR_STOP    ! Stop w/ error message
  USE PhysConstants                       ! Physical constants
  USE Precision_Mod                       ! Flexible precision definitions
  USE TIME_MOD                            ! Date & time routines
  USE TRANSFER_MOD                        ! Routines for casting 

  IMPLICIT NONE
  PRIVATE

# include "netcdf.inc"                    ! Include file for netCDF library
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Check_Dimensions
  PRIVATE :: Merra2_Read_A3cld
  PRIVATE :: Merra2_Read_A3dyn
  PRIVATE :: Merra2_Read_A3mstC
  PRIVATE :: Merra2_Read_A3mstE
  PRIVATE :: Get_Resolution_String
!
! !PUBLIC MEMBER FUNCTIONS:
! 
  PUBLIC  :: Merra2_Read_CN
  PUBLIC  :: Merra2_Read_A1
  PUBLIC  :: Merra2_Read_A3
  PUBLIC  :: Merra2_Read_I3_1
  PUBLIC  :: Merra2_Read_I3_2
  PUBLIC  :: Cleanup_Merra2_Read
!
! !REMARKS:
!  Assumes that you have a netCDF library (either v3 or v4) installed on 
!  your system. 
!
! !REVISION HISTORY:
!  12 Aug 2015 - R. Yantosca - Initial version, based on geosfp_read_mod.F90
!  03 Dec 2015 - R. Yantosca - Add file ID's as module variables
!  03 Dec 2015 - R. Yantosca - Add CLEANUP_MERRA2_READ to close any open 
!                              netCDF files left at the end of a simulation
!  02 Feb 2016 - E. Lundgren - Block of diagnostics with if defined BPCH
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  ! netCDF file ID's
  INTEGER :: fA1     = -1
  INTEGER :: fA3cld  = -1
  INTEGER :: fA3dyn  = -1
  INTEGER :: fA3mstC = -1
  INTEGER :: fA3mstE = -1
  INTEGER :: fI3_1   = -1
  INTEGER :: fI3_2   = -1

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Resolution_String
!
! !DESCRIPTION: Function Get\_Resolution\_String returns the proper filename 
!  extension for the GEOS-Chem horizontal grid resolution.  This is used to
!  construct the various file names.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_Resolution_String() RESULT( resString )
!
! !RETURN VALUE:
!
    CHARACTER(LEN=255) :: resString
! 
! !REVISION HISTORY:
!  12 Aug 2015 - R. Yantosca - Initial version, based on geosfp_read_mod.F90
!  13 Aug 2015 - R. Yantosca - MERRA2 data is now storead as netCDF4 (.nc4)
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( GRID4x5 )
    resString = '4x5.nc4'
     
#elif defined( GRID2x25 ) 
    resString = '2x25.nc4'

#elif defined( GRID1x125 )
    resString = '1x125.nc4'

#elif defined( GRID1x1 ) 
    resString = '1x1.nc4'

#elif defined( GRID05x0666 )
    resString = '05x0666.nc4'

#elif defined( GRID05x0625 ) && defined( NESTED_AS )
    resString = '05x0625.AS.nc4'

#elif defined( GRID05x0625 ) && defined( NESTED_EU )
    resString = '05x0625.EU.nc4'

#elif defined( GRID05x0625 ) && defined( NESTED_NA )
    resString = '05x0625.NA.nc4'

#elif defined( GRID05x0625 )
    resString = '05x0625.nc4'

#endif

  END FUNCTION Get_Resolution_String
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_Dimensions
!
! !DESCRIPTION: Subroutine CHECK\_DIMENSIONS checks to see if dimensions read 
!  from the netCDF file match the defined GEOS-Chem dimensions.  If not, then 
!  it will stop the GEOS-Chem simulation with an error message.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Check_Dimensions( lon, lat, lev, time, time_expected, caller )
!
! !INPUT PARAMETERS:
!
    INTEGER,          OPTIONAL, INTENT(IN)  :: lon            ! Lon dimension
    INTEGER,          OPTIONAL, INTENT(IN)  :: lat            ! Lat dimension
    INTEGER,          OPTIONAL, INTENT(IN)  :: lev            ! Alt dimension
    INTEGER,          OPTIONAL, INTENT(IN)  :: time           ! Time dimension
    INTEGER,          OPTIONAL, INTENT(IN)  :: time_expected  ! Expected # of 
                                                              !  time slots
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN)  :: caller         ! Name of caller
                                                              !  routine
! 
! !REMARKS:
!  Call this routine with keyword arguments, e.g
!     CALL CHECK_DIMENSION( lon=X,  lat=Y,           lev=Z,         &
!                           time=T, time_expected=8, caller=caller )
!
! !REVISION HISTORY:
!  02 Feb 2012 - R. Yantosca - Initial version
!  03 Feb 2012 - R. Yantosca - Now pass the caller routine name as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Error message string
    CHARACTER(LEN=255) :: errMsg                  

    ! Error check longitude dimension 
    IF ( PRESENT( lon ) ) THEN
       IF ( lon /= IIPAR ) THEN
          WRITE( errMsg, 100 ) lon, IIPAR
 100      FORMAT( 'Longitude dimension (', i5, &
                  ' ) does not match IIPAR ( ', i5, ')!' )
          CALL ERROR_STOP( errMsg, caller )
       ENDIF
    ENDIF


    ! Error check longitude dimension 
    IF ( PRESENT( lat ) ) THEN
       IF ( lat /= JJPAR ) THEN
          WRITE( errMsg, 110 ) lat, JJPAR
 110      FORMAT( 'Latitude dimension (', i5, &
                  ' ) does not match JJPAR ( ', i5, ')!' )
          CALL ERROR_STOP( errMsg, caller )
       ENDIF
    ENDIF


    ! Error check longitude dimension 
    IF ( PRESENT( lev ) ) THEN
       IF ( lev /= LGLOB ) THEN
          WRITE( errMsg, 120 ) lev, LGLOB
 120      FORMAT( 'Levels dimension (', i5, &
                  ' ) does not match LGLOB ( ', i5, ')!' )
          CALL ERROR_STOP( errMsg, caller )
       ENDIF
    ENDIF

    ! Error check time dimension 
    IF ( PRESENT( time ) .and. PRESENT( time_expected ) ) THEN
       IF ( time /= time_expected ) THEN
          WRITE( errMsg, 130 ) time, time_expected
 130      FORMAT( 'Time dimension (', i5, &
                  ' ) does not match expected # of times ( ', i5, ')!' )
          CALL ERROR_STOP( errMsg, caller )
       ENDIF
    ENDIF

  END SUBROUTINE Check_Dimensions
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Merra2_Read_cn
!
! !DESCRIPTION: Routine to read variables and attributes from a MERRA2
!  met fields file containing constant (CN) data.  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Merra2_Read_CN( Input_Opt, State_Met )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !REMARKS:
!  Even though the netCDF file is self-describing, the MERRA2 data, 
!  dimensions, and units are pre-specified according to the GMAO MERRA2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  12 Aug 2015 - R. Yantosca - Initial version, based on geosfp_read_mod.F90
!  13 Aug 2015 - R. Yantosca - Bug fix: change CN date to 2015/01/01
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fCN                ! netCDF file ID
    INTEGER            :: X, Y, T            ! netCDF file dimensions
    CHARACTER(LEN=16)  :: stamp              ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file            ! netCDF file name
    CHARACTER(LEN=255) :: v_name             ! netCDF variable name
    CHARACTER(LEN=255) :: dir                ! Data directory path
    CHARACTER(LEN=255) :: errMsg             ! Error message
    CHARACTER(LEN=255) :: caller             ! Name of this routine
                                
    ! Arrays                                 
    INTEGER            :: st3d(3), ct3d(3)   ! Start + count, for 3D arrays 
    REAL*4             :: Q(IIPAR,JJPAR)     ! Temporary data arrray

    !======================================================================
    ! Open the netCDF file
    !======================================================================

    ! Name of this routine (for error printout)
    caller  = "MERRA2_READ_CN (merra2_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( Input_Opt%MERRA2_DIR )
    CALL Expand_Date( dir, 20150101, 000000 )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String() 
    nc_file = 'MERRA2.YYYYMMDD.CN.' // TRIM( nc_file ) 
    CALL Expand_Date( nc_file, 20150101, 000000 )

    ! Construct complete file path
    nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
    ! Open netCDF file
    CALL NcOp_Rd( fCN, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fCN, 'lon',   X )
    CALL NcGet_DimLen( fCN, 'lat',   Y )
    CALL NcGet_DimLen( fCN, 'time',  T )

    ! Make sure the dimensions of the file are valid
    CALL Check_Dimensions( lon=X,           lat=Y,         time=T,  &
                           time_expected=1, caller=caller          )

    !======================================================================
    ! Read data from netCDF file
    !======================================================================

    ! Start and count indices
    st3d   = (/ 1,     1,     1 /)
    ct3d   = (/ IIPAR, JJPAR, 1 /)

    ! Read FRLAKE
    v_name = "FRLAKE"
    CALL NcRd( Q, fCN, TRIM(v_name), st3d, ct3d )
    State_Met%FRLAKE = Q

    ! Read FRLAND
    v_name = "FRLAND"
    CALL NcRd( Q, fCN, TRIM(v_name), st3d, ct3d )
    State_Met%FRLAND = Q

    ! Read FRLANDIC
    v_name = "FRLANDIC"
    CALL NcRd( Q, fCN, TRIM(v_name), st3d, ct3d )
    State_Met%FRLANDIC = Q
    
    ! Read FROCEAN
    v_name = "FROCEAN"
    CALL NcRd( Q, fCN, TRIM(v_name), st3d, ct3d )
    State_Met%FROCEAN = Q
    
    ! Read PHIS
    v_name = "PHIS"
    CALL NcRd( Q, fCN, TRIM(v_name), st3d, ct3d )
    State_Met%PHIS = Q

    ! Echo info
    stamp = TimeStamp_String( 20110101, 000000 )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 5  MERRA2 CN     met fields for ', a )

    !======================================================================
    ! Cleanup and quit
    !======================================================================

    ! Convert PHIS from [m2/s2] to [m]
    State_Met%PHIS = State_Met%PHIS / g0

#if defined( BPCH_DIAG )
    ! ND67 diagnostic 
    IF ( ND67 > 0 ) THEN
       AD67(:,:,15) = AD67(:,:,15) + State_Met%PHIS  ! Sfc geopotential [m]
    ENDIF
#endif

    ! Close netCDF file
    CALL NcCl( fCN )

  END SUBROUTINE Merra2_Read_CN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Merra2_Read_a1
!
! !DESCRIPTION: Routine to read variables and attributes from a MERRA2
!  met fields file containing 1-hr time-averaged (A1) data.  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Merra2_Read_A1( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
! 
    INTEGER,        INTENT(IN)    :: YYYYMMDD   ! GMT date in YYYY/MM/DD format
    INTEGER,        INTENT(IN)    :: HHMMSS     ! GMT time in hh:mm:ss   format
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
!
! !REMARKS:
!  Even though the netCDF file is self-describing, the MERRA2 data, 
!  dimensions, and units are pre-specified according to the GMAO MERRA2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!                                                                             .
!  Special handling for surface precipitation fields:
!  ---------------------------------------------------------------------------
!  In MERRA2 (and in MERRA), the PRECTOT etc. surface precipitation
!  met fields fields have units of [kg/m2/s].  In all other GEOS 
!  versions, PREACC and PRECON have units of [mm/day].  
!                                                                             .
!  Therefore, for backwards compatibility with existing code, apply 
!  the following unit conversion to the GEOS-5 PRECTOT and PRECCON 
!  fields:
!                                                                             .
!      kg  |    m3    | 86400 s | 1000 mm
!    ------+----------+---------+--------- = 86400 
!     m2 s |  1000 kg |  day    |   m
!               ^
!               |
!        1 / density of water 
!
! !REMARKS:
!  MERRA2 pressure quantities are stored on disk with units of Pa.  For now
!  we will convert to hPa for compatibility with GEOS-Chem.  But in the near
!  future we will probably recode GEOS-Chem to use the native units, which
!  facilitate GEOS-Chem HP development.
!
! !REVISION HISTORY:
!  12 Aug 2015 - R. Yantosca - Initial version, based on geosfp_read_mod.F90
!  23 Sep 2015 - E. Lundgren - Now assign SWGDN to State_Met SWGDN not RADSWG
!  03 Dec 2015 - R. Yantosca - Now open file only once per day
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: X, Y, T            ! netCDF file dimensions
    INTEGER            :: I, J               ! Do-loop indices
    INTEGER            :: time_index         ! Read this time slice of data
    CHARACTER(LEN=16)  :: stamp              ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file            ! netCDF file name
    CHARACTER(LEN=255) :: v_name             ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                ! Data directory path
    CHARACTER(LEN=255) :: errMsg             ! Error message
    CHARACTER(LEN=255) :: caller             ! Name of this routine

    ! Saved scalars
    INTEGER, SAVE      :: lastDate = -1      ! Stores last YYYYMMDD value
    INTEGER, SAVE      :: lastTime = -1      ! Stores last hhmmss value
    LOGICAL, SAVE      :: first    = .TRUE.  ! First time reading data?

    ! Arrays                                 
    INTEGER            :: st3d(3), ct3d(3)   ! Start + count, for 3D arrays 
    REAL*4             :: Q(IIPAR,JJPAR)     ! Temporary data arrray

    !======================================================================
    ! Skip if we have already read data for this date & time
    !======================================================================
    IF ( YYYYMMDD == lastDate .and. HHMMSS == lastTime ) THEN
       stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
       WRITE( 6, 20 ) stamp
 20    FORMAT( '     - MERRA2 A1 met fields for ', a,  &
               ' have been read already'                  ) 
       RETURN
    ENDIF

    !======================================================================
    ! Select the proper time slice
    !======================================================================
    
    ! Name of this routine (for error printout)
    caller  = "MERRA2_READ_A1 (merra2_read_mod.F90)"

    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 10000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > 24 ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 24!' )
       CALL Error_Stop( errMsg, caller )
    ENDIF
    
    !======================================================================
    ! Open the netCDF file only when necessary
    !======================================================================
    IF ( time_index == 1 .or. first ) THEN     
    
       ! Replace time & date tokens in the file name
       dir     = TRIM( Input_Opt%MERRA2_DIR )
       CALL Expand_Date( dir, YYYYMMDD, HHMMSS )

       ! Replace time & date tokens in the file name
       nc_file = Get_Resolution_String()
       nc_file = 'MERRA2.YYYYMMDD.A1.' // TRIM( nc_file )
       CALL Expand_Date( nc_file, YYYYMMDD, HHMMSS )

       ! Construct complete file path
       nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )

       ! Open netCDF file
       CALL NcOp_Rd( fA1, TRIM( nc_file ) )

       ! Read the dimensions from the netCDF file
       CALL NcGet_DimLen( fA1, 'lon',   X )
       CALL NcGet_DimLen( fA1, 'lat',   Y )
       CALL NcGet_DimLen( fA1, 'time',  T )

       ! Make sure the dimensions of the file are valid
       CALL Check_Dimensions( lon=X,            lat=Y,        time=T,  &
                              time_expected=24, caller=caller         )
    
       ! Reset first-time flag
       first = .FALSE.
    ENDIF

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================

    ! netCDF start & count indices
    st3d      = (/ 1,     1,     time_index /)      
    ct3d      = (/ IIPAR, JJPAR, 1          /)

    ! Read ALBEDO
    v_name = "ALBEDO"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%ALBD = Q

    ! Read CLDTOT
    v_name = "CLDTOT"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%CLDFRC = Q

    ! Read EFLUX
    v_name = "EFLUX"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%EFLUX = Q

    ! Read EVAP
    v_name = "EVAP"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%EVAP = Q

    ! Read FRSEAICE
    v_name = "FRSEAICE"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%FRSEAICE = Q

    ! Read FRSNO
    v_name = "FRSNO"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%FRSNO = Q

    ! Read GRN
    v_name = "GRN"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%GRN = Q

    ! Read GWETROOT
    v_name = "GWETROOT"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%GWETROOT = Q

    ! Read GWETTOP
    v_name = "GWETTOP"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%GWETTOP = Q

    ! Read HFLUX from file
    v_name = "HFLUX"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%HFLUX = Q

    ! Read LAI
    v_name = "LAI"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%LAI = Q

    ! Read LWI
    v_name = "LWI"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%LWI = Q

    ! Read LWGNT 
    v_name = "LWGNT"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%RADLWG = Q

    !-----------------------------------------------------------------------
    ! Comment this out for now, this field isn't needed (bmy, 2/2/12)
    !! Read LWTUP
    !v_name = "LWTUP"
    !CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    !State_Met%LWTUP = Q)
    !-----------------------------------------------------------------------

    ! Read PARDF
    v_name = "PARDF"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%PARDF = Q

    ! Read PARDR
    v_name = "PARDR"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%PARDR = Q

    ! Read PBLH
    v_name = "PBLH"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%PBLH = Q

    ! Read PRECANV
    v_name = "PRECANV"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%PRECANV = Q

    ! Read PRECCON
    v_name = "PRECCON"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%PRECCON = Q

    ! Read PRECLSC
    v_name = "PRECLSC"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%PRECLSC = Q

    ! Read PRECSNO
    v_name = "PRECSNO"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%PRECSNO = Q

    ! Read PRECTOT
    v_name = "PRECTOT"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%PRECTOT = Q

    !-----------------------------------------------------------------------
    ! Comment this out for now, this field isn't needed (bmy, 2/2/12)
    !! Read QV2M
    !v_name = "QV2M"
    !CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    !State_Met%QV2M = Q
    !-----------------------------------------------------------------------

    ! Read SEAICE00
    v_name = "SEAICE00"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%SEAICE00 = Q

    ! Read SEAICE10
    v_name = "SEAICE10"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%SEAICE10 = Q

    ! Read SEAICE20
    v_name = "SEAICE20"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%SEAICE20 = Q

    ! Read SEAICE30
    v_name = "SEAICE30"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%SEAICE30 = Q

    ! Read SEAICE40
    v_name = "SEAICE40"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%SEAICE40 = Q

    ! Read SEAICE50
    v_name = "SEAICE50"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%SEAICE50 = Q

    ! Read SEAICE60 
    v_name = "SEAICE60"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%SEAICE60 = Q

    ! Read SEAICE70
    v_name = "SEAICE70"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%SEAICE70 = Q

    ! Read SEAICE80
    v_name = "SEAICE80"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%SEAICE80 = Q

    ! Read SEAICE90
    v_name = "SEAICE90"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%SEAICE90 = Q

    ! Read SLP
    v_name = "SLP"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%SLP = Q

    ! Read SNODP
    v_name = "SNODP"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%SNODP = Q

    ! Read SNOMAS
    v_name = "SNOMAS"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%SNOMAS = Q

    ! Read SWGDN
    v_name = "SWGDN"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%SWGDN  = Q

    ! Read TO3
    v_name = "TO3"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%TO3 = Q

    ! Read TROPPT
    v_name = "TROPPT"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%TROPP = Q

    ! Read TS
    v_name = "TS"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%TSKIN = Q

    ! Read T2M
    v_name = "T2M"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%TS = Q

    ! Read U10M
    v_name = "U10M"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%U10M = Q

    ! Read USTAR
    v_name = "USTAR"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%USTAR = Q

    ! Read V10M
    v_name = "V10M"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met% V10M = Q

    ! Read Z0M
    v_name = "Z0M"
    CALL NcRd( Q, fA1, TRIM(v_name), st3d, ct3d )
    State_Met%Z0 = Q

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp                    
 10 FORMAT( '     - Found all 44 MERRA2 A1     met fields for ', a )

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! If it's the last time slice, then close the netCDF file
    ! and set the file ID to -1 to indicate that it's closed.
    IF ( time_index == 24 ) THEN
       CALL NcCl( fA1 )
       fA1 = -1
    ENDIF

    ! Increment the # of times A1 fields are read from disk
    CALL Set_Ct_A1( INCREMENT=.TRUE. )
      
    ! Convert surface precip fields from [kg/m2/s] --> [mm/day]
    State_Met%PRECANV = State_Met%PRECANV * 86400e0_fp
    State_Met%PRECCON = State_Met%PRECCON * 86400e0_fp
    State_Met%PRECLSC = State_Met%PRECLSC * 86400e0_fp
    State_Met%PRECTOT = State_Met%PRECTOT * 86400e0_fp
    
    ! Convert pressure quantities from [Pa] -> [hPa]
    State_Met%SLP     = State_Met%SLP     * 1e-2_fp
    State_Met%TROPP   = State_Met%TROPP   * 1e-2_fp

#if defined( BPCH_DIAG )
    ! ND67 diagnostic: surface fields
    IF ( ND67 > 0 ) THEN
       AD67(:,:,1 ) = AD67(:,:,1 ) + State_Met%HFLUX    ! Sens heat flux [W/m2]
       AD67(:,:,2 ) = AD67(:,:,2 ) + State_Met%SWGDN    ! SW rad @ sfc [W/m2]
       AD67(:,:,3 ) = AD67(:,:,3 ) + State_Met%PRECTOT  ! Tot prec [kg/m2/s]
       AD67(:,:,4 ) = AD67(:,:,4 ) + State_Met%PRECCON  ! Sfc conv prec[kg/m2/s]
       AD67(:,:,5 ) = AD67(:,:,5 ) + State_Met%TS       ! T @ 2m height [K]
       AD67(:,:,6 ) = AD67(:,:,6 ) + 0e0                !
       AD67(:,:,7 ) = AD67(:,:,7 ) + State_Met%USTAR    ! Friction vel [m/s]
       AD67(:,:,8 ) = AD67(:,:,8 ) + State_Met%Z0       ! Roughness height [m]
       AD67(:,:,9 ) = AD67(:,:,9 ) + State_Met%PBLH     ! PBL height [m]
       AD67(:,:,10) = AD67(:,:,10) + State_Met%CLDFRC   ! Column cld fraction
       AD67(:,:,11) = AD67(:,:,11) + State_Met%U10M     ! U-wind @ 10m [m/s]
       AD67(:,:,12) = AD67(:,:,12) + State_Met%V10M     ! V-wind @ 10m [m/s]
       AD67(:,:,14) = AD67(:,:,14) + State_Met%ALBD     ! Sfc albedo [unitless]
       AD67(:,:,17) = AD67(:,:,17) + State_Met%TROPP    ! T'pause pressure [hPa]
       AD67(:,:,18) = AD67(:,:,18) + State_Met%SLP      ! Sea level prs [hPa]
       AD67(:,:,19) = AD67(:,:,19) + State_Met%TSKIN    ! Sfc skin temp [K]
       AD67(:,:,20) = AD67(:,:,20) + State_Met%PARDF    ! Diffuse PAR [W/m2]
       AD67(:,:,21) = AD67(:,:,21) + State_Met%PARDR    ! Direct PAR [W/m2]
       AD67(:,:,22) = AD67(:,:,22) + State_Met%GWETTOP  ! Topsoil wetness [frac]
       AD67(:,:,23) = AD67(:,:,23) + State_Met%EFLUX    ! Latent heat flux [W/m2]
    ENDIF
#endif

    ! Save date & time for next iteration
    lastDate = YYYYMMDD
    lastTime = HHMMSS

  END SUBROUTINE Merra2_Read_A1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Merra2_Read_a3
!
! !DESCRIPTION: Convenience wrapper for the following routines which read
!  3-hour time averaged data from disk:
! \begin{itemize}
! \item Merra2\_Read\_A3cld
! \item Merra2\_Read\_A3dyn
! \item Merra2\_Read\_A3mstC
! \item Merra2\_Read\_A3mstE
! \end{itemize}
!
! !INTERFACE:
!
  SUBROUTINE Merra2_Read_A3( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
! 
    INTEGER,        INTENT(IN)    :: YYYYMMDD   ! GMT date in YYYY/MM/DD format
    INTEGER,        INTENT(IN)    :: HHMMSS     ! GMT time in hh:mm:ss   format
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
!
! !REVISION HISTORY:
!  12 Aug 2015 - R. Yantosca - Initial version, based on geosfp_read_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=16) :: stamp            ! Time and date stamp

    ! Saved scalars
    INTEGER, SAVE     :: lastDate = -1    ! Stores last YYYYMMDD value
    INTEGER, SAVE     :: lastTime = -1    ! Stores last hhmmss value

    !======================================================================
    ! Call individual routines for reading A3 data
    !======================================================================

    ! Test to see if we have already read this data in
    IF ( YYYYMMDD == lastDate .and. HHMMSS == lastTime ) THEN
       stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
       WRITE( 6, 20 ) stamp
 20    FORMAT( '     - MERRA2 A3 met fields for ', a,  &
               ' have been read already'                  ) 
       RETURN
    ENDIF

    ! Save date & time for next iteration
    lastDate = YYYYMMDD
    lastTime = HHMMSS

    ! Read all the diffeent A3 files
    CALL Merra2_Read_A3cld ( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
    CALL Merra2_Read_A3dyn ( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
    CALL Merra2_Read_A3mstC( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
    CALL Merra2_Read_A3mstE( YYYYMMDD, HHMMSS, Input_Opt, State_Met )

    !======================================================================
    ! Cleanup and quit
    !======================================================================

    ! Increment the # of times that A3 fields have been read
    CALL Set_Ct_A3( INCREMENT=.TRUE. )

    ! Save date & time for next iteration
    lastDate = YYYYMMDD
    lastTime = HHMMSS

  END SUBROUTINE Merra2_Read_A3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Merra2_Read_a3cld
!
! !DESCRIPTION: Routine to read variables and attributes from a MERRA2
!  met fields file containing 3-hr time-averaged (A3) data (cloud fields).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Merra2_Read_A3cld( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
! 
    INTEGER,        INTENT(IN)    :: YYYYMMDD   ! GMT date in YYYY/MM/DD format
    INTEGER,        INTENT(IN)    :: HHMMSS     ! GMT time in hh:mm:ss   format
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
!
! !REMARKS:
!  Even though the netCDF file is self-describing, the MERRA2 data, 
!  dimensions, and units are pre-specified according to the GMAO MERRA2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  12 Aug 2015 - R. Yantosca - Initial version, based on geosfp_read_mod.F90
!  03 Dec 2015 - R. Yantosca - Now open file only once per day
!  17 Mar 2016 - M. Sulprizio- Read optical depth into State_Met%OPTD instead of
!                              State_Met%OPTDEP (obsolete).
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine

    ! SAVEd scalars
    LOGICAL, SAVE      :: first = .TRUE.           ! First time reading data?

    ! Arrays                                 
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Q(IIPAR,JJPAR,LGLOB)     ! Temporary data arrray

    !======================================================================
    ! Select the proper time slice
    !======================================================================
    
    ! Name of this routine (for error printout)
    caller  = "MERRA2_READ_A3cld (Merra2_read_mod.F90)"

    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 030000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > 8 ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL ERROR_STOP( errMsg, caller )
    ENDIF
    
    !======================================================================
    ! Open the netCDF file only when necessary
    !======================================================================
    IF ( time_index == 1 .or. first ) THEN
       
       ! Replace time & date tokens in the file name
       dir     = TRIM( Input_Opt%MERRA2_DIR )
       CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

       ! Replace time & date tokens in the file name
       nc_file = Get_Resolution_String()
       nc_file = 'MERRA2.YYYYMMDD.A3cld.' // TRIM( nc_file )
       CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

       ! Construct complete file path
       nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )

       ! Open netCDF file
       CALL NcOp_Rd( fA3cld, TRIM( nc_file ) )
       
       ! Read the dimensions from the netCDF file
       CALL NcGet_DimLen( fA3cld, 'lon',   X )
       CALL NcGet_DimLen( fA3cld, 'lat',   Y )
       CALL NcGet_DimLen( fA3cld, 'lev',   Z )
       CALL NcGet_DimLen( fA3cld, 'time',  T )

       ! Make sure the dimensions of the file are valid
       CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                              time=T, time_expected=8, caller=caller )


       ! Reset first-time flag
       first = .FALSE.
    ENDIF

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! netCDF start & count indices
    st4d      = (/ 1,     1,     1,     time_index /)      
    ct4d      = (/ IIPAR, JJPAR, LGLOB, 1          /)

    ! Read CLOUD
    v_name = "CLOUD"
    CALL NcRd( Q, fA3cld, TRIM(v_name), st4d, ct4d )
    CALL TRANSFER_3d( Q, State_Met%CLDF )
    
    ! Read OPTDEPTH
    v_name = "OPTDEPTH"
    CALL NcRd( Q, fA3cld, TRIM(v_name), st4d, ct4d )
    CALL TRANSFER_3d( Q, State_Met%OPTD )

    ! Read QI
    v_name = "QI"
    CALL NcRd( Q, fA3cld, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%QI )

    ! Read QL
    v_name = "QL"
    CALL NcRd( Q, fA3cld, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%QL )

    ! Read TAUCLI
    v_name = "TAUCLI"
    CALL NcRd( Q, fA3cld, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%TAUCLI )

    ! Read TAUCLW
    v_name = "TAUCLW"
    CALL NcRd( Q, fA3cld, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%TAUCLW )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 6  MERRA2 A3cld  met fields for ', a )

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

#if defined( BPCH_DIAG )
    ! ND21 diagnostic: OPTD and CLDF
    IF ( ND21 > 0 ) THEN
       AD21(:,:,1:LD21,1) = AD21(:,:,1:LD21,1) + State_Met%OPTD(:,:,1:LD21)
       AD21(:,:,1:LD21,2) = AD21(:,:,1:LD21,2) + State_Met%CLDF(:,:,1:LD21)
    ENDIF
#endif

    ! If it's the last time slice, then close the netCDF file
    ! and set the file ID to -1 to indicate that it's closed.
    IF ( time_index == 8 ) THEN
       CALL NcCl( fA3cld )
       fA3cld = -1
    ENDIF

  END SUBROUTINE Merra2_Read_A3cld
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Merra2_Read_a3dyn
!
! !DESCRIPTION: Routine to read variables and attributes from a MERRA2
!  met fields file containing 3-hr time-averaged (A3) data (dynamics fields).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Merra2_Read_A3dyn( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput  
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
! 
    INTEGER,        INTENT(IN)    :: YYYYMMDD   ! GMT date in YYYY/MM/DD format
    INTEGER,        INTENT(IN)    :: HHMMSS     ! GMT time in hh:mm:ss   format
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
!
! !REMARKS:
!  Even though the netCDF file is self-describing, the MERRA2 data, 
!  dimensions, and units are pre-specified according to the GMAO MERRA2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  12 Aug 2015 - R. Yantosca - Initial version, based on geosfp_read_mod.F90
!  03 Dec 2015 - R. Yantosca - Now open file only once per day
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine

    ! SAVEd scalars
    LOGICAL, SAVE      :: first = .TRUE.           ! First time we read data?

    ! Arrays                                 
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Q (IIPAR,JJPAR,LGLOB  )  ! Temporary data arrray
    REAL*4             :: Qe(IIPAR,JJPAR,LGLOB+1)  ! Temporary data arrray

    !======================================================================
    ! Select the proper time slice
    !======================================================================
    
    ! Name of this routine (for error printout)
    caller  = "MERRA2_READ_A3dyn (merra2_read_mod.F90)"

    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 030000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > 8 ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL ERROR_STOP( errMsg, 'GEOS57_READ_A1 (geos57_read_mod.F90)' )
    ENDIF
    
    !======================================================================
    ! Open the netCDF file only when necessary
    !======================================================================
    IF ( time_index == 1 .or. first ) THEN
       
       ! Replace time & date tokens in the file name
       dir     = TRIM( Input_Opt%MERRA2_DIR )
       CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

       ! Replace time & date tokens in the file name
       nc_file = Get_Resolution_String()
       nc_file = 'MERRA2.YYYYMMDD.A3dyn.' // TRIM( nc_file )
       CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

       ! Construct complete file path
       nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
       ! Open netCDF file
       CALL NcOp_Rd( fA3dyn, TRIM( nc_file ) )

       ! Read the dimensions from the netCDF file
       CALL NcGet_DimLen( fA3dyn, 'lon',   X )
       CALL NcGet_DimLen( fA3dyn, 'lat',   Y )
       CALL NcGet_DimLen( fA3dyn, 'lev',   Z )
       CALL NcGet_DimLen( fA3dyn, 'time',  T )

       ! Make sure the dimensions of the file are valid
       CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                              time=T, time_expected=8, caller=caller )
       ! Reset first-time flag
       first = .FALSE.
    ENDIF

    !======================================================================
    ! Read data from the netCDF file (only when necessary)
    !======================================================================

    ! netCDF start & count indices
    st4d      = (/ 1,     1,     1,     time_index /)      
    ct4d      = (/ IIPAR, JJPAR, LGLOB, 1          /)

    ! Read DTRAIN
    v_name = "DTRAIN"
    CALL NcRd( Q, fA3dyn, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%DTRAIN )

    !! Read OMEGA  from file
    !v_name = "OMEGA"
    !CALL NcRd( Q, fA3dyn, TRIM(v_name), st4d, ct4d )
    !CALL Transfer_3d( Q, State_Met%OMEGA )

    ! Read RH
    v_name = "RH"
    CALL NcRd( Q, fA3dyn, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%RH )

    ! Read U
    v_name = "U"
    CALL NcRd( Q, fA3dyn, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%U )

    ! Read V
    v_name = "V"
    CALL NcRd( Q, fA3dyn, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%V )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 6  MERRA2 A3dyn  met fields for ', a )

    !======================================================================
    ! Unit conversions, diagnostics, cleanup, and quit
    !======================================================================
    
    ! Convert RH from [1] to [%]
    State_Met%RH = State_Met%RH * 100d0

#if defined( BPCH_DIAG )
    ! ND66 diagnostic: U, V, DTRAIN met fields
    IF ( ND66 > 0 ) THEN
       AD66(:,:,1:LD66,1) = AD66(:,:,1:LD66,1) + State_Met%U     (:,:,1:LD66)
       AD66(:,:,1:LD66,2) = AD66(:,:,1:LD66,2) + State_Met%V     (:,:,1:LD66)
       AD66(:,:,1:LD66,6) = AD66(:,:,1:LD66,6) + State_Met%DTRAIN(:,:,1:LD66)
    ENDIF

    ! ND67 diagnostic: CLDTOPS
    IF ( ND67 > 0 ) THEN
       AD67(:,:,16) = AD67(:,:,16) + State_Met%CLDTOPS         ! [levels]
    ENDIF
#endif

    ! If it's the last time slice, then close the netCDF file
    ! and set the file ID to -1 to indicate that it's closed.
    IF ( time_index == 8 ) THEN
       CALL NcCl( fA3dyn )
       fA3dyn = -1
    ENDIF

  END SUBROUTINE Merra2_Read_A3dyn
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Merra2_Read_a3mstc
!
! !DESCRIPTION: Routine to read variables and attributes from a MERRA2
!  met fields file containing 3-hr time-averaged (A3) data (moist fields,
!  saved on level centers).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Merra2_Read_A3mstC( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
! 
    INTEGER,        INTENT(IN)    :: YYYYMMDD   ! GMT date in YYYY/MM/DD format
    INTEGER,        INTENT(IN)    :: HHMMSS     ! GMT time in hh:mm:ss   format
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
!
! !REMARKS:
!  Even though the netCDF file is self-describing, the MERRA2 data, 
!  dimensions, and units are pre-specified according to the GMAO MERRA2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  12 Aug 2015 - R. Yantosca - Initial version, based on geosfp_read_mod.F90
!  03 Dec 2015 - R. Yantosca - Now open file only once per day
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine

    ! SAVEd scalars
    LOGICAL, SAVE      :: first = .TRUE.           ! First time reading data?
                                    
    ! Arrays                                 
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Q (IIPAR,JJPAR,LGLOB)    ! Temporary data arrray

    !======================================================================
    ! Select the proper time slice
    !======================================================================
    
    ! Name of this routine (for error printout)
    caller  = "MERRA2_READ_A3mstC (merra2_read_mod.F90)"

    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 030000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > 8 ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL Error_Stop( errMsg, caller )
    ENDIF

    !======================================================================
    ! Open the netCDF file only when necessary
    !======================================================================
    IF ( time_index == 1 .or. first ) THEN
       
       ! Replace time & date tokens in the file name
       dir     = TRIM( Input_Opt%MERRA2_DIR )
       CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

       ! Replace time & date tokens in the file name
       nc_file = Get_Resolution_String()
       nc_file = 'MERRA2.YYYYMMDD.A3mstC.' // TRIM( nc_file )
       CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

       ! Construct complete file path
       nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )

       ! Open netCDF file
       CALL NcOp_Rd( fA3mstC, TRIM( nc_file ) )

       ! Read the dimensions from the netCDF file
       CALL NcGet_DimLen( fA3mstC, 'lon',   X )
       CALL NcGet_DimLen( fA3mstC, 'lat',   Y )
       CALL NcGet_DimLen( fA3mstC, 'lev',   Z )
       CALL NcGet_DimLen( fA3mstC, 'time',  T )

       ! Make sure the dimensions of the file are valid
       CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                              time=T, time_expected=8, caller=caller )

       ! Reset first-time flag
       first = .FALSE.
    ENDIF

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! netCDF start & count indices
    st4d      = (/ 1,     1,     1,     time_index /)      
    ct4d      = (/ IIPAR, JJPAR, LGLOB, 1          /)

    ! Read DQRCU  from file
    v_name = "DQRCU"
    CALL NcRd( Q, fA3mstC, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%DQRCU )

    ! Read DQRLSAN
    v_name = "DQRLSAN"
    CALL NcRd( Q, fA3mstC, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%DQRLSAN )

    ! Read REEVAPCN
    v_name = "REEVAPCN"
    CALL NcRd( Q, fA3mstC, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%REEVAPCN )

    ! Read  from file
    v_name = "REEVAPLS"
    CALL NcRd( Q, fA3mstC, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%REEVAPLS )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 4  MERRA2 A3mstC met fields for ', a )

    !======================================================================
    ! Cleanup and quit
    !======================================================================

    ! If it's the last time slice, then close the netCDF file
    ! and set the file ID to -1 to indicate that it's closed.
    IF ( time_index == 8 ) THEN
       CALL NcCl( fA3mstC )
       fA3mstC = -1
    ENDIF

  END SUBROUTINE Merra2_read_A3mstC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Merra2_Read_a3mste
!
! !DESCRIPTION: Routine to read variables and attributes from a MERRA2
!  met fields file containing 3-hr time-averaged (A3) data (moist fields,
!  saved on level edges).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Merra2_read_A3mstE( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
! 
    INTEGER,        INTENT(IN)    :: YYYYMMDD   ! GMT date in YYYY/MM/DD format
    INTEGER,        INTENT(IN)    :: HHMMSS     ! GMT time in hh:mm:ss   format
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
!
! !REMARKS:
!  This routine was automatically generated by the Perl script ncCodeRead, 
!  and was subsequently hand-edited for compatibility with GEOS-Chem.
!                                                                             .
!  Even though the netCDF file is self-describing, the MERRA2 data, 
!  dimensions, and units are pre-specified according to the GMAO MERRA2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  12 Aug 2015 - R. Yantosca - Initial version, based on geosfp_read_mod.F90
!  03 Dec 2015 - R. Yantosca - Now open file only once per day
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine

    ! SAVEd scalars
    LOGICAL, SAVE      :: first = .TRUE.           ! First time reading data?
                                             
    ! Arrays                                 
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Qe(IIPAR,JJPAR,LGLOB+1)  ! Temporary data arrray

    !======================================================================
    ! Select the proper time slice
    !======================================================================
    
    ! Name of this routine (for error printout)
    caller  = "MERRA2_READ_A3mstE (merra2_read_mod.F90)"

    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 030000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > 8 ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL ERROR_STOP( errMsg, caller )
    ENDIF

    !======================================================================
    ! Open the netCDF file only when necessary
    !======================================================================
    IF ( time_index == 1 .or. first ) THEN
       
       ! Replace time & date tokens in the file name
       dir     = TRIM( Input_Opt%MERRA2_DIR )
       CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

       ! Replace time & date tokens in the file name
       nc_file = Get_Resolution_String()
       nc_file = 'MERRA2.YYYYMMDD.A3mstE.' // TRIM( nc_file )
       CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

       ! Construct complete file path
       nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
       ! Open netCDF file
       CALL NcOp_Rd( fA3mstE, TRIM( nc_file ) )

       ! Read the dimensions from the netCDF file
       CALL NcGet_DimLen( fA3mstE, 'lon',   X )
       CALL NcGet_DimLen( fA3mstE, 'lat',   Y )
       CALL NcGet_DimLen( fA3mstE, 'lev',   Z )
       CALL NcGet_DimLen( fA3mstE, 'time',  T )

       ! Make sure the dimensions of the file are valid
       CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z-1,       &
                              time=T, time_expected=8, caller=caller )

       ! Reset first-time flag
       first = .FALSE.
    ENDIF

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================

    ! netCDF start & count indices
    st4d      = (/ 1,     1,     1,       time_index /)      
    ct4d      = (/ IIPAR, JJPAR, LGLOB+1, 1          /)

    ! Read CMFMC (only in MERRA2*.nc files)
    v_name = "CMFMC"
    CALL NcRd( Qe, fA3mstE, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe, State_Met%CMFMC )

    ! Read PFICU
    v_name = "PFICU"
    CALL NcRd( Qe, fA3mstE, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe, State_Met%PFICU )

    ! Read PFILSAN
    v_name = "PFILSAN"
    CALL NcRd( Qe, fA3mstE, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe, State_Met%PFILSAN )

    ! Read PFLCU
    v_name = "PFLCU"
    CALL NcRd( Qe, fA3mstE, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe, State_Met%PFLCU )

    ! Read  from file
    v_name = "PFLLSAN"
    CALL NcRd( Qe, fA3mstE, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe, State_Met%PFLLSAN )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 4  MERRA2 A3mstE met fields for ', a )

    !=================================================================
    ! Diagnostics, cleanup and quit
    !=================================================================

    ! CLDTOPS = highest location of CMFMC in the column (I,J)
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       State_Met%CLDTOPS(I,J) = 1
       DO L = LLPAR, 1, -1
          IF ( State_Met%CMFMC(I,J,L) > 0d0 ) THEN
             State_Met%CLDTOPS(I,J) = L + 1
             EXIT
          ENDIF
       ENDDO
    ENDDO
    ENDDO

#if defined( BPCH_DIAG )
    ! ND66 diagnostic: CMFMC met field
    IF ( ND66 > 0 ) THEN
       AD66(:,:,1:LD66,5) = AD66(:,:,1:LD66,5) + State_Met%CMFMC(:,:,1:LD66)
    ENDIF
#endif

    ! If it's the last time slice, then close the netCDF file
    ! and set the file ID to -1 to indicate that it's closed.
    IF ( time_index == 8 ) THEN
       CALL NcCl( fA3mstE )
       fA3mstE = -1
    ENDIF

  END SUBROUTINE Merra2_Read_A3mstE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Merra2_Read_I3_1
!
! !DESCRIPTION: Routine to read variables and attributes from a MERRA2
!  met fields file containing 3-hr instantaneous (I3) data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Merra2_Read_I3_1( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
! 
    INTEGER,        INTENT(IN)    :: YYYYMMDD   ! GMT date in YYYY/MM/DD format
    INTEGER,        INTENT(IN)    :: HHMMSS     ! GMT time in hh:mm:ss   format
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
!
! !REMARKS:
!  Even though the netCDF file is self-describing, the MERRA2 data, 
!  dimensions, and units are pre-specified according to the GMAO MERRA2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
!  MERRA2 pressure quantities are stored on disk with units of Pa.  For now
!  we will convert to hPa for compatibility with GEOS-Chem.  But in the near
!  future we will probably recode GEOS-Chem to use the native units, which
!  facilitate GEOS-Chem HP development.
!
! !REVISION HISTORY:
!  12 Aug 2015 - R. Yantosca - Initial version, based on geosfp_read_mod.F90
!  03 Dec 2015 - R. Yantosca - Now open file only once per day
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine

    ! SAVEd scalars
    LOGICAL, SAVE      :: first = .TRUE.           ! First time reading data?
                                    
    ! Arrays                                 
    INTEGER            :: st3d(3), ct3d(3)         ! Start & count indices
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Q2(IIPAR,JJPAR      )    ! 2D temporary data arrray
    REAL*4             :: Q3(IIPAR,JJPAR,LGLOB)    ! 3D temporary data arrray

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! Name of this routine (for error printout)
    caller  = "MERRA2_READ_I3_1 (merra2_read_mod.F90)"

    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 030000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > 8 ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL Error_Stop( errMsg, caller )
    ENDIF

    !======================================================================
    ! Open the netCDF file only when necessary
    !======================================================================
    IF ( time_index == 1 .or. first ) THEN
       
       ! Replace time & date tokens in the file name
       dir     = TRIM( Input_Opt%MERRA2_DIR )
       CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

       ! Replace time & date tokens in the file name
       nc_file = Get_Resolution_String()
       nc_file = 'MERRA2.YYYYMMDD.I3.' // TRIM( nc_file )
       CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

       ! Construct complete file path
       nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )

       ! Open netCDF file
       CALL NcOp_Rd( fI3_1, TRIM( nc_file ) )

       ! Read the dimensions from the netCDF file
       CALL NcGet_DimLen( fI3_1, 'lon',   X )
       CALL NcGet_DimLen( fI3_1, 'lat',   Y )
       CALL NcGet_DimLen( fI3_1, 'lev',   Z )
       CALL NcGet_DimLen( fI3_1, 'time',  T )

       ! Make sure the dimensions of the file are valid
       CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                              time=T, time_expected=8, caller=caller )

       ! Reset first-time flag
       first = .FALSE.
    ENDIF

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    !-------------------------------------------------
    ! Read 3D data (2D spatial + 1D time )
    !-------------------------------------------------

    ! netCDF start & count indices
    st3d   = (/ 1,     1,     time_index /)
    ct3d   = (/ IIPAR, JJPAR, 1          /)

    ! Read PS
    v_name = "PS"
    CALL NcRd( Q2, fI3_1, TRIM(v_name), st3d, ct3d )
    State_Met%PS1_WET = Q2

    !-------------------------------------------------
    ! Read 4D data (3D spatial + 1D time)
    !-------------------------------------------------

    ! netCDF start + count indices
    st4d   = (/ 1,     1,     1,     time_index /)
    ct4d   = (/ IIPAR, JJPAR, LGLOB, 1          /)

    !----------------------------------------------------------------
    ! Prior to 2/3/12:
    ! For now, skip reading Potential Vorticity (bmy, 2/3/12)
    !! Read PV
    !v_name = "PV"
    !CALL NcRd( Q3, fI3_1, TRIM(v_name), st4d, ct4d )
    !CALL Transfer_3d( Q3, State_Met%PV )
    !----------------------------------------------------------------

    ! Read QV
    v_name = "QV"
    CALL NcRd( Q3, fI3_1, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q3, State_Met%SPHU1 )

    ! Read T
    v_name = "T"
    CALL NcRd( Q3, fI3_1, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q3, State_Met%TMPU1 )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 4  MERRA2 I3     met fields for ', a )

    !-------------------------------------------------
    ! Unit conversions & special handling
    !-------------------------------------------------
    WHERE ( State_Met%SPHU1 < 0e+0_fp ) 

       ! NOTE: Now set negative Q to a small positive # 
       ! instead of zero, so as not to blow up logarithms
       State_Met%SPHU1 = 1e-32_fp

    ELSEWHERE

       ! Convert MERRA2 specific humidity from [kg/kg] to [g/kg]
       State_Met%SPHU1 = State_Met%SPHU1 * 1000e+0_fp

    ENDWHERE

    ! Initialize State_Met%T to State_Met%TMPU1 and State_Met%SPHU to
    ! State_Met%SPHU1.  After all other met field reads (merra2_read_i3_2)
    ! we will interpolate State_Met%T from the values of State_Met vars 
    ! TMPU1 and TMPU2 and State_Met%SPHU from the values of State_Met vars 
    ! SPHU1 and SPHU2.
    State_Met%T         = State_Met%TMPU1
    State_Met%SPHU      = State_Met%SPHU1

    ! Convert PS1_WET from [Pa] to [hPa]
    State_Met%PS1_WET = State_Met%PS1_WET * 1e-2_fp

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! If it's the last time slice, then close the netCDF file
    ! and set the file ID to -1 to indicate that it's closed.
    IF ( time_index == 8 ) THEN
       CALL NcCl( fI3_1 )
       fI3_1 = -1
    ENDIF

    ! Increment the # of times I3 fields have been read
    CALL Set_Ct_I3( INCREMENT=.TRUE. )

#if defined( BPCH_DIAG )
    ! ND66 diagnostic: T1, QV1 met fields
    IF ( ND66 > 0 ) THEN
       AD66(:,:,1:LD66,3) = AD66(:,:,1:LD66,3) + State_Met%TMPU1(:,:,1:LD66)
       AD66(:,:,1:LD66,4) = AD66(:,:,1:LD66,4) + State_Met%SPHU1(:,:,1:LD66)
    ENDIF
#endif

  END SUBROUTINE Merra2_Read_I3_1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Merra2_Read_I3_2
!
! !DESCRIPTION: Routine to read variables and attributes from a MERRA2
!  met fields file containing 3-hr instantaneous (I3) data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Merra2_Read_I3_2( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
! 
    INTEGER,        INTENT(IN)    :: YYYYMMDD   ! GMT date in YYYY/MM/DD format
    INTEGER,        INTENT(IN)    :: HHMMSS     ! GMT time in hh:mm:ss   format
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
!
! !REMARKS:
!  Even though the netCDF file is self-describing, the MERRA2 data, 
!  dimensions, and units are pre-specified according to the GMAO MERRA2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
!  MERRA2 pressure quantities are stored on disk with units of Pa.  For now
!  we will convert to hPa for compatibility with GEOS-Chem.  But in the near
!  future we will probably recode GEOS-Chem to use the native units, which
!  facilitate GEOS-Chem HP development.
!
! !REVISION HISTORY:
!  12 Aug 2015 - R. Yantosca - Initial version, based on geosfp_read_mod.F90
!  03 Dec 2015 - R. Yantosca - Now open file only once per day
!  20 Sep 2016 - R. Yantosca - Bug fix: FIRST must be declared as LOGICAL

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine

    ! SAVEd scalars
    LOGICAL, SAVE      :: first = .TRUE.           ! First time reading data?
                                    
    ! Arrays                                 
    INTEGER            :: st3d(3), ct3d(3)         ! Start & count indices
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Q2(IIPAR,JJPAR      )    ! 2D temporary data arrray
    REAL*4             :: Q3(IIPAR,JJPAR,LGLOB)    ! 3D temporary data arrray

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! Name of this routine (for error printout)
    caller  = "MERRA2_READ_I3_2 (merra2_read_mod.F90)"

    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 030000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > 8 ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL Error_Stop( errMsg, caller )
    ENDIF

    !======================================================================
    ! Open the netCDF file only when necessary
    !======================================================================
    IF ( time_index == 1 .or. first ) THEN
       
       ! Replace time & date tokens in the file name
       dir     = TRIM( Input_Opt%MERRA2_DIR )
       CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

       ! Replace time & date tokens in the file name
       nc_file = Get_Resolution_String()
       nc_file = 'MERRA2.YYYYMMDD.I3.' // TRIM( nc_file )
       CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

       ! Construct complete file path
       nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )

       ! Open netCDF file
       CALL NcOp_Rd( fI3_2, TRIM( nc_file ) )

       ! Read the dimensions from the netCDF file
       CALL NcGet_DimLen( fI3_2, 'lon',   X )
       CALL NcGet_DimLen( fI3_2, 'lat',   Y )
       CALL NcGet_DimLen( fI3_2, 'lev',   Z )
       CALL NcGet_DimLen( fI3_2, 'time',  T )

       ! Make sure the dimensions of the file are valid
       CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                              time=T, time_expected=8, caller=caller )

       ! Reset the first-time flag
       first = .FALSE.
    ENDIF

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    !-------------------------------------------------
    ! Read 3D data (2D spatial + 1D time )
    !-------------------------------------------------

    ! netCDF start & count indices
    st3d   = (/ 1,     1,     time_index /)
    ct3d   = (/ IIPAR, JJPAR, 1          /)

    ! Read PS
    v_name = "PS"
    CALL NcRd( Q2, fI3_2, TRIM(v_name), st3d, ct3d )
    State_Met%PS2_WET = Q2

    !-------------------------------------------------
    ! Read 4D data (3D spatial + 1D time)
    !-------------------------------------------------

    ! netCDF start + count indices
    st4d   = (/ 1,     1,     1,     time_index /)
    ct4d   = (/ IIPAR, JJPAR, LGLOB, 1          /)

    !----------------------------------------------------------------
    ! Prior to 2/3/12:
    ! For now, skip reading Potential Vorticity (bmy, 2/3/12)
    !! Read PV
    !v_name = "PV"
    !CALL NcRd( Q3, fI3_2, TRIM(v_name), st4d, ct4d )
    !CALL Transfer_3d( Q3, State_Met%PV )
    !----------------------------------------------------------------

    ! Read QV
    v_name = "QV"
    CALL NcRd( Q3, fI3_2, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q3, State_Met%SPHU2 )

    ! Read T
    v_name = "T"
    CALL NcRd( Q3, fI3_2, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q3, State_Met%TMPU2 )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 4  MERRA2 I3     met fields for ', a )

    !-------------------------------------------------
    ! Unit conversions & special handling
    !-------------------------------------------------
    WHERE ( State_Met%SPHU2 < 0d0 ) 

       ! NOTE: Now set negative Q to a small positive # 
       ! instead of zero, so as not to blow up logarithms
       State_Met%SPHU2 = 1d-32

    ELSEWHERE

       ! Convert MERRA2 specific humidity from [kg/kg] to [g/kg]
       State_Met%SPHU2 = State_Met%SPHU2 * 1000d0

    ENDWHERE

    ! Convert PS2_WET from [Pa] to [hPa]
    State_Met%PS2_WET = State_Met%PS2_WET * 1e-2_fp

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! If it's the last time slice, then close the netCDF file
    ! and set the file ID to -1 to indicate that it's closed.
    IF ( time_index == 8 ) THEN
       CALL NcCl( fI3_2 )
       fI3_2 = -1
    ENDIF

    ! Increment the # of times I3 fields have been read
    CALL Set_Ct_I3( INCREMENT=.TRUE. )

#if defined( BPCH_DIAG )
    ! ND66 diagnostic: T2, QV2 met fields
    IF ( ND66 > 0 ) THEN
       AD66(:,:,1:LD66,3) = AD66(:,:,1:LD66,3) + State_Met%TMPU2(:,:,1:LD66)
       AD66(:,:,1:LD66,4) = AD66(:,:,1:LD66,4) + State_Met%SPHU2(:,:,1:LD66)
    ENDIF
#endif

  END SUBROUTINE Merra2_Read_I3_2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Merra2_Read
!
! !DESCRIPTION: Closes any open netCDF files at the end of a simulation.
!  This can occur if the simulation ends at a time other than 00:00 GMT.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Merra2_Read()
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! Close any netCDF file ID's that are still open
    !=======================================================================
    IF ( fA1 > 0 ) THEN
       CALL NcCl( fA1 )
       WRITE( 6, 100 ) 'A1    '
    ENDIF

    IF ( fA3cld > 0 ) THEN
       CALL NcCl( fA3cld  )
       WRITE( 6, 100 ) 'A3cld '
    ENDIF

    IF ( fA3dyn > 0 ) THEN
       CALL NcCl( fA3dyn  )
       WRITE( 6, 100 ) 'A3dyn '
    ENDIF

    IF ( fA3mstC > 0 ) THEN
       CALL NcCl( fA3mstC )
       WRITE( 6, 100 ) 'A3mstC'
    ENDIF

    IF ( fA3mstE > 0 ) THEN
       CALL NcCl( fA3mstE )
       WRITE( 6, 100 ) 'A3mstE'
    ENDIF

    IF ( fI3_1 > 0 ) THEN
       CALL NcCl( fI3_1   )
       WRITE( 6, 100 ) 'I3_1  '
    ENDIF

    IF ( fI3_2 > 0 ) THEN
       CALL NcCl( fI3_2   )
       WRITE( 6, 100 ) 'I3_2  '
    ENDIF

    ! Format strings
100 FORMAT( '     - Closing the MERRA2 ', a6, ' file at end of simulation' )

  END SUBROUTINE Cleanup_Merra2_Read
!EOC
END MODULE Merra2_Read_Mod
