!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: geosfp_read_mod
!
! !DESCRIPTION: Module GEOSFP\_READ\_MOD contains subroutines for reading the 
!  GEOS-FP data from disk (in netCDF format).
!\\
!\\
! !INTERFACE: 
!
MODULE GeosFp_Read_Mod
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
  USE CMN_GCTM_MOD                        ! Physical constants
  USE CMN_DIAG_MOD                        ! Diagnostic arrays & counters
  USE DIAG_MOD,      ONLY : AD66          ! Array for ND66 diagnostic  
  USE DIAG_MOD,      ONLY : AD67          ! Array for ND67 diagnostic
  USE ERROR_MOD,     ONLY : ERROR_STOP    ! Stop w/ error message
  USE TIME_MOD                            ! Date & time routines
  USE TRANSFER_MOD                        ! Routines for casting 

  IMPLICIT NONE
  PRIVATE

# include "netcdf.inc"                    ! Include file for netCDF library
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Check_Dimensions
  PRIVATE :: GeosFp_Read_A3cld
  PRIVATE :: GeosFp_Read_A3dyn
  PRIVATE :: GeosFp_Read_A3mstC
  PRIVATE :: GeosFp_Read_A3mstE
  PRIVATE :: Get_Resolution_String
!
! !PUBLIC MEMBER FUNCTIONS:
! 
  PUBLIC  :: GeosFp_Read_CN
  PUBLIC  :: GeosFp_Read_A1
  PUBLIC  :: GeosFp_Read_A3
  PUBLIC  :: GeosFp_Read_I3_1
  PUBLIC  :: GeosFp_Read_I3_2
!
! !REMARKS:
!  Assumes that you have a netCDF library (either v3 or v4) installed on 
!  your system. 
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  03 Feb 2012 - R. Yantosca - Add Geos57_Read_A3 wrapper function
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Add function Get_Resolution_String
!  05 Apr 2012 - R. Yantosca - Convert units for specific humidity properly
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  11 Apr 2013 - R. Yantosca - Now pass directory fields via Input_Opt
!  26 Sep 2013 - R. Yantosca - Renamed to geosfp_read_mod.F90
!  14 Jan 2014 - R. Yantosca - Remove "define GEOS572_FILES #ifdef blocks
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
! !IROUTINE: get_resolution_string
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
!  10 Feb 2012 - R. Yantosca - Initial version
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Remove SEAC4RS C-preprocssor switch
!  14 Jan 2014 - R. Yantosca - Now add NESTED_SE option
!EOP
!------------------------------------------------------------------------------
!BOC

#if   defined( GRID4x5 )
    resString = '4x5.nc'
     
#elif defined( GRID2x25 ) 
    resString = '2x25.nc'

#elif defined( GRID1x125 )
    resString = '1x125.nc'

#elif defined( GRID1x1 ) 
    resString = '1x1.nc'

#elif defined( GRID05x0666 )
    resString = '05x0666.nc'

#elif defined( GRID025x03125 ) && defined( NESTED_CH )
    resString = '025x03125.CH.nc'

#elif defined( GRID025x03125 ) && defined( NESTED_EU )
    resString = '025x03125.EU.nc'

#elif defined( GRID025x03125 ) && defined( NESTED_NA )
    resString = '025x03125.NA.nc'

#elif defined( GRID025x03125 ) && defined( NESTED_SE )
    resString = '025x03125.SE.nc'

#elif defined( GRID025x03125 )
    resString = '025x03125.nc'

#endif

  END FUNCTION Get_Resolution_String
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_dimensions
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
! !IROUTINE: geosfp_read_cn
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-FP
!  met fields file containing constant (CN) data.  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFp_Read_CN( Input_Opt, State_Met )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
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
!  This routine was automatically generated by the Perl script ncCodeRead, 
!  and was subsequently hand-edited for compatibility with GEOS-Chem.
!                                                                             .
!  Even though the netCDF file is self-describing, the GEOS-FP data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-FP
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  11 Apr 2013 - R. Yantosca - Now pass directory fields with Input_Opt
!  26 Sep 2013 - R. Yantosca - Renamed to GeosFp_Read_CN
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                ! netCDF file ID
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
    caller  = "GEOSFP_READ_CN (geosfp_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( Input_Opt%GEOS_FP_DIR )
    CALL Expand_Date( dir, 20110101, 000000 )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String() 
    nc_file = 'GEOSFP.YYYYMMDD.CN.' // TRIM( nc_file ) 
    CALL Expand_Date( nc_file, 20110101, 000000 )

    ! Construct complete file path
    nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'time',  T )

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
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%FRLAKE )

    ! Read FRLAND
    v_name = "FRLAND"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%FRLAND )

    ! Read FRLANDIC
    v_name = "FRLANDIC"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%FRLANDIC )
    
    ! Read FROCEAN
    v_name = "FROCEAN"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%FROCEAN )
    
    ! Read PHIS
    v_name = "PHIS"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%PHIS )

    ! Echo info
    stamp = TimeStamp_String( 20110101, 000000 )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 5  GEOS-FP CN     met fields for ', a )

    !======================================================================
    ! Cleanup and quit
    !======================================================================

    ! Convert PHIS from [m2/s2] to [m]
    State_Met%PHIS = State_Met%PHIS / g0

    ! ND67 diagnostic 
    IF ( ND67 > 0 ) THEN
       AD67(:,:,15) = AD67(:,:,15) + State_Met%PHIS  ! Sfc geopotential [m]
    ENDIF

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE GeosFp_Read_CN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geosfp_read_a1
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-FP
!  met fields file containing 1-hr time-averaged (A1) data.  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFp_Read_A1( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
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
!  Even though the netCDF file is self-describing, the GEOS-FP data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-FP
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!                                                                             .
!  Special handling for surface precipitation fields:
!  ---------------------------------------------------------------------------
!  In GEOS-FP (and in MERRA), the PRECTOT etc. surface precipitation
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
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  04 Jan 2013 - M. Payer    - Bug fix: Use State_Met%TSKIN for ND67 surface
!                              skin temperature diagnostic, not State_MET%TS
!  11 Apr 2013 - R. Yantosca - Now pass directory fields with Input_Opt
!  02 Dec 2013 - S. Philip   - Correction for GEOS-FP boundary layer height
!  04 Dec 2013 - R. Yantosca - Now comment out GEOS-FP BL height correction
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                ! netCDF file ID
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
                                             
    ! Arrays                                 
    INTEGER            :: st3d(3), ct3d(3)   ! Start + count, for 3D arrays 
    REAL*4             :: Q(IIPAR,JJPAR)     ! Temporary data arrray

    !======================================================================
    ! Open the netCDF file
    !======================================================================

    ! Skip if we have already read data for this date & time
    IF ( YYYYMMDD == lastDate .and. HHMMSS == lastTime ) THEN
       stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
       WRITE( 6, 20 ) stamp
 20    FORMAT( '     - GEOS-FP A1 met fields for ', a,  &
               ' have been read already'                  ) 
       RETURN
    ENDIF

    ! Name of this routine (for error printout)
    caller  = "GEOSFP_READ_A1 (geosfp_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( Input_Opt%GEOS_FP_DIR )
    CALL Expand_Date( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOSFP.YYYYMMDD.A1.' // TRIM( nc_file )
    CALL Expand_Date( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'time',  T )

    ! Make sure the dimensions of the file are valid
    CALL Check_Dimensions( lon=X,            lat=Y,        time=T,  &
                           time_expected=24, caller=caller         )
    
    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 10000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > T ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 24!' )
       CALL Error_Stop( errMsg, caller )
    ENDIF

    ! netCDF start & count indices
    st3d      = (/ 1,     1,     time_index /)      
    ct3d      = (/ IIPAR, JJPAR, 1          /)

    ! Read ALBEDO
    v_name = "ALBEDO"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%ALBD )

    ! Read CLDTOT
    v_name = "CLDTOT"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%CLDFRC )

    ! Read EFLUX
    v_name = "EFLUX"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%EFLUX )

    ! Read EVAP
    v_name = "EVAP"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%EVAP )

    ! Read FRSEAICE
    v_name = "FRSEAICE"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%FRSEAICE )

    ! Read FRSNO
    v_name = "FRSNO"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%FRSNO )

    ! Read GRN
    v_name = "GRN"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%GRN )

    ! Read GWETROOT
    v_name = "GWETROOT"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%GWETROOT )

    ! Read GWETTOP
    v_name = "GWETTOP"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%GWETTOP )

    ! Read HFLUX from file
    v_name = "HFLUX"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%HFLUX )

    ! Read LAI
    v_name = "LAI"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%LAI )

    ! Read LWI
    v_name = "LWI"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%LWI )

    ! Read LWGNT 
    v_name = "LWGNT"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%RADLWG )

    !-----------------------------------------------------------------------
    ! Comment this out for now, this field isn't needed (bmy, 2/2/12)
    !! Read LWTUP
    !v_name = "LWTUP"
    !CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    !CALL Transfer_2d( Q, State_Met%FRLAKE )
    !-----------------------------------------------------------------------

    ! Read PARDF
    v_name = "PARDF"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%PARDF )

    ! Read PARDR
    v_name = "PARDR"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%PARDR )

    ! Read PBLH
    v_name = "PBLH"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%PBLH )

    ! Read PRECANV
    v_name = "PRECANV"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%PRECANV )

    ! Read PRECCON
    v_name = "PRECCON"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%PRECCON )

    ! Read PRECLSC
    v_name = "PRECLSC"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%PRECLSC )

    ! Read PRECSNO
    v_name = "PRECSNO"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%PRECSNO )

    ! Read PRECTOT
    v_name = "PRECTOT"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%PRECTOT )

    !-----------------------------------------------------------------------
    ! Comment this out for now, this field isn't needed (bmy, 2/2/12)
    !! Read QV2M
    !v_name = "QV2M"
    !CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    !CALL Transfer_2d( Q, State_Met%QV2M )
    !-----------------------------------------------------------------------

    ! Read SEAICE00
    v_name = "SEAICE00"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%SEAICE00 )

    ! Read SEAICE10
    v_name = "SEAICE10"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%SEAICE10 )

    ! Read SEAICE20
    v_name = "SEAICE20"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%SEAICE20 )

    ! Read SEAICE30
    v_name = "SEAICE30"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%SEAICE30 )

    ! Read SEAICE40
    v_name = "SEAICE40"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%SEAICE40 )

    ! Read SEAICE50
    v_name = "SEAICE50"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%SEAICE50 )

    ! Read SEAICE60 
    v_name = "SEAICE60"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%SEAICE60 )

    ! Read SEAICE70
    v_name = "SEAICE70"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%SEAICE70 )

    ! Read SEAICE80
    v_name = "SEAICE80"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%SEAICE80 )

    ! Read SEAICE90
    v_name = "SEAICE90"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%SEAICE90 )

    ! Read SLP
    v_name = "SLP"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%SLP )

    ! Read SNODP
    v_name = "SNODP"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%SNODP )

    ! Read SNOMAS
    v_name = "SNOMAS"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%SNOMAS )

    ! Read SWGDN
    v_name = "SWGDN"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%RADSWG )

    ! Read TO3
    v_name = "TO3"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%TO3 )

    ! Read TROPPT
    v_name = "TROPPT"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%TROPP )

    ! Read TS
    v_name = "TS"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%TSKIN )

    ! Read T2M
    v_name = "T2M"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%TS )

    ! Read U10M
    v_name = "U10M"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%U10M )

    ! Read USTAR
    v_name = "USTAR"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%USTAR )

    ! Read V10M
    v_name = "V10M"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met% V10M )

    ! Read Z0M
    v_name = "Z0M"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, State_Met%Z0 )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp                    
 10 FORMAT( '     - Found all 44 GEOS-FP A1     met fields for ', a )

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! Close netCDF file
    CALL NcCl( fId )

    ! Increment the # of times A1 fields are read from disk
    CALL Set_Ct_A1( INCREMENT=.TRUE. )
      
    ! Convert surface precip fields from [kg/m2/s] --> [mm/day]
    State_Met%PRECANV = State_Met%PRECANV * 86400d0
    State_Met%PRECCON = State_Met%PRECCON * 86400d0
    State_Met%PRECLSC = State_Met%PRECLSC * 86400d0
    State_Met%PRECTOT = State_Met%PRECTOT * 86400d0
    
    ! ND67 diagnostic: surface fields
    IF ( ND67 > 0 ) THEN
       AD67(:,:,1 ) = AD67(:,:,1 ) + State_Met%HFLUX    ! Sens heat flux [W/m2]
       AD67(:,:,2 ) = AD67(:,:,2 ) + State_Met%RADSWG   ! SW rad @ sfc [W/m2]
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

    ! Save date & time for next iteration
    lastDate = YYYYMMDD
    lastTime = HHMMSS

  END SUBROUTINE GeosFp_Read_A1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geosfp_read_a3
!
! !DESCRIPTION: Convenience wrapper for the following routines which read
!  3-hour time averaged data from disk:
! \begin{itemize}
! \item GeosFp\_Read\_A3cld
! \item GeosFp\_Read\_A3dyn
! \item GeosFp\_Read\_A3mstC
! \item GeosFp\_Read\_A3mstE
! \end{itemize}
!
! !INTERFACE:
!
  SUBROUTINE GeosFp_Read_A3( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
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
!  30 Jan 2012 - R. Yantosca - Initial version
!  11 Apr 2013 - R. Yantosca - Now pass directory fields with Input_Opt
!  26 Sep 2013 - R. Yantosca - Renamed to GeosFp_Read_A3
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
 20    FORMAT( '     - GEOS-FP A3 met fields for ', a,  &
               ' have been read already'                  ) 
       RETURN
    ENDIF

    ! Save date & time for next iteration
    lastDate = YYYYMMDD
    lastTime = HHMMSS

    ! Read all the diffeent A3 files
    CALL GeosFp_Read_A3cld ( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
    CALL GeosFp_Read_A3dyn ( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
    CALL GeosFp_Read_A3mstC( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
    CALL GeosFp_Read_A3mstE( YYYYMMDD, HHMMSS, Input_Opt, State_Met )

    !======================================================================
    ! Cleanup and quit
    !======================================================================

    ! Increment the # of times that A3 fields have been read
    CALL Set_Ct_A3( INCREMENT=.TRUE. )

    ! Save date & time for next iteration
    lastDate = YYYYMMDD
    lastTime = HHMMSS

  END SUBROUTINE GeosFp_Read_A3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geosfp_read_a3cld
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-FP
!  met fields file containing 3-hr time-averaged (A3) data (cloud fields).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFp_Read_A3cld( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
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
!  Even though the netCDF file is self-describing, the GEOS-FP data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-FP
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
!  05 Apr 2012 - R. Yantosca - Fixed bug: TAUCLI was overwritten w/ TAUCLW
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  11 Apr 2013 - R. Yantosca - Now pass directory fields with Input_Opt
!  26 Sep 2013 - R. Yantosca - Renamed to GeosFp_Read_A3Cld
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                      ! netCDF file ID
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine
                                             
    ! Arrays                                 
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Q(IIPAR,JJPAR,LGLOB)     ! Temporary data arrray

    !======================================================================
    ! Open the netCDF file
    !======================================================================

    ! Name of this routine (for error printout)
    caller  = "GEOSFP_READ_A3cld (geosfp_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( Input_Opt%GEOS_FP_DIR )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOSFP.YYYYMMDD.A3cld.' // TRIM( nc_file )
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'lev',   Z )
    CALL NcGet_DimLen( fId, 'time',  T )

    ! Make sure the dimensions of the file are valid
    CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                           time=T, time_expected=8, caller=caller )

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 030000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > T ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL ERROR_STOP( errMsg, caller )
    ENDIF

    ! netCDF start & count indices
    st4d      = (/ 1,     1,     1,     time_index /)      
    ct4d      = (/ IIPAR, JJPAR, LGLOB, 1          /)

    ! Read CLOUD
    v_name = "CLOUD"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL TRANSFER_A6( Q, State_Met%CLDF )
    
    ! Read OPTDEPTH
    v_name = "OPTDEPTH"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL TRANSFER_A6( Q, State_Met%OPTDEP )

    ! Read QI
    v_name = "QI"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%QI )

    ! Read QL
    v_name = "QL"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%QL )

    ! Read TAUCLI
    v_name = "TAUCLI"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%TAUCLI )

    ! Read TAUCLW
    v_name = "TAUCLW"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%TAUCLW )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 6  GEOS-FP A3cld  met fields for ', a )

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE GeosFp_Read_A3cld
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geosfp_read_a3dyn
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-FP
!  met fields file containing 3-hr time-averaged (A3) data (dynamics fields).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFp_Read_A3dyn( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput  
    USE GIGC_State_Met_Mod, ONLY : MetState
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
!  Even though the netCDF file is self-describing, the GEOS-FP data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-FP
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  11 Apr 2013 - R. Yantosca - Now pass directories with Input_Opt
!  26 Sep 2013 - R. Yantosca - Renamed to GeosFp_Read_A3dyn
!  15 Nov 2013 - R. Yantosca - Now convert RH from [1] to [%], in order
!                              to be consistent with GEOS-Chem convention
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                      ! netCDF file ID
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine
                                             
    ! Arrays                                 
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Q (IIPAR,JJPAR,LGLOB  )  ! Temporary data arrray
    REAL*4             :: Qe(IIPAR,JJPAR,LGLOB+1)  ! Temporary data arrray

    !======================================================================
    ! Open the netCDF file
    !======================================================================

    ! Name of this routine (for error printout)
    caller  = "GEOSFP_READ_A3dyn (geosfp_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( Input_Opt%GEOS_FP_DIR )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOSFP.YYYYMMDD.A3dyn.' // TRIM( nc_file )
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'lev',   Z )
    CALL NcGet_DimLen( fId, 'time',  T )

    ! Make sure the dimensions of the file are valid
    CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                           time=T, time_expected=8, caller=caller )

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 030000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > T ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL ERROR_STOP( errMsg, 'GEOS57_READ_A1 (geos57_read_mod.F90)' )
    ENDIF

    !--------------------------------
    ! Read data on level centers
    !--------------------------------

    ! netCDF start & count indices
    st4d      = (/ 1,     1,     1,     time_index /)      
    ct4d      = (/ IIPAR, JJPAR, LGLOB, 1          /)

    ! Read DTRAIN
    v_name = "DTRAIN"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%DTRAIN )

    !! Read OMEGA  from file
    !v_name = "OMEGA"
    !CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    !CALL Transfer_3d( Q, State_Met%OMEGA )

    ! Read RH
    v_name = "RH"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%RH )

    ! Read U
    v_name = "U"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%U )

    ! Read V
    v_name = "V"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%V )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 6  GEOS-FP A3dyn  met fields for ', a )

    !======================================================================
    ! Unit conversions, diagnostics, cleanup, and quit
    !======================================================================
    
    ! Convert RH from [1] to [%]
    State_Met%RH = State_Met%RH * 100d0

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

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE GeosFp_Read_A3dyn
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geosfp_read_a3mstc
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-FP
!  met fields file containing 3-hr time-averaged (A3) data (moist fields,
!  saved on level centers).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFp_Read_A3mstC( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
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
!  Even though the netCDF file is self-describing, the GEOS-FP data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-FP
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  11 Apr 2013 - R. Yantosca - Now pass directory fields with Input_Opt
!  26 Sep 2013 - R. Yantosca - Renamed to GeosFp_Read_A3mstC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                      ! netCDF file ID
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine
                                    
    ! Arrays                                 
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Q (IIPAR,JJPAR,LGLOB)    ! Temporary data arrray

    !======================================================================
    ! Open the netCDF file
    !======================================================================

    ! Name of this routine (for error printout)
    caller  = "GEOSFP_READ_A3mstC (geosfp_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( Input_Opt%GEOS_FP_DIR )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOSFP.YYYYMMDD.A3mstC.' // TRIM( nc_file )
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'lev',   Z )
    CALL NcGet_DimLen( fId, 'time',  T )

    ! Make sure the dimensions of the file are valid
    CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                           time=T, time_expected=8, caller=caller )

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 030000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > T ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL Error_Stop( errMsg, caller )
    ENDIF

    ! netCDF start & count indices
    st4d      = (/ 1,     1,     1,     time_index /)      
    ct4d      = (/ IIPAR, JJPAR, LGLOB, 1          /)

    ! Read DQRCU  from file
    v_name = "DQRCU"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%DQRCU )

    ! Read DQRLSAN
    v_name = "DQRLSAN"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%DQRLSAN )

    ! Read REEVAPCN
    v_name = "REEVAPCN"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%REEVAPCN )

    ! Read  from file
    v_name = "REEVAPLS"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, State_Met%REEVAPLS )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 4  GEOS-FP A3mstC met fields for ', a )

    !======================================================================
    ! Cleanup and quit
    !======================================================================

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE GeosFp_read_A3mstC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geosfp_read_a3mste
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-FP
!  met fields file containing 3-hr time-averaged (A3) data (moist fields,
!  saved on level edges).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFp_read_A3mstE( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
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
!  Even though the netCDF file is self-describing, the GEOS-FP data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-FP
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  11 Apr 2013 - R. Yantosca - Now pass directory fields with Input_Opt
!  26 Sep 2013 - R. Yantosca - Renamed to GeosFp_Read_A3mstE
!  26 Sep 2013 - R. Yantosca - Now read CMFMC from GEOSFP*.nc files
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                      ! netCDF file ID
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine
                                             
    ! Arrays                                 
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Qe(IIPAR,JJPAR,LGLOB+1)  ! Temporary data arrray

    !======================================================================
    ! Open the netCDF file
    !======================================================================

    ! Name of this routine (for error printout)
    caller  = "GEOSFP_READ_A3mstE (geosfp_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( Input_Opt%GEOS_FP_DIR )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOSFP.YYYYMMDD.A3mstE.' // TRIM( nc_file )
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'lev',   Z )
    CALL NcGet_DimLen( fId, 'time',  T )

    ! Make sure the dimensions of the file are valid
    CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z-1,       &
                           time=T, time_expected=8, caller=caller )

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 030000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > T ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL ERROR_STOP( errMsg, caller )
    ENDIF

    ! netCDF start & count indices
    st4d      = (/ 1,     1,     1,       time_index /)      
    ct4d      = (/ IIPAR, JJPAR, LGLOB+1, 1          /)

    ! Read CMFMC (only in GEOSFP*.nc files)
    v_name = "CMFMC"
    CALL NcRd( Qe, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe, State_Met%CMFMC )

    ! Read PFICU
    v_name = "PFICU"
    CALL NcRd( Qe, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe, State_Met%PFICU )

    ! Read PFILSAN
    v_name = "PFILSAN"
    CALL NcRd( Qe, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe, State_Met%PFILSAN )

    ! Read PFLCU
    v_name = "PFLCU"
    CALL NcRd( Qe, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe, State_Met%PFLCU )

    ! Read  from file
    v_name = "PFLLSAN"
    CALL NcRd( Qe, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe, State_Met%PFLLSAN )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 4  GEOS-FP A3mstE met fields for ', a )

    !=================================================================
    ! Diagnostics, cleanup and quit
    !=================================================================

    ! ND66 diagnostic: CMFMC met field
    IF ( ND66 > 0 ) THEN
       AD66(:,:,1:LD66,5) = AD66(:,:,1:LD66,5) + State_Met%CMFMC(:,:,1:LD66)
    ENDIF

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE GeosFp_Read_A3mstE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geosfp_read_I3_1
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-FP
!  met fields file containing 3-hr instantaneous (I3) data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFp_Read_I3_1( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE DAO_MOD,            ONLY : T_FULLGRID
    USE DAO_MOD,            ONLY : T_FULLGRID_1
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
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
!  Even though the netCDF file is self-describing, the GEOS-FP data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-FP
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
!  05 Apr 2012 - R. Yantosca - Now convert QV1 from [kg/kg] to [g/kg]
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  11 Apr 2013 - R. Yantosca - Now pass directory fields with Input_Opt
!  06 Sep 2013 - R. Yantosca - Bug fix: we need to initialize State_Met%T
!                              with State_Met%TMPU1 to avoid errors.  The
!                              State_Met%T field will be set again in INTERP.
!  26 Sep 2013 - R. Yantosca - Renamed to GeosFp_Read_I3_1
!  29 Oct 2013 - R. Yantosca - Now read T_FULLGRID_1 for offline simulations
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                      ! netCDF file ID
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine
                                    
    ! Arrays                                 
    INTEGER            :: st3d(3), ct3d(3)         ! Start & count indices
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Q2(IIPAR,JJPAR      )    ! 2D temporary data arrray
    REAL*4             :: Q3(IIPAR,JJPAR,LGLOB)    ! 3D temporary data arrray

    !======================================================================
    ! Open the netCDF file
    !======================================================================

    ! Name of this routine (for error printout)
    caller  = "GEOSFP_READ_I3_1 (geosfp_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( Input_Opt%GEOS_FP_DIR )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOSFP.YYYYMMDD.I3.' // TRIM( nc_file )
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )

    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'lev',   Z )
    CALL NcGet_DimLen( fId, 'time',  T )

    ! Make sure the dimensions of the file are valid
    CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                           time=T, time_expected=8, caller=caller )

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 030000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > T ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL Error_Stop( errMsg, caller )
    ENDIF
    
    !-------------------------------------------------
    ! Read 3D data (2D spatial + 1D time )
    !-------------------------------------------------

    ! netCDF start & count indices
    st3d   = (/ 1,     1,     time_index /)
    ct3d   = (/ IIPAR, JJPAR, 1          /)

    ! Read PS
    v_name = "PS"
    CALL NcRd( Q2, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q2, State_Met%PS1 )

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
    !CALL NcRd( Q3, fId, TRIM(v_name), st4d, ct4d )
    !CALL Transfer_3d( Q3, !State_Met%PV )
    !----------------------------------------------------------------

    ! Read QV
    v_name = "QV"
    CALL NcRd( Q3, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q3, State_Met%SPHU1 )

    ! Read T
    v_name = "T"
    CALL NcRd( Q3, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q3, State_Met%TMPU1 )

    ! Also archive T_FULLGRID_1 for offline "specialty" simulations
    IF ( Input_Opt%ITS_A_SPECIALTY_SIM ) THEN
       T_FULLGRID_1 = Q3
       T_FULLGRID   = Q3   ! Also need to initialize T_FULLGRID
    ENDIF

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 4  GEOS-FP I3     met fields for ', a )

    !-------------------------------------------------
    ! Unit conversions & special handling
    !-------------------------------------------------
    WHERE ( State_Met%SPHU1 < 0d0 ) 

       ! NOTE: Now set negative Q to a small positive # 
       ! instead of zero, so as not to blow up logarithms
       State_Met%SPHU1 = 1d-32

    ELSEWHERE

       ! Convert GEOS-FP specific humidity from [kg/kg] to [g/kg]
       State_Met%SPHU1 = State_Met%SPHU1 * 1000d0

    ENDWHERE

    ! For now, copy State_Met%TMPU1 to State_Met%T.  At the next met field 
    ! read, we will State_Met%T from the values of State_Met%TMPU1 and
    ! State_Met%TMPU2. (bmy, 9/6/13)
    State_Met%T = State_Met%TMPU1

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! Close netCDF file
    CALL NcCl( fId )

    ! Increment the # of times I3 fields have been read
    CALL Set_Ct_I3( INCREMENT=.TRUE. )

    ! ND66 diagnostic: T1, QV1 met fields
    IF ( ND66 > 0 ) THEN
       AD66(:,:,1:LD66,3) = AD66(:,:,1:LD66,3) + State_Met%TMPU1(:,:,1:LD66)
       AD66(:,:,1:LD66,4) = AD66(:,:,1:LD66,4) + State_Met%SPHU1(:,:,1:LD66)
    ENDIF

  END SUBROUTINE GeosFp_Read_I3_1
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geosfp_read_I3_2
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-FP
!  met fields file containing 3-hr instantaneous (I3) data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GeosFp_Read_I3_2( YYYYMMDD, HHMMSS, Input_Opt, State_Met )
!
! !USES:
!
    USE DAO_MOD,            ONLY : T_FULLGRID_2
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
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
!  Even though the netCDF file is self-describing, the GEOS-FP data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-FP
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
!  05 Apr 2012 - R. Yantosca - Now convert QV2 from [kg/kg] to [g/kg]
!  09 Nov 2012 - M. Payer    - Copy all met fields to the State_Met derived type
!                              object
!  15 Nov 2012 - R. Yantosca - Now replace dao_mod.F arrays with State_Met
!  11 Apr 2013 - R. Yantosca - Now pass directory fields with Input_Opt
!  26 Sep 2013 - R. Yantosca - Rename to GeosFp_Read_I3_2
!  29 Oct 2013 - R. Yantosca - Now read T_FULLGRID_2 for offline simulations
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                      ! netCDF file ID
    INTEGER            :: I, J, L                  ! Loop indices
    INTEGER            :: X, Y, Z, T               ! netCDF file dimensions
    INTEGER            :: time_index               ! Read this slice of data
    CHARACTER(LEN=16)  :: stamp                    ! Time and date stamp
    CHARACTER(LEN=255) :: nc_file                  ! netCDF file name
    CHARACTER(LEN=255) :: v_name                   ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                      ! Data directory path
    CHARACTER(LEN=255) :: errMsg                   ! Error message
    CHARACTER(LEN=255) :: caller                   ! Name of this routine
                                    
    ! Arrays                                 
    INTEGER            :: st3d(3), ct3d(3)         ! Start & count indices
    INTEGER            :: st4d(4), ct4d(4)         ! Start & count indices
    REAL*4             :: Q2(IIPAR,JJPAR      )    ! 2D temporary data arrray
    REAL*4             :: Q3(IIPAR,JJPAR,LGLOB)    ! 3D temporary data arrray

    !======================================================================
    ! Open the netCDF file
    !======================================================================

    ! Name of this routine (for error printout)
    caller  = "GEOSFP_READ_I3_2 (geosfp_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( Input_Opt%GEOS_FP_DIR )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOSFP.YYYYMMDD.I3.' // TRIM( nc_file )
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( Input_Opt%DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )

    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'lev',   Z )
    CALL NcGet_DimLen( fId, 'time',  T )

    ! Make sure the dimensions of the file are valid
    CALL Check_Dimensions( lon=X,  lat=Y,           lev=Z,         &
                           time=T, time_expected=8, caller=caller )

    !======================================================================
    ! Read data from the netCDF file
    !======================================================================
    
    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 030000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > T ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
       CALL Error_Stop( errMsg, caller )
    ENDIF
    
    !-------------------------------------------------
    ! Read 3D data (2D spatial + 1D time )
    !-------------------------------------------------

    ! netCDF start & count indices
    st3d   = (/ 1,     1,     time_index /)
    ct3d   = (/ IIPAR, JJPAR, 1          /)

    ! Read PS
    v_name = "PS"
    CALL NcRd( Q2, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q2, State_Met%PS2 )

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
    !CALL NcRd( Q3, fId, TRIM(v_name), st4d, ct4d )
    !CALL Transfer_3d( Q3, State_Met%PV )
    !----------------------------------------------------------------

    ! Read QV
    v_name = "QV"
    CALL NcRd( Q3, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q3, State_Met%SPHU2 )

    ! Read T
    v_name = "T"
    CALL NcRd( Q3, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q3, State_Met%TMPU2 )

    ! Also archive T_FULLGRID_1 for offline "specialty" simulations
    IF ( Input_Opt%ITS_A_SPECIALTY_SIM ) THEN
       T_FULLGRID_2 = Q3
    ENDIF

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 4  GEOS-FP I3     met fields for ', a )

    !-------------------------------------------------
    ! Unit conversions & special handling
    !-------------------------------------------------
    WHERE ( State_Met%SPHU2 < 0d0 ) 

       ! NOTE: Now set negative Q to a small positive # 
       ! instead of zero, so as not to blow up logarithms
       State_Met%SPHU2 = 1d-32

    ELSEWHERE

       ! Convert GEOS-FP specific humidity from [kg/kg] to [g/kg]
       State_Met%SPHU2 = State_Met%SPHU2 * 1000d0

    ENDWHERE

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! Close netCDF file
    CALL NcCl( fId )

    ! Increment the # of times I3 fields have been read
    CALL Set_Ct_I3( INCREMENT=.TRUE. )

    ! ND66 diagnostic: T2, QV2 met fields
    IF ( ND66 > 0 ) THEN
       AD66(:,:,1:LD66,3) = AD66(:,:,1:LD66,3) + State_Met%TMPU2(:,:,1:LD66)
       AD66(:,:,1:LD66,4) = AD66(:,:,1:LD66,4) + State_Met%SPHU2(:,:,1:LD66)
    ENDIF

  END SUBROUTINE GeosFp_Read_I3_2
!EOC
END MODULE GeosFp_Read_Mod
