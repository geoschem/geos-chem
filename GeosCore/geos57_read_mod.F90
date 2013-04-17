!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: geos57_read_mod
!
! !DESCRIPTION: Module GEOS57\_READ\_MOD contains subroutines for reading the 
!  GEOS-5.7.2 data from disk (in netCDF format).
!\\
!\\
! !INTERFACE: 
!
MODULE Geos57_Read_Mod
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
  USE DIRECTORY_MOD                       ! Directory paths
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
  PRIVATE :: Geos57_Read_A3cld
  PRIVATE :: Geos57_Read_A3dyn
  PRIVATE :: Geos57_Read_A3mstC
  PRIVATE :: Geos57_Read_A3mstE
  PRIVATE :: Get_Resolution_String
!
! !PUBLIC MEMBER FUNCTIONS:
! 
  PUBLIC  :: Geos57_Read_CN
  PUBLIC  :: Geos57_Read_A1
  PUBLIC  :: Geos57_Read_A3
  PUBLIC  :: Geos57_Read_I3_1
  PUBLIC  :: Geos57_Read_I3_2
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
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
! !USES:
!
#   include "define.h"
!
! !RETURN VALUE:
!
    CHARACTER(LEN=255) :: resString
! 
! !REVISION HISTORY:
!  10 Feb 2012 - R. Yantosca - Initial version
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
    resString = '025x03125.SEA4CRS.nc'

#elif defined( GRID025x03125 ) && defined( SEAC4RS )
    resString = '025x03125.SEA4CRS.nc'

#elif defined( GRID025x03125 ) && defined( NESTED_EU )
    resString = '025x03125.EU.nc'

#elif defined( GRID025x03125 ) && defined( NESTED_NA )
    resString = '025x03125.NA.nc'

#elif defined( GRID025x03125 )
    resString = '025x03125.nc'

#endif

  END FUNCTION Get_Resolution_String
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geos57_read_cn
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-5.7.2
!  met fields file containing constant (CN) data.  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Geos57_Read_CN()
!
! !USES:
!
    USE DAO_MOD, ONLY : FRLAKE
    USE DAO_MOD, ONLY : FRLAND
    USE DAO_MOD, ONLY : FRLANDIC
    USE DAO_MOD, ONLY : FROCEAN
    USE DAO_MOD, ONLY : PHIS
!
! !REMARKS:
!  This routine was automatically generated by the Perl script ncCodeRead, 
!  and was subsequently hand-edited for compatibility with GEOS-Chem.
!                                                                             .
!  Even though the netCDF file is self-describing, the GEOS-5.7.2 data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-5.7.2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
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
    caller  = "GEOS57_READ_CN (geos57_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( GEOS_57_DIR )
    CALL Expand_Date( dir, 20110101, 000000 )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String() 
    nc_file = 'GEOS572.YYYYMMDD.CN.' // TRIM( nc_file ) 
    CALL Expand_Date( nc_file, 20110101, 000000 )

    ! Construct complete file path
    nc_file = TRIM( DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
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
    CALL Transfer_2d( Q, FRLAKE )

    ! Read FRLAND
    v_name = "FRLAND"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, FRLAND )

    ! Read FRLANDIC
    v_name = "FRLANDIC"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, FRLANDIC )
    
    ! Read FROCEAN
    v_name = "FROCEAN"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, FROCEAN )
    
    ! Read PHIS
    v_name = "PHIS"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, PHIS )

    ! Echo info
    stamp = TimeStamp_String( 20110101, 000000 )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 5  GEOS-5.7-x CN     met fields for ', a )

    !======================================================================
    ! Cleanup and quit
    !======================================================================

    ! Convert PHIS from [m2/s2] to [m]
    PHIS = PHIS / g0

    ! ND67 diagnostic 
    IF ( ND67 > 0 ) THEN
       AD67(:,:,15) = AD67(:,:,15) + PHIS  ! Sfc geopotential [m]
    ENDIF

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE Geos57_Read_CN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geos57_read_a1
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-5.7.2
!  met fields file containing 1-hr time-averaged (A1) data.  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Geos57_Read_A1( YYYYMMDD, HHMMSS )
!
! !USES:

    USE DAO_MOD, ONLY : ALBEDO   => ALBD
    USE DAO_MOD, ONLY : CLDTOT   => CLDFRC
    USE DAO_MOD, ONLY : EFLUX
    USE DAO_MOD, ONLY : EVAP
    USE DAO_MOD, ONLY : FRSEAICE
    USE DAO_MOD, ONLY : FRSNO
    USE DAO_MOD, ONLY : GRN
    USE DAO_MOD, ONLY : GWETROOT
    USE DAO_MOD, ONLY : GWETTOP
    USE DAO_MOD, ONLY : HFLUX
    USE DAO_MOD, ONLY : LAI
    USE DAO_MOD, ONLY : LWI
    USE DAO_MOD, ONLY : LWGNT    => RADLWG
    USE DAO_MOD, ONLY : PARDF
    USE DAO_MOD, ONLY : PARDR
    USE DAO_MOD, ONLY : PBLH     => PBL
    USE DAO_MOD, ONLY : PRECANV  => PREANV
    USE DAO_MOD, ONLY : PRECCON  => PRECON
    USE DAO_MOD, ONLY : PRECLSC  => PRELSC
    USE DAO_MOD, ONLY : PRECSNO
    USE DAO_MOD, ONLY : PRECTOT  => PREACC
    USE DAO_MOD, ONLY : SEAICE00  
    USE DAO_MOD, ONLY : SEAICE10
    USE DAO_MOD, ONLY : SEAICE20
    USE DAO_MOD, ONLY : SEAICE30
    USE DAO_MOD, ONLY : SEAICE40
    USE DAO_MOD, ONLY : SEAICE50
    USE DAO_MOD, ONLY : SEAICE60
    USE DAO_MOD, ONLY : SEAICE70
    USE DAO_MOD, ONLY : SEAICE80
    USE DAO_MOD, ONLY : SEAICE90
    USE DAO_MOD, ONLY : SLP
    USE DAO_MOD, ONLY : SNODP
    USE DAO_MOD, ONLY : SNOMAS
    USE DAO_MOD, ONLY : SWGDN    => RADSWG
    USE DAO_MOD, ONLY : SWGNT    => RADSWG
    USE DAO_MOD, ONLY : TROPPT   => TROPP
    USE DAO_MOD, ONLY : T2M      => TS
    USE DAO_MOD, ONLY : TS       => TSKIN
    USE DAO_MOD, ONLY : U10M
    USE DAO_MOD, ONLY : USTAR
    USE DAO_MOD, ONLY : V10M
    USE DAO_MOD, ONLY : Z0M      => Z0
!
! !INPUT PARAMETERS:
! 
    INTEGER, INTENT(IN) :: YYYYMMDD    ! GMT date in YYYY/MM/DD format
    INTEGER, INTENT(IN) :: HHMMSS      ! GMT time in hh:mm:ss   format
!
! !REMARKS:
!  This routine was automatically generated by the Perl script ncCodeRead, 
!  and was subsequently hand-edited for compatibility with GEOS-Chem.
!                                                                             .
!  Even though the netCDF file is self-describing, the GEOS-5.7.2 data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-5.7.2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!                                                                             .
!  Special handling for surface precipitation fields:
!  ---------------------------------------------------------------------------
!  In GEOS-5.7.x (and in MERRA), the PRECTOT etc. surface precipitation
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                ! netCDF file ID
    INTEGER            :: X, Y, T            ! netCDF file dimensions
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
 20    FORMAT( '     - GEOS-5.7.x A1 met fields for ', a,  &
               ' have been read already'                  ) 
       RETURN
    ENDIF

    ! Name of this routine (for error printout)
    caller  = "GEOS57_READ_A1 (geos57_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( GEOS_57_DIR )
    CALL Expand_Date( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOS572.YYYYMMDD.A1.' // TRIM( nc_file )
    CALL Expand_Date( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
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
    CALL Transfer_2d( Q, ALBEDO )

    ! Read CLDTOT
    v_name = "CLDTOT"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, CLDTOT )

    ! Read EFLUX
    v_name = "EFLUX"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, EFLUX )

    ! Read EVAP
    v_name = "EVAP"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, EVAP )

    ! Read FRSEAICE
    v_name = "FRSEAICE"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, FRSEAICE )

    ! Read FRSNO
    v_name = "FRSNO"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, FRSNO )

    ! Read GRN
    v_name = "GRN"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, GRN )

    ! Read GWETROOT
    v_name = "GWETROOT"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, GWETROOT )

    ! Read GWETTOP
    v_name = "GWETTOP"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, GWETTOP )

    ! Read HFLUX from file
    v_name = "HFLUX"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, HFLUX )

    ! Read LAI
    v_name = "LAI"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, LAI )

    ! Read LWI
    v_name = "LWI"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, LWI )

    ! Read LWGNT 
    v_name = "LWGNT"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, LWGNT )

    !-----------------------------------------------------------------------
    ! Comment this out for now, this field isn't needed (bmy, 2/2/12)
    !! Read LWTUP
    !v_name = "LWTUP"
    !CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    !CALL Transfer_2d( Q, FRLAKE )
    !-----------------------------------------------------------------------

    ! Read PARDF
    v_name = "PARDF"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, PARDF )

    ! Read PARDR
    v_name = "PARDR"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, PARDR )

    ! Read PBLH
    v_name = "PBLH"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, PBLH )

    ! Read PRECANV
    v_name = "PRECANV"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, PRECANV )

    ! Read PRECCON
    v_name = "PRECCON"
    CALL NcRd( PRECCON, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, PRECCON )

    ! Read PRECLSC
    v_name = "PRECLSC"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, PRECLSC )

    ! Read PRECSNO
    v_name = "PRECSNO"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, PRECSNO )

    ! Read PRECTOT
    v_name = "PRECTOT"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, PRECTOT )

    !-----------------------------------------------------------------------
    ! Comment this out for now, this field isn't needed (bmy, 2/2/12)
    !! Read QV2M
    !v_name = "QV2M"
    !CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    !CALL Transfer_2d( Q, QV2M )
    !-----------------------------------------------------------------------

    ! Read SEAICE00
    v_name = "SEAICE00"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, SEAICE00 )

    ! Read SEAICE10
    v_name = "SEAICE10"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, SEAICE10 )

    ! Read SEAICE20
    v_name = "SEAICE20"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, SEAICE20 )

    ! Read SEAICE30
    v_name = "SEAICE30"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, SEAICE30 )

    ! Read SEAICE40
    v_name = "SEAICE40"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, SEAICE40 )

    ! Read SEAICE50
    v_name = "SEAICE50"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, SEAICE50 )

    ! Read SEAICE60 
    v_name = "SEAICE60"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, SEAICE60 )

    ! Read SEAICE70
    v_name = "SEAICE70"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, SEAICE70 )

    ! Read SEAICE80
    v_name = "SEAICE80"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, SEAICE80 )

    ! Read SEAICE90
    v_name = "SEAICE90"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, SEAICE90 )

    ! Read SLP
    v_name = "SLP"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, SLP )

    ! Read SNODP
    v_name = "SNODP"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, SNODP )

    ! Read SNOMAS
    v_name = "SNOMAS"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, SNOMAS )

    ! Read SWGDN
    v_name = "SWGDN"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, SWGDN )

    ! Read TROPPT
    v_name = "TROPPT"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, TROPPT )

    ! Read TS
    v_name = "TS"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, TS )

    ! Read T2M
    v_name = "T2M"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, T2M )

    ! Read U10M
    v_name = "U10M"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, U10M )

    ! Read USTAR
    v_name = "USTAR"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, USTAR )

    ! Read V10M
    v_name = "V10M"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, V10M )

    ! Read Z0M
    v_name = "Z0M"
    CALL NcRd( Q, fId, TRIM(v_name), st3d, ct3d )
    CALL Transfer_2d( Q, Z0M )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp                    
 10 FORMAT( '     - Found all 44 GEOS-5.7-x A1     met fields for ', a )

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! Close netCDF file
    CALL NcCl( fId )

    ! Increment the # of times A1 fields are read from disk
    CALL Set_Ct_A1( INCREMENT=.TRUE. )
      
    ! Convert surface precip fields from [kg/m2/s] --> [mm/day]
    PRECANV = PRECANV * 86400d0
    PRECCON = PRECCON * 86400d0
    PRECLSC = PRECLSC * 86400d0
    PRECTOT = PRECTOT * 86400d0

    ! ND67 diagnostic: surface fields
    IF ( ND67 > 0 ) THEN
       AD67(:,:,1 ) = AD67(:,:,1 ) + HFLUX    ! Sensible heat flux [W/m2]
       AD67(:,:,2 ) = AD67(:,:,2 ) + SWGDN    ! Incident SW rad @ sfc [W/m2]
       AD67(:,:,3 ) = AD67(:,:,3 ) + PRECTOT  ! Tot prec @ sfc [kg/m2/s]
       AD67(:,:,4 ) = AD67(:,:,4 ) + PRECCON  ! CV prec @ sfc [kg/m2/s]
       AD67(:,:,5 ) = AD67(:,:,5 ) + T2M      ! T @ 2m height [K]
       AD67(:,:,6 ) = AD67(:,:,6 ) + 0e0      !
       AD67(:,:,7 ) = AD67(:,:,7 ) + USTAR    ! Friction velocity [m/s]
       AD67(:,:,8 ) = AD67(:,:,8 ) + Z0M      ! Roughness height [m]
       AD67(:,:,9 ) = AD67(:,:,9 ) + PBLH     ! PBL height [m]
       AD67(:,:,10) = AD67(:,:,10) + CLDTOT   ! Column cld fraction
       AD67(:,:,11) = AD67(:,:,11) + U10M     ! U-wind @ 10m [m/s]
       AD67(:,:,12) = AD67(:,:,12) + V10M     ! V-wind @ 10m [m/s]
       AD67(:,:,14) = AD67(:,:,14) + ALBEDO   ! Sfc albedo [unitless]
       AD67(:,:,17) = AD67(:,:,17) + TROPPT   ! T'pause pressure [hPa]
       AD67(:,:,18) = AD67(:,:,18) + SLP      ! Sea level pressure [hPa]
       AD67(:,:,19) = AD67(:,:,19) + TS       ! Sfc skin temperature [K]
       AD67(:,:,20) = AD67(:,:,20) + PARDF    ! Diffuse PAR [W/m2]
       AD67(:,:,21) = AD67(:,:,21) + PARDR    ! Direct PAR [W/m2]
       AD67(:,:,22) = AD67(:,:,22) + GWETTOP  ! Topsoil wetness [frac]
       AD67(:,:,23) = AD67(:,:,23) + EFLUX    ! Latent heat flux [W/m2]
    ENDIF

    ! Save date & time for next iteration
    lastDate = YYYYMMDD
    lastTime = HHMMSS

  END SUBROUTINE Geos57_Read_A1
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geos57_read_a3
!
! !DESCRIPTION: Convenience wrapper for the following routines which read
!  3-hour time averaged data from disk:
! \begin{itemize}
! \item Geos57\_Read\_A3cld
! \item Geos57\_Read\_A3dyn
! \item Geos57\_Read\_A3mstC
! \item Geos57\_Read\_A3mstE
! \end{itemize}
!
! !INTERFACE:
!
  SUBROUTINE Geos57_Read_A3( YYYYMMDD, HHMMSS )
!
! !INPUT PARAMETERS:
! 
    INTEGER, INTENT(IN) :: YYYYMMDD       ! GMT date in YYYY/MM/DD format
    INTEGER, INTENT(IN) :: HHMMSS         ! GMT time in hh:mm:ss   format
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
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
 20    FORMAT( '     - GEOS-5.7.x A3 met fields for ', a,  &
               ' have been read already'                  ) 
       RETURN
    ENDIF

    ! Save date & time for next iteration
    lastDate = YYYYMMDD
    lastTime = HHMMSS

    ! Read all the diffeent A3 files
    CALL Geos57_Read_A3cld ( YYYYMMDD, HHMMSS )
    CALL Geos57_Read_A3dyn ( YYYYMMDD, HHMMSS )
    CALL Geos57_Read_A3mstC( YYYYMMDD, HHMMSS )
    CALL Geos57_Read_A3mstE( YYYYMMDD, HHMMSS )

    !======================================================================
    ! Cleanup and quit
    !======================================================================

    ! Increment the # of times that A3 fields have been read
    CALL Set_Ct_A3( INCREMENT=.TRUE. )

    ! Save date & time for next iteration
    lastDate = YYYYMMDD
    lastTime = HHMMSS

  END SUBROUTINE Geos57_Read_A3
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geos57_read_a3cld
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-5.7.2
!  met fields file containing 3-hr time-averaged (A3) data (cloud fields).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Geos57_Read_A3cld( YYYYMMDD, HHMMSS )
!
! !USES:
!
    USE DAO_MOD, ONLY : CLOUD    => CLDF
    USE DAO_MOD, ONLY : OPTDEPTH => OPTDEP
    USE DAO_MOD, ONLY : QI
    USE DAO_MOD, ONLY : QL
    USE DAO_MOD, ONLY : TAUCLI
    USE DAO_MOD, ONLY : TAUCLW
!
! !INPUT PARAMETERS:
! 
    INTEGER, INTENT(IN) :: YYYYMMDD    ! GMT date in YYYY/MM/DD format
    INTEGER, INTENT(IN) :: HHMMSS      ! GMT time in hh:mm:ss   format
!
! !REMARKS:
!  This routine was automatically generated by the Perl script ncCodeRead, 
!  and was subsequently hand-edited for compatibility with GEOS-Chem.
!                                                                             .
!  Even though the netCDF file is self-describing, the GEOS-5.7.2 data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-5.7.2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
!  05 Apr 2012 - R. Yantosca - Fixed bug: TAUCLI was overwritten w/ TAUCLW
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
    caller  = "GEOS57_READ_A3cld (geos57_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( GEOS_57_DIR )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOS572.YYYYMMDD.A3cld.' // TRIM( nc_file )
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
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
    CALL TRANSFER_A6( Q, CLOUD )
    
    ! Read OPTDEPTH
    v_name = "OPTDEPTH"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL TRANSFER_A6( Q, OPTDEPTH )

    ! Read QI
    v_name = "QI"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, QI )

    ! Read QL
    v_name = "QL"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, QL )

    ! Read TAUCLI
    v_name = "TAUCLI"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, TAUCLI )

    ! Read TAUCLW
    v_name = "TAUCLW"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, TAUCLW )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 6  GEOS-5.7-x A3cld  met fields for ', a )

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE GEOS57_READ_A3cld
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geos57_read_a3dyn
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-5.7.2
!  met fields file containing 3-hr time-averaged (A3) data (dynamics fields).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS57_READ_A3dyn( YYYYMMDD, HHMMSS )
!
! !USES:
!
    USE DAO_MOD, ONLY : CLDTOPS 
    USE DAO_MOD, ONLY : CMFMC
    USE DAO_MOD, ONLY : DTRAIN
   !USE DAO_MOD, ONLY : OMEGA
    USE DAO_MOD, ONLY : RH
    USE DAO_MOD, ONLY : U      => UWND
    USE DAO_MOD, ONLY : V      => VWND
!
! !INPUT PARAMETERS:
! 
    INTEGER, INTENT(IN) :: YYYYMMDD    ! GMT date in YYYY/MM/DD format
    INTEGER, INTENT(IN) :: HHMMSS      ! GMT time in hh:mm:ss   format
!
! !REMARKS:
!  This routine was automatically generated by the Perl script ncCodeRead, 
!  and was subsequently hand-edited for compatibility with GEOS-Chem.
!                                                                             .
!  Even though the netCDF file is self-describing, the GEOS-5.7.2 data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-5.7.2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
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
    caller  = "GEOS57_READ_A3dyn (geos57_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( GEOS_57_DIR )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOS572.YYYYMMDD.A3dyn.' // TRIM( nc_file )
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
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
    ! Read data on level edges
    !--------------------------------

    ! netCDF start & count indices
    st4d      = (/ 1,     1,     1,       time_index /)      
    ct4d      = (/ IIPAR, JJPAR, LGLOB+1, 1          /)

    ! Read CMFMC
    v_name = "CMFMC"
    CALL NcRd( Qe, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe, CMFMC )

    !--------------------------------
    ! Read data on level centers
    !--------------------------------

    ! netCDF start & count indices
    st4d      = (/ 1,     1,     1,       time_index /)      
    ct4d      = (/ IIPAR, JJPAR, LGLOB,   1          /)

    ! Read DTRAIN
    v_name = "DTRAIN"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, DTRAIN )

    !! Read OMEGA  from file
    !v_name = "OMEGA"
    !CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    !CALL Transfer_3d( Q, OMEGA )

    ! Read RH
    v_name = "RH"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, RH )

    ! Read U
    v_name = "U"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, U )

    ! Read V
    v_name = "V"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, V )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 6  GEOS-5.7-x A3dyn  met fields for ', a )

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================
    
    ! CLDTOPS = highest location of CMFMC in the column (I,J)
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       CLDTOPS(I,J) = 1
       DO L = LLPAR, 1, -1
          IF ( CMFMC(I,J,L) > 0d0 ) THEN
             CLDTOPS(I,J) = L + 1
             EXIT
          ENDIF
       ENDDO
    ENDDO
    ENDDO

    ! ND66 diagnostic: U, V, CMFMC, DTRAIN met fields
    IF ( ND66 > 0 ) THEN
       AD66(:,:,1:LD66,1) = AD66(:,:,1:LD66,1) + U     (:,:,1:LD66) ! [m/s    ]
       AD66(:,:,1:LD66,2) = AD66(:,:,1:LD66,2) + V     (:,:,1:LD66) ! [m/s    ]
       AD66(:,:,1:LD66,5) = AD66(:,:,1:LD66,5) + CMFMC (:,:,1:LD66) ! [kg/m2/s]
       AD66(:,:,1:LD66,6) = AD66(:,:,1:LD66,6) + DTRAIN(:,:,1:LD66) ! [kg/m2/s]
    ENDIF

    ! ND67 diagnostic: CLDTOPS
    IF ( ND67 > 0 ) THEN
       AD67(:,:,16) = AD67(:,:,16) + CLDTOPS                        ! [levels]
    ENDIF

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE GEOS57_READ_A3dyn
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geos57_read_a3mstc
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-5.7.2
!  met fields file containing 3-hr time-averaged (A3) data (moist fields,
!  saved on level centers).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS57_READ_A3mstC( YYYYMMDD, HHMMSS )
!
! !USES:
!
    USE DAO_MOD, ONLY : DQRCU
    USE DAO_MOD, ONLY : DQRLSAN
    USE DAO_MOD, ONLY : REEVAPCN
    USE DAO_MOD, ONLY : REEVAPLS
!
! !INPUT PARAMETERS:
! 
    INTEGER, INTENT(IN) :: YYYYMMDD    ! GMT date in YYYY/MM/DD format
    INTEGER, INTENT(IN) :: HHMMSS      ! GMT time in hh:mm:ss   format
!
! !REMARKS:
!  This routine was automatically generated by the Perl script ncCodeRead, 
!  and was subsequently hand-edited for compatibility with GEOS-Chem.
!                                                                             .
!  Even though the netCDF file is self-describing, the GEOS-5.7.2 data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-5.7.2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
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
    caller  = "GEOS57_READ_A3mstC (geos57_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( GEOS_57_DIR )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOS572.YYYYMMDD.A3mstC.' // TRIM( nc_file )
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
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
    CALL Transfer_3d( Q, DQRCU )

    ! Read DQRLSAN
    v_name = "DQRLSAN"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, DQRLSAN )

    ! Read REEVAPCN
    v_name = "REEVAPCN"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, REEVAPCN )

    ! Read  from file
    v_name = "REEVAPLS"
    CALL NcRd( Q, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q, REEVAPLS )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 4  GEOS-5.7-x A3mstC met fields for ', a )

    !======================================================================
    ! Cleanup and quit
    !======================================================================

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE GEOS57_READ_A3mstC
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geos57_read_a3mste
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-5.7.2
!  met fields file containing 3-hr time-averaged (A3) data (moist fields,
!  saved on level edges).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOS57_READ_A3mstE( YYYYMMDD, HHMMSS )
!
! !USES:
!
    USE DAO_MOD, ONLY : PFICU    
    USE DAO_MOD, ONLY : PFILSAN
    USE DAO_MOD, ONLY : PFLCU
    USE DAO_MOD, ONLY : PFLLSAN
!
! !INPUT PARAMETERS:
! 
    INTEGER, INTENT(IN) :: YYYYMMDD    ! GMT date in YYYY/MM/DD format
    INTEGER, INTENT(IN) :: HHMMSS      ! GMT time in hh:mm:ss   format
!
! !REMARKS:
!  This routine was automatically generated by the Perl script ncCodeRead, 
!  and was subsequently hand-edited for compatibility with GEOS-Chem.
!                                                                             .
!  Even though the netCDF file is self-describing, the GEOS-5.7.2 data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-5.7.2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
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
    caller  = "GEOS57_READ_A3mstE (geos57_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( GEOS_57_DIR )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOS572.YYYYMMDD.A3mstE.' // TRIM( nc_file )
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
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

    ! Read PFICU
    v_name = "PFICU"
    CALL NcRd( Qe, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe, PFICU )

    ! Read PFILSAN
    v_name = "PFILSAN"
    CALL NcRd( Qe, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe, PFILSAN )

    ! Read PFLCU
    v_name = "PFLCU"
    CALL NcRd( Qe, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe, PFLCU )

    ! Read  from file
    v_name = "PFLLSAN"
    CALL NcRd( Qe, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d_Lp1( Qe,  PFLLSAN )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 4  GEOS-5.7-x A3mstE met fields for ', a )

    !=================================================================
    ! Cleanup and quit
    !=================================================================

    ! Close netCDF file
    CALL NcCl( fId )

  END SUBROUTINE Geos57_Read_A3mstE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geos57_read_I3_1
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-5.7.2
!  met fields file containing 3-hr instantaneous (I3) data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Geos57_Read_I3_1( YYYYMMDD, HHMMSS )
!
! !USES:
!
    USE DAO_MOD, ONLY : PS1
   !USE DAO_MOD, ONLY : PV1
    USE DAO_MOD, ONLY : QV1 => SPHU1
    USE DAO_MOD, ONLY : T1  => TMPU1
!
! !INPUT PARAMETERS:
! 
    INTEGER, INTENT(IN) :: YYYYMMDD    ! GMT date in YYYY/MM/DD format
    INTEGER, INTENT(IN) :: HHMMSS      ! GMT time in hh:mm:ss   format
!
! !REMARKS:
!  This routine was automatically generated by the Perl script ncCodeRead, 
!  and was subsequently hand-edited for compatibility with GEOS-Chem.
!                                                                             .
!  Even though the netCDF file is self-describing, the GEOS-5.7.2 data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-5.7.2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
!  05 Apr 2012 - R. Yantosca - Now convert QV1 from [kg/kg] to [g/kg]
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
    caller  = "GEOS57_READ_I3_1 (geos57_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( GEOS_57_DIR )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOS572.YYYYMMDD.I3.' // TRIM( nc_file )
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
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
    CALL Transfer_2d( Q2, PS1 )

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
    !CALL Transfer_3d( Q3, PV )
    !----------------------------------------------------------------

    ! Read QV
    v_name = "QV"
    CALL NcRd( Q3, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q3, QV1 )

    ! Read T
    v_name = "T"
    CALL NcRd( Q3, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q3, T1 )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 4  GEOS-5.7-x I3     met fields for ', a )

    !-------------------------------------------------
    ! Unit conversions & special handling
    !-------------------------------------------------
    WHERE ( QV1 < 0d0 ) 

       ! NOTE: Now set negative Q to a small positive # 
       ! instead of zero, so as not to blow up logarithms
       QV1 = 1d-32

    ELSEWHERE

       ! Convert GEOS-5.7.x specific humidity from [kg/kg] to [g/kg]
       QV1 = QV1 * 1000d0

    ENDWHERE

    !======================================================================
    ! Diagnostics, cleanup, and quit
    !======================================================================

    ! Close netCDF file
    CALL NcCl( fId )

    ! Increment the # of times I3 fields have been read
    CALL Set_Ct_I3( INCREMENT=.TRUE. )

    ! ND66 diagnostic: T1, QV1 met fields
    IF ( ND66 > 0 ) THEN
       AD66(:,:,1:LD66,3) = AD66(:,:,1:LD66,3) + T1 (:,:,1:LD66) ! [K   ]
       AD66(:,:,1:LD66,4) = AD66(:,:,1:LD66,4) + QV1(:,:,1:LD66) ! [g/kg]
    ENDIF

  END SUBROUTINE Geos57_Read_I3_1
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geos57_read_I3_2
!
! !DESCRIPTION: Routine to read variables and attributes from a GEOS-5.7.2
!  met fields file containing 3-hr instantaneous (I3) data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Geos57_Read_I3_2( YYYYMMDD, HHMMSS )
!
! !USES:
!
    USE DAO_MOD, ONLY : PS2
   !USE DAO_MOD, ONLY : PV2
    USE DAO_MOD, ONLY : QV2 => SPHU2
    USE DAO_MOD, ONLY : T2  => TMPU2
!
! !INPUT PARAMETERS:
! 
    INTEGER, INTENT(IN) :: YYYYMMDD    ! GMT date in YYYY/MM/DD format
    INTEGER, INTENT(IN) :: HHMMSS      ! GMT time in hh:mm:ss   format
!
! !REMARKS:
!  This routine was automatically generated by the Perl script ncCodeRead, 
!  and was subsequently hand-edited for compatibility with GEOS-Chem.
!                                                                             .
!  Even though the netCDF file is self-describing, the GEOS-5.7.2 data, 
!  dimensions, and units are pre-specified according to the GMAO GEOS-5.7.2
!  file specification.  Therefore we can "cheat" a little bit and not have
!  to read netCDF attributes to find what these values are.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  07 Feb 2012 - R. Yantosca - Now echo info after reading fields from disk
!  10 Feb 2012 - R. Yantosca - Now get a string for the model resolution
!  05 Apr 2012 - R. Yantosca - Now convert QV2 from [kg/kg] to [g/kg]
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
    caller  = "GEOS57_READ_I3_2 (geos57_read_mod.F90)"

    ! Replace time & date tokens in the file name
    dir     = TRIM( GEOS_57_DIR )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = Get_Resolution_String()
    nc_file = 'GEOS572.YYYYMMDD.I3.' // TRIM( nc_file )
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
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
    CALL Transfer_2d( Q2, PS2 )

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
    !CALL Transfer_3d( Q3, PV )
    !----------------------------------------------------------------

    ! Read QV
    v_name = "QV"
    CALL NcRd( Q3, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q3, QV2 )

    ! Read T
    v_name = "T"
    CALL NcRd( Q3, fId, TRIM(v_name), st4d, ct4d )
    CALL Transfer_3d( Q3, T2 )

    ! Echo info
    stamp = TimeStamp_String( YYYYMMDD, HHMMSS )
    WRITE( 6, 10 ) stamp
 10 FORMAT( '     - Found all 4  GEOS-5.7-x I3     met fields for ', a )

    !-------------------------------------------------
    ! Unit conversions & special handling
    !-------------------------------------------------
    WHERE ( QV2 < 0d0 ) 

       ! NOTE: Now set negative Q to a small positive # 
       ! instead of zero, so as not to blow up logarithms
       QV2 = 1d-32

    ELSEWHERE

       ! Convert GEOS-5.7.x specific humidity from [kg/kg] to [g/kg]
       QV2 = QV2 * 1000d0

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
       AD66(:,:,1:LD66,3) = AD66(:,:,1:LD66,3) + T2 (:,:,1:LD66) ! [K   ]
       AD66(:,:,1:LD66,4) = AD66(:,:,1:LD66,4) + QV2(:,:,1:LD66) ! [g/kg]
    ENDIF

  END SUBROUTINE Geos57_Read_I3_2
!EOC
END MODULE Geos57_Read_Mod
