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
MODULE GEOS57_READ_MOD
!
! !USES:
!
#if defined( USE_NETCDF )
  ! NcdfUtil modules for netCDF I/O
  USE m_netcdf_io_open                    ! netCDF open
  USE m_netcdf_io_get_dimlen              ! netCDF dimension queries
  USE m_netcdf_io_read                    ! netCDF data reads
  USE m_netcdf_io_readattr                ! netCDF attribute reads
  USE m_netcdf_io_close                   ! netCDF close
#endif

  ! GEOS-Chem modules
  USE CMN_SIZE_MOD                        ! Size parameters
  USE CMN_GCTM_MOD                        ! Physical constants
  USE CMN_DIAG_MOD                        ! Diagnostic arrays & counters
  USE DIAG_MOD,      ONLY : AD66          ! Array for ND66 diagnostic  
  USE DIAG_MOD,      ONLY : AD67          ! Array for ND67 diagnostic
  USE DIRECTORY_MOD                       ! Directory paths
  USE ERROR_MOD,     ONLY : ERROR_STOP    ! Stop w/ error message
  USE TIME_MOD,      ONLY : EXPAND_DATE   ! Routine for YMD token replace
  USE TRANSFER_MOD                        ! Routines for casting 

  IMPLICIT NONE
  PRIVATE

#if defined( USE_NETCDF )
# include "netcdf.inc"                    ! Include file for netCDF library
#endif
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: CHECK_DIMENSIONS
!
! !PUBLIC MEMBER FUNCTIONS:
! 
  PUBLIC  :: GEOS57_READ_CN
  PUBLIC  :: GEOS57_READ_A1
!  PUBLIC  :: GEOS57_READ_A3CLD
!  PUBLIC  :: GEOS57_READ_A3DYN
!  PUBLIC  :: GEOS57_READ_A3MSTC
!  PUBLIC  :: GEOS57_READ_A3MSTE
!  PUBLIC  :: GEOS57_READ_I3
!

! !REMARKS:
!  Assumes that you have a netCDF library (either v3 or v4) installed on your 
!  system. 
!
! !REVISION HISTORY:
!  19 Aug 2010 - R. Yantosca - Initial version, based on i6_read_mod.f
!  20 Aug 2010 - R. Yantosca - Moved include files to top of module
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
! !IROUTINE: check_dimensions
!
! !DESCRIPTION: Subroutine CHECK\_DIMENSIONS checks to see if dimensions read 
!  from the netCDF file match the defined GEOS-Chem dimensions.  If not, then 
!  it will stop the GEOS-Chem simulation with an error message.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHECK_DIMENSIONS( lon,          lat,  lev,           &
                               lev_expected, time, time_expected )
!
! !INPUT PARAMETERS:
!
    INTEGER, OPTIONAL :: lon             ! Longitude dimension
    INTEGER, OPTIONAL :: lat             ! Latitude dimension
    INTEGER, OPTIONAL :: lev             ! Altitude dimension
    INTEGER, OPTIONAL :: lev_expected    ! Expected # of levels
    INTEGER, OPTIONAL :: time            ! Time dimension
    INTEGER, OPTIONAL :: time_expected   ! Expected # of times
!
! !REMARKS:
!  Call this routine with keyword arguments, i.e.
!     CALL CHECK_DIMENSION( lon=X, lat=Y, lev=Z, time=T )
!
! !REVISION HISTORY:
!  02 Feb 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: errMsg, location

    ! routine 
    location = "CHECK_DIMENSION (geos57_read_mod.F90)"

    ! Error check longitude dimension 
    IF ( PRESENT( lon ) ) THEN
       IF ( lon /= IIPAR ) THEN
          WRITE( errMsg, 100 ) lon, IIPAR
 100      FORMAT( 'Longitude dimension (', i5, &
                  ' ) does not match IIPAR ( ', i5, ')!' )
          CALL ERROR_STOP( errMsg, location )
       ENDIF
    ENDIF


    ! Error check longitude dimension 
    IF ( PRESENT( lat ) ) THEN
       IF ( lat /= JJPAR ) THEN
          WRITE( errMsg, 110 ) lat, JJPAR
 110      FORMAT( 'Latitude dimension (', i5, &
                  ' ) does not match JJPAR ( ', i5, ')!' )
          CALL ERROR_STOP( errMsg, location )
       ENDIF
    ENDIF


    ! Error check longitude dimension 
    IF ( PRESENT( lev ) ) THEN
       IF ( lev /= lev_expected ) THEN
          WRITE( errMsg, 120 ) lev, lev_expected
 120      FORMAT( 'Levels dimension (', i5, &
                  ' ) does not match expected # of levels ( ', i5, ')!' )
          CALL ERROR_STOP( errMsg, location )
       ENDIF
    ENDIF

    ! Error check time dimension 
    IF ( PRESENT( time ) .and. PRESENT( time_expected ) ) THEN
       IF ( time /= time_expected ) THEN
          WRITE( errMsg, 130 ) time, time_expected
 130      FORMAT( 'Time dimension (', i5, &
                  ' ) does not match expected # of times ( ', i5, ')!' )
          CALL ERROR_STOP( errMsg, location )
       ENDIF
    ENDIF

  END SUBROUTINE CHECK_DIMENSIONS
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
  SUBROUTINE GEOS57_READ_CN()
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId                ! netCDF file ID
    INTEGER            :: X, Y, T            ! netCDF file dimensions
    CHARACTER(LEN=255) :: nc_file            ! netCDF file name
    CHARACTER(LEN=255) :: v_name             ! netCDF variable name
    CHARACTER(LEN=255) :: dir                ! Data directory path
    CHARACTER(LEN=255) :: errMsg             ! Error message
                                
    ! Arrays                                 
    INTEGER            :: st3d(3), ct3d(3)   ! Start + count, for 3D arrays 
    REAL*4             :: TEMP(IIPAR,JJPAR)  ! Temporary data arrray

#if defined( USE_NETCDF )

    !======================================================================
    ! Open the netCDF file
    !======================================================================

    ! Replace time & date tokens in the file name
    dir = TRIM( GEOS_57_DIR )
    CALL EXPAND_DATE( dir, 20110101, 000000 )

    ! Replace time & date tokens in the file name
    nc_file = 'GEOS572.YYYYMMDD.CN.4x5.nc'
    CALL EXPAND_DATE( nc_file, 20110101, 000000 )

    ! Construct complete file path
    nc_file = TRIM( DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'time',  T )

    ! Make sure the dimensions of the file are valid
    CALL Check_Dimensions( lon=X, lat=Y, time=T, time_expected=1 )

    !======================================================================
    ! Read data from netCDF file
    !======================================================================

    ! Start and count indices
    st3d   = (/ 1,     1,     1 /)
    ct3d   = (/ IIPAR, JJPAR, 1 /)

    ! Read FRLAKE
    v_name = "FRLAKE"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, FRLAKE )

    ! Read FRLAND
    v_name = "FRLAND"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, FRLAND )

    ! Read FRLANDIC
    v_name = "FRLANDIC"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, FRLANDIC )
    
    ! Read FROCEAN
    v_name = "FROCEAN"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, FROCEAN )
    
    ! Read PHIS
    v_name = "PHIS"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, PHIS )

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

#endif

  END SUBROUTINE GEOS57_READ_CN
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
  SUBROUTINE GEOS57_READ_A1( YYYYMMDD, HHMMSS )
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
    INTEGER            :: fId                ! netCDF file ID
    INTEGER            :: X, Y, T            ! netCDF file dimensions
    INTEGER            :: time_index         ! Read this time slice of data
    CHARACTER(LEN=255) :: nc_file            ! netCDF file name
    CHARACTER(LEN=255) :: v_name             ! netCDF variable name 
    CHARACTER(LEN=255) :: dir                ! Data directory path
    CHARACTER(LEN=255) :: errMsg             ! Error message
                                             
    ! Arrays                                 
    INTEGER            :: st3d(3), ct3d(3)   ! Start + count, for 3D arrays 
    REAL*4             :: TEMP(IIPAR,JJPAR)  ! Temporary data arrray

#if defined( USE_NETCDF )

    !======================================================================
    ! Open the netCDF file
    !======================================================================

    ! Replace time & date tokens in the file name
    dir = TRIM( GEOS_57_DIR )
    CALL EXPAND_DATE( dir, YYYYMMDD, HHMMSS )

    ! Replace time & date tokens in the file name
    nc_file = 'GEOS572.YYYYMMDD.A1.4x5.nc'
    CALL EXPAND_DATE( nc_file, YYYYMMDD, HHMMSS )

    ! Construct complete file path
    nc_file = TRIM( DATA_DIR ) // TRIM( dir ) // TRIM( nc_file )
    
    ! Open netCDF file
    CALL NcOp_Rd( fId, TRIM( nc_file ) )

    ! Read the dimensions from the netCDF file
    CALL NcGet_DimLen( fId, 'lon',   X )
    CALL NcGet_DimLen( fId, 'lat',   Y )
    CALL NcGet_DimLen( fId, 'time',  T )

    ! Make sure the dimensions of the file are valid
    CALL Check_Dimensions( lon=X, lat=Y, time=T, time_expected=24 )
    
    !======================================================================
    ! Read data from disk
    !======================================================================
    
    ! Find the proper time-slice to read from disk
    time_index = ( HHMMSS / 10000 ) + 1

    ! Stop w/ error if the time index is invalid
    IF ( time_index < 1 .or. time_index > T ) THEN
       WRITE( 6, 100 ) time_index
 100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 24!' )
       CALL ERROR_STOP( errMsg, 'GEOS57_READ_A1 (geos57_read_mod.F90)' )
    ENDIF

    ! netCDF start & count indices
    st3d      = (/ 1,     1,     time_index /)      
    ct3d      = (/ IIPAR, JJPAR, 1          /)

    ! Read ALBEDO
    v_name = "ALBEDO"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, ALBEDO )

    ! Read CLDTOT
    v_name = "CLDTOT"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, CLDTOT )

    ! Read EFLUX
    v_name = "EFLUX"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, EFLUX )

    ! Read EVAP
    v_name = "EVAP"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, EVAP )

    ! Read FRSEAICE
    v_name = "FRSEAICE"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, FRSEAICE )

    ! Read FRSNO
    v_name = "FRSNO"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, FRSNO )

    ! Read GRN
    v_name = "GRN"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, GRN )

    ! Read GWETROOT
    v_name = "GWETROOT"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, GWETROOT )

    ! Read GWETTOP
    v_name = "GWETTOP"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, GWETTOP )

    ! Read JFLUX from file
    v_name = "HFLUX"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, HFLUX )

    ! Read LAI
    v_name = "LAI"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, LAI )

    ! Read LWI
    v_name = "LWI"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, LWI )

    ! Read LWGNT 
    v_name = "LWGNT"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, LWGNT )

    !-----------------------------------------------------------------------
    ! Comment this out for now, this field isn't needed (bmy, 2/2/12)
    !! Read LWTUP
    !v_name = "LWTUP"
    !CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    !CALL TRANSFER_2D( TEMP, FRLAKE )
    !-----------------------------------------------------------------------

    ! Read PARDF
    v_name = "PARDF"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, PARDF )

    ! Read PARDR
    v_name = "PARDR"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, PARDR )

    ! Read PBLH
    v_name = "PBLH"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, PBLH )

    ! Read PRECANV
    v_name = "PRECANV"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, PRECANV )

    ! Read PRECCON
    v_name = "PRECCON"
    CALL NcRd( PRECCON, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, PRECCON )

    ! Read PRECLSC
    v_name = "PRECLSC"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, PRECLSC )

    ! Read PRECSNO
    v_name = "PRECSNO"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, PRECSNO )

    ! Read PRECTOT
    v_name = "PRECTOT"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, PRECTOT )

    !-----------------------------------------------------------------------
    ! Comment this out for now, this field isn't needed (bmy, 2/2/12)
    !! Read QV2M
    !v_name = "QV2M"
    !CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    !CALL TRANSFER_2D( TEMP, QV2M )
    !-----------------------------------------------------------------------


    ! Read SEAICE00
    v_name = "SEAICE00"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, SEAICE00 )

    ! Read SEAICE10
    v_name = "SEAICE10"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, SEAICE10 )

    ! Read SEAICE20
    v_name = "SEAICE20"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, SEAICE20 )

    ! Read SEAICE30
    v_name = "SEAICE30"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, SEAICE30 )

    ! Read SEAICE40
    v_name = "SEAICE40"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, SEAICE40 )

    ! Read SEAICE50
    v_name = "SEAICE50"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, SEAICE50 )

    ! Read SEAICE60 
    v_name = "SEAICE60"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, SEAICE60 )

    ! Read SEAICE70
    v_name = "SEAICE70"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, SEAICE70 )

    ! Read SEAICE80
    v_name = "SEAICE80"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, SEAICE80 )

    ! Read SEAICE90
    v_name = "SEAICE90"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, SEAICE90 )

    ! Read SLP
    v_name = "SLP"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, SLP )

    ! Read SNODP
    v_name = "SNODP"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, SNODP )

    ! Read SNOMAS
    v_name = "SNOMAS"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, SNOMAS )

    ! Read SWGDN
    v_name = "SWGDN"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, SWGDN )

    ! Read TROPPT
    v_name = "TROPPT"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, TROPPT )

    ! Read TS
    v_name = "TS"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, TS )

    ! Read T2M
    v_name = "T2M"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, T2M )

    ! Read U10M
    v_name = "U10M"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, U10M )

    ! Read USTAR
    v_name = "USTAR"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, USTAR )

    ! Read V10M
    v_name = "V10M"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, V10M )

    ! Read Z0M
    v_name = "Z0M"
    CALL NcRd( TEMP, fId, TRIM(v_name), st3d, ct3d )
    CALL TRANSFER_2D( TEMP, Z0M )

    !======================================================================
    !        %%%%% SPECIAL HANDLING FOR CERTAIN FIELDS %%%%% 
    !
    ! In GEOS-5.7.x (and in MERRA), the PRECTOT etc. surface precipi 
    ! met fields fields have units of [kg/m2/s].  In all other GEOS 
    ! versions, PREACC and PRECON have units of [mm/day].  
    !
    ! Therefore, for backwards compatibility with existing code, 
    ! apply the following unit conversion to the GEOS-5 PRECTOT and
    ! PRECCON fields:
    !
    !
    !     kg  |    m3    | 86400 s | 1000 mm
    !   ------+----------+---------+--------- = 86400 
    !    m2 s |  1000 kg |  day    |   m
    !              ^
    !              |
    !       1 / density of water 
    !===================================================================
      
    ! Convert from [kg/m2/s] --> [mm/day]
    PRECANV = PRECANV * 86400d0
    PRECCON = PRECCON * 86400d0
    PRECLSC = PRECLSC * 86400d0
    PRECTOT = PRECTOT * 86400d0

    !=================================================================
    ! ND67 diagnostic: A1 surface fields
    !=================================================================
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

    !=================================================================
    ! Cleanup and quit
    !=================================================================

    print*, '### st3d: ', st3d
    print*, '### t2m : ', minval(t2m), maxval(t2m)

    ! Close netCDF file
    CALL NcCl( fId )

#endif

  END SUBROUTINE GEOS57_READ_A1
!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_netcdf_file
!!
!! !DESCRIPTION: Routine to read variables and attributes from a netCDF
!!  file.  This routine was automatically generated by the Perl script
!!  NcdfUtilities/perl/ncCodeRead.
!!\
!!\
!! !INTERFACE:
!!
!      SUBROUTINE READ_FROM_NETCDF_FILE( fId )
!!
!! !USES:
!!
!      ! Modules for netCDF read
!      USE m_netcdf_io_open
!      USE m_netcdf_io_get_dimlen
!      USE m_netcdf_io_read
!      USE m_netcdf_io_readattr
!      USE m_netcdf_io_close
!
!      IMPLICIT NONE
!
!#     include "netcdf.inc"
!!
!! !OUTPUT PARAMETERS:
!!   
!      INTEGER, INTENT(INOUT) :: fId    ! netCDF file ID
!!
!! !REMARKS:
!!  Assumes that you have:
!!  (1) A netCDF library (either v3 or v4) installed on your system
!!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!!                                                                             .
!!  Although this routine was generated automatically, some further
!!  hand-editing may be required (i.e. to  specify the size of parameters, 
!!  and/or to assign values to variables.  Also, you can decide how to handle
!!  the variable attributes (or delete calls for reading attributes that you
!!  do not need).
!!
!! !REVISION HISTORY:
!!  30 Jan 2012 - R. Yantosca - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      !=================================================================
!      ! Variable declarations
!      !=================================================================
!
!      ! Data arrays
!      REAL*4             :: lon(IIPAR)
!      REAL*4             :: lat(JJPAR)
!      REAL*4             :: lev(LLPAR)
!      INTEGER            :: time(1)
!      REAL*4             :: CLOUD(IIPAR,JJPAR,LLPAR,1)
!      REAL*4             :: OPTDEPTH(IIPAR,JJPAR,LLPAR,1)
!      REAL*4             :: QI(IIPAR,JJPAR,LLPAR,1)
!      REAL*4             :: QL(IIPAR,JJPAR,LLPAR,1)
!      REAL*4             :: TAUCLI(IIPAR,JJPAR,LLPAR,1)
!      REAL*4             :: TAUCLW(IIPAR,JJPAR,LLPAR,1)
!
!      ! Character strings
!      CHARACTER(LEN=255) :: nc_file            ! netCDF file name
!      CHARACTER(LEN=255) :: v_name             ! netCDF variable name 
!      CHARACTER(LEN=255) :: a_name             ! netCDF attribute name
!      CHARACTER(LEN=255) :: a_val              ! netCDF attribute value
!
!      ! Arrays for netCDF start and count values
!      INTEGER            :: st1d(1), ct1d(1)   ! For 1D arrays    
!      INTEGER            :: st2d(2), ct2d(2)   ! For 2D arrays 
!      INTEGER            :: st3d(3), ct3d(3)   ! For 3D arrays 
!      INTEGER            :: st4d(4), ct4d(4)   ! For 4D arrays 
!      INTEGER            :: st5d(5), ct5d(5)   ! For 5D arrays 
!      INTEGER            :: st6d(6), ct6d(6)   ! For 6D arrays 
!
!      !=================================================================
!      ! Open and read data from the netCDF file
!      !=================================================================
!
!      ! Open netCDF file
!      nc_file = 'GEOS572.YYYYMMDD.A3cld.4x5.nc'
!      CALL Ncop_Rd( fId, TRIM(nc_file) )
!
!      !----------------------------------------
!      ! VARIABLE: lon
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lon"
!
!      ! Read lon from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ IIPAR /)
!      CALL NcRd( lon, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lon:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lon:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: lat
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lat"
!
!      ! Read lat from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ JJPAR /)
!      CALL NcRd( lat, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lat:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lat:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: lev
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lev"
!
!      ! Read lev from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ LLPAR /)
!      CALL NcRd( lev, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lev:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lev:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: time
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "time"
!
!      ! Read time from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ 1 /)
!      CALL NcRd( time, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the time:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:delta_t attribute
!      a_name = "delta_t"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:begin_date attribute
!      a_name = "begin_date"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:begin_time attribute
!      a_name = "begin_time"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:time_increment attribute
!      a_name = "time_increment"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: CLOUD
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "CLOUD"
!
!      ! Read CLOUD from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( CLOUD, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the CLOUD:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the CLOUD:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the CLOUD:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: OPTDEPTH
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "OPTDEPTH"
!
!      ! Read OPTDEPTH from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( OPTDEPTH, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the OPTDEPTH:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the OPTDEPTH:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the OPTDEPTH:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: QI
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "QI"
!
!      ! Read QI from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( QI, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the QI:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the QI:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the QI:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: QL
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "QL"
!
!      ! Read QL from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( QL, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the QL:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the QL:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the QL:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: TAUCLI
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "TAUCLI"
!
!      ! Read TAUCLI from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( TAUCLI, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the TAUCLI:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the TAUCLI:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the TAUCLI:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: TAUCLW
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "TAUCLW"
!
!      ! Read TAUCLW from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( TAUCLW, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the TAUCLW:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the TAUCLW:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the TAUCLW:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !=================================================================
!      ! Cleanup and quit
!      !=================================================================
!
!      ! Close netCDF file
!      CALL NcCl( fId )
!
!      END SUBROUTINE READ_FROM_NETCDF_FILE
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_netcdf_file
!!
!! !DESCRIPTION: Routine to read variables and attributes from a netCDF
!!  file.  This routine was automatically generated by the Perl script
!!  NcdfUtilities/perl/ncCodeRead.
!!\
!!\
!! !INTERFACE:
!!
!      SUBROUTINE READ_FROM_NETCDF_FILE( fId )
!!
!! !USES:
!!
!      ! Modules for netCDF read
!      USE m_netcdf_io_open
!      USE m_netcdf_io_get_dimlen
!      USE m_netcdf_io_read
!      USE m_netcdf_io_readattr
!      USE m_netcdf_io_close
!
!      IMPLICIT NONE
!
!#     include "netcdf.inc"
!!
!! !OUTPUT PARAMETERS:
!!   
!      INTEGER, INTENT(INOUT) :: fId    ! netCDF file ID
!!
!! !REMARKS:
!!  Assumes that you have:
!!  (1) A netCDF library (either v3 or v4) installed on your system
!!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!!                                                                             .
!!  Although this routine was generated automatically, some further
!!  hand-editing may be required (i.e. to  specify the size of parameters, 
!!  and/or to assign values to variables.  Also, you can decide how to handle
!!  the variable attributes (or delete calls for reading attributes that you
!!  do not need).
!!
!! !REVISION HISTORY:
!!  30 Jan 2012 - R. Yantosca - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      !=================================================================
!      ! Variable declarations
!      !=================================================================
!
!      ! Data arrays
!      REAL*4             :: lon(IIPAR)
!      REAL*4             :: lat(JJPAR)
!      REAL*4             :: lev(LLPAR)
!      INTEGER            :: time(1)
!      REAL*4             :: CMFMC(IIPAR,JJPAR,LLPAR+1,1)
!      REAL*4             :: DTRAIN(IIPAR,JJPAR,LLPAR,1)
!      REAL*4             :: OMEGA(IIPAR,JJPAR,LLPAR,1)
!      REAL*4             :: RH(IIPAR,JJPAR,LLPAR,1)
!      REAL*4             :: U(IIPAR,JJPAR,LLPAR,1)
!      REAL*4             :: V(IIPAR,JJPAR,LLPAR,1)
!
!      ! Character strings
!      CHARACTER(LEN=255) :: nc_file            ! netCDF file name
!      CHARACTER(LEN=255) :: v_name             ! netCDF variable name 
!      CHARACTER(LEN=255) :: a_name             ! netCDF attribute name
!      CHARACTER(LEN=255) :: a_val              ! netCDF attribute value
!
!      ! Arrays for netCDF start and count values
!      INTEGER            :: st1d(1), ct1d(1)   ! For 1D arrays    
!      INTEGER            :: st2d(2), ct2d(2)   ! For 2D arrays 
!      INTEGER            :: st3d(3), ct3d(3)   ! For 3D arrays 
!      INTEGER            :: st4d(4), ct4d(4)   ! For 4D arrays 
!      INTEGER            :: st5d(5), ct5d(5)   ! For 5D arrays 
!      INTEGER            :: st6d(6), ct6d(6)   ! For 6D arrays 
!
!      !=================================================================
!      ! Open and read data from the netCDF file
!      !=================================================================
!
!      ! Open netCDF file
!      nc_file = 'GEOS572.YYYYMMDD.A3dyn.4x5.nc'
!      CALL Ncop_Rd( fId, TRIM(nc_file) )
!
!      !----------------------------------------
!      ! VARIABLE: lon
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lon"
!
!      ! Read lon from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ IIPAR /)
!      CALL NcRd( lon, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lon:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lon:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: lat
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lat"
!
!      ! Read lat from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ JJPAR /)
!      CALL NcRd( lat, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lat:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lat:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: lev
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lev"
!
!      ! Read lev from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ LLPAR /)
!      CALL NcRd( lev, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lev:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lev:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: time
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "time"
!
!      ! Read time from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ 1 /)
!      CALL NcRd( time, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the time:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:delta_t attribute
!      a_name = "delta_t"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:begin_date attribute
!      a_name = "begin_date"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:begin_time attribute
!      a_name = "begin_time"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:time_increment attribute
!      a_name = "time_increment"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: CMFMC
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "CMFMC"
!
!      ! Read CMFMC from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR+1, 1 /)
!      CALL NcRd( CMFMC, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the CMFMC:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the CMFMC:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the CMFMC:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: DTRAIN
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "DTRAIN"
!
!      ! Read DTRAIN from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( DTRAIN, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the DTRAIN:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the DTRAIN:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the DTRAIN:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: OMEGA
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "OMEGA"
!
!      ! Read OMEGA from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( OMEGA, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the OMEGA:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the OMEGA:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the OMEGA:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: RH
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "RH"
!
!      ! Read RH from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( RH, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the RH:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the RH:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the RH:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: U
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "U"
!
!      ! Read U from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( U, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the U:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the U:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the U:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: V
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "V"
!
!      ! Read V from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( V, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the V:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the V:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the V:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !=================================================================
!      ! Cleanup and quit
!      !=================================================================
!
!      ! Close netCDF file
!      CALL NcCl( fId )
!
!      END SUBROUTINE READ_FROM_NETCDF_FILE
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_netcdf_file
!!
!! !DESCRIPTION: Routine to read variables and attributes from a netCDF
!!  file.  This routine was automatically generated by the Perl script
!!  NcdfUtilities/perl/ncCodeRead.
!!\
!!\
!! !INTERFACE:
!!
!      SUBROUTINE READ_FROM_NETCDF_FILE( fId )
!!
!! !USES:
!!
!      ! Modules for netCDF read
!      USE m_netcdf_io_open
!      USE m_netcdf_io_get_dimlen
!      USE m_netcdf_io_read
!      USE m_netcdf_io_readattr
!      USE m_netcdf_io_close
!
!      IMPLICIT NONE
!
!#     include "netcdf.inc"
!!
!! !OUTPUT PARAMETERS:
!!   
!      INTEGER, INTENT(INOUT) :: fId    ! netCDF file ID
!!
!! !REMARKS:
!!  Assumes that you have:
!!  (1) A netCDF library (either v3 or v4) installed on your system
!!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!!                                                                             .
!!  Although this routine was generated automatically, some further
!!  hand-editing may be required (i.e. to  specify the size of parameters, 
!!  and/or to assign values to variables.  Also, you can decide how to handle
!!  the variable attributes (or delete calls for reading attributes that you
!!  do not need).
!!
!! !REVISION HISTORY:
!!  30 Jan 2012 - R. Yantosca - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      !=================================================================
!      ! Variable declarations
!      !=================================================================
!
!      ! Data arrays
!      REAL*4             :: lon(IIPAR)
!      REAL*4             :: lat(JJPAR)
!      REAL*4             :: lev(LLPAR)
!      INTEGER            :: time(1)
!      REAL*4             :: DQRCU(IIPAR,JJPAR,LLPAR,1)
!      REAL*4             :: DQRLSAN(IIPAR,JJPAR,LLPAR,1)
!      REAL*4             :: REEVAPCN(IIPAR,JJPAR,LLPAR,1)
!      REAL*4             :: REEVAPLS(IIPAR,JJPAR,LLPAR,1)
!
!      ! Character strings
!      CHARACTER(LEN=255) :: nc_file            ! netCDF file name
!      CHARACTER(LEN=255) :: v_name             ! netCDF variable name 
!      CHARACTER(LEN=255) :: a_name             ! netCDF attribute name
!      CHARACTER(LEN=255) :: a_val              ! netCDF attribute value
!
!      ! Arrays for netCDF start and count values
!      INTEGER            :: st1d(1), ct1d(1)   ! For 1D arrays    
!      INTEGER            :: st2d(2), ct2d(2)   ! For 2D arrays 
!      INTEGER            :: st3d(3), ct3d(3)   ! For 3D arrays 
!      INTEGER            :: st4d(4), ct4d(4)   ! For 4D arrays 
!      INTEGER            :: st5d(5), ct5d(5)   ! For 5D arrays 
!      INTEGER            :: st6d(6), ct6d(6)   ! For 6D arrays 
!
!      !=================================================================
!      ! Open and read data from the netCDF file
!      !=================================================================
!
!      ! Open netCDF file
!      nc_file = 'GEOS572.YYYYMMDD.A3mstC.4x5.nc'
!      CALL Ncop_Rd( fId, TRIM(nc_file) )
!
!      !----------------------------------------
!      ! VARIABLE: lon
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lon"
!
!      ! Read lon from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ IIPAR /)
!      CALL NcRd( lon, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lon:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lon:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: lat
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lat"
!
!      ! Read lat from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ JJPAR /)
!      CALL NcRd( lat, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lat:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lat:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: lev
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lev"
!
!      ! Read lev from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ LLPAR /)
!      CALL NcRd( lev, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lev:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lev:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: time
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "time"
!
!      ! Read time from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ 1 /)
!      CALL NcRd( time, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the time:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:delta_t attribute
!      a_name = "delta_t"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:begin_date attribute
!      a_name = "begin_date"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:begin_time attribute
!      a_name = "begin_time"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:time_increment attribute
!      a_name = "time_increment"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: DQRCU
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "DQRCU"
!
!      ! Read DQRCU from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( DQRCU, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the DQRCU:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the DQRCU:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the DQRCU:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: DQRLSAN
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "DQRLSAN"
!
!      ! Read DQRLSAN from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( DQRLSAN, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the DQRLSAN:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the DQRLSAN:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the DQRLSAN:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: REEVAPCN
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "REEVAPCN"
!
!      ! Read REEVAPCN from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( REEVAPCN, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the REEVAPCN:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the REEVAPCN:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the REEVAPCN:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: REEVAPLS
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "REEVAPLS"
!
!      ! Read REEVAPLS from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( REEVAPLS, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the REEVAPLS:long_name attribute
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_netcdf_file
!!
!! !DESCRIPTION: Routine to read variables and attributes from a netCDF
!!  file.  This routine was automatically generated by the Perl script
!!  NcdfUtilities/perl/ncCodeRead.
!!\
!!\
!! !INTERFACE:
!!
!      SUBROUTINE READ_FROM_NETCDF_FILE( fId )
!!
!! !USES:
!!
!      ! Modules for netCDF read
!      USE m_netcdf_io_open
!      USE m_netcdf_io_get_dimlen
!      USE m_netcdf_io_read
!      USE m_netcdf_io_readattr
!      USE m_netcdf_io_close
!
!      IMPLICIT NONE
!
!#     include "netcdf.inc"
!!
!! !OUTPUT PARAMETERS:
!!   
!      INTEGER, INTENT(INOUT) :: fId    ! netCDF file ID
!!
!! !REMARKS:
!!  Assumes that you have:
!!  (1) A netCDF library (either v3 or v4) installed on your system
!!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!!                                                                             .
!!  Although this routine was generated automatically, some further
!!  hand-editing may be required (i.e. to  specify the size of parameters, 
!!  and/or to assign values to variables.  Also, you can decide how to handle
!!  the variable attributes (or delete calls for reading attributes that you
!!  do not need).
!!
!! !REVISION HISTORY:
!!  30 Jan 2012 - R. Yantosca - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      !=================================================================
!      ! Variable declarations
!      !=================================================================
!
!      ! Data arrays
!      REAL*4             :: lon(IIPAR)
!      REAL*4             :: lat(JJPAR)
!      REAL*4             :: lev(LLPAR+1)
!      INTEGER            :: time(1)
!      REAL*4             :: PFICU(IIPAR,JJPAR,)
!      REAL*4             :: PFILSAN(IIPAR,JJPAR,)
!      REAL*4             :: PFLCU(IIPAR,JJPAR,)
!      REAL*4             :: PFLLSAN(IIPAR,JJPAR,)
!
!      ! Character strings
!      CHARACTER(LEN=255) :: nc_file            ! netCDF file name
!      CHARACTER(LEN=255) :: v_name             ! netCDF variable name 
!      CHARACTER(LEN=255) :: a_name             ! netCDF attribute name
!      CHARACTER(LEN=255) :: a_val              ! netCDF attribute value
!
!      ! Arrays for netCDF start and count values
!      INTEGER            :: st1d(1), ct1d(1)   ! For 1D arrays    
!      INTEGER            :: st2d(2), ct2d(2)   ! For 2D arrays 
!      INTEGER            :: st3d(3), ct3d(3)   ! For 3D arrays 
!      INTEGER            :: st4d(4), ct4d(4)   ! For 4D arrays 
!      INTEGER            :: st5d(5), ct5d(5)   ! For 5D arrays 
!      INTEGER            :: st6d(6), ct6d(6)   ! For 6D arrays 
!
!      !=================================================================
!      ! Open and read data from the netCDF file
!      !=================================================================
!
!      ! Open netCDF file
!      nc_file = 'GEOS572.YYYYMMDD.A3mstE.4x5.nc'
!      CALL Ncop_Rd( fId, TRIM(nc_file) )
!
!      !----------------------------------------
!      ! VARIABLE: lon
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lon"
!
!      ! Read lon from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ IIPAR /)
!      CALL NcRd( lon, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lon:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lon:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: lat
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lat"
!
!      ! Read lat from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ JJPAR /)
!      CALL NcRd( lat, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lat:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lat:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: lev
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lev"
!
!      ! Read lev from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ LLPAR+1 /)
!      CALL NcRd( lev, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lev:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lev:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: time
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "time"
!
!      ! Read time from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ 1 /)
!      CALL NcRd( time, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the time:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:delta_t attribute
!      a_name = "delta_t"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:begin_date attribute
!      a_name = "begin_date"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:begin_time attribute
!      a_name = "begin_time"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:time_increment attribute
!      a_name = "time_increment"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: PFICU
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "PFICU"
!
!      ! Read PFICU from file
!      st3d   = (/ 1, 1, 1 /)
!      ct3d   = (/ IIPAR, JJPAR,  /)
!      CALL NcRd( PFICU, fId, TRIM(v_name), st3d, ct3d )
!
!      ! Read the PFICU:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the PFICU:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the PFICU:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: PFILSAN
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "PFILSAN"
!
!      ! Read PFILSAN from file
!      st3d   = (/ 1, 1, 1 /)
!      ct3d   = (/ IIPAR, JJPAR,  /)
!      CALL NcRd( PFILSAN, fId, TRIM(v_name), st3d, ct3d )
!
!      ! Read the PFILSAN:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the PFILSAN:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the PFILSAN:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: PFLCU
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "PFLCU"
!
!      ! Read PFLCU from file
!      st3d   = (/ 1, 1, 1 /)
!      ct3d   = (/ IIPAR, JJPAR,  /)
!      CALL NcRd( PFLCU, fId, TRIM(v_name), st3d, ct3d )
!
!      ! Read the PFLCU:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the PFLCU:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the PFLCU:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: PFLLSAN
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "PFLLSAN"
!
!      ! Read PFLLSAN from file
!      st3d   = (/ 1, 1, 1 /)
!      ct3d   = (/ IIPAR, JJPAR,  /)
!      CALL NcRd( PFLLSAN, fId, TRIM(v_name), st3d, ct3d )
!
!      ! Read the PFLLSAN:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the PFLLSAN:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the PFLLSAN:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !=================================================================
!      ! Cleanup and quit
!      !=================================================================
!
!      ! Close netCDF file
!      CALL NcCl( fId )
!
!      END SUBROUTINE READ_FROM_NETCDF_FILE
!!EOC
!
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the REEVAPLS:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the REEVAPLS:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !=================================================================
!      ! Cleanup and quit
!      !=================================================================
!
!      ! Close netCDF file
!      CALL NcCl( fId )
!
!      END SUBROUTINE READ_FROM_NETCDF_FILE
!!EOC
!!------------------------------------------------------------------------------
!!          Harvard University Atmospheric Chemistry Modeling Group            !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: read_netcdf_file
!!
!! !DESCRIPTION: Routine to read variables and attributes from a netCDF
!!  file.  This routine was automatically generated by the Perl script
!!  NcdfUtilities/perl/ncCodeRead.
!!\
!!\
!! !INTERFACE:
!!
!      SUBROUTINE READ_FROM_NETCDF_FILE( fId )
!!
!! !USES:
!!
!      ! Modules for netCDF read
!      USE m_netcdf_io_open
!      USE m_netcdf_io_get_dimlen
!      USE m_netcdf_io_read
!      USE m_netcdf_io_readattr
!      USE m_netcdf_io_close
!
!      IMPLICIT NONE
!
!#     include "netcdf.inc"
!!
!! !OUTPUT PARAMETERS:
!!   
!      INTEGER, INTENT(INOUT) :: fId    ! netCDF file ID
!!
!! !REMARKS:
!!  Assumes that you have:
!!  (1) A netCDF library (either v3 or v4) installed on your system
!!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!!                                                                             .
!!  Although this routine was generated automatically, some further
!!  hand-editing may be required (i.e. to  specify the size of parameters, 
!!  and/or to assign values to variables.  Also, you can decide how to handle
!!  the variable attributes (or delete calls for reading attributes that you
!!  do not need).
!!
!! !REVISION HISTORY:
!!  30 Jan 2012 - R. Yantosca - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      !=================================================================
!      ! Variable declarations
!      !=================================================================
!
!      ! Data arrays
!      REAL*4             :: lon(IIPAR)
!      REAL*4             :: lat(JJPAR)
!      REAL*4             :: lev(LLPAR)
!      INTEGER            :: time(1)
!      REAL*4             :: PS(IIPAR,JJPAR,1)
!      REAL*4             :: PV(IIPAR,JJPAR,LLPAR,1)
!
!      ! Character strings
!      CHARACTER(LEN=255) :: nc_file            ! netCDF file name
!      CHARACTER(LEN=255) :: v_name             ! netCDF variable name 
!      CHARACTER(LEN=255) :: a_name             ! netCDF attribute name
!      CHARACTER(LEN=255) :: a_val              ! netCDF attribute value
!
!      ! Arrays for netCDF start and count values
!      INTEGER            :: st1d(1), ct1d(1)   ! For 1D arrays    
!      INTEGER            :: st2d(2), ct2d(2)   ! For 2D arrays 
!      INTEGER            :: st3d(3), ct3d(3)   ! For 3D arrays 
!      INTEGER            :: st4d(4), ct4d(4)   ! For 4D arrays 
!      INTEGER            :: st5d(5), ct5d(5)   ! For 5D arrays 
!      INTEGER            :: st6d(6), ct6d(6)   ! For 6D arrays 
!
!      !=================================================================
!      ! Open and read data from the netCDF file
!      !=================================================================
!
!      ! Open netCDF file
!      nc_file = 'GEOS572.YYYYMMDD.I3.4x5.nc'
!      CALL Ncop_Rd( fId, TRIM(nc_file) )
!
!      !----------------------------------------
!      ! VARIABLE: lon
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lon"
!
!      ! Read lon from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ IIPAR /)
!      CALL NcRd( lon, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lon:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lon:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: lat
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lat"
!
!      ! Read lat from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ JJPAR /)
!      CALL NcRd( lat, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lat:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lat:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: lev
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "lev"
!
!      ! Read lev from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ LLPAR /)
!      CALL NcRd( lev, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the lev:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the lev:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: time
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "time"
!
!      ! Read time from file
!      st1d   = (/ 1 /)
!      ct1d   = (/ 1 /)
!      CALL NcRd( time, fId, TRIM(v_name), st1d, ct1d )
!
!      ! Read the time:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:delta_t attribute
!      a_name = "delta_t"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:begin_date attribute
!      a_name = "begin_date"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:begin_time attribute
!      a_name = "begin_time"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the time:time_increment attribute
!      a_name = "time_increment"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: PS
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "PS"
!
!      ! Read PS from file
!      st3d   = (/ 1, 1, 1 /)
!      ct3d   = (/ IIPAR, JJPAR, 1 /)
!      CALL NcRd( PS, fId, TRIM(v_name), st3d, ct3d )
!
!      ! Read the PS:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the PS:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the PS:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !----------------------------------------
!      ! VARIABLE: PV
!      !----------------------------------------
!
!      ! Read  from file
!      v_name = "PV"
!
!      ! Read PV from file
!      st4d   = (/ 1, 1, 1, 1 /)
!      ct4d   = (/ IIPAR, JJPAR, LLPAR, 1 /)
!      CALL NcRd( PV, fId, TRIM(v_name), st4d, ct4d )
!
!      ! Read the PV:long_name attribute
!      a_name = "long_name"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the PV:units attribute
!      a_name = "units"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      ! Read the PV:gamap_category attribute
!      a_name = "gamap_category"
!      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )
!
!      !=================================================================
!      ! Cleanup and quit
!      !=================================================================
!
!      ! Close netCDF file
!      CALL NcCl( fId )
!
!      END SUBROUTINE READ_FROM_NETCDF_FILE
!EOC
END MODULE GEOS57_READ_MOD
