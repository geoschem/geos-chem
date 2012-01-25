!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: TestNcdfUtil.F90
!
! !DESCRIPTION: Program TestNcdfUtilities.F90 is the standalone driver that
!  tests if the libNcUtils.a file was built correctly.
!\\
!\\
! !INTERFACE:
!
PROGRAM TestNcdfUtil
!
! !USES: 
!
!
! !USES:
!
  ! Modules for netCDF write
  USE m_netcdf_io_create
  USE m_netcdf_io_define
  USE m_netcdf_io_write
  USE m_netcdf_io_close
   
  ! Modules for netCDF read
  USE m_netcdf_io_open
  USE m_netcdf_io_close     
  USE m_netcdf_io_get_dimlen
  USE m_netcdf_io_read
  USE m_netcdf_io_readattr
  
  IMPLICIT NONE

  ! netCDF include files
# include "netcdf.inc"

!
! !BUGS:  
!  None known at this time
!
! !SEE ALSO: 
!  m_do_err_out.F90
!  m_netcdf_io_checks.F90
!  m_netcdf_io_close.F90
!  m_netcdf_io_create.F90
!  m_netcdf_io_define.F90
!  m_netcdf_io_get_dimlen.F90
!  m_netcdf_io_handle_err.F90
!  m_netcdf_io_open.F90
!  m_netcdf_io_read.F90
!  m_netcdf_io_write.F90
!
! !SYSTEM ROUTINES: 
!  None
!
! !REMARKS:
!  netCDF library modules originally written by Jules Kouatchou, GSFC
!  and re-packaged into NcdfUtilities by Bob Yantosca, Harvard Univ.

! !REVISION HISTORY: 
!  03 Jul 2008 - R. Yantosca (Harvard University) - Initial version
!  24 Jan 2012 - R. Yantosca - Modified to write COARDS-compliant output
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  ! Global private variables
  INTEGER, PARAMETER :: ILONG = 72               ! # of longitude grid points
  INTEGER, PARAMETER :: ILAT  = 46               ! # of latitude  grid points
  INTEGER, PARAMETER :: IVERT = 55               ! # of altitude  levels 
  INTEGER, PARAMETER :: ITIME = 1                ! # of times
  INTEGER            :: pCt                      ! # of passed tests
  INTEGER            :: tCt                      ! # of total tests
  INTEGER            :: I                        ! Loop index
  INTEGER            :: longdeg, latdeg          ! For longdat, latdat
  REAL*8             :: longDat(ILONG)           ! Longitude data
  REAL*8             :: latDat (ILONG)           ! Latitude data
  REAL*8             :: levDat (IVERT)           ! Altitude data
  INTEGER            :: timeDat(ITIME)           ! Time data

  ! Initialize
  pCt = 0
  tCt = 0

  ! Longitude data
  longdeg = 360.0 / REAL( ILONG )
  if ( mod( 360, ILONG) /= 0 ) longdeg = longdeg + 1
  do i = 1, ILONG
     longDat(i) = i*longdeg
  enddo

  ! Writing latitude data point
  latdeg  = 180.0 / REAL( ILAT )
  if ( mod( 180, ILAT ) /= 0 ) latdeg = latdeg + 1
  do i = 1, ilong
     latDat(i) = -90 + (i-0.5)*latdeg
  enddo

  ! Pressure
  do i = 1, IVERT
     levDat(i) = 1000.00 - (i-1)*(920.00/IVERT) 
  enddo

  ! Time data
  do i = 1, ITIME
     timeDat(i) = 0
  enddo

  ! Echo info
  WRITE( 6, '(a)' ) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  WRITE( 6, '(a)' ) '%%%  Testing libNcdfUtilities.a  %%%'
  WRITE( 6, '(a)' ) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

  ! Create a netCDF file
  CALL TestNcdfCreate

  ! And try to read it back
  CALL TestNcdfRead
!BOC

CONTAINS

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TestNcdfCreate
!
! !DESCRIPTION: Subroutine TestNcdfCreate creates a netCDF file
!  named \texttt{my\_filename.nc} with the following variables:
!
!  \begin{description}
!  \item[PSF] Surface pressure (2D variable)
!  \item[KEL] Temperature (3D variable)
!  \end{description}
!
!  Fake values are used for the data.  An unlimited dimension is employed
!  to write out several records of kel.  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TestNcdfCreate
!
! !REVISION HISTORY: 
!  03 Jul 2008 - R. Yantosca (Harvard University) - Initial version
!  24 Jan 2012 - R. Yantosca - Modified to provide COARDS-compliant output
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! For netCDF file I/O
    INTEGER             :: idLon,    idLat,   idLev,   idTime
    INTEGER             :: fId,      vId,     omode,   i
    INTEGER             :: ct1d(1),  ct3d(3), ct4d(4)
    INTEGER             :: st1d(1),  st3d(3), st4d(4)
    INTEGER             :: var1(1),  var3(3), var4(4)
    CHARACTER(LEN=255)  :: units,    delta_t, begin_d
    CHARACTER(LEN=255)  :: begin_t,  incr

    ! "Fake" data arrays
    REAL*4              :: PS( ILONG, ILAT,        ITIME )  ! surface pressure
    REAL*4              :: T ( ILONG, ILAT, IVERT, ITIME )  ! temperature
    
    !=========================================================================
    ! Create the netCDF file
    !=========================================================================

    ! Echo info
    WRITE( 6, '(a)' ) '=== Begin netCDF file creation test ==='

    CALL NcCr_Wr( fId, 'my_filename.nc' )
    
    ! Turn filling off
    CALL NcSetFill( fId, NF_NOFILL, omode )

    !=========================================================================
    ! Define the dimensions
    !=========================================================================
    
    ! Longitude dimension                      
    WRITE( 6, '(a)' ) 'Writing lon  (dim)      to netCDF file'    
    CALL NcDef_Dimension( fId, 'lon',  ILONG, idLon )
    
    ! Latitude dimension
    WRITE( 6, '(a)' ) 'Writing lat  (dim)      to netCDF file'    
    CALL NcDef_Dimension( fId, 'lat',  ILAT , idLat )
    
    ! Altitude dimension
    WRITE( 6, '(a)' ) 'Writing lev  (dim)       to netCDF file'    
    CALL NcDef_Dimension( fId, 'lev',  IVERT, idLev )

    ! Altitude dimension
    WRITE( 6, '(a)' ) 'Writing time (dim)       to netCDF file'    
    CALL NcDef_Dimension( fId, 'time', ITIME, idTime )
    
    !=========================================================================
    ! Define the variables and variable attributes 
    ! for COARDS compliance and GAMAP compliance
    !=========================================================================
    CALL NcDef_Glob_Attributes( fId, 'title',       'NcdfUtilities test file' )
    CALL NcDef_Glob_Attributes( fId, 'history',     'test file - 24 Jan 2011' )
    CALL NcDef_Glob_Attributes( fId, 'conventions', 'COARDS'                  )
    CALL NcDef_Glob_Attributes( fId, 'model',       'GEOS4'                   )
    CALL NcDef_Glob_Attributes( fId, 'nlayers',     '55'                      )
    CALL NcDef_Glob_Attributes( fId, 'start_date',  '20110101'                )
    CALL NcDef_Glob_Attributes( fId, 'start_time',  '00:00:00.0'              )
    CALL NcDef_Glob_Attributes( fId, 'end_date',    '20110101'                )
    CALL NcDef_Glob_Attributes( fId, 'end_time',    '23:59:59.0'              )
    CALL NcDef_Glob_Attributes( fId, 'delta_lon',   '5'                       )
    CALL NcDef_Glob_Attributes( fId, 'delta_lat',   '4'                       )
    CALL NcDef_Glob_Attributes( fId, 'delta_time',  '000000'                  )
    CALL NcDef_Glob_Attributes( fId, 'format',      'netCDF-3'                )

    !=========================================================================
    ! Define the variables and variable attributes
    !=========================================================================

    ! Define longitude variable
    var1 = (/ idLon /)
    CALL NcDef_Variable( fId, 'lon', NF_DOUBLE, 1, var1, vId )
    CALL NcDef_Var_Attributes( fId, vId,  'long_name', 'Longitude'   )
    CALL NcDef_Var_Attributes( fId, vId,  'units',     'degree_east' )
  
    ! Define latitude variable
    var1 = (/ idLat /)
    CALL NcDef_Variable( fId, 'lat', NF_DOUBLE, 1, var1, vId )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'Latitude'     )
    CALL NcDef_Var_Attributes( fId, vId, 'units',     'degree_north' )
    
    ! Define vertical (pressure) variable
    var1 = (/ idLev /)
    CALL NcDef_Variable( fId, 'lev', NF_DOUBLE, 1, var1, vId )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'Pressure' )
    CALL NcDef_Var_Attributes( fId, vId, 'units',     'hPa'      )
    
    ! Time index array (hardwire date to 2011/01/01)
    var1    = (/ idTime /)
    vId     = vId + 1
    units   = 'minutes since 2011-01-01 00:00:00 GMT'
    delta_t = '0000-00-00 00:00:00'
    begin_d = '20110101'
    begin_t = '000000'
    incr    = '000000'
    CALL NcDef_Variable      ( fId, 'time', NF_INT,  1, var1, vId           )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name',      'time'           )
    CALL NcDef_Var_Attributes( fId, vId, 'units',          TRIM( units   )  ) 
    CALL NcDef_Var_Attributes( fId, vId, 'delta_t',        TRIM( delta_t )  ) 
    CALL NcDef_Var_Attributes( fId, vId, 'begin_date',     TRIM( begin_d )  )
    CALL NcDef_Var_Attributes( fId, vId, 'begin_time',     TRIM( begin_t )  )
    CALL NcDef_Var_Attributes( fId, vId, 'time_increment', TRIM( incr    )  )


    ! Define surface pressure variable
    var3 = (/ idLon, idLat, idTime /)
    CALL NcDef_Variable      ( fId, 'PS', NF_FLOAT, 3, var3, vId )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name',      'Surface Pressure' )
    CALL NcDef_Var_Attributes( fId, vId, 'units',          'hPa'              )
    CALL NcDef_Var_Attributes( fId, vId, 'gamap_category', 'GMAO-2D'          )
    
    ! Define 
    var4 = (/ idLon, idLat, idLev, idTime /)
    CALL NcDef_Variable      ( fId, 'T', NF_FLOAT, 4, var4, vId )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name',      'Temperature'      )
    CALL NcDef_Var_Attributes( fId, vId, 'units',          'K'                )
    CALL NcDef_Var_Attributes( fId, vId, 'gamap_category', 'GMAO-3D$'         )

    
    !=========================================================================
    ! %%% END OF DEFINITION SECTION %%%
    ! %%% NOW WRITE DATA TO FILE    %%%
    !=========================================================================
    CALL NcEnd_def( fId )
    
    ! Write longitude
    WRITE( 6, '(a)' ) 'Writing lon  (1D array) to netCDF file'    
    st1d = (/ 1     /)
    ct1d = (/ ILONG /)
    CALL NcWr( longDat, fId, 'lon', st1d, ct1d )
    
    ! Write latitude
    WRITE( 6, '(a)' ) 'Writing lat  (1D array) to netCDF file'    
    st1d = (/ 1    /)
    ct1d = (/ ILAT /)
    CALL NcWr( latDat, fId, 'lat', st1d, ct1d )
    
    ! Write pressure levels
    WRITE( 6, '(a)' ) 'Writing lev  (1D array) to netCDF file'    
    st1d = (/ 1     /)
    ct1d = (/ IVERT /)
    CALL NcWr( levDat, fId, 'lev', st1d, ct1d )

    ! Write pressure levels
    WRITE( 6, '(a)' ) 'Writing time (1D array) to netCDF file'    
    st1d = (/ 1     /)
    ct1d = (/ ITIME /)
    CALL NcWr( timeDat, fId, 'time', st1d, ct1d )
    
    ! Write surface pressure (w/ fake values)
    WRITE( 6, '(a)' ) 'Writing PS   (3D array) to netCDF file'  
    PS    = 1e0
    st3d = (/ 1,     1,   1      /)
    ct3d = (/ ILONG, ILAT, ITIME /)
    CALL NcWr( PS, fId, 'PS', st3d, ct3d )
    
    ! Write temperature (w/ fake values)
    WRITE( 6, '(a)' ) 'Writing T    (4D array) to netCDF file'      
    T    = 1e0
    st4d = (/ 1,     1,    1,     1     /)
    ct4d = (/ ILONG, ILAT, IVERT, ITIME /)
    CALL NcWr( T, fId, 'T', st4d, ct4d )
    
    !=========================================================================
    ! Close the netCDF file
    !=========================================================================
    CALL NcCl( fId )

    ! Echo info
    WRITE( 6, '(a)' ) '=== End netCDF file creation test ==='

  END SUBROUTINE TestNcdfCreate
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TestNcdfRead
!
! !DESCRIPTION: Routine TestNcdfRead extracts the following fields from
!  the netCDF file \texttt{my\_filename.nc}:
!
!  \begin{description}
!  \item[PSF] Surface pressure (2D variable)
!  \item[KEL] Temperature (3D variable).
!  \end{description}
!
!  Note that the file \texttt{my\_filename.nc} was created with fake data
!  values by subroutine TestNcdfCreate.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TestNcdfRead
!
! !REVISION HISTORY: 
!  03 Jul 2008 - R. Yantosca (Harvard University) - Initial version
!  24 Jan 2012 - R. Yantosca - Modified to provide COARDS-compliant output
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    ! Scalars
    INTEGER              :: fId,     rc,      XDim
    INTEGER              :: YDim,    ZDim,    TDim
    INTEGER              :: ct1d(1), ct3d(3), ct4d(4)
    INTEGER              :: st1d(1), st3d(3), st4d(4)
    CHARACTER(LEN=255)   :: attValue

    ! Arrays
    REAL*8,  ALLOCATABLE :: lon(:),  lat(:),  lev(:)
    INTEGER, ALLOCATABLE :: time(:)
    REAL*4,  ALLOCATABLE :: PS(:,:,:)
    REAL*4,  ALLOCATABLE :: T(:,:,:,:)
  
    !=========================================================================
    ! Open the netCDF file
    !=========================================================================

    ! Echo info
    WRITE( 6, '(a)' ) '=== Begin netCDF file reading test ==='

    CALL Ncop_Rd( fId, 'my_filename.nc' )
    
    !=========================================================================
    ! Get the dimensions
    !=========================================================================
    CALL Ncget_Dimlen( fId, 'lon',  XDim )
    CALL Ncget_Dimlen( fId, 'lat',  YDim )
    CALL Ncget_Dimlen( fId, 'lev',  ZDim )
    CALL Ncget_Dimlen( fId, 'time', TDim )

    rc = XDim - ILONG
    CALL Check( 'Reading lon  (dim)   back from netCDF file', rc, pCt, tCt )

    rc = YDim - ILAT
    CALL Check( 'Reading lat  (dim)   back from netCDF file', rc, pCt, tCt )

    rc = ZDim - IVERT
    CALL Check( 'Reading lev  (dim)   back from netCDF file', rc, pCt, tCt ) 

    rc = TDim - ITIME
    CALL Check( 'Reading time (dim)   back from netCDF file', rc, pCt, tCt ) 

    !=========================================================================
    ! Read the LON variable
    !=========================================================================

    ! Read data
    ALLOCATE( lon( XDim ) )
    st1d = (/ 1    /)
    ct1d = (/ XDim /)
    CALL NcRd( lon, fId, 'lon', st1d, ct1d )

    ! Equality test
    rc = SUM( lon - longDat )
    CALL Check( 'Reading lon  (array) back from netCDF file', rc, pCt, tCt )

    !=========================================================================
    ! Read the LAT variable
    !=========================================================================

    ! Read data
    ALLOCATE( lat( YDim ) )
    st1d = (/ 1    /)
    ct1d = (/ YDim /)
    CALL NcRd( lat, fId, 'lat', st1d, ct1d )

    ! Equality test
    rc = SUM( lat - latDat )
    CALL Check( 'Reading lat  (array) back from netCDF file', rc, pCt, tCt )

    !=========================================================================
    ! Read the LEV variable
    !=========================================================================

    ! Read data
    ALLOCATE( lev( ZDim ) )
    st1d = (/ 1    /)
    ct1d = (/ ZDim /)
    CALL NcRd( lev, fId, 'lev', st1d, ct1d )

    ! Equality test
    rc = SUM( lev - levDat )
    CALL Check( 'Reading lev  (array) back from netCDF file', rc, pCt, tCt )

    !=========================================================================
    ! Read the TIME variable
    !=========================================================================

    ! Read data
    ALLOCATE( time( TDim ) )
    st1d = (/ 1    /)
    ct1d = (/ Tdim /)
    CALL NcRd( time, fId, 'time', st1d, ct1d )

    ! Equality test
    rc = SUM( time - timeDat )
    CALL Check( 'Reading time (array) back from netCDF file', rc, pCt, tCt )

    !=========================================================================
    ! Read the PS variable
    !=========================================================================

    ! Read data    
    ALLOCATE( ps( XDim, YDim, TDim ) )
    st3d = (/ 1,    1,    1    /)
    ct3d = (/ XDim, YDim, TDim /)
    CALL NcRd( ps, fId, 'PS', st3d, ct3d )

    ! Equality test
    rc = SUM( PS ) - SIZE( PS )
    CALL Check( 'Reading PS           back from netCDF file', rc, pCt, tCt )
    
    ! Read units attribute
    CALL NcGet_Var_Attributes( fId, 'PS', 'units', attValue )
    rc = ( TRIM( attValue ) == 'hPa' )
    CALL Check( 'Reading PS:units     back from netCDF file', rc, pCt, tCt )

    ! Read long_name attribute
    CALL NcGet_Var_Attributes( fId, 'PS', 'long_name', attValue )
    rc = ( TRIM( attValue ) == 'Surface Pressure' )
    CALL Check( 'Reading PS:long_name back from netCDF file', rc, pCt, tCt )

    !=========================================================================
    ! Read the T variable
    !=========================================================================

    ! Read data    
    ALLOCATE( T( XDim, YDim, ZDim, TDim ) )
    st4d = (/ 1,    1,    1,    1    /)
    ct4d = (/ XDim, YDim, ZDim, TDim /)
    CALL NcRd( T, fId, 'T', st4d, ct4d )

    ! Equality test
    rc = SUM( t ) - SIZE( t )
    CALL Check( 'Reading T            back from netCDF file', rc, pCt, tCt )

    ! Read units attribute
    CALL NcGet_Var_Attributes( fId, 'PS', 'units', attValue )
    rc = ( TRIM( attValue ) == 'hPa' )
    CALL Check( 'Reading T:units      back from netCDF file', rc, pCt, tCt )

    ! Read long_name attribute
    CALL NcGet_Var_Attributes( fId, 'PS', 'long_name', attValue )
    rc = ( TRIM( attValue ) == 'Surface Pressure' )
    CALL Check( 'Reading T:long_name  back from netCDF file', rc, pCt, tCt )

    !=========================================================================
    ! Read global attributes
    !=========================================================================

    ! Read title attribute
    CALL NcGet_Glob_Attributes( fId, 'title', attValue )
    rc = ( TRIM( attValue ) == 'NcdfUtilities test file' )
    CALL Check( 'Reading title        back from netCDF file', rc, pCt, tCt )

    ! Read start_date
    CALL NcGet_Glob_Attributes( fId, 'start_date', attValue )
    rc = ( TRIM( attValue ) == '20110101' )
    CALL Check( 'Reading start_date   back from netCDF file', rc, pCt, tCt )

    ! Read start_time
    CALL NcGet_Glob_Attributes( fId, 'start_time', attValue )
    rc = ( TRIM( attValue ) == '000000' )
    CALL Check( 'Reading start_time   back from netCDF file', rc, pCt, tCt )

    ! Close netCDF file
    CALL NcCl( fId )

    ! Cleanup
    IF ( ALLOCATED( lon  ) ) DEALLOCATE( lon  )
    IF ( ALLOCATED( lat  ) ) DEALLOCATE( lat  )
    IF ( ALLOCATED( lev  ) ) DEALLOCATE( lev  )
    IF ( ALLOCATED( time ) ) DEALLOCATE( time )
    IF ( ALLOCATED( PS   ) ) DEALLOCATE( PS   )
    IF ( ALLOCATED( T    ) ) DEALLOCATE( T    )

    ! Echo info
    WRITE( 6, '(a)' ) '=== End of netCDF file read test! ==='
    
  END SUBROUTINE TestNcdfRead
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check
!
! !DESCRIPTION: Subroutine that prints "PASSED" or "FAILED" after each test.  
!  Also increments the various counters of passed or failed tests.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Check( msg, rc, passCt, totCt )
!
! !INPUT PARAMETERS:
!    
    CHARACTER(LEN=*), INTENT(IN)    :: msg     ! message to print
    INTEGER,          INTENT(IN)    :: rc      ! Return code
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: passCt  ! # of passed tests
    INTEGER,          INTENT(INOUT) :: totCt   ! # of total tests
!
! !REVISION HISTORY: 
!  03 Jul 2008 - R. Yantosca (Harvard University) - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC  
!
! !LOCAL VARIABLES:
!
    INTEGER :: s

    ! length of message
    s = LEN( msg )

    IF ( rc == 0 ) THEN
       WRITE( 6, '(a)' ) msg // REPEAT( '.', 45-s ) // 'PASSED'
       passCt = passCt + 1
    ELSE
       WRITE( 6, '(a)' ) msg // REPEAT( '.', 45-s ) // 'FAILED'
    ENDIF

    totCt = totCt + 1

  END SUBROUTINE Check
!EOC

END PROGRAM TestNcdfUtil

