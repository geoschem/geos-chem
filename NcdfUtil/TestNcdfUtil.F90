!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  ! Modules for netCDF write
  USE m_netcdf_io_define
  USE m_netcdf_io_create
  USE m_netcdf_io_write

  ! Modules for netCDF read
  USE m_netcdf_io_open
  USE m_netcdf_io_get_dimlen
  USE m_netcdf_io_read
  USE m_netcdf_io_readattr
  USE m_netcdf_io_close

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
!  31 Jan 2012 - R. Yantosca - Bug fix in error checks for attributes
!  14 Jun 2012 - R. Yantosca - Now tests 2D character read/write
!  10 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2017 - R. Yantosca - Now write a test global attribute
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  ! Global private variables
  INTEGER, PARAMETER :: ILONG  = 72              ! # of longitude grid points
  INTEGER, PARAMETER :: ILAT   = 46              ! # of latitude  grid points
  INTEGER, PARAMETER :: IVERT  = 55              ! # of altitude  levels
  INTEGER, PARAMETER :: ITIME  = 1               ! # of times
  INTEGER, PARAMETER :: ICHAR1 = 2               ! # of times
  INTEGER, PARAMETER :: ICHAR2 = 20              ! # of times
  INTEGER            :: pCt                      ! # of passed tests
  INTEGER            :: tCt                      ! # of total tests
  INTEGER            :: I                        ! Loop index
  INTEGER            :: longdeg, latdeg          ! For longdat, latdat
  REAL*8             :: longDat(ILONG)           ! Longitude data
  REAL*8             :: latDat (ILAT )           ! Latitude data
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!  14 Jun 2012 - R. Yantosca - Now writes a 2-D character array
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! For netCDF file I/O
    INTEGER             :: idLon,    idLat,   idLev,   idTime
    INTEGER             :: idChar1,  idChar2
    INTEGER             :: fId,      vId,     omode,   i
    INTEGER             :: ct1d(1),  ct2d(2), ct3d(3), ct4d(4)
    INTEGER             :: st1d(1),  st2d(2), st3d(3), st4d(4)
    INTEGER             :: var1(1),  var2(2), var3(3), var4(4)
    CHARACTER(LEN=255)  :: units,    delta_t, begin_d
    CHARACTER(LEN=255)  :: begin_t,  incr

    ! "Fake" data arrays
    REAL*4              :: PS( ILONG, ILAT,        ITIME )  ! surface pressure
    REAL*4              :: T ( ILONG, ILAT, IVERT, ITIME )  ! temperature
    CHARACTER           :: DESC( ICHAR1, ICHAR2          )  ! Description
!
! !DEFINED PARAMETERS:
!
    LOGICAL, PARAMETER  :: COMPRESS = .TRUE.                ! Use compression

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

    ! Time dimension
    WRITE( 6, '(a)' ) 'Writing time  (dim     ) to netCDF file'
    CALL NcDef_Dimension( fId, 'time',     ITIME,  idTime )

    ! Altitude dimension
    WRITE( 6, '(a)' ) 'Writing lev   (dim     ) to netCDF file'
    CALL NcDef_Dimension( fId, 'lev',      IVERT,  idLev )

    ! Latitude dimension
    WRITE( 6, '(a)' ) 'Writing lat   (dim     ) to netCDF file'
    CALL NcDef_Dimension( fId, 'lat',      ILAT ,  idLat )

    ! Longitude dimension
    WRITE( 6, '(a)' ) 'Writing lon   (dim     ) to netCDF file'
    CALL NcDef_Dimension( fId, 'lon',      ILONG,  idLon )

    ! Character dimension 1
    WRITE( 6, '(a)' ) 'Writing cdim1 (dim     ) to netCDF file'
    CALL NcDef_Dimension( fId, 'cdim1', ICHAR1, idChar1 )

    ! Character dimension 1
    WRITE( 6, '(a)' ) 'Writing cdim2 (dim     ) to netCDF file'
    CALL NcDef_Dimension( fId, 'cdim2', ICHAR2, idChar2 )

    !=========================================================================
    ! Define the variables and variable attributes
    ! for COARDS compliance and GAMAP compliance
    !=========================================================================
    CALL NcDef_Glob_Attributes( fId, 'Title',       'NcdfUtilities test file' )
    CALL NcDef_Glob_Attributes( fId, 'History',     'test file - 24 Jan 2011' )
    CALL NcDef_Glob_Attributes( fId, 'Conventions', 'COARDS'                  )
    CALL NcDef_Glob_Attributes( fId, 'Model',       'GEOS4'                   )
    CALL NcDef_Glob_Attributes( fId, 'Nlayers',     '55'                      )
    CALL NcDef_Glob_Attributes( fId, 'Start_Date',  '20110101'                )
    CALL NcDef_Glob_Attributes( fId, 'Start_Time',  '00:00:00.0'              )
    CALL NcDef_Glob_Attributes( fId, 'End_Date',    '20110101'                )
    CALL NcDef_Glob_Attributes( fId, 'End_Time',    '23:59:59.0'              )
    CALL NcDef_Glob_Attributes( fId, 'Delta_Lon',   '5'                       )
    CALL NcDef_Glob_Attributes( fId, 'Delta_Lat',   '4'                       )
    CALL NcDef_Glob_Attributes( fId, 'Delta_time',  '000000'                  )
    CALL NcDef_Glob_Attributes( fId, 'Format',      'netCDF-3'                )
    CALL NcDef_Glob_Attributes( fId, 'valid_range', (/ -1e15, +1e15 /)        )
    CALL NcDef_Glob_Attributes( fId, 'id_number',   1                         )

    !=========================================================================
    ! Define the variables and variable attributes
    !=========================================================================

    ! Time index array (hardwire date to 2011/01/01)
    var1    = (/ idTime /)
    units   = 'minutes since 2011-01-01 00:00:00 GMT'
    delta_t = '0000-00-00 00:00:00'
    begin_d = '20110101'
    begin_t = '000000'
    incr    = '000000'
    CALL NcDef_Variable      ( fId, 'time', NF_INT,  1, var1, vId, COMPRESS )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name',      'time'           )
    CALL NcDef_Var_Attributes( fId, vId, 'units',          TRIM( units   )  )
    CALL NcDef_Var_Attributes( fId, vId, 'delta_t',        TRIM( delta_t )  )
    CALL NcDef_Var_Attributes( fId, vId, 'begin_date',     TRIM( begin_d )  )
    CALL NcDef_Var_Attributes( fId, vId, 'begin_time',     TRIM( begin_t )  )
    CALL NcDef_Var_Attributes( fId, vId, 'time_increment', TRIM( incr    )  )

    ! Define vertical (pressure) variable
    var1 = (/ idLev /)
    CALL NcDef_Variable( fId, 'lev', NF_DOUBLE, 1, var1, vId, COMPRESS )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'Pressure'       )
    CALL NcDef_Var_Attributes( fId, vId, 'units',     'hPa'            )

    ! Define latitude variable
    var1 = (/ idLat /)
    CALL NcDef_Variable( fId, 'lat', NF_DOUBLE, 1, var1, vId, COMPRESS )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name', 'Latitude'       )
    CALL NcDef_Var_Attributes( fId, vId, 'units',     'degrees_north'  )

    ! Define longitude variable
    var1 = (/ idLon /)
    CALL NcDef_Variable( fId, 'lon', NF_DOUBLE, 1, var1, vId, COMPRESS )
    CALL NcDef_Var_Attributes( fId, vId,  'long_name', 'Longitude'     )
    CALL NcDef_Var_Attributes( fId, vId,  'units',     'degrees_east'  )

    ! Define surface pressure variable
    var3 = (/ idLon, idLat, idTime /)
    CALL NcDef_Variable      ( fId, 'PS', NF_FLOAT, 3, var3, vId, COMPRESS    )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name',      'Surface Pressure' )
    CALL NcDef_Var_Attributes( fId, vId, 'units',          'hPa'              )
    CALL NcDef_Var_Attributes( fId, vId, 'gamap_category', 'GMAO-2D'          )
    CALL NcDef_Var_Attributes( fId, vId, 'missing_value',   1e15              )
    CALL NcDef_Var_Attributes( fId, vId, '_FillValue',      1e15              )
    CALL NcDef_Var_Attributes( fId, vId, 'valid_range',     (/-1e15, +1e15/)  )

    !=========================================================================
    ! %%% TEST RE-OPENING OF DEFINE MODE %%%
    !=========================================================================
    CALL NcEnd_Def( fId )
    WRITE( 6, '(a)' ) 'Testing re-opening of define mode'
    CALL NcBegin_Def( fId )

    ! Define temperature variable
    var4 = (/ idLon, idLat, idLev, idTime /)
    CALL NcDef_Variable      ( fId, 'T', NF_FLOAT, 4, var4, vId, COMPRESS     )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name',      'Temperature'      )
    CALL NcDef_Var_Attributes( fId, vId, 'units',          'K'                )
    CALL NcDef_Var_Attributes( fId, vId, 'gamap_category', 'GMAO-3D$'         )
    CALL NcDef_Var_Attributes( fId, vId, 'missing_value',   1e15              )
    CALL NcDef_Var_Attributes( fId, vId, '_FillValue',      1e15              )
    CALL NcDef_Var_Attributes( fId, vId, 'valid_range',     (/-1e15, +1e15/)  )


    ! Define description variable
    var2 = (/ idChar1, idChar2 /)
    CALL NcDef_Variable      ( fId, 'DESC', NF_CHAR, 2, var2, vId, COMPRESS   )
    CALL NcDef_Var_Attributes( fId, vId, 'long_name',      'Description'      )
    CALL NcDef_Var_Attributes( fId, vId, 'units',          '1'                )
    CALL NcDef_Var_Attributes( fId, vId, 'gamap_category', 'none'             )

    !=========================================================================
    ! %%% END OF DEFINITION SECTION %%%
    ! %%% NOW WRITE DATA TO FILE    %%%
    !=========================================================================
    CALL NcEnd_def( fId )

    ! Write longitude
    WRITE( 6, '(a)' ) 'Writing lon   (1D array) to netCDF file'
    st1d = (/ 1     /)
    ct1d = (/ ILONG /)
    CALL NcWr( longDat, fId, 'lon', st1d, ct1d )

    ! Write latitude
    WRITE( 6, '(a)' ) 'Writing lat   (1D array) to netCDF file'
    st1d = (/ 1    /)
    ct1d = (/ ILAT /)
    CALL NcWr( latDat, fId, 'lat', st1d, ct1d )

    ! Write pressure levels
    WRITE( 6, '(a)' ) 'Writing lev   (1D array) to netCDF file'
    st1d = (/ 1     /)
    ct1d = (/ IVERT /)
    CALL NcWr( levDat, fId, 'lev', st1d, ct1d )

    ! Write pressure levels
    WRITE( 6, '(a)' ) 'Writing time  (1D array) to netCDF file'
    st1d = (/ 1     /)
    ct1d = (/ ITIME /)
    CALL NcWr( timeDat, fId, 'time', st1d, ct1d )

    ! Write surface pressure (w/ fake values)
    WRITE( 6, '(a)' ) 'Writing PS    (3D array) to netCDF file'
    PS    = 1e0
    st3d = (/ 1,     1,   1      /)
    ct3d = (/ ILONG, ILAT, ITIME /)
    CALL NcWr( PS, fId, 'PS', st3d, ct3d )

    ! Write temperature (w/ fake values)
    WRITE( 6, '(a)' ) 'Writing T     (4D array) to netCDF file'
    T    = 1e0
    st4d = (/ 1,     1,    1,     1     /)
    ct4d = (/ ILONG, ILAT, IVERT, ITIME /)
    CALL NcWr( T, fId, 'T', st4d, ct4d )

    ! Initialzie the character array
    DO i = 1, ICHAR2
       DESC(1,i) = ACHAR(64+I)
       DESC(2,i) = ACHAR(96+I)
    ENDDO

    ! Write temperature (w/ fake values)
    WRITE( 6, '(a)' ) 'Writing DESC  (2D char ) to netCDF file'
    st2d = (/ 1,      1      /)
    ct2d = (/ ICHAR1, ICHAR2 /)
    CALL NcWr( DESC, fId, 'DESC', st2d, ct2d )

    !=========================================================================
    ! Close the netCDF file
    !=========================================================================
    CALL NcCl( fId )

    ! Echo info
    WRITE( 6, '(a)' ) '=== End netCDF file creation test ==='

  END SUBROUTINE TestNcdfCreate
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!  31 Jan 2012 - R. Yantosca - Bug fix in error checks for attributes
!  14 Jun 2012 - R. Yantosca - Now tests 2-D character read
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: fId,     rc,      XDim,    CDim1
    INTEGER                :: YDim,    ZDim,    TDim,    CDim2
    INTEGER                :: ct1d(1), ct2d(2), ct3d(3), ct4d(4)
    INTEGER                :: st1d(1), st2d(2), st3d(3), st4d(4)
    CHARACTER(LEN=255)     :: attValue
    REAL*4                 :: attValR4

    ! Arrays
    REAL*4                 :: valid(2)
    REAL*8,    ALLOCATABLE :: lon(:),  lat(:),  lev(:)
    INTEGER,   ALLOCATABLE :: time(:)
    REAL*4,    ALLOCATABLE :: PS(:,:,:)
    REAL*4,    ALLOCATABLE :: T(:,:,:,:)
    CHARACTER, ALLOCATABLE :: DESC(:,:)

    !=========================================================================
    ! Open the netCDF file
    !=========================================================================

    ! Echo info
    WRITE( 6, '(a)' ) '=== Begin netCDF file reading test ==='

    CALL Ncop_Rd( fId, 'my_filename.nc' )

    !=========================================================================
    ! Get the dimensions
    !=========================================================================
    CALL Ncget_Dimlen( fId, 'lon',   XDim  )
    CALL Ncget_Dimlen( fId, 'lat',   YDim  )
    CALL Ncget_Dimlen( fId, 'lev',   ZDim  )
    CALL Ncget_Dimlen( fId, 'time',  TDim  )
    CALL Ncget_Dimlen( fId, 'cdim1', CDim1 )
    CALL Ncget_Dimlen( fId, 'cdim2', CDim2 )

    rc = XDim - ILONG
    CALL Check( 'Reading lon   (dim  )  back from netCDF file', rc, pCt, tCt )

    rc = YDim - ILAT
    CALL Check( 'Reading lat   (dim  )  back from netCDF file', rc, pCt, tCt )

    rc = ZDim - IVERT
    CALL Check( 'Reading lev   (dim  )  back from netCDF file', rc, pCt, tCt )

    rc = TDim - ITIME
    CALL Check( 'Reading time  (dim  )  back from netCDF file', rc, pCt, tCt )

    rc = CDim1 - ICHAR1
    CALL Check( 'Reading cdim1 (dim  )  back from netCDF file', rc, pCt, tCt )

    rc = CDim2 - ICHAR2
    CALL Check( 'Reading cdim2 (dim  )  back from netCDF file', rc, pCt, tCt )

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
    CALL Check( 'Reading lon   (array)  back from netCDF file', rc, pCt, tCt )

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
    CALL Check( 'Reading lat   (array)  back from netCDF file', rc, pCt, tCt )

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
    CALL Check( 'Reading lev   (array)  back from netCDF file', rc, pCt, tCt )

    !=========================================================================
    ! Read the TIME variable
    !=========================================================================

    ! Read data
    ALLOCATE( time( TDim ) )
    st1d = (/ 1    /)
    ct1d = (/ TDim /)
    CALL NcRd( time, fId, 'time', st1d, ct1d )

    ! Equality test
    rc = SUM( time - timeDat )
    CALL Check( 'Reading time  (array)  back from netCDF file', rc, pCt, tCt )

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
    CALL Check( 'Reading PS             back from netCDF file', rc, pCt, tCt )

    ! Read units attribute
    CALL NcGet_Var_Attributes( fId, 'PS', 'units', attValue )
    IF ( TRIM( attValue ) == 'hPa' ) THEN
       rc = 0
    ELSE
       rc = -1
    ENDIF
    CALL Check( 'Reading PS:units       back from netCDF file', rc, pCt, tCt )

    ! Read long_name attribute
    CALL NcGet_Var_Attributes( fId, 'PS', 'long_name', attValue )
    IF ( TRIM( attValue ) == 'Surface Pressure' ) THEN
       rc = 0
    ELSE
       rc = -1
    ENDIF
    CALL Check( 'Reading PS:long_name   back from netCDF file', rc, pCt, tCt )

    ! Read _FillValue attribute
    CALL NcGet_Var_Attributes( fId, 'PS', '_FillValue', attValR4 )
    IF ( attValR4 == 1e15 ) THEN
       rc = 0
    ELSE
       rc = -1
    ENDIF
    CALL Check( 'Reading PS:_FillValue  back from netCDF file', rc, pCt, tCt )

    ! Read valid_range attribute
    CALL NcGet_Var_Attributes( fId, 'PS', 'valid_range', valid )
    IF ( valid(1) == -1e15 .and. valid(2) == 1e15 ) THEN
       rc = 0
    ELSE
       rc = -1
    ENDIF
    CALL Check( 'Reading PS:valid_range back from netCDF file', rc, pCt, tCt )

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
    CALL Check( 'Reading T              back from netCDF file', rc, pCt, tCt )

    ! Read units attribute
    CALL NcGet_Var_Attributes( fId, 'T', 'units', attValue )
    IF ( TRIM( attValue ) == 'K' ) THEN
       rc = 0
    ELSE
       rc = -1
    ENDIF
    CALL Check( 'Reading T:units        back from netCDF file', rc, pCt, tCt )

    ! Read long_name
    CALL NcGet_Var_Attributes( fId, 'T', 'long_name', attValue )
    IF ( TRIM( attValue ) == 'Temperature' ) THEN
       rc = 0
    ELSE
       rc = -1
    ENDIF
    CALL Check( 'Reading T:long_name    back from netCDF file', rc, pCt, tCt )

    ! Read _FillValue attribute
    CALL NcGet_Var_Attributes( fId, 'T', '_FillValue', attValR4 )
    IF ( attValR4 == 1e15 ) THEN
       rc = 0
    ELSE
       rc = -1
    ENDIF
    CALL Check( 'Reading T:_FillValue   back from netCDF file', rc, pCt, tCt )

    ! Read valid_range attribute
    CALL NcGet_Var_Attributes( fId, 'T', 'valid_range', valid )
    IF ( valid(1) == -1e15 .and. valid(2) == 1e15 ) THEN
       rc = 0
    ELSE
       rc = -1
    ENDIF
    CALL Check( 'Reading T:valid_range  back from netCDF file', rc, pCt, tCt )

    !=========================================================================
    ! Read the DESC variable
    !=========================================================================

    ! Read data
    ALLOCATE( DESC( CDim1, CDim2 ) )
    st2d = (/ 1,    1      /)
    ct2d = (/ CDim1, CDim2 /)
    CALL NcRd( DESC, fId, 'DESC', st2d, ct2d )

    ! Check that DESC was read properly
    rc = 0
    DO i = 1, ICHAR2
       IF ( ICHAR( DESC(1,i) ) - 64 /= I ) rc = 1
       IF ( ICHAR( DESC(2,i) ) - 96 /= I ) rc = 1
    ENDDO
    CALL Check( 'Reading DESC           back from netCDF file', rc, pCt, tCt )

    ! Read units attribute
    CALL NcGet_Var_Attributes( fId, 'DESC', 'units', attValue )
    IF ( TRIM( attValue ) == '1' ) THEN
       rc = 0
    ELSE
       rc = -1
    ENDIF
    CALL Check( 'Reading DESC:units     back from netCDF file', rc, pCt, tCt )

    ! Read long_name attribute
    CALL NcGet_Var_Attributes( fId, 'T', 'long_name', attValue )
    IF ( TRIM( attValue ) == 'Temperature' ) THEN
       rc = 0
    ELSE
       rc = -1
    ENDIF
    CALL Check( 'Reading DESC:long_name back from netCDF file', rc, pCt, tCt )

    !=========================================================================
    ! Read global attributes
    !=========================================================================

    ! Read title attribute
    CALL NcGet_Glob_Attributes( fId, 'Title', attValue )
    IF ( TRIM( attValue ) == 'NcdfUtilities test file' ) THEN
       rc = 0
    ELSE
       rc = -1
    ENDIF
    CALL Check( 'Reading title          back from netCDF file', rc, pCt, tCt )

    ! Read start_date
    CALL NcGet_Glob_Attributes( fId, 'Start_Date', attValue )
    IF ( TRIM( attValue ) == '20110101' ) THEN
       rc = 0
    ELSE
       rc = -1
    ENDIF
    CALL Check( 'Reading start_date     back from netCDF file', rc, pCt, tCt )

    ! Read start_time
    CALL NcGet_Glob_Attributes( fId, 'Start_Time', attValue )
    IF ( TRIM( attValue ) == '00:00:00.0' ) THEN
       rc = 0
    ELSE
       rc = -1
    ENDIF
    CALL Check( 'Reading start_time     back from netCDF file', rc, pCt, tCt )

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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
!  14 Jun 2012 - R. Yantosca - Now add 10 more . characters
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
       WRITE( 6, '(a)' ) msg // REPEAT( '.', 55-s ) // 'PASSED'
       passCt = passCt + 1
    ELSE
       WRITE( 6, '(a)' ) msg // REPEAT( '.', 55-s ) // 'FAILED'
    ENDIF

    totCt = totCt + 1

  END SUBROUTINE Check
!EOC

END PROGRAM TestNcdfUtil

