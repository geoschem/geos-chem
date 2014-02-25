!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: READCLIM
!
! !DESCRIPTION: Subroutine READCLIM reads in soiltype data, fertilizer data, 
!  and climate data.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READCLIM
!
! !USES:
!
      USE m_netcdf_io_open
      USE m_netcdf_io_get_dimlen
      USE m_netcdf_io_read
      USE m_netcdf_io_readattr
      USE m_netcdf_io_close
	! Regrid
      USE REGRID_A2A_MOD, ONLY : DO_REGRID_A2A
      USE DIRECTORY_MOD,  ONLY : DATA_DIR_1x1
	! CLIM and size
      USE CMN_SIZE_MOD   	! Size parameters
      USE COMMSOIL_MOD		! Soil variables

      IMPLICIT NONE
#     include "netcdf.inc"
      
! !REVISION HISTORY:
!  24 Aug 2012 - J.D. Maasakkers - Initial version
!  03 Jan 2013 - M. Payer        - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !Variables for NetCDF reading

      INTEGER :: fId    ! netCDF file ID

      ! Data arrays
      REAL*4             :: lon(720)
      REAL*4             :: lat(360)
      INTEGER            :: time(1)
      INTEGER            :: NON_ARID(720,360,1)
      INTEGER            :: ARID(720,360,1)

      ! Character strings
      CHARACTER(LEN=255) :: nc_dir             ! netCDF directory name
      CHARACTER(LEN=255) :: nc_file            ! netCDF file name
      CHARACTER(LEN=255) :: nc_path            ! netCDF path name
      CHARACTER(LEN=255) :: v_name             ! netCDF variable name 
      CHARACTER(LEN=255) :: a_name             ! netCDF attribute name
      CHARACTER(LEN=255) :: a_val              ! netCDF attribute value
      CHARACTER(LEN=255) :: LLFILENAME

      ! Arrays for netCDF start and count values
      INTEGER            :: st1d(1), ct1d(1)   ! For 1D arrays    
      INTEGER            :: st2d(2), ct2d(2)   ! For 2D arrays 
      INTEGER            :: st3d(3), ct3d(3)   ! For 3D arrays 
      INTEGER            :: st4d(4), ct4d(4)   ! For 4D arrays 
      INTEGER            :: st5d(5), ct5d(5)   ! For 5D arrays 
      INTEGER            :: st6d(6), ct6d(6)   ! For 6D arrays 
 
       ! Variables for regridding
       REAL*8 :: INGRID(720,360)
       REAL*8 :: OUTGRID(IIPAR,JJPAR)
      LOGICAL, SAVE      :: FIRST   = .TRUE.

      ! First-time only initialization
      IF ( FIRST ) THEN
         
         ! Reset First-time flag
         FIRST = .FALSE.

      !=================================================================
      ! Open and read data from the netCDF file
      !=================================================================

      ! Open netCDF file
      nc_path = TRIM( DATA_DIR_1x1 ) // 'soil_NOx_201208/' //
     &          'soilNOx.climate.generic.05x05.nc'
      CALL Ncop_Rd( fId, TRIM(nc_path) )

      ! Echo info to stdout
      nc_file = 'soilNOx.climate.generic.05x05.nc'  
      nc_dir  = TRIM( DATA_DIR_1x1 ) // 'soil_NOx_201208/' 
      WRITE( 6, 100 ) REPEAT( '%', 79 )
      WRITE( 6, 110 ) TRIM(nc_file)
      WRITE( 6, 120 ) TRIM(nc_dir)

      !----------------------------------------
      ! VARIABLE: NON_ARID
      !----------------------------------------

      ! Variable name
      v_name = "NON_ARID"

      ! Read NON_ARID from file
      st3d   = (/ 1, 1, 1 /)
      ct3d   = (/ 720, 360, 1 /)
      CALL NcRd( NON_ARID, fId, TRIM(v_name), st3d, ct3d )

      ! Read the NON_ARID:long_name attribute
      a_name = "long_name"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Read the NON_ARID:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)

      ! Read the NON_ARID:gamap_category attribute
      a_name = "gamap_category"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      !----------------------------------------
      ! VARIABLE: ARID
      !----------------------------------------

      ! Variable name
      v_name = "ARID"

      ! Read ARID from file
      st3d   = (/ 1, 1, 1 /)
      ct3d   = (/ 720, 360, 1 /)
      CALL NcRd( ARID, fId, TRIM(v_name), st3d, ct3d )

      ! Read the ARID:long_name attribute
      a_name = "long_name"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Read the ARID:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)

      ! Read the ARID:gamap_category attribute
      a_name = "gamap_category"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      !=================================================================
      ! Regrid
      !=================================================================

      ! 0.5x0.5 regrid - CLIM
      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/latlon_generic05x05.txt'

      ! Regrid non-arid layer [unitless]
      INGRID = NON_ARID(:,:,1)
      OUTGRID = 0.

      CALL DO_REGRID_A2A( LLFILENAME, 720,     360,
     &                    INGRID,     OUTGRID, IS_MASS=0 )

      CLIM(:,:,1) = OUTGRID 

      ! Regrid arid layer [unitless]
      INGRID = ARID(:,:,1)
      OUTGRID = 0.

      CALL DO_REGRID_A2A( LLFILENAME, 720,     360,
     &                    INGRID,     OUTGRID, IS_MASS=0 )

      CLIM(:,:,2) = OUTGRID 

      !=================================================================
      ! Cleanup and quit
      !=================================================================

      ! Close netCDF file
      CALL NcCl( fId )

      ! Echo info to stdout
      WRITE( 6, 140 )
      WRITE( 6, 100 ) REPEAT( '%', 79 )

      ! FORMAT statements
 100  FORMAT( a                                              )
 110  FORMAT( '%% Opening file  : ',         a               )
 120  FORMAT( '%%  in directory : ',         a, / , '%%'     )
 130  FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
 140  FORMAT( '%% Successfully closed file!'                 )

      ENDIF
!EOC
      END SUBROUTINE READCLIM
