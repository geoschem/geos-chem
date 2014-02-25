!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: READSOIL
!
! !DESCRIPTION: Subroutine READSOIL reads in soiltype data, fertilizer data, 
!  and climate data.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE READSOIL
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
	! LAND2 and size
      USE CMN_SIZE_MOD   	! Size parameters
      USE COMMSOIL_MOD		! Soil variables

      IMPLICIT NONE
#     include "netcdf.inc"
      
! !REVISION HISTORY:
!  24 Aug 2012 - J.D. Maasakkers - Initial version
!  03 Jan 2013 - M. Payer        - Renamed PERAREA to IS_MASS in DO_REGRID_A2A
!  20 Aug 2013 - R. Yantosca     - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !Variables for NetCDF reading

      INTEGER :: fId    ! netCDF file ID

      ! Data arrays
      REAL*4             :: lon(1440)
      REAL*4             :: lat(720)
      INTEGER            :: koppen(24)
      INTEGER            :: time(1)
      INTEGER            :: LANDFRAC(1440,720,24,1)

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
      INTEGER            :: K  
       REAL*8 :: INGRID(1440,720)
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
     &          'soilNOx.landtype.generic.025x025.nc'
      CALL Ncop_Rd( fId, TRIM(nc_path) )

      ! Echo info to stdout
      nc_file = 'soilNOx.landtype.generic.025x025.nc'  
      nc_dir  = TRIM( DATA_DIR_1x1 ) // 'soil_NOx_201208/' 
      WRITE( 6, 100 ) REPEAT( '%', 79 )
      WRITE( 6, 110 ) TRIM(nc_file)
      WRITE( 6, 120 ) TRIM(nc_dir)

      !----------------------------------------
      ! VARIABLE: LANDFRAC
      !----------------------------------------

      ! Variable name
      v_name = "LANDFRAC"

      ! Read LANDFRAC from file
      st4d   = (/ 1, 1, 1, 1 /)
      ct4d   = (/ 1440, 720, 24, 1 /)
      CALL NcRd( LANDFRAC, fId, TRIM(v_name), st4d, ct4d )

      ! Read the LANDFRAC:long_name attribute
      a_name = "long_name"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Read the LANDFRAC:units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      ! Echo info to stdout
      WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)

      ! Read the LANDFRAC:gamap_category attribute
      a_name = "gamap_category"
      CALL NcGet_Var_Attributes( fId,TRIM(v_name),TRIM(a_name),a_val )

      !=================================================================
      ! Regrid
      !=================================================================

      ! 0.25x0.25 regrid - LAND2
      ! File with lat/lon edges for regridding
      LLFILENAME = TRIM( DATA_DIR_1x1) //
     &             'MAP_A2A_Regrid_201203/latlon_generic025x025.txt'

      ! Regrid to current model resolution [unitless]
      DO K = 1,24
         INGRID = LANDFRAC(:,:,K,1)
         OUTGRID = 0.

         CALL DO_REGRID_A2A( LLFILENAME, 1440,    720,
     &                       INGRID,     OUTGRID, IS_MASS=0 )

         LAND2(:,:,K) = OUTGRID
      ENDDO 

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
      END SUBROUTINE READSOIL
