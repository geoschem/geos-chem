!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
!
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_mod
!
! !DESCRIPTION: Module NCDF\_MOD contains routines to read data from
! netCDF files.
!\\
!\\
! !INTERFACE:
!
      MODULE NCDF_MOD
!
! !USES:
!
      ! Modules for netCDF read
      USE m_netcdf_io_open
      USE m_netcdf_io_get_dimlen
      USE m_netcdf_io_read
      USE m_netcdf_io_readattr
      USE m_netcdf_io_close
      USE m_netcdf_io_create
      USE m_netcdf_io_define
      USE m_netcdf_io_write

      IMPLICIT NONE
#     include "define.h"
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: NC_READ
      PUBLIC  :: GET_TIME_IDX

      INTERFACE NC_READ
         MODULE PROCEDURE NC_READ_3D
         MODULE PROCEDURE NC_READ_4D
      END INTERFACE
!
! !REVISION HISTORY:
!  27 Jul 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
!
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nc_read_3d
!
! !DESCRIPTION: Routine to read a 3D-field and corresponding variables rom a
!  netCDF file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_READ_3D( NC_PATH, PARA   , ARRAY,   YEAR     ,
     &                       MONTH  , DAY    , HOUR  ,  VERBOSE  ,
     &                       UNITS  , TEMPRES, GRID  ,  SHORTNAME,
     &                       ILON   , ILAT   , ITIME  , LON      ,
     &                       LAT    , TIME )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : ERROR_STOP

      IMPLICIT NONE
#     include "netcdf.inc"
!
! !ARGUMENTS:
!
      ! Required input
      CHARACTER(LEN=*),  INTENT(IN)   :: NC_PATH  ! Entire file path
      CHARACTER(LEN=*),  INTENT(IN)   :: PARA     ! Paramater to read
      REAL*4,            POINTER      :: ARRAY(:,:,:) ! 3D output array

      ! Optional input
      INTEGER,          INTENT(IN),  OPTIONAL  :: YEAR  ! desired year
      INTEGER,          INTENT(IN),  OPTIONAL  :: MONTH ! desired month
      INTEGER,          INTENT(IN),  OPTIONAL  :: DAY   ! desired day
      INTEGER,          INTENT(IN),  OPTIONAL  :: HOUR  ! desired hour
      LOGICAL,          INTENT(IN),  OPTIONAL  :: VERBOSE

      ! Optional output
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: SHORTNAME ! shortname
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: UNITS    ! unit
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: TEMPRES  ! temp. res.
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: GRID     ! grid type
      INTEGER,          INTENT(OUT), OPTIONAL  :: ILON     ! lon extent
      INTEGER,          INTENT(OUT), OPTIONAL  :: ILAT     ! lat extent
      INTEGER,          INTENT(OUT), OPTIONAL  :: ITIME    ! time extent
      REAL*4,           POINTER,     OPTIONAL  :: LON(:)   ! lon values
      REAL*4,           POINTER,     OPTIONAL  :: LAT(:)   ! lat values
      INTEGER,          POINTER,     OPTIONAL  :: TIME(:)  ! time vals
!
! !REMARKS:
!
! !REVISION HISTORY:
!  27 Jul 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

!=================================================================
! Variable declarations
!=================================================================

      ! Data arrays
      REAL*4,  ALLOCATABLE   :: longitude(:)
      REAL*4,  ALLOCATABLE   :: latitude(:)
      INTEGER, ALLOCATABLE   :: timevec(:)
      INTEGER                :: fId       ! netCDF file ID
      CHARACTER(LEN=255)     :: v_name    ! netCDF variable name
      CHARACTER(LEN=255)     :: a_name    ! netCDF attribute name
      CHARACTER(LEN=255)     :: a_val     ! netCDF attribute value

      ! Arrays for netCDF start and count values
      INTEGER                :: st1d(1), ct1d(1)   ! For 1D arrays
      INTEGER                :: st3d(3), ct3d(3)   ! For 3D arrays

      ! Grid dimensions
      INTEGER                :: II, JJ, LL

      ! Verbose
      LOGICAL                :: VERB

      ! For time index
      INTEGER                :: USE_YEAR
      INTEGER                :: USE_MONTH
      INTEGER                :: USE_DAY
      INTEGER                :: USE_HOUR
      INTEGER                :: TIDX

      ! For error handling
      CHARACTER(LEN=255)     :: LOC
      CHARACTER(LEN=255)     :: MSG
      CHARACTER(LEN=4)       :: YYYY
      CHARACTER(LEN=2)       :: MM, DD, HH


!=================================================================
! Open and read data from the netCDF file
!=================================================================

      ! Init
      IF ( PRESENT ( VERBOSE ) ) THEN
         VERB = VERBOSE
      ELSE
         VERB = .TRUE.
      ENDIF

      ! Check for temporal resolution
      IF ( PRESENT ( YEAR ) ) THEN
         USE_YEAR = YEAR
      ELSE
         USE_YEAR = -1
      ENDIF
      IF ( PRESENT ( MONTH ) ) THEN
         USE_MONTH = MONTH
      ELSE
         USE_MONTH = -1
      ENDIF
      IF ( PRESENT ( DAY ) ) THEN
         USE_DAY = DAY
      ELSE
         USE_DAY = -1
      ENDIF
      IF ( PRESENT ( HOUR ) ) THEN
         USE_HOUR = HOUR
      ELSE
         USE_HOUR = -1
      ENDIF

      ! For error handling
      LOC = 'NC_READ_3D ("ncdf_mod.F")'

      ! Open netCDF file
      CALL Ncop_Rd( fId, TRIM(NC_PATH) )

      ! Prompt some output (testing only)
      WRITE( *, * ) '####################################'
      WRITE( *, * ) 'Reading ', TRIM(NC_PATH)

      ! Echo info to stdout
      IF ( VERB ) THEN
         WRITE( 6, 100 ) REPEAT( '%', 79 )
         WRITE( 6, 110 ) TRIM(NC_PATH)
      ENDIF

      !----------------------------------------
      ! VARIABLE: lon
      !----------------------------------------

      ! Variable name
      v_name = "lon"

      ! Get longitudinal dimension and allocate longitude vector
      CALL Ncget_Dimlen ( fId, v_name, II )
      ALLOCATE ( longitude(II) )

      ! Read lon from file
      st1d   = (/ 1  /)
      ct1d   = (/ II /)
      CALL NcRd( longitude, fId, TRIM(v_name), st1d, ct1d )

      IF ( VERB ) THEN
         ! Read the lon:long_name attribute
         a_name = "long_name"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Read the lon:units attribute
         a_name = "units"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Echo info to stdout
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: lat
      !----------------------------------------

      ! Variable name
      v_name = "lat"

      ! Get longitudinal dimension and allocate longitude vector
      CALL Ncget_Dimlen ( fId, v_name, JJ )
      ALLOCATE ( latitude(JJ) )

      ! Read lat from file
      st1d   = (/ 1  /)
      ct1d   = (/ JJ /)
      CALL NcRd( latitude, fId, TRIM(v_name), st1d, ct1d )

      IF ( VERB ) THEN
         ! Read the lat:long_name attribute
         a_name = "long_name"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Read the lat:units attribute
         a_name = "units"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Echo info to stdout
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: time
      !----------------------------------------

      ! Variable name
      v_name = "time"

      ! Get longitudinal dimension and allocate longitude vector
      CALL Ncget_Dimlen ( fId, v_name, LL )
      ALLOCATE ( timevec(LL) )

      ! Read time from file
      st1d   = (/ 1  /)
      ct1d   = (/ LL /)
      CALL NcRd( timevec, fId, TRIM(v_name), st1d, ct1d )

      IF ( VERB ) THEN
         ! Read the time:long_name attribute
         a_name = "long_name"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Read the time:units attribute
         a_name = "units"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Echo info to stdout
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: PARA
      !----------------------------------------

      ! Check for time index
      ! If year is given:
      IF ( USE_YEAR > 0 ) THEN

         ! Make sure that month and day have valid entries
         USE_MONTH = MIN ( MAX( 1, USE_MONTH ), 12 )
         USE_DAY   = MIN ( MAX( 1, USE_DAY   ), 31 )
         USE_HOUR  = MIN ( MAX( 0, USE_HOUR  ), 23 )

         ! Get the time index
         CALL GET_TIME_IDX ( LL, timevec, USE_YEAR, USE_MONTH,
     &      USE_DAY, USE_HOUR, TIDX )

         ! Error if no agreement was found
         IF ( TIDX == 0 ) THEN
            WRITE ( YYYY, '(i4)' ) USE_YEAR
            WRITE ( MM,   '(i2)' ) USE_MONTH
            WRITE ( DD,   '(i2)' ) USE_DAY
            WRITE ( HH,   '(i2)' ) USE_HOUR
            MSG = 'Cannot find time stamp ' // YYYY // '/'     //
     &         MM // '/' // DD // '; hour ' // HH // ' in file ' //
     &         TRIM(NC_PATH) // ' - netCDF reading failed!'
            CALL ERROR_STOP ( MSG, LOC )
         ENDIF

         ! Set time dimension to 1
         LL = 1

      ! If only month is given (for monthly data)
      ELSEIF ( USE_MONTH > 0 ) THEN

         ! Check if it's indeed monthly data:
         IF ( LL /= 12 ) THEN
            MSG = 'file ' // TRIM(NC_PATH) // ' has not monthly ' //
     &         'data, cannot read month!'
            CALL ERROR_STOP ( MSG, LOC )
         ENDIF

         ! Set time index to specified month
         TIDX = USE_MONTH

         ! Set time dimension to 1
         LL = 1

      ! Otherwise, read all temporal indices, start w/ index 1
      ELSE
         TIDX = 1

      ENDIF

      ! Variable name
      v_name = TRIM(PARA)

      ! Prompt some output (testing only)
      WRITE( *, * ) 'Reading para ', TRIM(PARA)
      WRITE( *, * ) '3D Array dimensions: ', II, ', ', JJ, ', ', LL
      WRITE( *, * ) 'TIDX ', TIDX

      ! Allocate the array
      ALLOCATE ( ARRAY(II,JJ,LL) )

      ! Read the 3D - field from file
      st3d   = (/ 1 , 1 , TIDX /)
      ct3d   = (/ II, JJ, LL   /)
      CALL NcRd( ARRAY, fId, TRIM(v_name), st3d, ct3d )

      ! Read some attributes for printing only
      IF ( VERB ) THEN
         ! Read the long_name attribute
         a_name = "long_name"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Read the units attribute
         a_name = "units"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)
         IF ( PRESENT ( UNITS ) ) UNITS = a_val

         ! Echo info to stdout
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)

         ! Read the gamap_category attribute
         a_name = "gamap_category"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)
      ENDIF

      ! ----------------------------
      ! Optional attributes
      ! ----------------------------

      IF ( PRESENT ( UNITS ) .AND. .NOT. VERB ) THEN
         a_name = "units"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)
         UNITS = TRIM(a_val)
      ENDIF

      ! Read the short_name attribute
      IF ( PRESENT ( SHORTNAME ) ) THEN
         a_name = "short_name"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)
         SHORTNAME = TRIM(a_val)
      ENDIF


!=================================================================
! Global attributes
!=================================================================
      IF ( PRESENT ( GRID ) ) THEN
         a_name = "Grid"
         CALL NcGet_Glob_Attributes(fId,TRIM(a_name),a_val)
         GRID = TRIM(a_val)
      ENDIF

      ! Read the temp_res attribute
      IF ( PRESENT ( TEMPRES ) ) THEN
         a_name = "Temp_Res"
         CALL NcGet_Glob_Attributes(fId,TRIM(a_name),a_val)
         TEMPRES = TRIM(a_val)
      ENDIF


!=================================================================
! Field dimensions
!=================================================================
      IF ( PRESENT ( ILON  ) ) ILON  = II
      IF ( PRESENT ( ILAT  ) ) ILAT  = JJ
      IF ( PRESENT ( ITIME ) ) ITIME = LL
      IF ( PRESENT ( LON   ) ) THEN
         ALLOCATE( LON(II) )
         LON   = longitude
      ENDIF
      IF ( PRESENT ( LAT   ) ) THEN
         ALLOCATE( LAT(JJ) )
         LAT   = latitude
      ENDIF
      IF ( PRESENT ( TIME  ) ) THEN
         ALLOCATE( TIME(LL) )
         IF ( LL == 1 ) THEN
            TIME = timevec(TIDX)
         ELSE
            TIME  = timevec
         ENDIF
      ENDIF


!=================================================================
! Cleanup and quit
!=================================================================

      ! Deallocate internal variables
      IF ( ALLOCATED ( longitude ) ) DEALLOCATE ( longitude )
      IF ( ALLOCATED ( latitude  ) ) DEALLOCATE ( latitude  )
      IF ( ALLOCATED ( timevec   ) ) DEALLOCATE ( timevec   )

      ! Close netCDF file
      CALL NcCl( fId )

      ! Prompt some output (testing only)
      WRITE( *, * ) '####################################'


      ! Echo info to stdout
      IF ( VERB ) THEN
         WRITE( 6, 140 )
         WRITE( 6, 100 ) REPEAT( '%', 79 )
      ENDIF

      ! FORMAT statements
 100  FORMAT( a                                              )
 110  FORMAT( '%% Opening file  : ',         a               )
 130  FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
 140  FORMAT( '%% Successfully closed file!'                 )

      END SUBROUTINE NC_READ_3D
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
!
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nc_read_4d
!
! !DESCRIPTION: Routine to read a 4D-field and corresponding variables from a
!  netCDF file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_READ_4D( NC_PATH, PARA   , ARRAY, YEAR     ,
     &                       MONTH  , DAY    , HOUR , VERBOSE  ,
     &                       UNITS  , TEMPRES, GRID , SHORTNAME,
     &                       ILON   , ILAT   , ILEV , ITIME    ,
     &                       LON    , LAT    , LEV  , TIME      )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : ERROR_STOP

      IMPLICIT NONE
#     include "netcdf.inc"
!
! !ARGUMENTS:
!
      ! Required
      CHARACTER(LEN=*), INTENT(IN)   :: NC_PATH  ! Entire file path
      CHARACTER(LEN=*), INTENT(IN)   :: PARA     ! Paramater to read
      REAL*4,           POINTER      :: ARRAY(:,:,:,:) ! 4D output array

      ! Optional
      INTEGER,          INTENT(IN) , OPTIONAL  :: YEAR
      INTEGER,          INTENT(IN) , OPTIONAL  :: MONTH
      INTEGER,          INTENT(IN) , OPTIONAL  :: DAY
      INTEGER,          INTENT(IN) , OPTIONAL  :: HOUR
      LOGICAL,          INTENT(IN) , OPTIONAL  :: VERBOSE

      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: SHORTNAME  ! Shortname
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: UNITS     ! Unit
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: TEMPRES   ! Temp. Res.
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL  :: GRID      ! Grid type
      INTEGER,          INTENT(OUT), OPTIONAL  :: ILON      ! lon ext.
      INTEGER,          INTENT(OUT), OPTIONAL  :: ILAT      ! lat ext.
      INTEGER,          INTENT(OUT), OPTIONAL  :: ILEV      ! vert. ext.
      INTEGER,          INTENT(OUT), OPTIONAL  :: ITIME     ! temp. ext.
      REAL*4,           POINTER,     OPTIONAL  :: LON(:)
      REAL*4,           POINTER,     OPTIONAL  :: LAT(:)
      REAL*4,           POINTER,     OPTIONAL  :: LEV(:)
      INTEGER,          POINTER,     OPTIONAL  :: TIME(:)
!
! !REMARKS:
!
! !REVISION HISTORY:
!  27 Jul 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

!=================================================================
! Variable declarations
!=================================================================

      ! Data arrays
      REAL*4,  ALLOCATABLE   :: longitude(:)
      REAL*4,  ALLOCATABLE   :: latitude(:)
      REAL*4,  ALLOCATABLE   :: levels(:)
      INTEGER, ALLOCATABLE   :: timevec(:)
      INTEGER                :: fId       ! netCDF file ID
      CHARACTER(LEN=255)     :: v_name    ! netCDF variable name
      CHARACTER(LEN=255)     :: a_name    ! netCDF attribute name
      CHARACTER(LEN=255)     :: a_val     ! netCDF attribute value

      ! Arrays for netCDF start and count values
      INTEGER                :: st1d(1), ct1d(1)   ! For 1D arrays
      INTEGER                :: st4d(4), ct4d(4)   ! For 3D arrays

      ! Grid dimensions
      INTEGER                :: II, JJ, KK, LL

      ! For time index
      INTEGER                :: USE_YEAR
      INTEGER                :: USE_MONTH
      INTEGER                :: USE_DAY
      INTEGER                :: USE_HOUR
      INTEGER                :: TIDX

      ! Verbose
      LOGICAL                :: VERB

      ! For error handling
      CHARACTER(LEN=255)     :: LOC
      CHARACTER(LEN=255)     :: MSG
      CHARACTER(LEN=4)       :: YYYY
      CHARACTER(LEN=2)       :: MM, DD, HH


!=================================================================
! NC_READ_4D begins here
!=================================================================

      ! Init
      IF ( PRESENT ( VERBOSE ) ) THEN
         VERB = VERBOSE
      ELSE
         VERB = .TRUE.
      ENDIF

      ! Check for temporal resolution
      IF ( PRESENT ( YEAR ) ) THEN
         USE_YEAR = YEAR
      ELSE
         USE_YEAR = -1
      ENDIF
      IF ( PRESENT ( MONTH ) ) THEN
         USE_MONTH = MONTH
      ELSE
         USE_MONTH = -1
      ENDIF
      IF ( PRESENT ( DAY ) ) THEN
         USE_DAY = DAY
      ELSE
         USE_DAY = -1
      ENDIF
      IF ( PRESENT ( HOUR ) ) THEN
         USE_HOUR = HOUR
      ELSE
         USE_HOUR = -1
      ENDIF

      ! For error handling
      LOC = 'NC_READ_4D ("ncdf_mod.F")'

      ! Open netCDF file
      CALL Ncop_Rd( fId, TRIM(NC_PATH) )

      IF ( VERB ) THEN
         ! Echo info to stdout
         WRITE( 6, 100 ) REPEAT( '%', 79 )
         WRITE( 6, 110 ) TRIM(NC_PATH)
      ENDIF

      !----------------------------------------
      ! VARIABLE: lon
      !----------------------------------------

      ! Variable name
      v_name = "lon"

      ! Get longitudinal dimension and allocate longitude vector
      CALL Ncget_Dimlen ( fId, v_name, II )
      ALLOCATE ( longitude(II) )

      ! Read lon from file
      st1d   = (/ 1 /)
      ct1d   = (/ II /)
      CALL NcRd( longitude, fId, TRIM(v_name), st1d, ct1d )

      IF ( VERB ) THEN
         ! Read the lon:long_name attribute
         a_name = "long_name"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Read the lon:units attribute
         a_name = "units"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Echo info to stdout
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: lat
      !----------------------------------------

      ! Variable name
      v_name = "lat"

      ! Get longitudinal dimension and allocate longitude vector
      CALL Ncget_Dimlen ( fId, v_name, JJ )
      ALLOCATE ( latitude(JJ) )

      ! Read lat from file
      st1d   = (/ 1 /)
      ct1d   = (/ JJ /)
      CALL NcRd( latitude, fId, TRIM(v_name), st1d, ct1d )

      IF ( VERB ) THEN
         ! Read the lat:long_name attribute
         a_name = "long_name"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Read the lat:units attribute
         a_name = "units"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Echo info to stdout
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: lev
      !----------------------------------------

      ! Variable name
      v_name = "lev"

      ! Get longitudinal dimension and allocate longitude vector
      CALL Ncget_Dimlen ( fId, v_name, KK )
      ALLOCATE ( levels(KK) )

      ! Read lon from file
      st1d   = (/ 1 /)
      ct1d   = (/ KK /)
      CALL NcRd( levels, fId, TRIM(v_name), st1d, ct1d )

      IF ( VERB ) THEN
         ! Read the lon:long_name attribute
         a_name = "long_name"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Read the lon:units attribute
         a_name = "units"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Echo info to stdout
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: time
      !----------------------------------------

      ! Variable name
      v_name = "time"

      ! Get longitudinal dimension and allocate longitude vector
      CALL Ncget_Dimlen ( fId, v_name, LL )
      ALLOCATE ( timevec(LL) )

      ! Read time from file
      st1d   = (/ 1 /)
      ct1d   = (/ LL /)
      CALL NcRd( timevec, fId, TRIM(v_name), st1d, ct1d )

      IF ( VERB ) THEN
         ! Read the time:long_name attribute
         a_name = "long_name"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Read the time:units attribute
         a_name = "units"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Echo info to stdout
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)
      ENDIF

      !----------------------------------------
      ! VARIABLE: PARA
      !----------------------------------------

      ! Check for time index
      ! If year is given:
      IF ( USE_YEAR > 0 ) THEN

         ! Make sure that month and day have valid entries
         USE_MONTH = MIN ( MAX( 1, USE_MONTH ), 12 )
         USE_DAY   = MIN ( MAX( 1, USE_DAY   ), 31 )
         USE_HOUR  = MIN ( MAX(01, USE_HOUR  ), 23 )

         ! Get the time index
         CALL GET_TIME_IDX ( LL, timevec, USE_YEAR, USE_MONTH,
     &      USE_DAY, USE_HOUR, TIDX )

         ! Error if no agreement was found
         IF ( TIDX == 0 ) THEN
            WRITE ( YYYY, '(i4)' ) USE_YEAR
            WRITE ( MM,   '(i2)' ) USE_MONTH
            WRITE ( DD,   '(i2)' ) USE_DAY
            WRITE ( HH,   '(i2)' ) USE_HOUR
            MSG = 'Cannot find time stamp ' // YYYY // '/'     //
     &         MM // '/' // DD // '; hour ' // HH // ' in file ' //
     &         TRIM(NC_PATH) // ' - netCDF reading failed!'
            CALL ERROR_STOP ( MSG, LOC )
         ENDIF

         ! Set time dimension to 1
         LL = 1

      ! If only month is given (for monthly data)
      ELSEIF ( USE_MONTH > 0 ) THEN

         ! Check if it's indeed monthly data:
         IF ( LL /= 12 ) THEN
            MSG = 'file ' // TRIM(NC_PATH) // ' has not monthly ' //
     &         'data, cannot read month!'
            CALL ERROR_STOP ( MSG, LOC )
         ENDIF

         ! Set time index to specified month
         TIDX = USE_MONTH

         ! Set time dimension to 1
         LL = 1

      ! Otherwise, read all temporal indices, start w/ index 1
      ELSE
         TIDX = 1

      ENDIF

      ! Variable name
      v_name = TRIM(PARA)

      ! Allocate the array
      ALLOCATE ( ARRAY(II,JJ,KK,LL) )

      ! Read the 3D - field from file
      st4d   = (/ 1,  1 , 1 , TIDX /)
      ct4d   = (/ II, JJ, KK, LL   /)
      CALL NcRd( ARRAY, fId, TRIM(v_name), st4d, ct4d )

      ! Read some attributes for printing only
      IF ( VERB ) THEN
         ! Read the long_name attribute
         a_name = "long_name"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)

         ! Read the units attribute
         a_name = "units"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)
         IF ( PRESENT ( UNITS ) ) UNITS = TRIM(a_val)

         ! Echo info to stdout
         WRITE( 6, 130 ) TRIM(v_name), TRIM(a_val)

         ! Read the gamap_category attribute
         a_name = "gamap_category"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)
      ENDIF

      ! ----------------------------
      ! Optional arguments
      ! ----------------------------

      IF ( PRESENT ( UNITS ) .AND. .NOT. VERB ) THEN
         a_name = "units"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)
         UNITS = TRIM(a_val)
      ENDIF

      ! Read the short_name attribute
      IF ( PRESENT ( SHORTNAME ) ) THEN
         a_name = "short_name"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)
         SHORTNAME = TRIM(a_val)
      ENDIF


!=================================================================
! Global attributes
!=================================================================
      IF ( PRESENT ( GRID ) ) THEN
         a_name = "Grid"
         CALL NcGet_Glob_Attributes(fId,TRIM(a_name),a_val)
         GRID = TRIM(a_val)
      ENDIF

      ! Read the temp_res attribute
      IF ( PRESENT ( TEMPRES ) ) THEN
         a_name = "Temp_Res"
         CALL NcGet_Glob_Attributes(fId,TRIM(a_name),a_val)
         TEMPRES = TRIM(a_val)
      ENDIF


!=================================================================
! Field dimensions
!=================================================================
      IF ( PRESENT ( ILON  ) ) ILON  = II
      IF ( PRESENT ( ILAT  ) ) ILAT  = JJ
      IF ( PRESENT ( ILEV  ) ) ILEV  = KK
      IF ( PRESENT ( ITIME ) ) ITIME = LL
      IF ( PRESENT ( LON   ) ) THEN
         ALLOCATE( LON(II) )
         LON   = longitude
      ENDIF
      IF ( PRESENT ( LAT   ) ) THEN
         ALLOCATE( LAT(JJ) )
         LAT   = latitude
      ENDIF
      IF ( PRESENT ( LEV   ) ) THEN
         ALLOCATE( LEV(KK) )
         LEV   = levels
      ENDIF
      IF ( PRESENT ( TIME  ) ) THEN
         ALLOCATE( TIME(LL) )
        IF ( LL == 1 ) THEN
            TIME = timevec(TIDX)
         ELSE
            TIME  = timevec
         ENDIF
      ENDIF


!=================================================================
! Cleanup and quit
!=================================================================

      ! Deallocate internal variables
      IF ( ALLOCATED ( longitude ) ) DEALLOCATE ( longitude )
      IF ( ALLOCATED ( latitude  ) ) DEALLOCATE ( latitude  )
      IF ( ALLOCATED ( levels    ) ) DEALLOCATE ( levels    )
      IF ( ALLOCATED ( timevec   ) ) DEALLOCATE ( timevec   )

      ! Close netCDF file
      CALL NcCl( fId )

      IF ( VERB ) THEN
         ! Echo info to stdout
         WRITE( 6, 140 )
         WRITE( 6, 100 ) REPEAT( '%', 79 )

         ! FORMAT statements
 100     FORMAT( a                                              )
 110     FORMAT( '%% Opening file  : ',         a               )
 130     FORMAT( '%% Successfully read ',       a, ' [', a, ']' )
 140     FORMAT( '%% Successfully closed file!'                 )
      ENDIF

      END SUBROUTINE NC_READ_4D
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group
!
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_time_idx
!
! !DESCRIPTION: Routine which returns the index with the specified time
! for a given year vector. All values assumed to be relative to
! 1985-01-01.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_TIME_IDX( ITIME, TIME, YEAR, MONTH, DAY, HOUR, IDX)
!
! !USES:
!
      USE BPCH2_MOD, ONLY : GET_TAU0
!
! !ARGUMENTS:
!
      ! Required
      INTEGER, INTENT( IN)           :: ITIME
      INTEGER, INTENT( IN)           :: TIME(ITIME)
      INTEGER, INTENT( IN)           :: YEAR
      INTEGER, INTENT( IN)           :: MONTH
      INTEGER, INTENT( IN)           :: DAY
      INTEGER, INTENT( IN)           :: HOUR
      INTEGER, INTENT(OUT)           :: IDX
!
! !REMARKS:
!
! !REVISION HISTORY:
!  04 Nov 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: II
      REAL*8                 :: TAU


!=================================================================
! GET_TIME_IDX starts here
!=================================================================

      ! Initialize output arguments
      IDX = 0

      ! Get time since 1985.
      TAU = GET_TAU0(MONTH,DAY,YEAR,HOUR)

      ! Loop over all time stamps and check which one is the one that
      ! corresponds to the specified year. The two time stamps have to
      ! match exactly.
      DO II = 1, ITIME
         IF ( ( TIME(II) - INT(TAU) ) == 0 ) THEN
            IDX = II
            EXIT
         ENDIF
      ENDDO

      END SUBROUTINE GET_TIME_IDX
!EOC
      END MODULE NCDF_MOD
!EOM
