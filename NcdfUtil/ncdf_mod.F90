!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GFSC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ncdf_mod.F90
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
  USE m_netcdf_io_checks

  IMPLICIT NONE
  PRIVATE
# include "netcdf.inc"
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: NC_OPEN
  PUBLIC  :: NC_CREATE
  PUBLIC  :: NC_SET_DEFMODE
  PUBLIC  :: NC_VAR_DEF
  PUBLIC  :: NC_VAR_CHUNK
  PUBLIC  :: NC_VAR_WRITE
  PUBLIC  :: NC_CLOSE
  PUBLIC  :: NC_READ_TIME
  PUBLIC  :: NC_READ_TIME_YYYYMMDDhhmm
  PUBLIC  :: NC_READ_VAR
  PUBLIC  :: NC_READ_ARR
  PUBLIC  :: NC_GET_REFDATETIME
  PUBLIC  :: NC_GET_GRID_EDGES
  PUBLIC  :: NC_GET_SIGMA_LEVELS
  PUBLIC  :: NC_WRITE
  PUBLIC  :: NC_ISMODELLEVEL
  PUBLIC  :: GET_TAU0
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: GET_TIDX
  PRIVATE :: TIMEUNIT_CHECK
  PRIVATE :: NC_WRITE_3D
  PRIVATE :: NC_WRITE_4D
  PRIVATE :: NC_VAR_WRITE_INT_1D
  PRIVATE :: NC_VAR_WRITE_INT_2D
  PRIVATE :: NC_VAR_WRITE_INT_3D
  PRIVATE :: NC_VAR_WRITE_INT_4D
  PRIVATE :: NC_VAR_WRITE_R4_1D
  PRIVATE :: NC_VAR_WRITE_R4_2D
  PRIVATE :: NC_VAR_WRITE_R4_3D
  PRIVATE :: NC_VAR_WRITE_R4_4D
  PRIVATE :: NC_VAR_WRITE_R8_0D
  PRIVATE :: NC_VAR_WRITE_R8_1D
  PRIVATE :: NC_VAR_WRITE_R8_2D
  PRIVATE :: NC_VAR_WRITE_R8_3D
  PRIVATE :: NC_VAR_WRITE_R8_4D
  PRIVATE :: NC_READ_VAR_SP
  PRIVATE :: NC_READ_VAR_DP
  PRIVATE :: NC_GET_GRID_EDGES_SP
  PRIVATE :: NC_GET_GRID_EDGES_DP
  PRIVATE :: NC_GET_GRID_EDGES_C
  PRIVATE :: NC_GET_SIGMA_LEVELS_SP
  PRIVATE :: NC_GET_SIGMA_LEVELS_DP
  PRIVATE :: NC_GET_SIGMA_LEVELS_C
  PRIVATE :: NC_GET_SIG_FROM_HYBRID
  PRIVATE :: NC_READ_VAR_CORE
!
! !REVISION HISTORY:
!  27 Jul 2012 - C. Keller   - Initial version
!  13 Jun 2014 - R. Yantosca - Now use F90 free-format indentation
!  13 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  10 Jul 2014 - R. Yantosca - Add GET_TAU0 as a PRIVATE local routine
!  12 Dec 2014 - C. Keller   - Added NC_ISMODELLEVEL
!  19 Sep 2016 - R. Yantosca - Rewrite NC_VAR_WRITE overloaded functions to
!                              remove optional args (which chokes Gfortran)
!  19 Sep 2016 - R. Yantosca - Now include netcdf.inc once at top of module
!  19 Sep 2016 - R. Yantosca - Remove extra IMPLICIT NONE statements, we only
!                              need to declare it once at the top of module
!  10 Apr 2017 - R. Yantosca - Renamed routine NC_READ_TIME_YYYYMMDDhh to
!                              NC_READ_TIME_YYYYMMDDhhmm, to indicate that
!                              it will now uses YYYYYMMDDhhmm format
!  09 Aug 2017 - R. Yantosca - Add public routine NC_SET_DEFMODE
!  25 Aug 2017 - R. Yantosca - Add NC_Var_Write_*_0D routines
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES:
!
  INTERFACE NC_WRITE
     MODULE PROCEDURE NC_WRITE_3D
     MODULE PROCEDURE NC_WRITE_4D
  END INTERFACE NC_WRITE

  INTERFACE NC_READ_VAR
     MODULE PROCEDURE NC_READ_VAR_SP
     MODULE PROCEDURE NC_READ_VAR_DP
  END INTERFACE NC_READ_VAR

  INTERFACE NC_GET_GRID_EDGES
     MODULE PROCEDURE NC_GET_GRID_EDGES_SP
     MODULE PROCEDURE NC_GET_GRID_EDGES_DP
  END INTERFACE NC_GET_GRID_EDGES

  INTERFACE NC_GET_SIGMA_LEVELS
     MODULE PROCEDURE NC_GET_SIGMA_LEVELS_SP
     MODULE PROCEDURE NC_GET_SIGMA_LEVELS_DP
  END INTERFACE NC_GET_SIGMA_LEVELS

  INTERFACE NC_VAR_WRITE
     MODULE PROCEDURE NC_VAR_WRITE_INT_0D
     MODULE PROCEDURE NC_VAR_WRITE_INT_1D
     MODULE PROCEDURE NC_VAR_WRITE_INT_2D
     MODULE PROCEDURE NC_VAR_WRITE_INT_3D
     MODULE PROCEDURE NC_VAR_WRITE_INT_4D
     MODULE PROCEDURE NC_VAR_WRITE_R4_0D
     MODULE PROCEDURE NC_VAR_WRITE_R4_1D
     MODULE PROCEDURE NC_VAR_WRITE_R4_2D
     MODULE PROCEDURE NC_VAR_WRITE_R4_3D
     MODULE PROCEDURE NC_VAR_WRITE_R4_4D
     MODULE PROCEDURE NC_VAR_WRITE_R8_0D
     MODULE PROCEDURE NC_VAR_WRITE_R8_1D
     MODULE PROCEDURE NC_VAR_WRITE_R8_2D
     MODULE PROCEDURE NC_VAR_WRITE_R8_3D
     MODULE PROCEDURE NC_VAR_WRITE_R8_4D
  END INTERFACE NC_VAR_WRITE

CONTAINS
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Open
!
! !DESCRIPTION: Simple wrapper routine to open the given netCDF file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_OPEN( FileName, fID )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   )  :: FileName
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)  :: fID
!
! !REVISION HISTORY:
!  04 Nov 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! NC_OPEN begins here
    !=================================================================

    ! Open netCDF file
    CALL Ncop_Rd( fId, TRIM(FileName) )


  END SUBROUTINE NC_OPEN
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Close
!
! !DESCRIPTION: Simple wrapper routine to close the given lun.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_CLOSE( fID )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN   ) :: fID
!
! !REVISION HISTORY:
!  04 Nov 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! NC_CLOSE begins here
    !=================================================================

    CALL NcCl( fID )

  END SUBROUTINE NC_CLOSE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Set_DefMode
!
! !DESCRIPTION: Toggles netCDF define mode on or off.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Nc_Set_DefMode( fId, On, Off )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: fId   ! netCDF file ID
    LOGICAL, OPTIONAL   :: On    ! On=T   will turn on  netCDF define mode
    LOGICAL, OPTIONAL   :: Off   ! Off=T  will turn off netCDF define mdoe
!
! !REMARKS:
!  This is a convenience wrapper for routines NcBegin_Def and NcEnd_Def in
!  NcdfUtil module m_netcdf_define_mod.F90.
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! If the ON switch is passed then ...
    IF ( PRESENT( On ) ) THEN
       IF ( On ) THEN
          CALL NcBegin_Def( fId )  ! Turn define mode on
          RETURN
       ELSE
          CALL NcEnd_Def( fId )    ! Turn define mode off
          RETURN
       ENDIF
    ENDIF

    ! If the OFF switch is passed then ,,,
    IF ( PRESENT( Off ) ) THEN
       IF ( Off ) THEN
          CALL NcEnd_Def( fId )      ! Turn define mode off
          RETURN
       ELSE
          CALL NcBegin_Def( fId )    ! Turn define mode on
          RETURN
       ENDIF
    ENDIF


  END SUBROUTINE Nc_Set_DefMode
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Read_Time
!
! !DESCRIPTION: Subroutine NC\_READ\_TIME reads the time variable of the
! given fID and returns the time slices and unit.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_READ_TIME( fID,     nTime,        timeUnit, &
                           timeVec, timeCalendar, RC       )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   )            :: fID
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)            :: nTime
    CHARACTER(LEN=*), INTENT(  OUT)            :: timeUnit
    REAL*8,           POINTER,       OPTIONAL  :: timeVec(:)
    CHARACTER(LEN=*), INTENT(  OUT), OPTIONAL  :: timeCalendar
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)            :: RC
!
! !REVISION HISTORY:
!  04 Nov 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                :: hasTime
    CHARACTER(LEN=255)     :: v_name             ! netCDF variable name
    CHARACTER(LEN=255)     :: a_name             ! netCDF attribute name
    CHARACTER(LEN=255)     :: a_val              ! netCDF attribute value
    INTEGER                :: st1d(1), ct1d(1)   ! For 1D arrays

    ! Arrays
    REAL*8 , ALLOCATABLE   :: tmpTime(:)

    !=================================================================
    ! NC_READ_TIME begins here
    !=================================================================

    ! Init
    RC      = 0
    nTime   = 0
    hasTime = .FALSE.

    ! Variable name
    v_name = "time"

    ! Check if dimension "time" exist
    hasTime = Ncdoes_Dim_Exist ( fID, TRIM(v_name) )

    ! If time dim not found, also check for dimension "date"
    IF ( .NOT. hasTime ) THEN
       v_name   = "date"
       hasTime = Ncdoes_Dim_Exist ( fID, TRIM(v_name) )
    ENDIF

    ! Return here if no time variable defined
    IF ( .NOT. hasTime ) RETURN

    ! Get dimension length
    CALL Ncget_Dimlen ( fID, TRIM(v_name), nTime )

    ! Read time/date units attribute
    a_name = "units"
    CALL NcGet_Var_Attributes( fID,          TRIM(v_name), &
                               TRIM(a_name), timeUnit     )

    ! Read time vector from file.
    IF ( PRESENT(timeVec) ) THEN
       IF ( ASSOCIATED(timeVec) ) DEALLOCATE ( timeVec)
       ALLOCATE ( tmpTime(nTime) )
       ALLOCATE ( timeVec(nTime) )
       st1d = (/ 1     /)
       ct1d = (/ nTime /)
       CALL NcRd( tmpTime, fID, TRIM(v_name), st1d, ct1d )
       timevec(:) = tmpTime
       DEALLOCATE(tmpTime)
    ENDIF

    ! Read calendar attribute
    IF ( PRESENT( timeCalendar ) ) THEN
       CALL NcGet_Var_Attributes( fId, v_name, 'calendar', timeCalendar )
    ENDIF

  END SUBROUTINE NC_READ_TIME
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Read_Var_Sp
!
! !DESCRIPTION: Subroutine NC\_READ\_VAR\_SP reads the given variable from the
! given fID and returns the corresponding variable values and units.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_READ_VAR_SP( fID, Var, nVar, varUnit, varVec, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   )            :: fID
    CHARACTER(LEN=*), INTENT(IN   )            :: var
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)            :: nVar
    CHARACTER(LEN=*), INTENT(  OUT)            :: varUnit
    REAL*4,           POINTER                  :: varVec(:)
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)            :: RC
!
! !REVISION HISTORY:
!  04 Nov 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    CALL NC_READ_VAR_CORE( fID, Var, nVar, varUnit, varVecSp=varVec, RC=RC )

  END SUBROUTINE NC_READ_VAR_SP
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Read_Var_Dp
!
! !DESCRIPTION: Subroutine NC\_READ\_VAR\_DP reads the given variable from the
! given fID and returns the corresponding variable values and units.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_READ_VAR_DP( fID, Var, nVar, varUnit, varVec, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   )            :: fID
    CHARACTER(LEN=*), INTENT(IN   )            :: var
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)            :: nVar
    CHARACTER(LEN=*), INTENT(  OUT)            :: varUnit
    REAL*8,           POINTER                  :: varVec(:)
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)            :: RC
!
! !REVISION HISTORY:
!  04 Nov 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    CALL NC_READ_VAR_CORE( fID, Var, nVar, varUnit, varVecDp=varVec, RC=RC )

  END SUBROUTINE NC_READ_VAR_DP
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Read_Var_Core
!
! !DESCRIPTION: Subroutine NC\_READ\_VAR\_CORE reads the given variable from the
! given fID and returns the corresponding variable values and units.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_READ_VAR_CORE( fID, Var, nVar, varUnit, varVecDp, varVecSp, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   )            :: fID
    CHARACTER(LEN=*), INTENT(IN   )            :: var
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)            :: nVar
    CHARACTER(LEN=*), INTENT(  OUT)            :: varUnit
    REAL*4,           POINTER,       OPTIONAL  :: varVecSp(:)
    REAL*8,           POINTER,       OPTIONAL  :: varVecDp(:)
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)            :: RC
!
! !REVISION HISTORY:
!  04 Nov 2012 - C. Keller   - Initial version
!  20 Feb 2015 - R. Yantosca - Need to add attType to Ncdoes_Attr_Exist
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL                :: hasVar
    CHARACTER(LEN=255)     :: v_name             ! netCDF variable name
    CHARACTER(LEN=255)     :: a_name             ! netCDF attribute name
    CHARACTER(LEN=255)     :: a_val              ! netCDF attribute value
    INTEGER                :: a_type             ! netCDF attribute type
    INTEGER                :: st1d(1), ct1d(1)   ! For 1D arrays
    INTEGER                :: I

    !=================================================================
    ! NC_READ_VAR_CORE begins here
    !=================================================================

    ! Init
    RC      = 0
    nVar    = 0
    hasVar  = .FALSE.

    ! Variable name
    v_name = var

    ! Check if variable exists
    hasVar = Ncdoes_Dim_Exist ( fID, TRIM(v_name) )

    ! Return here if variable not defined
    IF ( .NOT. hasVar ) RETURN

    ! Get dimension length
    CALL Ncget_Dimlen ( fID, TRIM(v_name), nVar )

    ! Read vector from file.
    IF ( PRESENT(VarVecSp) ) THEN
       IF ( ASSOCIATED( VarVecSp ) ) DEALLOCATE(VarVecSp)
       ALLOCATE ( VarVecSp(nVar) )
       st1d = (/ 1    /)
       ct1d = (/ nVar /)
       CALL NcRd( VarVecSp, fID, TRIM(v_name), st1d, ct1d )
    ENDIF
    IF ( PRESENT(VarVecDp) ) THEN
       IF ( ASSOCIATED( VarVecDp ) ) DEALLOCATE(VarVecDp)
       ALLOCATE ( VarVecDp(nVar) )
       st1d = (/ 1    /)
       ct1d = (/ nVar /)
       CALL NcRd( VarVecDp, fID, TRIM(v_name), st1d, ct1d )
    ENDIF

    ! Read units attribute. If unit attribute does not exist, return
    ! empty string (dimensionless vertical coordinates do not require
    ! a units attribute).
    a_name  = "units"
    hasVar  = Ncdoes_Attr_Exist ( fId, TRIM(v_name), TRIM(a_name), a_type )
    IF ( .NOT. hasVar ) THEN
       varUnit = ''
    ELSE
       CALL NcGet_Var_Attributes( fID,          TRIM(v_name), &
                                  TRIM(a_name), varUnit     )

       ! Check if the last character of VarUnit is the ASCII null character
       ! ("\0", ASCII value = 0), which is used to denote the end of a string.
       ! The ASCII null character may be introduced if the netCDF file was
       ! written using a language other than Fortran.  The compiler might
       ! interpret the null character as part of the string instead of as
       ! an empty space.  If the null space is there, then replace it with
       ! a Fortran empty string value (''). (bmy, 7/17/18)
       I = LEN_TRIM( VarUnit )
       IF ( ICHAR( VarUnit(I:I) ) == 0 ) THEN
          VarUnit(I:I) = ''
       ENDIF
    ENDIF

  END SUBROUTINE NC_READ_VAR_CORE
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Read_Arr
!
! !DESCRIPTION: Routine NC\_READ\_ARR reads variable ncVar into a 4-D array
! (lon,lat,lev,time). Domain boundaries can be provided by input arguments
! lon1,lon2, lat1,lat2, lev1,lev2, and time1,time2. The level and time bounds
! are optional and can be set to zero (lev1=0 and/or time1=0) for data with
! undefined level/time coordinates.
!\\
!\\
! The default behavior for time slices is to read all slices (time1:time2),
! and pass all of them to the output array. It is also possible to assign
! specific weights (wgt1 and wgt2) to the two time slices time1 and time2,
! respectively. In this case, only those two slices will be read and merged
! using the given weights. The output array will then contain only one time
! dimension. Negative weights are currently not supported and will be ignored,
! e.g. providing negative weights has the same effect as providing no weights
! at all.
!\\
!\\
! If the passed variable contains attribute names `offset` and/or
! `scale\_factor`, those operations will be applied to the data array
! before returning it.
!\\
!\\
! Missing values in the netCDF file are replaced with value 'MissVal'
! (default = 0). Currently, the routine identifies attributes 'missing\_value'
! and '\_FillValue' as missing values.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_READ_ARR( fID,    ncVar,   lon1,    lon2,  lat1,  &
                          lat2,   lev1,    lev2,    time1, time2, &
                          ncArr,  VarUnit, MissVal, wgt1,  wgt2,  &
                          ArbIdx, RC                               )
!
! !USES:
!
    USE CHARPAK_MOD, ONLY : TRANLC
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)            :: fID
    CHARACTER(LEN=*), INTENT(IN)            :: ncVar        ! variable to read
    INTEGER,          INTENT(IN)            :: lon1,  lon2
    INTEGER,          INTENT(IN)            :: lat1,  lat2
    INTEGER,          INTENT(IN)            :: lev1,  lev2
    INTEGER,          INTENT(IN)            :: time1, time2
    REAL*4,           INTENT(IN ), OPTIONAL :: MissVal
    REAL*4,           INTENT(IN ), OPTIONAL :: wgt1
    REAL*4,           INTENT(IN ), OPTIONAL :: wgt2
    INTEGER,          INTENT(IN ), OPTIONAL :: ArbIdx      ! Index of arbitrary additional dimension (-1 if none)
!
! !OUTPUT PARAMETERS:
!
    ! Array to write data
    REAL*4,           POINTER               :: ncArr(:,:,:,:)

    ! Optional output
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: VarUnit
!
! !INPUT/OUTPUT PARAMETERS:
!
    ! Error handling
    INTEGER,          INTENT(INOUT)         :: RC
!
! !REVISION HISTORY:
!  27 Jul 2012 - C. Keller - Initial version
!  18 Jan 2012 - C. Keller - Now reads 4D, 3D, and 2D arrays, with
!                            optional dimensions level and time.
!  18 Apr 2012 - C. Keller - Now also read & apply offset and scale factor
!  27 Feb 2015 - C. Keller - Added weights.
!  22 Sep 2015 - C. Keller - Added arbitrary dimension index.
!  20 Nov 2015 - C. Keller - Bug fix: now read times if weights need be applied.
!  23 Nov 2015 - C. Keller - Initialize all temporary arrays to 0.0 when allocating
!  09 Jan 2017 - C. Keller - Bug fix: store time-weighted arrays in temporary array
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
    CHARACTER(LEN=255)     :: v_name    ! netCDF variable name
    CHARACTER(LEN=255)     :: a_name    ! netCDF attribute name
    CHARACTER(LEN=255)     :: a_val     ! netCDF attribute value
    INTEGER                :: a_type    ! netCDF attribute type
    REAL*8                 :: corr      ! netCDF attribute value

    ! Arrays for netCDF start and count values
    INTEGER                :: I, nRead, l1, l2
    INTEGER                :: ndims
    INTEGER                :: nlon,  nlat, nlev, ntime, arbdim
    INTEGER                :: nclev, nctime
    INTEGER                :: s1, s2, s3, s4, s5
    INTEGER                :: n1, n2, n3, n4, n5
    INTEGER                :: nt, st, tdim, sti, nti
    INTEGER                :: st2d(2), ct2d(2)   ! For 2D arrays
    INTEGER                :: st3d(3), ct3d(3)   ! For 3D arrays
    INTEGER                :: st4d(4), ct4d(4)   ! For 4D arrays
    INTEGER                :: st5d(5), ct5d(5)   ! For 5D arrays

    ! Temporary arrays
    REAL*4, ALLOCATABLE    :: TMPARR_5D(:,:,:,:,:)
    REAL*4, ALLOCATABLE    :: WGTARR_5D(:,:,:,:,:)
    REAL*4, ALLOCATABLE    :: TMPARR_4D(:,:,:,:)
    REAL*4, ALLOCATABLE    :: WGTARR_4D(:,:,:,:)
    REAL*4, ALLOCATABLE    :: TMPARR_3D(:,:,:)
    REAL*4, ALLOCATABLE    :: WGTARR_3D(:,:,:)
    REAL*4, ALLOCATABLE    :: TMPARR_2D(:,:)

    ! Logicals
    LOGICAL                :: FlipZ
    LOGICAL                :: ReadAtt

    ! Missing value
    REAL*8                 :: miss8
    REAL*4                 :: miss4
    REAL*4                 :: MissValue

    ! Weights
    LOGICAL                :: ApplyWeights
    REAL*4                 :: weight1, weight2

    ! For error handling
    CHARACTER(LEN=255)     :: LOC, MSG

    !=================================================================
    ! NC_READ_ARR begins here
    !=================================================================

    !-----------------------------------------------------------------
    ! Initialize
    !-----------------------------------------------------------------

    ! For error handling
    LOC = 'NC_READ_ARR ("ncdf_mod.F")'

    ! Eventually deallocate output array
    IF ( ASSOCIATED ( ncArr ) ) DEALLOCATE ( ncArr )

    ! weights to be applied to time1 and time2 (if any):
    weight1 = -999.0
    weight2 = -999.0
    IF(PRESENT(wgt1)) weight1 = wgt1
    IF(PRESENT(wgt2)) weight2 = wgt2

    ! apply weights?
    IF ( time1 > 0 .AND. weight1 >= 0.0 ) THEN
       ApplyWeights = .TRUE.
    ELSE
       ApplyWeights = .FALSE.
    ENDIF

    ! # of horizontal dimensions to read
    nLon = lon2 - lon1 + 1
    nLat = lat2 - lat1 + 1

    ! # of vertical levels
    FlipZ = .FALSE. ! Flip z-axis?
    l1    = lev1    ! Lower level to be read
    l2    = lev2    ! Upper level to be read
    IF ( lev1 > 0 ) THEN

       ! Check if we need to flip the vertical axis
       IF ( lev1 > lev2 ) THEN
          FlipZ = .TRUE.
          l1    = lev2
          l2    = lev1
       ENDIF

       ! Number of levels to be read
       nLev = l2 - l1 + 1

    ! no vertical levels:
    ELSE
       nLev = 0
    ENDIF

    ! # of time slices
    ! read all time slices time1:time2:
    IF ( time1 > 0 .AND. weight1 < 0.0 ) THEN
       ntime = time2 - time1 + 1
    ! Interpolate amongs 2 time slices:
    ELSEIF ( ApplyWeights ) THEN
       ntime = 1
    ! no time dimension:
    ELSE
       ntime = 0
    ENDIF

    ! # of arbitrary other dimensions
    arbdim = -1
    IF ( PRESENT(ArbIdx) ) THEN
       IF ( ArbIdx > 0 ) THEN
          arbdim = ArbIdx
       ENDIF
    ENDIF

    ! Set dimensions of output array
    ! --> must have at least dimension 1
    nclev  = max(nlev ,1)
    nctime = max(ntime,1)

    ! set total number of dimensions to be read. This is at least 2 and
    ! at most 5.
    ndims = 2
    if ( nlev   > 0 ) ndims = ndims + 1
    if ( ntime  > 0 ) ndims = ndims + 1
    if ( arbdim > 0 ) ndims = ndims + 1

    !----------------------------------------
    ! Read array
    !----------------------------------------

    ! Variable name
    v_name = TRIM(ncVar)

    ! Allocate the output array
    ALLOCATE ( ncArr( nLon, nLat, ncLev, ncTime ) )
    ncArr = 0.0

    ! Define number of required reads and time dimension on temporary array
    nRead = 1
    IF ( ntime > 0 ) THEN
       IF ( ApplyWeights ) THEN
          nRead = 2
          nt    = 2
       ELSE
          nRead = 1
          nt    = ntime
       ENDIF
    ENDIF

    !----------------------------------------
    ! Read 5D array:
    IF ( ndims == 5 ) THEN

       ! Allocate array. If time weights are applied, the two
       ! time slices are read into TMPARR_5D and then temporarily
       ! stored in WGTARR_5D. Same applies to 4D and 3D below.
       ! (ckeller, 01/09/17)
       IF ( ApplyWeights ) THEN
          ALLOCATE ( TMPARR_5D( nlon, nlat, nlev, 1, 1 ) )
          TMPARR_5D = 0.0
          ALLOCATE ( WGTARR_5D( nlon, nlat, nlev, nt, 1 ) )
          WGTARR_5D = 0.0
       ELSE
          ALLOCATE ( TMPARR_5D( nlon, nlat, nlev, nt, 1 ) )
          TMPARR_5D = 0.0
       ENDIF

       ! Set default start/end indeces
       s1 = lon1
       n1 = nlon
       s2 = lat1
       n2 = nlat
       s3 = l1
       n3 = nlev
       s5 = arbdim
       n5 = 1

       ! Read arrays from file
       DO I = 1, nRead

          ! time index
          IF ( .NOT. ApplyWeights ) THEN
             s4 = time1
             n4 = ntime
          ELSE
             IF ( I == 1 ) THEN
                s4 = time1
             ELSE
                s4 = time2
             ENDIF
             n4 = 1
          ENDIF

          st5d = (/ s1, s2, s3, s4, s5 /)
          ct5d = (/ n1, n2, n3, n4, n5 /)
          CALL NcRd( TMPARR_5D, fId, TRIM(v_name), st5d, ct5d )

          ! Eventually pass time weighted arrays to temporary array
          IF ( ApplyWeights ) THEN
             WGTARR_5D(:,:,:,I,:) = TMPARR_5D(:,:,:,1,:)
          ENDIF

       ENDDO

       ! Pass to output array. Eventually apply time weights.
       IF ( ApplyWeights ) THEN
          ncArr(:,:,:,1) = WGTARR_5D(:,:,:,1,1) * weight1 &
                         + WGTARR_5D(:,:,:,2,1) * weight2
       ELSE
          ncArr(:,:,:,:) = TMPARR_5D(:,:,:,:,1)
       ENDIF

       ! Cleanup
       DEALLOCATE(TMPARR_5D)
       IF(ALLOCATED(WGTARR_5D)) DEALLOCATE(WGTARR_5D)
    ENDIF

    !----------------------------------------
    ! Read 4D array:
    ! This can be:
    ! - lon,lat,lev,time
    ! - lon,lat,lev,arb
    ! - lon,lat,time,arb
    IF ( ndims == 4 ) THEN

       ! Allocate temporary array
       s1    = lon1
       n1    = nlon
       s2    = lat1
       n2    = nlat
       tdim  = -1

       ! 3rd and 4th dim

       ! lev is defined
       IF ( nlev > 0 ) THEN
          s3   = l1
          n3   = nlev
          ! plus time...
          IF ( ntime > 0 ) THEN
             n4   = nt
             tdim = 4
          ! ... or plus arbitrary dim
          ELSE
             s4 = arbdim
             n4 = 1
          ENDIF

       ! lev not defined: time + arbitrary dim
       ELSE
          n3 = nt
          tdim = 3
          s4 = arbdim
          n4 = 1
       ENDIF

       IF ( ApplyWeights ) THEN
          ALLOCATE ( WGTARR_4D(n1,n2,n3,n4) )
          WGTARR_4D = 0.0
          IF ( tdim == 3 ) THEN
             ALLOCATE ( TMPARR_4D(n1,n2,1,n4) )
             TMPARR_4D = 0.0
          ELSEIF ( tdim == 4 ) THEN
             ALLOCATE ( TMPARR_4D(n1,n2,n3,1) )
             TMPARR_4D = 0.0
          ENDIF

       ELSE
          ALLOCATE ( TMPARR_4D(n1,n2,n3,n4) )
          TMPARR_4D = 0.0
       ENDIF

       ! Read arrays from file
       DO I = 1, nRead

          ! time index
          IF ( .NOT. ApplyWeights ) THEN
             sti = time1
             nti = ntime
          ELSE
             IF ( I == 1 ) THEN
                sti = time1
             ELSE
                sti = time2
             ENDIF
             nti = 1
          ENDIF

          ! need to adjust time index: this is either 3rd or 4th dimension:
          IF ( tdim == 3 ) THEN
             s3 = sti
             n3 = nti
          ELSEIF ( tdim == 4 ) THEN
             s4 = sti
             n4 = nti
          ENDIF

          st4d = (/ s1, s2, s3, s4 /)
          ct4d = (/ n1, n2, n3, n4 /)

          ! Read data from disk
          CALL NcRd( TMPARR_4D, fId, TRIM(v_name), st4d, ct4d )

          ! Eventually pass time weighted arrays to temporary array
          IF ( ApplyWeights ) THEN
             IF ( tdim == 3 ) THEN
                WGTARR_4D(:,:,I,:) = TMPARR_4D(:,:,1,:)
             ELSEIF ( tdim == 4 ) THEN
                WGTARR_4D(:,:,:,I) = TMPARR_4D(:,:,:,1)
             ENDIF
          ENDIF
       ENDDO

       ! Pass to output array. Eventually apply time weights.
       IF ( ApplyWeights ) THEN
          IF ( tdim == 3 ) THEN
             ncArr(:,:,:,1) = WGTARR_4D(:,:,1,:) * weight1 &
                            + WGTARR_4D(:,:,2,:) * weight2
          ELSEIF ( tdim == 4 ) THEN
             ncArr(:,:,:,1) = WGTARR_4D(:,:,:,1) * weight1 &
                            + WGTARR_4D(:,:,:,2) * weight2
          ENDIF
       ELSE
          ncArr(:,:,:,:) = TMPARR_4D(:,:,:,:)
       ENDIF

       ! Cleanup
       DEALLOCATE(TMPARR_4D)
       IF(ALLOCATED(WGTARR_4D)) DEALLOCATE(WGTARR_4D)
    ENDIF

    !----------------------------------------
    ! Read 3D array:
    ! This can be:
    ! - lon,lat,lev
    ! - lon,lat,time
    ! - lon,lat,arb
    IF ( ndims == 3 ) THEN

       ! Allocate temporary array
       s1    = lon1
       n1    = nlon
       s2    = lat1
       n2    = nlat
       tdim  = -1

       ! 3rd dim:
       ! - lev is defined:
       IF ( nlev > 0 ) THEN
          s3   = l1
          n3   = nlev
       ! - time is defined:
       ELSEIF ( ntime > 0 ) THEN
          n3   = nt
          tdim = 3
       ! - arbitrary dimension is defined:
       ELSEIF ( arbdim > 0 ) THEN
          s3   = arbdim
          n3   = 1
       ENDIF

       IF ( ApplyWeights ) THEN
          ALLOCATE ( TMPARR_3D(n1,n2,1) )
          TMPARR_3D = 0.0
          ALLOCATE ( WGTARR_3D(n1,n2,n3) )
          WGTARR_3D = 0.0
       ELSE
          ALLOCATE ( TMPARR_3D(n1,n2,n3) )
          TMPARR_3D = 0.0
       ENDIF

       ! Read arrays from file
       DO I = 1, nRead

          ! time index
          IF ( tdim  == 3 ) THEN
             IF ( .NOT. ApplyWeights ) THEN
                s3 = time1
                n3 = ntime
             ELSE
                IF ( I == 1 ) THEN
                   s3 = time1
                ELSE
                   s3 = time2
                ENDIF
                n3 = 1
             ENDIF
          ENDIF

          st3d = (/ s1, s2, s3 /)
          ct3d = (/ n1, n2, n3 /)
          CALL NcRd( TMPARR_3D, fId, TRIM(v_name), st3d, ct3d )

          ! Eventually pass time weighted arrays to temporary array
          IF ( ApplyWeights ) THEN
           WGTARR_3D(:,:,I) = TMPARR_3D(:,:,1)
          ENDIF

       ENDDO

       ! Pass to output array. Eventually apply time weights.
       IF ( ApplyWeights ) THEN
          ncArr(:,:,1,1) = WGTARR_3D(:,:,1) * weight1 &
                         + WGTARR_3D(:,:,2) * weight2
       ELSE
          IF ( tdim == 3 ) THEN
             ncArr(:,:,1,:) = TMPARR_3D(:,:,:)
          ELSE
             ncArr(:,:,:,1) = TMPARR_3D(:,:,:)
          ENDIF
       ENDIF

       ! Cleanup
       IF(ALLOCATED(TMPARR_3D)) DEALLOCATE(TMPARR_3D)
       IF(ALLOCATED(WGTARR_3D)) DEALLOCATE(WGTARR_3D)
    ENDIF

    !----------------------------------------
    ! Read a 2D array (lon and lat only):
    IF ( ndims == 2 ) THEN
       ALLOCATE ( TMPARR_2D( nLon, nLat ) )
       TMPARR_2D = 0.0
       st2d      = (/ lon1, lat1 /)
       ct2d      = (/ nlon, nlat /)
       CALL NcRd( TMPARR_2D, fId, TRIM(v_name), st2d, ct2d )
       ncArr(:,:,1,1) = TMPARR_2D(:,:)
       DEALLOCATE(TMPARR_2D)
    ENDIF

    ! ------------------------------------------
    ! Eventually apply scale / offset factors
    ! ------------------------------------------

    ! Check for scale factor
    a_name  = "scale_factor"
    ReadAtt = Ncdoes_Attr_Exist ( fId, TRIM(v_name), TRIM(a_name), a_type )

    IF ( ReadAtt ) THEN
       CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),corr)
       ncArr(:,:,:,:) = ncArr(:,:,:,:) * corr
    ENDIF

    ! Check for offset factor
    a_name  = "add_offset"
    ReadAtt = Ncdoes_Attr_Exist ( fId, TRIM(v_name), TRIM(a_name), a_type )

    IF ( ReadAtt ) THEN
       CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),corr)
       ncArr(:,:,:,:) = ncArr(:,:,:,:) + corr
    ENDIF

    ! ------------------------------------------
    ! Check for filling values
    ! NOTE: Test for REAL*4 and REAL*8
    ! ------------------------------------------

    ! Define missing value
    IF ( PRESENT(MissVal) ) THEN
       MissValue = MissVal
    ELSE
       MissValue = 0.0
    ENDIF

    ! 1: 'missing_value'
    a_name  = "missing_value"
    ReadAtt = Ncdoes_Attr_Exist ( fId, TRIM(v_name), TRIM(a_name), a_type )
    IF ( ReadAtt ) THEN
       IF ( a_type == NF_REAL ) THEN
          CALL NcGet_Var_Attributes( fId, TRIM(v_name), TRIM(a_name), miss4 )
          WHERE ( ncArr == miss4 )
             ncArr = MissValue
          END WHERE
       ELSE IF ( a_type == NF_DOUBLE ) THEN
          CALL NcGet_Var_Attributes( fId, TRIM(v_name), TRIM(a_name), miss8 )
          miss4 = REAL( miss8 )
          WHERE ( ncArr == miss4 )
             ncArr = MissValue
          END WHERE
       ENDIF
    ENDIF

    ! 2: '_FillValue'
    a_name  = "_FillValue"
    ReadAtt = Ncdoes_Attr_Exist ( fId, TRIM(v_name), TRIM(a_name), a_type )
    IF ( ReadAtt ) THEN
       IF ( a_type == NF_REAL ) THEN
          CALL NcGet_Var_Attributes( fId, TRIM(v_name), TRIM(a_name), miss4 )
          WHERE ( ncArr == miss4 )
             ncArr = MissValue
          END WHERE
       ELSE IF ( a_type == NF_DOUBLE ) THEN
          CALL NcGet_Var_Attributes( fId, TRIM(v_name), TRIM(a_name), miss8 )
          miss4 = REAL( miss8 )
          WHERE ( ncArr == miss4 )
             ncArr = MissValue
          END WHERE
       ENDIF
    ENDIF

    ! ------------------------------------------
    ! Flip z-axis if needed
    ! ------------------------------------------
    IF ( FlipZ ) THEN
       ncArr(:,:,:,:) = ncArr(:,:,ncLev:1:-1,:)
    ENDIF

    ! ----------------------------
    ! Read optional arguments
    ! ----------------------------

    ! Read units
    IF ( PRESENT(VarUnit) )THEN
       a_name = "units"
       CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)
       VarUnit = TRIM(a_val)

       ! Check if the last character of VarUnit is the ASCII null character
       ! ("\0", ASCII value = 0), which is used to denote the end of a string.
       ! The ASCII null character may be introduced if the netCDF file was
       ! written using a language other than Fortran.  The compiler might
       ! interpret the null character as part of the string instead of as
       ! an empty space.  If the null space is there, then replace it with
       ! a Fortran empty string value (''). (bmy, 7/17/18)
       I = LEN_TRIM( VarUnit )
       IF ( ICHAR( VarUnit(I:I) ) == 0 ) THEN
          VarUnit(I:I) = ''
       ENDIF
    ENDIF

    !=================================================================
    ! Cleanup and quit
    !=================================================================

    ! Return w/ success
    RC = 0

  END SUBROUTINE NC_READ_ARR
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Read_Time_yyyymmddhhmm
!
! !DESCRIPTION: Returns a vector containing the datetimes (YYYYMMDDhhmm) of
! all time slices in the netCDF file.
!\\
! !INTERFACE:
!
  SUBROUTINE NC_READ_TIME_YYYYMMDDhhmm( fID,              nTime,    &
                                        all_YYYYMMDDhhmm, timeUnit, &
                                        refYear,          RC )
!
! !USES:
!
    USE JULDAY_MOD, ONLY : JULDAY, CALDATE
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   )           :: fID
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*8,           POINTER                 :: all_YYYYMMDDhhmm(:)
    CHARACTER(LEN=*), INTENT(  OUT), OPTIONAL :: timeUnit
    INTEGER,          INTENT(  OUT), OPTIONAL :: refYear
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: nTime
    INTEGER,          INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  27 Jul 2012 - C. Keller   - Initial version
!  09 Oct 2014 - C. Keller   - Now also support 'minutes since ...'
!  05 Nov 2014 - C. Keller   - Bug fix if reference datetime is in minutes.
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  05 Apr 2017 - C. Keller   - Now also support 'seconds since ...'
!  10 Apr 2017 - R. Yantosca - Now return times in YYYYMMDDhhmm
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    CHARACTER(LEN=255)  :: ncUnit
    INTEGER             :: refYr, refMt, refDy, refHr, refMn, refSc
    INTEGER             :: T, YYYYMMDD, hhmmss
    REAL*8              :: realrefDy, refJulday, tJulday

    ! Pointers
    REAL*8,   POINTER   :: tVec(:)

    !=================================================================
    ! NC_READ_TIME_YYYYMMDDhhmm begins here
    !=================================================================

    ! Init values
    RC = 0
    tVec => NULL()
    IF ( PRESENT(TimeUnit) ) TimeUnit = ''
    IF ( PRESENT(refYear ) ) refYear  = 0

    ! Read time vector
    CALL NC_READ_TIME ( fID, nTime, ncUnit, timeVec=tVec, RC=RC )
    IF ( RC/=0 ) RETURN

    ! If nTime is zero, return here!
    IF ( nTime == 0 ) RETURN

    ! Get reference date in julian days
    CALL NC_GET_REFDATETIME ( ncUnit, refYr, refMt, &
                              refDy,  refHr, refMn, refSc, RC )
    IF ( RC /= 0 ) RETURN
    realrefDy =         refDy              &
              + ( MAX(0,refHr) / 24d0    ) &
              + ( MAX(0,refMn) / 1440d0  ) &
              + ( MAX(0,refSc) / 86400d0 )
    refJulday = JULDAY ( refYr, refMt, realrefDy )

    ! NOTE: It seems that there is an issue with reference dates
    ! between 1800 and 1901: the respective time stamps all seem to
    ! be off by one day (this problem doesn't appear for netCDF files
    ! with reference date zero, i.e. hours since 1-1-1)!
    ! I'm not sure what causes this problem, but adding one day to
    ! reference dates that lie between 1600 and 1900 seems to fix the
    ! problem.
    ! TODO: requires more testing!
    IF ( refYr <= 1900 .AND. refYr >= 1600 ) THEN
       refJulday = refJulday + 1.0
       !PRINT *, 'Reference julian day increased by one day!!!'
    ENDIF

    ! Get calendar dates
    IF ( ASSOCIATED ( all_YYYYMMDDhhmm ) ) DEALLOCATE( all_YYYYMMDDhhmm )
    ALLOCATE( all_YYYYMMDDhhmm(nTime) )
    all_YYYYMMDDhhmm = 0.0d0

    ! Construct julian date for every available time slice. Make sure it is
    ! in the proper 'units', e.g. in days, hours or minutes, depending on
    ! the reference unit.
    DO T = 1, nTime
       tJulDay = tVec(T)
       IF ( refHr >= 0 ) tJulday = tJulday / 24.d0
       IF ( refMn >= 0 ) tJulday = tJulday / 60.d0
       IF ( refSc >= 0 ) tJulday = tJulday / 60.d0
       tJulday = tJulday + refJulday
       CALL CALDATE ( tJulday, YYYYMMDD, hhmmss )
       all_YYYYMMDDhhmm(T) = ( DBLE( YYYYMMDD ) * 1d4   ) + &
                             ( DBLE( hhmmss     / 100 ) )
    ENDDO

    ! Cleanup
    IF ( ASSOCIATED( tVec ) ) DEALLOCATE( tVec )

    ! Return
    IF ( PRESENT(timeUnit) ) timeUnit = ncUnit
    IF ( PRESENT(refYear ) ) refYear  = refYr
    RC = 0

  END SUBROUTINE NC_READ_TIME_YYYYMMDDhhmm
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Get_RefDateTime
!
! !DESCRIPTION: Returns the reference datetime (tYr / tMt / tDy / tHr /
! tMn ) of the provided time unit. For now, supported formats are
! "days since YYYY-MM-DD", "hours since YYYY-MM-DD HH:MM:SS", and
! "minutes since YYYY-MM-DD HH:NN:SS". For times in days since refdate,
! the returned reference hour rHr is set to -1. The same applies for the
! reference minute for units in days / hours since XXX.
!\\
! !INTERFACE:
!
  SUBROUTINE NC_GET_REFDATETIME( tUnit, tYr, tMt, tDy, tHr, tMn, tSc, RC )
!
! !USES:
!
    USE CHARPAK_MOD, ONLY : TRANLC
!
! !INPUT PARAMETERS:
!
    ! Required
    CHARACTER(LEN=*), INTENT( IN)    :: tUnit
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)    :: tYr
    INTEGER,          INTENT(OUT)    :: tMt
    INTEGER,          INTENT(OUT)    :: tDy
    INTEGER,          INTENT(OUT)    :: tHr
    INTEGER,          INTENT(OUT)    :: tMn
    INTEGER,          INTENT(OUT)    :: tSc
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  18 Jan 2012 - C. Keller - Initial version
!  09 Oct 2014 - C. Keller - Now also support 'minutes since ...'
!  20 Nov 2015 - C. Keller - Now also support 'seconds since ...'
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)    :: LOC, MSG
    CHARACTER(LEN=255)    :: MIRRUNIT
    INTEGER               :: TTYPE, STAT, L1, L2
    INTEGER               :: MINLEN, STRLEN, I

    !=================================================================
    ! NC_GET_REFDATETIME starts here
    !=================================================================

    ! Init
    LOC = 'NC_GET_REFDATETIME (ncdf_mod.F)'

    ! ----------------------------------------------------------------------
    ! Determine time unit type
    ! ----------------------------------------------------------------------

    ! Mirror time unit and convert to lower case
    MIRRUNIT = tUnit
    CALL TRANLC( MIRRUNIT )

    ! Check for reference time unit '(days, hours, minutes) since ...'
    ! Set beginning of reference date according to the unit and define
    ! minimum string length required by unit.

    ! 'days since YYYY-M-D'
    IF ( MIRRUNIT(1:10) == 'days since' ) THEN
       TTYPE  = 1
       L1     = 12
       MINLEN = 19

    ! 'hours since YYYY-M-D h:m:s'
    ELSEIF ( MIRRUNIT(1:11) == 'hours since' ) THEN
       TTYPE  = 2
       L1     = 13
       MINLEN = 26

    ! 'minutes since YYYY-M-D h:m:s'
    ELSEIF ( MIRRUNIT(1:13) == 'minutes since' ) THEN
       TTYPE  = 3
       L1     = 15
       MINLEN = 28

    ! 'seconds since YYYY-M-D h:m:s'
    ELSEIF ( MIRRUNIT(1:13) == 'seconds since' ) THEN
       TTYPE  = 4
       L1     = 15
       MINLEN = 28

    ! Return w/ error otherwise
    ELSE
       PRINT *, 'Invalid time unit: ' // TRIM(tUnit)
       RC = -999; RETURN
    ENDIF

    ! Check if time string is long enough or not
    STRLEN = LEN(tUnit)
    IF ( STRLEN < MINLEN ) THEN
       PRINT *, 'Time unit string too short: ' // TRIM(tUnit)
       RC = -999; RETURN
    ENDIF

    ! ----------------------------------------------------------------------
    ! Determine reference time/date
    ! Get the year, month, day and hour from the string
    ! '... since YYYY-MM-DD hh:mm:ss

    ! Read reference year, i.e. from beginning of date string until
    ! first separator sign (-).
    DO I=L1,STRLEN
       IF(tUnit(I:I) == '-') EXIT
    ENDDO
    L2 = I-1

    READ( tUnit(L1:L2),'(i4)', IOSTAT=STAT ) tYr
    IF ( STAT /= 0 ) THEN
       PRINT *, 'Invalid year in ' // TRIM(tUnit)
       RC = -999; RETURN
    ENDIF

    ! Advance in date string: now read reference month.
    L1 = L2 + 2
    DO I=L1,STRLEN
       IF(tUnit(I:I) == '-') EXIT
    ENDDO
    L2 = I-1
    READ( tUnit(L1:L2), '(i2)', IOSTAT=STAT ) tMt
    IF ( STAT /= 0 ) THEN
       PRINT *, 'Invalid month in ' // TRIM(tUnit)
       RC = -999; RETURN
    ENDIF

    ! Advance in date string: now read reference day.
    L1 = L2 + 2
    DO I=L1,STRLEN
       IF(tUnit(I:I) == ' ') EXIT
    ENDDO
    L2 = I-1
    READ( tUnit(L1:L2), '(i2)', IOSTAT=STAT ) tDy
    IF ( STAT /= 0 ) THEN
       PRINT *, 'Invalid day in ' // TRIM(tUnit)
       RC = -999; RETURN
    ENDIF

    ! Get reference hour only if 'hours/minutes/seconds since'.
    IF ( TTYPE > 1 ) THEN

       ! Reference hour
       L1 = L2 + 2
       DO I=L1,STRLEN
          IF(tUnit(I:I) == ':') EXIT
       ENDDO
       L2 = I-1
       READ( tUnit(L1:L2), '(i2)', IOSTAT=STAT ) tHr
       IF ( STAT /= 0 ) THEN
          PRINT *, 'Invalid hour in ', TRIM(tUnit)
          RC = -999; RETURN
       ENDIF

    ELSE
       ! Set reference hour to -1
       tHr = -1
    ENDIF

    ! Get reference minute only if 'minutes since...'
    IF ( TTYPE>2 ) THEN

       ! Reference minute
       L1 = L2 + 2
       DO I=L1,STRLEN
          IF(tUnit(I:I) == ':') EXIT
       ENDDO
       L2 = I-1
       READ( tUnit(L1:L2), '(i2)', IOSTAT=STAT ) tMn
       IF ( STAT /= 0 ) THEN
          PRINT *, 'Invalid minute in ', TRIM(tUnit)
          RC = -999; RETURN
       ENDIF

    ELSE
       ! Set reference minute to -1
       tMn = -1
    ENDIF

    ! Get reference minute only if 'seconds since...'
    IF ( TTYPE>3 ) THEN

       ! Reference second
       L1 = L2 + 2
       DO I=L1,STRLEN
          IF(tUnit(I:I) == ':') EXIT
       ENDDO
       L2 = I-1
       READ( tUnit(L1:L2), '(i2)', IOSTAT=STAT ) tSc
       IF ( STAT /= 0 ) THEN
          PRINT *, 'Invalid second in ', TRIM(tUnit)
          RC = -999; RETURN
       ENDIF

    ELSE
       ! Set reference second to -1
       tSc = -1
    ENDIF

    ! Return w/ success
    RC = 0

  END SUBROUTINE NC_GET_REFDATETIME
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Tidx
!
! !DESCRIPTION: Routine GET\_TIDX returns the index with the specified time
! for a given time vector.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_TIDX( TDIM, TIMEVEC,  TTYPE, TOFFSET, &
                       YEAR, MONTH,    DAY,   HOUR,    &
                       TIDX, TDIMREAD, RC )
!
! !INPUT PARAMETERS:
!
    ! Required
    INTEGER, INTENT(   IN)           :: TDIM
    INTEGER, INTENT(   IN)           :: TTYPE
    REAL*8,  INTENT(   IN)           :: TOFFSET
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT)           :: TIMEVEC(TDIM)
    INTEGER, INTENT(INOUT)           :: YEAR
    INTEGER, INTENT(INOUT)           :: MONTH
    INTEGER, INTENT(INOUT)           :: DAY
    INTEGER, INTENT(INOUT)           :: HOUR
    INTEGER, INTENT(INOUT)           :: RC
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(  OUT)           :: TIDX
    INTEGER, INTENT(  OUT)           :: TDIMREAD
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
    INTEGER                :: II, iiDiff, minDiff
    REAL*8                 :: TAU
    CHARACTER(LEN=255)     :: MSG, LOC

    !=================================================================
    ! GET_TIDX starts here
    !=================================================================

    ! Init
    LOC = 'GET_TIDX (ncdf_mod.F)'
    TIDX = 0
    minDiff = -999

    !-----------------------------------------------------------------
    ! If year is given, compare netcdf-tau against desired tau
    !-----------------------------------------------------------------
    IF ( YEAR > 0 ) THEN

       ! Restrict month, day and hour to valid values
       MONTH = MIN ( MAX( 1, MONTH ), 12 )
       DAY   = MIN ( MAX( 1, DAY   ), 31 )
       HOUR  = MIN ( MAX( 0, HOUR  ), 23 )

       ! Read desired tau => hours relative to G-C reference time
       TAU = GET_TAU0( MONTH, DAY, YEAR, HOUR )

       ! Convert to 'hours since ...' if unit is 'days since ...'
       IF ( TTYPE == 2 ) THEN
          TIMEVEC(:) = TIMEVEC(:) * 24
       ENDIF

       ! Convert time stamps to hours since G-C reference time
       TIMEVEC(:) = TIMEVEC(:) + INT(TOFFSET)

       ! Compare wanted tau to tau's of ncdf-file.
       ! Loop over all time stamps and check which one is closest
       ! to the specified one. Print a warning if time stamps don't
       ! match!
       DO II = 1, TDIM

          ! Difference between time stamps
          iiDiff = ABS( TIMEVEC(II) - INT(TAU) )

          ! Check if this is closest time stamp so far, and save this
          ! index and difference
          IF ( iiDiff < minDiff .OR. II == 1 ) THEN
             minDiff = iiDiff
             TIDX    = II
          ENDIF

          ! Exit loop if difference is zero
          IF ( minDiff == 0 ) EXIT

       ENDDO

       ! Warning if time stamps did not match
       IF ( minDiff /= 0 ) THEN
          PRINT *, 'In NCDF_MOD: Time stamp not found ' // &
                   'take closest timestamp!'
       ENDIF

       ! Set number of time stamps to be read to 1
       TDIMREAD = 1

    !-----------------------------------------------------------------
    ! If only month is given, assume netCDF file to contain monthly
    ! data and pick the desired month.
    !-----------------------------------------------------------------
    ELSEIF ( MONTH > 0 ) THEN

       ! Check if it's indeed monthly data:
       IF ( TDIM /= 12 ) THEN
          PRINT *, 'Array is not monthly '
          RC = -999; RETURN
       ENDIF

       ! Set time index to specified month
       TIDX = MONTH

       ! Set number of time stamps to be read to 1
       TDIMREAD = 1

    !-----------------------------------------------------------------
    ! If hour is given, assume netCDF file to contain hourly data
    ! and pick the desired hour.
    !-----------------------------------------------------------------
    ELSEIF ( HOUR >= 0 ) THEN

       ! Check if it's indeed hourly data:
       IF ( TDIM /= 24 ) THEN
          PRINT *, 'Array is not hourly'
          RC = -999; RETURN
       ENDIF

       ! Set time index to specified hour (+1 since hour 0 is idx 1)
       TIDX = HOUR + 1

       ! Set number of time stamps to be read to 1
       TDIMREAD = 1

    !-----------------------------------------------------------------
    ! Otherwise, read all time dimensions
    !-----------------------------------------------------------------
    ELSE
       TIDX     = 1
       TDIMREAD = TDIM
    ENDIF

    ! Return w/ success
    RC = 0

  END SUBROUTINE GET_TIDX
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TimeUnit_Check
!
! !DESCRIPTION: Makes a validity check of the passed unit string.
! Supported formats are "days since YYYY-MM-DD" (TIMETYPE=1) and
! "hours since YYYY-MM-DD HH:MM:SS" (TIMETYPE=2).
!\\
!\\
! The output argument TOFFSET gives the offset of the ncdf reference
! time relative to Geos-Chem reference time (in hours).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TIMEUNIT_CHECK( TIMEUNIT, TIMETYPE, TOFFSET, FILENAME, RC )
!
! !USES:
!
    USE CHARPAK_MOD, ONLY : TRANLC
!
! !INPUT PARAMETERS:
!
    ! Required
    CHARACTER(LEN=*), INTENT(IN  )  :: TIMEUNIT
    CHARACTER(LEN=*), INTENT(IN  )  :: FILENAME
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT) :: TIMETYPE
    REAL*8,           INTENT(  OUT) :: TOFFSET
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  18 Jan 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)    :: LOC, MSG
    CHARACTER(LEN=255)    :: MIRRUNIT
    INTEGER               :: STAT, L1, L2
    INTEGER               :: TTYPE
    INTEGER               :: YYYY, MM, DD, HH
    INTEGER               :: STRLEN

    !=================================================================
    ! TIMEUNIT_CHECK starts here
    !=================================================================

    ! Init
    LOC = 'TIMEUNIT_CHECK (ncdf_mod.F)'

    ! Check length of time unit string. This must be at least 21
    ! ("days since YYYY:MM:DD" is of length 21)
    STRLEN = LEN(TIMEUNIT)
    IF ( STRLEN < 21 ) THEN
       PRINT *, 'Time unit string too short: ' // TRIM(FILENAME)
       RC = -999; RETURN
    ENDIF

    ! ----------------------------------------------------------------------
    ! Determine time unit type
    ! ----------------------------------------------------------------------

    ! Mirror time unit and convert to lower case
    MIRRUNIT = TIMEUNIT
    CALL TRANLC( MIRRUNIT )

    ! Check for 'hours since'. If true, set TTYPE to 1 and set the
    ! begin of the effective date string to 12. Also check if the time
    ! string is at least of length 25, which is required for this
    ! unit.
    IF ( MIRRUNIT(1:11) == 'hours since' ) THEN
       TTYPE = 1
       L1    = 13
       IF ( STRLEN < 25 ) THEN
          PRINT *, 'Time unit string too short: ' // TRIM(FILENAME)
          RC = -999; RETURN
       ENDIF

    ! Check for 'days since'. If true, set TTYPE to 2 and set the
    ! begin of the effective date string to 11.
    ELSEIF ( MIRRUNIT(1:10) == 'days since' ) THEN
       TTYPE = 2
       L1    = 12
    ELSE
       ! Return w/ error
       PRINT *, 'Invalid time unit in', TRIM(FILENAME)
       RC = -999; RETURN
    ENDIF

    ! ----------------------------------------------------------------------
    ! Determine reference time/date
    ! Get the year, month, day and hour from the string
    ! '... since YYYY-MM-DD hh:mm:ss
    ! ----------------------------------------------------------------------

    ! Read reference year, i.e. first four integers
    L2 = L1 + 3
    READ( TIMEUNIT(L1:L2),'(i4)', IOSTAT=STAT ) YYYY
    IF ( STAT /= 0 ) THEN
       PRINT *, 'Invalid year in ', TRIM(TIMEUNIT), &
            ' in file'             , TRIM(FILENAME)
       RC = -999; RETURN
    ENDIF

    ! Read reference month. Typically, the month is represented by
    ! two characters, i.e. 1 is 01, etc.
    L1 = L2 + 2
    L2 = L1 + 1
    READ( TIMEUNIT(L1:L2), '(i2)', IOSTAT=STAT ) MM
    ! Also check for the case where the month is only one character:
    IF ( STAT /= 0 ) THEN
       L2 = L1
       READ( TIMEUNIT(L1:L2), '(i2)', IOSTAT=STAT ) MM
       IF ( STAT /= 0 ) THEN
          PRINT *, 'Invalid month in ', TRIM(TIMEUNIT), &
                   ' in file'         , TRIM(FILENAME)
          RC = -999; RETURN
       ENDIF
    ENDIF

    ! Reference day. Typically, the day is represented by two
    ! characters, i.e. 1 is 01, etc.
    L1 = L2 + 2
    L2 = L1 + 1
    READ( TIMEUNIT(L1:L2), '(i2)', IOSTAT=STAT ) DD
    ! Also check for the case where the day is only one character:
    IF ( STAT /= 0 ) THEN
       L2 = L1
       READ( TIMEUNIT(L1:L2), '(i2)', IOSTAT=STAT ) DD
       IF ( STAT /= 0 ) THEN
          PRINT *, 'Invalid day in ', TRIM(TIMEUNIT), &
                   ' in file'       , TRIM(FILENAME)
          RC = -999; RETURN
       ENDIF
    ENDIF

    ! Get reference hour only if 'hours since...'
    IF ( TTYPE == 1 ) THEN

       ! Reference hour
       L1 = L2 + 2
       L2 = L1 + 1
       READ( TIMEUNIT(L1:L2), '(i2)', IOSTAT=STAT ) HH
       IF ( STAT /= 0 ) THEN
          L2 = L1
          READ( TIMEUNIT(L1:L2), '(i2)', IOSTAT=STAT ) HH
          IF ( STAT /= 0 ) THEN
             PRINT *, 'Invalid hour in ', TRIM(TIMEUNIT), &
                      ' in file'            , TRIM(FILENAME)
             RC = -999; RETURN
          ENDIF
       ENDIF

    ELSE
       ! Set reference hour to 0
       HH = 0

    ENDIF

    ! Get reference tau relative to G-C reference time, i.e. the
    ! offset of the netCDF reference time to the G-C reference time.
    ! This is hours since G-C reftime.
    TOFFSET = GET_TAU0( MM, DD, YYYY, HH )

    ! Remove one day if TOFFSET is negative, i.e. if the netCDF
    ! reference time is older than G-C reference time. We have to do
    ! this because GET_TAU0 does count the last day in this case!
    IF ( TOFFSET < 0d0 ) THEN
       TOFFSET = TOFFSET + 24d0
    ENDIF

    ! Output argument
    TIMETYPE = TTYPE

    ! Return w/ success
    RC = 0

  END SUBROUTINE TIMEUNIT_CHECK
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Get_Grid_Edges_Sp
!
! !DESCRIPTION: Routine to get the longitude or latitude edges. If the edge
! cannot be read from the netCDF file, they are calculated from the provided
! grid midpoints. Use the axis input argument to discern between longitude
! (axis 1) and latitude (axis 2).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_GET_GRID_EDGES_SP( fID, AXIS, MID, NMID, EDGE, NEDGE, RC )
!
! !USES:
!
    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   ) :: fID             ! Ncdf File ID
    INTEGER,          INTENT(IN   ) :: AXIS            ! 1=lon, 2=lat
    REAL*4,           INTENT(IN   ) :: MID(NMID)       ! midpoints
    INTEGER,          INTENT(IN   ) :: NMID            ! # of midpoints
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*4,           POINTER       :: EDGE(:)         ! edges
    INTEGER,          INTENT(INOUT) :: NEDGE           ! # of edges
    INTEGER,          INTENT(INOUT) :: RC              ! Return code
!
! !REVISION HISTORY:
!  16 Jul 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! NC_GET_GRID_EDGES_SP begins here
    !======================================================================

    CALL NC_GET_GRID_EDGES_C( fID, AXIS, NMID, NEDGE, RC, &
                              MID4=MID,  EDGE4=EDGE )

  END SUBROUTINE NC_GET_GRID_EDGES_SP
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Get_Grid_Edges_Dp
!
! !DESCRIPTION: Routine to get the longitude or latitude edges. If the edge
! cannot be read from the netCDF file, they are calculated from the provided
! grid midpoints. Use the axis input argument to discern between longitude
! (axis 1) and latitude (axis 2).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_GET_GRID_EDGES_DP( fID, AXIS, MID, NMID, EDGE, NEDGE, RC )
!
! !USES:
!
    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   ) :: fID             ! Ncdf File ID
    INTEGER,          INTENT(IN   ) :: AXIS            ! 1=lon, 2=lat
    REAL*8,           INTENT(IN   ) :: MID(NMID)       ! midpoints
    INTEGER,          INTENT(IN   ) :: NMID            ! # of midpoints
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*8,           POINTER       :: EDGE(:)         ! edges
    INTEGER,          INTENT(INOUT) :: NEDGE           ! # of edges
    INTEGER,          INTENT(INOUT) :: RC              ! Return code
!
! !REVISION HISTORY:
!  16 Jul 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! NC_GET_GRID_EDGES_DP begins here
    !======================================================================

    CALL NC_GET_GRID_EDGES_C( fID, AXIS, NMID, NEDGE, RC, &
                              MID8=MID,  EDGE8=EDGE )

  END SUBROUTINE NC_GET_GRID_EDGES_DP
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Get_Grid_Edges_C
!
! !DESCRIPTION: Routine to get the longitude or latitude edges. If the edge
! cannot be read from the netCDF file, they are calculated from the provided
! grid midpoints. Use the axis input argument to discern between longitude
! (axis 1) and latitude (axis 2).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_GET_GRID_EDGES_C( fID, AXIS, NMID, NEDGE, RC, &
                                  MID4, MID8, EDGE4, EDGE8 )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   ) :: fID             ! Ncdf File ID
    INTEGER,          INTENT(IN   ) :: AXIS            ! 1=lon, 2=lat
    REAL*4, OPTIONAL, INTENT(IN   ) :: MID4(NMID)       ! midpoints
    REAL*8, OPTIONAL, INTENT(IN   ) :: MID8(NMID)       ! midpoints
    INTEGER,          INTENT(IN   ) :: NMID            ! # of midpoints
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*4, OPTIONAL, POINTER       :: EDGE4(:)         ! edges
    REAL*8, OPTIONAL, POINTER       :: EDGE8(:)         ! edges
    INTEGER,          INTENT(INOUT) :: NEDGE           ! # of edges
    INTEGER,          INTENT(INOUT) :: RC              ! Return code
!
! !REVISION HISTORY:
!  16 Jul 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL              :: PoleMid
    INTEGER              :: I, AS
    CHARACTER(LEN=255)   :: ncVar, ThisUnit

    !======================================================================
    ! NC_GET_GRID_EDGES_DP begins here
    !======================================================================

    ! Error trap: edge and mid must be same kind
    IF ( PRESENT(EDGE4) ) THEN
       IF ( .NOT. PRESENT(MID4) ) THEN
          PRINT *, 'If you provide EDGE4, you must also provide MID4'
          RC = -999
          RETURN
       ENDIF
    ELSEIF ( PRESENT(EDGE8) ) THEN
       IF ( .NOT. PRESENT(MID8) ) THEN
          PRINT *, 'If you provide EDGE8, you must also provide MID8'
          RC = -999
          RETURN
       ENDIF
    ELSE
       PRINT *, 'EDGE4 or EDGE8 must be given'
       RC = -999
       RETURN
    ENDIF

    ! Try to read edges from ncdf file
    IF ( AXIS == 1 ) THEN
       ncVar = 'lon_edge'
    ELSEIF ( AXIS == 2 ) THEN
       ncVar = 'lat_edge'
    ENDIF

    IF ( PRESENT(EDGE4) ) THEN
       CALL NC_READ_VAR( fID, TRIM(ncVar), nEdge, ThisUnit, Edge4, RC )
    ELSE
       CALL NC_READ_VAR( fID, TRIM(ncVar), nEdge, ThisUnit, Edge8, RC )
    ENDIF
    IF ( RC /= 0 ) RETURN

    ! Also try 'XXX_edges'
    IF ( nEdge == 0 ) THEN
       IF ( AXIS == 1 ) THEN
          ncVar = 'lon_edges'
       ELSEIF ( AXIS == 2 ) THEN
          ncVar = 'lat_edges'
       ENDIF
       IF ( PRESENT(EDGE4) ) THEN
          CALL NC_READ_VAR( fID, 'lon_edges', nEdge, ThisUnit, Edge4, RC )
       ELSE
          CALL NC_READ_VAR( fID, 'lon_edges', nEdge, ThisUnit, Edge8, RC )
       ENDIF
       IF ( RC /= 0 ) RETURN
    ENDIF

    ! Sanity check if edges are read from files: dimension must be nlon + 1!
    IF ( nEdge > 0 ) THEN
       IF ( nEdge /= (nMid + 1) ) THEN
          PRINT *, 'Edge has incorrect length!'
          RC = -999; RETURN
       ENDIF

    ! If not read from file, calculate from provided lon midpoints.
    ELSE

       nEdge = nMid + 1
       IF ( PRESENT(EDGE4) ) THEN
          IF ( ASSOCIATED ( Edge4 ) ) DEALLOCATE( Edge4 )
          ALLOCATE ( Edge4(nEdge), STAT=AS )
          IF ( AS /= 0 ) THEN
             PRINT *, 'Edge alloc. error in NC_GET_LON_EDGES (ncdf_mod.F90)'
             RC = -999; RETURN
          ENDIF
          Edge4 = 0.0
       ELSE
          IF ( ASSOCIATED ( Edge8 ) ) DEALLOCATE( Edge8 )
          ALLOCATE ( Edge8(nEdge), STAT=AS )
          IF ( AS /= 0 ) THEN
             PRINT *, 'Edge alloc. error in NC_GET_LON_EDGES (ncdf_mod.F90)'
             RC = -999; RETURN
          ENDIF
          Edge8 = 0.0d0
       ENDIF

       ! Get leftmost edge by extrapolating from first two midpoints.
       ! Error trap: for latitude axis, first edge must not be below -90!
       IF ( PRESENT(EDGE4) ) THEN
          Edge4(1) = Mid4(1) - ( (Mid4(2) - Mid4(1) ) / 2.0 )
          IF ( Edge4(1) < -90.0 .AND. AXIS == 2 ) Edge4(1) = -90.0
       ELSE
          Edge8(1) = Mid8(1) - ( (Mid8(2) - Mid8(1) ) / 2.0d0 )
          IF ( Edge8(1) < -90.0d0 .AND. AXIS == 2 ) Edge8(1) = -90.0d0
       ENDIF

       ! Calculate second edge. We need to catch the case where the first
       ! latitude mid-point is -90 (this is the case for GEOS-5 generic
       ! grids...). In that case, the second edge is put in the middle of
       ! the first two mid points (e.g. between -90 and -89). In all other
       ! case, we calculate it from the previously calculated left edge.
       IF ( PRESENT(EDGE4) ) THEN
          IF ( Mid4(1) == Edge4(1) ) THEN
             Edge4(2) = Mid4(1) + ( Mid4(2) - Mid4(1) ) / 2.0
             PoleMid  = .TRUE.
          ELSE
             Edge4(2) = Mid4(1) + Mid4(1) - Edge4(1)
             PoleMid  = .FALSE.
          ENDIF

          ! Sequentially calculate the right edge from the previously
          ! calculated left edge.
          DO I = 2, nMid
             Edge4(I+1) = Mid4(I) + Mid4(I) - Edge4(I)
          ENDDO

          ! Error check: max. lat edge must not exceed +90!
          IF ( Edge4(nMId+1) > 90.01 .AND. AXIS == 2 ) THEN
             IF ( PoleMid ) THEN
                Edge4(nMid+1) = 90.0
             ELSE
                PRINT *, 'Uppermost latitude edge above 90 deg north!'
                PRINT *, Edge4
                RC = -999; RETURN
             ENDIF
          ENDIF

       ! Real8
       ELSE
          IF ( Mid8(1) == Edge8(1) ) THEN
             Edge8(2) = Mid8(1) + ( Mid8(2) - Mid8(1) ) / 2.0d0
             PoleMid  = .TRUE.
          ELSE
             Edge8(2) = Mid8(1) + Mid8(1) - Edge8(1)
             PoleMid  = .FALSE.
          ENDIF

          ! Sequentially calculate the right edge from the previously
          ! calculated left edge.
          DO I = 2, nMid
             Edge8(I+1) = Mid8(I) + Mid8(I) - Edge8(I)
          ENDDO

          ! Error check: max. lat edge must not exceed +90!
          IF ( Edge8(nMId+1) > 90.01d0 .AND. AXIS == 2 ) THEN
             IF ( PoleMid ) THEN
                Edge8(nMid+1) = 90.0d0
             ELSE
                PRINT *, 'Uppermost latitude edge above 90 deg north!'
                PRINT *, Edge8
                RC = -999; RETURN
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    ! Return w/ success
    RC = 0

  END SUBROUTINE NC_GET_GRID_EDGES_C
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Get_Sigma_Levels_Sp
!
! !DESCRIPTION: Wrapper routine to get the sigma levels in single precision.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_GET_SIGMA_LEVELS_SP( fID,  ncFile, levName, lon1, lon2, lat1, &
                                     lat2, lev1,   lev2,    time, SigLev, dir, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   ) :: fID             ! Ncdf File ID
    CHARACTER(LEN=*), INTENT(IN   ) :: ncFile          ! ncFile
    CHARACTER(LEN=*), INTENT(IN   ) :: levName         ! variable name
    INTEGER,          INTENT(IN   ) :: lon1            ! lon lower bound
    INTEGER,          INTENT(IN   ) :: lon2            ! lon upper bound
    INTEGER,          INTENT(IN   ) :: lat1            ! lat lower bound
    INTEGER,          INTENT(IN   ) :: lat2            ! lat upper bound
    INTEGER,          INTENT(IN   ) :: lev1            ! lev lower bound
    INTEGER,          INTENT(IN   ) :: lev2            ! lev upper bound
    INTEGER,          INTENT(IN   ) :: time            ! time index
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*4,           POINTER       :: SigLev(:,:,:)   ! sigma levels
    INTEGER,          INTENT(INOUT) :: dir             ! axis direction (1=up;-1=down)
    INTEGER,          INTENT(INOUT) :: RC              ! Return code
!
! !REVISION HISTORY:
!  03 Oct 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

  CALL NC_GET_SIGMA_LEVELS_C( fID,  ncFile, levName, lon1, lon2, lat1, &
                              lat2, lev1,   lev2,    time, dir,  RC,   &
                              SigLev4=SigLev )

  END SUBROUTINE NC_GET_SIGMA_LEVELS_SP
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Get_Sigma_Levels_Dp
!
! !DESCRIPTION: Wrapper routine to get the sigma levels in double precision.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_GET_SIGMA_LEVELS_DP( fID,  ncFile, levName, lon1, lon2, lat1, &
                                     lat2, lev1,   lev2,    time, SigLev, dir, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   ) :: fID             ! Ncdf File ID
    CHARACTER(LEN=*), INTENT(IN   ) :: ncFile          ! ncFile
    CHARACTER(LEN=*), INTENT(IN   ) :: levName         ! variable name
    INTEGER,          INTENT(IN   ) :: lon1            ! lon lower bound
    INTEGER,          INTENT(IN   ) :: lon2            ! lon upper bound
    INTEGER,          INTENT(IN   ) :: lat1            ! lat lower bound
    INTEGER,          INTENT(IN   ) :: lat2            ! lat upper bound
    INTEGER,          INTENT(IN   ) :: lev1            ! lev lower bound
    INTEGER,          INTENT(IN   ) :: lev2            ! lev upper bound
    INTEGER,          INTENT(IN   ) :: time            ! time index
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*8,           POINTER       :: SigLev(:,:,:)   ! sigma levels
    INTEGER,          INTENT(INOUT) :: dir             ! axis direction (1=up;-1=down)
    INTEGER,          INTENT(INOUT) :: RC              ! Return code
!
! !REVISION HISTORY:
!  03 Oct 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

  CALL NC_GET_SIGMA_LEVELS_C( fID,  ncFile, levName, lon1, lon2, lat1, &
                              lat2, lev1,   lev2,    time, dir,  RC,   &
                              SigLev8=SigLev )

  END SUBROUTINE NC_GET_SIGMA_LEVELS_DP
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Get_Sigma_Levels_C
!
! !DESCRIPTION: Routine to get the sigma levels from the netCDF file
! within the given grid bounds and for the given time index. This routine
! attempts to construct the 3D sigma values from provided variable levName.
! The vertical coordinate system is determined based upon the variable
! attribute "standard\_name".
!\\
!\\
! For now, only hybrid sigma coordinate systems are supported, and the
! standard\_name attribute must follow CF conventions and be set to
! "atmosphere\_hybrid\_sigma\_pressure\_coordinate".
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_GET_SIGMA_LEVELS_C( fID,  ncFile, levName, lon1, lon2, lat1, &
                                    lat2, lev1,   lev2,    time, dir,  RC,   &
                                    SigLev4, SigLev8 )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   ) :: fID             ! Ncdf File ID
    CHARACTER(LEN=*), INTENT(IN   ) :: ncFile          ! ncFile
    CHARACTER(LEN=*), INTENT(IN   ) :: levName         ! variable name
    INTEGER,          INTENT(IN   ) :: lon1            ! lon lower bound
    INTEGER,          INTENT(IN   ) :: lon2            ! lon upper bound
    INTEGER,          INTENT(IN   ) :: lat1            ! lat lower bound
    INTEGER,          INTENT(IN   ) :: lat2            ! lat upper bound
    INTEGER,          INTENT(IN   ) :: lev1            ! lev lower bound
    INTEGER,          INTENT(IN   ) :: lev2            ! lev upper bound
    INTEGER,          INTENT(IN   ) :: time            ! time index
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT) :: dir             ! axis direction (1=up;-1=down)
    INTEGER,          INTENT(INOUT) :: RC              ! Return code
    REAL*4, OPTIONAL, POINTER       :: SigLev4(:,:,:)  ! sigma levels w/in
    REAL*8, OPTIONAL, POINTER       :: SigLev8(:,:,:)  ! specified boundaries
!
! !REVISION HISTORY:
!  03 Oct 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)   :: stdname
    CHARACTER(LEN=255)   :: a_name    ! netCDF attribute name
    INTEGER              :: a_type    ! netCDF attribute type
    LOGICAL              :: ok

    !======================================================================
    ! NC_GET_SIGMA_LEVELS begins here
    !======================================================================

    !------------------------------------------------------------------------
    ! Get level standard name. This attribute will be used to identify
    ! the coordinate system
    !------------------------------------------------------------------------
    ok = Ncdoes_Var_Exist( fID, TRIM(levName) )
    IF ( .NOT. ok ) THEN
       WRITE(*,*) 'Cannot find level variable ', TRIM(levName), ' in ', TRIM(ncFile), '!'
       RC = -999
       RETURN
    ENDIF

    ! Get standard name
    a_name = "standard_name"
    IF ( .NOT. NcDoes_Attr_Exist ( fID,          TRIM(levName),     &
                                   TRIM(a_Name), a_type         ) ) THEN
       WRITE(*,*) 'Cannot find level attribute ', TRIM(a_name), ' in variable ', &
                  TRIM(levName), ' - File: ', TRIM(ncFile), '!'
       RC = -999
       RETURN
    ENDIF
    CALL NcGet_Var_Attributes( fID, TRIM(levName), TRIM(a_name), stdname )

    !------------------------------------------------------------------------
    ! Call functions to calculate sigma levels depending on the coordinate
    ! system.
    !------------------------------------------------------------------------

    IF ( TRIM(stdname) == 'atmosphere_hybrid_sigma_pressure_coordinate' ) THEN

       IF ( PRESENT(SigLev4) ) THEN
          CALL NC_GET_SIG_FROM_HYBRID ( fID,  levName, lon1, lon2, lat1, lat2, &
                                        lev1, lev2,    time, dir, RC, SigLev4=SigLev4 )
       ELSEIF ( PRESENT(SigLev8) ) THEN
          CALL NC_GET_SIG_FROM_HYBRID ( fID,  levName, lon1, lon2, lat1, lat2, &
                                        lev1, lev2,    time, dir, RC, SigLev8=SigLev8 )
       ELSE
          WRITE(*,*) 'SigLev array is missing!'
          RC = -999
          RETURN
       ENDIF
       IF ( RC /= 0 ) RETURN

    ! NOTE: for now, only hybrid sigma coordinates are supported!
    ELSE
       WRITE(*,*) 'Invalid level standard name: ', TRIM(stdname), ' in ', TRIM(ncFile)
       RC = -999
       RETURN
    ENDIF

    ! Return w/ success
    RC = 0

  END SUBROUTINE NC_GET_SIGMA_LEVELS_C
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Get_Sig_From_Hybrid
!
! !DESCRIPTION: Calculates the sigma level field for a hybrid sigma coordinate
! system:
!
! sigma(i,j,l,t) = ( a(l) * p0 + b(l) * ps(i,j,t) ) / ps(i,j,t)
!
! or (p0=1):
!
! sigma(i,j,l,t) = ( ap(l) + b(l) * ps(i,j,t) ) / ps(i,j,t)
!
! where sigma are the sigma levels, ap and bp are the hybrid sigma coordinates,
! p0 is the constant reference pressure, and ps is the surface pressure. The
! variable names of ap, p0, bp, and ps are taken from level attribute
! `formula\_terms`.
!\\
!\\
! The direction of the vertical coordinate system is determined from attribute
! `positive` (up or down) or - if not found - from the b values, whereby it is
! assumed that the higher b value is found at the surface. The return argument
! dir is set to 1 for upward coordinates (level 1 is surface level) and -1 for
! downward coordinates (level 1 is top of atmosphere).
!\\
!\\
! !REMARKS:
! Example of valid netCDF meta-data: The attributes `standard\_name` and
! `formula\_terms` are required, as is the 3D surface pressure field.
!
! double lev(lev) ;\\
!        lev:standard_name = "atmosphere_hybrid_sigma_pressure_coordinate" ;\\
!        lev:units = "level" ;\\
!        lev:positive = "down" ;\\
!        lev:formula_terms = "ap: hyam b: hybm ps: PS" ;\\
! double hyam(nhym) ;\\
!        hyam:long_name = "hybrid A coefficient at layer midpoints" ;\\
!        hyam:units = "hPa" ;\\
! double hybm(nhym) ;\\
!        hybm:long_name = "hybrid B coefficient at layer midpoints" ;\\
!        hybm:units = "1" ;\\
! double time(time) ;\\
!        time:standard_name = "time" ;\\
!        time:units = "days since 2000-01-01 00:00:00" ;\\
!        time:calendar = "standard" ;\\
! double PS(time, lat, lon) ;\\
!        PS:long_name = "surface pressure" ;\\
!        PS:units = "hPa" ;\\
!
! !INTERFACE:
!
  SUBROUTINE NC_GET_SIG_FROM_HYBRID ( fID,  levName, lon1, lon2, lat1, lat2, &
                                      lev1, lev2,    time, dir,  RC,   sigLev4, sigLev8 )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   ) :: fID             ! Ncdf File ID
    CHARACTER(LEN=*), INTENT(IN   ) :: levName         ! variable name
    INTEGER,          INTENT(IN   ) :: lon1            ! lon lower bound
    INTEGER,          INTENT(IN   ) :: lon2            ! lon upper bound
    INTEGER,          INTENT(IN   ) :: lat1            ! lat lower bound
    INTEGER,          INTENT(IN   ) :: lat2            ! lat upper bound
    INTEGER,          INTENT(IN   ) :: lev1            ! lev lower bound
    INTEGER,          INTENT(IN   ) :: lev2            ! lev upper bound
    INTEGER,          INTENT(IN   ) :: time            ! time index
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*4, OPTIONAL, POINTER       :: SigLev4(:,:,:)  ! sigma levels w/in
    REAL*8, OPTIONAL, POINTER       :: SigLev8(:,:,:)  ! specified boundaries
    INTEGER,          INTENT(  OUT) :: dir             ! axis direction (1=up;-1=down)
    INTEGER,          INTENT(INOUT) :: RC              ! Return code
!
! !REVISION HISTORY:
!  03 Oct 2014 - C. Keller   - Initial version
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: I, J, l1, l2, AS
    INTEGER              :: nlev, nlat, nlon
    INTEGER              :: nlevs
    INTEGER              :: st1d(1), ct1d(1)
    LOGICAL              :: ok
    REAL*4, POINTER      :: a(:)
    REAL*4, POINTER      :: b(:)
    REAL*4, POINTER      :: ps(:,:,:,:)
    REAL*8               :: p0
    CHARACTER(LEN=255)   :: formula, ThisUnit
    CHARACTER(LEN=255)   :: aname, bname, psname, p0name
    CHARACTER(LEN=255)   :: a_name    ! netCDF attribute name
    INTEGER              :: a_type    ! netCDF attribute type

    !======================================================================
    ! NC_GET_SIG_FROM_HYBRID begins here
    !======================================================================

    ! Init
    p0 = -999.d0
    a  => NULL()
    b  => NULL()
    ps => NULL()

    ! Get desired grid dimensions.
    nlon = lon2 - lon1 + 1
    nlat = lat2 - lat1 + 1
    nlev = lev2 - lev1 + 1

    ! Get dimension length
    CALL Ncget_Dimlen ( fID, TRIM(LevName), nlevs )

    ! Sanity check
    IF ( nlevs < nlev ) THEN
       WRITE(*,*) TRIM(LevName), ' is only of length ', nlevs, ' - required is: ', nlev
       RC = -999
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Get formula and parse variable names (ap, bp, p0, ps)
    !------------------------------------------------------------------------

    ! Get formula
    a_name = "formula_terms"
    IF ( .NOT. NcDoes_Attr_Exist ( fID,          TRIM(levName),            &
                                   TRIM(a_name), a_type         ) ) THEN
       WRITE(*,*) 'Cannot find attribute ', TRIM(a_name), ' in variable ', &
                  TRIM(levName)
       RC = -999
       RETURN
    ENDIF
    CALL NcGet_Var_Attributes( fID, TRIM(levName), TRIM(a_name), formula )

    ! Get variable names
    !-------------------
    I = INDEX( formula, 'a:' )
    IF ( I > 0 ) THEN
       CALL GetVarFromFormula( formula, 'a:',  aname, RC )
       IF ( RC /= 0 ) RETURN
       CALL GetVarFromFormula( formula, 'p0:', p0name, RC )
       IF ( RC /= 0 ) RETURN
    ELSE
       CALL GetVarFromFormula( formula, 'ap:', aname, RC )
       IF ( RC /= 0 ) RETURN
       p0 = 1.0d0
    ENDIF
    IF ( RC /= 0 ) RETURN

    CALL GetVarFromFormula( formula, 'b:', bname, RC )
    IF ( RC /= 0 ) RETURN

    CALL GetVarFromFormula( formula, 'ps:', psname, RC )
    IF ( RC /= 0 ) RETURN

    !------------------------------------------------------------------------
    ! Read variables from file.
    !------------------------------------------------------------------------

    ALLOCATE ( a(nlevs), b(nlevs) )
    st1d = (/ 1     /)
    ct1d = (/ nlevs /)

    ! read a
    !-------
    IF ( .NOT. Ncdoes_Var_Exist( fID, TRIM(aname) ) ) THEN
       WRITE(*,*) 'Cannot find variable ', TRIM(aname), '!'
       RC = -999
       RETURN
    ENDIF
    CALL NcRd( a, fID, TRIM(aname), st1d, ct1d )

    ! eventually read p0
    !-------------------
    IF ( p0 < 0.0d0 ) THEN
       IF ( .NOT. Ncdoes_Var_Exist( fID, TRIM(p0name) ) ) THEN
          WRITE(*,*) 'Cannot find variable ', TRIM(p0name), '!'
          RC = -999
          RETURN
       ENDIF
    CALL NcRd( p0, fID, TRIM(p0name) )
    ENDIF

    ! read b
    !-------
    IF ( .NOT. Ncdoes_Var_Exist( fID, TRIM(bname) ) ) THEN
       WRITE(*,*) 'Cannot find variable ', TRIM(bname), '!'
       RC = -999
       RETURN
    ENDIF
    CALL NcRd( b, fID, TRIM(bname), st1d, ct1d )

    ! Read ps
    !--------
    CALL NC_READ_ARR( fID, TRIM(psname), lon1, lon2, lat1, &
                      lat2, 0, 0, time,  time, ps, VarUnit=thisUnit, RC=RC )
    IF ( RC /= 0 ) RETURN

    !------------------------------------------------------------------------
    ! Determine positive axis ('up' or 'down')
    ! Try to read it from the netCDF meta data (attribute `positive`). If not
    ! found, determine it from b values (b value at surface higher than at
    ! top of atmosphere).
    !------------------------------------------------------------------------
    a_name = "positive"
    IF ( NcDoes_Attr_Exist( fID, TRIM(levName), TRIM(a_name), a_type ) ) THEN
       CALL NcGet_Var_Attributes( fID, TRIM(levName), TRIM(a_name), formula )
       IF ( TRIM(formula) == 'up' ) THEN
          dir = 1
       ELSEIF ( TRIM(formula) == 'down' ) THEN
          dir = -1
       ELSE
          WRITE(*,*) 'level attribute `positive` must be `up` ', &
                     'or `down`, instead: ', TRIM(formula)
          RC = -999
          RETURN
       ENDIF

    ! determine direction from b values.
    ELSE

       IF ( b(1) > b(nlevs) ) THEN
          dir = 1
       ELSE
          dir = -1
       ENDIF
    ENDIF

    !------------------------------------------------------------------------
    ! Determine vertical indeces to be used. It is possible to calculate
    ! the pressure only for a given number of layers (as specified by input
    ! arguments lev1 and lev2). Assume those are always from bottom to top,
    ! i.e. counting `upwards`.
    !------------------------------------------------------------------------

    IF ( dir == -1 ) THEN
       l1 = nlevs - lev2 + 1
       l2 = nlevs - lev1 + 1
    ELSE
       l1 = lev1
       l2 = lev2
    ENDIF

    !------------------------------------------------------------------------
    ! Calculate sigma values at grid edges
    !------------------------------------------------------------------------

    IF ( PRESENT(SigLev4) ) THEN
       IF ( ASSOCIATED(SigLev4) ) DEALLOCATE(SigLev4)
       ALLOCATE(SigLev4(nlon,nlat,nlev),STAT=AS)
    ELSEIF ( PRESENT(SigLev8) ) THEN
       IF ( ASSOCIATED(SigLev8) ) DEALLOCATE(SigLev8)
       ALLOCATE(SigLev8(nlon,nlat,nlev),STAT=AS)
    ELSE
       WRITE(*,*) 'SigLev must be provided!'
       RC = -999
       RETURN
    ENDIF
    IF ( AS /= 0 ) THEN
       WRITE(*,*) 'Cannot allocate SigLev!'
       RC = -999
       RETURN
    ENDIF

    DO J=1,nlat
    DO I=1,nlon
       IF ( PRESENT(SigLev4) ) THEN
          SigLev4(i,j,:) = ( ( a(l1:l2) * p0 ) + ( b(l1:l2) * ps(i,j,1,1) ) ) &
                        / ps(i,j,1,1)
       ELSE
          SigLev8(i,j,:) = ( ( a(l1:l2) * p0 ) + ( b(l1:l2) * ps(i,j,1,1) ) ) &
                        / ps(i,j,1,1)
       ENDIF
    ENDDO
    ENDDO

    ! Cleanup
    IF ( ASSOCIATED(a ) ) DEALLOCATE(a )
    IF ( ASSOCIATED(b ) ) DEALLOCATE(b )
    IF ( ASSOCIATED(ps) ) DEALLOCATE(ps)

    ! Return w/ success
    RC = 0

  END SUBROUTINE NC_GET_SIG_FROM_HYBRID
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetVarFromFormula
!
! !DESCRIPTION: helper function to extract the variable name from a vertical
! coordinate formula.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetVarFromFormula ( formula, inname, outname, RC )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: formula
    CHARACTER(LEN=*), INTENT(IN   ) :: inname
!
! !INPUT/OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(  OUT) :: outname
    INTEGER,          INTENT(INOUT) :: RC              ! Return code
!
! !REVISION HISTORY:
!  03 Oct 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: I, J, IDX, LN

    !======================================================================
    ! GetVarFromFormula begins here
    !======================================================================

    ! maximum length
    LN = LEN(TRIM(formula))

    ! Get start index of string
    !--------------------------
    I = INDEX( TRIM(formula), TRIM(inname) )
    IF ( I <= 0 ) THEN
       WRITE(*,*) 'Cannot extract ', TRIM(inname), ' from ', TRIM(formula)
       RC = -999
       RETURN
    ENDIF

    ! The variable name follows the formula string plus one space!
    I = I + LEN(inname) + 1

    outname = ''
    IDX = 1
    DO J = I, LN
       IF ( formula(J:J) == ' ' ) EXIT
       outname(IDX:IDX) = formula(J:J)
       IDX = IDX + 1
    ENDDO

    ! Return w/ success
    RC = 0

  END SUBROUTINE GetVarFromFormula
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Write_3d
!
! !DESCRIPTION: Routine to write time slices of 2D fields into netCDF.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_WRITE_3D( ncFile,  I,  J,    T,  N,   lon, lat, &
                          time,    timeUnit, ncVars,  ncUnits,  &
                          ncLongs, ncShorts, ncArrays            )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)  :: ncFile             ! file path+name
    INTEGER,          INTENT(IN)  :: I                  ! # of lons
    INTEGER,          INTENT(IN)  :: J                  ! # of lats
    INTEGER,          INTENT(IN)  :: T                  ! # of time slices
    INTEGER,          INTENT(IN)  :: N                  ! # of vars
    REAL*4,           INTENT(IN)  :: lon(I)             ! longitude
    REAL*4,           INTENT(IN)  :: lat(J)             ! latitude
    REAL*4,           INTENT(IN)  :: time(T)            ! time
    CHARACTER(LEN=*), INTENT(IN)  :: timeUnit           ! time unit
    CHARACTER(LEN=*), INTENT(IN)  :: ncVars(N)          ! nc variables
    CHARACTER(LEN=*), INTENT(IN)  :: ncUnits(N)         ! var units
    CHARACTER(LEN=*), INTENT(IN)  :: ncLongs(N)         ! var long names
    CHARACTER(LEN=*), INTENT(IN)  :: ncShorts(N)        ! var short names
    REAL*4, TARGET,   INTENT(IN)  :: ncArrays(I,J,T,N)  ! var arrays
!
! !REMARKS:
!  Created with the ncCodeRead script of the NcdfUtilities package,
!  with subsequent hand-editing.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: fId, II
    REAL*4, POINTER    :: tmpArr(:,:,:) => NULL()

    !======================================================================
    ! NC_WRITE_3D begins here
    !======================================================================

    CALL NC_DEFINE(ncFile=ncFile, nLon=I,            nLat=J,         &
                   nTime=T,       timeUnit=timeUnit, ncVars=ncVars,  &
                   ncUnits=ncUnits,ncLongs=ncLongs,ncShorts=ncShorts,&
                   fId=fId )

    CALL NC_WRITE_DIMS( fID=fId, lon=lon, lat=lat, time=time )

    DO II = 1, N
       tmpArr => ncArrays(:,:,:,II)
       CALL NC_WRITE_DATA_3D ( fId, ncVars(II), tmpArr )
       tmpArr => NULL()
    ENDDO

    CALL NcCl( fId )

  END SUBROUTINE NC_WRITE_3D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Write_4d
!
! !DESCRIPTION: Routine to write time slices of 3D fields into netCDF.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_WRITE_4D (ncFile,  I, J, L, T, N, lon, lat, lev, &
                          time,    timeUnit, ncVars,  ncUnits,   &
                          ncLongs, ncShorts, ncArrays             )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)  :: ncFile   ! file path+name
    INTEGER,          INTENT(IN)  :: I        ! # of lons
    INTEGER,          INTENT(IN)  :: J        ! # of lats
    INTEGER,          INTENT(IN)  :: L        ! # of levs
    INTEGER,          INTENT(IN)  :: T        ! # of time slices
    INTEGER,          INTENT(IN)  :: N        ! # of vars
    REAL*4,           INTENT(IN)  :: lon(:)   ! longitude
    REAL*4,           INTENT(IN)  :: lat(:)   ! latitude
    REAL*4,           INTENT(IN)  :: lev(:)   ! levels
    REAL*4,           INTENT(IN)  :: time(:)  ! time
    CHARACTER(LEN=*), INTENT(IN)  :: timeUnit ! time unit
    CHARACTER(LEN=*), INTENT(IN)  :: ncVars(:)    ! nc variables
    CHARACTER(LEN=*), INTENT(IN)  :: ncUnits(:)   ! var units
    CHARACTER(LEN=*), INTENT(IN)  :: ncLongs(:)   ! var long names
    CHARACTER(LEN=*), INTENT(IN)  :: ncShorts(:)  ! var short names
    REAL*4, TARGET,   INTENT(IN)  :: ncArrays(:,:,:,:,:)  ! var arrays
!
! !REMARKS:
!  Created with the ncCodeRead script of the NcdfUtilities package,
!  with subsequent hand-editing.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: II, fID
    REAL*4, POINTER    :: tmpArr(:,:,:,:) => NULL()

    !======================================================================
    ! NC_WRITE begins here
    !======================================================================

    CALL NC_DEFINE(ncFile=ncFile, nLon=I, nLat=J, nLev=L,            &
                   nTime=T,  timeUnit=timeUnit, ncVars=ncVars,       &
                   ncUnits=ncUnits,ncLongs=ncLongs,ncShorts=ncShorts,&
                   fId=fId )

    CALL NC_WRITE_DIMS( fID=fId, lon=lon, lat=lat, time=time, lev=lev)

    DO II = 1, size(ncVars)
       tmpArr => ncArrays(:,:,:,:,II)
       CALL NC_WRITE_DATA_4D ( fId, ncVars(II), tmpArr )
       tmpArr => NULL()
    ENDDO

    CALL NcCl( fId )

  END SUBROUTINE NC_WRITE_4D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Define
!
! !DESCRIPTION: Routine to define the variables and attributes of a netCDF
!  file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_DEFINE ( ncFile,  nLon,    nLat,    nLev,    nTime,&
                         timeUnit, ncVars,  ncUnits, ncLongs, ncShorts, fId )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),   INTENT(IN   )  :: ncFile      ! ncdf file path + name
    INTEGER,            INTENT(IN   )  :: nLon        ! # of lons
    INTEGER,            INTENT(IN   )  :: nLat        ! # of lats
    INTEGER, OPTIONAL,  INTENT(IN   )  :: nLev        ! # of levels
    INTEGER,            INTENT(IN   )  :: nTime       ! # of time stamps
    CHARACTER(LEN=*),   INTENT(IN   )  :: timeUnit    ! time unit
    CHARACTER(LEN=*),   INTENT(IN   )  :: ncVars(:)   ! ncdf variables
    CHARACTER(LEN=*),   INTENT(IN   )  :: ncUnits(:)  ! var units
    CHARACTER(LEN=*),   INTENT(IN   )  :: ncLongs(:)  ! var long names
    CHARACTER(LEN=*),   INTENT(IN   )  :: ncShorts(:) ! var short names
!
! !OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(  OUT)  :: fId      ! netCDF file ID
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  10 May 2017 - R. Yantosca - Don't manually increment vId, it's returned
!                              as an output from NCDEF_VARIABLE
!  18 May 2018 - C. Holmes   - Define time as an unlimited dimension
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Declare netCDF variable ID and fill mode
    INTEGER            :: vId
    INTEGER            :: omode

    ! Variables for netCDF dimensions
    INTEGER            :: id_lon
    INTEGER            :: id_lat
    INTEGER            :: id_time
    INTEGER            :: id_lev

    ! Character strings
    CHARACTER(LEN=255) :: v_name      ! netCDF variable name
    CHARACTER(LEN=255) :: a_name      ! netCDF attribute name
    CHARACTER(LEN=255) :: a_val       ! netCDF attribute value
    CHARACTER(LEN=3  ) :: idstr       ! tracer ID string

    ! Arrays for netCDF dimension IDs
    INTEGER            :: var1d(1)    ! For 1D arrays
    INTEGER            :: var3d(3)    ! For 3D arrays
    INTEGER            :: var4d(4)    ! For 4D arrays

    ! Other variables
    INTEGER            :: I

    !=================================================================
    ! %%%%% NETCDF DEFINITION SECTION %%%%%
    !=================================================================

    ! Initialize the variable ID counter
    vId = 0

    ! Open filename
    CALL NcCr_Wr( fId, TRIM(ncFile) )

    ! Turn filling off
    CALL NcSetFill( fId, NF_NOFILL, omode )

    !--------------------------------
    ! GLOBAL ATTRIBUTES
    !--------------------------------

    ! Define the title global attribute
    a_name = "Title"
    a_val  = "Field generated by ncdf_util.F"
    CALL NcDef_Glob_Attributes( fId, TRIM(a_name), TRIM(a_val) )

    ! Define the history global attribute
    a_name = "History"
    a_val  = "Initial version"
    CALL NcDef_Glob_Attributes( fId, TRIM(a_name), TRIM(a_val) )

    ! Define the conventions global attribute
    a_name = "Conventions"
    a_val  = "COARDS"
    CALL NcDef_Glob_Attributes( fId, TRIM(a_name), TRIM(a_val) )

    ! Define the format global attribute
    a_name = "Format"
    a_val  = "netCDF-3"
    CALL NcDef_Glob_Attributes( fId, TRIM(a_name), TRIM(a_val) )

    !--------------------------------
    ! DIMENSIONS
    !--------------------------------

    ! Define lon dimension
    v_name = "lon"
    CALL NcDef_Dimension( fId, TRIM(v_name), nlon, id_lon )

    ! Define lat dimension
    v_name = "lat"
    CALL NcDef_Dimension( fId, TRIM(v_name), nlat, id_lat )

    ! Define lev dimension
    IF ( PRESENT(nlev) ) THEN
       v_name = "lev"
       CALL NcDef_Dimension( fId, TRIM(v_name), nlev, id_lev )
    ENDIF

    ! Define time dimension
    v_name = "time"
    CALL NcDef_Dimension( fId, TRIM(v_name), ntime, id_time, unlimited=.true. )

    !--------------------------------
    ! VARIABLE: lon
    !--------------------------------

    ! Define the "lon" variable
    v_name = "lon"
    var1d = (/ id_lon /)
    CALL NcDef_Variable( fId, TRIM(v_name), NF_FLOAT, 1, var1d, vId )

    ! Define the "lon:long_name" attribute
    a_name = "long_name"
    a_val  = "Longitude"
    CALL NcDef_Var_Attributes( fId, vId, TRIM(a_name), TRIM(a_val) )

    ! Define the "lon:units" attribute
    a_name = "units"
    a_val  = "degrees_east"
    CALL NcDef_Var_Attributes( fId, vId, TRIM(a_name), TRIM(a_val) )

    !--------------------------------
    ! VARIABLE: lat
    !--------------------------------

    ! Define the "lat" variable
    v_name = "lat"
    var1d = (/ id_lat /)
    CALL NcDef_Variable( fId, TRIM(v_name), NF_FLOAT, 1, var1d, vId )

    ! Define the "lat:long_name" attribute
    a_name = "long_name"
    a_val  = "Latitude"
    CALL NcDef_Var_Attributes( fId, vId, TRIM(a_name), TRIM(a_val) )

    ! Define the "lat:units" attribute
    a_name = "units"
    a_val  = "degrees_north"
    CALL NcDef_Var_Attributes( fId, vId, TRIM(a_name), TRIM(a_val) )

    !--------------------------------
    ! VARIABLE: lev
    !--------------------------------

    IF ( PRESENT(nlev) ) THEN

       ! Define the "levels" variable
       v_name = "lev"
       var1d = (/ id_lev /)
       CALL NcDef_Variable( fId, TRIM(v_name), NF_INT, 1, var1d, vId )

       ! Define the "time:long_name" attribute
       a_name = "long_name"
       a_val  = "Levels"
       CALL NcDef_Var_Attributes( fId, vId, TRIM(a_name), TRIM(a_val))

       ! Define the "time:units" attribute
       a_name = "units"
       a_val  = "unitless"
       CALL NcDef_Var_Attributes( fId, vId, TRIM(a_name), TRIM(a_val))
    ENDIF

    !--------------------------------
    ! VARIABLE: time
    !--------------------------------

    ! Define the "time" variable
    v_name = "time"
    var1d = (/ id_time /)
    CALL NcDef_Variable( fId, TRIM(v_name), NF_INT, 1, var1d, vId )

    ! Define the "time:long_name" attribute
    a_name = "long_name"
    a_val  = "Time"
    CALL NcDef_Var_Attributes( fId, vId, TRIM(a_name), TRIM(a_val) )

    ! Define the "time:units" attribute
    a_name = "units"
    a_val  = trim(timeUnit)
    CALL NcDef_Var_Attributes( fId, vId, TRIM(a_name), TRIM(a_val) )

    !--------------------------------
    ! Define variables
    !--------------------------------

    DO I = 1, SIZE(ncVars)

       v_name = TRIM(ncVars(I))
       IF ( PRESENT(nlev) ) THEN
          var4d = (/ id_lon, id_lat, id_lev, id_time /)
          CALL NcDef_Variable(fId,TRIM(v_name),NF_DOUBLE,4,var4d,vId)
       ELSE
          var3d = (/ id_lon, id_lat, id_time /)
          CALL NcDef_Variable(fId,TRIM(v_name),NF_DOUBLE,3,var3d,vId)
       ENDIF

       ! Define the long_name attribute
       a_name = "long_name"
       a_val  = TRIM(ncLongs(I))
       CALL NcDef_Var_Attributes(fId, vId, TRIM(a_name), TRIM(a_val) )

       ! Define the short_name attribute
       a_name = "short_name"
       a_val  = TRIM(ncShorts(I))
       CALL NcDef_Var_Attributes(fId, vId, TRIM(a_name), TRIM(a_val) )

       ! Define the units attribute
       a_name = "units"
       a_val  = TRIM(ncUnits(I))
       CALL NcDef_Var_Attributes(fId, vId, TRIM(a_name), TRIM(a_val) )
    ENDDO

    !=================================================================
    ! %%%%% END OF NETCDF DEFINITION SECTION %%%%%
    !=================================================================
    CALL NcEnd_Def( fId )

  END SUBROUTINE NC_DEFINE
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Write_Dims
!
! !DESCRIPTION: Routine to write dimension arrays to a netCDF file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_WRITE_DIMS( fID, lon, lat, time, lev )
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: fId
!
! !INPUT PARAMETERS:
!
    REAL*4,           INTENT(IN   ) :: lon(:)
    REAL*4,           INTENT(IN   ) :: lat(:)
    REAL*4,           INTENT(IN   ) :: time(:)
    REAL*4, OPTIONAL, INTENT(IN   ) :: lev(:)
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!  13 Jun 2014 - R. Yantosca - Avoid array temporaries
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Character strings
    CHARACTER(LEN=255) :: v_name             ! netCDF variable name

    ! Arrays for netCDF start and count values
    INTEGER            :: st1d(1), ct1d(1)   ! For 1D arrays
    INTEGER            :: v_size

    !=================================================================
    ! Define lon/lat
    !=================================================================

    ! Write lon to netCDF file
    v_name = "lon"
    v_size = size( lon, 1 )
    st1d   = (/ 1      /)
    ct1d   = (/ v_size /)
    CALL NcWr( lon, fId, TRIM(v_name), st1d, ct1d )

    ! Write lat to netCDF file
    v_name = "lat"
    v_size = size( lat, 1 )
    st1d   = (/ 1      /)
    ct1d   = (/ v_size /)
    CALL NcWr( lat, fId, TRIM(v_name), st1d, ct1d )

    ! Write lev to netCDF file
    IF ( PRESENT(lev) ) THEN
       v_name = "lev"
       v_size = size( lev, 1 )
       st1d   = (/ 1      /)
       ct1d   = (/ v_size /)
       CALL NcWr( lev, fId, TRIM(v_name), st1d, ct1d )
    ENDIF

    ! Write passed time integer to netCDF file
    v_name = "time"
    v_size = size( time, 1 )
    st1d   = (/ 1      /)
    ct1d   = (/ v_size /)
    CALL NcWr( time, fId, TRIM(v_name), st1d, ct1d )

  END SUBROUTINE NC_WRITE_DIMS
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Nrite_Data_3d
!
! !DESCRIPTION: Routine to write a 3-D array to a netCDF file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_WRITE_DATA_3D ( fID, ncVar, Array )
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: fId
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: ncVar
    REAL*4,           POINTER       :: Array(:,:,:)
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays for netCDF start and count values
    INTEGER            :: st3d(3), ct3d(3)   ! For 3D arrays

    !=================================================================
    ! Write data to netCDF file
    !=================================================================

    st3d = (/ 1, 1, 1 /)
    ct3d = (/ size(array,1), size(array,2), size(array,3) /)
    CALL NcWr( ARRAY, fId, TRIM(ncVar), st3d, ct3d )

  END SUBROUTINE NC_WRITE_DATA_3D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Write_Data_4d
!
! !DESCRIPTION: Routine to write a 4-D array to a netCDF file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_WRITE_DATA_4D ( fID, ncVar, Array )
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: fId
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: ncVar
    REAL*4,           POINTER       :: Array(:,:,:,:)
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  30 Jan 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays for netCDF start and count values
    INTEGER            :: st4d(4), ct4d(4)   ! For 4D arrays

    !=================================================================
    ! Write data to netCDF file
    !=================================================================

    st4d   = (/ 1, 1, 1, 1 /)
    ct4d   = (/ size(array,1), size(array,2), &
                  size(array,3), size(array,4)   /)
    CALL NcWr( ARRAY, fId, TRIM(ncVar), st4d, ct4d )

  END SUBROUTINE NC_WRITE_DATA_4D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Create
!
! !DESCRIPTION: Creates a new netCDF file and defines several global
!  attributes.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Nc_Create( NcFile,      Title,          nLon,                   &
                        nLat,        nLev,           nTime,                  &
                        fId,         lonID,          latId,                  &
                        levId,       timeId,         VarCt,                  &
                        Create_NC4,  KeepDefMode,    NcFormat,               &
                        Conventions, History,        ProdDateTime,           &
                        Reference,   Contact,        nIlev,                  &
                        iLevId,      StartTimeStamp, EndTimeStamp           )
!
! !INPUT PARAMETERS:
!
    ! Required arguments
    CHARACTER(LEN=*), INTENT(IN   )  :: ncFile         ! ncdf file path + name
    CHARACTER(LEN=*), INTENT(IN   )  :: title          ! ncdf file title
    INTEGER,          INTENT(IN   )  :: nLon           ! # of lons
    INTEGER,          INTENT(IN   )  :: nLat           ! # of lats
    INTEGER,          INTENT(IN   )  :: nLev           ! # of level midpoints
    INTEGER,          INTENT(IN   )  :: nTime          ! # of times
    INTEGER,          OPTIONAL       :: nILev          ! # of level interfaces

    ! Optional arguments (mostly global attributes)
    LOGICAL,          OPTIONAL       :: Create_Nc4     ! Save as netCDF-4
    LOGICAL,          OPTIONAL       :: KeepDefMode    ! If = T, then don't
                                                       !  exit define mode
    CHARACTER(LEN=*), OPTIONAL       :: NcFormat       ! e.g. netCDF-4
    CHARACTER(LEN=*), OPTIONAL       :: Conventions    ! e.g. COARDS, CF, etc.
    CHARACTER(LEN=*), OPTIONAL       :: History        ! History glob attribute
    CHARACTER(LEN=*), OPTIONAL       :: ProdDateTime   ! Time/date of production
    CHARACTER(LEN=*), OPTIONAL       :: Reference      ! Reference string
    CHARACTER(LEN=*), OPTIONAL       :: Contact        ! People to contact
    CHARACTER(LEN=*), OPTIONAL       :: StartTimeStamp ! Timestamps at start
    CHARACTER(LEN=*), OPTIONAL       :: EndTimeStamp   !  and end of simulation
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)  :: fId            ! file id
    INTEGER,          INTENT(  OUT)  :: lonId          ! lon  dimension id
    INTEGER,          INTENT(  OUT)  :: latId          ! lat  dimension id
    INTEGER,          INTENT(  OUT)  :: levId          ! lev  dimension id
    INTEGER,          INTENT(  OUT)  :: timeId         ! time dimension id
    INTEGER,          INTENT(  OUT)  :: VarCt          ! variable counter
    INTEGER,          OPTIONAL       :: ilevId         ! ilev dimension id
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  11 Jan 2016 - R. Yantosca - Added optional CREATE_NC4 to save as netCDF-4
!  14 Jan 2016 - E. Lundgren - Pass title string for netcdf metadata
!  08 Aug 2017 - R. Yantosca - Add more optional arguments (mostly global atts)
!  08 Aug 2017 - R. Yantosca - Now define in dims in order: time,lev,lat,lon
!  08 Aug 2017 - R. Yantosca - Add optional KeepDefMode argument so that we can
!                              stay in netCDF define mode upon leaving this
!                              routine (i.e. to define variables afterwards)
!  24 Aug 2017 - R. Yantosca - Added nIlev and iLevId variables so that we can
!                               create the iLev dimension (level interfaces)
!  24 Jan 2018 - R. Yantosca - Add update frequency as an optional global attr
!  31 Jan 2018 - R. Yantosca - Add StartTimeStamp, EndTimeStamp arguments
!  18 May 2018 - C. Holmes   - Define time as an unlimited dimension
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: omode
    LOGICAL            :: Save_As_Nc4
    LOGICAL            :: QuitDefMode

    ! Strings
    CHARACTER(LEN=255) :: ThisHistory
    CHARACTER(LEN=255) :: ThisNcFormat
    CHARACTER(LEN=255) :: ThisConv
    CHARACTER(LEN=255) :: ThisPdt
    CHARACTER(LEN=255) :: ThisReference
    CHARACTER(LEN=255) :: ThisContact
    CHARACTER(LEN=255) :: ThisStartTimeStamp
    CHARACTER(LEN=255) :: ThisEndTimeStamp

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Create file as NetCDF4?
    IF ( PRESENT( Create_Nc4 ) ) THEN
       Save_As_Nc4 = Create_Nc4
    ELSE
       Save_As_Nc4 = .FALSE.
    ENDIF

    ! Should we exit netCDF define mode before leaving this routine?
    IF ( PRESENT( KeepDefMode ) ) THEN
       QuitDefMode = ( .not. KeepDefMode )
    ELSE
       QuitDefMode = .TRUE.
    ENDIF

    ! History global attribute
    IF ( PRESENT( History ) ) THEN
       ThisHistory = TRIM( History )
    ELSE
       ThisHistory = 'Created by routine NC_CREATE (in ncdf_mod.F90)'
    ENDIF

    ! NetCDF format global attribute
    IF ( PRESENT( NcFormat ) ) Then
       ThisNcFormat = NcFormat
    ELSE
       IF ( Save_As_Nc4 ) THEN
          ThisNcFormat = 'NetCDF-4'
       ELSE
          ThisNcFormat = 'NetCDF-3'
       ENDIF
    ENDIF

    ! Conventions global attribute (assume COARDS)
    IF ( PRESENT( Conventions ) ) THEN
       ThisConv = TRIM( Conventions )
    ELSE
       ThisConv = 'COARDS'
    ENDIF

    ! Conventions global attribute (assume COARDS)
    IF ( PRESENT( ProdDateTime ) ) THEN
       ThisPdt= TRIM( ProdDateTime )
    ENDIF

    ! Conventions global attribute (assume COARDS)
    IF ( PRESENT( Reference ) ) THEN
       ThisReference = TRIM( Reference )
    ELSE
       ThisReference = ''
    ENDIF

    ! Contact
    IF ( PRESENT( Contact ) ) THEN
       ThisContact = TRIM( Contact )
    ELSE
       ThisContact = ''
    ENDIF

    ! Starting date and time of the simulation
    IF ( PRESENT( StartTimeStamp ) ) THEN
       ThisStartTimeStamp = TRIM( StartTimeStamp )
    ELSE
       ThisStartTimeStamp = ''
    ENDIF

    ! Ending date and time of the simulation
    IF ( PRESENT( EndTimeStamp ) ) THEN
       ThisEndTimeStamp = TRIM( EndTimeStamp )
    ELSE
       ThisEndTimeStamp = ''
    ENDIF

    !=======================================================================
    ! Open the file
    !=======================================================================

    ! Open filename.  Save file in netCDF-4 format if requested by user.
    CALL NcCr_Wr( fId, TRIM( ncFile ), Save_As_Nc4 )

    ! Turn filling off
    CALL NcSetFill( fId, NF_NOFILL, omode )

    !=======================================================================
    ! Set global attributes
    !=======================================================================

    ! These attributes are required for COARDS or CF conventions
    CALL NcDef_Glob_Attributes(  fId, 'title',        TRIM( Title         ) )
    CALL NcDef_Glob_Attributes(  fId, 'history',      TRIM( ThisHistory   ) )
    CALL NcDef_Glob_Attributes(  fId, 'format',       TRIM( ThisNcFormat  ) )
    CALL NcDef_Glob_Attributes(  fId, 'conventions',  TRIM( ThisConv      ) )

    ! These attributes are optional
    IF ( PRESENT( ProdDateTime ) ) THEN
     CALL NcDef_Glob_Attributes( fId, 'ProdDateTime', TRIM( ThisPdt       ) )
    ENDIF

    IF ( PRESENT( Reference ) ) THEN
     CALL NcDef_Glob_Attributes( fId, 'reference',    TRIM( ThisReference ) )
    ENDIF

    IF ( PRESENT( Contact ) ) THEN
     CALL NcDef_Glob_Attributes( fId, 'contact',      TRIM( ThisContact   ) )
    ENDIF

    IF ( PRESENT( StartTimeStamp ) ) THEN
     CALL NcDef_Glob_Attributes( fId, 'simulation_start_date_and_time',      &
                                       TRIM( ThisStartTimeStamp   )         )
    ENDIF

    IF ( PRESENT( EndTimeStamp ) ) THEN
     CALL NcDef_Glob_Attributes( fId, 'simulation_end_date_and_time',        &
                                       TRIM( ThisEndTimeStamp )             )
    ENDIF

    !=======================================================================
    ! Set dimensions
    !=======================================================================

    ! Time
    CALL NcDef_Dimension( fId, 'time', nTime, TimeId, unlimited=.true. )

    ! Level midpoints
    IF ( nLev > 0 ) THEN
       CALL NcDef_Dimension( fId, 'lev',  nLev,  levId  )
    ELSE
       levId = -1
    ENDIF

    ! Optional ILev dimension: level interfaces
    IF ( PRESENT( nIlev ) .and. PRESENT( iLevId ) ) THEN
       IF ( nILev > 0 ) THEN
          CALL NcDef_Dimension( fId, 'ilev', nIlev, iLevId )
       ELSE
          iLevId = -1
       ENDIF
    ENDIF

    ! Lat and lon
    CALL NcDef_Dimension( fId, 'lat',  nLat,  latId  )
    CALL NcDef_Dimension( fId, 'lon',  nLon,  lonId  )

    ! Close definition section
    IF ( QuitDefMode ) THEN
       CALL NcEnd_Def( fId )
    ENDIF

    ! Initialize variable counter
    VarCt = -1

  END SUBROUTINE Nc_Create
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Def
!
! !DESCRIPTION: Defines a new netCDF variable along with its attributes.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_Var_Def( fId,       lonId,        latId,        levId,       &
                         TimeId,    VarName,      VarLongName,  VarUnit,     &
                         DataType,  VarCt,        DefMode,      Compress,    &
                         AddOffset, MissingValue, ScaleFactor,  Calendar,    &
                         Axis,      StandardName, FormulaTerms, AvgMethod,   &
                         Positive,  iLevId,       nUpdates                  )
!
! !INPUT PARAMETERS:
!
    ! Required inputs
    INTEGER,          INTENT(IN   ) :: fId          ! file ID
    INTEGER,          INTENT(IN   ) :: lonId        ! ID of lon      (X) dim
    INTEGER,          INTENT(IN   ) :: latId        ! ID of lat      (Y) dim
    INTEGER,          INTENT(IN   ) :: levId        ! ID of lev ctr  (Z) dim
    INTEGER,          OPTIONAL      :: iLevId       ! ID of lev edge (I) dim
    INTEGER,          INTENT(IN   ) :: TimeId       ! ID of time     (T) dim
    CHARACTER(LEN=*), INTENT(IN   ) :: VarName      ! Variable name
    CHARACTER(LEN=*), INTENT(IN   ) :: VarLongName  ! Long name description
    CHARACTER(LEN=*), INTENT(IN   ) :: VarUnit      ! Units
    INTEGER,          INTENT(IN   ) :: DataType     ! 1=Int, 4=float, 8=double

    ! Optional inputs
    LOGICAL,          OPTIONAL      :: DefMode      ! Toggles define mode
    LOGICAL,          OPTIONAL      :: Compress     ! Toggles compression
    REAL*4,           OPTIONAL      :: AddOffset    ! Add offset attribute
    REAL*4,           OPTIONAL      :: MissingValue ! Missing value attribute
    REAL*4,           OPTIONAL      :: ScaleFactor  ! Scale factor attribute
    CHARACTER(LEN=*), OPTIONAL      :: Calendar     ! Calendar for time var
    CHARACTER(LEN=*), OPTIONAL      :: Axis         ! Axis for index vars
    CHARACTER(LEN=*), OPTIONAL      :: StandardName ! Standard name attribute
    CHARACTER(LEN=*), OPTIONAL      :: FormulaTerms ! Formula for vert coords
    CHARACTER(LEN=*), OPTIONAL      :: AvgMethod    ! Averaging method
    CHARACTER(LEN=*), OPTIONAL      :: Positive     ! Positive dir (up or down)
    REAL*4,           OPTIONAL      :: nUpdates     ! # of updates (for time-
                                                    !  averaged fields only)
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: VarCt        ! variable counter
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  21 Jan 2017 - C. Holmes   - Added optional DefMode argument to avoid
!                              excessive switching between define & data modes
!  18 Feb 2017 - C. Holmes   - Enable netCDF-4 compression
!  08 Aug 2017 - R. Yantosca - Add more optional arguments for variable atts
!  24 Aug 2017 - R. Yantosca - Added StandardName, FormulaTerms arguments
!  24 Aug 2017 - R. Yantosca - Added optional Ilev dimension so that we can
!                               define variables on level interfaces
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays
    INTEGER, ALLOCATABLE :: VarDims(:)

    ! Scalars
    INTEGER              :: nDim,     Pos
    INTEGER              :: NF_TYPE,  tmpIlevId
    LOGICAL              :: isDefMode

    ! Strings
    CHARACTER(LEN=80)    :: Att

    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume file is not in define mode unless explicitly told otherwise
    IF ( PRESENT( DefMode ) ) THEN
       isDefMode = DefMode
    ELSE
       isDefMode = .FALSE.
    ENDIF

    ! Test if iLevId (dimension for level interfaces) is present
    IF ( PRESENT( iLevId ) ) THEN
       tmpIlevId = iLevId
    ELSE
       tmpIlevId = -1
    ENDIF

    !=======================================================================
    ! DEFINE VARIABLE
    !=======================================================================

    ! Reopen definition section, if necessary
    IF ( .not. isDefMode ) CALL NcBegin_Def( fId )

    VarCt = VarCt + 1

    ! number of dimensions
    nDim = 0
    IF ( lonId     >= 0 ) nDim = nDim + 1
    IF ( latId     >= 0 ) nDim = nDim + 1
    IF ( levId     >= 0 ) nDim = nDim + 1
    IF ( tmpIlevId >= 0 ) nDim = nDim + 1
    if ( timeId    >= 0 ) nDim = nDim + 1

    ! write dimensions
    ALLOCATE( VarDims(nDim) )
    Pos = 1
    IF ( lonId >= 0 ) THEN
       VarDims(Pos) = lonId
       Pos          = Pos + 1
    ENDIF
    IF ( latId >= 0 ) THEN
       VarDims(Pos) = latId
       Pos          = Pos + 1
    ENDIF
    IF ( levId >= 0 ) THEN
       VarDims(Pos) = levId
       Pos          = Pos + 1
    ENDIF
    IF ( tmpIlevId >= 0 ) THEN
       VarDims(Pos) = tmpIlevId
       Pos          = Pos + 1
    ENDIF
    IF ( timeId >= 0 ) THEN
       VarDims(Pos) = timeId
       Pos          = Pos + 1
    ENDIF

    ! Set data type
    IF ( DataType == 1 ) THEN
       NF_TYPE = NF_INT
    ELSEIF ( DataType == 4 ) THEN
       NF_TYPE = NF_FLOAT
    ELSEIF ( DataType == 8 ) THEN
       NF_TYPE = NF_DOUBLE
    ELSE
       NF_TYPE = NF_FLOAT
    ENDIF

    !-----------------------------------------------------------------------
    ! Define variable
    !-----------------------------------------------------------------------
    CALL NcDef_Variable( fId,  TRIM(VarName), NF_TYPE,              &
                         nDim, VarDims,       VarCt,     Compress  )
    DEALLOCATE( VarDims )

    !-----------------------------------------------------------------------
    ! Define variable atttibutes (some are optional)
    !-----------------------------------------------------------------------

    ! long_name (reuired)
    Att = 'long_name'
    CALL NcDef_Var_Attributes(  fId, VarCt, TRIM(Att), TRIM(VarLongName) )

    ! units (requited)
    Att = 'units'
    CALL NcDef_Var_Attributes(  fId, VarCt, TRIM(Att),  TRIM(VarUnit) )

    ! add_offset (optional)
    IF ( PRESENT( AddOffset ) ) THEN
       Att = 'add_offset'
       CALL NcDef_Var_Attributes( fId, VarCt, TRIM(Att), AddOffset )
    ENDIF

    ! scale_factor (optional)
    IF ( PRESENT( ScaleFactor ) ) THEN
       Att = 'scale_factor'
       CALL NcDef_Var_Attributes( fId, VarCt, TRIM(Att), ScaleFactor )
    ENDIF

    ! missing_value (optional but recommended)
    IF ( PRESENT( MissingValue ) ) THEN
       Att = '_FillValue'
       CALL NcDef_Var_Attributes( fId, VarCt, TRIM(Att),  MissingValue )
    ENDIF

    ! calendar (only used for time) -- skip if null string
    IF ( PRESENT( Calendar ) ) THEN
       IF ( LEN_TRIM( Calendar ) > 0 ) THEN
          Att = 'calendar'
          CALL NcDef_Var_Attributes( fId, VarCt, TRIM(Att), TRIM(Calendar) )
       ENDIF
    ENDIF

    ! axis (only used for index variables) -- skip if null string
    IF ( PRESENT( Axis ) ) THEN
       IF ( LEN_TRIM( Axis ) > 0 ) THEN
          Att = 'axis'
          CALL NcDef_Var_Attributes( fId, VarCt, TRIM(Att), TRIM(Axis) )
       ENDIF
    ENDIF

    ! averaging_method (optional) -- skip if null string
    IF ( PRESENT( AvgMethod ) ) THEN
       IF ( LEN_TRIM( AvgMethod ) > 0 ) THEN
          Att = 'averaging_method'
          CALL NcDef_Var_Attributes( fId, VarCt, TRIM(Att), TRIM(AvgMethod) )
       ENDIF
    ENDIF

    ! averaging_method (optional) -- skip if null string
    IF ( PRESENT( Positive ) ) THEN
       IF ( LEN_TRIM( Positive ) > 0 ) THEN
          Att = 'positive'
          CALL NcDef_Var_Attributes( fId, VarCt, TRIM(Att), TRIM(Positive) )
       ENDIF
    ENDIF

    ! Standard name (optional) -- skip if null string
    IF ( PRESENT( StandardName ) ) THEN
       IF ( LEN_TRIM( StandardName ) > 0 ) THEN
          Att = 'standard_name'
          CALL NcDef_Var_Attributes( fId, VarCt, TRIM(Att), TRIM(StandardName))
       ENDIF
    ENDIF

    ! Formula terms (optional) -- skip if null string
    IF ( PRESENT( FormulaTerms ) ) THEN
       IF ( LEN_TRIM( FormulaTerms ) > 0 ) THEN
          Att = 'formula_terms'
          CALL NcDef_Var_Attributes( fId, VarCt, TRIM(Att), TRIM(FormulaTerms))
       ENDIF
    ENDIF

    ! Number of updates
    IF ( PRESENT( nUpdates ) ) THEN
       IF ( nUpdates > 0.0 ) THEN
          Att = 'number_of_updates'
          CALL NcDef_Var_Attributes( fId, VarCt, TRIM(Att), nUpdates )
       ENDIF
    ENDIF

    ! Close definition section, if necessary
    IF ( .not. isDefMode ) CALL NcEnd_Def( fId )

  END SUBROUTINE NC_Var_Def
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Chunk
!
! !DESCRIPTION: Turns on chunking for a netCDF variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Nc_Var_Chunk( fId, vId, ChunkSizes, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)  :: fId            ! NetCDF file ID
    INTEGER, INTENT(IN)  :: vId            ! NetCDF variable ID
    INTEGER, INTENT(IN)  :: ChunkSizes(:)  ! NetCDF chunk sizes for each dim
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC             ! Success or failure?
!
! !REMARKS:
!  RC will return an error (nonzero) status if chunking cannot be activated.
!  Most often, this is because support for netCDF-4 compression is disabled,
!  or if the netCDF file is not a netCDF-4 file.  In this case, RC will have
!  an error code of -111.
!
! !REVISION HISTORY:
!  28 Aug 2017 - R. Yantosca - Initial version
!  11 Sep 2017 - R. Yantosca - Do not call NF_DEF_VAR_CHUNKING if the netCDF
!                               library was built w/o compression enabled
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( NC_HAS_COMPRESSION )

    ! Turn on chunking for this variable
    ! But only if the netCDF library supports it
    RC = NF_Def_Var_Chunking( fId, vId, NF_CHUNKED, ChunkSizes )

#else

    ! Otherwise return success
    RC = 0

#endif

  END SUBROUTINE Nc_Var_Chunk
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_R8_0d
!
! !DESCRIPTION: Writes data of a 0-D double precision variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_R8_0D( fId, VarName, Var )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId           ! file ID
    CHARACTER(LEN=*), INTENT(IN)  :: VarName       ! variable name
    REAL(kind=8)                  :: Var           ! Variable to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  25 Aug 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    ! Write to netCDF file
    CALL NcWr( Var, fId, VarName )

  END SUBROUTINE NC_VAR_WRITE_R8_0d
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_R8_1d
!
! !DESCRIPTION: Writes data of a 1-D double precision variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_R8_1D( fId, VarName, Arr1D )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId           ! file ID
    CHARACTER(LEN=*), INTENT(IN)  :: VarName       ! variable name
    REAL(kind=8),     POINTER     :: Arr1D(:)      ! array to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  16 Jun 2014 - R. Yantosca - Now use simple arrays instead of allocating
!  19 Sep 2016 - R. Yantosca - Renamed to NC_VAR_WRITE_R8_1D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays
    INTEGER :: St1d(1), Ct1d(1)

    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    ! Set start & count arrays
    St1d(1) = 1
    Ct1d(1) = SIZE( Arr1d, 1 )

    ! Write to netCDF file
    CALL NcWr( Arr1d, fId, VarName, St1d, Ct1d )

  END SUBROUTINE NC_VAR_WRITE_R8_1D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_R8_2d
!
! !DESCRIPTION: Writes data of a 2-D double precision variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_R8_2D( fId, VarName, Arr2D )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: fId            ! file ID
    CHARACTER(LEN=*), INTENT(IN) :: VarName        ! variable name
    REAL(kind=8),     POINTER    :: Arr2D(:,:)     ! array to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  16 Jun 2014 - R. Yantosca - Now use simple arrays instead of allocating
!  19 Sep 2016 - R. Yantosca - Renamed to NC_VAR_WRITE_R8_2D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays
    INTEGER :: St2d(2), Ct2d(2)

    ! Scalars
    INTEGER :: I,       nDim

    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    ! Set start & count arrays
    nDim = 2
    DO I =1, nDim
       St2d(I) = 1
       Ct2d(I) = SIZE( Arr2d, I )
    ENDDO

    ! Write to netCDF file
    CALL NcWr( Arr2d, fId, VarName, St2d, Ct2d )

  END SUBROUTINE NC_VAR_WRITE_R8_2D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_R8_3D
!
! !DESCRIPTION: Writes data of a 3-D double precision variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_R8_3D( fId, VarName, Arr3D )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: fId            ! file ID
    CHARACTER(LEN=*), INTENT(IN) :: VarName        ! variable name
    REAL(kind=8),     POINTER    :: Arr3D(:,:,:)   ! array to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  16 Jun 2014 - R. Yantosca - Now use simple arrays instead of allocating
!  19 Sep 2016 - R. Yantosca - Renamed to NC_VAR_WRITE_R8_3D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays
    INTEGER :: St3d(3), Ct3d(3)

    ! Scalars
    INTEGER :: I,       nDim

    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    ! Set start & count arrays
    nDim = 3
    DO I = 1, nDim
       St3d(I) = 1
       Ct3d(I) = SIZE( Arr3d, I )
    ENDDO

    ! Write data to netCDF file
    CALL NcWr( Arr3d, fId, VarName, St3d, Ct3d )

  END SUBROUTINE NC_VAR_WRITE_R8_3D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_r8_4d
!
! !DESCRIPTION: Writes data of a 4-D double precision variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_R8_4D( fId, VarName, Arr4D )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: fId            ! file ID
    CHARACTER(LEN=*), INTENT(IN) :: VarName        ! variable name
    REAL(kind=8),     POINTER    :: Arr4D(:,:,:,:) ! array to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  16 Jun 2014 - R. Yantosca - Now use simple arrays instead of allocating
!  19 Sep 2016 - R. Yantosca - Renamed to NC_VAR_WRITE_R8_4D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays
    INTEGER :: St4d(4), Ct4d(4)

    ! Scalars
    INTEGER :: I,       nDim

    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    ! Set start & count arrays
    nDim = 4
    DO I = 1, nDim
       St4d(I) = 1
       Ct4d(I) = SIZE( Arr4d, I )
    ENDDO

    ! Write to netCDF file
    CALL NcWr( Arr4d, fId, VarName, St4d, Ct4d )

  END SUBROUTINE NC_VAR_WRITE_R8_4D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_R4_0d
!
! !DESCRIPTION: Writes data of a 0-D single-precision variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_R4_0d( fId, VarName, Var )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId           ! file ID
    CHARACTER(LEN=*), INTENT(IN)  :: VarName       ! variable name
    REAL(kind=4)                  :: Var           ! Variable to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  25 Aug 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    ! Write to netCDF file
    CALL NcWr( Var, fId, VarName )

  END SUBROUTINE NC_VAR_WRITE_R4_0D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_r4_1d
!
! !DESCRIPTION: Writes data of a single precision variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_R4_1D( fId, VarName, Arr1D )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: fId            ! file ID
    CHARACTER(LEN=*), INTENT(IN) :: VarName        ! variable name
    REAL(kind=4),     POINTER    :: Arr1D(:)       ! array to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  16 Jun 2014 - R. Yantosca - Now use simple arrays instead of allocating
!  19 Sep 2016 - R. Yantosca - Renamed to NC_VAR_WRITE_R4_1D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays
    INTEGER :: St1d(1), Ct1d(1)

    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    ! Set start & count arrays
    St1d(1) = 1
    Ct1d(1) = SIZE( Arr1d, 1 )

    ! Write to netCDF file
    CALL NcWr( Arr1d, fId, VarName, St1d, Ct1d )

  END SUBROUTINE NC_VAR_WRITE_R4_1D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_r4_2D
!
! !DESCRIPTION: Writes data of a 2-D single precision variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_R4_2D( fId, VarName, Arr2D )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: fId            ! file ID
    CHARACTER(LEN=*), INTENT(IN) :: VarName        ! variable name
    REAL(kind=4),     POINTER    :: Arr2D(:,:)     ! array to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  16 Jun 2014 - R. Yantosca - Now use simple arrays instead of allocating
!  19 Sep 2016 - R. Yantosca - Renamed to NC_VAR_WRITE_R4_2D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays
    INTEGER :: St2d(2), Ct2d(2)

    ! Scalars
    INTEGER :: I,       nDim

    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    ! Set start & count arrays
    nDim = 2
    DO I = 1, nDim
       St2d(I) = 1
       Ct2d(I) = SIZE( Arr2d, I )
    ENDDO

    ! Write to netCDF file
    CALL NcWr( Arr2d, fId, VarName, St2d, Ct2d )

  END SUBROUTINE NC_VAR_WRITE_R4_2D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_r4_3d
!
! !DESCRIPTION: Writes data of a 3-D single precision variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_R4_3D( fId, VarName, Arr3D )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId            ! file ID
    CHARACTER(LEN=*), INTENT(IN)  :: VarName        ! variable name
    REAL(kind=4),     POINTER     :: Arr3D(:,:,:)   ! array to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  16 Jun 2014 - R. Yantosca - Now use simple arrays instead of allocating
!  19 Sep 2016 - R. Yantosca - Renamed to NC_VAR_WRITE_R4_3D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays
    INTEGER :: St3d(3), Ct3d(3)

    ! Scalars
    INTEGER :: I, nDim

    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    ! Set start & count arrays
    nDim = 3
    DO I = 1, nDim
       St3d(I) = 1
       Ct3d(I) = SIZE( Arr3d, I )
    ENDDO

    ! Write to netCDF file
    CALL NcWr( Arr3d, fId, VarName, St3d, Ct3d )

  END SUBROUTINE NC_VAR_WRITE_R4_3D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_r4_4d
!
! !DESCRIPTION: Writes data of a 4-D single precision variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_R4_4D( fId, VarName, Arr4D )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: fId            ! file ID
    CHARACTER(LEN=*), INTENT(IN) :: VarName        ! variable name
    REAL(kind=4),     POINTER    :: Arr4D(:,:,:,:) ! array to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  16 Jun 2014 - R. Yantosca - Now use simple arrays instead of allocating
!  19 Sep 2016 - R. Yantosca - Renamed to NC_VAR_WRITE_R4_1D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays
    INTEGER :: St4d(4), Ct4d(4)

    ! Scalars
    INTEGER :: I,       nDim

    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    nDim = 4
    DO I = 1, nDim
       St4d(I) = 1
       Ct4d(I) = SIZE( Arr4d, I )
    ENDDO

    ! Write to netCDF file
    CALL NcWr( Arr4d, fId, VarName, St4d, Ct4d )

  END SUBROUTINE NC_VAR_WRITE_R4_4D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_Int_0d
!
! !DESCRIPTION: Writes data of a 0-D integer variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_INT_0d( fId, VarName, Var )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)  :: fId           ! file ID
    CHARACTER(LEN=*), INTENT(IN)  :: VarName       ! variable name
    INTEGER                       :: Var           ! Variable to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  25 Aug 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    ! Write to netCDF file
    CALL NcWr( Var, fId, VarName )

  END SUBROUTINE NC_VAR_WRITE_INT_0D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_int_1d
!
! !DESCRIPTION: Writes data of an 1-D integer variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_INT_1D( fId, VarName, Arr1D )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: fId            ! file ID
    CHARACTER(LEN=*), INTENT(IN) :: VarName        ! variable name
    INTEGER,          POINTER    :: Arr1D(:)       ! array to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  16 Jun 2014 - R. Yantosca - Now use simple arrays instead of allocating
!  19 Sep 2016 - R. Yantosca - Renamed to NC_VAR_WRITE_INT_1D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays
    INTEGER :: St1d(1), Ct1d(1)

    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    ! Set start & count arrays
    St1d(1) = 1
    Ct1d(1) = SIZE( Arr1d, 1 )

    ! Write to netCDF file
    CALL NcWr( Arr1d, fId, VarName, St1d, Ct1d )

  END SUBROUTINE NC_VAR_WRITE_INT_1D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_int_2d
!
! !DESCRIPTION: writes data of an 2-D integer variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_INT_2D( fId, VarName, Arr2D )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: fId            ! file ID
    CHARACTER(LEN=*), INTENT(IN) :: VarName        ! variable name
    INTEGER,          POINTER    :: Arr2D(:,:)     ! array to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  16 Jun 2014 - R. Yantosca - Now use simple arrays instead of allocating
!  19 Sep 2016 - R. Yantosca - Renamed to NC_VAR_WRITE_INT_2D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays
    INTEGER :: St2d(2), Ct2d(2)

    ! Scalars
    INTEGER :: I,       nDim

    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    ! Set start & count arrays
    nDim = 2
    DO I = 1, nDim
       St2d(I) = 1
       Ct2d(I) = SIZE( Arr2d, I )
    ENDDO

    ! Write to netCDF file
    CALL NcWr( Arr2d, fId, VarName, St2d, Ct2d )

  END SUBROUTINE NC_VAR_WRITE_INT_2D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_int_3d
!
! !DESCRIPTION: writes data of an 3-D integer variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_INT_3D( fId, VarName, Arr3D )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: fId            ! file ID
    CHARACTER(LEN=*), INTENT(IN) :: VarName        ! variable name
    INTEGER,          POINTER    :: Arr3D(:,:,:)   ! array to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  16 Jun 2014 - R. Yantosca - Now use simple arrays instead of allocating
!  19 Sep 2016 - R. Yantosca - Renamed to NC_VAR_WRITE_INT_3D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays
    INTEGER :: St3d(3), Ct3d(3)

    ! Scalars
    INTEGER :: I,       nDim

    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    ! Set start & count arrays
    nDim = 3
    DO I = 1, nDim
       St3d(I) = 1
       Ct3d(I) = SIZE( Arr3d, I )
    ENDDO

    ! Write to netCDF file
    CALL NcWr( Arr3d, fId, trim(VarName), St3d, Ct3d )

  END SUBROUTINE NC_VAR_WRITE_INT_3D
!EOC
!------------------------------------------------------------------------------
!       NcdfUtilities: by Harvard Atmospheric Chemistry Modeling Group        !
!                      and NASA/GSFC, SIVO, Code 610.3                        !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_Var_Write_int_4d
!
! !DESCRIPTION: writes data of an 4-Dinteger variable.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE NC_VAR_WRITE_INT_4D( fId, VarName, Arr4D )
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: fId            ! file ID
    CHARACTER(LEN=*), INTENT(IN) :: VarName        ! variable name
    INTEGER,          POINTER    :: Arr4D(:,:,:,:) ! array to be written
!
! !REMARKS:
!  Assumes that you have:
!  (1) A netCDF library (either v3 or v4) installed on your system
!  (2) The NcdfUtilities package (from Bob Yantosca) source code
!                                                                             .
!  Although this routine was generated automatically, some further
!  hand-editing may be required.
!
! !REVISION HISTORY:
!  15 Jun 2012 - C. Keller   - Initial version
!  16 Jun 2014 - R. Yantosca - Now use simple arrays instead of allocating
!  19 Sep 2016 - R. Yantosca - Renamed to NC_VAR_WRITE_INT_1D
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Arrays
    INTEGER :: St4d(4), Ct4d(4)

    ! Scalars
    INTEGER :: I, nDim

    !--------------------------------
    ! WRITE DATA
    !--------------------------------

    ! Set start & count arrays
    nDim = 4
    DO I = 1, nDim
       St4d(I) = 1
       Ct4d(I) = SIZE( Arr4d, I )
    ENDDO

    ! Write to netCDF file
    CALL NcWr( Arr4d, fId, VarName, St4d, Ct4d )

  END SUBROUTINE NC_VAR_WRITE_INT_4D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Tau0
!
! !DESCRIPTION: Function GET\_TAU0\_6A returns the corresponding TAU0 value
!  for the first day of a given MONTH of a given YEAR.  This is necessary to
!  index monthly mean binary punch files, which are used as input to GEOS-Chem.
!\\
!\\
!  This function takes 3 mandatory arguments (MONTH, DAY, YEAR) and 3
!  optional arguments (HOUR, MIN, SEC).  It is intended to replace the current
!  2-argument version of GET\_TAU0.  The advantage being that GET\_TAU0\_6A
!  can compute a TAU0 for any date and time in the GEOS-Chem epoch, rather
!  than just the first day of each month.  Overload this w/ an interface so
!  that the user can also choose the version of GET\_TAU0 w/ 2 arguments
!  (MONTH, YEAR), which is the prior version.
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_TAU0( MONTH, DAY, YEAR, HOUR, MIN, SEC ) RESULT( THIS_TAU0 )
!
! !USES:
!
    USE JULDAY_MOD, ONLY : JULDAY
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)           :: MONTH
    INTEGER, INTENT(IN)           :: DAY
    INTEGER, INTENT(IN)           :: YEAR
    INTEGER, INTENT(IN), OPTIONAL :: HOUR
    INTEGER, INTENT(IN), OPTIONAL :: MIN
    INTEGER, INTENT(IN), OPTIONAL :: SEC
!
! !RETURN VALUE:
!
    REAL*8                        :: THIS_TAU0   ! TAU0 timestamp
!
! !REMARKS:
!  TAU0 is hours elapsed since 00:00 GMT on 01 Jan 1985.
!
! !REVISION HISTORY:
!  (1 ) 1985 is the first year of the GEOS epoch.
!  (2 ) Add TAU0 values for years 1985-2001 (bmy, 8/1/00)
!  (3 ) Correct error for 1991 TAU values.  Also added 2002 and 2003.
!        (bnd, bmy, 1/4/01)
!  (4 ) Updated comments  (bmy, 9/26/01)
!  (5 ) Now references JULDAY from "julday_mod.f" (bmy, 11/20/01)
!  (6 ) Now references ERROR_STOP from "error_mod.f"  (bmy, 10/15/02)
!  20 Nov 2009 - R. Yantosca - Added ProTeX header
!  10 Jul 2014 - R. Yantosca - Add this routine as a PRIVATE module variable
!                              to prevent ncdf_mod.F90 from using bpch2_mod.F
!  10 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: TMP_HOUR, TMP_MIN, TMP_SEC
    REAL*8  :: DAYS

    !=======================================================================
    ! GET_TAU0 begins here!
    !=======================================================================

    ! Error checking
    IF ( MONTH < 1 .or. MONTH > 12 ) THEN
       WRITE( 6, 100 )
100    FORMAT( 'Invalid MONTH selection!  STOP in GET_TAU0 (ncdf_mod.F90)!' )
       STOP
    ENDIF

    ! Error checking
    IF ( DAY < 1 .or. DAY > 31 ) THEN
       WRITE( 6, 110 )
110    FORMAT( 'Invalid DAY selection!  STOP in GET_TAU0 (ncdf_mod.F90)!' )
       STOP
    ENDIF

    ! If HOUR isn't passed, default to 0
    IF ( PRESENT( HOUR ) ) THEN
       TMP_HOUR = HOUR
    ELSE
       TMP_HOUR = 0
    ENDIF

    ! If MIN isn't passed, default to 0
    IF ( PRESENT( MIN ) ) THEN
       TMP_MIN = MIN
    ELSE
       TMP_MIN = 0
    ENDIF

    ! If SEC isn't passed, default to 0
    IF ( PRESENT( SEC ) ) THEN
       TMP_SEC = SEC
    ELSE
       TMP_SEC = 0
    ENDIF

    ! Number of days since midnight on 1/1/1985
    THIS_TAU0 = JULDAY( YEAR, MONTH, DBLE( DAY ) ) - 2446066.5d0

    ! Multiply by 24 to get hours since 1/1/1985
    ! Also add in the hours elapsed since midnight on this date
    THIS_TAU0 = ( THIS_TAU0 * 24d0 ) + ( TMP_HOUR         ) + &
                ( TMP_MIN   / 60d0 ) + ( TMP_SEC / 3600d0 )

  END FUNCTION GET_TAU0
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nc_IsModelLevels
!
! !DESCRIPTION: Function NC\_ISMODELLEVELS returns true if (and only if) the
!  long name of the level variable name of the given file ID contains the
!  character "GEOS-Chem level".
!\\
!\\
! !INTERFACE:
!
  FUNCTION NC_ISMODELLEVEL( fID, lev_name ) RESULT ( IsModelLevel )
!
! !USES:
!
#   include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN) :: fID        ! file ID
    CHARACTER(LEN=*), INTENT(IN) :: lev_name   ! level variable name
!
! !RETURN VALUE:
!
    LOGICAL                      :: IsModelLevel
!
! !REVISION HISTORY:
!  12 Dec 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL                :: HasLngN
    CHARACTER(LEN=255)     :: a_name, LngName
    INTEGER                :: a_type

    !=======================================================================
    ! NC_ISMODELLEVEL begins here!
    !=======================================================================

    ! Init
    IsModelLevel = .FALSE.

    ! Check if there is a long_name attribute
    a_name = "long_name"
    HasLngN = Ncdoes_Attr_Exist ( fId, TRIM(lev_name), TRIM(a_name), a_type )

    ! Only if attribute exists...
    IF ( HasLngN ) THEN
       ! Read attribute
       CALL NcGet_Var_Attributes( fID, TRIM(lev_name), TRIM(a_name), LngName )

       ! See if this is a GEOS-Chem model level
       IF ( INDEX(TRIM(LngName),"GEOS-Chem level") > 0 ) THEN
          IsModelLevel = .TRUE.
       ENDIF
    ENDIF

  END FUNCTION NC_ISMODELLEVEL
!EOC
END MODULE NCDF_MOD
