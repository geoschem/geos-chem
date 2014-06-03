!-----------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      USE m_netcdf_io_checks

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: NC_OPEN
      PUBLIC  :: NC_CREATE
      PUBLIC  :: NC_VAR_DEF
      PUBLIC  :: NC_VAR_WRITE
      PUBLIC  :: NC_CLOSE
      PUBLIC  :: NC_READ_TIME
      PUBLIC  :: NC_READ_TIME_YYYYMMDDhh
      PUBLIC  :: NC_READ_VAR
      PUBLIC  :: NC_READ_ARR
      PUBLIC  :: NC_GET_REFDATETIME

      PUBLIC  :: NC_READ
      PUBLIC  :: NC_READ_GRID
      PUBLIC  :: NC_WRITE
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: GET_TIDX
      PRIVATE :: TIMEUNIT_CHECK
!
! !REVISION HISTORY:
!  27 Jul 2012 - C. Keller - Initial version
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES:
!
      INTERFACE NC_WRITE
         MODULE PROCEDURE NC_WRITE_3D
         MODULE PROCEDURE NC_WRITE_4D
      END INTERFACE

      INTERFACE NC_VAR_WRITE
         MODULE PROCEDURE NC_VAR_WRITE_INT
         MODULE PROCEDURE NC_VAR_WRITE_R4
         MODULE PROCEDURE NC_VAR_WRITE_R8
      END INTERFACE

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NC_OPEN
!
! !DESCRIPTION: Simple wrapper routine to open the given netCDF file. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_OPEN ( FileName, fID ) 
!
! !ARGUMENTS:
!   
      CHARACTER(LEN=*), INTENT(IN   )  :: FileName 
      INTEGER,          INTENT(  OUT)  :: fID 
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NC_CLOSE
!
! !DESCRIPTION: Simple wrapper routine to close the given lun. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_CLOSE ( fID ) 
!
! !ARGUMENTS:
!   
      INTEGER,          INTENT(IN   )  :: fID
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NC_READ_TIME
!
! !DESCRIPTION: Subroutine NC\_READ\_TIME reads the time variable of the
! given fID and returns the time slices and unit. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_READ_TIME ( fID,     nTime,        timeUnit, &
                                timeVec, timeCalendar, RC       ) 
!
! !ARGUMENTS:
!   
      INTEGER,          INTENT(IN   )            :: fID
      INTEGER,          INTENT(  OUT)            :: nTime 
      CHARACTER(LEN=*), INTENT(  OUT)            :: timeUnit 
      INTEGER,          POINTER,       OPTIONAL  :: timeVec(:)
      CHARACTER(LEN=*), INTENT(  OUT), OPTIONAL  :: timeCalendar
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
      LOGICAL                :: hasTime
      CHARACTER(LEN=255)     :: v_name             ! netCDF variable name 
      CHARACTER(LEN=255)     :: a_name             ! netCDF attribute name
      CHARACTER(LEN=255)     :: a_val              ! netCDF attribute value
      INTEGER                :: st1d(1), ct1d(1)   ! For 1D arrays    
      INTEGER, ALLOCATABLE   :: tmpTime(:)

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
         CALL NcGet_Var_Attributes( fId, v_name, 'calendar', & 
                                    timeCalendar              )
      ENDIF

      END SUBROUTINE NC_READ_TIME
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NC_READ_VAR
!
! !DESCRIPTION: Subroutine NC\_READ\_VAR reads the given variable from the
! given fID and returns the corresponding variable values and units. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_READ_VAR ( fID,    Var, nVar, varUnit,  &
                               varVec, RC                    ) 
!
! !ARGUMENTS:
!   
      INTEGER,          INTENT(IN   )            :: fID
      CHARACTER(LEN=*), INTENT(IN   )            :: var 
      INTEGER,          INTENT(  OUT)            :: nVar
      CHARACTER(LEN=*), INTENT(  OUT)            :: varUnit 
      REAL,             POINTER,       OPTIONAL  :: varVec(:)
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
      LOGICAL                :: hasVar
      CHARACTER(LEN=255)     :: v_name             ! netCDF variable name 
      CHARACTER(LEN=255)     :: a_name             ! netCDF attribute name
      CHARACTER(LEN=255)     :: a_val              ! netCDF attribute value
      INTEGER                :: st1d(1), ct1d(1)   ! For 1D arrays    
      INTEGER, ALLOCATABLE   :: tmpVar(:)

      !=================================================================
      ! NC_READ_VAR begins here
      !=================================================================

      ! Init
      RC      = 0
      nVar    = 0
      hasVar  = .FALSE.

      ! Variable name
      v_name = var

      ! Check if variable exists 
      hasVar = Ncdoes_Dim_Exist ( fID, TRIM(v_name) ) 

      ! Return here if no lon variable defined 
      IF ( .NOT. hasVar ) RETURN 
      
      ! Get dimension length
      CALL Ncget_Dimlen ( fID, TRIM(v_name), nVar )

      ! Read time vector from file. 
      ! Pass to output
      IF ( PRESENT(VarVec) ) THEN
         IF ( ASSOCIATED( VarVec ) ) DEALLOCATE(VarVec)
         ALLOCATE ( varvec(nVar) )
         ALLOCATE ( tmpVar(nVar) )
         st1d = (/ 1    /)
         ct1d = (/ nVar /)
         CALL NcRd( tmpVar, fID, TRIM(v_name), st1d, ct1d )
         varvec(:) = tmpVar
         DEALLOCATE(tmpVar)
      ENDIF

      ! Read units attribute
      a_name = "units"
      CALL NcGet_Var_Attributes( fID,          TRIM(v_name), &
                                 TRIM(a_name), varUnit     )

      END SUBROUTINE NC_READ_VAR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nc_read_arr
!
! !DESCRIPTION: Routine to read data from a netCDF file. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_READ_ARR( fID,   ncVar,   lon1, lon2,  lat1,  &
                              lat2,  lev1,    lev2, time1, time2, & 
                              ncArr, VarUnit, RC                   ) 
!
! !USES:
!
      USE CHARPAK_MOD,        ONLY : TRANLC

      IMPLICIT NONE
#     include "netcdf.inc"
!
! !ARGUMENTS:
!   
      ! Required input
      INTEGER,          INTENT(IN)            :: fID 
      CHARACTER(LEN=*), INTENT(IN)            :: ncVar        ! variable to read
      INTEGER,          INTENT(IN)            :: lon1,  lon2 
      INTEGER,          INTENT(IN)            :: lat1,  lat2
      INTEGER,          INTENT(IN)            :: lev1,  lev2
      INTEGER,          INTENT(IN)            :: time1, time2

      ! Required output: Array to write data
      REAL*4,           POINTER               :: ncArr(:,:,:,:)

      ! Optional output
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: VarUnit 

      ! Error handling
      INTEGER,          INTENT(INOUT)         :: RC
!
! !REVISION HISTORY:
!  27 Jul 2012 - C. Keller - Initial version
!  18 Jan 2012 - C. Keller - Now reads 4D, 3D, and 2D arrays, with
!                            optional dimensions level and time.
!  18 Apr 2012 - C. Keller - Now also read & apply offset and scale factor
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
      REAL*8                 :: corr      ! netCDF attribute value 

      ! Arrays for netCDF start and count values
      INTEGER                :: nlon,  nlat, nlev, ntime
      INTEGER                :: nclev, nctime 
      INTEGER                :: st2d(2), ct2d(2)   ! For 2D arrays 
      INTEGER                :: st3d(3), ct3d(3)   ! For 3D arrays 
      INTEGER                :: st4d(4), ct4d(4)   ! For 4D arrays 

      ! Temporary arrays
      REAL*4, ALLOCATABLE    :: TMPARR_3D(:,:,:)
      REAL*4, ALLOCATABLE    :: TMPARR_2D(:,:)

      ! Logicals
      LOGICAL                :: ReadAtt

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

      ! Extract dimensions to read 
      nLon = lon2 - lon1 + 1
      nLat = lat2 - lat1 + 1
      IF ( lev1 > 0 ) THEN
         nLev = lev2 - lev1 + 1
      ELSE
         nLev = 0
      ENDIF
      IF ( time1 > 0 ) THEN
         ntime = time2 - time1 + 1
      ELSE
         ntime = 0
      ENDIF

      ! Set dimensions of output array
      ! --> must have at least dimension 1
      nclev  = max(nlev ,1)
      nctime = max(ntime,1)

      !----------------------------------------
      ! Read array
      !----------------------------------------

      ! Variable name
      v_name = TRIM(ncVar)
   
      ! Allocate the array
      ALLOCATE ( ncArr( nLon, nLat, ncLev, ncTime ) )
      ncArr = 0d0  

      ! To read a 4D array:
      IF ( (nLev>0) .AND. (nTime>0) ) THEN
         st4d = (/ lon1, lat1, lev1, time1 /)
         ct4d = (/ nlon, nlat, nlev, ntime /)
         CALL NcRd( ncArr, fId, TRIM(v_name), st4d, ct4d )
   
      ! To read a 3D array:
      ! ==> Level defined but not time:
      ELSEIF ( nLev>0 ) THEN
         ALLOCATE ( TMPARR_3D( nLon, nLat, nLev ) )
         st3d = (/ lon1, lat1, lev1 /)
         ct3d = (/ nlon, nlat, nlev /)
         CALL NcRd( TMPARR_3D, fId, TRIM(v_name), st3d, ct3d )
         ncArr(:,:,:,1) = TMPARR_3D(:,:,:)

      ! ==> Time defined but not level:
      ELSEIF ( nTime>0 ) THEN

        ALLOCATE ( TMPARR_3D( nLon, nLat, nTime ) ) 
         st3d   = (/ lon1, lat1, time1 /)
         ct3d   = (/ nlon, nlat, ntime /)
         CALL NcRd( TMPARR_3D, fId, TRIM(v_name), st3d, ct3d )
         ncArr(:,:,1,:) = TMPARR_3D(:,:,:)

      ! Otherwise, read a 2D array (lon and lat only):
      ELSE
         ALLOCATE ( TMPARR_2D( nLon, nLat ) )
         st2d   = (/ lon1, lat1 /)
         ct2d   = (/ nlon, nlat /)
         CALL NcRd( TMPARR_2D, fId, TRIM(v_name), st2d, ct2d )
         ncArr(:,:,1,1) = TMPARR_2D(:,:)
      ENDIF

      ! ------------------------------------------
      ! Eventually apply scale / offset factors
      ! ------------------------------------------

      ! Check for scale factor
      a_name  = "scale_factor"
      ReadAtt = Ncdoes_Attr_Exist ( fId, TRIM(v_name), TRIM(a_name) ) 

      IF ( ReadAtt ) THEN
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),corr)
         ncArr(:,:,:,:) = ncArr(:,:,:,:) * corr
      ENDIF

      ! Check for offset factor
      a_name  = "add_offset"
      ReadAtt = Ncdoes_Attr_Exist ( fId, TRIM(v_name), TRIM(a_name) ) 

      IF ( ReadAtt ) THEN
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),corr)
         ncArr(:,:,:,:) = ncArr(:,:,:,:) + corr
      ENDIF

      ! ----------------------------
      ! Read optional arguments
      ! ----------------------------

      ! Read units
      IF ( PRESENT(VarUnit) )THEN
         a_name = "units"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)
         VarUnit = TRIM(a_val)
      ENDIF

      !=================================================================
      ! Cleanup and quit
      !=================================================================

      ! Deallocate internal variables
      IF ( ALLOCATED ( TMPARR_3D ) ) DEALLOCATE ( TMPARR_3D )
      IF ( ALLOCATED ( TMPARR_2D ) ) DEALLOCATE ( TMPARR_2D )

      ! Return w/ success
      RC = 0 

      END SUBROUTINE NC_READ_ARR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NC_READ_TIME_YYYYMMDDhh
!
! !DESCRIPTION: Returns a vector containing the datetimes (YYYYMMDDhh) of 
! all time slices in the netCDF file.
!\\
! !INTERFACE:
!
      SUBROUTINE NC_READ_TIME_YYYYMMDDhh( fID, nTime, all_YYYYMMDDhh, &
                                          timeUnit,   refYear, RC )
!
! !USES:
!
      USE JULDAY_MOD, ONLY : JULDAY, CALDATE
!
! !ARGUMENTS:
!   
      INTEGER,          INTENT(IN   )           :: fID
      INTEGER,          INTENT(INOUT)           :: nTime
      INTEGER,          POINTER                 :: all_YYYYMMDDhh(:)
      CHARACTER(LEN=*), INTENT(  OUT), OPTIONAL :: timeUnit
      INTEGER,          INTENT(  OUT), OPTIONAL :: refYear 
      INTEGER,          INTENT(INOUT)           :: RC 
!
!
! !REVISION HISTORY:
!  27 Jul 2012 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=255)  :: ncUnit
      INTEGER, POINTER    :: tVec(:) => NULL()
      INTEGER             :: refYr, refMt, refDy, refHr
      INTEGER             :: T, YYYYMMDD, hhmmss 
      REAL*8              :: realrefDy, refJulday, tJulday

      !=================================================================
      ! NC_READ_TIME_YYYYMMDDhh begins here 
      !=================================================================

      ! Read time vector
      CALL NC_READ_TIME ( fID, nTime, ncUnit, timeVec=tVec, RC=RC ) 
      IF ( RC/=0 ) RETURN 

      ! If nTime is zero, return here!
      IF ( nTime == 0 ) RETURN

      ! Get reference date in julian days
      CALL NC_GET_REFDATETIME ( ncUnit, refYr, refMt, &
                                refDy,  refHr, RC      ) 
      IF ( RC /= 0 ) RETURN
      realrefDy = refDy + ( MAX(0,refHr) / 24d0 )
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
!         PRINT *, 'Reference julian day increased by one day!!!'
      ENDIF

      ! Get calendar dates
      IF ( ASSOCIATED ( all_YYYYMMDDhh ) ) DEALLOCATE( all_YYYYMMDDhh )
      ALLOCATE( all_YYYYMMDDhh(nTime) )
      all_YYYYMMDDhh = 0 
 
      ! Construct julian date for every available time slice
      DO T = 1, nTime
         tJulday = real(tVec(T), kind=8)
         IF ( refHr >= 0 ) tJulday = tJulday / 24.d0
         tJulday = tJulday + refJulday
         CALL CALDATE ( tJulday, YYYYMMDD, hhmmss )
         all_YYYYMMDDhh(T) = YYYYMMDD * 100 & 
                           + FLOOR( MOD(hhmmss,1000000) / 1.0d4 )
      ENDDO

      ! Cleanup
      IF ( ASSOCIATED( tVec ) ) DEALLOCATE( tVec )

      ! Return
      IF ( PRESENT(timeUnit) ) timeUnit = ncUnit
      IF ( PRESENT(refYear ) ) refYear  = refYr 
      RC = 0

      END SUBROUTINE NC_READ_TIME_YYYYMMDDhh 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NC_GET_REFDATETIME 
!
! !DESCRIPTION: Returns the reference datetime (tYr / tMt / tDy / tHr )
! of the provided time unit. For now, supported formats are 
! "days since YYYY-MM-DD" and "hours since YYYY-MM-DD HH:MM:SS". For
! times in days since refdate, the returned reference hour rHr is set 
! to -1.
!\\
! !INTERFACE:
!
      SUBROUTINE NC_GET_REFDATETIME( tUnit, tYr, tMt, tDy, tHr, RC ) 
!
! !USES:
!
      USE CHARPAK_MOD,        ONLY : TRANLC
!
! !ARGUMENTS:
!   
      ! Required
      CHARACTER(LEN=*), INTENT( IN)    :: tUnit 
      INTEGER,          INTENT(OUT)    :: tYr
      INTEGER,          INTENT(OUT)    :: tMt
      INTEGER,          INTENT(OUT)    :: tDy
      INTEGER,          INTENT(OUT)    :: tHr
      INTEGER,          INTENT(INOUT)  :: RC 
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
      INTEGER               :: TTYPE, STAT, L1, L2
      INTEGER               :: STRLEN, I

      !=================================================================
      ! NC_GET_REFDATETIME starts here 
      !=================================================================

      ! Init
      LOC = 'NC_GET_REFDATETIME (ncdf_mod.F)'

      ! Check length of time unit string. This must be at least 21 
      ! (==> "days since YYYY:MM:DD" is of length 21)
      STRLEN = LEN(tUnit)
      IF ( STRLEN < 21 ) THEN
         PRINT *, 'Time unit string too short!' 
         RC = -999; RETURN
      ENDIF

      ! ----------------------------------------------------------------------
      ! Determine time unit type
      ! ----------------------------------------------------------------------

      ! Mirror time unit and convert to lower case
      MIRRUNIT = tUnit
      CALL TRANLC( MIRRUNIT )

      ! Check for 'hours since'. If true, set TTYPE to 1 and set the
      ! begin of the effective date string to 12. Also check if the time
      ! string is at least of length 25, which is required for this
      ! unit.
      IF ( MIRRUNIT(1:11) == 'hours since' ) THEN
         TTYPE = 1
         L1    = 13
         IF ( STRLEN < 25 ) THEN
            PRINT *, 'Time unit string too short: ' // TRIM(tUnit)
            RC = -999; RETURN
         ENDIF 

      ! Check for 'days since'. If true, set TTYPE to 2 and set the
      ! begin of the effective date string to 11.
      ELSEIF ( MIRRUNIT(1:10) == 'days since' ) THEN
         TTYPE = 2
         L1    = 12
      ELSE
         ! Return w/ error
         PRINT *, 'Invalid time unit: ' // TRIM(tUnit) 
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

      READ( tUnit(L1:L2),'(i)', IOSTAT=STAT ) tYr 
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
      READ( tUnit(L1:L2), '(i)', IOSTAT=STAT ) tMt 
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
      READ( tUnit(L1:L2), '(i)', IOSTAT=STAT ) tDy
      IF ( STAT /= 0 ) THEN
         PRINT *, 'Invalid day in ' // TRIM(tUnit)
         RC = -999; RETURN
      ENDIF

      ! Get reference hour only if 'hours since...'
      IF ( TTYPE == 1 ) THEN

         ! Reference hour
         L1 = L2 + 2
         DO I=L1,STRLEN
            IF(tUnit(I:I) == ':') EXIT 
         ENDDO
         L2 = I-1
         READ( tUnit(L1:L2), '(i)', IOSTAT=STAT ) tHr 
         IF ( STAT /= 0 ) THEN
            PRINT *, 'Invalid hour in ', TRIM(tUnit)
            RC = -999; RETURN
         ENDIF
 
      ELSE 
         ! Set reference hour to -1
         tHr = -1

      ENDIF

      ! Return w/ success
      RC = 0 

      END SUBROUTINE NC_GET_REFDATETIME 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nc_read
!
! !DESCRIPTION: Routine to read emission netCDF data. 
!  The netCDF file must have the variables "lon" and "lat" and,
!  optionally, the variables "lev" and "time" or "date" (in this order).  
!  The (optional) input arguments YEAR, MONTH, DAY and HOUR determine 
!  the time stamp to be read. If not specified, the whole data array
!  (i.e. all time stamps) is read. \\
!  It is possible to only specify some of the time variables. If the year
!  is provided, the other time variables are set to the respective input
!  argument, if provided, and to the begin of the year (Jan 01, 00:00)
!  otherwise.
!  If only the month is provided, the netCDF emission array is assumed to 
!  contain monthly data, and the respective array is used (an error is 
!  returned if the array is not of dimension 12). The same applies for 
!  the time variable HOUR.\\
!  Note that the longitude unit must be "degrees_east", and the 
!  latitude unit must be "degrees_north", otherwise the routine stops 
!  with an error.\\
!  The time format is assumed to be in relative time, i.e. time since a
!  reference date. Currently supported time formats are "days since 
!  YYYY-MM-DD" and "hours since YYYY-MM-DD HH:MN:SS" (with MN and SS 
!  being 00). The reference time doesn't have to be the Geos-Chem
!  reference time (which is, 01/01/1985 00:00:00). However, the
!  calendar is always assumed to be standard, i.e. gregorian calendar 
!  with leapyears being considered. We may want to change that in 
!  future so that other calendars are also supported...
!  The (optional) output argument TIME is always given as hours since 
!  G-C reference time.
!\\
!
! !NOTES: 
!  - For lon, lat and lev, always the whole array is read. 
!    Eventually change in future so that we can also read subarrays. 
!  - Currently, only standard calendar format is supported, i.e. 
!    gregorian calendar (which includes leapyears). Eventually add 
!    support for other calendars in future!
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_READ( FILENAME, PARA,  ARRAY, YEAR, MONTH, &
                          DAY,      HOUR,  UNITS, NLON, NLAT,  & 
                          NLEV,     NTIME, LON,   LAT,  LEV,   &
                          TIME,     RC )
!
! !USES:
!
      USE CHARPAK_MOD,        ONLY : TRANLC

      IMPLICIT NONE
#     include "netcdf.inc"
!
! !ARGUMENTS:
!   
      ! Required input
      CHARACTER(LEN=*), INTENT(IN)            :: FILENAME ! File path 
      CHARACTER(LEN=*), INTENT(IN)            :: PARA     ! Paramater

      ! Optional input
      INTEGER,          INTENT(IN), OPTIONAL  :: YEAR    ! Year to read
      INTEGER,          INTENT(IN), OPTIONAL  :: MONTH   ! Month to read
      INTEGER,          INTENT(IN), OPTIONAL  :: DAY     ! Day to read
      INTEGER,          INTENT(IN), OPTIONAL  :: HOUR    ! Hour to read

      ! Required output: Array to write data
      REAL*4,           POINTER               :: ARRAY(:,:,:,:)

      ! Optional output
      CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: UNITS     ! Unit
      INTEGER,          INTENT(OUT), OPTIONAL :: NLON      ! lon ext.
      INTEGER,          INTENT(OUT), OPTIONAL :: NLAT      ! lat ext.
      INTEGER,          INTENT(OUT), OPTIONAL :: NTIME     ! time ext.
      INTEGER,          INTENT(OUT), OPTIONAL :: NLEV      ! vert. ext.
      REAL*4,           POINTER,     OPTIONAL :: LON(:)    ! Grid lons
      REAL*4,           POINTER,     OPTIONAL :: LAT(:)    ! Grid lats
      REAL*4,           POINTER,     OPTIONAL :: LEV(:)    ! Grid levs
      INTEGER,          POINTER,     OPTIONAL :: TIME(:)   ! Grid time

      ! Error handling
      INTEGER,          INTENT(INOUT)         :: RC

!
! !REMARKS:
!
! !REVISION HISTORY:
!  27 Jul 2012 - C. Keller - Initial version
!  18 Jan 2012 - C. Keller - Now reads 4D, 3D, and 2D arrays, with
!                            optional dimensions level and time.
!  18 Apr 2012 - C. Keller - Now also read offset and scale factor
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
      REAL*8                 :: corr      ! netCDF attribute value 

      ! Arrays for netCDF start and count values
      INTEGER                :: st1d(1), ct1d(1)   ! For 1D arrays    
      INTEGER                :: st2d(2), ct2d(2)   ! For 2D arrays 
      INTEGER                :: st3d(3), ct3d(3)   ! For 3D arrays 
      INTEGER                :: st4d(4), ct4d(4)   ! For 4D arrays 

      ! Grid dimensions
      INTEGER                :: XDIM, YDIM, ZDIM, TDIM

      ! For time index
      CHARACTER(LEN=255)     :: TIMEUNIT, CALENDAR
      INTEGER                :: USEYEAR, USEMONTH, USEDAY, USEHOUR
      INTEGER                :: TIDX, TDIMREAD
      INTEGER                :: TTYPE
      REAL*8                 :: TOFFSET

      ! Temporary 3D array
      REAL*4, ALLOCATABLE    :: TMPARR_3D(:,:,:)
      REAL*4, ALLOCATABLE    :: TMPARR_2D(:,:)

      ! For error handling
      CHARACTER(LEN=255)     :: LOC, MSG
      CHARACTER(LEN=4)       :: YYYY
      CHARACTER(LEN=2)       :: MM, DD, HH
      LOGICAL                :: ReadAtt, ReadTime, ReadLev

      !=================================================================
      ! NC_READ begins here
      !=================================================================

      !-----------------------------------------------------------------
      ! Initialize
      !-----------------------------------------------------------------

      ! For error handling
      LOC = 'NC_READ ("ncdf_mod.F")'

      ! Eventually deallocate output array
      IF ( ASSOCIATED ( ARRAY ) ) DEALLOCATE ( ARRAY )

      ! Default time stamps to read
      USEYEAR  = -1
      USEMONTH = -1
      USEDAY   = -1
      USEHOUR  = -1

      ! Check if time stamps are provided
      IF ( PRESENT(YEAR) ) THEN
         USEYEAR = YEAR
      ENDIF
      IF ( PRESENT(MONTH) ) THEN
         USEMONTH = MONTH
      ENDIF
      IF ( PRESENT(DAY) ) THEN
         USEDAY = DAY
      ENDIF
      IF ( PRESENT(HOUR) ) THEN
         USEHOUR = HOUR
      ENDIF

      !=================================================================
      ! Read netCDF
      !=================================================================

      ! Open netCDF file
      CALL Ncop_Rd( fId, TRIM(FILENAME) )

      !----------------------------------------
      ! VARIABLE: lon
      !----------------------------------------

      ! Variable name
      v_name = "lon"
      a_name = "units"

      ! Read lon units attribute. Longitudes must be in degrees_east!
      CALL NcGet_Var_Attributes( fId, TRIM(v_name), TRIM(a_name),a_val)
      IF ( INDEX ( TRIM(a_val), 'degrees_east' ) == 0 ) THEN
         PRINT *, 'Wrong longitude unit in file ' // TRIM(FILENAME)
         PRINT *, 'units: ' // TRIM(a_val)
         RC = -999; RETURN
      ENDIF

      ! Get longitudinal dimension and allocate longitude vector
      CALL Ncget_Dimlen ( fId, TRIM(v_name), XDIM )

      ! Read lon from file
      ALLOCATE ( longitude(XDIM) )
      st1d   = (/ 1    /)
      ct1d   = (/ XDIM /)
      CALL NcRd( longitude, fId, TRIM(v_name), st1d, ct1d )

      !----------------------------------------
      ! VARIABLE: lat
      !----------------------------------------

      ! Variable name
      v_name = "lat"
      a_name = "units"

      ! Read lat units attribute. Latitudes must be in degrees_north!
      CALL NcGet_Var_Attributes( fId, TRIM(v_name), TRIM(a_name),a_val)
      IF ( INDEX ( TRIM(a_val), 'degrees_north' ) == 0 ) THEN
         PRINT *, 'Wrong latitude unit in file ' // TRIM(FILENAME)
         PRINT *, 'units: ' // TRIM(a_val)
         RC = -999; RETURN
      ENDIF

      ! Get latitudinal dimension and allocate latitude vector
      CALL Ncget_Dimlen ( fId, TRIM(v_name), YDIM )

      ! Read lat from file
      ALLOCATE ( latitude(YDIM) )
      st1d   = (/ 1    /)
      ct1d   = (/ YDIM /)
      CALL NcRd( latitude, fId, TRIM(v_name), st1d, ct1d )

      !----------------------------------------
      ! VARIABLE: lev
      !----------------------------------------

      ! Variable name
      v_name = "lev"

      ! Check if dimension exist
      ReadLev = Ncdoes_Dim_Exist ( fId, TRIM(v_name) ) 

      ! If 'lev' not found, check for 'height'
      IF ( .NOT. ReadLev ) THEN
         v_name = "height"
         ReadLev = Ncdoes_Dim_Exist ( fId, TRIM(v_name) ) 
      ENDIF

      ! Get level dimension and allocate level vector.
      IF ( ReadLev ) THEN
         CALL Ncget_Dimlen ( fId, TRIM(v_name), ZDIM )

         ! Read lev from file
         ALLOCATE ( levels(ZDIM) )
         st1d   = (/ 1    /)
         ct1d   = (/ ZDIM /)
         CALL NcRd( levels, fId, TRIM(v_name), st1d, ct1d )

      ! Set number of levels to one if ncdf has no level variable
      ELSE
         ZDIM = 1
      ENDIF

      !----------------------------------------
      ! VARIABLE: time
      !----------------------------------------

      ! Variable name
      v_name = "time"
      a_name = "units"

      ! Check if dimension "time" exist
      ReadTime = Ncdoes_Dim_Exist ( fId, TRIM(v_name) ) 

      ! If time dim not found, also check for dimension "date"
      IF ( .NOT. ReadTime ) THEN
         v_name   = "date"
         ReadTime = Ncdoes_Dim_Exist ( fId, TRIM(v_name) ) 
      ENDIF

      ! Do the following only if time variable was found in netCDF:
      IF ( ReadTime ) THEN

         ! Get dimension length
         CALL Ncget_Dimlen ( fId, TRIM(v_name), TDIM )

         ! Read time vector from file. 
         ALLOCATE ( timevec(TDIM) )
         st1d = (/ 1    /)
         ct1d = (/ TDIM /)
         CALL NcRd( timevec, fId, TRIM(v_name), st1d, ct1d )

!         ! Read calendar attribute
!         CALL NcGet_Var_Attributes( fId, v_name, 'calendar', 
!                                    CALENDAR )
!    
!         ! Currently, only the standard calendar (i.e. gregorian) is
!         ! supported:
!         CALL TRANLC ( CALENDAR ) ! lower case
!         IF ( ( TRIM(CALENDAR) /= 'gregorian' ) .AND. 
!     &        ( TRIM(CALENDAR) /= 'standard'  )      ) THEN
!            PRINT *, 'Wrong time calendar in ' // TRIM(FILENAME)
!            RC = -999; RETURN
!         ENDIF

         ! Read time/date units attribute
         CALL NcGet_Var_Attributes( fId,          TRIM(v_name), & 
                                    TRIM(a_name), TIMEUNIT )

         ! Check time unit and determine offset of ncdf reference time 
         ! compared to the Geos-Chem reference time.
         ! Time unit must be 'hours since ...' or 'days since ...'
         CALL TIMEUNIT_CHECK( TRIM(TIMEUNIT), TTYPE, TOFFSET, &
              FILENAME, RC )
         IF ( RC /= 0 ) RETURN

         ! Get time index where to start reading and number of time 
         ! indices to be read. Also, convert timevec to 'hours since
         ! G-C reference time'.
         CALL GET_TIDX ( TDIM,    timevec,  TTYPE,  TOFFSET, &
                         USEYEAR, USEMONTH, USEDAY, USEHOUR, &
                         TIDX,    TDIMREAD, RC )
         IF ( RC /= 0 ) RETURN

      ! If time is not defined, set number of time steps to be read
      ! (TDIMREAD) and starting point (TIDX) to 1.
      ELSE
         TIDX     = 1
         TDIMREAD = 1
      ENDIF

      !----------------------------------------
      ! Read array
      !----------------------------------------

      ! Variable name
      v_name = TRIM(PARA)
   
      ! Allocate the array
      ALLOCATE ( ARRAY( XDIM, YDIM, ZDIM, TDIMREAD ) )
   
      ! Read the field from file
   
      ! To read a 4D array:
      IF ( ReadTime .AND. ReadLev ) THEN
         st4d = (/ 1,    1,    1,    TIDX     /)
         ct4d = (/ XDIM, YDIM, ZDIM, TDIMREAD /)
         CALL NcRd( ARRAY, fId, TRIM(v_name), st4d, ct4d )
   
      ! To read a 3D array:
      ! Level defined but not time:
      ELSEIF ( ReadLev ) THEN
         ALLOCATE ( TMPARR_3D( XDIM, YDIM, ZDIM ) )
         st3d = (/ 1,    1,    1    /)
         ct3d = (/ XDIM, YDIM, ZDIM /)
         CALL NcRd( TMPARR_3D, fId, TRIM(v_name), st3d, ct3d )
         ARRAY(:,:,:,1) = TMPARR_3D(:,:,:)

      ! Time defined but not level:
      ELSEIF ( ReadTime ) THEN
         ALLOCATE ( TMPARR_3D( XDIM, YDIM, TDIMREAD) )
         st3d   = (/ 1,    1,    TIDX     /)
         ct3d   = (/ XDIM, YDIM, TDIMREAD /)
         CALL NcRd( TMPARR_3D, fId, TRIM(v_name), st3d, ct3d )
         ARRAY(:,:,1,:) = TMPARR_3D(:,:,:)

      ! Otherwise, we read a 2D array (lon and lat only):
      ELSE
         ALLOCATE ( TMPARR_2D( XDIM, YDIM ) )
         st2d   = (/ 1,    1    /)
         ct2d   = (/ XDIM, YDIM /)
         CALL NcRd( TMPARR_2D, fId, TRIM(v_name), st2d, ct2d )
         ARRAY(:,:,1,1) = TMPARR_2D(:,:)

      ENDIF

      ! ----------------------------
      ! Read optional arguments
      ! ----------------------------

      ! Read units
      IF ( PRESENT(UNITS)  )THEN
         v_name = TRIM(PARA)
         a_name = "units"
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),a_val)
         UNITS = TRIM(a_val)
      ENDIF

      ! ------------------------------------------
      ! Eventually apply scale / offset factors
      ! ------------------------------------------

      ! Check for scale factor
      v_name  = TRIM(PARA)
      a_name  = "scale_factor"
      ReadAtt = Ncdoes_Attr_Exist ( fId, TRIM(v_name), TRIM(a_name) ) 

      IF ( ReadAtt ) THEN
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),corr)
         ARRAY(:,:,:,:) = ARRAY(:,:,:,:) * corr
      ENDIF

      ! Check for offset factor
      v_name  = TRIM(PARA)
      a_name  = "add_offset"
      ReadAtt = Ncdoes_Attr_Exist ( fId, TRIM(v_name), TRIM(a_name) ) 

      IF ( ReadAtt ) THEN
         CALL NcGet_Var_Attributes(fId,TRIM(v_name),TRIM(a_name),corr)
         ARRAY(:,:,:,:) = ARRAY(:,:,:,:) + corr
      ENDIF

      !=================================================================
      ! Prepare (optional) output
      !=================================================================

      ! (Output) array dimensions.
      IF ( PRESENT(NLON ) ) NLON  = XDIM
      IF ( PRESENT(NLAT ) ) NLAT  = YDIM
      IF ( PRESENT(NLEV ) ) NLEV  = ZDIM
      IF ( PRESENT(NTIME) ) NTIME = TDIMREAD

      ! Dimension values

      ! Longitude values
      IF ( PRESENT(LON ) ) THEN
         ALLOCATE( LON( XDIM ) )
         LON = longitude
      ENDIF

      ! Latitude values
      IF ( PRESENT(LAT) ) THEN
         ALLOCATE( LAT( YDIM ) )
         LAT = latitude
      ENDIF

      ! Horizontal level values
      IF ( PRESENT(LEV) ) THEN

         ! If defined
         IF ( ReadLev ) THEN
            ALLOCATE( LEV( ZDIM ) )
            LEV = levels

         ! Default       
         ELSE
            ALLOCATE( LEV(1) )
            LEV = -999d0
         ENDIF
      ENDIF

      ! Time values
      IF ( PRESENT(TIME) ) THEN

         ! If defined
         IF ( ReadTime ) THEN
            ALLOCATE( TIME( TDIMREAD ) )

            ! Get the specific time stamp that has been used:
            IF ( TDIMREAD == 1 ) THEN
               TIME = timevec( TIDX )

            ! Get entire time vector
            ELSE
               TIME  = timevec
            ENDIF

         ! Default
         ELSE
            ALLOCATE( TIME(1) )
            TIME = -999
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
      IF ( ALLOCATED ( TMPARR_3D ) ) DEALLOCATE ( TMPARR_3D )
      IF ( ALLOCATED ( TMPARR_2D ) ) DEALLOCATE ( TMPARR_2D )

      ! Close netCDF file
      CALL NcCl( fId )

      ! Return w/ success
      RC = 0 

      END SUBROUTINE NC_READ
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_tidx
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
! !USES:
!
      USE BPCH2_MOD, ONLY : GET_TAU0
!
! !ARGUMENTS:
!   
      ! Required
      INTEGER, INTENT(   IN)           :: TDIM
      INTEGER, INTENT(INOUT)           :: TIMEVEC(TDIM)
      INTEGER, INTENT(   IN)           :: TTYPE
      REAL*8,  INTENT(   IN)           :: TOFFSET
      INTEGER, INTENT(INOUT)           :: YEAR
      INTEGER, INTENT(INOUT)           :: MONTH
      INTEGER, INTENT(INOUT)           :: DAY
      INTEGER, INTENT(INOUT)           :: HOUR
      INTEGER, INTENT(  OUT)           :: TIDX
      INTEGER, INTENT(  OUT)           :: TDIMREAD
      INTEGER, INTENT(INOUT)           :: RC
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: timeunit_check
!
! !DESCRIPTION: Makes a validity check of the passed unit string.
! Supported formats are "days since YYYY-MM-DD" (TIMETYPE=1) and 
! "hours since YYYY-MM-DD HH:MM:SS" (TIMETYPE=2).\\
! The output argument TOFFSET gives the offset of the ncdf reference 
! time relative to Geos-Chem reference time (in hours).
!
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE TIMEUNIT_CHECK( TIMEUNIT, TIMETYPE, TOFFSET, & 
                                 FILENAME, RC )
!
! !USES:
!
      USE BPCH2_MOD,          ONLY : GET_TAU0
      USE CHARPAK_MOD,        ONLY : TRANLC
!
! !ARGUMENTS:
!   
      ! Required
      CHARACTER(LEN=*), INTENT( IN)    :: TIMEUNIT
      INTEGER,          INTENT(OUT)    :: TIMETYPE
      REAL*8,           INTENT(OUT)    :: TOFFSET
      CHARACTER(LEN=*), INTENT( IN)    :: FILENAME
      INTEGER,          INTENT(INOUT)  :: RC 
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
      READ( TIMEUNIT(L1:L2),'(i)', IOSTAT=STAT ) YYYY
      IF ( STAT /= 0 ) THEN
         PRINT *, 'Invalid year in ', TRIM(TIMEUNIT), &
               ' in file'             , TRIM(FILENAME)
         RC = -999; RETURN
      ENDIF 

      ! Read reference month. Typically, the month is represented by
      ! two characters, i.e. 1 is 01, etc.
      L1 = L2 + 2
      L2 = L1 + 1
      READ( TIMEUNIT(L1:L2), '(i)', IOSTAT=STAT ) MM
      ! Also check for the case where the month is only one character:
      IF ( STAT /= 0 ) THEN
         L2 = L1
         READ( TIMEUNIT(L1:L2), '(i)', IOSTAT=STAT ) MM
         IF ( STAT /= 0 ) THEN
            PRINT *, 'Invalid month in ', TRIM(TIMEUNIT), &
                   ' in file'             , TRIM(FILENAME)
            RC = -999; RETURN
         ENDIF 
      ENDIF

      ! Reference day. Typically, the day is represented by two 
      ! characters, i.e. 1 is 01, etc.
      L1 = L2 + 2
      L2 = L1 + 1
      READ( TIMEUNIT(L1:L2), '(i)', IOSTAT=STAT ) DD
      ! Also check for the case where the day is only one character:
      IF ( STAT /= 0 ) THEN
         L2 = L1
         READ( TIMEUNIT(L1:L2), '(i)', IOSTAT=STAT ) DD
         IF ( STAT /= 0 ) THEN
            PRINT *, 'Invalid day in ', TRIM(TIMEUNIT), &
                   ' in file'           , TRIM(FILENAME)
            RC = -999; RETURN
         ENDIF 
      ENDIF

      ! Get reference hour only if 'hours since...'
      IF ( TTYPE == 1 ) THEN

         ! Reference hour
         L1 = L2 + 2
         L2 = L1 + 1
         READ( TIMEUNIT(L1:L2), '(i)', IOSTAT=STAT ) HH
         IF ( STAT /= 0 ) THEN
            L2 = L1
            READ( TIMEUNIT(L1:L2), '(i)', IOSTAT=STAT ) HH
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nc_read_grid
!
! !DESCRIPTION: Routine to read variables and attributes from a netCDF
!  file.  This routine was automatically generated by the Perl script
!  NcdfUtilities/perl/ncCodeRead.
!  If the grid edges cannot be read from the netCDF file, the grid
!  centers are read instead and the grid edges are calculated from 
!  these centers.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_READ_GRID( IM,        JM,        fileName, &
                               lon_edges, lat_sines, RC         ) 
!
! !USES:
!
      USE CMN_GCTM_MOD, ONLY : PI_180

      IMPLICIT NONE

#   include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
      INTEGER,          INTENT(IN)  :: IM                ! # of longitudes
      INTEGER,          INTENT(IN)  :: JM                ! # of latitudes
      CHARACTER(LEN=*), INTENT(IN)  :: fileName          ! File w/ grid info
!
! !OUTPUT PARAMETERS:
!   
      REAL*4,           INTENT(OUT) :: lon_edges(IM+1)   ! Lon edges [degrees]
      REAL*4,           INTENT(OUT) :: lat_sines(JM+1)   ! SIN( latitude edges )
      INTEGER,          INTENT(INOUT)  :: RC 
!
! !REMARKS:
!  Created with the ncCodeRead script of the NcdfUtilities package,
!  with subsequent hand-editing.
!
! !REVISION HISTORY:
!  23 Aug 2012 - R. Yantosca - Initial version
!  24 Jan 2012 - C. Keller   - Added grid edge calculation from mid points
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER            :: fId                          ! netCDF file ID
      INTEGER            :: STAT                         ! error status          
      INTEGER            :: II, JJ 

      ! Arrays
      INTEGER            :: st1d(1), ct1d(1)             ! netCDF start & count

      ! For calculating the edges from the mid-points:
      REAL*4             :: lon_mid(IM)                  ! longitude midoints
      REAL*4             :: lat_mid(JM)                  ! latitude midpoints
      LOGICAL            :: REDO

      ! For error handling
      CHARACTER(LEN=255) :: MSG, LOC

      !======================================================================
      ! NC_READ_GRID begins here
      !======================================================================

      ! For error handling
      LOC = 'NC_READ_GRID (ncdf_mod.F90)' 

      ! Open file for reading
      CALL Ncop_Rd( fId, TRIM( fileName ) )

      !----------------------------------------
      ! Longitude edges
      !----------------------------------------

      ! Try to read lon_edges directly from netCDF
      st1d = (/ 1    /)
      ct1d = (/ IM+1 /)
      CALL NcRd( lon_edges, fId,  "lon_edges", st1d, ct1d, & 
                 err_stop=.FALSE., stat=STAT )

      ! If reading of lon_edges failed, calculate them from the longitude 
      ! centers. These values are assumed to be in degrees_east.
      IF ( STAT /= 0 ) THEN
         lon_mid(:) = 0d0
         st1d = (/ 1  /)
         ct1d = (/ IM /)
         CALL NcRd( lon_mid, fId,  "lon", st1d, ct1d, &           
                    err_stop=.FALSE., stat=STAT )
       
         ! Return w/ error if reading failed
         IF ( STAT /= 0 ) THEN
            PRINT *, 'Cannot calculate longitude edges!', TRIM(fileName)
            RC = -999; RETURN
         ENDIF

         ! Longitude values must be steadily increasing.
         REDO = .TRUE.
         DO JJ = 1, 10 ! Try max. 10 times

            ! Leave if data is steadily increasing
            IF ( .NOT. REDO ) EXIT

            ! Exit after 10 iterations
            IF ( JJ == 10 ) THEN
               PRINT *, 'Longitudes are not steadily increasing!', &
                  TRIM(fileName)
               RC = -999; RETURN
            ENDIF

            ! Loop over all lon values
            DO II = 1, IM

               ! If we reach the last grid box, all values are steadily
               ! increasing (otherwise, the II loop would have been
               ! interrupted and the check would have been repeaded).
               IF ( II == IM ) THEN
                  REDO = .FALSE.
                  EXIT
               ENDIF

               ! Check if next longitude value is lower, in which case 
               ! we subtract 360 degrees from all longitude values up
               ! to this point. Exit the II loop and repeat the check.
               IF ( lon_mid(II+1) < lon_mid(II) ) THEN
                  lon_mid(1:II) = lon_mid(1:II) - 360d0
                  REDO = .TRUE.
                  EXIT
               ENDIF

            ENDDO !II
         ENDDO !JJ

         ! Get leftmost lon edge by extrapolating from first two 
         ! midpoints.
         lon_edges(1) = lon_mid(1) - & 
                        ( (lon_mid(2) - lon_mid(1) ) / 2d0 )
       
         ! Sequentially calculate the right edge from the previously 
         ! calculated left edge.
         DO II = 1, IM
            lon_edges(II+1) = lon_mid(II) + lon_mid(II) - lon_edges(II)
         ENDDO 
      ENDIF

      !----------------------------------------
      ! Latitude edges
      !----------------------------------------

      ! Try to read lat_sines directly from netCDF
      st1d = (/ 1    /)
      ct1d = (/ JM+1 /)
      CALL NcRd( lat_sines, fId,  "lat_sines", st1d, ct1d, & 
                 err_stop=.FALSE., stat=STAT )

      ! If reading of lat_sines failed, try to read lat_edges and
      ! calculate the sin of them
      IF ( STAT /= 0 ) THEN

         ! Try to read lat_edges
         CALL NcRd( lat_sines, fId,  "lat_edges", st1d, ct1d, &
                    err_stop=.FALSE., stat=STAT )

         ! If successful, calculate the sines of the edges
         IF ( STAT == 0 ) THEN
            DO II = 1, JM+1
               lat_sines(II) = SIN( lat_sines(II) * PI_180 )
            ENDDO
         ENDIF
      ENDIF

      ! If both of the above failed, read latitude centers and calculate
      ! lat edges from it
      IF ( STAT /= 0 ) THEN
         st1d = (/ 1  /)
         ct1d = (/ JM /)
         CALL NcRd( lat_mid, fId,  "lat", st1d, ct1d, &           
                    err_stop=.FALSE., stat=STAT )

         ! Return w/ error if reading failed
         IF ( STAT /= 0 ) THEN
            PRINT *, 'Cannot calculate latitude edges!', TRIM(fileName)
            RC = -999; RETURN
         ENDIF

         ! Get first edge by extrapolating from first two midpoints.
         lat_sines(1) = lat_mid(1) - & 
                        ( (lat_mid(2) - lat_mid(1) ) / 2d0 )
       
         ! If first edge is less than -90 degrees, manually adjust it
         IF ( lat_sines(1) < -90d0 ) THEN
            lat_sines(1) = -90d0
         ENDIF 

         ! Sequentially calculate the upper lat edge from the previously 
         ! calculated lower edge. Also calculate the sine of the latitude. 
         DO II = 1, JM
            lat_sines(II+1) = lat_mid(II) + lat_mid(II) - lat_sines(II)
            lat_sines(II)   = SIN( lat_sines(II) * PI_180 )
         ENDDO 

         ! The last edge must not exceed 90 deg north
         IF ( lat_sines(JM+1) > 90.01d0 ) THEN
            PRINT *, 'Uppermost latitude edge above 90 deg north!', &
               TRIM(fileName)
            PRINT *, lat_sines
            RC = -999; RETURN
         ELSE
            lat_sines(JM+1) = SIN( lat_sines(JM+1) * PI_180 )
         ENDIF
      ENDIF

      ! Close netCDF file
      CALL NcCl( fId )

      ! Return w/ success
      RC = 0 

      END SUBROUTINE NC_READ_GRID
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nc_write_3d
!
! !DESCRIPTION: Routine to write time slices of 2D fields into netCDF. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_WRITE_3D (ncFile,  I,  J,    T,  N,   lon, lat, &
                              time,    timeUnit, ncVars,  ncUnits,  &
                              ncLongs, ncShorts, ncArrays            )
!
! !USES:
!
      IMPLICIT NONE

#   include "netcdf.inc"
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: nc_write_4d
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
! !USES:
!
      IMPLICIT NONE

#   include "netcdf.inc"
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NC_DEFINE 
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
! !USES:
!
      IMPLICIT NONE

#     include "netcdf.inc"
!
! !OUTPUT PARAMETERS:
! 
      CHARACTER(LEN=*),   INTENT(IN   )  :: ncFile   ! ncdf file path + name 
      INTEGER,            INTENT(IN   )  :: nLon     ! # of lons 
      INTEGER,            INTENT(IN   )  :: nLat     ! # of lats
      INTEGER, OPTIONAL,  INTENT(IN   )  :: nLev     ! # of levels
      INTEGER,            INTENT(IN   )  :: nTime    ! # of time stamps
      CHARACTER(LEN=*),   INTENT(IN   )  :: timeUnit ! time unit 
      CHARACTER(LEN=*),   INTENT(IN   )  :: ncVars(:)    ! ncdf variables
      CHARACTER(LEN=*),   INTENT(IN   )  :: ncUnits(:)   ! var units
      CHARACTER(LEN=*),   INTENT(IN   )  :: ncLongs(:)   ! var long names
      CHARACTER(LEN=*),   INTENT(IN   )  :: ncShorts(:)  ! var short names
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
      CALL NcDef_Dimension( fId, TRIM(v_name), ntime, id_time )

      !--------------------------------
      ! VARIABLE: lon
      !--------------------------------

      ! Define the "lon" variable
      v_name = "lon"
      vId    = vId + 1
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
      vId    = vId + 1
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
         vId    = vId + 1
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
      vId    = vId + 1
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
         vId    = vId + 1
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_to_netcdf_file
!
! !DESCRIPTION: Routine to write data to a netCDF file.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_WRITE_DIMS( fID, lon, lat, time, lev ) 
!
! !USES:
!
      IMPLICIT NONE

#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
      INTEGER,          INTENT(INOUT )           :: fId
      REAL*4,           INTENT(IN    )           :: lon(:) 
      REAL*4,           INTENT(IN    )           :: lat(:) 
      REAL*4,           INTENT(IN    )           :: time(:) 
      REAL*4, OPTIONAL, INTENT(IN    )           :: lev(:) 
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
      ! Character strings
      CHARACTER(LEN=255) :: v_name             ! netCDF variable name
      
      ! Arrays for netCDF start and count values
      INTEGER            :: st1d(1), ct1d(1)   ! For 1D arrays

      !=================================================================
      ! Define lon/lat
      !=================================================================

      ! Write lon to netCDF file
      st1d   = (/ 1 /)
      ct1d   = (/ size(lon,1) /)
      v_name = "lon"
      CALL NcWr( lon, fId, TRIM(v_name), st1d, ct1d )

      ! Write lat to netCDF file
      st1d   = (/ 1 /)
      ct1d   = (/ size(lat,1) /)
      v_name = "lat"
      CALL NcWr( lat, fId, TRIM(v_name), st1d, ct1d )

      ! Write lev to netCDF file
      IF ( PRESENT(lev) ) THEN
         st1d   = (/ 1 /)
         ct1d   = (/ size(lev,1) /)
         v_name = "lev"
         CALL NcWr( lev, fId, TRIM(v_name), st1d, ct1d )
      ENDIF

      ! Write passed time integer to netCDF file
      st1d   = (/ 1 /)
      ct1d   = (/ size(time,1) /)
      v_name = "time"
      CALL NcWr( time, fId, TRIM(v_name), st1d, ct1d )

      END SUBROUTINE NC_WRITE_DIMS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_to_netcdf_file
!
! !DESCRIPTION: Routine to write data to a netCDF file.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_WRITE_DATA_3D ( fID, ncVar, Array )
!
! !USES:
!
      IMPLICIT NONE

#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
      INTEGER,          INTENT(INOUT)           :: fId
      CHARACTER(LEN=*), INTENT(IN   )           :: ncVar
      REAL*4,           POINTER                 :: Array(:,:,:)
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

      st3d   = (/ 1, 1, 1 /)
      ct3d   = (/ size(array,1), size(array,2), size(array,3) /)
      CALL NcWr( ARRAY, fId, TRIM(ncVar), st3d, ct3d )

      END SUBROUTINE NC_WRITE_DATA_3D
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_to_netcdf_file
!
! !DESCRIPTION: Routine to write data to a netCDF file.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_WRITE_DATA_4D ( fID, ncVar, Array )
!
! !USES:
!
      IMPLICIT NONE

#     include "netcdf.inc"
!
! !INPUT PARAMETERS:
!
      INTEGER,          INTENT(INOUT)           :: fId
      CHARACTER(LEN=*), INTENT(IN   )           :: ncVar
      REAL*4,           POINTER                 :: Array(:,:,:,:)
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
!BOP
!
! !IROUTINE: NC_CREATE
!
! !DESCRIPTION: creates a new netCDF file. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_CREATE( NcFile, nLon,  nLat,  nLev,  nTime, &
                            fId,    lonID, latId, levId, timeId, VarCt )
!
! !USES:
!
      IMPLICIT NONE

#     include "netcdf.inc"
!
! !OUTPUT PARAMETERS:
! 
      CHARACTER(LEN=*),   INTENT(IN   )  :: ncFile   ! ncdf file path + name 
      INTEGER,            INTENT(IN   )  :: nLon     ! # of lons 
      INTEGER,            INTENT(IN   )  :: nLat     ! # of lats 
      INTEGER,            INTENT(IN   )  :: nLev     ! # of levs 
      INTEGER,            INTENT(IN   )  :: nTime    ! # of times 
      INTEGER,            INTENT(  OUT)  :: fId      ! file id 
      INTEGER,            INTENT(  OUT)  :: lonId    ! lon id 
      INTEGER,            INTENT(  OUT)  :: latId    ! lat id 
      INTEGER,            INTENT(  OUT)  :: levId    ! lev id 
      INTEGER,            INTENT(  OUT)  :: timeId   ! time id 
      INTEGER,            INTENT(  OUT)  :: VarCt    ! variable counter 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: omode

      ! Open filename
      CALL NcCr_Wr( fId, TRIM(ncFile) )

      ! Turn filling off
      CALL NcSetFill( fId, NF_NOFILL, omode )     

      !--------------------------------
      ! SET GLOBAL ATTRIBUTES
      !--------------------------------
      CALL NcDef_Glob_Attributes( fId, 'title',   'HEMCO diagnostics' ) 
      CALL NcDef_Glob_Attributes( fId, 'history', 'NC_CREATE.F90'     ) 
      CALL NcDef_Glob_Attributes( fId, 'format',  'netCDF-3'          )

      !--------------------------------
      ! SET DIMENSIONS
      !--------------------------------
      CALL NcDef_Dimension ( fId, 'lon',  nLon,  lonId  ) 
      CALL NcDef_Dimension ( fId, 'lat',  nLat,  latId  ) 
      IF ( nLev > 0 ) THEN
         CALL NcDef_Dimension ( fId, 'lev',  nLev,  levId  )
      ELSE
         levId = -1
      ENDIF 
      CALL NcDef_Dimension ( fId, 'time', nTime, TimeId ) 

      ! Close definition section
      CALL NcEnd_Def( fId )

      ! Initialize variable counter
      VarCt = -1

      END SUBROUTINE NC_CREATE
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NC_VAR_DEF
!
! !DESCRIPTION: defines a new variable. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_VAR_DEF ( fId, lonId, latId, levId, TimeId, &
                              VarName, VarLongName, VarUnit,    &
                              DataType, VarCt )
!
! !USES:
!
      IMPLICIT NONE

#     include "netcdf.inc"
!
! !OUTPUT PARAMETERS:
! 
      INTEGER,            INTENT(IN   )  :: fId      ! file ID 
      INTEGER,            INTENT(IN   )  :: lonId 
      INTEGER,            INTENT(IN   )  :: latId 
      INTEGER,            INTENT(IN   )  :: levId 
      INTEGER,            INTENT(IN   )  :: TimeId 
      CHARACTER(LEN=*),   INTENT(IN   )  :: VarName
      CHARACTER(LEN=*),   INTENT(IN   )  :: VarLongName
      CHARACTER(LEN=*),   INTENT(IN   )  :: VarUnit
      INTEGER,            INTENT(IN   )  :: DataType ! 1=Int, 4=float, 8=double 
      INTEGER,            INTENT(INOUT)  :: VarCt    ! variable counter 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, ALLOCATABLE :: VarDims(:) 
      INTEGER              :: nDim, Pos
      INTEGER              :: NF_TYPE

      !--------------------------------
      ! DEFINE VARIABLE 
      !--------------------------------

      ! Reopen definition section
      CALL NcBegin_Def( fId )

      ! Increate variable counter
      VarCt = VarCt + 1
    
      ! number of dimensions
      nDim = 0
      if ( lonId  >= 0 ) nDim = nDim + 1
      if ( latId  >= 0 ) nDim = nDim + 1
      if ( levId  >= 0 ) nDim = nDim + 1
      if ( timeId >= 0 ) nDim = nDim + 1

      ! write dimensions
      allocate( VarDims(nDim) )
      Pos = 1
      if ( lonId >= 0 ) then
         VarDims(Pos) = lonId
         Pos          = Pos + 1
      endif
      if ( latId >= 0 ) then
         VarDims(Pos) = latId
         Pos          = Pos + 1
      endif
      if ( levId >= 0 ) then
         VarDims(Pos) = levId
         Pos          = Pos + 1
      endif
      if ( timeId >= 0 ) then
         VarDims(Pos) = timeId
         Pos          = Pos + 1
      endif

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
 
      ! Define variable
      CALL NcDef_Variable( fId,  TRIM(VarName), NF_TYPE,   &
                           nDim, VarDims,       VarCt       )
      deallocate( VarDims )

      ! Define variable atttibutes
      CALL NcDef_Var_Attributes( fId, VarCt, 'long_name', TRIM(VarLongName) )
      CALL NcDef_Var_Attributes( fId, VarCt, 'units',     TRIM(VarUnit    ) )

      ! Close definition section
      CALL NcEnd_Def( fId )

      END SUBROUTINE NC_VAR_DEF
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NC_VAR_WRITE_R8
!
! !DESCRIPTION: writes data of a double precision variable. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_VAR_WRITE_R8 ( fId, VarName, Arr1D, Arr2D, Arr3D, Arr4D )
!
! !USES:
!
      IMPLICIT NONE

#     include "netcdf.inc"
!
! !OUTPUT PARAMETERS:
! 
      INTEGER,          INTENT(IN)           :: fId            ! file ID 
      CHARACTER(LEN=*), INTENT(IN)           :: VarName        ! variable name      
      REAL(kind=8),     POINTER,   OPTIONAL  :: Arr1D(:)       ! array to be written 
      REAL(kind=8),     POINTER,   OPTIONAL  :: Arr2D(:,:)     ! array to be written 
      REAL(kind=8),     POINTER,   OPTIONAL  :: Arr3D(:,:,:)   ! array to be written 
      REAL(kind=8),     POINTER,   OPTIONAL  :: Arr4D(:,:,:,:) ! array to be written 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, ALLOCATABLE :: St(:), Ct(:)
      INTEGER              :: I, nDim

      !--------------------------------
      ! WRITE DATA 
      !--------------------------------
   
      ! 1D data 
      if ( present(Arr1d) ) then
         nDim = 1
         allocate( St(nDim), Ct(nDim) ) 
         St(1) = 1
         Ct(1) = size(Arr1d,1)
         CALL NcWr( Arr1d, fId, trim(VarName), St, Ct )

      ! 2D data 
      elseif ( present(arr2d) ) then
         nDim = 2
         allocate( St(nDim), Ct(nDim) ) 
         do i=1,nDim
            St(i) = 1
            Ct(i) = size(Arr2d,i)
         enddo
         CALL NcWr( Arr2d, fId, trim(VarName), St, Ct )

      ! 3D data
      elseif ( present(arr3d) ) then
         nDim = 3
         allocate( St(nDim), Ct(nDim) ) 
         do i=1,nDim
            St(i) = 1
            Ct(i) = size(Arr3d,i)
         enddo
         CALL NcWr( Arr3d, fId, trim(VarName), St, Ct )

      ! 4D data
      elseif ( present(arr4d) ) then
         nDim = 4
         allocate( St(nDim), Ct(nDim) ) 
         do i=1,nDim
            St(i) = 1
            Ct(i) = size(Arr4d,i)
         enddo
         CALL NcWr( Arr4d, fId, trim(VarName), St, Ct )
      endif

      ! cleanup
      if(allocated(St)) deallocate( St )
      if(allocated(Ct)) deallocate( Ct )
 
      END SUBROUTINE NC_VAR_WRITE_R8
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NC_VAR_WRITE_R4
!
! !DESCRIPTION: writes data of a single precision variable. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_VAR_WRITE_R4 ( fId, VarName, Arr1D, Arr2D, Arr3D, Arr4D )
!
! !USES:
!
      IMPLICIT NONE

#     include "netcdf.inc"
!
! !OUTPUT PARAMETERS:
! 
      INTEGER,          INTENT(IN)           :: fId            ! file ID 
      CHARACTER(LEN=*), INTENT(IN)           :: VarName        ! variable name      
      REAL(kind=4),     POINTER,   OPTIONAL  :: Arr1D(:)       ! array to be written 
      REAL(kind=4),     POINTER,   OPTIONAL  :: Arr2D(:,:)     ! array to be written 
      REAL(kind=4),     POINTER,   OPTIONAL  :: Arr3D(:,:,:)   ! array to be written 
      REAL(kind=4),     POINTER,   OPTIONAL  :: Arr4D(:,:,:,:) ! array to be written 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, ALLOCATABLE :: St(:), Ct(:)
      INTEGER              :: I, nDim

      !--------------------------------
      ! WRITE DATA 
      !--------------------------------
   
      ! 1D data 
      if ( present(Arr1d) ) then
         nDim = 1
         allocate( St(nDim), Ct(nDim) ) 
         St(1) = 1
         Ct(1) = size(Arr1d,1)
         CALL NcWr( Arr1d, fId, trim(VarName), St, Ct )

      ! 2D data 
      elseif ( present(arr2d) ) then
         nDim = 2
         allocate( St(nDim), Ct(nDim) ) 
         do i=1,nDim
            St(i) = 1
            Ct(i) = size(Arr2d,i)
         enddo
         CALL NcWr( Arr2d, fId, trim(VarName), St, Ct )

      ! 3D data
      elseif ( present(arr3d) ) then
         nDim = 3
         allocate( St(nDim), Ct(nDim) ) 
         do i=1,nDim
            St(i) = 1
            Ct(i) = size(Arr3d,i)
         enddo
         CALL NcWr( Arr3d, fId, trim(VarName), St, Ct )

      ! 4D data
      elseif ( present(arr4d) ) then
         nDim = 4
         allocate( St(nDim), Ct(nDim) ) 
         do i=1,nDim
            St(i) = 1
            Ct(i) = size(Arr4d,i)
         enddo
         CALL NcWr( Arr4d, fId, trim(VarName), St, Ct )
      endif

      ! cleanup
      if(allocated(St)) deallocate( St )
      if(allocated(Ct)) deallocate( Ct )
 
      END SUBROUTINE NC_VAR_WRITE_R4
!EOC
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: NC_VAR_WRITE_INT
!
! !DESCRIPTION: writes data of an integer variable. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NC_VAR_WRITE_INT ( fId, VarName, Arr1D, Arr2D, Arr3D, Arr4D )
!
! !USES:
!
      IMPLICIT NONE

#     include "netcdf.inc"
!
! !OUTPUT PARAMETERS:
! 
      INTEGER,          INTENT(IN)           :: fId            ! file ID 
      CHARACTER(LEN=*), INTENT(IN)           :: VarName        ! variable name      
      INTEGER,          POINTER,   OPTIONAL  :: Arr1D(:)       ! array to be written 
      INTEGER,          POINTER,   OPTIONAL  :: Arr2D(:,:)     ! array to be written 
      INTEGER,          POINTER,   OPTIONAL  :: Arr3D(:,:,:)   ! array to be written 
      INTEGER,          POINTER,   OPTIONAL  :: Arr4D(:,:,:,:) ! array to be written 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER, ALLOCATABLE :: St(:), Ct(:)
      INTEGER              :: I, nDim

      !--------------------------------
      ! WRITE DATA 
      !--------------------------------
   
      ! 1D data 
      if ( present(Arr1d) ) then
         nDim = 1
         allocate( St(nDim), Ct(nDim) ) 
         St(1) = 1
         Ct(1) = size(Arr1d,1)
         CALL NcWr( Arr1d, fId, trim(VarName), St, Ct )

      ! 2D data 
      elseif ( present(arr2d) ) then
         nDim = 2
         allocate( St(nDim), Ct(nDim) ) 
         do i=1,nDim
            St(i) = 1
            Ct(i) = size(Arr2d,i)
         enddo
         CALL NcWr( Arr2d, fId, trim(VarName), St, Ct )

      ! 3D data
      elseif ( present(arr3d) ) then
         nDim = 3
         allocate( St(nDim), Ct(nDim) ) 
         do i=1,nDim
            St(i) = 1
            Ct(i) = size(Arr3d,i)
         enddo
         CALL NcWr( Arr3d, fId, trim(VarName), St, Ct )

      ! 4D data
      elseif ( present(arr4d) ) then
         nDim = 4
         allocate( St(nDim), Ct(nDim) ) 
         do i=1,nDim
            St(i) = 1
            Ct(i) = size(Arr4d,i)
         enddo
         CALL NcWr( Arr4d, fId, trim(VarName), St, Ct )
      endif

      ! cleanup
      if(allocated(St)) deallocate( St )
      if(allocated(Ct)) deallocate( Ct )
 
      END SUBROUTINE NC_VAR_WRITE_INT
!EOC
      END MODULE NCDF_MOD
!EOM
