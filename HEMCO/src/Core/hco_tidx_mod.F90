!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_tidx_mod.F90
!
! !DESCRIPTION: Module HCO\_tIdx\_Mod contains routines and variables
! to organize and index data array time slices.
!\\
!\\
! The HEMCO data containers can hold multiple 2D or 3D data arrays
! (aligned in a vector), providing a 4th dimension (time). During
! emission calculation, only the vector element ('time slice')
! representative for the given time will be used. Currently, the
! following time slice intervals are supported:
!
! \begin{itemize}
! \item Constant: Only one time slice (default)
! \item Hourly: Between 2-24 time slices. These slices will be split
!         into even day bins, and are cycled through accordingly.
!         the local time is used to obtain the current valid index at
!         each longitude.
! \item Hourly\_gridded: As hourly, but uses the same time slice index
!         across all longitudes (based upon UTC time).
! \item Weekdaily: Seven time slices, representing the days of the
!         week: Sun, Mon, ..., Sat. Uses local time.
! \item Monthly: 12 time slices, representing the months of the year:
!         Jan, ..., Dec. Uses local time.
! \end{itemize}
!
! The time slice cycling frequency is automatically determined during
! creation of a data container - based upon the time stamp settings in
! the HEMCO configuration file and the time information read from the
! data (netCDF) file.
!\\
!\\
! Spatial uniform data (i.e. nx = ny = 1) is always assigned the local
! time cycle intervals (hourly, weekdaily, or monthly), e.g. the local
! time is used at every grid box when picking the time slice at a given
! time. For gridded data, it's assumed that local-time effects are already
! taken into account and UTC time is used at all locations to select the
! currently valid time slice. The exception is weekdaily data, which is
! always assumed to be in local time.
!\\
!\\
! Structure AlltIDx organizes the indexing of the vector arrays. It
! contains the current valid time slice index for each of the above
! defined time slice intervals. Each data container points to one of the
! elements of AlltIDx, according to the temporal dimension of its array.
! The values of AlltIDx become update on every HEMCO time step.
!\\
!\\
!
! !INTERFACE:
!
MODULE HCO_tIdx_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Types_Mod,  ONLY : TimeIdx
  USE HCO_Types_Mod,  ONLY : TimeIdxCollection

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: tIDx_Assign
  PUBLIC :: tIDx_Init
  PUBLIC :: tIDx_GetIndx
  PUBLIC :: tIDx_Cleanup
  PUBLIC :: tIDx_IsInRange
  PUBLIC :: HCO_GetPrefTimeAttr
  PUBLIC :: HCO_ExtractTime
!
! !REMARKS:
!  The current local time implementation assumes a regular grid,
!  i.e. local time does not change with latitude!
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller   - Initialization
!  22 Aug 2013 - C. Keller   - Some time slice updates.
!  08 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  08 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  03 Dec 2014 - C. Keller   - Major update: now calculate the time slice
!                              indeces on the fly instead of storing them in
!                              precalculated vectors.
!  25 Feb 2015 - R. Yantosca - Comment out WEEKDAY_GRID, it is not used
!                              anymore.  This avoids seg faults.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: tIDx_Init
!
! !DESCRIPTION: Subroutine tIDx\_Init initializes the time slice index
! collection.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE tIDx_Init( HcoState, RC )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState  ! Hemco state
    INTEGER,         INTENT(INOUT)  :: RC        ! Return code
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! tIDx_Init begins here!
    !======================================================================

    ! Enter
    CALL HCO_ENTER ( HcoState%Config%Err, 'tIDx_Init (hco_tidx_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Allocate collection of time indeces
    ALLOCATE ( HcoState%AlltIDx )

    ! Initialize the vectors holding the currently valid time slice
    ! indices for the various cycle intervals. Only create longitude-
    ! dependent time slice vector if cycle interval changes with
    ! longitude, i.e. if it is time-zone dependent.

    ! ----------------------------------------------------------------
    ! "CONSTANT" => only one time slice, no cycling.
    ! ----------------------------------------------------------------
    ALLOCATE ( HcoState%AlltIDx%CONSTANT )
    HcoState%AlltIDx%CONSTANT%TypeID       = 0
    HcoState%AlltIDx%CONSTANT%TempRes      = "Constant"

    ! ----------------------------------------------------------------
    ! "HOURLY" => changes every hour, longitude-dependent
    ! ----------------------------------------------------------------
    ALLOCATE ( HcoState%AlltIDx%HOURLY )
    HcoState%AlltIDx%HOURLY%TypeID       = 24
    HcoState%AlltIDx%HOURLY%TempRes      = "Hourly"

    ! ----------------------------------------------------------------
    ! "HOURLY_GRID" => changes every hour, longitude-independent
    ! ----------------------------------------------------------------
    ALLOCATE ( HcoState%AlltIDx%HOURLY_GRID )
    HcoState%AlltIDx%HOURLY_GRID%TypeID       = 241
    HcoState%AlltIDx%HOURLY_GRID%TempRes      = "Hourly_Grid"

    ! ----------------------------------------------------------------
    ! "WEEKDAY" => changes every weekday, longitude-dependent
    ! ----------------------------------------------------------------
    ALLOCATE ( HcoState%AlltIDx%WEEKDAY )
    HcoState%AlltIDx%WEEKDAY%TypeID       = 7
    HcoState%AlltIDx%WEEKDAY%TempRes      = "Weekday"

    ! ----------------------------------------------------------------
    ! "MONTHLY" => changes every month, longitude-dependent
    ! ----------------------------------------------------------------
    ALLOCATE ( HcoState%AlltIDx%MONTHLY )
    HcoState%AlltIDx%MONTHLY%TypeID       = 12
    HcoState%AlltIDx%MONTHLY%TempRes      = "Monthly"

    ! Return w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

  END SUBROUTINE tIDx_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: tIDx_Set
!
! !DESCRIPTION: Subroutine tIDx\_Set linkes the passed TimeIDx type to the
! corresponding element in the TimeIdx collection, according to the given
! TypeID. TypeID can be one of the following:
!
! \begin{itemize}
! \item 1:   Constant
! \item 24:  Hourly
! \item 241: Hourly\_Grid
! \item 7:   Weekday
! \item 71:  Weekday\_Grid
! \item 12:  Monthly
! \end{itemize}
!
! !INTERFACE:
!
  SUBROUTINE tIDx_Set( HcoState, ctIDx, TypeID )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER     :: HcoState ! HEMCO state
    TYPE(TimeIdx),   POINTER     :: ctIDx   ! container TimeIDx
    INTEGER,         INTENT(IN)  :: TypeID  ! type ID
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! tIDx_Set begins here!
    !======================================================================

    ! Set pointer according to TempRes
    SELECT CASE ( TypeID )

       CASE ( 1 )
          ctIDx => HcoState%AlltIDx%CONSTANT

       CASE ( 24 )
          ctIDx => HcoState%AlltIDx%HOURLY

       CASE ( 241 )
          ctIDx => HcoState%AlltIDx%HOURLY_GRID

       CASE ( 7 )
          ctIDx => HcoState%AlltIDx%WEEKDAY

       CASE ( 12 )
          ctIDx => HcoState%AlltIDx%MONTHLY

       CASE DEFAULT
          ctIDx => HcoState%AlltIDx%CONSTANT

    END SELECT

  END SUBROUTINE tIDx_Set
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: tIDx_Cleanup
!
! !DESCRIPTION: Subroutine tIDx\_Cleanup deallocates the time slice
! index collection.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE tIDx_Cleanup( AlltIDx )
!
! !Input/output arguments:
!
     TYPE(TimeIdxCollection), POINTER :: AlltIDx
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! tIDx_Cleanup begins here!
    !======================================================================

    IF ( ASSOCIATED( AlltIDx ) ) THEN

       IF ( ASSOCIATED(AlltIDx%CONSTANT) ) THEN
          DEALLOCATE(AlltIDx%CONSTANT)
       ENDIF

       IF ( ASSOCIATED(AlltIDx%HOURLY) ) THEN
          DEALLOCATE(AlltIDx%HOURLY)
       ENDIF

       IF ( ASSOCIATED(AlltIDx%HOURLY_GRID) ) THEN
          DEALLOCATE(AlltIDx%HOURLY_GRID)
       ENDIF

       IF ( ASSOCIATED(AlltIDx%WEEKDAY) ) THEN
          DEALLOCATE(AlltIDx%WEEKDAY)
       ENDIF

       IF ( ASSOCIATED(AlltIDx%MONTHLY) ) THEN
          DEALLOCATE(AlltIDx%MONTHLY)
       ENDIF

       ! Also deallocate AlltIDx pointer
       DEALLOCATE( AlltIDx )
    ENDIF
    AlltIDx => NULL()

  END SUBROUTINE tIDx_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: tIDx_GetIndx
!
! !DESCRIPTION: Subroutine tIDx\_GetIndx calculates the current valid
! index values for the given file data time slice type and longitude
! location
!\\
! !INTERFACE:
!
  FUNCTION tIDx_GetIndx ( HcoState, Dta, I, J ) RESULT ( Indx )
!
! !USES:
!
    USE HCO_State_Mod,    ONLY : HCO_State
    USE HCO_Types_Mod,    ONLY : FileData
    USE HCO_CLOCK_MOD,    ONLY : HcoClock_Get, HcoClock_GetLocal
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER    :: HcoState  ! Hemco state
    TYPE(FileData),  POINTER    :: Dta       ! File data object
    INTEGER,         INTENT(IN) :: I         ! Longitude index of interest
    INTEGER,         INTENT(IN) :: J         ! Latitude  index of interest
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER                     :: Indx      ! Index
!
! !REVISION HISTORY:
!  02 Dec 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: HH, WD, MM, RC
    REAL(hp)           :: LonHH
    REAL(dp)           :: frac

    !======================================================================
    ! tIDx_GetIndx begins here!
    !======================================================================

    ! Default value (=> This will cause the code to crash!)
    Indx = -1

    ! ----------------------------------------------------------------
    ! Calculate time slice index for the given time slice type
    ! ----------------------------------------------------------------
    SELECT CASE ( Dta%tIDx%TypeID )

       ! Constant: there is only one time slice
       CASE ( 1 )
          Indx = 1

       ! Hourly data (local time)
       ! Indx returns the time slice representative for the LOCAL time
       ! at longitude Lon.
       CASE ( 24 )
          CALL HcoClock_GetLocal( HcoState, I, J, cH=LonHH, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          Indx = FLOOR(LonHH) + 1

       ! Hourly data (already gridded)
       ! Gridded hourly data is assumed to be already adjusted for
       ! local time effects, hence just point to the time slice
       ! of current UTC time. Add one since hour starts at 0.
       CASE ( 241 )
          CALL HcoClock_Get( HcoState%Clock, cH=HH, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          Indx = HH + 1

       ! Weekday data (local time)
       ! Indx returns the time slice representative for the LOCAL
       ! weekday at longitude Lon.
       CASE ( 7 )

          CALL HcoClock_GetLocal( HcoState, I, J, cWeekday=WD, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          Indx = WD + 1

       ! Monthly data (local time)
       ! Monthly data is always in local time.
       ! For gridded monthly data, only the current valid time slice
       ! is kept in memory (and updated whenever a new month is entered).
       CASE ( 12 )
          CALL HcoClock_GetLocal( HcoState, I, J, cMM = MM, RC=RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          Indx = MM

       ! Default: assume it's constant
       CASE DEFAULT
          Indx = 1

    END SELECT

    ! For hourly data with less than 24 time slices, i.e. time
    ! intervals of more than 1 hour, map 24 hour index onto the
    ! reduced time slice.
    ! For example, for 3-hourly data (8 vector elements), this will
    ! return index 1 for hours 0-2am, index 2 for 3-5am, etc.
    ! Dta%DeltaT denotes the time difference between between two time
    ! slices (in hours).
    IF ( Dta%DeltaT > 1 .AND. Dta%DeltaT < 24 ) THEN
       frac = DBLE(Indx) / DBLE(Dta%DeltaT)
       Indx = CEILING( frac )
    ENDIF

    ! Sanity check: index must not exceed time dimension
    IF ( Indx > Dta%nt ) THEN
       WRITE(*,*) 'Indx exceeds # of time slices! ', Indx, Dta%nt, TRIM(Dta%ncFile)
       Indx = -1
    ENDIF

  END FUNCTION tIDx_GetIndx
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tIDx_Assign
!
! !DESCRIPTION: Subroutine tIDx\_Assign assigns the time index
! pointer of the file data object of the passed list container
! to the corresponding time index element of the time index
! collection. The time slice cycle interval is determined from
! the number of time slices (length of the data array vector) and
! the time difference (deltaT) between them.
!\\
!\\
! Note: Scale factors read directly from the HEMCO configuration file
! become their delta T assigned upon reading from file (see
! subroutine Register\_Scal in hco\_config\_mod.F90).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE tIDx_Assign( HcoState, Dct, RC )
!
! !USES:
!
    USE HCO_State_Mod,    ONLY : HCO_State
    USE HCO_Types_Mod,    ONLY : DataCont
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState  ! Hemco state
    TYPE(DataCont),  POINTER       :: Dct       ! Data container
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  13 Jan 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: nx, ny, nt, ntexp, dt
    INTEGER            :: cTypeID
    CHARACTER(LEN=255) :: MSG

    !-----------------------------------
    ! tIDx_Assign begins here!
    !-----------------------------------

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, &
                   'tIDx_Assign (hco_tidx_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Check if already done
    IF ( ASSOCIATED( Dct%Dta%tIDx ) ) THEN
       CALL HCO_LEAVE( HcoState%Config%Err, RC )
       RETURN
    ENDIF

    ! Get array dimensions
    nt = Dct%Dta%nt

    ! Do the following only for defined data arrays
    IF ( nt > 0 ) THEN

       ! -------------------------------------------------------------
       ! Extract data array dimensions and delta t [h] between time
       ! slices
       ! -------------------------------------------------------------
       dt = Dct%Dta%DeltaT
       IF ( Dct%Dta%SpaceDim <= 2 ) THEN
          nx = SIZE(Dct%Dta%V2(1)%Val,1)
          ny = SIZE(Dct%Dta%V2(1)%Val,2)
       ELSE
          nx = SIZE(Dct%Dta%V3(1)%Val,1)
          ny = SIZE(Dct%Dta%V3(1)%Val,2)
       ENDIF

       ! -------------------------------------------------------------
       ! Auto detection of temporal resolution.
       ! -------------------------------------------------------------

       ! No time variation
       IF ( nt == 1 ) THEN
          cTypeID = 1

       ! -------------------------------------------------------------
       ! Daily varying data.
       ! For now, these arrays must contain seven time slices, i.e.
       ! one for every weekday (starting w/ sunday!).
       ! -------------------------------------------------------------
       ELSEIF ( nt == 7 ) THEN

          ! Sanity check if dt = 24 hours
          IF ( dt /= 24 ) THEN
             MSG = '7 time slices but delta t is not 24 hours!' // &
                  TRIM(Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF

          ! Check if data is in local time or not, set TempRes attribute
          ! accordingly. Data will only be in local time if data is
          ! spatially uniform or country-specific data. The IsLocTime
          ! flag is set in hco_config_mod.F90 (for data read from other
          ! sources than netCDF) or in hcoio_dataread_mod.F90 (for
          ! weekdaily data read from netCDF).
          IF ( .NOT. Dct%Dta%IsLocTime ) THEN
             MSG = 'Weekday data must be in local time!' // &
                  TRIM(Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN

          ELSE
             cTypeID = 7
!             cTypeID = 71
          ENDIF

       ! -------------------------------------------------------------
       ! Monthly varying data
       ! This must be a scalar. Gridded monthly data should be read
       ! slice by slice!
       ! -------------------------------------------------------------
       ELSEIF ( nt == 12 .AND. dt /= 2 ) THEN
          IF ( nx == 1 .AND. ny == 1 ) THEN
             cTypeID = 12
          ELSE
             MSG = 'Monthly data must not be gridded:' // &
                  TRIM(Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF

       ! -------------------------------------------------------------
       ! Hourly varying data.
       ! Hourly time slices can represent more than just one hour,
       ! e.g. an array can hold only 8 time slices, each representing
       ! 3-hourly intervals of the day. The only constraint is that
       ! the number of time slices have to be a factor of 24!
       ! -------------------------------------------------------------
       ELSEIF ( dt < 24 .AND. dt > 0 ) THEN

          ! Check if we can evenly split up data within a day
          IF ( MOD(24,dt) /= 0 ) THEN
             MSG = 'Cannot properly split up hourly data!' // &
                  TRIM(Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF

          ! Also check if nt is indeed of expected size
          ntexp = 24/dt
          IF ( ntexp /= nt ) THEN
             MSG = 'Wrong delta t and/or number of time slices!' // &
                  TRIM(Dct%cName)
             CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
             RETURN
          ENDIF

          ! Check if data is in local time or not, set TempRes attribute
          ! accordingly. Data will only be in local time if data is
          ! read from other sources than netCDF. The corresponding
          ! IsLocTime flag is set in hco_config_mod.F90.
          IF ( Dct%Dta%IsLocTime ) THEN
             cTypeID = 24
          ELSE
             cTypeID = 241
          ENDIF

       ! -------------------------------------------------------------
       ! Return w/ error otherwise
       ! -------------------------------------------------------------
       ELSE
          MSG = 'Invalid time slice for field ' // &
               TRIM(Dct%cName)
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC )
          RETURN
       ENDIF

       ! -------------------------------------------------------------
       ! Establish the appropriate pointer for the
       ! 4th dimension (temporal resolution) of the field array
       ! -------------------------------------------------------------
       CALL tIDx_Set( HcoState, Dct%Dta%tIDx, cTypeID )

    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE tIDx_Assign
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: tIDx_IsInRange
!
! !DESCRIPTION: Subroutine tIDx\_IsInRange returns true if the passed datetime
! is within the range of the date ranges of the data container.
!\\
! !INTERFACE:
!
  FUNCTION tIDx_IsInRange ( Lct, Yr, Mt, Dy, Hr ) RESULT ( InRange )
!
! !USES:
!
    USE HCO_TYPES_MOD,    ONLY : ListCont
!
! !INPUT PARAMETERS:
!
    TYPE(ListCont),  POINTER    :: Lct       ! File data object
    INTEGER,         INTENT(IN) :: Yr
    INTEGER,         INTENT(IN) :: Mt
    INTEGER,         INTENT(IN) :: Dy
    INTEGER,         INTENT(IN) :: Hr
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL                     :: InRange
!
! !REVISION HISTORY:
!  04 Mar 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Init
    InRange = .TRUE.

    IF ( Lct%Dct%Dta%ncYrs(1) /= Lct%Dct%Dta%ncYrs(2) ) THEN
       IF ( Yr < Lct%Dct%Dta%ncYrs(1) .OR. &
            Yr > Lct%Dct%Dta%ncYrs(2)       ) THEN
          InRange = .FALSE.
       ENDIF
    ENDIF

    IF ( Lct%Dct%Dta%ncMts(1) /= Lct%Dct%Dta%ncMts(2) ) THEN
       IF ( Mt < Lct%Dct%Dta%ncMts(1) .OR. &
            Mt > Lct%Dct%Dta%ncMts(2)       ) THEN
          InRange = .FALSE.
       ENDIF
    ENDIF

    IF ( Lct%Dct%Dta%ncDys(1) /= Lct%Dct%Dta%ncDys(2) ) THEN
       IF ( Dy < Lct%Dct%Dta%ncDys(1) .OR. &
            Dy > Lct%Dct%Dta%ncDys(2)       ) THEN
          InRange = .FALSE.
       ENDIF
    ENDIF

    IF ( Lct%Dct%Dta%ncHrs(1) /= Lct%Dct%Dta%ncHrs(2) ) THEN
       IF ( Hr < Lct%Dct%Dta%ncHrs(1) .OR. &
            Hr > Lct%Dct%Dta%ncHrs(2)       ) THEN
          InRange = .FALSE.
       ENDIF
    ENDIF

  END FUNCTION tIDx_IsInRange
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_GetPrefTimeAttr
!
! !DESCRIPTION: Subroutine HCO\_GetPrefTimeAttr returns the preferred time
! reading attributes for a given field, based upon the specs set in the HEMCO
! configuration file and the current simulation date. This routine is used
! to select the proper time slice to be read at the given time.
!\\
!\\
! The time reading attributes are set to the value that is closest (or equal)
! to the current datetime but within the range provided in the configuration
! file. If the year, month, or day is not specified, the current simulation
! date values are taken. For unspecified hours, a value of -1 is returned
! (this allows to read all hourly slices at once).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_GetPrefTimeAttr( HcoState, Lct,            &
                                  readYr,   readMt, readDy, &
                                  readHr,   readMn, RC       )
!
! !USES:
!
    USE HCO_STATE_MOD,     ONLY : HCO_State
    USE HCO_TYPES_MOD,     ONLY : ListCont
    USE HCO_TYPES_MOD,     ONLY : HCO_CFLAG_RANGE, HCO_CFLAG_EXACT
    USE HCO_CLOCK_MOD,     ONLY : HcoClock_Get
    USE HCO_TIMESHIFT_MOD, ONLY : TimeShift_Apply
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState  ! List container
    TYPE(ListCont),  POINTER       :: Lct       ! List container
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(  OUT) :: readYr    ! preferred year
    INTEGER,         INTENT(  OUT) :: readMt    ! preferred month
    INTEGER,         INTENT(  OUT) :: readDy    ! preferred day
    INTEGER,         INTENT(  OUT) :: readHr    ! preferred hour
    INTEGER,         INTENT(  OUT) :: readMn    ! preferred minute
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  13 Jan 2014 - C. Keller - Initial version
!  29 Feb 2016 - C. Keller - Added time shift option
!  03 Mar 2017 - C. Keller - Added option to deal with UTC weekdays
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER   :: cYr, cMt, cDy, cWd, cHr, cMn, sYr
    LOGICAL   :: InRange

    !-----------------------------------
    ! HCO_GetPrefTimeAttr begins here!
    !-----------------------------------

    ! Init
    RC = HCO_SUCCESS

    ! Get current time
    CALL HcoClock_Get( HcoState%Clock,                              &
                       cYYYY    = cYr, cMM = cMt, cDD = cDy,        &
                       cWEEKDAY = cWd, cH  = cHr, cM  = cMn,        &
                       sYYYY    = sYr, RC = RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Always use simulation year when specified
    IF ( Lct%Dct%Dta%UseSimYear ) cYr = sYr

    ! preferred minute is always current one
    readMn = cMn

    ! If data is in weekdays, set day to weekday: 1=Sun, ..., 7=Sat
    IF ( Lct%Dct%Dta%ncDys(1) == 1 .AND. &
         Lct%Dct%Dta%ncDys(2) == 7        ) THEN
       cDy = cWd + 1
    ENDIF

    ! -------------------------------------------------------------
    ! If CycleFlag is set to range, the preferred datetime is
    ! the current date if we are within the provided range, and
    ! invalid otherwise.
    ! -------------------------------------------------------------
    IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGE ) THEN

       ! Are we in range?
       InRange = tIDx_IsInRange ( Lct, cYr, cMt, cDy, cHr )

       IF ( InRange ) THEN
          readYr = cYr
          readMt = cMt
          readDy = cDy
          IF ( Lct%Dct%Dta%ncHrs(1) == -1 ) THEN
             readHr = -1
          ELSE
             readHr = cHr
          ENDIF
       ELSE
          readYr = -1
          readMt = -1
          readDy = -1
          readHr = -1
       ENDIF

       ! Eventually shift reference time by amount specified
       ! in HEMCO configuration file
       CALL TimeShift_Apply ( HcoState, Lct, &
                              readYr, readMt, readDy, readHr, readMn, RC )

       ! Don't need below
       RETURN
    ENDIF

    ! -------------------------------------------------------------
    ! If CycleFlag is set to exact, the preferred datetime is
    ! always the current date.
    ! -------------------------------------------------------------
    IF ( Lct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT ) THEN
       readYr = cYr
       readMt = cMt
       readDy = cDy
       IF ( Lct%Dct%Dta%ncHrs(1) == -1 ) THEN
          readHr = -1
       ELSE
          readHr = cHr
       ENDIF

       ! Eventually shift reference time by amount specified
       ! in HEMCO configuration file
       CALL TimeShift_Apply ( HcoState, Lct, &
                              readYr, readMt, readDy, readHr, readMn, RC )

       ! Don't need below
       RETURN
    ENDIF

    ! -------------------------------------------------------------
    ! If CycleFlag is anything else, select the time attributes
    ! as specified in the configuration file and closest to
    ! current simulation date.
    ! -------------------------------------------------------------

    ! Year
    IF ( Lct%Dct%Dta%ncYrs(1) == Lct%Dct%Dta%ncYrs(2) ) THEN
       readYr = Lct%Dct%Dta%ncYrs(1)
    ELSE
       readYr = MIN( MAX(cYr,Lct%Dct%Dta%ncYrs(1)), &
                         Lct%Dct%Dta%ncYrs(2)          )
    ENDIF

    ! Month
    IF ( Lct%Dct%Dta%ncMts(1) == Lct%Dct%Dta%ncMts(2) ) THEN
       readMt = Lct%Dct%Dta%ncMts(1)
    ELSE
       readMt = MIN( MAX(cMt,Lct%Dct%Dta%ncMts(1)), &
                         Lct%Dct%Dta%ncMts(2)          )
    ENDIF

    ! Day
    IF ( Lct%Dct%Dta%ncDys(1) == Lct%Dct%Dta%ncDys(2) ) THEN
       readDy = Lct%Dct%Dta%ncDys(1)
    ELSE
       readDy = MIN( MAX(cDy,Lct%Dct%Dta%ncDys(1)), &
                         Lct%Dct%Dta%ncDys(2)          )
    ENDIF

    ! Hour
    IF ( Lct%Dct%Dta%ncHrs(1) == Lct%Dct%Dta%ncHrs(2) ) THEN
       readHr = Lct%Dct%Dta%ncHrs(1)
    ELSE
       readHr = MIN( MAX(cHr,Lct%Dct%Dta%ncHrs(1)), &
                         Lct%Dct%Dta%ncHrs(2)          )
    ENDIF

    ! For weekday data, set day to 1. The seven day slices to be
    ! read will be detected based upon the current year/month.
    IF ( Lct%Dct%Dta%ncDys(1) == -10 ) readDy = 1

    ! Don't allow invalid entries for years, months or days, i.e.
    ! force readYr, readMt and readDy to be positive!
    ! readHr can be negative, in which case all hourly fields will
    ! become read!
    IF ( readYr < 0 ) readYr = cYr
    IF ( readMt < 0 ) readMt = cMt
    IF ( readDy < 0 ) readDy = cDy

    ! Eventually shift reference time by amount specified
    ! in HEMCO configuration file
    CALL TimeShift_Apply ( HcoState, Lct, &
                           readYr, readMt, readDy, readHr, readMn, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_GetPrefTimeAttr
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_ExtractTime
!
! !DESCRIPTION: Subroutine HCO\_ExtractTime extracts the time stamp (ranges)
! from the passed character string. The string is expected to be in
! format Yr/Mt/Dy/Hr. Valid values for Yr, Mt, Dy and Hr are:
!
! \begin{enumerate}
! \item Range of values, separated by - sign: e.g. 2000-2010.
! \item Single value: 2000
! \item Wildcard character (default = *). In this case, the data interval
!  is determined automatically by HEMCO based on the number of time slices
!  found in the data set.
! \item Time tokens: $YYYY, $MM, $DD, $HH. When reading the data, these values
!  will be substituted by the current simulation date.
! \item String 'WD'. Denotes that the data contains weekday data. It is
!  expected that the first slice represents Sunday. Weekday data can be used
!  in combination with annual or monthly data. In that case, there need to be
!  seven entries for every year and/or month, respectively.
! \end{enumerate}
!
! The extracted time stamp is written into the arrays ncYrs, ncMts,
! ncDys and ncHrs of the passed data structure Dta.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ExtractTime ( HcoConfig, CharStr, Dta, RC )
!
! !USES:
!
    USE CHARPAK_MOD,        ONLY : STRSPLIT
    USE HCO_TYPES_MOD,      ONLY : FileData
    USE HCO_TYPES_MOD,      ONLY : ConfigObj
    USE HCO_TYPES_MOD,      ONLY : HCO_UFLAG_ALWAYS
    USE HCO_EXTLIST_MOD,    ONLY : HCO_GetOpt
    USE HCO_TIMESHIFT_MOD,  ONLY : TimeShift_Set
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER       :: HcoConfig  ! config obj
    CHARACTER(LEN=*), INTENT(IN   ) :: CharStr
    TYPE(FileData),   POINTER       :: Dta
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  18 Sep 2013 - C. Keller - Initial version (update)
!  29 Feb 2016 - C. Keller - Added time shift option
!  03 Mar 2017 - C. Keller - Added option to deal with UTC weekdays
!  08 Aug 2018 - C. Keller - Don't set hours to -1 for local time,
!                            this is obsolete.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: I, I0, I1, N, M
    INTEGER               :: TimeVec(8)
    CHARACTER(LEN=255)    :: SUBSTR(255), DATERNG(255), MSG
    CHARACTER(LEN=255)    :: LOC = 'HCO_ExtractTime (hco_tidx_mod.F90)'

    !=================================================================
    ! HCO_ExtractTime begins here!
    !=================================================================

    ! Check for case where time flag is set to just one wildcard
    ! character. In this case, we want to update the file on every
    ! HEMCO time step, i.e. it will be added to readlist 'Always'
    ! (in hco_readlist_mod.F90).
    IF ( TRIM(CharStr) == HCO_GetOpt(HcoConfig%ExtList,'Wildcard') ) THEN
       Dta%UpdtFlag = HCO_UFLAG_ALWAYS
       Dta%ncYrs    = -999
       Dta%ncMts    = -999
       Dta%ncDys    = -999
       Dta%ncHrs    = -999

       RC = HCO_SUCCESS
       RETURN
    ENDIF

    ! Init
    TimeVec(:) = -1

    ! Extract strings to be translated into integers
    CALL STRSPLIT( CharStr, HCO_GetOpt(HcoConfig%ExtList,'Separator'), SUBSTR, N )
    IF ( N < 4 ) THEN
       MSG = 'Time stamp must have at least 4 elements: ' // TRIM(CharStr)
       CALL HCO_ERROR( HcoConfig%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Extract year, month, day, and hour range from the four
    ! substrings

    ! Do for each substring:
    DO I = 1,4

       ! Indices in TimeVec:
       I0 = (I-1) * 2 + 1
       I1 = I0 + 1

       ! For wildcard character, set lower and upper limit both to -1.
       ! In this case, the whole time slice will be read into file!
       IF ( TRIM(SUBSTR(I)) == TRIM(HCO_GetOpt(HcoConfig%ExtList,'Wildcard') ) ) THEN
          TimeVec(I0:I1) = -1

       ! Characters YYYY, MM, DD, and/or HH can be used to ensure that
       ! the current simulation time is being used.
       ELSEIF ( INDEX( TRIM(SUBSTR(I)), 'YYYY' ) > 0 .OR. &
                INDEX( TRIM(SUBSTR(I)), 'MM'   ) > 0 .OR. &
                INDEX( TRIM(SUBSTR(I)), 'DD'   ) > 0 .OR. &
                INDEX( TRIM(SUBSTR(I)), 'HH'   ) > 0       ) THEN
          TimeVec(I0:I1) = -999

       ! For weekdaily data in UTC coordinates, set day flags to 1-7.
       ELSEIF ( I==3 .AND. INDEX( TRIM(SUBSTR(I)), 'UTCWD' ) > 0 ) THEN
          TimeVec(I0) = 1
          TimeVec(I1) = 7
          Dta%IsLocTime = .FALSE.

       ! For the daily index, value 'WD' is also supported. This
       ! indicates weekdays (Sun, Mon, ..., Sat). Use a special
       ! flag here to expliclity state that these are weekday data.
       ELSEIF ( I==3 .AND. INDEX( TRIM(SUBSTR(I)), 'WD' ) > 0 ) THEN
          TimeVec(I0:I1) = -10

       ! For the hourly index, value 'LH' is also supported. This
       ! indicates local hours. Use a special flag here to expliclity
       ! state that these data are local hours.
       ELSEIF ( I==4 .AND. INDEX( TRIM(SUBSTR(I)), 'LH' ) > 0 ) THEN
          TimeVec(I0:I1) = -10
          Dta%IsLocTime = .TRUE.

       ! Otherwise, check for date range and set lower and upper bound
       ! accordingly.
       ELSE
          CALL STRSPLIT( SUBSTR(I), '-', DATERNG, M )

          ! If range is given:
          IF ( M == 2 ) THEN
             READ( DATERNG(1), * ) TimeVec(I0)
             READ( DATERNG(2), * ) TimeVec(I1)

          ! Use same value if only one value is given:
          ELSEIF ( M == 1 ) THEN
             READ( DATERNG(1), * ) TimeVec(I0)
             TimeVec(I1) = TimeVec(I0)
          ELSE
             MSG = 'Cannot extract time stamp: ' // TRIM(CharStr)
             CALL HCO_ERROR( HcoConfig%Err, MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
       ENDIF

    ENDDO !I

    ! Pass to list container
    Dta%ncYrs = TimeVec(1:2)
    Dta%ncMts = TimeVec(3:4)
    Dta%ncDys = TimeVec(5:6)
    Dta%ncHrs = TimeVec(7:8)

    ! Check for local times.
    ! Hourly and daily data that is in local time shall be completely
    ! read into memory. This will ensure that for every time zone, the
    ! correct values can be selected.
!    !IF ( Dta%IsLocTime ) THEN
!       !Dta%ncDys = -1
!       !Dta%ncHrs = -1
!    ENDIF

    ! Weekdaily data is always in local time. All seven time slices will
    ! be read into memory.
    IF ( Dta%ncDys(1) == -10 ) Dta%IsLocTime  = .TRUE.

    ! If time shift is specified, archive it in attribute 'tShift'.
    IF ( N > 4 ) THEN
       CALL TimeShift_Set( HcoConfig, Dta, SUBSTR(5), RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ExtractTime
!EOC
END MODULE HCO_TIDX_MOD
