!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_clock_mod.F90
!
! !DESCRIPTION: Module HCO\_Clock\_Mod contains routines and variables
! to handle the HEMCO time and calendar settings through the HEMCO clock
! object. The HEMCO clock carries information of the current UTC/local
! time as well as the UTC times from the previous time step. These
! values should be updated on every time step (--> HcoClock\_Set). It
! also contains separate variables for the current and previous emission
! datetime (year, month, day, hour, min, sec). This allows us to keep
! track of emission time steps even if the emission time steps are less
! frequent than the regular time step.
!\\
!\\
! Subroutine HcoClock\_EmissionsDone indicates that emisisons have been
! completely calculated for the current time step. Any calls to
! HcoClock\_Get will return IsEmisTime FALSE until the clock has been
! advanced to the next emission time step (via HcoClock\_Set).
!\\
!\\
! The HEMCO clock object HcoClock is a private object and cannot be
! accessed directly from outside of this module. The HcoClock\_Get
! routine should be used instead. There are also some wrapper routines
! for frequently used checks, i.e. if this is a new year, month, etc.
!\\
!\\
! Local times are calculated for 26 time zones, ranging from UTC-12 hours
! to UTC+13 hours. The time zone to be used at a given grid point is
! based on its geographical position. By default, the time zone is picked
! according to the longitude, with each time zone spanning 15 degrees.
! More detailed time zones can be provided through an external input file,
! specified in the HEMCO configuration file. The field name must be
! `TIMEZONES`, and the file must contain UTC offsets in hours. If such a
! file is provided, the time zones are determined based on these values.
! Minute offsets are ignored, e.g. UTC+9hr30min is treated as UTC+9hr. If
! the input field contains any invalid values (e.g. outside the range of
! UTC-12 - UTC+13 hours), the default algorithm is applied.
!\\
!\\
! The HEMCO clock object also controls cases where the emission dates shall
! be held constant, e.g. for simulations where emission year 2000 shall be
! used irrespective of the simulation date. Fixed simulation dates can be
! set in the settings section of the HEMCO configuration file via settings
! `Emission year`, `Emission month`, `Emission day`, and `Emission hour`.
! Only a subset of those settings can be provided, in which case all other
! time attributes will be taken from the simulation datetime.
!\\
!\\
! !INTERFACE:
!
MODULE HCO_CLOCK_MOD
!
! !USES:
!
  USE HCO_Error_Mod
  USE Julday_Mod
  USE HCO_TYPES_MOD, ONLY : HcoClock

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  ! HEMCO Clock object:
  PUBLIC :: HcoClock_Init
  PUBLIC :: HcoClock_InitTzPtr
  PUBLIC :: HcoClock_Set
  PUBLIC :: HcoClock_Get
  PUBLIC :: HcoClock_GetLocal
  PUBLIC :: HcoClock_Cleanup
  PUBLIC :: HcoClock_NewYear
  PUBLIC :: HcoClock_NewMonth
  PUBLIC :: HcoClock_NewDay
  PUBLIC :: HcoClock_NewHour
  PUBLIC :: HcoClock_New3Hour
  PUBLIC :: HcoClock_First
  PUBLIC :: HcoClock_Rewind
  PUBLIC :: HcoClock_CalcDOY
  PUBLIC :: HcoClock_Increase
  PUBLIC :: HcoClock_EmissionsDone
  PUBLIC :: HcoClock_SetLast
  PUBLIC :: Get_LastDayOfMonth
!
! !REMARKS:
!  The current local time implementation assumes a regular grid,
!  i.e. local time does not change with latitude
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller   - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  08 Oct 2014 - C. Keller   - Added mid-month day calculation
!  03 Dec 2014 - C. Keller   - Now use fixed number of time zones (24)
!  12 Jan 2015 - C. Keller   - Added emission time variables.
!  02 Feb 2015 - C. Keller   - Added option to get time zones from file
!  23 Feb 2015 - R. Yantosca - Added routine HcoClock_InitTzPtr
!  11 Jun 2015 - C. Keller   - Added simulation times and option to fix
!                              emission year, month, day, and/or hour.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
!
! !PRIVATE TYPES:
!
!
! !LOCAL VARIABLES:
!
  ! HcoClock is the variable for the HEMCO clock object
!  TYPE(HcoClock),    POINTER :: HcoClock => NULL()

  ! Midmonth days for a regular year.
  ! These can be used to obtain the mid-month day of the current month.
  INTEGER, PARAMETER   :: MidMon(13) = (/  15,  45,  74, 105,      &
                                          135, 166, 196, 227,      &
                                          258, 288, 319, 349, 380/)

  ! Number of time zones. Time zone index 1 is UTC-12. Time zone index
  ! 25 is UTC+12. Add one more to account for UTC+13.
  INTEGER, PARAMETER   :: nTimeZones = 26

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_Init
!
! !DESCRIPTION: Subroutine HcoClock\_Init initializes the HEMCO clock.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoClock_Init ( HcoState, RC )
!
! !USES:
!
    USE HCO_ARR_MOD,     ONLY : HCO_ArrInit
    USE HCO_STATE_MOD,   ONLY : HCO_State
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState  ! HEMCO state obj
    INTEGER,         INTENT(INOUT)  :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  10 Sep 2013 - C. Keller - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                  :: AS

    !======================================================================
    ! HcoClock_Init begins here!
    !======================================================================

    ! Enter
    CALL HCO_ENTER ( HcoState%Config%Err, 'HcoClock_Init (HCO_CLOCK_MOD.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Eventually allocate clock object and set all values to -1
    IF ( .NOT. ASSOCIATED ( HcoState%Clock ) ) THEN
       ALLOCATE ( HcoState%Clock )
       HcoState%Clock%PrevYear     = -1
       HcoState%Clock%PrevMonth    = -1
       HcoState%Clock%PrevDay      = -1
       HcoState%Clock%PrevHour     = -1
       HcoState%Clock%PrevMin      = -1
       HcoState%Clock%PrevSec      = -1
       HcoState%Clock%PrevDOY      = -1
       HcoState%Clock%PrevWD       = -1

       HcoState%Clock%ThisYear     = -1
       HcoState%Clock%ThisMonth    = -1
       HcoState%Clock%ThisDay      = -1
       HcoState%Clock%ThisHour     = -1
       HcoState%Clock%ThisMin      = -1
       HcoState%Clock%ThisSec      = -1
       HcoState%Clock%ThisDOY      = -1
       HcoState%Clock%ThisWD       = -1
       HcoState%Clock%MonthLastDay = -1

       HcoState%Clock%ThisEYear    = -1
       HcoState%Clock%ThisEMonth   = -1
       HcoState%Clock%ThisEDay     = -1
       HcoState%Clock%ThisEHour    = -1
       HcoState%Clock%ThisEMin     = -1
       HcoState%Clock%ThisESec     = -1

       HcoState%Clock%PrevEYear    = -1
       HcoState%Clock%PrevEMonth   = -1
       HcoState%Clock%PrevEDay     = -1
       HcoState%Clock%PrevEHour    = -1
       HcoState%Clock%PrevEMin     = -1
       HcoState%Clock%PrevESec     = -1

       HcoState%Clock%SimYear      = -1
       HcoState%Clock%SimMonth     = -1
       HcoState%Clock%SimDay       = -1
       HcoState%Clock%SimHour      = -1
       HcoState%Clock%SimMin       = -1
       HcoState%Clock%SimSec       = -1

       HcoState%Clock%isLast       = .FALSE.

       ! local time vectors
       HcoState%Clock%ntz = nTimeZones

       ALLOCATE ( HcoState%Clock%ThisLocYear(HcoState%Clock%ntz), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( HcoState%Config%Err, 'ThisLocYear', RC )
          RETURN
       ENDIF
       HcoState%Clock%ThisLocYear(:) = -1

       ALLOCATE ( HcoState%Clock%ThisLocMonth(HcoState%Clock%ntz), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'ThisLocMonth', RC )
          RETURN
       ENDIF
       HcoState%Clock%ThisLocMonth(:) = -1

       ALLOCATE ( HcoState%Clock%ThisLocDay(HcoState%Clock%ntz), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'ThisLocDay', RC )
          RETURN
       ENDIF
       HcoState%Clock%ThisLocDay(:) = -1

       ALLOCATE ( HcoState%Clock%ThisLocWD(HcoState%Clock%ntz), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'ThisLocWD', RC )
          RETURN
       ENDIF
       HcoState%Clock%ThisLocWD(:) = -1

       ALLOCATE ( HcoState%Clock%ThisLocHour(HcoState%Clock%ntz), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, 'ThisLocHour', RC )
          RETURN
       ENDIF
       HcoState%Clock%ThisLocHour(:) = -1.0_sp

       HcoState%Clock%nSteps     = 0
       HcoState%Clock%nEmisSteps = 0
       HcoState%Clock%LastEStep  = 0

       ! Initialize TIMEZONES array. Initialize as pointer (dims=0)
       CALL HCO_ArrInit( HcoState%Clock%TIMEZONES, 0, 0, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

    ENDIF

    ! Return w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

  END SUBROUTINE HcoClock_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_InitTzPtr
!
! !DESCRIPTION: Subroutine HcoClock\_InitTzPtr initializes the TIMEZONES
!  module variable.  TIMEZONES points to the timezones data (i.e. offsets
!  from UTC in hours) as read from disk.  If the timezones data file is not
!  being used, then the TIMEZONES pointer will be left unassociated.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoClock_InitTzPtr( HcoState, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
    USE HCO_EMISLIST_MOD, ONLY : HCO_GetPtr
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState  ! HcoState object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT)          :: RC         ! Success or failure?
!
! !REMARKS:
!  This routine has to be called in the HCO_Run routine, immediately after
!  the call to ReadList_Read.   The HEMCO configuration file has to be read
!  first in order to determine if we are getting our timezone information from
!  a file, or if we are computing it just based on longitude in the default
!  manner.
!
! !REVISION HISTORY:
!  23 Feb 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL :: FOUND

    ! Make sure HcoClock obj. is associated
    IF ( .NOT. ASSOCIATED(HcoState%Clock) ) THEN
       CALL HCO_WARNING( HcoState%Config%Err, &
        'CANNOT SET TIMEZONES - HEMCO CLOCK IS NOT DEFINED', &
        RC, WARNLEV=1, THISLOC='HcoClock_InitTzPtr (hco_clock_mod.F90)' )
        RETURN
    ENDIF

    ! Look for the time zone pointer
    CALL HCO_GetPtr ( HcoState, 'TIMEZONES', &
                      HcoState%Clock%TIMEZONES%Val, RC, FOUND=FOUND )

    ! Print a message
    IF ( HcoState%amIRoot ) THEN
       IF ( FOUND ) THEN
          CALL HCO_MSG( HcoState%Config%Err, &
           'TIMEZONES (i.e. OFFSETS FROM UTC) WERE READ FROM A FILE' )
       ELSE
          CALL HCO_MSG( HcoState%Config%Err, &
           'TIMEZONES (i.e. OFFSETS FROM UTC) WERE COMPUTED FROM LONGITUDE' )
       ENDIF
    ENDIF

  END SUBROUTINE HcoClock_InitTzPtr
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_Set
!
! !DESCRIPTION: Subroutine HcoClock\_Set updates the HEMCO clock. These
! routine should be called at the beginning of every emission time step!
! If the current day of year (cDoy) is not provided, it is automatically
! calculated from the current date.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoClock_Set ( HcoState, cYr, cMt, cDy, cHr, &
                            cMin, cSec, cDOY, IsEmisTime, RC    )
!
! !USES:
!
    USE HCO_TYPES_MOD,   ONLY : ConfigObj, Ext
    USE HCO_STATE_MOD,   ONLY : HCO_State
    USE HCO_EXTLIST_MOD, ONLY : GetExtOpt, CoreNr
!
! !INPUT PARAMETERS:
!
    INTEGER,         INTENT(IN   )           :: cYr       ! Current year
    INTEGER,         INTENT(IN   )           :: cMt       ! Current month
    INTEGER,         INTENT(IN   )           :: cDy       ! Current day
    INTEGER,         INTENT(IN   )           :: cHr       ! Current hour
    INTEGER,         INTENT(IN   )           :: cMin      ! Current minute
    INTEGER,         INTENT(IN   )           :: cSec      ! Current second
    INTEGER,         INTENT(IN   ), OPTIONAL :: cDoy      ! Current day of year
    LOGICAL,         INTENT(IN   ), OPTIONAL :: IsEmisTime! Is it time for emissions?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER                 :: HcoState  ! HcoState object
    INTEGER,         INTENT(INOUT)           :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller   - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  08 Jul 2014 - C. Keller   - Now calculate DOY if not provided
!  12 Jan 2015 - C. Keller   - Added IsEmisTime
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
    TYPE(HcoClock), POINTER  :: Clock
    TYPE(ConfigObj),POINTER  :: CF
    REAL(sp)                 :: UTC
    INTEGER                  :: DUM, DOY
    INTEGER                  :: UseYr, UseMt, UseDy, UseHr
    CHARACTER(LEN=255)       :: MSG
    LOGICAL                  :: FND, NewStep, EmisTime, WasEmisTime

    !======================================================================
    ! Clock_Set begins here!
    !======================================================================

    ! Assume success until otherwise
    RC = HCO_SUCCESS

    ! Get HEMCO clock object
    Clock => HcoState%Clock
    CF    => HcoState%Config

    ! ----------------------------------------------------------------
    ! On first call, check if fixed emission dates are to be used.
    ! Those can be set in the HEMCO configuration file.
    ! ----------------------------------------------------------------
    IF ( Clock%nSteps == 0 ) THEN
       CALL GetExtOpt( CF, CoreNr, 'Emission year', OptValInt=DUM, &
                       FOUND=FND, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( FND ) THEN
          Clock%FixYY = DUM
          WRITE(MSG,*) 'Emission year will be fixed to day ', Clock%FixYY
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       CALL GetExtOpt( CF, CoreNr, 'Emission month', OptValInt=DUM, &
                       FOUND=FND, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( FND ) THEN
          Clock%FixMM = DUM
          WRITE(MSG,*) 'Emission month will be fixed to day ', Clock%FixMM
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       CALL GetExtOpt( CF, CoreNr, 'Emission day', OptValInt=DUM, &
                       FOUND=FND, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( FND ) THEN
          Clock%Fixdd = DUM
          WRITE(MSG,*) 'Emission day will be fixed to day ', Clock%Fixdd
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       CALL GetExtOpt( CF, CoreNr, 'Emission hour', OptValInt=DUM, &
                       FOUND=FND, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( FND ) THEN
          Clock%Fixhh = DUM
          WRITE(MSG,*) 'Emission hour will be fixed to day ', Clock%Fixhh
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------
    ! Is this a new time step comparted to the most current one in
    ! memory?
    ! ----------------------------------------------------------------
    NewStep = .TRUE.
    IF ( Clock%SimYear==cYr  .AND. Clock%SimMonth==cMt .AND. &
         Clock%SimDay ==cDy  .AND. Clock%SimHour ==cHr .AND. &
         Clock%SimMin ==cMin .AND. Clock%SimSec  ==cSec ) THEN
       NewStep = .FALSE.
    ENDIF

    ! ----------------------------------------------------------------
    ! Update current and previous time stamps
    ! ----------------------------------------------------------------
    IF ( NewStep ) THEN
       Clock%PrevYear   = Clock%ThisYear
       Clock%PrevMonth  = Clock%ThisMonth
       Clock%PrevDay    = Clock%ThisDay
       Clock%PrevHour   = Clock%ThisHour
       Clock%PrevMin    = Clock%ThisMin
       Clock%PrevSec    = Clock%ThisSec
       Clock%PrevDOY    = Clock%ThisDOY
       Clock%PrevWD     = Clock%ThisWD

       ! Set simulation date
       Clock%SimYear    = cYr
       Clock%SimMonth   = cMt
       Clock%SimDay     = cDy
       Clock%SimHour    = cHr
       Clock%SimMin     = cMin
       Clock%SimSec     = cSec

       ! Check for fixed dates
       UseYr = cYr
       UseMt = cMt
       UseDy = cDy
       UseHr = cHr
       IF ( Clock%FixYY > 0 ) UseYr = Clock%FixYY
       IF ( Clock%FixMM > 0 ) UseMt = Clock%FixMM
       IF ( Clock%Fixdd > 0 ) UseDy = Clock%Fixdd
       IF ( Clock%Fixhh > 0 ) UseHr = Clock%Fixhh

       ! Set day of year: calculate if not specified
       IF ( PRESENT(cDOY)    .AND. Clock%FixYY<0 .AND. &
            Clock%FixMM<0 .AND. Clock%Fixdd<0        ) THEN
          DOY = cDOY
       ELSE
          DOY = HcoClock_CalcDOY( UseYr, UseMt, UseDy )
       ENDIF

       Clock%ThisYear   = UseYr
       Clock%ThisMonth  = UseMt
       Clock%ThisDay    = UseDy
       Clock%ThisHour   = UseHr
       Clock%ThisMin    = cMin
       Clock%ThisSec    = cSec
       Clock%ThisDOY    = DOY

       ! UTC decimal time
       UTC = ( REAL( Clock%ThisHour, sp )             ) + &
             ( REAL( Clock%ThisMin , sp ) / 60.0_sp   ) + &
             ( REAL( Clock%ThisSec , sp ) / 3600.0_sp )
       Clock%ThisWD = HCO_GetWeekday ( UseYr, UseMt, UseDy, UTC )

       ! ----------------------------------------------------------------
       ! Get last day of this month (only if month has changed)
       ! ----------------------------------------------------------------
       IF ( Clock%ThisMonth /= Clock%PrevMonth ) THEN
          Clock%MonthLastDay = &
               Get_LastDayOfMonth( Clock%ThisMonth, Clock%ThisYear )
       ENDIF

       ! ----------------------------------------------------------------
       ! Set local times
       ! ----------------------------------------------------------------
       CALL Set_LocalTime ( HcoState, Clock, UTC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! ----------------------------------------------------------------
       ! Update counter
       ! ----------------------------------------------------------------
       Clock%nSteps = Clock%nSteps + 1

    ENDIF !New time step

    ! ----------------------------------------------------------------
    ! Emission time steps
    ! ----------------------------------------------------------------
    EmisTime = .FALSE.
    IF ( PRESENT(IsEmisTime) ) EmisTime = IsEmisTime

    ! If this is an emission time step, force current values to be in
    ! sync with the other values.
    IF ( EmisTime ) THEN

       ! Check if previous emission time step is different
       IF ( ( Clock%ThisEYear   /= Clock%ThisYear   ) .OR. &
            ( Clock%ThisEMonth  /= Clock%ThisMonth  ) .OR. &
            ( Clock%ThisEDay    /= Clock%ThisDay    ) .OR. &
            ( Clock%ThisEHour   /= Clock%ThisHour   ) .OR. &
            ( Clock%ThisEMin    /= Clock%ThisMin    ) .OR. &
            ( Clock%ThisESec    /= Clock%ThisSec    )       ) THEN

          ! Set previous values
          Clock%PrevEYear  = Clock%ThisEYear
          Clock%PrevEMonth = Clock%ThisEMonth
          Clock%PrevEDay   = Clock%ThisEDay
          Clock%PrevEHour  = Clock%ThisEHour
          Clock%PrevEMin   = Clock%ThisEMin
          Clock%PrevESec   = Clock%ThisESec

          ! Update current values
          Clock%ThisEYear  = Clock%ThisYear
          Clock%ThisEMonth = Clock%ThisMonth
          Clock%ThisEDay   = Clock%ThisDay
          Clock%ThisEHour  = Clock%ThisHour
          Clock%ThisEMin   = Clock%ThisMin
          Clock%ThisESec   = Clock%ThisSec

          ! Increase counter
          Clock%nEmisSteps = Clock%nEmisSteps + 1

       ! Set EmisTime to false to make sure that the verbose message
       ! below won't be printed.
       ELSE
          EmisTime = .FALSE.
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------
    ! Verbose mode
    ! ----------------------------------------------------------------
    IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
       IF ( NewStep ) THEN
          WRITE(MSG,110) Clock%ThisYear, Clock%ThisMonth, &
                         Clock%ThisDay,  Clock%ThisHour,  &
                         Clock%ThisMin,  Clock%ThisSec
          CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1='=')
          WRITE(MSG,120) Clock%ThisWD
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(MSG,130) EmisTime
          CALL HCO_MSG(HcoState%Config%Err,MSG,SEP2=' ')
       ELSEIF ( EmisTime ) THEN
          WRITE(MSG,140) Clock%ThisYear, Clock%ThisMonth, &
                         Clock%ThisDay,  Clock%ThisHour,  &
                         Clock%ThisMin,  Clock%ThisSec
          CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1=' ', SEP2=' ')
       ENDIF
    ENDIF

    ! Cleanup
    Clock => NULL()
    CF    => NULL()

110 FORMAT( 'Set HEMCO clock to ', i4,'-',i2.2,'-',i2.2,' ', &
            i2.2,':',i2.2,':',i2.2 )
120 FORMAT( 'The weekday is (0=Sun,...,6=Sat): ', i1.1 )
130 FORMAT( 'Is this an emission time step   : ', L1   )
140 FORMAT( 'Declared emission time: ', i4,'-',i2.2,'-',i2.2,' ', &
            i2.2,':',i2.2,':',i2.2 )

  END SUBROUTINE HcoClock_Set
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_Get
!
! !DESCRIPTION: Subroutine HcoClock\_Get returns the selected UTC variables
! from the HEMCO clock object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoClock_Get ( Clock,                          &
                            cYYYY,   cMM, cDD,  cH,         &
                            cM,      cS,  cDOY, cWEEKDAY,   &
                            pYYYY,   pMM, pDD,  pH,         &
                            pM,      pS,  pDOY, pWEEKDAY,   &
                            sYYYY,   sMM, sDD,  sH,         &
                            sM,      sS,                    &
                            LMD,     nSteps,    cMidMon,    &
                            dslmm,   dbtwmm,    IsEmisTime, &
                            IsLast,  RC )
!
! !INPUT PARAMETERS:
!
    TYPE(HcoClock), POINTER              :: Clock   ! HEMCO clock obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(  OUT), OPTIONAL     :: cYYYY      ! Current year
    INTEGER, INTENT(  OUT), OPTIONAL     :: cMM        ! Current month
    INTEGER, INTENT(  OUT), OPTIONAL     :: cDD        ! Current day
    INTEGER, INTENT(  OUT), OPTIONAL     :: cH         ! Current hour
    INTEGER, INTENT(  OUT), OPTIONAL     :: cM         ! Current minute
    INTEGER, INTENT(  OUT), OPTIONAL     :: cS         ! Current second
    INTEGER, INTENT(  OUT), OPTIONAL     :: cDOY       ! Current day of year
    INTEGER, INTENT(  OUT), OPTIONAL     :: cWEEKDAY   ! Current weekday
    INTEGER, INTENT(  OUT), OPTIONAL     :: pYYYY      ! Previous year
    INTEGER, INTENT(  OUT), OPTIONAL     :: pMM        ! Previous month
    INTEGER, INTENT(  OUT), OPTIONAL     :: pDD        ! Previous day
    INTEGER, INTENT(  OUT), OPTIONAL     :: pH         ! Previous hour
    INTEGER, INTENT(  OUT), OPTIONAL     :: pM         ! Previous minute
    INTEGER, INTENT(  OUT), OPTIONAL     :: pS         ! Previous second
    INTEGER, INTENT(  OUT), OPTIONAL     :: sYYYY      ! Simulation year
    INTEGER, INTENT(  OUT), OPTIONAL     :: sMM        ! Simulation month
    INTEGER, INTENT(  OUT), OPTIONAL     :: sDD        ! Simulation day
    INTEGER, INTENT(  OUT), OPTIONAL     :: sH         ! Simulation hour
    INTEGER, INTENT(  OUT), OPTIONAL     :: sM         ! Simulation minute
    INTEGER, INTENT(  OUT), OPTIONAL     :: sS         ! Simulation second
    INTEGER, INTENT(  OUT), OPTIONAL     :: pDOY       ! Previous day of year
    INTEGER, INTENT(  OUT), OPTIONAL     :: pWEEKDAY   ! Previous weekday
    INTEGER, INTENT(  OUT), OPTIONAL     :: LMD        ! Last day of month
    INTEGER, INTENT(  OUT), OPTIONAL     :: nSteps     ! # of passed steps
    INTEGER, INTENT(  OUT), OPTIONAL     :: cMidMon    ! Mid-month day of curr. month
    INTEGER, INTENT(  OUT), OPTIONAL     :: dslmm      ! days since last mid-month
    INTEGER, INTENT(  OUT), OPTIONAL     :: dbtwmm     ! days between mid-month
    LOGICAL, INTENT(  OUT), OPTIONAL     :: IsEmisTime ! days between mid-month
    LOGICAL, INTENT(  OUT), OPTIONAL     :: IsLast     ! last time step?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT)               :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  08 Oct 2014 - C. Keller   - Added mid-month day arguments
!  12 Jan 2015 - C. Keller   - Added IsEmisTime
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    !======================================================================
    ! HcoClock_Get begins here!
    !======================================================================

    ! Set selected datetime variables
    IF ( PRESENT(cYYYY     ) ) cYYYY      = Clock%ThisYear
    IF ( PRESENT(cMM       ) ) cMM        = Clock%ThisMonth
    IF ( PRESENT(cDD       ) ) cDD        = Clock%ThisDay
    IF ( PRESENT(cH        ) ) cH         = Clock%ThisHour
    IF ( PRESENT(cM        ) ) cM         = Clock%ThisMin
    IF ( PRESENT(cS        ) ) cS         = Clock%ThisSec
    IF ( PRESENT(cDOY      ) ) cDOY       = Clock%ThisDOY
    IF ( PRESENT(cWEEKDAY  ) ) cWEEKDAY   = Clock%ThisWD

    IF ( PRESENT(pYYYY     ) ) pYYYY      = Clock%PrevYear
    IF ( PRESENT(pMM       ) ) pMM        = Clock%PrevMonth
    IF ( PRESENT(pDD       ) ) pDD        = Clock%PrevDay
    IF ( PRESENT(pH        ) ) pH         = Clock%PrevHour
    IF ( PRESENT(pM        ) ) pM         = Clock%PrevMin
    IF ( PRESENT(pS        ) ) pS         = Clock%PrevSec
    IF ( PRESENT(pDOY      ) ) pDOY       = Clock%PrevDOY
    IF ( PRESENT(pWEEKDAY  ) ) pWEEKDAY   = Clock%PrevWD

    IF ( PRESENT(sYYYY     ) ) sYYYY      = Clock%SimYear
    IF ( PRESENT(sMM       ) ) sMM        = Clock%SimMonth
    IF ( PRESENT(sDD       ) ) sDD        = Clock%SimDay
    IF ( PRESENT(sH        ) ) sH         = Clock%SimHour
    IF ( PRESENT(sM        ) ) sM         = Clock%SimMin
    IF ( PRESENT(sS        ) ) sS         = Clock%SimSec

    IF ( PRESENT(LMD       ) ) LMD        = Clock%MonthLastDay

    IF ( PRESENT(nSteps    ) ) nSteps     = Clock%nSteps

    ! Mid-month day related variables
    IF ( PRESENT(cMidMon   ) ) cMidMon    = MidMon(Clock%ThisMonth)

    ! Days since passing the most recent mid-month day. From modis_lai_mod.F90
    IF ( PRESENT(dslmm ) ) THEN
       IF ( Clock%ThisDOY < MidMon(1) ) THEN
          dslmm = 365 + Clock%ThisDoy - MidMon(12)
       ELSE
          dslmm = MidMon(Clock%ThisMonth+1) - MidMon(Clock%ThisMonth)
       ENDIF
    ENDIF

    ! Days between most recently passed mid-month day and next one.
    IF ( PRESENT(dbtwmm ) ) THEN

       ! If day of year is earlier than first mid-month day, we are between
       ! December and January.
       IF ( Clock%ThisDOY < MidMon(1) ) THEN
          dbtwmm = MidMon(13) - MidMon(12)

       ! If day of year is earlier than mid-month day of current month, the
       ! day difference has to be taken relative to previous month' mid-day
       ELSEIF ( Clock%ThisDOY < MidMon(Clock%ThisMonth) ) THEN
          dbtwmm = MidMon(Clock%ThisMonth) - MidMon(Clock%ThisMonth-1)

       ! If day of year is after than mid-month day of current month, the
       ! day difference has to be taken relative to current month' mid-day
       ELSE
          dbtwmm = MidMon(Clock%ThisMonth+1) - MidMon(Clock%ThisMonth)
       ENDIF
    ENDIF

    ! Is it time for emissions?
    IF ( PRESENT(IsEmisTime) ) THEN
       IsEmisTime = .FALSE.
       IF ( ( Clock%ThisEYear  == Clock%ThisYear   ) .AND. &
            ( Clock%ThisEMonth == Clock%ThisMonth  ) .AND. &
            ( Clock%ThisEDay   == Clock%ThisDay    ) .AND. &
            ( Clock%ThisEHour  == Clock%ThisHour   ) .AND. &
            ( Clock%ThisEMin   == Clock%ThisMin    ) .AND. &
            ( Clock%ThisESec   == Clock%ThisSec    ) .AND. &
            ( Clock%LastEStep  /= Clock%nEmisSteps )        ) THEN
          IsEmisTime = .TRUE.
       ENDIF
    ENDIF

    ! Last step
    IF ( PRESENT(IsLast) ) IsLast = Clock%isLast

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HcoClock_Get
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_GetLocal
!
! !DESCRIPTION: Subroutine HcoClock\_GetLocal returns the selected local
! time variables from the HEMCO clock object for the given longitude and
! latitude. At the moment, the time zone is selected purely on the given
! longitude and the passed latitude is not evaluated.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoClock_GetLocal ( HcoState,   I,     J, cYYYY, cMM, &
                                 cDD, cH,    CWEEKDAY, RC, verb     )
!
! !USES:
!
    USE HCO_STATE_MOD,    ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER          :: HcoState  ! Hemco state
    INTEGER,  INTENT(IN   )           :: I         ! Longitude index
    INTEGER,  INTENT(IN   )           :: J         ! Latitude  index
!
! !OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(  OUT), OPTIONAL :: cYYYY     ! Current year
    INTEGER,  INTENT(  OUT), OPTIONAL :: cMM       ! Current month
    INTEGER,  INTENT(  OUT), OPTIONAL :: cDD       ! Current day
    REAL(hp), INTENT(  OUT), OPTIONAL :: cH        ! Current hour
    INTEGER,  INTENT(  OUT), OPTIONAL :: cWEEKDAY  ! Current weekday
    INTEGER,  INTENT(IN   ), OPTIONAL :: verb      ! verbose
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(INOUT)           :: RC        ! Success or failure?
!
! !REMARKS:
!  Module variable TIMEZONES points to the timezone data (i.e. offsets in
!  hours from UTC) as read from disk.  The data file containing UTC offsets
!  is specified in the "NON-EMISSIONS DATA" section of the HEMCO configuraiton
!  file, under the container name "TIMEZONES".
!
!  The TIMEZONES module variable is initialized by calling HcoClock_InitTzPtr.
!  in the HEMCO run method HCO_Run (in module hco_driver_mod.F90).  The call
!  to HcoClock_InitTzPtr immediately follows the call to ReadList_Read, and
!  is only done on the very first emissions timestep.  The reason we have to
!  initialize the TIMEZONES module variable in the run method (instead of in
!  the init method) is because the HEMCO configuration file has to be read
!  before the timezones  data can be loaded into a HEMCO data container.
!
!  If we are not reading timezone data from a file, then the TIMEZONES
!  module variable will remain unassociated.
!
!  This fix was necessary in order to avoid segmentation faults when running
!  with OpenMP parallelization turned on.
!
!     -- Bob Yantosca (23 Feb 2015)
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller   - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  23 Feb 2015 - R. Yantosca - Remove call to Hco_GetPtr: this was causing
!                              errors on runs with OpenMP parallelization
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Local variables
    REAL(hp)                 :: LON
    INTEGER                  :: IX, OFFSET
    LOGICAL                  :: FOUND

    !======================================================================
    ! HcoClock_GetLocal begins here!
    !======================================================================

    ! Get longitude (degrees east)
    LON = HcoState%Grid%XMID%Val(I,J)

    ! Longitude must be between -180 and +180
    IF ( LON < -180_hp ) THEN
       DO WHILE ( LON < 180_hp )
          LON = LON + 360_hp
       ENDDO
    ENDIF
    IF ( LON > 180_hp ) THEN
       DO WHILE ( LON > 180_hp )
          LON = LON - 360_hp
       ENDDO
    ENDIF

    ! Get time zone index for the given position (longitude and latitude).
    ! Use gridded time zones if available, and default 15 degrees time zone
    ! bins otherwise.

    ! Init
    IX = -1

    ! First try to get time zone index from gridded data
    IF ( ASSOCIATED(HcoState%Clock%TIMEZONES) ) THEN
       IF ( ASSOCIATED(HcoState%Clock%TIMEZONES%Val) ) THEN

          ! Offset from UTC in hours
          OFFSET = FLOOR(HcoState%Clock%TIMEZONES%Val(I,J))

          ! Extract time zone index from offset. Index 13 is UTC=0.
          ! Valid offset is between -12 and +13
          IF ( OFFSET >= -12 .AND. OFFSET <= 13 ) THEN
             IX = 13 + OFFSET
          ENDIF

       ENDIF
    ENDIF

    ! Use default approach if (a) time zone file is not provided; (b) no valid
    ! time zone was found for this grid box.
    IF ( IX < 0 ) THEN
       ! Get time zone index for the given longitude, i.e. see into which time
       ! zone the given longitude falls.
       LON = ( LON + 180_hp ) / 15_hp
       IX  = FLOOR(LON) + 1

       ! Avoid ix=25 if longitude is exactly 180
       IF ( IX==25 ) IX = 1
    ENDIF

    ! Check time zone index
    IF ( IX > HcoState%Clock%ntz ) THEN
       CALL HCO_ERROR ( HcoState%Config%Err, 'time zone index too large!', RC )
       RETURN
    ENDIF

    ! Set defined variables
    IF ( PRESENT(cYYYY   ) ) cYYYY    = HcoState%Clock%ThisLocYear(IX)
    IF ( PRESENT(cMM     ) ) cMM      = HcoState%Clock%ThisLocMonth(IX)
    IF ( PRESENT(cDD     ) ) cDD      = HcoState%Clock%ThisLocDay(IX)
    IF ( PRESENT(cH      ) ) cH       = HcoState%Clock%ThisLocHour(IX)
    IF ( PRESENT(cWEEKDAY) ) cWEEKDAY = HcoState%Clock%ThisLocWD(IX)

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HcoClock_GetLocal
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_First
!
! !DESCRIPTION: Function HcoClock\_First returns TRUE on the first HEMCO
! time step, FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HcoClock_First( Clock, EmisTime ) RESULT ( First )
!
! !INPUT ARGUMENTS:
!
    TYPE(HcoClock),  POINTER    :: Clock
    LOGICAL,         INTENT(IN) :: EmisTime
!
! !RETURN VALUE:
!
    LOGICAL :: First
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller   - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  06 Apr 2015 - C. Keller   - Now use nEmisSteps and nSteps instead of
!                              previous years.
!EOP
!------------------------------------------------------------------------------
!BOC
    IF ( EmisTime ) THEN
       !First = ( Clock%nEmisSteps == 1 )
       First = ( Clock%nEmisSteps <= 1 )
    ELSE
       !First = ( Clock%nSteps     == 1 )
       First = ( Clock%nSteps     <= 1 )
    ENDIF

  END FUNCTION HcoClock_First
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_Rewind
!
! !DESCRIPTION: Function HcoClock\_Rewind returns TRUE if the last archived
! HEMCO time step is not in the past.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HcoClock_Rewind( Clock, EmisTime ) RESULT ( Rwnd )
!
! !INPUT ARGUMENTS:
!
    TYPE(HcoClock),  POINTER    :: Clock
    LOGICAL,         INTENT(IN) :: EmisTime
!
! !RETURN VALUE:
!
    LOGICAL :: Rwnd
!
! !REVISION HISTORY:
!  08 May 2015 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    INTEGER :: YYYYMMDD,  HHMMSS
    INTEGER :: pYYYYMMDD, pHHMMSS

    ! Init
    Rwnd = .FALSE.

    ! Get current and previous date & time
    IF ( EmisTime ) THEN
       YYYYMMDD  = Clock%ThisEYear  * 10000 + &
                   Clock%ThisEMonth * 100   + &
                   Clock%ThisEDay
       HHMMSS    = Clock%ThisEHour  * 10000 + &
                   Clock%ThisEMin   * 100   + &
                   Clock%ThisESec

       pYYYYMMDD = Clock%PrevEYear  * 10000 + &
                   Clock%PrevEMonth * 100   + &
                   Clock%PrevEDay
       pHHMMSS   = Clock%PrevEHour  * 10000 + &
                   Clock%PrevEMin   * 100   + &
                   Clock%PrevESec

    ELSE
       YYYYMMDD  = Clock%ThisYear  * 10000 + &
                   Clock%ThisMonth * 100   + &
                   Clock%ThisDay
       HHMMSS    = Clock%ThisHour  * 10000 + &
                   Clock%ThisMin   * 100   + &
                   Clock%ThisSec

       pYYYYMMDD = Clock%PrevYear  * 10000 + &
                   Clock%PrevMonth * 100   + &
                   Clock%PrevDay
       pHHMMSS   = Clock%PrevHour  * 10000 + &
                   Clock%PrevMin   * 100   + &
                   Clock%PrevSec
    ENDIF

    ! Check if current date & time is in future
    IF ( ( pHHMMSS   >  HHMMSS ) .AND. &
         ( pYYYYMMDD >= YYYYMMDD )      ) THEN
       Rwnd = .TRUE.
    ELSEIF ( pYYYYMMDD > YYYYMMDD ) THEN
       Rwnd = .TRUE.
    ENDIF

  END FUNCTION HcoClock_Rewind
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_NewYear
!
! !DESCRIPTION: Function HcoClock\_NewYear returns TRUE if this is a new
! year (compared to the previous emission time step), FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HcoClock_NewYear( Clock, EmisTime ) RESULT ( NewYear )
!
! !INPUT ARGUMENTS:
!
    TYPE(HcoClock),  POINTER    :: Clock
    LOGICAL,         INTENT(IN) :: EmisTime
!
! !RETURN VALUE:
!
    LOGICAL  :: NewYear
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( EmisTime ) THEN
       NewYear = ( Clock%ThisEYear /= Clock%PrevEYear )
    ELSE
       NewYear = ( Clock%ThisYear /= Clock%PrevYear )
    ENDIF

  END FUNCTION HcoClock_NewYear
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_NewMonth
!
! !DESCRIPTION: Function HcoClock\_NewMonth returns TRUE if this is a new
! month (compared to the previous emission time step), FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HcoClock_NewMonth( Clock, EmisTime ) RESULT ( NewMonth )
!
! !INPUT ARGUMENTS:
!
    TYPE(HcoClock),  POINTER    :: Clock
    LOGICAL,         INTENT(IN) :: EmisTime
!
! !RETURN VALUE:
!
    LOGICAL  :: NewMonth
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( EmisTime ) THEN
    NewMonth = ( Clock%ThisEMonth /= Clock%PrevEMonth )
    ELSE
       NewMonth = ( Clock%ThisMonth /= Clock%PrevMonth )
    ENDIF

  END FUNCTION HcoClock_NewMonth
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_NewDay
!
! !DESCRIPTION: Function HcoClock\_NewDay returns TRUE if this is a new
! day (compared to the previous emission time step), FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HcoClock_NewDay( Clock, EmisTime ) RESULT ( NewDay )
!
! !INPUT ARGUMENTS:
!
    TYPE(HcoClock),  POINTER    :: Clock
    LOGICAL,         INTENT(IN) :: EmisTime
!
! !RETURN VALUE:
!
    LOGICAL  :: NewDay
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( EmisTime ) THEN
       NewDay = ( Clock%ThisEDay /= Clock%PrevEDay )
    ELSE
       NewDay = ( Clock%ThisDay /= Clock%PrevDay )
    ENDIF

  END FUNCTION HcoClock_NewDay
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_NewHour
!
! !DESCRIPTION: Function HcoClock\_NewHour returns TRUE if this is a new
! hour (compared to the previous emission time step), FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HcoClock_NewHour( Clock, EmisTime ) RESULT ( NewHour )
!
! !INPUT ARGUMENTS:
!
    TYPE(HcoClock),  POINTER    :: Clock
    LOGICAL,         INTENT(IN) :: EmisTime
!
! !RETURN VALUE:
!
    LOGICAL  :: NewHour
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( EmisTime ) THEN
       NewHour = ( Clock%ThisEHour /= Clock%PrevEHour )
    ELSE
       NewHour = ( Clock%ThisHour /= Clock%PrevHour )
    ENDIF

  END FUNCTION HcoClock_NewHour
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_New3Hour
!
! !DESCRIPTION: Function HcoClock\_New3Hour returns TRUE if this is a new
! 3-hour timestep, FALSE otherwise.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HcoClock_New3Hour( Clock, EmisTime ) RESULT ( New3Hour )
!
! !INPUT ARGUMENTS:
!
    TYPE(HcoClock),  POINTER    :: Clock
    LOGICAL,         INTENT(IN) :: EmisTime
!
! !RETURN VALUE:
!
    LOGICAL  :: New3Hour
!
! !REVISION HISTORY:
!  08 Dec 2019 - M. Sulprizio- Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: HHMMSS

    ! Compute hour-minute-second variable
    HHMMSS = Clock%ThisHour * 10000 + Clock%ThisMin * 100 + Clock%ThisSec

    ! Read 3-hourly fields as hour 0, 3, 6, 9, 12, 15, 18 UTC
    New3Hour = ( MOD( HHMMSS, 030000 ) == 0 )

  END FUNCTION HcoClock_New3Hour
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_Cleanup
!
! !DESCRIPTION: Subroutine HcoClock\_Cleanup removes the given HcoHcoClock
! type.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoClock_Cleanup ( Clock )
!
! !USES:
!
    USE HCO_ARR_MOD,  ONLY : HCO_ArrCleanup
!
! !INPUT ARGUMENTS:
!
    TYPE(HcoClock), POINTER    :: Clock
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! HcoClock_Cleanup begins here!
    !======================================================================

    IF ( ASSOCIATED( Clock ) ) THEN
       ! Make sure TIMEZONES array does not point to any content any more.
       CALL HCO_ArrCleanup( Clock%TIMEZONES, DeepClean=.FALSE.)

       DEALLOCATE ( Clock )
    ENDIF
    Clock => NULL()

  END SUBROUTINE HcoClock_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: HCO_GetWeekday
!
! !DESCRIPTION: Function HCO\_GetWeekday returns the weekday for the
! given date (year, month, day).
! 0 = Sunday, 1 = Monday, ..., 6 = Saturday.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_GetWeekday( year, month, day, gmt ) RESULT ( weekday )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)   :: year
    INTEGER,  INTENT(IN)   :: month
    INTEGER,  INTENT(IN)   :: day
    REAL(sp), INTENT(IN)   :: gmt
!
! !RETURN VALUE:
!
    INTEGER               :: weekday
!
! ! NOTES: This function is largely based on the GEOS-Chem functions
! in time_mod.F.
!
! !REVISION HISTORY:
!  18 Dec 2013 - C. Keller - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(dp) :: A, B, JD, THISDAY, TMP

    !--------------------------
    ! HCO_GetWeekday begins here
    !--------------------------

    ! Day of week w/r/t the GMT date
    ! Use same algorithm as in routine SET_CURRENT_TIME
    THISDAY     = DAY + ( GMT / 24.0_dp )
    JD          = JULDAY( YEAR, MONTH, THISDAY )
    A           = ( JD + 1.5_dp ) / 7_dp
    B           = ( A - INT( A ) ) * 7_dp
    B           = INT( NINT( B* 1e5_dp + SIGN(5.0_dp,B) ) / 10_dp ) / 1e4_dp
    weekday     = INT( B )

  END FUNCTION HCO_GetWeekday
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_lastdayofmonth
!
! !DESCRIPTION: Function GET\_LASTDAYOFMONTH returns the last day of MONTH.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Get_LastDayOfMonth( Month, Year ) RESULT ( LastDay )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: Month
    INTEGER, INTENT(IN) :: Year
!
! !RETURN VALUE:
!
    INTEGER             :: LastDay
!
! !REVISION HISTORY:
!  13 Jan 2014 - C. Keller - Initial version
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC

    !-----------------------------------
    ! GET_LASTDAYOFMONTH begins here!
    !-----------------------------------

    ! Set default value (MSL: 11/20/14)
    LastDay = 31

    ! Select month
    SELECT CASE ( Month )

       ! Months with 31 days
       CASE (1,3,5,7,8,10,12)
          LastDay = 31

       ! Months with 30 days
       CASE (4,6,9,11)
          LastDay = 30

       ! February
       CASE (2)
          LastDay = 28

          ! Check for leap years:
          IF ( (MOD(Year,4  ) == 0) .AND. &
               (MOD(Year,400) /= 0)        ) THEN
             LastDay = 29
          ENDIF

    END SELECT

  END FUNCTION Get_LastDayOfMonth
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_LocalTime
!
! !DESCRIPTION: Subroutine Set\_LocalTime sets the local time vectors in
! the HEMCO clock object. Local time is calculated for each of the 24
! defined time zones.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_LocalTime ( HcoState, Clock, UTC, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,   ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_STATE), POINTER       :: HcoState
    TYPE(HcoClock),  POINTER       :: Clock  ! Clock object
    REAL(sp),        INTENT(IN   ) :: UTC       ! UTC time
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  13 Jan 2014 - C. Keller   - Initial version
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  03 Dec 2014 - C. Keller   - Now use fixed number of time zones (24)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I,      MtLastDay
    REAL(sp) :: LocDt, DECloc
    REAL(sp) :: ThisLocHour
    INTEGER  :: ThisLocYear, ThisLocMonth
    INTEGER  :: ThisLocDay,  ThisLocWD

    !-----------------------------------
    ! SET_LOCALTIME begins here!
    !-----------------------------------

    ! Enter
    CALL HCO_ENTER ( HcoState%Config%Err, 'SET_LOCALTIME (HCO_CLOCK_MOD.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Loop over all time zones to account for different local times.
    DO I = 1, Clock%ntz

       ! local time shift relative to UTC
       LocDt = -12_sp + real(I-1,sp)

       ! local decimal time
       DECloc = UTC + LocDt

       ! Extract local dates

       ! defaults
       ThisLocYear  = Clock%ThisYear
       ThisLocMonth = Clock%ThisMonth

       ! Case 1: Local time is one day behind UTC.
       IF ( DECloc < 0.0_sp ) THEN
          ThisLocHour = DECloc + 24_sp

          ! Adjust local weekday
          ThisLocWD = Clock%ThisWD - 1
          IF ( ThisLocWD < 0 ) ThisLocWD = 6

          ! Adjust local day. Also correct local
          ! month/year if needed!
          ThisLocDay = Clock%ThisDay - 1
          IF ( ThisLocDay == 0 ) THEN
             ThisLocMonth = ThisLocMonth - 1
             IF ( ThisLocMonth == 0 ) THEN
                ThisLocMonth = 12
                ThisLocYear  = Clock%ThisYear - 1
             ENDIF
             ThisLocDay = Get_LastDayOfMonth( ThisLocMonth, &
                                              ThisLocYear  )
          ENDIF

       ! Case 2: Local time is one day ahead UTC.
       ELSE IF ( DECloc >= 24.0_sp ) THEN
          ThisLocHour = DECloc - 24_sp

          ! Adjust local weekday
          ThisLocWD  = Clock%ThisWD + 1
          IF ( ThisLocWD > 6 ) ThisLocWD = 0

          ! Adjust local day. Also correct local
          ! month/year if needed!
          ThisLocDay = Clock%ThisDay + 1
          IF ( ThisLocDay > Clock%MonthLastDay ) THEN
             ThisLocMonth = ThisLocMonth + 1
             IF ( ThisLocMonth == 13 ) THEN
                ThisLocMonth = 1
                ThisLocYear  = Clock%ThisYear + 1
             ENDIF
             ThisLocDay = Get_LastDayOfMonth( ThisLocMonth, &
                                                ThisLocYear  )
          ENDIF

       ! Case 3: Local time is same day as UTC.
       ELSE
          ThisLocHour = DECloc

          ! local day same as utc day
          ThisLocWD   = Clock%ThisWD
          ThisLocDay  = Clock%ThisDay
       ENDIF

       ! Error trap: prevent local time from being 24
       ! (can occur due to rounding errors)
       IF ( ThisLocHour == 24.0_sp ) THEN
          ThisLocHour = 0.0_sp
       ENDIF

       ! Pass to Clock
       Clock%ThisLocYear(I)  = ThisLocYear
       Clock%ThisLocMonth(I) = ThisLocMonth
       Clock%ThisLocDay(I)   = ThisLocDay
       Clock%ThisLocWD(I)    = ThisLocWD
       Clock%ThisLocHour(I)  = ThisLocHour
    ENDDO !I

    ! Leave w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

  END SUBROUTINE Set_LocalTime
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_CalcDOY
!
! !DESCRIPTION: FUNCTION HcoClock\_CalcDOY calculates the day of year
! for the given year, month, and day.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HcoClock_CalcDOY( YYYY, MM, DD ) RESULT ( DOY )
!
! !INPUT ARGUMENTS:
!
    INTEGER, INTENT(IN) :: YYYY  ! Year
    INTEGER, INTENT(IN) :: MM    ! Month
    INTEGER, INTENT(IN) :: DD    ! Day
!
! !RETURN VALUE:
!
    INTEGER             :: DOY   ! Day of year
!
! !REVISION HISTORY:
!  08 Jul 2014 - C. Keller - Initial version

!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER :: TMP, N

    !-----------------------------------
    ! HcoClock_CalcDOY begins here
    !-----------------------------------

    ! Init
    DOY = 0

    ! Add total days of all month up to current month MM
    DO N = 1, MM-1
       TMP = Get_LastDayOfMonth( N, YYYY )
       DOY = DOY + TMP
    ENDDO

    ! Add all days of current month
    DOY = DOY + DD

  END FUNCTION HcoClock_CalcDOY
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_Increase
!
! !DESCRIPTION: Subroutine HcoClock\_Increase increases the HEMCO clock by the
! specified time.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoClock_Increase ( HcoState, TimeStep, EmisTime, RC )
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState  ! Hemco state
    REAL(sp),        INTENT(IN   ) :: TimeStep  ! Time step increase [s]
    LOGICAL,         INTENT(IN   ) :: EmisTime  ! Is new time step emission time?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  29 Jul 2014 - C. Keller - Initial version
!  08 Sep 2014 - C. Keller - Bug fix: now calculate UTC as fraction of day.
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(HcoClock),  POINTER  :: Clock
    INTEGER                   :: YYYYMMDD, HHMMSS
    INTEGER                   :: Yr, Mt, Dy, Hr, Mn, Sc
    REAL(dp)                  :: DAY, UTC, JD

    !-----------------------------------
    ! HcoClock_Increase begins here!
    !-----------------------------------

    ! Get pointer to HEMCO clock
    Clock => HcoState%Clock

    ! Get current date as Julian day.
    UTC = ( REAL(Clock%ThisHour,dp) / 24.0_dp    ) + &
          ( REAL(Clock%ThisMin ,dp) / 1440.0_dp  ) + &
          ( REAL(Clock%ThisSec ,dp) / 86400.0_dp )
    DAY = REAL(Clock%ThisDay,dp) + UTC
    JD  = JULDAY( Clock%ThisYear, Clock%ThisMonth, DAY )

    ! Add time step
    JD = JD + ( REAL(TimeStep,dp) / 86400.0_dp )

    ! Translate back into dates.
    CALL CALDATE( JD, YYYYMMDD, HHMMSS )
    Yr = FLOOR ( MOD( YYYYMMDD, 100000000) / 1.0e4_dp )
    Mt = FLOOR ( MOD( YYYYMMDD, 10000    ) / 1.0e2_dp )
    Dy = FLOOR ( MOD( YYYYMMDD, 100      ) / 1.0e0_dp )

    Hr = FLOOR ( MOD(   HHMMSS, 1000000  ) / 1.0e4_dp )
    Mn = FLOOR ( MOD(   HHMMSS, 10000    ) / 1.0e2_dp )
    Sc = FLOOR ( MOD(   HHMMSS, 100      ) / 1.0e0_dp )

    ! Update HEMCO clock to new values
    CALL HcoClock_Set ( HcoState, Yr, Mt, Dy, Hr, Mn, Sc, &
                        IsEmisTime=EmisTime, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Cleanup
    Clock => NULL()

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HcoClock_Increase
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_EmissionsDone
!
! !DESCRIPTION: Subroutine HcoClock\_EmissionsDone marks the current (emission)
! time step as having emissions completed. This is useful if the HEMCO core
! routines are called multiple times on the same time step, e.g. if there are
! two run phases.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoClock_EmissionsDone( Clock, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(HcoClock),  POINTER       :: Clock     ! HEMCO clock obj
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  13 Jan 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Update flag
    Clock%LastEStep = Clock%nEmisSteps

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HcoClock_EmissionsDone
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_SetLast
!
! !DESCRIPTION: Subroutine HcoClock\_SetLast sets the IsLast flag.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoClock_SetLast( Clock, IsLast, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(HcoClock),  POINTER       :: Clock     ! HEMCO clock obj
    LOGICAL,         INTENT(IN   ) :: IsLast    ! Is last time step?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  01 Nov 2016 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Update flag
    Clock%IsLast = IsLast

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HcoClock_SetLast
!EOC
END MODULE HCO_CLOCK_MOD
!EOM
