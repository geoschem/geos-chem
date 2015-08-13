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
! !INTERFACE: 
!
MODULE HCO_Clock_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE Julday_Mod

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
  PUBLIC :: HcoClock_First
  PUBLIC :: HcoClock_Rewind
  PUBLIC :: HcoClock_GetMinResetFlag
  PUBLIC :: HcoClock_CalcDOY
  PUBLIC :: HcoClock_Increase
  PUBLIC :: HcoClock_EmissionsDone
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  ! Reset flags used by diagnostics. These are used to identify
  ! the diagnostics that are at the end of their averaging interval.
  INTEGER, PARAMETER, PUBLIC  :: ResetFlagManually = -1
  INTEGER, PARAMETER, PUBLIC  :: ResetFlagEnd      =  0
  INTEGER, PARAMETER, PUBLIC  :: ResetFlagAnnually =  1
  INTEGER, PARAMETER, PUBLIC  :: ResetFlagMonthly  =  10
  INTEGER, PARAMETER, PUBLIC  :: ResetFlagDaily    =  100
  INTEGER, PARAMETER, PUBLIC  :: ResetFlagHourly   =  1000
  INTEGER, PARAMETER, PUBLIC  :: ResetFlagAlways   =  10000
!
! !PRIVATE TYPES:
!
  !-------------------------------------------------------------------------
  ! HCO_Clock: Derived type definition for the HEMCO clock object
  !-------------------------------------------------------------------------
  TYPE HCO_Clock

     ! Current time stamp (UTC)
     INTEGER            :: ThisYear        ! year 
     INTEGER            :: ThisMonth       ! month
     INTEGER            :: ThisDay         ! day
     INTEGER            :: ThisHour        ! hour
     INTEGER            :: ThisMin         ! minute
     INTEGER            :: ThisSec         ! second
     INTEGER            :: ThisDOY         ! day of year
     INTEGER            :: ThisWD          ! weekday (0=Sun,...,6=Sat)
     INTEGER            :: MonthLastDay    ! Last day of month: 28,29,30,31

     ! Current local times
     INTEGER            :: ntz             ! number of time zones 
     INTEGER,  POINTER  :: ThisLocYear(:)  ! local year 
     INTEGER,  POINTER  :: ThisLocMonth(:) ! local month 
     INTEGER,  POINTER  :: ThisLocDay(:)   ! local day
     INTEGER,  POINTER  :: ThisLocWD(:)    ! local weekday
     REAL(sp), POINTER  :: ThisLocHour(:)  ! local hour

     ! Previous time stamp (UTC)
     INTEGER            :: PrevYear 
     INTEGER            :: PrevMonth   
     INTEGER            :: PrevDay  
     INTEGER            :: PrevHour    
     INTEGER            :: PrevMin 
     INTEGER            :: PrevSec 
     INTEGER            :: PrevDOY 
     INTEGER            :: PrevWD

     ! Current emission time stamp
     INTEGER            :: ThisEYear
     INTEGER            :: ThisEMonth
     INTEGER            :: ThisEDay
     INTEGER            :: ThisEHour
     INTEGER            :: ThisEMin
     INTEGER            :: ThisESec

     ! Previous emission time stamp
     INTEGER            :: PrevEYear
     INTEGER            :: PrevEMonth
     INTEGER            :: PrevEDay
     INTEGER            :: PrevEHour
     INTEGER            :: PrevEMin
     INTEGER            :: PrevESec

     ! total number of elapsed time steps and emission time steps
     ! LastEStep denotes the last nEmisSteps values for which
     ! emissions have been calculated.
     INTEGER            :: nSteps
     INTEGER            :: nEmisSteps
     INTEGER            :: LastEStep 

  END TYPE HCO_Clock
!
! !LOCAL VARIABLES:
!
  ! HcoClock is the variable for the HEMCO clock object 
  TYPE(HCO_Clock),    POINTER :: HcoClock => NULL()

  ! Current minimum reset flag. This value will be reevaluated on 
  ! every time step. If we enter a new month, for instance, it will
  ! be set to ResetFlagMonthly, so that all diagnostics with a reset
  ! flag equal or higher than this value will be written to disk.
  INTEGER                     :: CurrMinResetFlag  = ResetFlagAlways + 1

  ! Midmonth days for a regular year.
  ! These can be used to obtain the mid-month day of the current month.
  INTEGER, PARAMETER   :: MidMon(13) = (/  15,  45,  74, 105,      &
                                          135, 166, 196, 227,      &
                                          258, 288, 319, 349, 380/)

  ! Number of time zones. Time zone index 1 is UTC-12. Time zone index
  ! 25 is UTC+12. Add one more to account for UTC+13. 
  INTEGER, PARAMETER   :: nTimeZones = 26

  ! Pointer to gridded time zones. Will only hold data if a field `TIMEZONES`
  ! is provided in the HEMCO configuration file.
  !
  ! NOTE: This pointer is initialized by a call to HcoClock_InitTzPtr
  ! from the HEMCO run routine (HCO_Run, in hco_driver_mod.F90).
  ! This is necessary to avoid segmentation faults when running with
  ! OpenMP turned on. (bmy, 2/23/15)
  REAL(sp), POINTER    :: TIMEZONES(:,:) => NULL()

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
    USE HCO_STATE_MOD, ONLY : HCO_State
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER        :: HcoState  ! HcoState object
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
    INTEGER :: AS

    !======================================================================
    ! HcoClock_Init begins here!
    !======================================================================

    ! Enter
    CALL HCO_ENTER ( 'HcoClock_Init (HCO_CLOCK_MOD.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Eventually allocate clock object and set all values to -1
    IF ( .NOT. ASSOCIATED ( HcoClock ) ) THEN
       ALLOCATE ( HcoClock )
       HcoClock%PrevYear     = -1
       HcoClock%PrevMonth    = -1
       HcoClock%PrevDay      = -1
       HcoClock%PrevHour     = -1
       HcoClock%PrevMin      = -1
       HcoClock%PrevSec      = -1
       HcoClock%PrevDOY      = -1
       HcoClock%PrevWD       = -1

       HcoClock%ThisYear     = -1
       HcoClock%ThisMonth    = -1
       HcoClock%ThisDay      = -1
       HcoClock%ThisHour     = -1
       HcoClock%ThisMin      = -1
       HcoClock%ThisSec      = -1
       HcoClock%ThisDOY      = -1
       HcoClock%ThisWD       = -1
       HcoClock%MonthLastDay = -1

       HcoClock%ThisEYear    = -1
       HcoClock%ThisEMonth   = -1
       HcoClock%ThisEDay     = -1
       HcoClock%ThisEHour    = -1
       HcoClock%ThisEMin     = -1
       HcoClock%ThisESec     = -1

       HcoClock%PrevEYear    = -1
       HcoClock%PrevEMonth   = -1
       HcoClock%PrevEDay     = -1
       HcoClock%PrevEHour    = -1
       HcoClock%PrevEMin     = -1
       HcoClock%PrevESec     = -1

       ! local time vectors
       HcoClock%ntz = nTimeZones 

       ALLOCATE ( HcoClock%ThisLocYear(HcoClock%ntz), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( 'ThisLocYear', RC )
          RETURN
       ENDIF
       HcoClock%ThisLocYear(:) = -1

       ALLOCATE ( HcoClock%ThisLocMonth(HcoClock%ntz), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( 'ThisLocMonth', RC )
          RETURN
       ENDIF
       HcoClock%ThisLocMonth(:) = -1

       ALLOCATE ( HcoClock%ThisLocDay(HcoClock%ntz), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( 'ThisLocDay', RC )
          RETURN
       ENDIF
       HcoClock%ThisLocDay(:) = -1

       ALLOCATE ( HcoClock%ThisLocWD(HcoClock%ntz), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( 'ThisLocWD', RC )
          RETURN
       ENDIF
       HcoClock%ThisLocWD(:) = -1

       ALLOCATE ( HcoClock%ThisLocHour(HcoClock%ntz), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( 'ThisLocHour', RC )
          RETURN
       ENDIF
       HcoClock%ThisLocHour(:) = -1.0_sp

       HcoClock%nSteps     = 0
       HcoClock%nEmisSteps = 0
       HcoClock%LastEStep  = 0

    ENDIF

    ! Return w/ success
    CALL HCO_LEAVE ( RC )

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
  SUBROUTINE HcoClock_InitTzPtr( am_I_Root, RC )
!
! !USES:
!
    USE HCO_EMISLIST_MOD, ONLY : HCO_GetPtr
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN   )  :: am_I_Root  ! Root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT)  :: RC         ! Success or failure?
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
!  10 Mar 2015 - C. Keller   - Packed message into am_I_Root statement
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL :: FOUND

    ! Look for the time zone pointer 
    CALL HCO_GetPtr ( .FALSE., 'TIMEZONES', TIMEZONES, RC, FOUND=FOUND )

    ! Print a message
    IF ( am_I_Root ) THEN
       IF ( FOUND ) THEN
          CALL HCO_MSG( &
           'TIMEZONES (i.e. OFFSETS FROM UTC) WERE READ FROM A FILE' )
       ELSE
          CALL HCO_MSG( &
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
  SUBROUTINE HcoClock_Set ( am_I_Root, HcoState, cYr,  cMt,   cDy, cHr, &
                            cMin,      cSec,     cDOY, IsEmisTime, RC    )
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   )           :: am_I_Root ! Root CPU?
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
    REAL(sp)            :: UTC
    INTEGER             :: DOY
    CHARACTER(LEN=255)  :: MSG
    LOGICAL             :: NewStep, EmisTime, WasEmisTime
   
    !======================================================================
    ! HcoClock_Set begins here!
    !======================================================================

    ! Assume success until otherwise
    RC = HCO_SUCCESS

    ! ----------------------------------------------------------------
    ! Is this a new time step comparted to the most current one in
    ! memory? 
    ! ----------------------------------------------------------------
    NewStep = .TRUE.
    IF ( HcoClock%ThisYear==cYr  .AND. HcoClock%ThisMonth==cMt .AND. &
         HcoClock%ThisDay ==cDy  .AND. HcoClock%ThisHour ==cHr .AND. &
         HcoClock%ThisMin ==cMin .AND. HcoClock%ThisSec  ==cSec ) THEN
       NewStep = .FALSE.
    ENDIF

    ! ----------------------------------------------------------------
    ! Update current and previous time stamps 
    ! ----------------------------------------------------------------
    IF ( NewStep ) THEN
       HcoClock%PrevYear   = HcoClock%ThisYear
       HcoClock%PrevMonth  = HcoClock%ThisMonth
       HcoClock%PrevDay    = HcoClock%ThisDay
       HcoClock%PrevHour   = HcoClock%ThisHour
       HcoClock%PrevMin    = HcoClock%ThisMin
       HcoClock%PrevSec    = HcoClock%ThisSec
       HcoClock%PrevDOY    = HcoClock%ThisDOY
       HcoClock%PrevWD     = HcoClock%ThisWD
   
       ! Set day of year: calculate if not specified
       IF ( PRESENT(cDOY) ) THEN
          DOY = cDOY
       ELSE
          DOY = HcoClock_CalcDOY( cYr, cMt, cDy )
       ENDIF
   
       HcoClock%ThisYear   = cYr 
       HcoClock%ThisMonth  = cMt 
       HcoClock%ThisDay    = cDy 
       HcoClock%ThisHour   = cHr 
       HcoClock%ThisMin    = cMin 
       HcoClock%ThisSec    = cSec 
       HcoClock%ThisDOY    = DOY 
   
       ! UTC decimal time
       UTC = ( REAL( HcoClock%ThisHour, sp )             ) + &
             ( REAL( HcoClock%ThisMin , sp ) / 60.0_sp   ) + &
             ( REAL( HcoClock%ThisSec , sp ) / 3600.0_sp )
       HcoClock%ThisWD = HCO_GetWeekday ( cYr, cMt, cDy, UTC ) 
   
       ! ----------------------------------------------------------------
       ! Get last day of this month (only if month has changed)
       ! ----------------------------------------------------------------
       IF ( HcoClock%ThisMonth /= HcoClock%PrevMonth ) THEN 
          HcoClock%MonthLastDay = &
               Get_LastDayOfMonth( HcoClock%ThisMonth, HcoClock%ThisYear )
       ENDIF
   
       ! ----------------------------------------------------------------
       ! Set local times 
       ! ----------------------------------------------------------------
       CALL Set_LocalTime ( UTC, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
   
       ! ----------------------------------------------------------------
       ! Update counter
       ! ----------------------------------------------------------------
       HcoClock%nSteps = HcoClock%nSteps + 1

       ! ----------------------------------------------------------------
       ! Update diagnostics reset flag
       ! Needs to be done after updating the counter 
       ! ----------------------------------------------------------------
       CurrMinResetFlag = HcoClock_SetMinResetFlag()

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
       IF ( ( HcoClock%ThisEYear   /= HcoClock%ThisYear   ) .OR. &
            ( HcoClock%ThisEMonth  /= HcoClock%ThisMonth  ) .OR. &
            ( HcoClock%ThisEDay    /= HcoClock%ThisDay    ) .OR. &
            ( HcoClock%ThisEHour   /= HcoClock%ThisHour   ) .OR. &
            ( HcoClock%ThisEMin    /= HcoClock%ThisMin    ) .OR. &
            ( HcoClock%ThisESec    /= HcoClock%ThisSec    )       ) THEN 

          ! Set previous values
          HcoClock%PrevEYear  = HcoClock%ThisEYear
          HcoClock%PrevEMonth = HcoClock%ThisEMonth
          HcoClock%PrevEDay   = HcoClock%ThisEDay
          HcoClock%PrevEHour  = HcoClock%ThisEHour
          HcoClock%PrevEMin   = HcoClock%ThisEMin
          HcoClock%PrevESec   = HcoClock%ThisESec

          ! Update current values
          HcoClock%ThisEYear  = HcoClock%ThisYear
          HcoClock%ThisEMonth = HcoClock%ThisMonth
          HcoClock%ThisEDay   = HcoClock%ThisDay
          HcoClock%ThisEHour  = HcoClock%ThisHour
          HcoClock%ThisEMin   = HcoClock%ThisMin
          HcoClock%ThisESec   = HcoClock%ThisSec

          ! Increase counter
          HcoClock%nEmisSteps = HcoClock%nEmisSteps + 1

       ! Set EmisTime to false to make sure that the verbose message
       ! below won't be printed.
       ELSE
          EmisTime = .FALSE.
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------
    ! Verbose mode
    ! ----------------------------------------------------------------
    IF ( HCO_IsVerb( 1 ) ) THEN
       IF ( NewStep ) THEN
          WRITE(MSG,110) HcoClock%ThisYear, HcoClock%ThisMonth, &
                         HcoClock%ThisDay,  HcoClock%ThisHour,  &
                         HcoClock%ThisMin,  HcoClock%ThisSec
          CALL HCO_MSG(MSG,SEP1=' ')
          WRITE(MSG,120) HcoClock%ThisWD 
          CALL HCO_MSG(MSG)
          WRITE(MSG,130) EmisTime
          CALL HCO_MSG(MSG,SEP2=' ')
       ELSEIF ( EmisTime ) THEN
          WRITE(MSG,140) HcoClock%ThisYear, HcoClock%ThisMonth, &
                         HcoClock%ThisDay,  HcoClock%ThisHour,  &
                         HcoClock%ThisMin,  HcoClock%ThisSec
          CALL HCO_MSG(MSG,SEP1=' ', SEP2=' ')
       ENDIF
    ENDIF

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
  SUBROUTINE HcoClock_Get ( cYYYY,   cMM, cDD,  cH,       & 
                            cM,      cS,  cDOY, cWEEKDAY, &
                            pYYYY,   pMM, pDD,  pH,       &
                            pM,      pS,  pDOY, pWEEKDAY, &
                            LMD,     nSteps,    cMidMon,  &
                            dslmm,   dbtwmm,    IsEmisTime, RC ) 
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
    INTEGER, INTENT(  OUT), OPTIONAL     :: pDOY       ! Previous day of year
    INTEGER, INTENT(  OUT), OPTIONAL     :: pWEEKDAY   ! Previous weekday
    INTEGER, INTENT(  OUT), OPTIONAL     :: LMD        ! Last day of month 
    INTEGER, INTENT(  OUT), OPTIONAL     :: nSteps     ! # of passed steps 
    INTEGER, INTENT(  OUT), OPTIONAL     :: cMidMon    ! Mid-month day of curr. month
    INTEGER, INTENT(  OUT), OPTIONAL     :: dslmm      ! days since last mid-month
    INTEGER, INTENT(  OUT), OPTIONAL     :: dbtwmm     ! days between mid-month
    LOGICAL, INTENT(  OUT), OPTIONAL     :: IsEmisTime ! days between mid-month
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

    !======================================================================
    ! HcoClock_Get begins here!
    !======================================================================

    IF ( PRESENT(cYYYY     ) ) cYYYY      = HcoClock%ThisYear
    IF ( PRESENT(cMM       ) ) cMM        = HcoClock%ThisMonth
    IF ( PRESENT(cDD       ) ) cDD        = HcoClock%ThisDay 
    IF ( PRESENT(cH        ) ) cH         = HcoClock%ThisHour
    IF ( PRESENT(cM        ) ) cM         = HcoClock%ThisMin 
    IF ( PRESENT(cS        ) ) cS         = HcoClock%ThisSec 
    IF ( PRESENT(cDOY      ) ) cDOY       = HcoClock%ThisDOY
    IF ( PRESENT(cWEEKDAY  ) ) cWEEKDAY   = HcoClock%ThisWD

    IF ( PRESENT(pYYYY     ) ) pYYYY      = HcoClock%PrevYear
    IF ( PRESENT(pMM       ) ) pMM        = HcoClock%PrevMonth
    IF ( PRESENT(pDD       ) ) pDD        = HcoClock%PrevDay 
    IF ( PRESENT(pH        ) ) pH         = HcoClock%PrevHour
    IF ( PRESENT(pM        ) ) pM         = HcoClock%PrevMin 
    IF ( PRESENT(pS        ) ) pS         = HcoClock%PrevSec 
    IF ( PRESENT(pDOY      ) ) pDOY       = HcoClock%PrevDOY
    IF ( PRESENT(pWEEKDAY  ) ) pWEEKDAY   = HcoClock%PrevWD
    
    IF ( PRESENT(LMD       ) ) LMD        = HcoClock%MonthLastDay

    IF ( PRESENT(nSteps    ) ) nSteps     = HcoClock%nSteps

    ! Mid-month day related variables
    IF ( PRESENT(cMidMon   ) ) cMidMon    = MidMon(HcoClock%ThisMonth)

    ! Days since passing the most recent mid-month day. From modis_lai_mod.F90
    IF ( PRESENT(dslmm ) ) THEN
       IF ( HcoClock%ThisDOY < MidMon(1) ) THEN
          dslmm = 365 + HcoClock%ThisDoy - MidMon(12)
       ELSE
          dslmm = MidMon(HcoClock%ThisMonth+1) - MidMon(HcoClock%ThisMonth)
       ENDIF
    ENDIF

    ! Days between most recently passed mid-month day and next one.
    IF ( PRESENT(dbtwmm ) ) THEN

       ! If day of year is earlier than first mid-month day, we are between
       ! December and January.
       IF ( HcoClock%ThisDOY < MidMon(1) ) THEN
          dbtwmm = MidMon(13) - MidMon(12)

       ! If day of year is earlier than mid-month day of current month, the
       ! day difference has to be taken relative to previous month' mid-day
       ELSEIF ( HcoClock%ThisDOY < MidMon(HcoClock%ThisMonth) ) THEN
          dbtwmm = MidMon(HcoClock%ThisMonth) - MidMon(HcoClock%ThisMonth-1)

       ! If day of year is after than mid-month day of current month, the
       ! day difference has to be taken relative to current month' mid-day
       ELSE
          dbtwmm = MidMon(HcoClock%ThisMonth+1) - MidMon(HcoClock%ThisMonth)
       ENDIF
    ENDIF

    ! Is it time for emissions?
    IF ( PRESENT(IsEmisTime) ) THEN
       IsEmisTime = .FALSE.
       IF ( ( HcoClock%ThisEYear  == HcoClock%ThisYear   ) .AND. & 
            ( HcoClock%ThisEMonth == HcoClock%ThisMonth  ) .AND. & 
            ( HcoClock%ThisEDay   == HcoClock%ThisDay    ) .AND. & 
            ( HcoClock%ThisEHour  == HcoClock%ThisHour   ) .AND. & 
            ( HcoClock%ThisEMin   == HcoClock%ThisMin    ) .AND. & 
            ( HcoClock%ThisESec   == HcoClock%ThisSec    ) .AND. &
            ( HcoClock%LastEStep  /= HcoClock%nEmisSteps )        ) THEN
          IsEmisTime = .TRUE.
       ENDIF 
    ENDIF

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
                                 cDD, cH,    CWEEKDAY, RC           )
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
    REAL(hp)      :: LON
    INTEGER       :: IX, OFFSET
    LOGICAL       :: FOUND
    LOGICAL, SAVE :: FIRST = .TRUE.

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
    IF ( ASSOCIATED(TIMEZONES) ) THEN

       ! Offset from UTC in hours
       OFFSET = FLOOR(TIMEZONES(I,J))

       ! Extract time zone index from offset. Index 13 is UTC=0.
       ! Valid offset is between -12 and +13
       IF ( OFFSET >= -12 .AND. OFFSET <= 13 ) THEN
          IX = 13 + OFFSET 
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
    IF ( IX > HcoClock%ntz ) THEN
       CALL HCO_ERROR ( 'time zone index too large!', RC )
       RETURN
    ENDIF

    ! Set defined variables
    IF ( PRESENT(cYYYY   ) ) cYYYY    = HcoClock%ThisLocYear(IX)
    IF ( PRESENT(cMM     ) ) cMM      = HcoClock%ThisLocMonth(IX)
    IF ( PRESENT(cDD     ) ) cDD      = HcoClock%ThisLocDay(IX)
    IF ( PRESENT(cH      ) ) cH       = HcoClock%ThisLocHour(IX)
    IF ( PRESENT(cWEEKDAY) ) cWEEKDAY = HcoClock%ThisLocWD(IX)

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
  FUNCTION HcoClock_First( EmisTime ) RESULT ( First )
!
! !INPUT ARGUMENTS:
!
    LOGICAL, INTENT(IN) :: EmisTime
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
       First = ( HcoClock%nEmisSteps == 1 )      
    ELSE
       First = ( HcoClock%nSteps     == 1 )
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
  FUNCTION HcoClock_Rewind( EmisTime ) RESULT ( Rwnd )
!
! !INPUT ARGUMENTS:
!
    LOGICAL, INTENT(IN) :: EmisTime
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
       YYYYMMDD  = HcoClock%ThisEYear  * 10000 + &
                   HcoClock%ThisEMonth * 100   + &
                   HcoClock%ThisEDay
       HHMMSS    = HcoClock%ThisEHour  * 10000 + &
                   HcoClock%ThisEMin   * 100   + &
                   HcoClock%ThisESec

       pYYYYMMDD = HcoClock%PrevEYear  * 10000 + &
                   HcoClock%PrevEMonth * 100   + &
                   HcoClock%PrevEDay
       pHHMMSS   = HcoClock%PrevEHour  * 10000 + &
                   HcoClock%PrevEMin   * 100   + &
                   HcoClock%PrevESec

    ELSE
       YYYYMMDD  = HcoClock%ThisYear  * 10000 + &
                   HcoClock%ThisMonth * 100   + &
                   HcoClock%ThisDay
       HHMMSS    = HcoClock%ThisHour  * 10000 + &
                   HcoClock%ThisMin   * 100   + &
                   HcoClock%ThisSec

       pYYYYMMDD = HcoClock%PrevYear  * 10000 + &
                   HcoClock%PrevMonth * 100   + &
                   HcoClock%PrevDay
       pHHMMSS   = HcoClock%PrevHour  * 10000 + &
                   HcoClock%PrevMin   * 100   + &
                   HcoClock%PrevSec
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
  FUNCTION HcoClock_NewYear( EmisTime ) RESULT ( NewYear )
!
! !INPUT ARGUMENTS:
!
    LOGICAL, INTENT(IN) :: EmisTime
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
       NewYear = ( HcoClock%ThisEYear /= HcoClock%PrevEYear )      
    ELSE
       NewYear = ( HcoClock%ThisYear /= HcoClock%PrevYear )      
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
  FUNCTION HcoClock_NewMonth( EmisTime ) RESULT ( NewMonth )
!
! !INPUT ARGUMENTS:
!
    LOGICAL, INTENT(IN) :: EmisTime
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
    NewMonth = ( HcoClock%ThisEMonth /= HcoClock%PrevEMonth )      
    ELSE
       NewMonth = ( HcoClock%ThisMonth /= HcoClock%PrevMonth )      
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
  FUNCTION HcoClock_NewDay( EmisTime ) RESULT ( NewDay )
!
! !INPUT ARGUMENTS:
!
    LOGICAL, INTENT(IN) :: EmisTime
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
       NewDay = ( HcoClock%ThisEDay /= HcoClock%PrevEDay )      
    ELSE
       NewDay = ( HcoClock%ThisDay /= HcoClock%PrevDay )      
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
  FUNCTION HcoClock_NewHour( EmisTime ) RESULT ( NewHour )
!
! !INPUT ARGUMENTS:
!
    LOGICAL, INTENT(IN) :: EmisTime
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
       NewHour = ( HcoClock%ThisEHour /= HcoClock%PrevEHour )      
    ELSE
       NewHour = ( HcoClock%ThisHour /= HcoClock%PrevHour )      
    ENDIF

  END FUNCTION HcoClock_NewHour
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
  SUBROUTINE HcoClock_Cleanup
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

    IF ( ASSOCIATED( HcoClock ) ) THEN
       DEALLOCATE ( HcoClock )
    ENDIF
    HcoClock => NULL()

    ! Reset current minimum reset flag to default (initial) value
    CurrMinResetFlag  = ResetFlagAlways + 1

    ! Make sure TIMEZONES array does not point to any content any more.
    TIMEZONES => NULL()

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
! !IROUTINE: HcoClock_SetMinResetFlag
!
! !DESCRIPTION: Function HcoClock\_SetMinResetFlag sets the minimum ResetFlag
! for the current HEMCO time, as used by the HEMCO diagnostics. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HcoClock_SetMinResetFlag() RESULT ( MinResetFlag ) 
!
! !RETURN VALUE:
!
    INTEGER :: MinResetFlag 
!
! !REVISION HISTORY: 
!  13 Jan 2014 - C. Keller   - Initial version 
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  12 May 2015 - R. Yantosca - Bug fix: PGI expects routine name to end w/ ()
!EOP
!------------------------------------------------------------------------------
!BOC

    !-----------------------------------
    ! HCOCLOCK_SETMINRESETFLAG begins here! 
    !-----------------------------------

    MinResetFlag = ResetFlagAlways

    ! MinResetFlag is the smallest ResetFlag number that has to
    ! be considered for the current time stamp. ResetFlag increases
    ! with reset frequency (1 = annually, 2 = monthly, etc.), so 
    ! if MinResetFlag is 2 (new month), all containers with 
    ! ResetFlag 2, 3, and 4 (monthly, daily, hourly) should be 
    ! considered.
    IF ( HcoClock_First( .FALSE. ) ) THEN
       ! MinResetFlag should be default on first HEMCO call! 
    ELSEIF ( HcoClock_NewYear( .FALSE. ) ) THEN
       MinResetFlag = ResetFlagAnnually 
    ELSEIF ( HcoClock_NewMonth( .FALSE. ) ) THEN
       MinResetFlag = ResetFlagMonthly
    ELSEIF ( HcoClock_NewDay( .FALSE. ) ) THEN
       MinResetFlag = ResetFlagDaily
    ELSEIF ( HcoClock_NewHour( .FALSE. ) ) THEN
       MinResetFlag = ResetFlagHourly
    ENDIF

  END FUNCTION HcoClock_SetMinResetFlag
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoClock_GetMinResetFlag
!
! !DESCRIPTION: Function HcoClock\_GetMinResetFlag returns the minimum 
! ResetFlag for the current HEMCO time, as used by the HEMCO diagnostics. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HcoClock_GetMinResetFlag() RESULT ( MinResetFlag ) 
!
! !RETURN VALUE:
!
    INTEGER :: MinResetFlag 
!
! !REVISION HISTORY: 
!  13 Jan 2014 - C. Keller   - Initial version 
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  12 May 2015 - R. Yantosca - Bug fix: PGI expects routine name to end w/ ()

!EOP
!------------------------------------------------------------------------------
!BOC

    !-----------------------------------
    ! HCOCLOCK_GETMINRESETFLAG begins here! 
    !-----------------------------------

    MinResetFlag = CurrMinResetFlag 

  END FUNCTION HcoClock_GetMinResetFlag
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
  SUBROUTINE Set_LocalTime ( UTC, RC ) 
!
! !USES:
!
!
! !INPUT PARAMETERS: 
!
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
    CALL HCO_ENTER ( 'SET_LOCALTIME (HCO_CLOCK_MOD.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Loop over all time zones to account for different local times.
    DO I = 1, HcoClock%ntz

       ! local time shift relative to UTC 
       LocDt = -12_sp + real(I-1,sp)

       ! local decimal time 
       DECloc = UTC + LocDt

       ! Extract local dates
         
       ! defaults
       ThisLocYear  = HcoClock%ThisYear
       ThisLocMonth = HcoClock%ThisMonth

       ! Case 1: Local time is one day behind UTC.  
       IF ( DECloc < 0.0_sp ) THEN
          ThisLocHour = DECloc + 24_sp

          ! Adjust local weekday
          ThisLocWD = HcoClock%ThisWD - 1
          IF ( ThisLocWD < 0 ) ThisLocWD = 6

          ! Adjust local day. Also correct local
          ! month/year if needed!
          ThisLocDay = HcoClock%ThisDay - 1
          IF ( ThisLocDay == 0 ) THEN
             ThisLocMonth = ThisLocMonth - 1
             IF ( ThisLocMonth == 0 ) THEN
                ThisLocMonth = 12
                ThisLocYear  = HcoClock%ThisYear - 1
             ENDIF
             ThisLocDay = Get_LastDayOfMonth( ThisLocMonth, &
                                              ThisLocYear  ) 
          ENDIF

       ! Case 2: Local time is one day ahead UTC.
       ELSE IF ( DECloc >= 24.0_sp ) THEN
          ThisLocHour = DECloc - 24_sp

          ! Adjust local weekday 
          ThisLocWD  = HcoClock%ThisWD + 1
          IF ( ThisLocWD > 6 ) ThisLocWD = 0

          ! Adjust local day. Also correct local
          ! month/year if needed!
          ThisLocDay = HcoClock%ThisDay + 1
          IF ( ThisLocDay > HcoClock%MonthLastDay ) THEN
             ThisLocMonth = ThisLocMonth + 1
             IF ( ThisLocMonth == 13 ) THEN
                ThisLocMonth = 1
                ThisLocYear  = HcoClock%ThisYear + 1
             ENDIF
             ThisLocDay = Get_LastDayOfMonth( ThisLocMonth, &
                                                ThisLocYear  ) 
          ENDIF

       ! Case 3: Local time is same day as UTC.
       ELSE
          ThisLocHour = DECloc

          ! local day same as utc day
          ThisLocWD   = HcoClock%ThisWD
          ThisLocDay  = HcoClock%ThisDay
       ENDIF
 
       ! Error trap: prevent local time from being 24
       ! (can occur due to rounding errors)
       IF ( ThisLocHour == 24.0_sp ) THEN
          ThisLocHour = 0.0_sp
       ENDIF
 
       ! Pass to HcoClock
       HcoClock%ThisLocYear(I)  = ThisLocYear
       HcoClock%ThisLocMonth(I) = ThisLocMonth
       HcoClock%ThisLocDay(I)   = ThisLocDay
       HcoClock%ThisLocWD(I)    = ThisLocWD
       HcoClock%ThisLocHour(I)  = ThisLocHour
    ENDDO !I

    ! Leave w/ success
    CALL HCO_LEAVE ( RC )

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
  SUBROUTINE HcoClock_Increase ( am_I_Root, HcoState, TimeStep, EmisTime, RC ) 
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_State
!
! !INPUT PARAMETERS: 
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root ! Root CPU 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: YYYYMMDD, HHMMSS
    INTEGER             :: Yr, Mt, Dy, Hr, Mn, Sc
    REAL(dp)            :: DAY, UTC, JD

    !-----------------------------------
    ! HcoClock_Increase begins here! 
    !-----------------------------------

    ! Get current date as Julian day.
    UTC = ( REAL(HcoClock%ThisHour,dp) / 24.0_dp    ) + &
          ( REAL(HcoClock%ThisMin ,dp) / 1440.0_dp  ) + &
          ( REAL(HcoClock%ThisSec ,dp) / 86400.0_dp )
    DAY = REAL(HcoClock%ThisDay,dp) + UTC
    JD  = JULDAY( HcoClock%ThisYear, HcoClock%ThisMonth, DAY )

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
    CALL HcoClock_Set ( am_I_Root, HcoState, Yr, Mt, Dy, Hr, Mn, Sc, &
                        IsEmisTime=EmisTime, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

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
  SUBROUTINE HcoClock_EmissionsDone( am_I_Root, RC ) 
!
! !INPUT PARAMETERS: 
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root ! Root CPU 
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
    HcoClock%LastEStep = HcoClock%nEmisSteps

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HcoClock_EmissionsDone
!EOC
END MODULE HCO_CLOCK_MOD
!EOM
