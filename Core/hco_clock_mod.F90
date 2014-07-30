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
! values should be updated on every time step (--> HcoClock\_Set).
!\\
!\\
! The HEMCO clock object HcoClock is a private object and cannot be
! accessed directly from outside of this module. The HcoClock\_Get 
! routine should be used instead. There are also some wrapper routines
! for frequently used checks, i.e. if this is a new year, month, etc.
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
  PUBLIC :: HcoClock_Set
  PUBLIC :: HcoClock_Get
  PUBLIC :: HcoClock_GetLocal
  PUBLIC :: HcoClock_Cleanup
  PUBLIC :: HcoClock_NewYear
  PUBLIC :: HcoClock_NewMonth
  PUBLIC :: HcoClock_NewDay
  PUBLIC :: HcoClock_NewHour
  PUBLIC :: HcoClock_First
  PUBLIC :: HcoClock_GetMinResetFlag
  PUBLIC :: HcoClock_CalcDOY
  PUBLIC :: HcoClock_Increase
!
! !REMARKS:
!  The current local time implementation assumes a regular grid,
!  i.e. local time does not change with latitude
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller   - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  ! Reset flags used by diagnostics. These are used to identify
  ! the diagnostics that are at the end of their averaging interval.
  INTEGER, PARAMETER, PUBLIC  :: ResetFlagAnnually = 1
  INTEGER, PARAMETER, PUBLIC  :: ResetFlagMonthly  = 10
  INTEGER, PARAMETER, PUBLIC  :: ResetFlagDaily    = 100
  INTEGER, PARAMETER, PUBLIC  :: ResetFlagHourly   = 1000
!
! !PRIVATE TYPES:
!
  !-------------------------------------------------------------------------
  ! HCO_Clock: Derived type definition for the HEMCO clock object
  !-------------------------------------------------------------------------
  TYPE HCO_Clock

     ! Current emission time stamp (UTC)
     INTEGER            :: ThisYear     ! year 
     INTEGER            :: ThisMonth    ! month
     INTEGER            :: ThisDay      ! day
     INTEGER            :: ThisHour     ! hour
     INTEGER            :: ThisMin      ! minute
     INTEGER            :: ThisSec      ! second
     INTEGER            :: ThisDOY      ! day of year
     INTEGER            :: ThisWD       ! weekday (0=Sun,...,6=Sat)
     INTEGER            :: MonthLastDay ! Last day of month: 28,29,30,31

     ! Current local times
     INTEGER            :: nx              ! vector length == # of lons
     INTEGER,  POINTER  :: ThisLocYear(:)  ! local year 
     INTEGER,  POINTER  :: ThisLocMonth(:) ! local month 
     INTEGER,  POINTER  :: ThisLocDay(:)   ! local day
     INTEGER,  POINTER  :: ThisLocWD(:)    ! local weekday
     REAL(sp), POINTER  :: ThisLocHour(:)  ! local hour

     ! Previous emission time stamp (UTC)
     INTEGER            :: PrevYear 
     INTEGER            :: PrevMonth   
     INTEGER            :: PrevDay  
     INTEGER            :: PrevHour    
     INTEGER            :: PrevMin 
     INTEGER            :: PrevSec 
     INTEGER            :: PrevDOY 
     INTEGER            :: PrevWD

     ! total number of emission time steps 
     INTEGER            :: nSteps

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
  INTEGER                     :: CurrMinResetFlag  = ResetFlagHourly + 1

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

       ! local time vectors
       HcoClock%NX = HcoState%NX

       ALLOCATE ( HcoClock%ThisLocYear(HcoClock%NX), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( 'ThisLocYear', RC )
          RETURN
       ENDIF
       HcoClock%ThisLocYear(:) = -1

       ALLOCATE ( HcoClock%ThisLocMonth(HcoClock%NX), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( 'ThisLocMonth', RC )
          RETURN
       ENDIF
       HcoClock%ThisLocMonth(:) = -1

       ALLOCATE ( HcoClock%ThisLocDay(HcoClock%NX), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( 'ThisLocDay', RC )
          RETURN
       ENDIF
       HcoClock%ThisLocDay(:) = -1

       ALLOCATE ( HcoClock%ThisLocWD(HcoClock%NX), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( 'ThisLocWD', RC )
          RETURN
       ENDIF
       HcoClock%ThisLocWD(:) = -1

       ALLOCATE ( HcoClock%ThisLocHour(HcoClock%NX), STAT=AS )
       IF ( AS /= 0 ) THEN
          CALL HCO_ERROR ( 'ThisLocHour', RC )
          RETURN
       ENDIF
       HcoClock%ThisLocHour(:) = -1.0_sp

       HcoClock%nSteps = 0
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
  SUBROUTINE HcoClock_Set ( am_I_Root, HcoState, cYr,  cMt,  cDy, &
                            cHr,       cMin,     cSec, cDOY, RC    )
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
    REAL(sp)            :: UTC
    INTEGER             :: DOY
    CHARACTER(LEN=255)  :: MSG

    !======================================================================
    ! HcoClock_Set begins here!
    !======================================================================

    ! Assume success until otherwise
    RC = HCO_SUCCESS

    ! ----------------------------------------------------------------
    ! Don't update if this date is already in memory 
    ! ----------------------------------------------------------------
    IF ( HcoClock%ThisYear==cYr  .AND. HcoClock%ThisMonth==cMt .AND. &
         HcoClock%ThisDay ==cDy  .AND. HcoClock%ThisHour ==cHr .AND. &
         HcoClock%ThisMin ==cMin .AND. HcoClock%ThisSec  ==cSec ) THEN
       RETURN
    ENDIF

    ! ----------------------------------------------------------------
    ! Update current and previous time stamps 
    ! ----------------------------------------------------------------
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
    CALL Set_LocalTime ( HcoState, UTC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! ----------------------------------------------------------------
    ! Update diagnostics reset flag 
    ! ----------------------------------------------------------------
    CurrMinResetFlag = HcoClock_SetMinResetFlag()

    ! ----------------------------------------------------------------
    ! Verbose mode
    ! ----------------------------------------------------------------
    IF ( am_I_Root .AND. HCO_VERBOSE_CHECK() ) THEN
       WRITE(MSG,100) HcoClock%ThisYear, HcoClock%ThisMonth, &
                      HcoClock%ThisDay,  HcoClock%ThisHour,  &
                      HcoClock%ThisMin,  HcoClock%ThisSec
       CALL HCO_MSG(MSG,SEP1=' ')
       WRITE(MSG,110) HcoClock%ThisWD 
       CALL HCO_MSG(MSG,SEP2=' ')
    ENDIF

    ! ----------------------------------------------------------------
    ! Update / reset counters
    ! ----------------------------------------------------------------
    HcoClock%nSteps = HcoClock%nSteps + 1

100 FORMAT( 'Set HEMCO clock to ', i4,'-',i2.2,'-',i2.2,' ', &
            i2.2,':',i2.2,':',i2.2 )
110 FORMAT( 'The weekday is (0=Sun,...,6=Sat): ', i1.1 )

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
  SUBROUTINE HcoClock_Get ( cYYYY, cMM, cDD,  cH,           & 
                            cM,    cS,  cDOY, cWEEKDAY,     &
                            pYYYY, pMM, pDD,  pH,           &
                            pM,    pS,  pDOY, pWEEKDAY,     &
                            LMD,   nSteps,    RC             ) 
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(  OUT), OPTIONAL     :: cYYYY     ! Current year   
    INTEGER, INTENT(  OUT), OPTIONAL     :: cMM       ! Current month 
    INTEGER, INTENT(  OUT), OPTIONAL     :: cDD       ! Current day     
    INTEGER, INTENT(  OUT), OPTIONAL     :: cH        ! Current hour   
    INTEGER, INTENT(  OUT), OPTIONAL     :: cM        ! Current minute 
    INTEGER, INTENT(  OUT), OPTIONAL     :: cS        ! Current second 
    INTEGER, INTENT(  OUT), OPTIONAL     :: cDOY      ! Current day of year
    INTEGER, INTENT(  OUT), OPTIONAL     :: cWEEKDAY  ! Current weekday
    INTEGER, INTENT(  OUT), OPTIONAL     :: pYYYY     ! Previous year   
    INTEGER, INTENT(  OUT), OPTIONAL     :: pMM       ! Previous month  
    INTEGER, INTENT(  OUT), OPTIONAL     :: pDD       ! Previous day    
    INTEGER, INTENT(  OUT), OPTIONAL     :: pH        ! Previous hour   
    INTEGER, INTENT(  OUT), OPTIONAL     :: pM        ! Previous minute 
    INTEGER, INTENT(  OUT), OPTIONAL     :: pS        ! Previous second 
    INTEGER, INTENT(  OUT), OPTIONAL     :: pDOY      ! Previous day of year
    INTEGER, INTENT(  OUT), OPTIONAL     :: pWEEKDAY  ! Previous weekday
    INTEGER, INTENT(  OUT), OPTIONAL     :: LMD       ! Last day of month 
    INTEGER, INTENT(  OUT), OPTIONAL     :: nSteps    ! # of passed steps 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT)               :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! HcoClock_Get begins here!
    !======================================================================

    IF ( PRESENT(cYYYY   ) ) cYYYY    = HcoClock%ThisYear
    IF ( PRESENT(cMM     ) ) cMM      = HcoClock%ThisMonth
    IF ( PRESENT(cDD     ) ) cDD      = HcoClock%ThisDay 
    IF ( PRESENT(cH      ) ) cH       = HcoClock%ThisHour
    IF ( PRESENT(cM      ) ) cM       = HcoClock%ThisMin 
    IF ( PRESENT(cS      ) ) cS       = HcoClock%ThisSec 
    IF ( PRESENT(cDOY    ) ) cDOY     = HcoClock%ThisDOY
    IF ( PRESENT(cWEEKDAY) ) cWEEKDAY = HcoClock%ThisWD

    IF ( PRESENT(pYYYY   ) ) pYYYY    = HcoClock%PrevYear
    IF ( PRESENT(pMM     ) ) pMM      = HcoClock%PrevMonth
    IF ( PRESENT(pDD     ) ) pDD      = HcoClock%PrevDay 
    IF ( PRESENT(pH      ) ) pH       = HcoClock%PrevHour
    IF ( PRESENT(pM      ) ) pM       = HcoClock%PrevMin 
    IF ( PRESENT(pS      ) ) pS       = HcoClock%PrevSec 
    IF ( PRESENT(pDOY    ) ) pDOY     = HcoClock%PrevDOY
    IF ( PRESENT(pWEEKDAY) ) pWEEKDAY = HcoClock%PrevWD
    
    IF ( PRESENT(LMD     ) ) LMD      = HcoClock%MonthLastDay

    IF ( PRESENT(nSteps  ) ) nSteps   = HcoClock%nSteps

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
! time variables from the HEMCO clock object. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoClock_GetLocal ( IX,  cYYYY, cMM,          &
                                 cDD, cH,    CWEEKDAY, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN   )           :: IX        ! Latitude index
!
! !OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(  OUT), OPTIONAL :: cYYYY     ! Current year   
    INTEGER,  INTENT(  OUT), OPTIONAL :: cMM       ! Current month 
    INTEGER,  INTENT(  OUT), OPTIONAL :: cDD       ! Current day     
    REAL(sp), INTENT(  OUT), OPTIONAL :: cH        ! Current hour   
    INTEGER,  INTENT(  OUT), OPTIONAL :: cWEEKDAY  ! Current weekday
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(INOUT)           :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! HcoClock_GetLocal begins here!
    !======================================================================
      
    ! Check passed index
    IF ( IX > HcoClock%NX ) THEN
       CALL HCO_ERROR ( 'longitude index too big!', RC )
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
  FUNCTION HcoClock_First() RESULT ( First )
!
! !RETURN VALUE:
!
    LOGICAL :: First
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC

    First = ( HcoClock%PrevYear < 0 )      

  END FUNCTION HcoClock_First
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
  FUNCTION HcoClock_NewYear() RESULT ( NewYear )
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

    NewYear = ( HcoClock%ThisYear /= HcoClock%PrevYear )      

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
  FUNCTION HcoClock_NewMonth() RESULT ( NewMonth )
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

    NewMonth = ( HcoClock%ThisMonth /= HcoClock%PrevMonth )      

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
  FUNCTION HcoClock_NewDay() RESULT ( NewDay )
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

    NewDay = ( HcoClock%ThisDay /= HcoClock%PrevDay )      

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
  FUNCTION HcoClock_NewHour() RESULT ( NewHour )
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

    NewHour = ( HcoClock%ThisHour /= HcoClock%PrevHour )      

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
    CurrMinResetFlag  = ResetFlagHourly + 1

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
  FUNCTION HcoClock_SetMinResetFlag RESULT ( MinResetFlag ) 
!
! !RETURN VALUE:
!
    INTEGER :: MinResetFlag 
!
! !REVISION HISTORY: 
!  13 Jan 2014 - C. Keller - Initial version 
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC

    !-----------------------------------
    ! HCOCLOCK_SETMINRESETFLAG begins here! 
    !-----------------------------------

    MinResetFlag = ResetFlagHourly + 1

    ! MinResetFlag is the smallest ResetFlag number that has to
    ! be considered for the current time stamp. ResetFlag increases
    ! with reset frequency (1 = annually, 2 = monthly, etc.), so 
    ! if MinResetFlag is 2 (new month), all containers with 
    ! ResetFlag 2, 3, and 4 (monthly, daily, hourly) should be 
    ! considered.
    IF ( HcoClock_First() ) THEN
       ! MinResetFlag should be default on first HEMCO call! 
    ELSEIF ( HcoClock_NewYear() ) THEN
       MinResetFlag = ResetFlagAnnually 
    ELSEIF ( HcoClock_NewMonth() ) THEN
       MinResetFlag = ResetFlagMonthly
    ELSEIF ( HcoClock_NewDay() ) THEN
       MinResetFlag = ResetFlagDaily
    ELSEIF ( HcoClock_NewHour() ) THEN
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
  FUNCTION HcoClock_GetMinResetFlag RESULT ( MinResetFlag ) 
!
! !RETURN VALUE:
!
    INTEGER :: MinResetFlag 
!
! !REVISION HISTORY: 
!  13 Jan 2014 - C. Keller - Initial version 
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
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
! the HEMCO clock object. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_LocalTime ( HcoState, UTC, RC ) 
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_State
!
! !INPUT PARAMETERS: 
!
    REAL(sp),        INTENT(IN   ) :: UTC       ! UTC 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState  ! Hemco state
    INTEGER,         INTENT(INOUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY: 
!  13 Jan 2014 - C. Keller - Initial version 
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: I,      MtLastDay
    REAL(sp) :: DECloc
    REAL(sp) :: ThisLocHour
    INTEGER  :: ThisLocYear, ThisLocMonth
    INTEGER  :: ThisLocDay,  ThisLocWD

    !-----------------------------------
    ! SET_LOCALTIME begins here! 
    !-----------------------------------

    ! Enter
    CALL HCO_ENTER ( 'SET_LOCALTIME (HCO_CLOCK_MOD.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Loop over longitude to account for different local times.
    DO I = 1, HcoState%NX

       ! local decimal time 
       DECloc = UTC + ( HcoState%Grid%XMID( I, 1 ) / 15_sp )

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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE HcoClock_Increase ( am_I_Root, HcoState, TimeStep, RC ) 
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
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY: 
!  29 Jul 2014 - C. Keller - Initial version 
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
    UTC = ( REAL( HcoClock%ThisHour, dp )             ) + &
          ( REAL( HcoClock%ThisMin , dp ) / 60.0_dp   ) + &
          ( REAL( HcoClock%ThisSec , dp ) / 3600.0_dp )
    DAY = REAL(HcoClock%ThisDay, dp) + UTC
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
    CALL HcoClock_Set ( am_I_Root, HcoState, Yr, Mt, Dy, Hr, Mn, Sc, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HcoClock_Increase
!EOC
END MODULE HCO_CLOCK_MOD
!EOM
