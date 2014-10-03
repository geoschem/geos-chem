!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_time_mod 
!
! !DESCRIPTION: Module HCO\_TIME\_MOD contains all routines, derived
! types and variables to handle the HEMCO time and calendar settings.
! \\
!
! Module Routines:
! ============================================================================
!
! - PUBLIC:
!
! - PRIVATE:
!
! \\
! !INTERFACE: 
!
      MODULE HCO_TIME_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD

      IMPLICIT NONE

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      ! HEMCO Clock object:
      PUBLIC :: INIT_CLOCK
      PUBLIC :: SET_CLOCK
      PUBLIC :: CLEANUP_CLOCK

      ! array indices 
      PUBLIC :: Init_tSlc
      PUBLIC :: Set_tSlc
      PUBLIC :: Update_tSlc
      PUBLIC :: Cleanup_tSlc
!
! !PUBLIC TYPES
!
      PUBLIC :: HcoClock
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!  22 Aug 2013 - C. Keller - Renamed from ng_emis_time_mod to
!                            hco_time_mod
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE MEMBER FUNCTIONS:
!
      ! HEMCO time clock
      TYPE HcoClock

         ! Current emission time stamp (UTC)
         INTEGER           :: ThisYear     ! year 
         INTEGER           :: ThisMonth    ! month
         INTEGER           :: ThisDay      ! day
         INTEGER           :: ThisHour     ! hour
         INTEGER           :: ThisMin      ! minute
         INTEGER           :: ThisSec      ! second
         INTEGER           :: ThisDOY      ! day of year
         INTEGER           :: ThisWD       ! weekday (0=Sun,...,6=Sat)
         INTEGER           :: MonthLastDay ! Last day of month: 28,29,30,31

         ! Previous emission time stamp (UTC)
         INTEGER           :: PrevYear 
         INTEGER           :: PrevMonth   
         INTEGER           :: PrevDay  
         INTEGER           :: PrevHour    
         INTEGER           :: PrevMin 
         INTEGER           :: PrevSec 
         INTEGER           :: PrevDOY 
         INTEGER           :: PrevWD
      ENDTYPE HcoClock

      ! tSlc contains the pointers to the current time slices 
      TYPE TSLC 
         INTEGER, POINTER  :: CONSTANT   (:)
         INTEGER, POINTER  :: HOURLY     (:)
         INTEGER, POINTER  :: HOURLY_GRID(:)
         INTEGER, POINTER  :: WEEKDY     (:)
         INTEGER, POINTER  :: WEEKDY_GRID(:)
         INTEGER, POINTER  :: MONTHLY    (:)
      END TYPE TSLC

      ! Types
      TYPE(tSlc),   POINTER       :: Active_tSlc => NULL()
!
      CONTAINS
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Init_tSlc
!
! !DESCRIPTION: Subroutine Init\_tSlc initializes the index list for the 
!  time slices
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_tSlc ( ISIZE, RC )
!
! !ARGUMENTS:
!
      INTEGER, INTENT(IN   )  :: ISIZE
      INTEGER, INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! Init_tSlc begins here!
      !======================================================================

      ! Allocate the Active_tSlc structure
      ALLOCATE ( Active_tSlc ) 

      ! Initialize the indices for the 4th (i.e. the temporal) array dimension
      ALLOCATE ( Active_tSlc%CONSTANT(ISIZE) )
      Active_tSlc%CONSTANT(:)    = 1

      ALLOCATE ( Active_tSlc%HOURLY(ISIZE) )
      Active_tSlc%HOURLY(:)      = 1

      ALLOCATE ( Active_tSlc%HOURLY_GRID(ISIZE) )
      Active_tSlc%HOURLY_GRID(:) = 1

      ALLOCATE ( Active_tSlc%WEEKDY(ISIZE) )
      Active_tSlc%WEEKDY(:)      = 1

      ALLOCATE ( Active_tSlc%WEEKDY_GRID(ISIZE) )
      Active_tSlc%WEEKDY_GRID(:) = 1

      ALLOCATE ( Active_tSlc%MONTHLY(ISIZE) )
      Active_tSlc%MONTHLY(:) = 1

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Init_tSlc
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Set_tSlc
!
! !DESCRIPTION: Subroutine Set\_tSlc defines the appropriate pointer for 
!  the given temporal resolution
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_tSlc ( TimeSlc, TempRes, ISIZE )
!
! !USES:
!
! !ARGUMENTS:
!
      INTEGER,          POINTER     :: TimeSlc(:)
      CHARACTER(LEN=*), INTENT(IN)  :: TempRes
      INTEGER,          INTENT(IN)  :: ISIZE
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! Set_tSlc begins here!
      !======================================================================

      ! Eventually allocate the variables
      IF ( .NOT. ASSOCIATED ( TimeSlc ) ) THEN
         ALLOCATE ( TimeSlc(ISIZE) )
      ENDIF

      ! Set pointer according to TempRes
      SELECT CASE ( TRIM(TempRes) )

         CASE ( 'CONSTANT' )
            TimeSlc => Active_tSlc%CONSTANT

         CASE ( 'HOURLY' )
            TimeSlc => Active_tSlc%HOURLY

         CASE ( 'HOURLY_GRID' )
            TimeSlc => Active_tSlc%HOURLY_GRID

         CASE ( 'WEEKDY' )
            TimeSlc => Active_tSlc%WEEKDY

         CASE ( 'WEEKDY_GRID' )
            TimeSlc => Active_tSlc%WEEKDY_GRID

         CASE ( 'MONTHLY' )
            TimeSlc => Active_tSlc%MONTHLY

         CASE DEFAULT
            TimeSlc => Active_tSlc%CONSTANT

      END SELECT

      END SUBROUTINE Set_tSlc
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: INIT_CLOCK
!
! !DESCRIPTION: Subroutine INIT\_CLOCK initializes the given clock. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_CLOCK ( Clock, RC )
!
! !USES:
!
!
! !ARGUMENTS:
!
      TYPE(HcoClock), POINTER        :: Clock
      INTEGER,        INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  10 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! Init_CLOCK begins here!
      !======================================================================
 
      ! Eventually allocate clock object and set all values to -1
      IF ( .NOT. ASSOCIATED ( Clock ) ) THEN
         ALLOCATE ( Clock )
         Clock%PrevYear     = -1
         Clock%PrevMonth    = -1
         Clock%PrevDay      = -1
         Clock%PrevHour     = -1
         Clock%PrevMin      = -1
         Clock%PrevSec      = -1
         Clock%PrevDOY      = -1
         Clock%PrevWD       = -1

         Clock%ThisYear     = -1
         Clock%ThisMonth    = -1
         Clock%ThisDay      = -1
         Clock%ThisHour     = -1
         Clock%ThisMin      = -1
         Clock%ThisSec      = -1
         Clock%ThisDOY      = -1
         Clock%ThisWD       = -1
         Clock%MonthLastDay = -1
      ENDIF

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE INIT_CLOCK
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: SET_CLOCK
!
! !DESCRIPTION: Subroutine SET\_CLOCK updates the given clock. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_CLOCK ( Clock, HcoState, RC ) 
!
! !USES:
!
      USE HCO_TYPE_MOD, ONLY   : HCO_State
!
! !ARGUMENTS:
!
      TYPE(HcoClock),  POINTER       :: Clock     ! Hemco clock object
      TYPE(HCO_State), POINTER       :: HcoState  ! Hemco state 
      INTEGER,         INTENT(INOUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      INTEGER :: LastDay

      !======================================================================
      ! SET_CLOCK begins here!
      !======================================================================

      ! First set previous time stamps 
      Clock%PrevYear   = Clock%ThisYear
      Clock%PrevMonth  = Clock%ThisMonth
      Clock%PrevDay    = Clock%ThisDay
      Clock%PrevHour   = Clock%ThisHour
      Clock%PrevMin    = Clock%ThisMin
      Clock%PrevSec    = Clock%ThisSec
      Clock%PrevDOY    = Clock%ThisDOY
      Clock%PrevWD     = Clock%ThisWD

      ! Update current time stamps 
      Clock%ThisYear   = HcoState%sYear 
      Clock%ThisMonth  = HcoState%sMonth
      Clock%ThisDay    = HcoState%sDay
      Clock%ThisHour   = HcoState%sHour
      Clock%ThisMin    = HcoState%sMin
      Clock%ThisSec    = HcoState%sSec
      Clock%ThisDOY    = HcoState%sDayOfYear
      Clock%ThisWD     = HcoState%sWeekDay

      ! Get last day of this month (only if month has changed)
      IF ( Clock%ThisMonth /= Clock%PrevMonth ) THEN 

         LastDay = 31 ! TEMPORARY FIX

         ! Select month
         SELECT CASE ( HcoState%sMonth ) 
      
            ! Months with 30 days
            CASE (4,6,9,11)
               LastDay = 30
   
            ! February
            CASE (2)
               LastDay = 28
   
               ! Check for leap years:
               IF ( (MOD(HcoState%sYear,4  ) == 0) .AND. &
                    (MOD(HcoState%sYear,400) /= 0)        ) THEN 
                  LastDay = 29
               ENDIF

            ! Months with 31 days
            CASE DEFAULT
               LastDay = 31                  

         END SELECT 

         ! Update MonthLastDay attribute
         Clock%MonthLastDay = LastDay
      ENDIF

      END SUBROUTINE SET_CLOCK
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Cleanup_tSlc
!
! !DESCRIPTION: Subroutine Cleanup\_tSlc deallocates the time slice
! index list. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_tSlc
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! Cleanup_tSlc begins here!
      !======================================================================

      IF ( ASSOCIATED( Active_tSlc ) ) THEN

         IF ( ASSOCIATED(Active_tSlc%CONSTANT) ) THEN
            DEALLOCATE(Active_tSlc%CONSTANT) 
         ENDIF

         IF ( ASSOCIATED(Active_tSlc%HOURLY) ) THEN
            DEALLOCATE(Active_tSlc%HOURLY) 
         ENDIF

         IF ( ASSOCIATED(Active_tSlc%HOURLY_GRID) ) THEN
            DEALLOCATE(Active_tSlc%HOURLY_GRID) 
         ENDIF

         IF ( ASSOCIATED(Active_tSlc%WEEKDY) ) THEN
            DEALLOCATE(Active_tSlc%WEEKDY) 
         ENDIF

         IF ( ASSOCIATED(Active_tSlc%WEEKDY_GRID) ) THEN
            DEALLOCATE(Active_tSlc%WEEKDY_GRID) 
         ENDIF
  
         IF ( ASSOCIATED(Active_tSlc%MONTHLY) ) THEN
            DEALLOCATE(Active_tSlc%MONTHLY) 
         ENDIF
  
         ! Also deallocate Active_tSlc pointer 
         DEALLOCATE( Active_tSlc )
 
      ENDIF

      END SUBROUTINE Cleanup_tSlc
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: CLEANUP_Clock
!
! !DESCRIPTION: Subroutine CLEANUP\_CLOCK removes the given HcoClock
! type.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_CLOCK ( Clock )
!
! !ARGUMENTS
!
      TYPE(HcoClock), POINTER   :: Clock

!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! CLEANUP_Clock begins here!
      !======================================================================

      IF ( ASSOCIATED( Clock ) ) THEN
         DEALLOCATE ( Clock )
      ENDIF

      END SUBROUTINE CLEANUP_CLOCK
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Update_tSlc
!
! !DESCRIPTION: Subroutine Update\_tSlc makes sure that the correct time
! slices will be used (at every given longitude).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Update_tSlc ( Clock, HcoState, RC )
!
! !USES:
!
      USE HCO_TYPE_MOD,  ONLY : HCO_State
!
! !ARGUMENTS:
!
      TYPE(HcoClock),  POINTER       :: Clock
      TYPE(HCO_State), POINTER       :: HcoState
      INTEGER,         INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  29 Dec 2012 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I,      MtLastDay
      INTEGER            :: Hutc,   Mutc,  Sutc,  WDutc, Dyutc, Mtutc
      INTEGER            :: Hloc,   WDloc, Dyloc, Mtloc 
      REAL*8             :: DECutc, DECloc
      CHARACTER(LEN=255) :: MSG, LOC

      ! testing only
      LOGICAL, SAVE      :: FIRST = .TRUE.

      !======================================================================
      ! Update_tSlc begins here!
      !======================================================================

      ! For error handling
      LOC = 'Update_tScl (HCO_TIME_MOD.F90)'

      ! Get all times in UTC from Clock object
      Mtutc     = Clock%ThisMonth
      Dyutc     = Clock%ThisDay
      Hutc      = Clock%ThisHour
      Mutc      = Clock%ThisMin
      Sutc      = Clock%ThisSec
      WDutc     = Clock%ThisWD
      MtLastDay = Clock%MonthLastDay

      ! UTC decimal time
      DECutc = ( DBLE( Hutc )          ) + &
               ( DBLE( Mutc ) / 60d0   ) + &
               ( DBLE( Sutc ) / 3600d0 )

      ! Loop over longitudes to account for different local times.
      DO I = 1, SIZE(Active_tSlc%HOURLY)

         ! local decimal time 
         DECloc = DECutc + ( HcoState%XMID( I, 1, 1 ) / 15d0 )

         ! Extract local hour and weekday
         ! Case 1: Local time is one day behind UTC.  
         IF ( DECloc < 0d0 ) THEN
            WDloc = WDutc - 1
            IF ( WDloc < 0 ) WDloc = 6
            Hloc  = FLOOR( DECloc + 24d0 )
 
            ! local day one day behind:
            Dyloc = Dyutc - 1 

         ! Case 2: Local time is one day ahead UTC.
         ELSEIF ( DECloc >= 24d0 ) THEN
            WDloc = WDutc + 1
            IF ( WDloc > 6 ) WDloc = 0
            Hloc  = FLOOR( DECloc - 24d0 )

            ! local day one day ahead:
            Dyloc = Dyutc + 1

         ! Case 3: Local time is same day as UTC.
         ELSE
            WDloc = WDutc
            Hloc  = FLOOR( DECloc )

            ! local day as utc day
            Dyloc = Dyutc

         ENDIF

         ! HOURLY:
         ! These are the pointers to the hourly indices. Hourly points
         ! to the hourly time slice representative for the LOCAL time at
         ! longitude ii. 
         ! Add one hour to local hour to make sure that hour 0 becomes 
         ! index 1. 
         Active_tSlc%HOURLY(I) = Hloc + 1

         ! HOURLY_GRID:
         ! Gridded hourly data is assumed to be already adjusted for
         ! local time effects, hence just point to the time slice
         ! of current UTC time.
         Active_tSlc%HOURLY_GRID(I) = Hutc + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!          TEMPORARY FIX FOR CONSISTENCY W/ FORMER VERSION
!
         ! WEEKDY:
         ! For non-gridded factors, take into account local time:
!         Active_tSlc%WEEKDY(I) = WDloc + 1
         Active_tSlc%WEEKDY(I) = WDutc + 1
         IF ( FIRST .AND. I == 1 ) THEN
            MSG =  'Constant weekday used, needs to be fixed!'
            CALL HCO_WARNING( MSG, LOC, RC )
         ENDIF 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         ! WEEKDY_GRID
         ! For gridded weekday factors, just use the UTC slice:
         Active_tSlc%WEEKDY_GRID(I) = WDutc + 1

         ! MONTHLY data
         ! If local day = utc day, then local month is equal to utc
         ! month
         IF ( Dyloc == Dyutc ) THEN
            Mtloc = Mtutc
         
         ! If behind, check if we have to adjust month
         ELSEIF ( Dyloc < Dyutc ) THEN
            IF ( Dyloc == 0 ) THEN
               Mtloc = Mtutc - 1
            ELSE
               Mtloc = Mtutc
            ENDIF
        
            ! If we fall back from Jan to Dec:
            IF ( Mtloc == 0 ) Mtloc = 12
       
         ! If ahead, check if the local time is in a new month
         ELSEIF ( Dyloc > Dyutc ) THEN

            IF ( Dyutc == MtLastDay ) THEN
               Mtloc = Mtutc + 1
            ELSE
               Mtloc = Mtutc
            ENDIF
            
            ! If we enter a new year: 
            IF ( Mtloc > 12 ) Mtloc = 1
         ENDIF

         ! Set month index accordingly
         Active_tSlc%MONTHLY(I) = Mtloc 

      ENDDO !I

      ! Adjust first flag
      FIRST = .FALSE. 

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Update_tSlc
!EOC
      END MODULE HCO_TIME_MOD
!EOM
