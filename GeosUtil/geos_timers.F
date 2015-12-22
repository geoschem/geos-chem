!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: geos_timers
!
! !DESCRIPTION: Module geos\_timers is used to track and time how long
! specified parts of GEOS-Chem take to run.
!\\
!\\
! !INTERFACE:
!
      MODULE GEOS_TIMERS
!
! !USES:
! 
      USE GIGC_Errcode_Mod

      IMPLICIT NONE
      PRIVATE
! 
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC   :: GEOS_Timer_Setup     ! Init Method
      PUBLIC   :: GEOS_Timer_Add       ! Adds a timer.
      PUBLIC   :: GEOS_Timer_Start     ! Starts a timer ticking.
      PUBLIC   :: GEOS_Timer_End       ! Stops a timer ticking.
      PUBLIC   :: GEOS_Timer_Print     ! Prints the specified timer.
      PUBLIC   :: GEOS_Timer_PrintAll  ! Prints all timers.
      PUBLIC   :: GEOS_Timer_StopAll   ! Stops all currently running timers.
      
! 
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE  :: GEOS_Timer_Find      ! Finds the specified timer.
      PRIVATE  :: GEOS_Timer_PrintNum  ! Prints the timer by number.
      PRIVATE  :: GEOS_Timer_TheTime   ! Returns the current time in MS.
      PRIVATE  :: GEOS_Timer_TimePrint ! Formats the seconds when printing.
!
! !DEFINED PARAMETERS:
! 

!
! !REMARKS:
!  This module helps track valuable timing information.
!
! !REVISION HISTORY:
!  (1 ) Created. (myannetti, 07/23/15)
!EOP
!-----------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
#if defined( USE_TIMERS )

      ! Enable the timers boolean.
      LOGICAL, PARAMETER         :: DoTimers = .true.

#else

      ! Disable the timers boolean.
      LOGICAL, PARAMETER         :: DoTimers = .false.

#endif

      ! What mode the timers should be in. Defaults to 1.
      ! 1: CPU Time
      ! 2: Real Time
      ! 3: MPI time
      INTEGER                    :: TimerMode = 1

      ! Current number of timers. Starts at 0.
      INTEGER                    :: TimerCurrentSize = 0
      ! Maximum Supported Timers. Increasing will increase memory
      ! footprint.
      INTEGER, PARAMETER         :: TimerMaxSize = 15

      ! The definition of the GEOS_Timer type.
      TYPE GEOS_Timer
        LOGICAL                  :: ENABLED
        CHARACTER(LEN=30)        :: TIMER_NAME
        REAL*8                   :: TOTAL_TIME
        REAL*8                   :: START_TIME
        REAL*8                   :: END_TIME
      END TYPE GEOS_Timer

      ! The array of timers. Determined by TimerMaxSize.
      TYPE(GEOS_Timer), DIMENSION(TimerMaxSize) :: SavedTimers


      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_Timer_Setup
!
! !DESCRIPTION: Set up the GEOS\_Timer for first use.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GEOS_Timer_Setup( TheMode )
!
! !INPUT PARAMETERS: 
!
      INTEGER,     INTENT(IN)  :: TheMode      ! The mode the timers should be.
!
! !OUTPUT PARAMETERS:
!

!
! !REMARKS:
!  This currently only needs to run if you want to manually set the mode.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
        RETURN                       ! If no timers, do nothing.
      ENDIF

      IF (TheMode .lt. 1) THEN
        PRINT*, "** WARNING: INVALID TIMER TYPE!"
        RETURN
      ENDIF

      TimerMode = TheMode
      PRINT*, "** GEOS TIMERS SET UP!"
      END SUBROUTINE GEOS_Timer_Setup
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_Timer_Add
!
! !DESCRIPTION: Adds a new timer to the timer list.
!  Returns status of success.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GEOS_Timer_Add( TimerName ) RESULT ( RC )
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*),   INTENT(IN) :: TimerName    ! Name for timer.
!
! !OUTPUT PARAMETERS:
!
      INTEGER    :: RC          ! Success / Failure

!
! !REMARKS:
!  This only fails if the timers are full.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
        RC = GIGC_SUCCESS
        RETURN                       ! If no timers, do nothing.
      ENDIF

      ! Now we are sure that timers are enabled.
      ! We need to check if the timers are full.
      IF (TimerCurrentSize < TimerMaxSize) THEN         ! There's room.
        ! Increase the timer current size by one.
        TimerCurrentSize = TimerCurrentSize + 1

        ! Set the defaults of the new Timer.
        SavedTimers(TimerCurrentSize) % ENABLED    = .false.
        SavedTimers(TimerCurrentSize) % TIMER_NAME = TimerName
        SavedTimers(TimerCurrentSize) % TOTAL_TIME = 0
        SavedTimers(TimerCurrentSize) % START_TIME = 0
        SavedTimers(TimerCurrentSize) % END_TIME   = 0

        ! Success.
!        PRINT*, "Timer added at slot ", TimerCurrentSize
        RC = GIGC_SUCCESS
      ELSE                                              ! There's not room.
        PRINT*," ** WARNING: Max timers reached!"
        RC = GIGC_FAILURE                               ! Failure
      ENDIF

      END FUNCTION GEOS_Timer_Add
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_Timer_Start
!
! !DESCRIPTION: Starts a timer ticking.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GEOS_Timer_Start( TimerName ) RESULT ( RC )
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*),   INTENT(IN) :: TimerName    ! Name for timer.
!
! !OUTPUT PARAMETERS:
!
      INTEGER   :: RC          ! Success / Failure

!
! !REMARKS:
!  This must be called to start a timer ticking.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER    ::  TimerLoc           ! Timer number
      REAL*8     ::  TheTime            ! Returned Time from method
      CHARACTER(LEN=30)   :: TempTimerName

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
        RC = GIGC_SUCCESS
        RETURN                       ! If no timers, do nothing.
      ENDIF

      TempTimerName = TimerName

      ! First we must find the specified timer.
      TimerLoc = GEOS_Timer_Find( TempTimerName )

      IF (TimerLoc .eq. 0) THEN
        PRINT*, "** WARNING: Timer not found: ", TimerName
        RC = GIGC_FAILURE
        RETURN
      ENDIF

      ! Now we do some minor error checking
      IF ( (SavedTimers(TimerLoc) % ENABLED) .eq. .true.) THEN
        PRINT*, "** WARNING: Timer already running: ", TimerName
        RC = GIGC_FAILURE
        RETURN
      ENDIF
      
      ! Timer isn't enabled, it's been found, so we enable it
      SavedTimers(TimerLoc)%ENABLED = .true.

      ! And we note the current time
      ! 1: CPU Time
      ! 2: Real Time
      ! 3: MPI time
      IF ( TimerMode .eq. 1 ) THEN
       TheTime = GEOS_Timer_TheTime()
      ENDIF
      
!     PRINT*, "** RETURNED TIME (START): ", TheTime

      SavedTimers(TimerLoc)%START_TIME = TheTime

      ! Success.
      RC = GIGC_SUCCESS

      END FUNCTION GEOS_Timer_Start
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_Timer_End
!
! !DESCRIPTION: Stops a timer ticking. Adds elapsed time to total.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GEOS_Timer_End( TimerName ) RESULT ( RC )
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*),   INTENT(IN) :: TimerName    ! Name for timer.
!
! !OUTPUT PARAMETERS:
!
      INTEGER    :: RC          ! Success / Failure

!
! !REMARKS:
!  Without this routine being called, a timer will not add to its total.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER    ::  TimerLoc           ! Timer number
      REAL*8     ::  TheTime            ! Returned Time from method
      REAL*8     ::  Diff               ! Difference in times
      CHARACTER(LEN=30)   :: TempTimerName

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
        RC = GIGC_SUCCESS
        RETURN                       ! If no timers, do nothing.
      ENDIF

      TempTimerName = TimerName

      TimerLoc = GEOS_Timer_Find( TempTimerName )

      IF (TimerLoc .eq. 0) THEN
        PRINT*, "** WARNING: Timer not found: ", TimerName
        RC = GIGC_FAILURE
        RETURN
      ENDIF

      ! Now we do some minor error checking
      IF ( (SavedTimers(TimerLoc) % ENABLED) .eq. .false.) THEN
        PRINT*, "** WARNING: Timer not running: ", TimerName
        RC = GIGC_FAILURE
        RETURN
      ENDIF

      ! Timer is enabled, it's been found, so we disable it
      SavedTimers(TimerLoc)%ENABLED = .false.

      ! And we note the current time
      ! 1: CPU Time
      ! 2: Real Time
      ! 3: MPI time
      IF ( TimerMode .eq. 1 ) THEN
       TheTime = GEOS_Timer_TheTime()
      ENDIF

!      PRINT*, "** RETURNED TIME (END): ", TheTime

      SavedTimers(TimerLoc)%END_TIME = TheTime

      ! Get the difference to the times
      Diff = SavedTimers(TimerLoc)%END_TIME -
     &       SavedTimers(TimerLoc)%START_TIME

      ! Error check...
      IF ( Diff .lt. 0 ) THEN
       PRINT*, "** WARNING: Timer returned invalid value: ", TimerName
       Diff = 0
      ENDIF

      ! And add difference to current value of total time
      SavedTimers(TimerLoc)%TOTAL_TIME =
     &   (SavedTimers(TimerLoc)%TOTAL_TIME) + Diff

      ! Success.
      RC = GIGC_SUCCESS

!
      END FUNCTION GEOS_Timer_End
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_Timer_Print
!
! !DESCRIPTION: Prints the specified GEOS\_Timer by name.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GEOS_Timer_Print( TimerName ) RESULT ( RC )
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*),   INTENT(IN) :: TimerName    ! Name for timer.
!
! !OUTPUT PARAMETERS:
!
      INTEGER    :: RC          ! Success / Failure
!
! !REMARKS:
!  This is useful if you only want to print a single timer.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=30)   :: TempTimerName
      INTEGER    ::  TimerLoc           ! Timer number

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
        RC = GIGC_SUCCESS
        RETURN                       ! If no timers, do nothing.
      ENDIF

      TempTimerName = TimerName

      TimerLoc = GEOS_Timer_Find( TempTimerName )

      IF (TimerLoc .eq. 0) THEN
        PRINT*, "** WARNING: Timer not found: ", TimerName
        RC = GIGC_FAILURE
        RETURN
      ENDIF      

      CALL GEOS_Timer_PrintNum(TimerLoc)

      ! Success.
      RC = GIGC_SUCCESS

      END FUNCTION GEOS_Timer_Print
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_Timer_PrintAll
!
! !DESCRIPTION: Prints all GEOS\_Timers.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GEOS_Timer_PrintAll() RESULT ( RC )
!
! !INPUT PARAMETERS: 
!

!
! !OUTPUT PARAMETERS:
!
      INTEGER           :: RC          ! Success / Failure


!
! !REMARKS:
!  This prints all timers in the order added.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I            ! For the for loop

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
        RC = GIGC_SUCCESS
        RETURN                       ! If no timers, do nothing.
      ENDIF

      IF(TimerCurrentSize < 1) THEN  ! Return 0 if it's empty
        RC = GIGC_FAILURE
        RETURN
      ENDIF

      DO I = 1, TimerCurrentSize
        CALL GEOS_Timer_PrintNum(I)
      ENDDO

      ! Success.
      RC = GIGC_SUCCESS

      END FUNCTION GEOS_Timer_PrintAll
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_Timer_StopAll
!
! !DESCRIPTION: Stops all GEOS\_Timers.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GEOS_Timer_StopAll() RESULT ( RC )
!
! !INPUT PARAMETERS: 
!

!
! !OUTPUT PARAMETERS:
!
      INTEGER           :: RC          ! Success / Failure

!
! !REMARKS:
!  This stops all currently running timers. Used during crashes.
! 
! !REVISION HISTORY: 
!  11 Aug 2015 - M. Yannetti - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER           :: I            ! For the for loop

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
        RC = GIGC_SUCCESS
        RETURN                       ! If no timers, do nothing.
      ENDIF

      IF(TimerCurrentSize < 1) THEN  ! Return 0 if it's empty
        RC = GIGC_FAILURE
        RETURN
      ENDIF

      DO I = 1, TimerCurrentSize
        IF ( (SavedTimers(I) % ENABLED) ) THEN
         ! Note it
         PRINT*, "Timer forced to stop due to error: ",
     &        SavedTimers(I) % TIMER_NAME
         ! Yes, this is inefficient. Should have another function
         ! written eventually to replace using the normal one.
         RC = GEOS_Timer_End( SavedTimers(I) % TIMER_NAME )
        ENDIF
      ENDDO

      ! Success.
      RC = GIGC_SUCCESS

      END FUNCTION GEOS_Timer_StopAll

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!

! !IROUTINE: GEOS_Timer_PrintNum
!
! !DESCRIPTION: Prints GEOS\_Timer by number.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GEOS_Timer_PrintNum( SlotNumber )
!
! !INPUT PARAMETERS: 
!
      INTEGER,              INTENT(IN) :: SlotNumber  ! The slot of the timer.

!
! !OUTPUT PARAMETERS:
!

!
! !REMARKS:
!  This actually does the printing, and is called by other print
!  routines.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
        RETURN                       ! If no timers, do nothing.
      ENDIF

      IF(TimerCurrentSize < 1) THEN  ! Return if it's empty
        RETURN
      ENDIF
      
      PRINT*, "Timer data for timer name: ",
     &  SavedTimers(SlotNumber) % TIMER_NAME

      IF ( (SavedTimers(SlotNumber) % ENABLED) ) THEN
       PRINT*, "** WARNING: Timer still enabled! "
      ENDIF

      CALL GEOS_Timer_TimePrint (
     &      SavedTimers(SlotNumber) % TOTAL_TIME )

      END SUBROUTINE GEOS_Timer_PrintNum


!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_Timer_Find
!
! !DESCRIPTION: Finds the number of the specified GEOS\_Timer.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GEOS_Timer_Find( TimerName ) RESULT ( SlotNumber ) 
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=30),   INTENT(IN) :: TimerName    ! Name for timer.
!
! !RETURN VALUE:
!
      INTEGER           :: SlotNumber  ! The slot of the timer.

!
! !REMARKS:
!  This is a private routine.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I            ! For the for loop

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
        RETURN                       ! If no timers, do nothing.
      ENDIF

      SlotNumber = 0

      IF(TimerCurrentSize .lt. 1) THEN  ! Return 0 if it's empty
        RETURN
      ENDIF

      DO I = 1, TimerCurrentSize, 1
       IF((SavedTimers(I) % TIMER_NAME) .eq. TimerName) THEN
           SlotNumber = I
       ENDIF
      ENDDO

      END FUNCTION GEOS_Timer_Find
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_Timer_TheTime
!
! !DESCRIPTION: Returns the current time in MS.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GEOS_Timer_TheTime() RESULT ( TotalTime )
!
! !INPUT PARAMETERS: 
!

!
! !RETURN VALUE:
!
      REAL*8           :: TotalTime  ! The current calculated time.

!
! !REMARKS:
!  This is a private routine.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER        :: TIME_VALUE            ! For the function
      INTEGER        :: TIME_CLOCK            ! For the function

      ! First it is important to check if timers are enabled.
!      IF (.not. DoTimers) THEN       ! Assume False for optimization.
!        RETURN                       ! If no timers, do nothing.
!      ENDIF

      ! Let's call the intrinsic function...
      CALL SYSTEM_CLOCK(TIME_VALUE, TIME_CLOCK)
!      PRINT*, "TIME_VALUE: ", TIME_VALUE
!      PRINT*, "TIME_CLOCK: ", TIME_CLOCK
!      CALL FLUSH(6)
      TotalTime = REAL(TIME_VALUE) / REAL(TIME_CLOCK)

      END FUNCTION GEOS_Timer_TheTime

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!

! !IROUTINE: GEOS_Timer_TimePrint
!
! !DESCRIPTION: Formats the printing of the seconds.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GEOS_Timer_TimePrint( InputSecs )
!
! !INPUT PARAMETERS: 
!
      REAL*8,              INTENT(IN) :: InputSecs  ! Real value of seconds

!
! !OUTPUT PARAMETERS:
!

!
! !REMARKS:
!  This is a private subroutine.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8    :: SecsLeft    ! How many seconds are 'left'
      INTEGER   :: IntSecs     ! Same as above, but integer
      INTEGER   :: TotalMS    = 0  ! Total Milliseconds
      INTEGER   :: TotalSecs  = 0  ! Total Seconds
      INTEGER   :: TotalMins  = 0  ! Total Minutes
      INTEGER   :: TotalHours = 0  ! Total Hours
      CHARACTER(LEN=100) :: OutputStr ! Combined output string
      CHARACTER(LEN=10)  :: TempStr   ! This is needed to remove whitespace.

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
        RETURN                       ! If no timers, do nothing.
      ENDIF

      ! Copy the value over
      IntSecs = INT(InputSecs)
      SecsLeft = InputSecs - REAL(IntSecs)

      IF (IntSecs < 0) THEN ! Invalid time
        PRINT*, "Invalid running time: negative value!"
        RETURN
      ELSEIF (SecsLeft .eq. 0) THEN ! Zero timer
        PRINT*, "The timer did not run."
        RETURN
      ENDIF

!      PRINT*, "INT: ", IntSecs
!      PRINT*, "REAL: ", InputSecs
!      PRINT*, "REMAIN: ", SecsLeft

      ! Start the calculation loop
      ! 1 hour = 60 minutes = 3600 seconds
      TotalHours = FLOOR(REAL(IntSecs)/3600.0)
      IntSecs = IntSecs - (TotalHours*3600)
      TotalMins = FLOOR(REAL(IntSecs)/60.0)
      IntSecs = IntSecs - (TotalMins*60)
      TotalSecs = IntSecs
      ! and the milliseconds
      SecsLeft = SecsLeft * 1000
      TotalMS = INT(SecsLeft)

      ! This big ugly part is entirely to make the output look nicer.
      OutputStr = "Total Running time: "

      IF (TotalHours > 0) THEN
       WRITE(TempStr, "(I10)") TotalHours
       OutputStr = TRIM(OutputStr) // " " // TRIM(ADJUSTL(TempStr))
       OutputStr = TRIM(OutputStr) // " hours, "
      ENDIF
      IF (TotalMins > 0) THEN
       WRITE(TempStr, "(I10)") TotalMins
       OutputStr = TRIM(OutputStr) // " " // TRIM(ADJUSTL(TempStr)) 
       OutputStr = TRIM(OutputStr) // " minutes, "
      ENDIF
      IF (TotalSecs > 0) THEN
       WRITE(TempStr, "(I10)") TotalSecs
       OutputStr = TRIM(OutputStr) // " " // TRIM(ADJUSTL(TempStr)) 
       OutputStr = TRIM(OutputStr) // " seconds, and "
      ENDIF
      WRITE(TempStr, "(I10)") TotalMS
      OutputStr = TRIM(OutputStr) // " " //  TRIM(ADJUSTL(TempStr)) 
      OutputStr = TRIM(OutputStr) // " ms."

      PRINT*, TRIM(OutputStr)

      END SUBROUTINE GEOS_Timer_TimePrint

!EOC
      END MODULE GEOS_TIMERS
