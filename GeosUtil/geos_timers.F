!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: geos_timers
!
! !DESCRIPTION: Module GEOS\_TIMERS is used to track and time how long
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
      PUBLIC  :: GEOS_Timer_Setup     ! Init Method
      PUBLIC  :: GEOS_Timer_Add       ! Adds a timer.
      PUBLIC  :: GEOS_Timer_Start     ! Starts a timer ticking.
      PUBLIC  :: GEOS_Timer_End       ! Stops a timer ticking.
      PUBLIC  :: GEOS_Timer_Print     ! Prints the specified timer.
      PUBLIC  :: GEOS_Timer_PrintAll  ! Prints all timers.
      PUBLIC  :: GEOS_Timer_StopAll   ! Stops all currently running timers.
! 
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: GEOS_Timer_Find      ! Finds the specified timer.
      PRIVATE :: GEOS_Timer_PrintNum  ! Prints the timer by number.
      PRIVATE :: GEOS_Timer_TheTime   ! Returns the current time in MS.
      PRIVATE :: GEOS_Timer_TimePrint ! Formats the seconds when printing.
!
! !REMARKS:
!  This module helps track valuable timing information.
!
! !REVISION HISTORY:
!  23 Jul 2015 - M. Yannetti - Initial version.
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
      INTEGER, INTENT(IN)  :: TheMode     ! Timer mode
                                          ! 1:CPU time, 2:Real time, 3:MPI time
!
! !REMARKS:
!  This currently only needs to run if you want to manually set the mode.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!  27 Oct 2015 - M. Sulprizio- Change from a function to a subroutine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
         RETURN                      ! If no timers, do nothing.
      ENDIF

      IF (TheMode .lt. 1) THEN
         PRINT*, "GEOS_Timer_Setup: WARNING: INVALID TIMER TYPE! "
     &           "The following timer modes are supported: "
     &           "(1) CPU time, (2) Real time, or (3) MPI time."
         RETURN
      ENDIF

      TimerMode = TheMode

      ! Debug
      !PRINT*, "GEOS_Timer_Setup: Done setting up GEOS-Chem timers"

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
      SUBROUTINE GEOS_Timer_Add( TimerName, RC )
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN)    :: TimerName   ! Name for timer.
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(INOUT) :: RC          ! Success / Failure

!
! !REMARKS:
!  This only fails if the timers are full.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!  27 Oct 2015 - M. Sulprizio- Change from a function to a subroutine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
         RC = GIGC_SUCCESS
         RETURN                      ! If no timers, do nothing.
      ENDIF

      ! Now we are sure that timers are enabled.
      ! We need to check if the timers are full.
      IF (TimerCurrentSize < TimerMaxSize) THEN         ! There's room.

         ! Increase the timer current size by one.
         TimerCurrentSize = TimerCurrentSize + 1

         ! Set the defaults of the new Timer.
         SavedTimers(TimerCurrentSize)%ENABLED    = .false.
         SavedTimers(TimerCurrentSize)%TIMER_NAME = TimerName
         SavedTimers(TimerCurrentSize)%TOTAL_TIME = 0
         SavedTimers(TimerCurrentSize)%START_TIME = 0
         SavedTimers(TimerCurrentSize)%END_TIME   = 0

         ! Debug
         !PRINT*, TimerName, "timer added at slot ", TimerCurrentSize

         ! Success.
         RC = GIGC_SUCCESS

      ELSE                                             ! There's not room.

         PRINT*," ** WARNING: Max timers reached! **"
         PRINT*,"    TimerCurrentSize = ", TimerCurrentSize
         PRINT*,"    TimerMaxSize     = ", TimerMaxSize

         ! Failure
         RC = GIGC_FAILURE

      ENDIF

      END SUBROUTINE GEOS_Timer_Add
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
      SUBROUTINE GEOS_Timer_Start( TimerName, RC )
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN)    :: TimerName   ! Name for timer.
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(INOUT) :: RC          ! Success / Failure
!
! !REMARKS:
!  This must be called to start a timer ticking.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!  27 Oct 2015 - M. Sulprizio- Change from a function to a subroutine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: TimerLoc           ! Timer number
      REAL*8            :: TheTime            ! Returned Time from method
      CHARACTER(LEN=30) :: TempTimerName

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
         RC = GIGC_SUCCESS
         RETURN                      ! If no timers, do nothing.
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
      IF ( (SavedTimers(TimerLoc)%ENABLED) .eq. .true.) THEN
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
      
      ! Debug
      !PRINT*, "** RETURNED TIME (START): ", TheTime

      SavedTimers(TimerLoc)%START_TIME = TheTime

      ! Success.
      RC = GIGC_SUCCESS

      END SUBROUTINE GEOS_Timer_Start
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
      SUBROUTINE GEOS_Timer_End( TimerName, RC )
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN)    :: TimerName   ! Name for timer.
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(INOUT) :: RC          ! Success / Failure
!
! !REMARKS:
!  Without this routine being called, a timer will not add to its total.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!  27 Oct 2015 - M. Sulprizio- Change from a function to a subroutine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      INTEGER           :: TimerLoc           ! Timer number
      REAL*8            :: TheTime            ! Returned Time from method
      REAL*8            :: Diff               ! Difference in times
      CHARACTER(LEN=30) :: TempTimerName

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
         RC = GIGC_SUCCESS
         RETURN                      ! If no timers, do nothing.
      ENDIF

      TempTimerName = TimerName

      TimerLoc = GEOS_Timer_Find( TempTimerName )

      IF (TimerLoc .eq. 0) THEN
         PRINT*, "** WARNING: Timer not found: ", TimerName
         RC = GIGC_FAILURE
         RETURN
      ENDIF

      ! Now we do some minor error checking
      IF ( (SavedTimers(TimerLoc)%ENABLED) .eq. .false.) THEN
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

      ! Debug
      !PRINT*, "** RETURNED TIME (END): ", TheTime

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

      END SUBROUTINE GEOS_Timer_End
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
      SUBROUTINE GEOS_Timer_Print( TimerName, am_I_Root, RC )
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN)    :: TimerName   ! Name for timer.
      LOGICAL,          INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(INOUT) :: RC          ! Success / Failure
!
! !REMARKS:
!  This is useful if you only want to print a single timer.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!  27 Oct 2015 - M. Sulprizio- Change from a function to a subroutine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      CHARACTER(LEN=30) :: TempTimerName
      INTEGER           :: TimerLoc           ! Timer number

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
         RC = GIGC_SUCCESS
         RETURN                      ! If no timers, do nothing.
      ENDIF

      TempTimerName = TimerName

      TimerLoc = GEOS_Timer_Find( TempTimerName )

      IF (TimerLoc .eq. 0) THEN
         PRINT*, "** WARNING: Timer not found: ", TimerName
         RC = GIGC_FAILURE
         RETURN
      ENDIF      

      CALL GEOS_Timer_PrintNum( TimerLoc, am_I_Root )

      ! Success.
      RC = GIGC_SUCCESS

      END SUBROUTINE GEOS_Timer_Print
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_Timer_PrintAll
!
! !DESCRIPTION: Prints all GEOS\_Timers to log file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GEOS_Timer_PrintAll( am_I_Root, RC )
!
! !INPUT PARAMETERS: 
!
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Is this the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success / Failure
!
! !REMARKS:
!  This prints all timers in the order added.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!  27 Oct 2015 - M. Sulprizio- Change from a function to a subroutine and
!                              modify to print timers out in a table
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
         RC = GIGC_SUCCESS
         RETURN                      ! If no timers, do nothing.
      ENDIF

      IF(TimerCurrentSize < 1) THEN  ! Return 0 if it's empty
         RC = GIGC_FAILURE
         RETURN
      ENDIF

      IF ( am_I_Root ) THEN
         WRITE( 6, *     ) ''
         WRITE( 6, '(a)' ) REPEAT( '=', 79 )
         WRITE( 6, '(a)' ) 'G E O S - C H E M   T I M E R S'
         WRITE( 6, *     ) ''
         WRITE( 6, 100   ) 'Timer name','DD-hh:mm:ss.SSS'
         WRITE( 6, '(a)' ) REPEAT( '-', 55 )
 100     FORMAT( 2x,a10,23x,a15 )
      ENDIF


      DO I = 1, TimerCurrentSize
         CALL GEOS_Timer_PrintNum( I, am_I_Root )
      ENDDO

      ! Success.
      RC = GIGC_SUCCESS

      END SUBROUTINE GEOS_Timer_PrintAll
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
      SUBROUTINE GEOS_Timer_StopAll( RC )
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success / Failure
!
! !REMARKS:
!  This stops all currently running timers. Used during crashes.
! 
! !REVISION HISTORY: 
!  11 Aug 2015 - M. Yannetti - Initial version.
!  27 Oct 2015 - M. Sulprizio- Change from a function to a subroutine
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
         RC = GIGC_SUCCESS
         RETURN                      ! If no timers, do nothing.
      ENDIF

      IF(TimerCurrentSize < 1) THEN  ! Return 0 if it's empty
         RC = GIGC_FAILURE
         RETURN
      ENDIF

      DO I = 1, TimerCurrentSize
         IF ( (SavedTimers(I)%ENABLED) ) THEN
            PRINT*, "Timer forced to stop due to error: ",
     &               SavedTimers(I)%TIMER_NAME

            ! Yes, this is inefficient. Should have another function
            ! written eventually to replace using the normal one.
            CALL GEOS_Timer_End( SavedTimers(I)%TIMER_NAME, RC )
         ENDIF
      ENDDO

      ! Success.
      RC = GIGC_SUCCESS

      END SUBROUTINE GEOS_Timer_StopAll
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
      SUBROUTINE GEOS_Timer_PrintNum( SlotNumber, am_I_Root )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: SlotNumber  ! The slot of the timer.
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !REMARKS:
!  This actually does the printing, and is called by other print
!  routines.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!  27 Oct 2015 - M. Sulprizio- Change from a function to a subroutine
!EOP
!------------------------------------------------------------------------------
!BOC
!
      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
         RETURN                      ! If no timers, do nothing.
      ENDIF

      IF(TimerCurrentSize < 1) THEN  ! Return if it's empty
         RETURN
      ENDIF
      
      IF ( (SavedTimers(SlotNumber)%ENABLED) ) THEN
         PRINT*, "** WARNING: Timer still enabled! "
      ENDIF

      CALL GEOS_Timer_TimePrint( SlotNumber, am_I_Root )

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
      INTEGER :: I

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
         RETURN                      ! If no timers, do nothing.
      ENDIF

      SlotNumber = 0

      IF(TimerCurrentSize .lt. 1) THEN  ! Return 0 if it's empty
         RETURN
      ENDIF

      DO I = 1, TimerCurrentSize, 1
         IF((SavedTimers(I)%TIMER_NAME) .eq. TimerName) THEN
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
! !RETURN VALUE:
!
      REAL*8 :: TotalTime  ! The current calculated time.

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

      ! Debug
      !PRINT*, "TIME_VALUE: ", TIME_VALUE
      !PRINT*, "TIME_CLOCK: ", TIME_CLOCK
      !CALL FLUSH(6)

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
! !DESCRIPTION: Formats the time and prints it out to the log file.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GEOS_Timer_TimePrint( SlotNumber, am_I_Root )
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: SlotNumber  ! The slot of the timer.
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
!
! !REMARKS:
!  This is a private subroutine.
! 
! !REVISION HISTORY: 
!  24 Jul 2015 - M. Yannetti - Initial version.
!  27 Oct 2015 - M. Sulprizio- Change from a function to a subroutine and
!                              modify to print timers out in a table in the
!                              DD-hh:mm:ss.SSS format 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL*8             :: InputSecs       ! Real value of secones
      REAL*8             :: SecsLeft        ! How many seconds are 'left'
      INTEGER            :: IntSecs         ! Same as above, but integer
      INTEGER            :: TotalMS    = 0  ! Total Milliseconds
      INTEGER            :: TotalSecs  = 0  ! Total Seconds
      INTEGER            :: TotalMins  = 0  ! Total Minutes
      INTEGER            :: TotalHours = 0  ! Total Hours
      INTEGER            :: TotalDays  = 0  ! Total Days
      CHARACTER(LEN=100) :: OutputStr       ! Combined output string
      CHARACTER(LEN=10)  :: TempStr         ! Needed to remove whitespace.
      CHARACTER(LEN=2)   :: DD, HH, MM, SS
      CHARACTER(LEN=3)   :: MS

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
        RETURN                       ! If no timers, do nothing.
      ENDIF

      ! Copy the timer value
      InputSecs = SavedTimers(SlotNumber)%TOTAL_TIME
      IntSecs   = INT(InputSecs)
      SecsLeft  = InputSecs - REAL(IntSecs)

      IF ( InputSecs < 0 ) THEN ! Invalid time
         WRITE( 6, 110 ), SavedTimers(SlotNumber)%TIMER_NAME
 110     FORMAT(2x,a30,': Invalid run time - negative value')
         RETURN
      ELSEIF ( InputSecs .eq. 0 ) THEN ! Zero timer
         WRITE( 6, 120 ), SavedTimers(SlotNumber)%TIMER_NAME
 120     FORMAT(2x,a30,': The timer did not run')
         RETURN
      ENDIF

      ! Debug
      !PRINT*, "INT   : ", IntSecs
      !PRINT*, "REAL  : ", InputSecs
      !PRINT*, "REMAIN: ", SecsLeft

      !----------------------------
      ! Calculate hours
      !----------------------------
      TotalHours = FLOOR(REAL(IntSecs)/3600.0)
      IntSecs    = IntSecs - (TotalHours*3600)

      !----------------------------
      ! Calculate days (if needed)
      !----------------------------
      IF ( TotalHours > 24 ) THEN
         TotalDays  = FLOOR(REAL(TotalHours)/24.0)
         TotalHours = TotalHours - (TotalDays*24)
      ENDIF

      !----------------------------
      ! Calculate minutes
      !----------------------------
      TotalMins  = FLOOR(REAL(IntSecs)/60.0)
      IntSecs    = IntSecs - (TotalMins*60)

      !----------------------------
      ! Calculate seconds
      !----------------------------
      TotalSecs  = IntSecs

      !----------------------------
      ! Calculate milliseconds
      !----------------------------
      SecsLeft = SecsLeft * 1000
      TotalMS = INT(SecsLeft)

      WRITE( DD, '(i2.2)' ) TotalDays
      WRITE( HH, '(i2.2)' ) TotalHours
      WRITE( MM, '(i2.2)' ) TotalMins
      WRITE( SS, '(i2.2)' ) TotalSecs
      WRITE( MS, '(i3.3)' ) TotalMS

      WRITE( 6, 130 ) SavedTimers(SlotNumber)%TIMER_NAME,
     &                DD, HH, MM, SS, MS
 130  FORMAT(2x,a30,':',2x,a2,'-',a2,':',a2,':',a2,'.',a3)

      END SUBROUTINE GEOS_Timer_TimePrint

!EOC
      END MODULE GEOS_TIMERS
