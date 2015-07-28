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
      
! 
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE  :: GEOS_Timer_Find      ! Finds the specified timer.
      PRIVATE  :: GEOS_Timer_PrintNum  ! Prints the timer by number.

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
#if defined( TIMERS )

      ! Enable the timers boolean.
      LOGICAL, PARAMETER         :: DoTimers = .true.

#else

      ! Disable the timers boolean.
      LOGICAL, PARAMETER         :: DoTimers = .false.

#endif

      ! What mode the timers should be in. Defaults to 0.
      ! 0: CPU Time
      ! 1: Real Time
      ! 2: MPI time
      INTEGER                    :: TimerMode = 0

      ! Current number of timers. Starts at 0.
      INTEGER                    :: TimerCurrentSize = 0
      ! Maximum Supported Timers. Increasing will increase memory
      ! footprint.
      INTEGER, PARAMETER         :: TimerMaxSize = 15

      ! The definition of the GEOS_Timer type.
      TYPE GEOS_Timer
        LOGICAL                  :: ENABLED
        CHARACTER(LEN=255)       :: TIMER_NAME
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
! !DESCRIPTION: Set up the GEOS_Timer for first use.
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

      TimerMode = TheMode

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
      CHARACTER(LEN=255),   INTENT(IN) :: TimerName    ! Name for timer.
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
      CHARACTER(LEN=255),   INTENT(IN) :: TimerName    ! Name for timer.
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

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
        RETURN                       ! If no timers, do nothing.
      ENDIF

      ! First we must find the specified timer.
      TimerLoc = GEOS_Timer_Find( TimerName )

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
      ! 0: CPU Time
      ! 1: Real Time
      ! 2: MPI time
      IF ( TimerMode .eq. 0 ) THEN
       CALL TIMER( TheTime )
      ENDIF

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
      CHARACTER(LEN=255),   INTENT(IN) :: TimerName    ! Name for timer.
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

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
        RETURN                       ! If no timers, do nothing.
      ENDIF

      TimerLoc = GEOS_Timer_Find( TimerName )

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
      ! 0: CPU Time
      ! 1: Real Time
      ! 2: MPI time
      IF ( TimerMode .eq. 0 ) THEN
       CALL TIMER( TheTime )
      ENDIF

      SavedTimers(TimerLoc)%END_TIME = TheTime

      ! Get the difference to the times
      Diff = SavedTimers(TimerLoc)%END_TIME -
     &       SavedTimers(TimerLoc)%START_TIME

      ! And add difference to current value of total time
      SavedTimers(TimerLoc)%TOTAL_TIME =
     &   SavedTimers(TimerLoc)%END_TIME + Diff

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
! !DESCRIPTION: Prints the specified GEOS_Timer by name.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GEOS_Timer_Print( TimerName ) RESULT ( RC )
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=255),   INTENT(IN) :: TimerName    ! Name for timer.
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
      INTEGER    ::  TimerLoc           ! Timer number

      ! First it is important to check if timers are enabled.
      IF (.not. DoTimers) THEN       ! Assume False for optimization.
        RETURN                       ! If no timers, do nothing.
      ENDIF

      TimerLoc = GEOS_Timer_Find( TimerName )

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
! !DESCRIPTION: Prints all GEOS_Timers.
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

! !IROUTINE: GEOS_Timer_PrintNum
!
! !DESCRIPTION: Prints GEOS_Timer by number.
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

      PRINT*, "Total Running time: ",
     &  SavedTimers(SlotNumber) % TOTAL_TIME, " seconds"

      END SUBROUTINE GEOS_Timer_PrintNum


!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOS_Timer_Find
!
! !DESCRIPTION: Finds the number of the specified GEOS_Timer.
!\\
!\\
! !INTERFACE:
!
      FUNCTION GEOS_Timer_Find( TimerName ) RESULT ( SlotNumber ) 
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=255),   INTENT(IN) :: TimerName    ! Name for timer.
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

      IF(TimerCurrentSize < 1) THEN  ! Return 0 if it's empty
        RETURN
      ENDIF

      DO I = 1, TimerCurrentSize
       IF(SavedTimers(TimerCurrentSize) % TIMER_NAME 
     &  .eq. TimerName) THEN
           SlotNumber = I
       ENDIF
      ENDDO

      END FUNCTION GEOS_Timer_Find
!EOC
      END MODULE GEOS_TIMERS
