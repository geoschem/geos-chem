!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_error_mod 
!
! !DESCRIPTION: Module HCO\_ERROR\_MOD contains routines and variables 
! for error handling in HEMCO. It also contains definitions of some 
! globally used parameter, such as the single/double precision 
! definitions.
! The error settings are specified in the HEMCO configuration file and
! error handling is performed according to these settings. They include: 
! (1) HEMCO logfile: all HEMCO information is written into the specified
!     logfile. The logfile can be set to the wildcard character, in which
!     case the standard output will be used (this may be another opened
!     logfile!).
! (2) Verbose: if true, this will run HEMCO in verbose mode. Note that
!     this may significantly slow down the model!
! (3) Track: if true, this will print the current location in the code
!     into the logfile. This is primarily for debugging.
! (4) Wildcard: wildcard character used in the HEMCO input file, e.g.
!     for time stamps.
! (5) Separator: separator character used in the HEMCO input file, e.g.
!     to separate scale factors and time stamp elements. 
! (6) Show warnings: if TRUE, prompt all warnings to the HEMCO logfile.
!     Otherwise, no warnings will be prompted but the total number of
!     warnings occurred will still be shown at the end of the run.
! (7) Only unitless scale factors: If set to TRUE, this will force all
!     scale factors to be 'unitless', as specified in module 
!     HCO\_UNIT\_MOD (code will return with error if scale factor is 
!     not unitless).
! \\
! !INTERFACE: 
!
MODULE HCO_ERROR_MOD
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC           :: HCO_ERROR
  PUBLIC           :: HCO_WARNING 
  PUBLIC           :: HCO_MSG
  PUBLIC           :: HCO_ENTER
  PUBLIC           :: HCO_LEAVE
  PUBLIC           :: HCO_ERROR_SET
  PUBLIC           :: HCO_ERROR_FINAL
  PUBLIC           :: HCO_VERBOSE_SET
  PUBLIC           :: HCO_VERBOSE_CHECK
  PUBLIC           :: HCO_FORCESCAL_CHECK
  PUBLIC           :: HCO_WILDCARD
  PUBLIC           :: HCO_SEP
  PUBLIC           :: HCO_LOGFILE_OPEN
  PUBLIC           :: HCO_LOGFILE_CLOSE
!
! !DEFINED PARAMETERS:
!
  ! Double and single precision definitions
  INTEGER, PARAMETER, PUBLIC  :: df = kind(0.d0)     ! default 
  INTEGER, PARAMETER, PUBLIC  :: dp = kind(0.d0)     ! double precision
  INTEGER, PARAMETER, PUBLIC  :: sp = 4              ! single precision
  INTEGER, PARAMETER, PUBLIC  :: hp = 8              ! HEMCO precision 

  ! Error success/failure definitions
  INTEGER, PARAMETER, PUBLIC  :: HCO_SUCCESS = 0
  INTEGER, PARAMETER, PUBLIC  :: HCO_FAIL    = -999
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE VARIABLES:
!
  TYPE :: HcoErr
     LOGICAL                     :: Track
     LOGICAL                     :: Verbose
     LOGICAL                     :: LogIsOpen
     LOGICAL                     :: ForceScal 
     LOGICAL                     :: ShowWarnings 
     INTEGER                     :: nWarnings
     INTEGER                     :: CurrLoc
     CHARACTER(LEN=255), POINTER :: Loc(:)
     CHARACTER(LEN=255)          :: LogFile
     CHARACTER(LEN=1)            :: Wildcard
     CHARACTER(LEN=1)            :: Separator
     INTEGER                     :: Lun
  END TYPE HcoErr
!
! !DEFINED PARAMETERS:
!
  ! MAXNEST is the maximum accepted subroutines nesting level.
  ! This only applies to routines with activated error tracking,
  ! i.e. which use HCO_ENTER/HCO_LEAVE statements.
  INTEGER, PARAMETER       :: MAXNEST =  10

  ! Err is the (internal) error type holding information on the
  ! logfile and which part of the code is executed. 
  TYPE(HcoErr), POINTER    :: Err     => NULL()

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_error
!
! !DESCRIPTION: Subroutine HCO_ERROR promts an error message and sets RC to 
! HCO_FAIL. Note that this routine does not stop a run, but it will cause a 
! stop at higher level (when RC gets evaluated). 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ERROR ( ErrMsg, RC, THISLOC )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   )            :: ErrMsg 
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL  :: THISLOC 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)            :: RC 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: I, J
    CHARACTER(LEN=255) :: MSG

    !======================================================================
    ! HCO_ERROR begins here 
    !======================================================================

    ! Print error message
    MSG =  'HEMCO ERROR: ' // TRIM(ErrMsg)
    CALL HCO_MSG ( MSG, SEP1='!' ) 

    ! Print error location
    IF ( PRESENT(THISLOC) ) THEN
       MSG = 'ERROR LOCATION: ' // TRIM( THISLOC )
       CALL HCO_MSG ( MSG )

    ! Traceback
    ELSE
       DO I = 0, Err%CurrLoc-1
          J = Err%CurrLoc-I
          MSG =  'ERROR LOCATION: ' // TRIM( Err%Loc(J) )
          CALL HCO_MSG ( MSG ) 
       ENDDO
    ENDIF

    MSG = ''
    CALL HCO_MSG ( MSG, SEP2='!' ) 

    ! Return w/ error
    RC = HCO_FAIL 

  END SUBROUTINE HCO_ERROR
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_warning
!
! !DESCRIPTION: Subroutine HCO_WARNING promts a warning message without 
! forcing HEMCO to stop, i.e. return code is set to HCO_SUCCESS. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_WARNING ( ErrMsg, RC, THISLOC )
!
! !INPUT PARAMETERS"
!
    CHARACTER(LEN=*), INTENT(IN   )            :: ErrMsg 
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL  :: THISLOC 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)            :: RC 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255) :: MSG

    !======================================================================
    ! HCO_WARNING begins here 
    !======================================================================

    IF ( Err%ShowWarnings ) THEN

       ! Print warning
       MSG = 'HEMCO WARNING: ' // TRIM( ErrMsg )
       CALL HCO_MSG ( MSG ) 

       ! Print location
       IF ( PRESENT(THISLOC) ) THEN
          MSG = '--> LOCATION: ' // TRIM(THISLOC)
          CALL HCO_MSG ( MSG ) 
       ELSEIF ( Err%CurrLoc > 0 ) THEN
          MSG = '--> LOCATION: ' // TRIM(Err%Loc(Err%CurrLoc))
          CALL HCO_MSG ( MSG ) 
       ENDIF
    ENDIF

    ! Return w/ success
    Err%nWarnings = Err%nWarnings + 1
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_WARNING
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_msg
!
! !DESCRIPTION: Subroutine HCO_MSG passes message msg to the HEMCO
! logfile (or to standard output if the logfile is not open).
! Sep1 and Sep2 denote line delimiters before and after the message,
! respectively.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_MSG ( Msg, Sep1, Sep2 )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL  :: Msg
    CHARACTER(LEN=1), INTENT(IN   ), OPTIONAL  :: Sep1
    CHARACTER(LEN=1), INTENT(IN   ), OPTIONAL  :: Sep2
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
    LOGICAL  :: IsOpen
    INTEGER  :: LUN

    !======================================================================
    ! HCO_MSG begins here 
    !======================================================================

    ! Check if Err object is indeed defined
    IF ( .NOT. ASSOCIATED(Err) ) THEN
       IsOpen = .FALSE.
    ELSE
       IsOpen = Err%LogIsOpen
    ENDIF

    ! Use standard output if file not open
    IF ( .NOT. IsOpen ) THEN
       IF ( PRESENT(MSG) ) PRINT *, TRIM(MSG)

    ! Print message to error file 
    ELSE
       LUN = Err%LUN

       IF (LUN > 0 ) THEN
          IF ( PRESENT(SEP1) ) THEN
             WRITE(LUN,'(a)') REPEAT( SEP1, 79) 
          ENDIF
          IF ( PRESENT(MSG) ) THEN
             WRITE(LUN,*) TRIM(MSG)
          ENDIF
          IF ( PRESENT(SEP2) ) THEN
             WRITE(LUN,'(a)') REPEAT( SEP2, 79) 
          ENDIF
       ELSE
          IF ( PRESENT(SEP1) ) THEN
             WRITE(*,'(a)') REPEAT( SEP1, 79) 
          ENDIF
          IF ( PRESENT(MSG) ) THEN
             WRITE(*,*) TRIM(MSG)
          ENDIF
          IF ( PRESENT(SEP2) ) THEN
             WRITE(*,'(a)') REPEAT( SEP2, 79) 
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE HCO_MSG
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_enter
!
! !DESCRIPTION: Subroutine HCO\_LEAVE is called upon entering a routine. 
! It organizes the traceback handling. It is recommended to call this
! routine for 'big' routines but NOT for routines/functions that are 
! frequently called, e.g. inside of loops!
! Note that all subroutines calling HCO_ENTER must also call HCO_LEAVE!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ENTER ( thisLoc, RC )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: thisLoc 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255)  :: Msg, Loc

    !======================================================================
    ! HCO_ENTER begins here 
    !======================================================================

    IF ( .NOT. ASSOCIATED(Err) ) RETURN 

    ! Increment position 
    Err%CurrLoc = Err%CurrLoc + 1
    IF ( Err%CurrLoc > MaxNest ) THEN 
       Msg = 'MaxNest too low, cannot enter ' // TRIM(thisLoc)
       CALL HCO_ERROR ( Msg, RC )
       RETURN 
    ENDIF

    ! Error trap
    IF ( Err%CurrLoc <= 0 ) THEN
       Msg = 'CurrLoc is zero, cannot enter: ' // TRIM(thisLoc)
       CALL HCO_ERROR ( Msg, RC )
       RETURN 
    ENDIF

    ! Register current routine
    Err%Loc(Err%CurrLoc) = thisLoc

    ! Track location if enabled 
    IF ( Err%Track ) THEN
       WRITE(MSG,100) TRIM(thisLoc), Err%CurrLoc
       CALL HCO_MSG( MSG )
    ENDIF

    ! Set RC to success
    RC = HCO_SUCCESS

100 FORMAT( 'HEMCO: Entering ', a, ' (', i2, ')' )  

  END SUBROUTINE HCO_ENTER
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_leave
!
! !DESCRIPTION: Subroutine HCO\_LEAVE is called upon leaving a routine. 
! It organizes the traceback handling. It is recommended to call this
! routine for 'big' routines but NOT for routines/functions that are 
! frequently called, e.g. inside of loops!
! Note that all subroutines calling HCO_LEAVE must also call HCO_ENTER!
!\\
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_LEAVE ( RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255)  :: MSG, LOC

    !======================================================================
    ! HCO_LEAVE begins here 
    !======================================================================

    IF ( .NOT. ASSOCIATED(Err) ) RETURN 
    
    ! Track location if enabled 
    IF ( Err%Track ) THEN
       WRITE(MSG,110) TRIM(Err%Loc(Err%CurrLoc)), Err%CurrLoc
       CALL HCO_MSG( MSG )
    ENDIF

    ! Remove from list 
    Err%Loc(Err%CurrLoc) = ''

    ! Remove current position
    Err%CurrLoc = Err%CurrLoc - 1

    ! Error trap
    IF ( Err%CurrLoc < 0 ) THEN
       Msg = 'CurrLoc is below zero, this should never happen!!' 
       CALL HCO_ERROR ( Msg, RC )
       RETURN 
    ENDIF

    ! Return w/ success 
    RC = HCO_SUCCESS

110 FORMAT( 'HEMCO: Leaving ', a, ' (', i2, ')' )  

  END SUBROUTINE HCO_LEAVE
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_error_set
!
! !DESCRIPTION: Subroutine HCO\_ERROR_SET defines the HEMCO error
! settings. This routine is called at the beginning of a HEMCO
! simulation. Its input parameter are directly taken from the
! HEMCO configuration file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ERROR_SET( LogFile, Verbose, Wildcard, Separator, &
                            ForceScalUnit, ShowWarnings, Track, RC )
!
!  !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)     :: LogFile        ! logfile path+name
    LOGICAL,          INTENT(IN)     :: Verbose        ! run in verbose mode?
    CHARACTER(LEN=1), INTENT(IN)     :: Wildcard       ! wildcard character
    CHARACTER(LEN=1), INTENT(IN)     :: Separator      ! separator character
    LOGICAL,          INTENT(IN)     :: ForceScalUnit  ! allow only unitless scale factors?
    LOGICAL,          INTENT(IN)     :: ShowWarnings   ! prompt warnings?
    LOGICAL,          INTENT(IN)     :: Track          ! track code location?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

    INTEGER :: Lun

    !======================================================================
    ! HCO_ERROR_SET begins here 
    !======================================================================

    ! Nothing to do if already defined
    RC = HCO_SUCCESS
    IF ( ASSOCIATED(Err) ) RETURN 

    ! Allocate error type
    ALLOCATE(Err)
    ALLOCATE(Err%Loc(MAXNEST))
    Err%Loc(:) = ''

    ! Pass values
    Err%LogFile      = TRIM(LogFile)
    Err%Verbose      = Verbose
    Err%ForceScal    = ForceScalUnit 
    Err%Track        = Track
    Err%ShowWarnings = ShowWarnings
    Err%Wildcard     = Wildcard
    Err%Separator    = Separator 

    ! Init misc. values
    Err%LogIsOpen = .FALSE.
    Err%nWarnings = 0
    Err%CurrLoc   = 0

    ! Set lun to -1 (write into default file) or 0 (write into specified
    ! logfile)
    IF ( INDEX(TRIM(Err%LogFile),TRIM(Wildcard)) > 0 ) THEN
       LUN = -1
    ELSE
       LUN = 0
    ENDIF
    Err%Lun = LUN

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_ERROR_SET
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_error_final
!
! !DESCRIPTION: Subroutine HCO\_ERROR_FINAL finalizes the error type.
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ERROR_FINAL 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

    INTEGER :: STAT 

    !======================================================================
    ! HCO_ERROR_FINAL begins here 
    !======================================================================

    ! Eventually close logfile
    CALL HCO_LOGFILE_CLOSE ( ShowSummary=.TRUE. ) 

    IF ( ASSOCIATED(Err) ) THEN
       IF ( ASSOCIATED(Err%Loc) ) DEALLOCATE(Err%Loc)
       DEALLOCATE(Err)
       Err=>NULL()
    ENDIF

  END SUBROUTINE HCO_ERROR_FINAL
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_verbose_set
!
! !DESCRIPTION: Subroutine HCO\_VERBOSE\_SET sets the verbose flag. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_VERBOSE_SET ( isVerbose )
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN) :: isVerbose 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! HCO_VERBOSE_SET begins here 
    !======================================================================

    IF ( ASSOCIATED(Err) ) THEN
       Err%Verbose = isVerbose 
    ENDIF

  END SUBROUTINE HCO_VERBOSE_SET
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_verbose_check
!
! !DESCRIPTION: Function HCO\_VERBOSE\_CHECK returns the verbose flag. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_VERBOSE_CHECK() RESULT( isVerbose )
!
! !RETURN VALUE:
!
    LOGICAL  :: isVerbose
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! HCO_VERBOSE_CHECK begins here 
    !======================================================================

    IF ( ASSOCIATED(Err) ) THEN
       isVerbose = Err%Verbose
    ELSE
       isVerbose = .FALSE.
    ENDIF

  END FUNCTION HCO_VERBOSE_CHECK
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_forcescal_check
!
! !DESCRIPTION: Function HCO\_FORCESCAL\_CHECK returns TRUE if scale 
! factors shall be forced to be unitless, FALSE otherwise. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_FORCESCAL_CHECK() RESULT( ForceScalUnit )
!
! !RETURN VALUE::
!
    LOGICAL  :: ForceScalUnit
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! HCO_FORCESCAL_CHECK begins here 
    !======================================================================

    IF ( ASSOCIATED(Err) ) THEN
       ForceScalUnit = Err%ForceScal
    ELSE
       ForceScalUnit = .TRUE.
    ENDIF

  END FUNCTION HCO_FORCESCAL_CHECK
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_wildcard
!
! !DESCRIPTION: Function HCO\_WILDCARD returns the WILDCARD character. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_WILDCARD() RESULT( WILDCARD ) 
!
! !RETURN VALUE:
!
    CHARACTER(LEN=1) :: WILDCARD 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! HCO_WILDCARD begins here 
    !======================================================================

    IF ( ASSOCIATED(Err) ) THEN
       WILDCARD = Err%WILDCARD
    ELSE
       WILDCARD = '*' 
    ENDIF

  END FUNCTION HCO_WILDCARD
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_SEP
!
! !DESCRIPTION: Function HCO\_SEP returns the separator character. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_SEP() RESULT( SEP ) 
!
! !RETURN VALUE::
!
    CHARACTER(LEN=1) :: SEP 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! HCO_SEP begins here 
    !======================================================================

    IF ( ASSOCIATED(Err) ) THEN
       SEP = Err%Separator
    ELSE
       SEP = '/' 
    ENDIF

  END FUNCTION HCO_SEP
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_logfile_open
!
! !DESCRIPTION: Subroutine HCO\_LOGFILE\_OPEN opens the HEMCO logfile
! (if not yet open). 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_LOGFILE_OPEN ( RC )
!
! !USES:
! 
    USE inquireMod,   ONLY : findFreeLUN
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,   INTENT(INOUT)  :: RC 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255) :: MSG
    INTEGER            :: IOS, LUN, FREELUN
    LOGICAL            :: isopen, exists
    LOGICAL, SAVE      :: FIRST = .TRUE.

    !======================================================================
    ! HCO_LOGFILE_OPEN begins here 
    !======================================================================

    ! Init
    RC = HCO_SUCCESS

    ! Check if object exists
    IF ( .NOT. ASSOCIATED(Err)) THEN
       PRINT *, 'Cannot open logfile - Err object not defined!'
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Don't do anything if we write into standard output!
    IF ( Err%LUN < 0 ) THEN
       Err%LogIsOpen = .TRUE.

    ! Explicit HEMCO logfile:
    ELSE
    
       ! Find free LUN just in case we need it!
       FREELUN = findFreeLun()
   
       ! Inquire if file is already open
       INQUIRE( FILE=TRIM(Err%LogFile), OPENED=isOpen, EXIST=exists, NUMBER=LUN )
   
       ! File exists and is opened ==> nothing to do
       IF ( exists .AND. isOpen ) THEN
          Err%LUN       = LUN
          Err%LogIsOpen = .TRUE.
   
       ! File exists but not opened ==> reopen
       ELSEIF (exists .AND. .NOT. isOpen ) THEN
          OPEN ( UNIT=FREELUN,   FILE=TRIM(Err%LogFile), STATUS='OLD', &
                 ACTION='WRITE', ACCESS='APPEND',     IOSTAT=IOS   )
          IF ( IOS /= 0 ) THEN
             PRINT *, 'Cannot reopen logfile: ' // TRIM(Err%LogFile)
             RC = HCO_FAIL
             RETURN
          ENDIF
          Err%LUN       = FREELUN
          Err%LogIsOpen = .TRUE.

       ! File does not yet exist ==> open new file
       ELSE
          OPEN ( UNIT=FREELUN, FILE=TRIM(Err%LogFile),     & 
               STATUS='NEW',  ACTION='WRITE', IOSTAT=IOS )
          IF ( IOS /= 0 ) THEN
             PRINT *, 'Cannot create logfile: ' // TRIM(Err%LogFile)
             RC = HCO_FAIL
             RETURN
          ENDIF
          Err%LUN       = FREELUN
          Err%LogIsOpen = .TRUE.
       ENDIF
    ENDIF

    ! Write header on first call
    IF ( FIRST ) THEN
       IF ( Err%LUN < 0 ) THEN
          WRITE(*,'(a)') REPEAT( '-', 79) 
          WRITE(*,*    ) 'Using HEMCO v1.0'
          WRITE(*,'(a)') REPEAT( '-', 79) 
       ELSE
          LUN = Err%LUN
          WRITE(LUN,'(a)') REPEAT( '-', 79) 
          WRITE(LUN,*    ) 'Using HEMCO v1.0'
          WRITE(LUN,'(a)') REPEAT( '-', 79) 
       ENDIF

       FIRST = .FALSE.
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_LOGFILE_OPEN
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_logfile_close
!
! !DESCRIPTION: Subroutine HCO\_LOGFILE\_CLOSE closes the HEMCO logfile.
! If argument ShowSummary is enabled, it will prompt a summary of the 
! HEMCO run up to this point (# of warnings, etc.). 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_LOGFILE_CLOSE ( ShowSummary ) 
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN), OPTIONAL :: ShowSummary
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER            :: IOS
    LOGICAL            :: Summary
    CHARACTER(LEN=255) :: MSG

    !======================================================================
    ! HCO_LOGFILE_CLOSE begins here 
    !======================================================================

    ! Check if object exists
    IF ( .NOT. ASSOCIATED(Err)) RETURN
    IF ( .NOT. Err%LogIsOpen  ) RETURN

    ! Show summary?
    IF ( PRESENT(ShowSummary) ) THEN
       Summary = ShowSummary 
    ELSE
       Summary = .FALSE.
    ENDIF

    ! Eventually print summary
    IF ( Summary ) THEN
       MSG = ' '
       CALL HCO_MSG ( MSG )
       MSG = 'HEMCO FINISHED'
       CALL HCO_MSG ( MSG, SEP1='-' )
  
       WRITE(MSG,'(A20,I6)') 'Number of warnings: ', Err%nWarnings 
       CALL HCO_MSG ( MSG, SEP2='-' )
    ENDIF

    ! Close logfile only if lun is defined 
    IF ( Err%Lun>0 ) THEN
       CLOSE ( UNIT=Err%Lun, IOSTAT=IOS )
       IF ( IOS/= 0 ) THEN
          PRINT *, 'Cannot close logfile: ' // TRIM(Err%LogFile)
       ENDIF
    ENDIF
    Err%LogIsOpen = .FALSE.

  END SUBROUTINE HCO_LOGFILE_CLOSE
!EOC
      END MODULE HCO_ERROR_MOD

