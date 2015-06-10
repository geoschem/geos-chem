!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_error_mod.F90
!
! !DESCRIPTION: Module HCO\_Error\_Mod contains routines and variables 
! for error handling and logfile messages in HEMCO. It also contains 
! definitions of some globally used parameter, such as the single/double 
! precision as well as the HEMCO precision definitions. The HEMCO precision 
! is used for almost all HEMCO internal data arrays and can be changed below 
! if required. 
!\\
!\\
! The error settings are specified in the HEMCO configuration file and
! error handling is performed according to these settings. They include: 
!
! \begin{enumerate}
! \item HEMCO logfile: all HEMCO information is written into the specified
!     logfile. The logfile can be set to the wildcard character, in which
!     case the standard output will be used (this may be another opened
!     logfile).
! \item Verbose: Number indicating the verbose level to be used. 
!     0 = no verbose, 3 = very verbose. The verbose level can be set in
!     the HEMCO configuration file. The default value is 0.
! \item Show warnings: if TRUE, prompt all warnings to the HEMCO logfile.
!     Otherwise, no warnings will be prompted but the total number of
!     warnings occurred will still be shown at the end of the run.
! \end{enumerate}
! 
! The error settings are set via subroutine HCO\_ERROR\_SET, called when
! reading section 'settings' of the HEMCO configuration file (subroutine
! Config\_ReadFile in hco\_config\_mod.F90). The currently active verbose
! settings can be checked using subroutines HCO\_IsVerb and
! HCO\_VERBOSE\_INQ. Messages can be written into the logfile using 
! subroutine HCO\_MSG. Note that the logfile actively need to be opened 
! (HCO\_LOGFILE\_OPEN) before writing to it.
!\\
!\\
! The verbose and warning settings are all set to false if it's not the 
! root CPU. 
!\\
!\\
! !INTERFACE: 
!
MODULE HCO_Error_Mod
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
  PUBLIC           :: HCO_IsVerb
  PUBLIC           :: HCO_VERBOSE_INQ
  PUBLIC           :: HCO_LOGFILE_OPEN
  PUBLIC           :: HCO_LOGFILE_CLOSE
!
! !MODULE VARIABLES:
!
  ! Double and single precision definitions
  INTEGER, PARAMETER, PUBLIC  :: dp = kind(0.d0)     ! double precision
  INTEGER, PARAMETER, PUBLIC  :: sp = kind(0.0)      ! single precision
#if defined( USE_REAL8 )
  INTEGER, PARAMETER, PUBLIC  :: hp = kind(0.d0)     ! HEMCO precision r8
#else
  INTEGER, PARAMETER, PUBLIC  :: hp = kind(0.0)      ! HEMCO precision r4
#endif
  INTEGER, PARAMETER, PUBLIC  :: i4 = 4              ! FourByteInt 
  INTEGER, PARAMETER, PUBLIC  :: i8 = 8              ! EightByteInt

  ! Error success/failure definitions
  INTEGER, PARAMETER, PUBLIC  :: HCO_SUCCESS = 0
  INTEGER, PARAMETER, PUBLIC  :: HCO_FAIL    = -999

  ! Tiny value for math operations: 
  ! --> deprecated. Use TINY(1.0_hp) instead!
  REAL(hp), PARAMETER, PUBLIC :: HCO_TINY    = 1.0e-32_hp

  ! Missing value
  ! Note: define missing value as single precision because all data arrays
  ! are read/stored in single precision.
  REAL(sp), PARAMETER, PUBLIC :: HCO_MISSVAL = -1.e-31_sp

  ! Cycle flags. Used to determine the temporal behavior of data
  ! fields. For example, data set to 'cycle' will be recycled if
  ! the simulation date is outside of the datetime range of the
  ! source data. See data reading routine (hcoio_dataread_mod.F90)
  ! for more details. 
  INTEGER, PARAMETER, PUBLIC  :: HCO_CFLAG_CYCLE = 1
  INTEGER, PARAMETER, PUBLIC  :: HCO_CFLAG_RANGE = 2
  INTEGER, PARAMETER, PUBLIC  :: HCO_CFLAG_EXACT = 3
  INTEGER, PARAMETER, PUBLIC  :: HCO_CFLAG_INTER = 4

  ! Data container update flags. At the moment, those only indicate
  ! if a container gets updated every time step or based upon the 
  ! details specifications ('srcTime') in the HEMCO configuration file.
  INTEGER, PARAMETER, PUBLIC  :: HCO_UFLAG_FROMFILE = 1
  INTEGER, PARAMETER, PUBLIC  :: HCO_UFLAG_ALWAYS   = 2

  ! Data container types. These are used to distinguish between
  ! base emissions, scale factors and masks.
  INTEGER, PARAMETER, PUBLIC  :: HCO_DCTTYPE_BASE = 1
  INTEGER, PARAMETER, PUBLIC  :: HCO_DCTTYPE_SCAL = 2
  INTEGER, PARAMETER, PUBLIC  :: HCO_DCTTYPE_MASK = 3

  ! HEMCO version number. Only increase after significant changes
  CHARACTER(LEN=12), PARAMETER, PUBLIC :: HCO_VERSION = 'v1.1.005'
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller   - Initialization
!  12 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  12 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!  03 Mar 2015 - C. Keller   - Added HCO_CFLAG_* and HCO_DCTTYPE_*
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE VARIABLES:
!
  TYPE :: HcoErr
     LOGICAL                     :: IsRoot
     LOGICAL                     :: LogIsOpen
     INTEGER                     :: Warnings 
     INTEGER                     :: Verbose
     INTEGER                     :: nWarnings
     INTEGER                     :: CurrLoc
     CHARACTER(LEN=255), POINTER :: Loc(:)
     CHARACTER(LEN=255)          :: LogFile
     INTEGER                     :: Lun
  END TYPE HcoErr

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
! !IROUTINE: HCO_Error
!
! !DESCRIPTION: Subroutine HCO\_Error promts an error message and sets RC to 
! HCO\_FAIL. Note that this routine does not stop a run, but it will cause a 
! stop at higher level (when RC gets evaluated). 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Error( ErrMsg, RC, THISLOC )
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

  END SUBROUTINE HCO_Error
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Warning
!
! !DESCRIPTION: Subroutine HCO\_Warning promts a warning message without 
! forcing HEMCO to stop, i.e. return code is set to HCO\_SUCCESS. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Warning( ErrMsg, RC, WARNLEV, THISLOC )
!
! !INPUT PARAMETERS"
!
    CHARACTER(LEN=*), INTENT(IN   )            :: ErrMsg 
    INTEGER         , INTENT(IN   ), OPTIONAL  :: WARNLEV 
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL  :: THISLOC 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)            :: RC 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!  26 Mar 2015 - C. Keller - Added warning levels
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER            :: WLEV
    CHARACTER(LEN=255) :: MSG

    !======================================================================
    ! HCO_WARNING begins here 
    !======================================================================

    IF ( PRESENT(WARNLEV) ) THEN
       WLEV = WARNLEV
    ELSE
       WLEV = 3
    ENDIF

    IF ( Err%Warnings >= WLEV ) THEN

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

       ! Increase # of warnings
       Err%nWarnings = Err%nWarnings + 1
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_Warning
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_MSG
!
! !DESCRIPTION: Subroutine HCO\_MSG passes message msg to the HEMCO
! logfile (or to standard output if the logfile is not open).
! Sep1 and Sep2 denote line delimiters before and after the message,
! respectively.
! The optional argument Verb denotes the minimum verbose level associated
! with this message. The message will only be prompted if the verbose level
! on this CPU (e.g. of this Err object) is at least as high as Verb.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_MSG( Msg, Sep1, Sep2, Verb )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL  :: Msg
    CHARACTER(LEN=1), INTENT(IN   ), OPTIONAL  :: Sep1
    CHARACTER(LEN=1), INTENT(IN   ), OPTIONAL  :: Sep2
    INTEGER,          INTENT(IN   ), OPTIONAL  :: Verb
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller   - Initialization
!  20 May 2015 - R. Yantosca - Minor formatting fix: use '(a)' instead of *
!                              to avoid line wrapping around at 80 columns.
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
       
       ! Don't print if this is not the root CPU
       IF ( .NOT. Err%IsRoot ) RETURN

       ! Don't print if verbose level is smaller than verbose level of this
       ! CPU.
       IF ( PRESENT( Verb ) ) THEN
          IF ( Verb < Err%Verbose ) RETURN
       ENDIF
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
!             WRITE(LUN,*) TRIM(MSG)
             WRITE(LUN,'(a)') TRIM(MSG)
          ENDIF
          IF ( PRESENT(SEP2) ) THEN
             WRITE(LUN,'(a)') REPEAT( SEP2, 79) 
          ENDIF
       ELSE
          IF ( PRESENT(SEP1) ) THEN
             WRITE(*,'(a)') REPEAT( SEP1, 79) 
          ENDIF
          IF ( PRESENT(MSG) ) THEN
!             WRITE(*,*) TRIM(MSG)
             WRITE(*,'(a)') TRIM(MSG)
          ENDIF
          IF ( PRESENT(SEP2) ) THEN
             WRITE(*,'(a)') REPEAT( SEP2, 79) 
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE HCO_Msg
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Enter
!
! !DESCRIPTION: Subroutine HCO\_Enter is called upon entering a routine. 
! It organizes the traceback handling. It is recommended to call this
! routine for 'big' routines but NOT for routines/functions that are 
! frequently called, e.g. inside of loops!
!\\
!\\
! Note that all subroutines calling HCO\_Enter must also call HCO\_Leave!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Enter( thisLoc, RC )
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
       CALL HCO_Error( Msg, RC )
       RETURN 
    ENDIF

    ! Error trap
    IF ( Err%CurrLoc <= 0 ) THEN
       Msg = 'CurrLoc is zero, cannot enter: ' // TRIM(thisLoc)
       CALL HCO_Error( Msg, RC )
       RETURN 
    ENDIF

    ! Register current routine
    Err%Loc(Err%CurrLoc) = thisLoc

    ! Track location if enabled 
    IF ( Err%Verbose >= 3 ) THEN
       WRITE(MSG,100) TRIM(thisLoc), Err%CurrLoc
       CALL HCO_Msg( MSG )
    ENDIF

    ! Set RC to success
    RC = HCO_SUCCESS

100 FORMAT( 'HEMCO: Entering ', a, ' (', i2, ')' )  

  END SUBROUTINE HCO_Enter
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Leave
!
! !DESCRIPTION: Subroutine HCO\_Leave is called upon leaving a routine. 
! It organizes the traceback handling. It is recommended to call this
! routine for 'big' routines but NOT for routines/functions that are 
! frequently called, e.g. inside of loops!
!\\
!\\
! Note that all subroutines calling HCO\_Leave must also call HCO\_Enter!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Leave( RC )
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
    CHARACTER(LEN=255)  :: MSG, LOC

    !======================================================================
    ! HCO_LEAVE begins here 
    !======================================================================

    IF ( .NOT. ASSOCIATED(Err) ) RETURN 
    
    ! Track location if enabled 
    IF ( Err%Verbose >= 3 ) THEN
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

  END SUBROUTINE HCO_Leave
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Error_Set
!
! !DESCRIPTION: Subroutine HCO\_Error\_Set defines the HEMCO error
! settings. This routine is called at the beginning of a HEMCO
! simulation. Its input parameter are directly taken from the
! HEMCO configuration file. If LogFile is set to '*' (asterik), 
! all output is directed to the standard output.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_ERROR_SET( am_I_Root, LogFile,      & 
                            Verbose,   WarningLevel, RC )
!
!  !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)     :: am_I_Root      ! Root CPU? 
    CHARACTER(LEN=*), INTENT(IN)     :: LogFile        ! logfile path+name
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: Verbose        ! verbose level 
    INTEGER,          INTENT(INOUT)  :: WarningLevel   ! warning level 
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

    ! Set verbose to -1 if this is not the root CPU. This will disable any
    ! log-file messages 
    IF ( .NOT. am_I_Root ) THEN
       Verbose      = -1
       WarningLevel =  0
    ENDIF

    ! Pass values
    Err%IsRoot       = am_I_Root
    Err%LogFile      = TRIM(LogFile)
    Err%Verbose      = Verbose
    Err%Warnings     = WarningLevel

    ! Init misc. values
    Err%LogIsOpen = .FALSE.
    Err%nWarnings = 0
    Err%CurrLoc   = 0

    ! If Logfile is set to '*', set lun to -1 (--> write into default file). 
    ! Otherwise, set lun to 0 (--> write into specified logfile)
    IF ( TRIM(Err%LogFile) == '*' ) THEN
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
! !IROUTINE: HCO_Error_Final
!
! !DESCRIPTION: Subroutine HCO\_Error\_Final finalizes the error type.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Error_Final 
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
    CALL HCO_Logfile_Close( ShowSummary=.TRUE. ) 

    IF ( ASSOCIATED(Err) ) THEN
       IF ( ASSOCIATED(Err%Loc) ) DEALLOCATE(Err%Loc)
       DEALLOCATE(Err)
    ENDIF
    Err => NULL()

  END SUBROUTINE HCO_Error_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_Verbose_Inq
!
! !DESCRIPTION: Function HCO\_Verbose\_Inq returns the HEMCO verbose number. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_VERBOSE_INQ ( ) RESULT ( VerbNr )
!
! !OUTPUT PARAMETERS:
!  
    INTEGER :: VerbNr 
!
! !REVISION HISTORY:
!  15 Mar 2015 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! HCO_VERBOSE_INQ begins here 
    !======================================================================

    IF ( .NOT. ASSOCIATED(Err) ) THEN
       VerbNr = -1
    ELSE
       VerbNr = Err%Verbose
    ENDIF

  END FUNCTION HCO_VERBOSE_INQ
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_IsVerb
!
! !DESCRIPTION: Function HCO\_IsVerb returns true if the HEMCO verbose number
!  is equal to or larger than the passed number. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION HCO_IsVerb ( VerbNr ) RESULT ( IsVerb )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: VerbNr 
!
! !OUTPUT PARAMETERS:
!  
    LOGICAL             :: IsVerb
!
! !REVISION HISTORY:
!  15 Mar 2015 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! HCO_IsVerb begins here 
    !======================================================================

    IF ( .NOT. ASSOCIATED(Err) ) THEN
       IsVerb = .FALSE. 
    ELSE
       IsVerb = ( Err%Verbose >= VerbNr ) 
    ENDIF

  END FUNCTION HCO_IsVerb
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_LOGFILE_OPEN
!
! !DESCRIPTION: Subroutine HCO\_LOGFILE\_OPEN opens the HEMCO logfile
! (if not yet open). 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_LogFile_Open( RC )
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
!  23 Sep 2013 - C. Keller   - Initialization
!  14 Aug 2014 - R. Yantosca - Add FORM='FORMATTED' to the OPEN statement
!                              so that the HEMCO log will be a text file
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

    ! Never open if we are not on the root CPU
    IF ( .NOT. Err%IsRoot ) RETURN

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

          ! Replace existing file on first call
          IF ( FIRST ) THEN
             OPEN ( UNIT=FREELUN,   FILE=TRIM(Err%LogFile), STATUS='REPLACE', &
                    ACTION='WRITE', FORM='FORMATTED',       IOSTAT=IOS         )
             IF ( IOS /= 0 ) THEN
                PRINT *, 'Cannot create logfile: ' // TRIM(Err%LogFile)
                RC = HCO_FAIL
                RETURN
             ENDIF

          ! Reopen otherwise
          ELSE
             OPEN ( UNIT=FREELUN,   FILE=TRIM(Err%LogFile), STATUS='OLD',     &
                    ACTION='WRITE', ACCESS='APPEND',        FORM='FORMATTED', &
                    IOSTAT=IOS   )
             IF ( IOS /= 0 ) THEN
                PRINT *, 'Cannot reopen logfile: ' // TRIM(Err%LogFile)
                RC = HCO_FAIL
                RETURN
             ENDIF
          ENDIF

          Err%LUN       = FREELUN
          Err%LogIsOpen = .TRUE.

       ! File does not yet exist ==> open new file
       ELSE
          OPEN ( UNIT=FREELUN,    FILE=TRIM(Err%LogFile),      & 
                 STATUS='NEW',    ACTION='WRITE', IOSTAT=IOS,  &
                 FORM='FORMATTED'                             )
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
          WRITE(*,'(A12,A12)') 'Using HEMCO ', HCO_VERSION
          WRITE(*,'(a)') REPEAT( '-', 79) 
       ELSE
          LUN = Err%LUN
          WRITE(LUN,'(a)') REPEAT( '-', 79) 
          WRITE(LUN,'(A12,A12)') 'Using HEMCO ', HCO_VERSION
          WRITE(LUN,'(a)') REPEAT( '-', 79) 
       ENDIF

       FIRST = .FALSE.
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCO_Logfile_Open
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_LogFile_Close
!
! !DESCRIPTION: Subroutine HCO\_LOGFILE\_CLOSE closes the HEMCO logfile.
! If argument ShowSummary is enabled, it will prompt a summary of the 
! HEMCO run up to this point (number of warnings, etc.). 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_LogFile_Close( ShowSummary ) 
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
    IF ( .NOT. Err%IsRoot     ) RETURN

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
       MSG = 'HEMCO ' // TRIM(HCO_VERSION) // ' FINISHED.'
       CALL HCO_MSG ( MSG, SEP1='-' )
 
       WRITE(MSG,'(A16,I1,A12,I6)') &
          'Warnings (level ', Err%Warnings, ' or lower): ', Err%nWarnings
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

  END SUBROUTINE HCO_LogFile_Close
!EOC
END MODULE HCO_Error_Mod

