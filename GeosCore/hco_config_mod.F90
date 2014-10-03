!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_config_mod
!
! !DESCRIPTION: Module HCO\_CONFIG\_MOD contains routines related 
! to the HEMCO configuration file. It reads the content of the 
! configuration file, checks which entires therein are actually used 
! for this simulation run, and adds this data to the list of files to 
! be read.
! \\
! !DETAILS: The configuration file provides all source file information of the
! emission fields and scale factors to be used. It is read at the
! beginning of a simulation run, and the information of all files to be used 
! in this simulation are written into one of the following reading
! lists: OneTimeList, YearList, MonthList, DayList, HourList. The files
! in the respective lists are read once, every year, every month, every
! day, and every hour. The time stamp column of the configuration file
! ('Yr/Mt/Dy/Hr') determines in which list information is written. For
! example, if Mt is set to -999 (e.g. 2000/-999/1/0), this entry is
! added to the month list, and the current simulation month is used.
! 
! Spatially homogenous scale factors ('scalars') can directly be defined in 
! the configuration file (column 'Scalar'). All file information (sourceFile,
! sourceVar, and time stamp) will be ignored in this case. More than one
! scalar can be specified, in which case these are interpreted as
! temporal variations. For 2, 3, 4, or 7 entries, the provided scalars
! are assumed to be day-of-week scale factors. For 6, 8, 12, or 24
! entries, the scalars are interpreted as diurnal scale factors.
! \\
! !INTERFACE: 
!
      MODULE HCO_CONFIG_MOD
!
! !USES:
!
      USE HCO_ERROR_MOD
      USE HCO_EMISLL_MOD,         ONLY : SclMax
      USE HCO_TYPE_MOD,           ONLY : HCO_State

      IMPLICIT NONE

      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      ! Public routines
      PUBLIC :: SetReadList 
      PUBLIC :: Read2Buffer
      PUBLIC :: GetNextLine
!
! !PRIVATE:
!
      PRIVATE :: Base2Buffer
      PRIVATE :: Scal2Buffer
      PRIVATE :: Mask2Buffer
      PRIVATE :: Ext2Buffer
      PRIVATE :: Settings2Buffer
      PRIVATE :: Set_Settings
      PRIVATE :: AddNewLine
      PRIVATE :: ConfigRead_Line
      PRIVATE :: RegisterPrepare
      PRIVATE :: RgstCheck
      PRIVATE :: Calc_Coverage
      PRIVATE :: Add_ScalIDs
      PRIVATE :: Register_Base
      PRIVATE :: Register_Scal
      PRIVATE :: DoCharSplit_R8
      PRIVATE :: DoCharSplit_Int
      PRIVATE :: ReadAndSplit_Line
      PRIVATE :: ExtractTime
      PRIVATE :: CleanBuffer 
!
! !REVISION HISTORY:
!  18 Jun 2013 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE TYPES/ARGUMENTS:
!
      INTERFACE DoCharSplit
         MODULE PROCEDURE DoCharSplit_R8
         MODULE PROCEDURE DoCharSplit_INT
      END INTERFACE
!
! !MODULE TYPES/ARGUMENTS:
!
      ! Container that carries one line of Input file information
      TYPE, PUBLIC :: HcoCfgLin
         INTEGER                   :: DataType
         CHARACTER(LEN= 31)        :: cName
         CHARACTER(LEN=255)        :: srcFile
         CHARACTER(LEN= 31)        :: srcVar
         CHARACTER(LEN= 31)        :: srcTime
         CHARACTER(LEN= 31)        :: srcDim
         CHARACTER(LEN= 31)        :: srcUnit
         CHARACTER(LEN=255)        :: SplitChar
         CHARACTER(LEN= 31)        :: Char1
         INTEGER                   :: Int1
         INTEGER                   :: Int2
         INTEGER                   :: Int3 
         INTEGER                   :: ExtNr 
         INTEGER                   :: Line 
         INTEGER                   :: ScalIDs(SclMax) 
         TYPE(HcoCfgLin), POINTER  :: NextLine
      ENDTYPE HcoCfgLin

      ! Linked list w/ all input file information
      TYPE(HcoCfgLin), POINTER     :: HcoCfgFil => NULL()

      ! Hemco settings (these initial values will be overwritten 
      ! by the values specified in the input file) 
      LOGICAL                      :: verbose = .FALSE.

      ! # of lines in buffer 
      INTEGER              :: LinesInBuffer = 0

      ! Set TRUE for debugging
      LOGICAL              :: DEBUG_ = .false.

      ! Define tab, space and comment character
      CHARACTER(LEN=1)     :: TAB     = ACHAR(9)
      CHARACTER(LEN=1)     :: SPACE   = ' '
      CHARACTER(LEN=1)     :: COMMENT = '#'

      ! nnScals
      INTEGER              :: nnScals = 200

      !----------------------------------------------------------------
      ! MODULE ROUTINES follow below
      !----------------------------------------------------------------

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: SetReadList
!
! !DESCRIPTION: Subroutine SetReadList writes data to the data reading
! lists. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SetReadList ( am_I_Root, HcoState, XRNG, YRNG, RC )
!
! !ARGUMENTS:
!
      LOGICAL,         INTENT(IN)    :: am_I_Root   ! root CPU?
      TYPE(HCO_State), POINTER       :: HcoState    ! HEMCO state
      REAL*8,          INTENT(IN   ) :: XRNG(2)     ! Grid lon range
      REAL*8,          INTENT(IN   ) :: YRNG(2)     ! Grid lat range
      INTEGER,         INTENT(INOUT) :: RC          ! Error stat
!
! !REVISION HISTORY:
!  18 Jun 2013 - C. Keller: Initialization
!  17 Sep 2013 - C. Keller: Now get data from buffer
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      INTEGER              :: RgstScls(nnScals)

      !======================================================================
      ! SetReadList begins here
      !======================================================================

      ! Init
      RgstScls(:)  = -1

      ! Eventually read HEMCO config file into buffer (if not yet done so)
      IF ( LinesInBuffer == 0 ) THEN
         CALL Read2Buffer( am_I_Root, HcoState%ConfigFile, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN
      ENDIF

      IF ( DEBUG_ .and. am_I_Root ) WRITE(*,*) 'File read!'

      ! Define HEMCO settings and write them into HEMCO state
      CALL Set_Settings ( am_I_Root, HcoState, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      IF ( DEBUG_ .and. am_I_Root ) WRITE(*,*) 'Settings set!'
      IF ( HcoState%verbose ) WRITE(*,*) 'HEMCO verbose active!'

      ! Prepare data in buffer 
      CALL RegisterPrepare ( am_I_Root, HcoState, XRNG, YRNG, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      IF ( HcoState%verbose ) WRITE(*,*) 'Cover and Tracer ID set!'

      ! Register base emissions and write out scale factor IDs to be
      ! registered in the next step
      CALL Register_Base ( am_I_Root, HcoState, RgstScls, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Register scale factors (incl. masks)
      CALL Register_Scal ( am_I_Root, HcoState, RgstScls, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Cleanup buffer
      CALL CleanBuffer

      ! Leave w/ success
      RC = HCO_SUCCESS 

      END SUBROUTINE SetReadList 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Read2Buffer
!
! !DESCRIPTION: Subroutine Read2Buffer reads the HEMCO configuration file into
! the internal HcoCfgFil derived type.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Read2Buffer( am_I_Root, ConfigFile, RC )
!
! !USES:
!
      USE inquireMod,         ONLY : findFreeLUN
      USE CHARPAK_MOD,        ONLY : STRREPL
!
! !ARGUMENTS:
!
      LOGICAL,            INTENT(IN)    :: am_I_Root   ! root CPU?
      CHARACTER(LEN=*),   INTENT(IN)    :: ConfigFile  ! Full file name
      INTEGER,            INTENT(INOUT) :: RC          ! Success?
!
! !REVISION HISTORY:
!  17 Sep 2013 - C. Keller: Initialization (update)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      INTEGER              :: IU_HCO, IOS
      LOGICAL              :: AIR, EOF
      CHARACTER(LEN=255)   :: LINE, LOC

      !======================================================================
      ! Read2Buffer begins here
      !======================================================================

      ! Enter
      LOC = 'Read2Buffer (HCO_CONFIG_MOD.F90)' 

      ! Leave here if data already in buffer 
      IF ( LinesInBuffer > 0 ) THEN
         RC = HCO_SUCCESS 
         RETURN
      ENDIF

      ! For convenience only
      AIR = am_I_Root

      ! debugging
      IF ( AIR ) WRITE( 6, * ) 'reading HEMCO configuration file: ', &
                               TRIM(ConfigFile)

      ! Find free LUN
      IU_HCO = findFreeLUN()

      ! Open input file
      OPEN ( IU_HCO, FILE=TRIM( ConfigFile ), STATUS='OLD', IOSTAT=IOS )
      IF ( IOS /= 0 ) THEN
         LINE = 'Error reading ' // TRIM(ConfigFile)
         CALL HCO_ERROR( LINE, LOC, RC ); RETURN 
      ENDIF 

      ! Loop until EOF 
      DO

         ! Read a line from the file, exit if EOF
         LINE = ConfigRead_Line ( IU_HCO, EOF )
         IF ( EOF ) EXIT

         ! Replace tab characters in LINE (if any) w/ spaces
         CALL STRREPL( LINE, TAB, SPACE )

         ! Start reading base emissions if LINE contains the character
         ! string 'BEGIN BASE EMISSIONS'
         IF ( INDEX ( LINE, 'BEGIN SECTION BASE EMISSIONS' ) > 0 ) THEN
            IF ( DEBUG_ ) WRITE( 6, * ) 'reading base emissions' 
            CALL Base2Buffer ( AIR, IU_HCO, EOF, RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
            IF ( EOF ) EXIT

         ! Start reading scale factors if LINE contains the character
         ! string 'BEGIN SCALE FACTORS'
         ELSE IF ( INDEX ( LINE, 'BEGIN SECTION SCALE FACTORS' ) > 0 ) THEN
            IF ( DEBUG_ ) WRITE( 6, * ) 'reading scale factors' 
            CALL Scal2Buffer ( AIR, IU_HCO, EOF, RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
            IF ( EOF ) EXIT

         ! Start reading masks if LINE contains the character
         ! string 'BEGIN MASKS'
         ELSE IF ( INDEX ( LINE, 'BEGIN SECTION MASKS' ) > 0 ) THEN
            IF ( DEBUG_ ) WRITE( 6, * ) 'reading masks' 
            CALL Mask2Buffer ( AIR, IU_HCO, EOF, RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
            IF ( EOF ) EXIT

         ! Read HEMCO settings. Their values are directly evaluated
         ! and saved outside of the HcoCfgFil object.
         ELSE IF ( INDEX ( LINE, 'BEGIN SECTION SETTINGS' ) > 0 ) THEN
            IF ( DEBUG_ ) WRITE( 6, * ) 'reading settings'
            CALL Settings2Buffer( AIR, IU_HCO, EOF, RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
            IF ( EOF ) EXIT

         ! Start reading extension data if LINE contains the character
         ! string 'BEGIN EXTENSIONS'
         ELSE IF ( INDEX ( LINE, 'BEGIN SECTION EXTENSIONS' ) > 0 ) THEN
            IF ( DEBUG_ ) WRITE( 6, * ) 'reading extensions' 
            CALL Ext2Buffer ( AIR, IU_HCO, EOF, RC )
            IF ( RC /= HCO_SUCCESS ) RETURN
            IF ( EOF ) EXIT

         ENDIF
      ENDDO

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Read2Buffer
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: GetNextLine 
!
! !DESCRIPTION: Subroutine GetNextLine points the passed HcoCfgLin pointer
! to the next container of HcoCfgFil. If the pointer is not associated yet, 
! it is pointed to the head of HcoCfgFil. The returned HcoCfgLin pointer 
! becomes nullified if it already contains the last line, and the error 
! status changes to HCO_FAIL.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GetNextLine( ThisLine, RC ) 
!
! !ARGUMENTS:
!
      TYPE(HcoCfgLin),    POINTER       :: ThisLine
      INTEGER,            INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  21 Nov 2013 - C. Keller: Initialization 
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! GetNextLine begins here
      !======================================================================

      ! Point to head of HcoCfgFil if passed HcoCfgLin is not yet defined.
      IF ( .NOT. ASSOCIATED (ThisLine) ) THEN
         ThisLine => HcoCfgFil
      ELSE
         ThisLine => ThisLine%NextLine
      ENDIF

      ! Set return flag 
      IF ( .NOT. ASSOCIATED ( ThisLine ) ) THEN
         RC = HCO_FAIL
      ELSE
         RC = HCO_SUCCESS
      ENDIF

      END SUBROUTINE GetNextLine 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Base2Buffer
!
! !DESCRIPTION: Subroutine Base2Buffer reads the HEMCO base emissions
! information and writes them into the internal HcoCfgFil derived type.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Base2Buffer( am_I_Root, IU_HCO, EOF, RC )
!
! !USES:
!
      USE HCO_EXTS_MOD,   ONLY : AddExt
!
! !ARGUMENTS:
!
      LOGICAL,            INTENT(IN)    :: am_I_Root   ! root CPU?
      INTEGER,            INTENT(IN)    :: IU_HCO
      LOGICAL,            INTENT(INOUT) :: EOF 
      INTEGER,            INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  17 Sep 2013 - C. Keller: Initialization (update)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      INTEGER                   :: STAT
      CHARACTER(LEN= 31)        :: cName
      CHARACTER(LEN=255)        :: srcFile
      CHARACTER(LEN= 31)        :: srcVar
      CHARACTER(LEN= 31)        :: srcTime
      CHARACTER(LEN= 31)        :: srcDim
      CHARACTER(LEN= 31)        :: srcUnit
      CHARACTER(LEN= 31)        :: modSpec
      CHARACTER(LEN=255)        :: scalIDs
      INTEGER                   :: Cat
      INTEGER                   :: Hier
      INTEGER                   :: ExtNr 
      CHARACTER(LEN=255)        :: LOC

      TYPE(HcoCfgLin), POINTER  :: NewLine => NULL()

      !======================================================================
      ! Base2Buffer begins here
      !======================================================================

      ! Enter
      LOC = 'Base2Buffer'

      ! Register HEMCO core as extension Nr. 0.
      CALL AddExt ( 'CORE', 0 )

      ! Repeat until end of base emissions section is found 
      DO

         !-------------------------------------------------------------------- 
         ! Read line and get desired character strings
         !-------------------------------------------------------------------- 
         CALL ReadAndSplit_Line ( am_I_Root, IU_HCO, cName,   2,  &
                                  srcFile,   3,      srcVar,  4,  &
                                  srcTime,   5,      srcDim,  6,  & 
                                  srcUnit,   7,      modSpec, 8,  &
                                  scalIDs,   9,      Cat,    10,  &
                                  Hier,      11,     ExtNr,   1,  &
                                  STAT                             )

         !-------------------------------------------------------------------- 
         ! Error checks
         !-------------------------------------------------------------------- 

         ! Check for end of file            
         IF ( STAT < 0 ) THEN
            EOF = .TRUE.
            EXIT
         ENDIF

         ! Skip this entry if commented line
         IF ( STAT == 1 ) THEN
            CYCLE
         ENDIF

         ! Leave routine here if end of section encountered
         IF ( STAT == 10 ) THEN
            EXIT
         ENDIF
               
         ! Error if not enough entries found 
         IF ( STAT == 100 ) THEN
            CALL HCO_ERROR ( 'STAT == 100', LOC, RC ); RETURN 
         ENDIF

         ! Output status should be 0 if none of the statuses above applies 
         IF ( STAT /= 0 ) THEN
            CALL HCO_ERROR ( 'STAT /= 0', LOC, RC ); RETURN
         ENDIF

         !-------------------------------------------------------------------- 
         ! Add to HcoCfgFil
         !-------------------------------------------------------------------- 
        
         ! Add blank line to HcoCfgFil linked list 
         CALL AddNewLine ( NewLine )

         ! Set values in new entry
         NewLine%DataType  = 1 
         NewLine%cName     = cName
         NewLine%srcFile   = srcFile
         NewLine%srcVar    = srcVar
         NewLine%srcTime   = srcTime
         NewLine%srcDim    = srcDim
         NewLine%srcUnit   = srcUnit
         NewLine%splitChar = scalIDs
         NewLine%Char1     = modSpec
         NewLine%Int1      = Cat
         NewLine%Int2      = Hier
         NewLine%ExtNr     = 0

         ! Free pointer
         NewLine => NULL()

      ENDDO

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Base2Buffer
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Scal2Buffer
!
! !DESCRIPTION: Subroutine Scal2Buffer reads the HEMCO scale factor 
! information and writes them into the internal HcoCfgFil derived type.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Scal2Buffer( am_I_Root, IU_HCO, EOF, RC )
!
! !USES:
!
!
! !ARGUMENTS:
!
      LOGICAL,            INTENT(IN)    :: am_I_Root   ! root CPU?
      INTEGER,            INTENT(IN)    :: IU_HCO
      LOGICAL,            INTENT(INOUT) :: EOF 
      INTEGER,            INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  17 Sep 2013 - C. Keller: Initialization (update)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      INTEGER                   :: STAT
      CHARACTER(LEN= 31)        :: cName
      CHARACTER(LEN=255)        :: srcFile
      CHARACTER(LEN= 31)        :: srcVar
      CHARACTER(LEN= 31)        :: srcTime
      CHARACTER(LEN= 31)        :: srcDim
      CHARACTER(LEN= 31)        :: srcUnit
      CHARACTER(LEN=255)        :: Scalar
      CHARACTER(LEN= 31)        :: notUsedC
      INTEGER                   :: scalID
      INTEGER                   :: Oper
      INTEGER                   :: notUsedI
      CHARACTER(LEN=255)        :: LOC

      TYPE(HcoCfgLin), POINTER  :: NewLine => NULL()

      !======================================================================
      ! Scal2Buffer begins here
      !======================================================================

      ! Enter
      LOC = 'Scal2Buffer'

      ! Repeat until end of scale factor section is found 
      DO

         !-------------------------------------------------------------------- 
         ! Read line and get desired character strings
         !-------------------------------------------------------------------- 
         CALL ReadAndSplit_Line ( am_I_Root, IU_HCO, cName,     2, &
                                  srcFile,   3,      srcVar,    4, &
                                  srcTime,   5,      srcDim,    6, & 
                                  srcUnit,   7,      Scalar,    9, &
                                  notUsedC, -1,      Oper,      8, &
                                  ScalID,    1,      NotUsedI, -1, &
                                  STAT                              )

         !-------------------------------------------------------------------- 
         ! Error checks
         !-------------------------------------------------------------------- 

         ! Check for end of file            
         IF ( STAT < 0 ) THEN
            EOF = .TRUE.
            EXIT
         ENDIF

         ! Skip this entry if commented line
         IF ( STAT == 1 ) THEN
            CYCLE
         ENDIF

         ! Leave routine here if end of section encountered
         IF ( STAT == 10 ) THEN
            EXIT
         ENDIF
               
         ! Error if not enough entries found 
         IF ( STAT == 100 ) THEN
            CALL HCO_ERROR ( 'STAT == 100', LOC, RC ); RETURN 
         ENDIF

         ! Output status should be 0 if none of the statuses above applies 
         IF ( STAT /= 0 ) THEN
            CALL HCO_ERROR ( 'STAT /= 0', LOC, RC ); RETURN
         ENDIF

         !-------------------------------------------------------------------- 
         ! Add to HcoCfgFil
         !-------------------------------------------------------------------- 
        
         ! Add blank line to HcoCfgFil linked list 
         CALL AddNewLine ( NewLine )

         ! Set values in new entry
         NewLine%DataType  = 2 
         NewLine%Int1      = ScalID
         NewLine%cName     = cName
         NewLine%srcFile   = srcFile
         NewLine%srcVar    = srcVar
         NewLine%srcTime   = srcTime
         NewLine%srcDim    = srcDim
         NewLine%srcUnit   = srcUnit
         NewLine%splitChar = Scalar
         NewLine%Int2      = Oper

         ! Free pointer
         NewLine => NULL()

      ENDDO

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Scal2Buffer
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Mask2Buffer
!
! !DESCRIPTION: Subroutine Mask2Buffer reads the HEMCO mask data 
! information and writes them into the internal HcoCfgFil derived type.  
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Mask2Buffer( am_I_Root, IU_HCO, EOF, RC )
!
! !USES:
!
!
! !ARGUMENTS:
!
      LOGICAL,            INTENT(IN)    :: am_I_Root   ! root CPU?
      INTEGER,            INTENT(IN)    :: IU_HCO
      LOGICAL,            INTENT(INOUT) :: EOF 
      INTEGER,            INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  17 Sep 2013 - C. Keller: Initialization (update)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      INTEGER                   :: STAT
      CHARACTER(LEN= 31)        :: cName
      CHARACTER(LEN=255)        :: srcFile
      CHARACTER(LEN= 31)        :: srcVar
      CHARACTER(LEN= 31)        :: srcTime
      CHARACTER(LEN= 31)        :: srcDim
      CHARACTER(LEN= 31)        :: srcUnit
      CHARACTER(LEN= 31)        :: notUsedC
      CHARACTER(LEN=255)        :: MaskBnd
      INTEGER                   :: scalID
      INTEGER                   :: Oper
      INTEGER                   :: NotUsedI
      CHARACTER(LEN=255)        :: LOC

      TYPE(HcoCfgLin), POINTER  :: NewLine => NULL()

      !======================================================================
      ! Mask2Buffer begins here
      !======================================================================

      ! Enter
      LOC = 'Mask2Buffer'

      ! Repeat until end of mask section is found 
      DO

         !-------------------------------------------------------------------- 
         ! Read line and get desired integers and character 
         !-------------------------------------------------------------------- 
         CALL ReadAndSplit_Line ( am_I_Root, IU_HCO, cName,     2, &
                                  srcFile,   3,      srcVar,    4, &
                                  srcTime,   5,      srcDim,    6, & 
                                  srcUnit,   7,      MaskBnd,   9, &
                                  notUsedC,  -1,     ScalID,    1, &
                                  Oper,      8,      NotUsedI, -1, &
                                  STAT                              )

         !-------------------------------------------------------------------- 
         ! Error checks
         !-------------------------------------------------------------------- 

         ! Check for end of file            
         IF ( STAT < 0 ) THEN
            EOF = .TRUE.
            EXIT
         ENDIF

         ! Skip this entry if commented line
         IF ( STAT == 1 ) THEN
            CYCLE
         ENDIF

         ! Leave routine here if end of section encountered
         IF ( STAT == 10 ) THEN
            EXIT
         ENDIF
               
         ! Error if not enough entries found 
         IF ( STAT == 100 ) THEN
            CALL HCO_ERROR ( 'STAT == 100', LOC, RC ); RETURN 
         ENDIF

         ! Output status should be 0 if none of the statuses above applies 
         IF ( STAT /= 0 ) THEN
            CALL HCO_ERROR ( 'STAT /= 0', LOC, RC ); RETURN
         ENDIF

         !-------------------------------------------------------------------- 
         ! Add to HcoCfgFil
         !-------------------------------------------------------------------- 
        
         ! Add blank line to HcoCfgFil linked list 
         CALL AddNewLine ( NewLine )

         ! Set values in new entry
         NewLine%DataType  = 3 
         NewLine%Int1      = ScalID
         NewLine%cName     = cName
         NewLine%srcFile   = srcFile
         NewLine%srcVar    = srcVar
         NewLine%srcTime   = srcTime
         NewLine%srcDim    = srcDim
         NewLine%srcUnit   = srcUnit
         NewLine%splitChar = MaskBnd
         NewLine%Int2      = Oper

         ! Free pointer
         NewLine => NULL()

      ENDDO

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Mask2Buffer
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Ext2Buffer
!
! !DESCRIPTION: Subroutine Ext2Buffer reads the extension switches and
! registers all activated extensions along with the corresponding
! extension number.\\
! In addition, information of all data files linked to any of the extensions 
! is read and stored in the internal derived type HcoCfgFil. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Ext2Buffer( am_I_Root, IU_HCO, EOF, RC )
!
! !USES:
!
      USE HCO_EXTS_MOD,     ONLY : AddExt
      USE CHARPAK_MOD,      ONLY : TRANLC, STRREPL
!
! !ARGUMENTS:
!
      LOGICAL,            INTENT(IN)    :: am_I_Root   ! root CPU?
      INTEGER,            INTENT(IN)    :: IU_HCO
      LOGICAL,            INTENT(INOUT) :: EOF 
      INTEGER,            INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  17 Sep 2013 - C. Keller: Initialization (update)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      INTEGER                   :: STAT
      CHARACTER(LEN= 31)        :: Char1 
      CHARACTER(LEN=255)        :: Char2 
      CHARACTER(LEN= 31)        :: srcVar
      CHARACTER(LEN= 31)        :: srcTime
      CHARACTER(LEN= 31)        :: srcDim
      CHARACTER(LEN= 31)        :: srcUnit
      CHARACTER(LEN= 31)        :: modSpec
      CHARACTER(LEN=255)        :: scalIDs
      INTEGER                   :: Cat
      INTEGER                   :: Hier
      INTEGER                   :: ExtNr 
      CHARACTER(LEN=255)        :: LOC

      LOGICAL                   :: IsDataFile
      CHARACTER(LEN=255)        :: iLine

      TYPE(HcoCfgLin), POINTER  :: NewLine => NULL()

      !======================================================================
      ! Ext2Buffer begins here
      !======================================================================

      ! Enter
      LOC = 'Ext2Buffer'

      !======================================================================= 
      ! The extension section consists of two parts: 
      !
      ! - all extension switches: 
      ! ExtName  on/off  ExtNr
      ! SeaFlux  on      101
      !
      ! - information on data files used by the extensions but managed
      ! through HEMCO. These information is provided the same way as base 
      ! emissions (in fact, extension data is treated in HEMCO as base
      ! emission data with ExtNr different than 0).
      !
      ! In the following, it is first assumed that the line of the config 
      ! file contains data file information. If this assumption fails 
      ! because the config line does not contain enough columns, the same 
      ! line is treated as line containing an extension switch.
      !======================================================================= 

      ! Repeat until end of section or end of file
      DO 

         ! Set data file switch
         IsDataFile = .TRUE.

         ! Try reading data file information from the next line 
         ! of the config file. Write this line out. 
         CALL ReadAndSplit_Line ( am_I_Root, IU_HCO, Char1,    2, &
                                  Char2,     3,      srcVar,   4, &
                                  srcTime,   5,      srcDim,   6, & 
                                  srcUnit,   7,      modSpec,  8, &
                                  scalIDs,   9,      Cat,     10, &
                                  Hier,      11,     ExtNr,    1, &
                                  STAT,      OutLine = iLine       )

         ! If reading of data file information failed because the 
         ! line in the config file did not contain enough columns
         ! (==> STAT = 100), this line is probably a extension switch.
         ! Reread line and save out extension name, on/off toggle and 
         ! extension number. 
         IF ( STAT == 100 ) THEN

            ! Now interpret line read above as extension toggle 
            ! line (ExtName on/off ExtNr). 
            CALL ReadAndSplit_Line ( am_I_Root, IU_HCO, Char1,    1, &
                                     Char2,      2,     srcVar,  -1, &
                                     srcTime,   -1,     srcDim,  -1, &
                                     srcUnit,   -1,     modSpec, -1, &
                                     scalIDs,   -1,     Cat,     -1, &
                                     Hier,      -1,     ExtNr,    3, &
                                     STAT,      InLine = iLine        )
            
            ! Adjust data file switch
            IsDataFile = .FALSE.

         ENDIF

         !-------------------------------------------------------------------- 
         ! Error checks 
         !-------------------------------------------------------------------- 

         ! Leave loop here if end of line encountered
         IF ( STAT < 0 ) THEN
            EOF = .TRUE.
            EXIT
         ENDIF

         ! Jump to next line if commented line encountered
         IF ( STAT == 1 ) THEN
            CYCLE 
         ENDIF

         ! Leave routine here if end of section encountered
         IF ( STAT == 10 ) THEN
            EXIT
         ENDIF

         ! Error if not enough entries found 
         IF ( STAT == 100 ) THEN
            CALL HCO_ERROR ( 'STAT == 100', LOC, RC ); RETURN
         ENDIF

         ! Output status should be 0 if none of the statuses above applies 
         IF ( STAT /= 0 ) THEN
            CALL HCO_ERROR ( 'STAT /= 0', LOC, RC ); RETURN
         ENDIF

         !-------------------------------------------------------------------- 
         ! Register read information 
         !-------------------------------------------------------------------- 

         ! For data files, save all information in a new line of the internal
         ! config file type HcoCfgFil, in analogy to base emissions data.
         IF ( IsDataFile ) THEN

            ! Add blank line to HcoCfgFil linked list 
            CALL AddNewLine ( NewLine )

            ! Set values in new entry
            NewLine%DataType  = 1 
            NewLine%cName     = Char1 
            NewLine%srcFile   = Char2 
            NewLine%srcVar    = srcVar
            NewLine%srcTime   = srcTime
            NewLine%srcDim    = srcDim
            NewLine%srcUnit   = srcUnit
            NewLine%splitChar = scalIDs
            NewLine%Char1     = modSpec
            NewLine%Int1      = Cat
            NewLine%Int2      = Hier
            NewLine%ExtNr     = ExtNr 

            ! Free pointer
            NewLine => NULL()

         ! For extension switches, check if this extension is activated. If 
         ! so, the extension is registered along with its extension number.
         ELSE

            ! The second column is the true/false toggle. So register
            ! extension only if it is used.
            CALL TRANLC(TRIM(Char2))
         
            IF ( DEBUG_ ) THEN
               write(*,*) 'checking extension ', TRIM(Char1)
               write(*,*) 'use? ', TRIM(Char2)
            ENDIF

            IF ( TRIM(Char2) == 'on' ) THEN

               ! Remove colon after extension name 
               CALL STRREPL ( Char1, ':', SPACE )
            
               ! Register extension
               CALL AddExt ( TRIM(Char1), ExtNr )
            ENDIF
         ENDIF
      ENDDO

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Ext2Buffer
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Settings2Buffer
!
! !DESCRIPTION: Subroutine Settings2Buffer reads the HEMCO settings data 
! information, evaluates the values and saved them into internal
! variables (NOT in HcoCfgFil!). 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Settings2Buffer( am_I_Root, IU_HCO, EOF, RC )
!
! !USES:
!
      USE CHARPAK_MOD,        ONLY : STRREPL, STRSPLIT, TRANLC
!
! !ARGUMENTS:
!
      LOGICAL,            INTENT(IN)    :: am_I_Root   ! root CPU?
      INTEGER,            INTENT(IN)    :: IU_HCO
      LOGICAL,            INTENT(INOUT) :: EOF 
      INTEGER,            INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  17 Sep 2013 - C. Keller: Initialization (update)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      INTEGER               :: I, N
      CHARACTER(LEN=255)    :: LINE
      CHARACTER(LEN=255)    :: LOC
      CHARACTER(LEN=255)    :: SUBSTR(255) 

      !======================================================================
      ! Settings2Buffer begins here
      !======================================================================

      ! Enter
      LOC = 'Mask2Buffer'

      ! Do until exit 
      DO 

         ! Read line 
         LINE = ConfigRead_Line ( IU_HCO, EOF )

         ! Return if EOF
         IF ( EOF ) THEN
            RC = HCO_SUCCESS; RETURN 
         ENDIF

         ! Jump to next line if line is commented out
         IF ( LINE(1:1) == COMMENT ) CYCLE

         ! Split character string
         CALL STRREPL ( LINE, TAB,   SPACE     )
         CALL STRSPLIT( LINE, SPACE, SUBSTR, N ) 

         ! Leave here if empty line encountered 
         IF ( N <= 1 ) THEN
            RC = HCO_SUCCESS; RETURN 
         ENDIF

         ! ---------------------------------------------------------------------
         ! Check for verbose switch 
         ! ---------------------------------------------------------------------
         IF ( TRIM(SUBSTR(1)) == 'Verbose:' ) THEN

            ! Write into integer I and set flag in HcoState
            CALL TRANLC( TRIM(SUBSTR(2)) )
            IF ( TRIM(SUBSTR(2)) == 'true' ) THEN
               Verbose = .TRUE.
            ELSE
               Verbose = .FALSE.
            ENDIF
         ENDIF

      ENDDO

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Settings2Buffer
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: AddNewLine 
!
! !DESCRIPTION: Subroutine AddNewLine adds a new (blank) container to
! the HcoCfgFil object. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE AddNewLine ( NewLine ) 
!
! !USES:
!
!
! !ARGUMENTS:
!
      TYPE(HcoCfgLin), POINTER    :: NewLine
!
! !REVISION HISTORY:
!  17 Sep 2013 - C. Keller: Initialization (update)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      TYPE(HcoCfgLin), POINTER    :: BlankLine => NULL()

      !======================================================================
      ! AddNewLine begins here
      !======================================================================

      ! Allocate blank line
      ALLOCATE ( BlankLine ) 

      ! Init values
      BlankLine%DataType   = 0
      BlankLine%cName      = ''
      BlankLine%srcFile    = ''
      BlankLine%srcVar     = ''
      BlankLine%srcTime    = ''
      BlankLine%srcDim     = ''
      BlankLine%srcUnit    = ''
      BlankLine%SplitChar  = ''
      BlankLine%Char1      = ''
      BlankLine%Int1       = -999
      BlankLine%Int2       = -999
      BlankLine%Int3       = -999
      BlankLine%ExtNr      = -999
      BlankLine%ScalIDs(:) = -1

      ! Connect blank line with Input file object
      BlankLine%NextLine => HcoCfgFil
      HcoCfgFil          => BlankLine

      ! Update number of lines in buffer and assign correct line number 
      ! to this container
      LinesInBuffer  = LinesInBuffer + 1
      BlankLine%Line = LinesInBuffer

      ! Output pointer points to the blank line
      NewLine => BlankLine

      END SUBROUTINE AddNewLine 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ConfigRead_Line 
!
! !DESCRIPTION: Subroutine ConfigRead\_Line reads a line from the
! emissions configuration file.  
! Adapted from function READ\_ONE\_LINE in input\_mod.F.
!\\
!\\
! !INTERFACE:
!
      FUNCTION ConfigRead_Line( IU_HCO, EOF ) RESULT( LINE )
!
! !USES:
!
      USE FILE_MOD, ONLY : IOERROR
!
! !INPUT PARAMETERS: 
!
      INTEGER,          INTENT(IN)           :: IU_HCO
!
! !OUTPUT PARAMETERS:
!
      LOGICAL,          INTENT(OUT)          :: EOF   ! Denotes EOF 
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  27 Aug 2010 - R. Yantosca - Added ProTeX headers
!  03 Aug 2012 - R. Yantosca - Now make IU_GEOS a global module variable
!                              so that we can define it with findFreeLun
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: IOS
      CHARACTER(LEN=255) :: LINE, MSG

      !=================================================================
      ! ConfigRead_Line begins here!
      !=================================================================

      ! Initialize
      EOF = .FALSE.

      ! Read a line from the file
      READ( IU_HCO, '(a)', IOSTAT=IOS ) LINE

      ! IO Status < 0: EOF condition
      IF ( IOS < 0 ) THEN
         EOF = .TRUE.
         IF ( DEBUG_ ) WRITE( 6, * ) 'End of file encountered'
         RETURN
      ENDIF

      ! IO Status > 0: true I/O error condition
      IF ( IOS > 0 ) CALL IOERROR( IOS, IU_HCO, 'ConfigRead_Line:1' )

      ! Print the line (if necessary)
      IF ( DEBUG_ ) WRITE( 6, '(a)' ) TRIM( LINE )

      END FUNCTION ConfigRead_Line 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RegisterPrepare
!
! !DESCRIPTION: Subroutine SetCoverAndID extracts the coverages of
! all defined mask fields and the tracer IDs of all base emissions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RegisterPrepare( am_I_Root, HcoState, XRNG, YRNG, RC)
!
! !USES:
!
      USE HCO_Exts_MOD,     ONLY : ExtNrInUse
      USE HCO_Type_Mod,     ONLY : HCO_GetSpecID
!
! !INPUT PARAMETERS: 
!
      LOGICAL,          INTENT(IN    )       :: am_I_Root
      TYPE(HCO_State),  POINTER              :: HcoState   
      REAL*8,           INTENT(IN    )       :: XRNG(2)
      REAL*8,           INTENT(IN    )       :: YRNG(2)
!
! !OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(INOUT)        :: RC
! 
! !REVISION HISTORY: 
!  18 Sep 2013 - C. Keller - Initial version (update) 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE(HcoCfgLin), POINTER     :: ThisLine => NULL()
      INTEGER                      :: ThisCover, ThisID, N, FLAG
      REAL*8                       :: lon1, lon2, lat1, lat2
      INTEGER                      :: SplitInts(SclMax)
      REAL*8                       :: SplitReals(4)
      CHARACTER(LEN=255)           :: MSG, LOC
      LOGICAL                      :: verb

      !=================================================================
      ! RegisterPrepare begins here!
      !=================================================================

      ! Enter
      LOC = 'RegisterPrepare'

      ! Check for verbose flag
      verb = HcoState%verbose .and. am_I_Root

      ! verbose
      IF ( verb ) THEN
         WRITE(6,*) 'Start to prepare fields for registering!'
         WRITE(6,*) 'XRNG: ', XRNG
         WRITE(6,*) 'YRNG: ', YRNG
      ENDIF

      ! Get next (first) line of HcoCfgFil 
      CALL GetNextLine ( ThisLine, FLAG ) 

      ! Loop over all lines
      DO WHILE ( FLAG == HCO_SUCCESS ) 

         ! verbose
         IF ( verb ) WRITE(6,*) 'Prepare ', TRIM(ThisLine%cName) 

         ! For base fields or data fields used in one of the HEMCO
         ! extensions:
         IF ( ThisLine%DataType == 1 ) THEN

            ! Only do for entries that will be used! 
            IF ( ExtNrInUse( ThisLine%ExtNr ) ) THEN

               ! Extract species ID. This will return -1 for undefined
               ! species and 0 for wildcard character.
               ThisID = HCO_GetSpecID ( ThisLine%Char1, HcoState )
               ThisLine%Int3 = ThisID

               ! verbose
               IF ( verb ) WRITE(6,*) 'Assigned ID: ', ThisLine%Int3

               ! Extract scale factors to be applied (if any). Skip for
               ! undefined species
               IF ( ThisID >= 0 ) THEN
                  CALL DoCharSplit ( ThisLine%SplitChar, HcoState, &
                                     SplitInts, N, RC )
                  IF ( RC /= HCO_SUCCESS ) RETURN

                  IF ( N > 0 ) THEN
                     ThisLine%ScalIDs(1:N) = SplitInts(1:N)
                  ENDIF
               ENDIF
            ENDIF

         ! Calculate coverage for masks 
         ELSEIF ( ThisLine%DataType == 3 ) THEN

            ! Extract grid box edges 
            CALL DoCharSplit ( ThisLine%SplitChar, HcoState, &
                               SplitReals, N, RC )
            IF ( RC /= HCO_SUCCESS ) RETURN

            ! Error if less/more than four variables extracted
            IF ( N /= 4 ) THEN
               MSG = 'Cannot properly read mask coverage: ' // &
                  TRIM(ThisLine%cName)
               CALL HCO_ERROR ( MSG, LOC, RC ); RETURN
            ENDIF

            ! Get mask edges
            lon1 = SplitReals(1)
            lat1 = SplitReals(2)
            lon2 = SplitReals(3)
            lat2 = SplitReals(4)

            ! Get Coverage
            ThisCover = CALC_COVERAGE( lon1,    lon2,    &
                                       lat1,    lat2,    &
                                       XRNG(1), XRNG(2), &
                                       YRNG(1), YRNG(2)   )
        
            ! Update information in HcoCfgFil
            ThisLine%Int3 = ThisCover 

            IF ( verb ) THEN
               WRITE(6,*) 'Coverage: ', ThisLine%Int3 
            ENDIF 
         ENDIF

         ! Advance to next line
         CALL GetNextLine ( ThisLine, FLAG ) 
      ENDDO

      ! Cleanup
      ThisLine => NULL()

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE RegisterPrepare 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Set_Settings 
!
! !DESCRIPTION: Subroutine Set\_Settings writes the HEMCO settings into
! State\_HEMCO. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Set_Settings ( am_I_Root, HcoState, RC )
!
! !ARGUMENTS:
!
      LOGICAL,          INTENT(IN)    :: am_I_Root 
      TYPE(HCO_State),  POINTER       :: HcoState 
      INTEGER,          INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  18 Jun 2013 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Verbose?
      HcoState%verbose = verbose

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Set_Settings 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Register_Base
!
! !DESCRIPTION: Subroutine Register\_Base registers all base emission
! data and writes out all associated scale factor IDs. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Register_Base ( am_I_Root, HcoState, RgstScls, RC )
!
! !USES:
!
      USE HCO_Exts_MOD,        ONLY : ExtNrInUse
      USE HCO_READLIST_MOD,    ONLY : ReadList_Set
!
! !ARGUMENTS:
!
      LOGICAL,          INTENT(IN)    :: am_I_Root 
      TYPE(HCO_State),  POINTER       :: HcoState 
      INTEGER,          INTENT(INOUT) :: RgstScls(nnScals) 
      INTEGER,          INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  18 Jun 2013 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      TYPE(HcoCfgLin), POINTER   :: ThisLine => NULL()

      INTEGER               :: Scls(SclMax)
      INTEGER               :: SplitInts(SclMax)
      INTEGER               :: N, ThisID
      INTEGER               :: TIMEVEC(8) 
      INTEGER               :: RgstFlag, FLAG 
      LOGICAL               :: Add, Verb

      !======================================================================
      ! Register_Base begins here
      !======================================================================

      ! Enter 
      Verb = HcoState%verbose .and. am_I_Root

      ! Debugging only
      IF ( Verb ) WRITE(*,*) 'starting registering base emissions'

      ! Point to next (first) line in HcoCfgFil 
      CALL GetNextLine ( ThisLine, FLAG )

      ! Loop over temporary arrays
      DO WHILE ( FLAG == HCO_SUCCESS )  

         ! Skip entry if it's not a base field 
         IF ( (ThisLine%DataType /= 1) ) THEN
            CALL GetNextLine ( ThisLine, FLAG ); CYCLE
         ENDIF

         ! Skip entry if this data is not used by any of the extensions
         IF ( .NOT. ExtNrInUse(ThisLine%ExtNr) ) THEN
            IF ( Verb ) WRITE(6,*) 'Skip because extension turned off!'
            CALL GetNextLine ( ThisLine, FLAG ); CYCLE
         ENDIF 

         ! Debugging
         IF ( Verb ) WRITE(6,*) 'Checking ', TRIM(ThisLine%cName)

         ! Don't register base fields for undefined species, unless it's
         ! data used for one of the HEMCO extensions and the species
         ! flag was set to * (= always read). In the latter case, the
         ! tracer ID search performed in RegisterPrepare returned 0.
         IF ( ThisLine%Int3 < 0 ) THEN
            IF ( Verb ) WRITE(6,*) 'Skip because species not defined!'
            CALL GetNextLine ( ThisLine, FLAG ); CYCLE
         ENDIF

         ! Check if we really need to add this inventory or if there
         ! exists another inventory for this tracer with higher
         ! hierarchy but same category and same extension nr..
         ! Also extract vector of scale factor IDs to be applied to this
         ! base field (Scls).
         ! We can skip this step for extension data fields
         IF ( ThisLine%DataType == 1 ) THEN
            CALL RgstCheck( am_I_Root, HcoState, ThisLine, RgstFlag, RC )
            IF ( RC /= HCO_SUCCESS ) RETURN

            IF ( RgstFlag == 0 ) THEN
               CALL GetNextLine ( ThisLine, FLAG ); CYCLE
            ENDIF
        
            IF ( Verb ) write(*,*) 'RgstFlag = ', RgstFlag
         ELSE
            RgstFlag = 1
         ENDIF

         ! Set Add flag
         IF ( RgstFlag == 2 ) THEN
            Add = .TRUE.
         ELSE
            Add = .FALSE.
         ENDIF

         ! For debugging
         IF ( Verb ) WRITE( 6, * ) 'Registering base field ', &
                                   TRIM(ThisLine%cName)

         ! Extract time stamp to be read
         CALL ExtractTime( ThisLine%srcTime, HcoState, TimeVec, RC ) 
         IF ( RC /= HCO_SUCCESS ) RETURN


         ! Register information in ReadList
         CALL ReadList_Set( am_I_Root  = am_I_Root,               &
                            HcoState   = HcoState,                &
                            DataType   = ThisLine%DataType,       &
                            cName      = TRIM( ThisLine%cName),   &
                            ncFile     = TRIM( ThisLine%srcFile), &
                            ncPara     = TRIM( ThisLine%srcVar),  &
                            ncYrs      = TimeVec(1:2),            &
                            ncMts      = TimeVec(3:4),            &
                            ncDys      = TimeVec(5:6),            &
                            ncHrs      = TimeVec(7:8),            &
                            TrcID      = ThisLine%Int3,           &
                            UseScalIDs = ThisLine%ScalIDs,        &
                            Cat        = ThisLine%Int1,           &
                            Hier       = ThisLine%Int2,           &
                            Oper       = 1,                       &
                            ExtNr      = ThisLine%ExtNr,          &
                            ncRead     = .TRUE.,                  &
                            ESMF_Dim   = ThisLine%srcDim,         &
                            ESMF_Unit  = ThisLine%srcUnit,        &
                            Add        = Add,                     &
                            RC         = RC                        )
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! Collect all scale factor IDs in ScalIDs vector
         CALL Add_ScalIDs ( ThisLine%ScalIDs, RgstScls, RC )
         IF ( RC /= HCO_SUCCESS ) RETURN

         ! Advance to next line
         CALL GetNextLine ( ThisLine, FLAG )

      ENDDO 

      ! Cleanup 
      ThisLine => NULL()

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Register_Base
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Register_Scal
!
! !DESCRIPTION: Subroutine Register\_Scal registers all scale factors.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Register_Scal ( am_I_Root, HcoState, RgstScls, RC )
!
! !USES:
!
      USE HCO_READLIST_MOD,    ONLY : ReadList_Set
!
! !ARGUMENTS:
!
      LOGICAL,          INTENT(IN)    :: am_I_Root      ! am I root? 
      TYPE(HCO_State),  POINTER       :: HcoState 
      INTEGER,          INTENT(IN   ) :: RgstScls(:) 
      INTEGER,          INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  18 Jun 2013 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      TYPE(HcoCfgLin), POINTER  :: ThisLine => NULL()
      INTEGER               :: I, N, FLAG
      INTEGER               :: TimeVec(8)
      CHARACTER(LEN=255)    :: MSG, LOC
      CHARACTER(LEN=  5)    :: strID
      REAL*8, POINTER       :: Array(:,:,:,:) => NULL()
      REAL*8                :: kkVal
      REAL*8                :: ScalFacts(24)
      INTEGER               :: KK
      LOGICAL               :: Verb

      !======================================================================
      ! Register_Scal begins here
      !======================================================================

      ! Enter
      LOC = 'Register_Scal' 
      Verb = HcoState%verbose .and. am_I_Root

      IF ( Verb ) WRITE( 6, * ) 'Registering IDs ', RgstScls

      ! Loop over all scale factor ids
      DO I = 1, SIZE(RgstScls)

         ! RgstScls are filled starting with the first element, so leave
         ! routine as soon as we reach the first non-defined (i.e.
         ! negative) slot
         IF ( RgstScls(I) < 0 ) EXIT

         ! Make ThisLine point to first element of HcoCfgFil
         ThisLine => NULL()
         CALL GetNextLine ( ThisLine, FLAG )         

         ! Loop over all lines in Input file and find the one with the
         ! correct scale factor ID 
         DO WHILE ( FLAG == HCO_SUCCESS )  

            ! Leave if this is the scale field (or mask) with the wanted
            ! scale factor ID. 
            IF ( (ThisLine%DataType == 2) .OR. &
                 (ThisLine%DataType == 3)        ) THEN

                 IF ( ThisLine%Int1 == RgstScls(I) ) EXIT
            ENDIF

            ! Advance to next line otherwise
            CALL GetNextLine ( ThisLine, FLAG )         
         ENDDO      

         ! Return error if scale factor ID not found
         IF ( .NOT. ASSOCIATED(ThisLine) ) THEN
            WRITE ( strID, * ) RgstScls(I) 
            MSG = 'Scale factor not found: ' // strID
            CALL HCO_ERROR ( MSG, LOC, RC); RETURN 
         ENDIF
  
         ! Check for scalar data
         IF ( TRIM(ThisLine%srcFile) == '-' ) THEN 

            IF ( Verb ) WRITE( 6, * ) 'Reading factors from input file!'

            ! Read data into array
            CALL  DoCharSplit ( ThisLine%splitChar, HcoState, &
                                ScalFacts, N, RC ) 
            IF ( RC /= HCO_SUCCESS ) RETURN

            ! Return w/ error if no scale factor defined
            IF ( N == 0 ) THEN
               MSG = 'No scale factors found: ' // TRIM(ThisLine%cName)
               CALL HCO_ERROR ( MSG, LOC, RC); RETURN 
            ENDIF

            ! Copy data into array. Assume all data is temporal
            ! dimension!
            ALLOCATE ( Array(1,1,1,N) )
            Array(1,1,1,1:N) = ScalFacts(1:N)

            ! Set dummy variables for time ranges
            TimeVec(:) = -1
  
            CALL ReadList_Set( am_I_Root = am_I_Root,             &
                               HcoState = HcoState,               &
                               DataType  = ThisLine%DataType,     &
                               cName  = TRIM( ThisLine%cName   ), &
                               ncFile = TRIM( ThisLine%srcFile ), &
                               ncPara = '-',                      &
                               ncYrs  = TimeVec(1:2),             &
                               ncMts  = TimeVec(3:4),             &
                               ncDys  = TimeVec(5:6),             &
                               ncHrs  = TimeVec(7:8),             &
                               ScalID = ThisLine%Int1,            &
                               Oper   = ThisLine%Int2,            &
                               ncRead    = .FALSE.,               &
                               ESMF_Dim  = ThisLine%srcDim,       &
                               ESMF_Unit = ThisLine%srcUnit,      &
                               Add       = .FALSE.,               &
                               RC        = RC,                    &
                               ARRAY     = ARRAY                   )
            IF ( RC /= HCO_SUCCESS ) RETURN

            ! Deallocate array
            IF ( ASSOCIATED ( Array ) ) DEALLOCATE ( Array )

         ! If not scalar data
         ELSE

            ! Extract time stamp to be read
            CALL ExtractTime ( ThisLine%srcTime, HcoState, TimeVec, RC ) 
            IF ( RC /= HCO_SUCCESS ) RETURN

            CALL ReadList_Set( am_I_Root = am_I_Root,                &
                               HcoState  = HcoState,                 &
                               DataType  = ThisLine%DataType,        &
                               cName     = TRIM( ThisLine%cName   ), &
                               ncFile    = TRIM( ThisLine%srcFile ), &
                               ncPara    = TRIM( ThisLine%srcVar  ), &
                               ncYrs     = TimeVec(1:2),             &
                               ncMts     = TimeVec(3:4),             &
                               ncDys     = TimeVec(5:6),             &
                               ncHrs     = TimeVec(7:8),             &
                               ScalID    = ThisLine%Int1,            &
                               Oper      = ThisLine%Int2,            &
                               ncRead    = .TRUE.,                   &
                               ESMF_Dim  = ThisLine%srcDim,          &
                               ESMF_Unit = ThisLine%srcUnit,         &
                               Add       = .FALSE.,                  &
                               RC        = RC                         )
            IF ( RC /= HCO_SUCCESS ) RETURN

         ENDIF

         ! Print some information if verbose mode is on 
         IF ( Verb ) THEN
            WRITE(*,*) 'Scale field registered: ', TRIM(ThisLine%cName)
            WRITE(*,*) 'ncFile                : ', TRIM(ThisLine%srcFile)
            WRITE(*,*) 'ScalID                : ', ThisLine%Int1
            WRITE(*,*) 'DataType              : ', ThisLine%DataType
         ENDIF

      ENDDO !I

      ! Cleanup
      ThisLine => NULL()

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Register_Scal
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtractTime 
!
! !DESCRIPTION: Subroutine ExtractTime extracts the time stamp (ranges)
! from the passed character string. The string is expected to be in
! format Yr/Mt/Dy/Hr. Valid values for Yr, Mt, Dy and Hr are:
! (1) Range of values, separated by - sign: e.g. 2000-2010.
! (2) Single value: 2000
! (3) Asterik as wildcard character: *
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ExtractTime ( CharStr, HcoState, TimeVec, RC ) 
!
! !USES:
!
      USE CHARPAK_MOD,        ONLY : STRSPLIT
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN   )       :: CharStr 
      TYPE(HCO_State),  POINTER             :: HcoState 
!
! !OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(  OUT)       :: TimeVec(8)
      INTEGER,          INTENT(INOUT)       :: RC 
! 
! !REVISION HISTORY: 
!  18 Sep 2013 - C. Keller - Initial version (update) 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      INTEGER               :: I, I0, I1, N
      CHARACTER(LEN=255)    :: SUBSTR(255), DATERNG(255)
      CHARACTER(LEN=255)    :: LOC 

      !=================================================================
      ! ExtractTime begins here!
      !=================================================================

      ! Enter
      LOC = 'ExtractTime'

      ! Init
      TimeVec(:) = -1

      ! Extract strings to be translated into integers 
      CALL STRSPLIT( CharStr, '/', SUBSTR, N )
      IF ( N > 4 ) THEN
         CALL HCO_ERROR( 'Too many substrings!', LOC, RC )
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
         IF ( TRIM(SUBSTR(I)) == TRIM(HcoState%WILDCARD) ) THEN
            TimeVec(I0:I1) = -1 

         ! Otherwise, check if date range if given and set lower and
         ! upper bound accordingly.
         ELSE
            CALL STRSPLIT( SUBSTR(I), '-', DATERNG, N )

            ! If range is given:
            IF ( N == 2 ) THEN
               READ( DATERNG(1), * ) TimeVec(I0)
               READ( DATERNG(2), * ) TimeVec(I1)

            ! Use same value if only one value is given:
            ELSEIF ( N == 1 ) THEN
               READ( DATERNG(1), * ) TimeVec(I0)
               TimeVec(I1) = TimeVec(I0)
            ELSE
               CALL HCO_ERROR( 'Cannot extract year!', LOC, RC )
               RETURN
            ENDIF
         ENDIF
      ENDDO

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE ExtractTime 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DoCharSplit_R8 
!
! !DESCRIPTION: Subroutine DoCharSplit\_R8 splits the passed character
! string into N real8 values. All values are separated by a forward slash
! (/) in the original string.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DoCharSplit_R8 ( CharStr, HcoState, Reals, N, RC ) 
!
! !USES:
!
      USE CHARPAK_MOD,        ONLY : STRSPLIT
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN   )       :: CharStr 
      TYPE(HCO_State),  POINTER             :: HcoState 
!
! !OUTPUT PARAMETERS:
!
      REAL*8,           INTENT(  OUT)       :: Reals(:)
      INTEGER,          INTENT(  OUT)       :: N
      INTEGER,          INTENT(INOUT)       :: RC
! 
! !REVISION HISTORY: 
!  18 Sep 2013 - C. Keller - Initial version (update) 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      INTEGER               :: I
      CHARACTER(LEN=255)    :: SUBSTR(255)
      CHARACTER(LEN=255)    :: LOC

      !=================================================================
      ! DoCharSplit_R8 begins here!
      !=================================================================

      ! Enter
      LOC = 'DoCharSplit_R8'

      ! Init
      Reals(:) = -999d0

      ! Extract strings to be translated into integers 
      CALL STRSPLIT( CharStr, '/', SUBSTR, N )
      IF ( N > SIZE(Reals,1) ) THEN
         CALL HCO_ERROR( 'Too many substrings!', LOC, RC )
         RETURN
      ENDIF   

      ! Return here if no entry found
      IF ( N == 0 ) THEN
         RC = HCO_SUCCESS; RETURN
      ENDIF

      ! Pass all extracted strings to integer vector. Replace wildcard
      ! character with -999!
      DO I = 1, N
         IF ( TRIM(SUBSTR(I)) == TRIM(HcoState%WILDCARD) ) THEN
            Reals(I) = -999d0
         ELSEIF ( TRIM(SUBSTR(I)) == '-' ) THEN
            Reals(I) = -999d0
         ELSE
            READ( SUBSTR(I), * ) Reals(I) 
         ENDIF
      ENDDO

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE DoCharSplit_R8 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DoCharSplit_INT 
!
! !DESCRIPTION: Subroutine DoCharSplit\_INT splits the passed character
! string into N integers. All integers are separated by a forward slash
! (/) in the original string.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DoCharSplit_INT ( CharStr, HcoState, Ints, N, RC ) 
!
! !USES:
!
      USE CHARPAK_MOD,        ONLY : STRSPLIT
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN   )       :: CharStr 
      TYPE(HCO_State),  POINTER             :: HcoState 
!
! !OUTPUT PARAMETERS:
!
      INTEGER,          INTENT(  OUT)       :: Ints(:)
      INTEGER,          INTENT(  OUT)       :: N
      INTEGER,          INTENT(INOUT)       :: RC
! 
! !REVISION HISTORY: 
!  18 Sep 2013 - C. Keller - Initial version (update) 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      INTEGER               :: I
      CHARACTER(LEN=255)    :: SUBSTR(255)
      CHARACTER(LEN=255)    :: LOC

      !=================================================================
      ! DoCharSplit_INT begins here!
      !=================================================================

      ! Enter
      LOC = 'DoCharSplit_Int'

      ! Init
      Ints(:) = -999

      ! Extract strings to be translated into integers 
      CALL STRSPLIT( CharStr, '/', SUBSTR, N )
      IF ( N > SIZE(Ints,1) ) THEN
         CALL HCO_ERROR( 'Too many substrings!', LOC, RC )
         RETURN
      ENDIF   

      ! Return here if no entry found
      IF ( N == 0 ) THEN
         RC = HCO_SUCCESS; RETURN 
      ENDIF

      ! Pass all extracted strings to integer vector. Replace wildcard
      ! character with -999!
      DO I = 1, N
         IF ( TRIM(SUBSTR(I)) == TRIM(HcoState%WILDCARD) ) THEN
            Ints(I) = -999
         ELSEIF ( TRIM(SUBSTR(I)) == '-' ) THEN
            Ints(I) = -999
         ELSE
            READ( SUBSTR(I), * ) Ints(I) 
         ENDIF
      ENDDO

      ! Leave w/ success
      RC = HCO_SUCCESS 

      END SUBROUTINE DoCharSplit_INT 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: RgstCheck 
!
! !DESCRIPTION: Function RgstCheck checks if there is another inventory
! with full coverage for this tracer but with a higher hierarchy, in
! which case we don't have to register the current data. It also checks
! if this is sectorial data that has to be added to data of another
! inventory, e.g. EDGAR sectorial data. 
! The output flags are as follows: 
! 0 = don't need to read the data 
! 1 = read data, don't add
! 2 = read data, add
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RgstCheck ( am_I_Root, HcoState,   &
                             CheckLine, RgstFlag,  RC ) 
!
! !ARGUMENTS:
!
      LOGICAL,         INTENT(IN   )   :: am_I_Root
      TYPE(HCO_State), POINTER         :: HcoState
      TYPE(HcoCfgLin), POINTER         :: CheckLine
      INTEGER,         INTENT(  OUT)   :: RgstFlag 
      INTEGER,         INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      INTEGER                   :: ThisTrc, ThisCat, ThisHier, ThisScal
      INTEGER                   :: ThisExt
      TYPE(HcoCfgLin), POINTER  :: ThisLine => NULL()
      TYPE(HcoCfgLin), POINTER  :: MaskLine => NULL()
      INTEGER                   :: I, J, FLAG1, FLAG2
      LOGICAL                   :: partCov, verb
      CHARACTER(LEN=255)        :: MSG, LOC
      CHARACTER(LEN=  7)        :: strID

      !======================================================================
      ! RgstCheck begins here
      !======================================================================

      ! Enter
      LOC = 'RgstCheck'

      ! Get verbose flag from HCO_State
      verb = HcoState%verbose .AND. am_I_Root

      ! Set output flag to 1 to start with (read without adding)
      RgstFlag = 1 

      ! Get Tracer ID, category and hierarchy of entry to be checked
      ThisExt  = CheckLine%ExtNr
      ThisCat  = CheckLine%Int1
      ThisHier = CheckLine%Int2
      ThisTrc  = CheckLine%Int3

      ! verbose
      IF ( verb ) THEN
         WRITE(6,*) 'RgstCheck: ', TRIM(CheckLine%cName)
         WRITE(6,*) 'ThisExt  : ', ThisExt
         WRITE(6,*) 'ThisCat  : ', ThisCat 
         WRITE(6,*) 'ThisHier : ', ThisHier 
         WRITE(6,*) 'ThisTrc  : ', ThisTrc
      ENDIF

      ! Check all scale factors of the original field to see if one of
      ! them is a mask that has no valid entries over the domain of this
      ! CPU. In this case we don't have to consider this field at all!
      DO I = 1, SclMax

         ! Check if it's a valid scale factor
         IF ( CheckLine%ScalIDs(I) < 0 ) CYCLE

         ! Loop over all lines
         MaskLine => NULL()
         CALL GetNextLine ( MaskLine, FLAG1 ) 
         DO WHILE ( FLAG1 == HCO_SUCCESS ) 

            ! Leave loop if this is the line of interest
            IF ( ( MaskLine%DataType > 1                     ) .AND. & 
                 ( MaskLine%Int1     == CheckLine%ScalIDs(I) ) ) EXIT

            ! Advance to next line otherwise
            CALL GetNextLine ( MaskLine, FLAG1 ) 
         ENDDO

         ! Error if scale factor not found
         IF ( .NOT. ASSOCIATED ( MaskLine ) ) THEN
            WRITE ( strID, * ) CheckLine%ScalIDs(I) 
            MSG = 'Scale factor not found: ' // TRIM(strID) 
            CALL HCO_ERROR ( MSG, LOC, RC); RETURN
         ENDIF

         ! Check if this is a mask with zero coverage over this CPU, in
         ! which case we don't need to consider the base field at all! 
         IF ( MaskLine%DataType == 3 .AND. MaskLine%Int3 == 0 ) THEN
            RgstFlag = 0  
            IF ( am_I_Root .and. verbose ) THEN 
               WRITE(6,*) 'Skip field for this CPU: ', &
                          TRIM(CheckLine%cName)
            ENDIF

            ! Return
            RC = HCO_SUCCESS
            RETURN 
         ENDIF
      ENDDO ! I 

      ! Now find out if there is another base field in the input file  
      ! for the same tracer and with same emission category (and
      ! extension number), but higher hierarchy and full coverage over 
      ! this CPU (in which case we can ignore the base field that is 
      ! checked). 

      ! Initialize looping pointer
      ThisLine => NULL()
      CALL GetNextLine ( ThisLine, FLAG1 ) 

      ! Loop over Inputfile
      DO WHILE ( FLAG1 == HCO_SUCCESS )
        
         ! Advance to next line if this is the line to be checked 
         IF ( CheckLine%Line == ThisLine%Line ) THEN
            CALL GetNextLine ( ThisLine, FLAG1 ); CYCLE
         ENDIF 

         ! Advance to next line if this is not a base field
         IF ( ThisLine%DataType /= 1 ) THEN
            CALL GetNextLine ( ThisLine, FLAG1 ); CYCLE
         ENDIF 

         ! Advance to next line if this is not the same extension nr
         IF ( ThisLine%ExtNr /= ThisExt ) THEN
            CALL GetNextLine ( ThisLine, FLAG1 ); CYCLE
         ENDIF 

         ! Advance to next line if this is not the same tracer
         IF ( ThisLine%Int3 /= ThisTrc ) THEN
            CALL GetNextLine ( ThisLine, FLAG1 ); CYCLE
         ENDIF 

         ! Advance to next line if this is not the same category
         IF ( ThisLine%Int1 /= ThisCat ) THEN
            CALL GetNextLine ( ThisLine, FLAG1 ); CYCLE
         ENDIF 

         ! Advance to next line if this is lower hierarchy
         IF ( ThisLine%Int2 < ThisHier ) THEN
            CALL GetNextLine ( ThisLine, FLAG1 ); CYCLE
         ENDIF 

         ! Check for partial coverage. In this case, we cannot say for 
         ! sure whether or not the field to be checked will be entirely 
         ! overwritten by this new field.
         partCov = .FALSE.

         ! Check all scale factors of both base field to check if:
         ! (1) The original base field has a mask that does not cover
         ! the current domain, in which case we don't need to read this
         ! base field anyway
         DO I = 1, SclMax

            ! Check if it's a valid scale factor
            IF ( ThisLine%ScalIDs(I) < 0 ) CYCLE

            ! Find line with wanted scale ID 
            MaskLine => NULL()
            CALL GetNextLine ( MaskLine, FLAG2 )  
            DO WHILE ( FLAG2 == HCO_SUCCESS ) 

               ! Leave loop if this is the line of interest
               IF ( ( MaskLine%DataType > 1                    ) .AND. &
                    ( MaskLine%Int1     == ThisLine%ScalIDs(I) ) ) EXIT

               ! Advance to next line otherwise
               CALL GetNextLine ( MaskLine, FLAG2 )  
            ENDDO

            ! Error if scale factor not found
            IF ( .NOT. ASSOCIATED ( MaskLine ) ) THEN
            MSG = 'Scale factor error: ' // TRIM(ThisLine%cName)
            CALL HCO_ERROR ( MSG, LOC, RC); RETURN
            ENDIF

            ! Check if this is a mask with only partial overlap.  
            IF ( MaskLine%DataType == 3 .AND. MaskLine%Int3 < 1 ) THEN
               partCov = .TRUE.
               EXIT 
            ENDIF
         ENDDO ! I 

         ! If we made it up to here and partCov is still FALSE, then
         ! the found field has the same tracer ID, the same category, 
         ! and a higher (or the same) hierarchy as the ckecked field I.

         ! If hierarchy of found field is higher than checked field and 
         ! if it has total coverage over this CPU, it will always replace 
         ! all values of the checked field, hence set RgstFlag to 0 
         ! (=don't read field I) and return here.
         IF ( (ThisLine%Int2 > ThisHier) .AND. (.NOT. partCov) ) THEN

            IF ( verb ) THEN
               WRITE(6,*) 'Skip field ',  TRIM(CheckLine%cName), &
                       ' because of ', TRIM(ThisLine%cName)

               RgstFlag = 0
            ENDIF
            
            ! Return
            RC = HCO_SUCCESS 
            RETURN
         ENDIF

         ! If the found and checked field have same hierarchy, same
         ! field name, same scale factors, and same update frequencies, 
         ! we may add the two fields together.
         ! Don't set the ADD flag
         ! to true to the field with the lowest line number! This will
         ! be the field who is stored first in the reading list (and
         ! hence will be read first). 
         ! Don't return here, because it is still possible that there is
         ! another inventory in the list coming up which overwrites this
         ! inventory.
        
         IF ( (TRIM(CheckLine%cName) == TRIM(ThisLine%cName)) .AND. & 
              (ThisHier              == ThisLine%Int2       ) ) THEN 

            ! If the two fields have same container name, they must also
            ! have the same scale factors and time update frequencies!!
            IF ( TRIM(CheckLine%SplitChar) /= &
                 TRIM(ThisLine%SplitChar ) ) THEN
                 MSG = 'Fields must have same scale factors: ' // &
                       TRIM(CheckLine%cName) // ' and ' //        &
                       TRIM(ThisLine%cName)
               CALL HCO_ERROR ( MSG, LOC, RC); RETURN
            ENDIF
            IF ( TRIM(CheckLine%srcTime) /= &
                 TRIM(ThisLine%srcTime ) ) THEN
                 MSG = 'Fields must have same update freqs: ' // &
                       TRIM(CheckLine%cName) // ' and ' //       &
                       TRIM(ThisLine%cName)
               CALL HCO_ERROR ( MSG, LOC, RC); RETURN 
            ENDIF

            IF ( CheckLine%Line > ThisLine%Line ) THEN
               RgstFlag = 2
            ELSE
               RgstFlag = 1
            ENDIF
         ENDIF

         ! Advance to next line
         CALL GetNextLine ( ThisLine, FLAG1 )  

      ENDDO !Loop over all entries in HcoCfgFil 

      ! Free pointers 
      ThisLine => NULL()
      MaskLine => NULL()

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE RgstCheck
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: CALC_COVERAGE
!
! !DESCRIPTION: Function CALC_COVERAGE calculates the coverage of
! the specified lon/lat box with the area covered by the inventory.
! Returns 0 if no overlap, 1 if complete overlap, and -1 for partial
! overlap.
!\\
!\\
! !INTERFACE:
!
      FUNCTION CALC_COVERAGE ( INVI1, INVI2, INVJ1, INVJ2,  &
                               SIMI1, SIMI2, SIMJ1, SIMJ2 ) &
      RESULT ( COVERAGE ) 
!
! !ARGUMENTS:
!
      REAL*8,           INTENT(IN   )  :: INVI1
      REAL*8,           INTENT(IN   )  :: INVI2
      REAL*8,           INTENT(IN   )  :: INVJ1
      REAL*8,           INTENT(IN   )  :: INVJ2
      REAL*8,           INTENT(IN   )  :: SIMI1
      REAL*8,           INTENT(IN   )  :: SIMI2
      REAL*8,           INTENT(IN   )  :: SIMJ1
      REAL*8,           INTENT(IN   )  :: SIMJ2
      INTEGER                          :: COVERAGE
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!

      !======================================================================
      ! CALC_COVERAGE begins here
      !======================================================================

      ! Check if specified area does not overlap with inventory
      COVERAGE = 1
      IF ( INVI1 > SIMI2 .OR. INVI2 < SIMI1 .OR. &
           INVJ1 > SIMJ2 .OR. INVJ2 < SIMJ1        ) THEN
         COVERAGE = 0
         RETURN
      ENDIF

      ! Check for partial coverage
      IF ( INVI1 > SIMI1 .OR. INVI2 < SIMI2 .OR. &
           INVJ1 > SIMJ1 .OR. INVJ2 < SIMJ2        ) THEN
         COVERAGE = -1
      ENDIF

      END FUNCTION CALC_COVERAGE
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Add_ScalIDs
!
! !DESCRIPTION: Subroutine Add\_ScalIDs adds the scale factor IDs to the
! ScalIDs vector. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Add_ScalIDs ( UseScalIDs, RgstScls, RC ) 
!
! !USES:
!
      USE HCO_EMISLL_MOD, ONLY  : SclMax
!
! !ARGUMENTS:
!
      INTEGER,         INTENT(IN   )  :: UseScalIDs(SclMax)
      INTEGER,         INTENT(INOUT)  :: RgstScls(nnScals) 
      INTEGER,         INTENT(INOUT)  :: RC
!
!
! !REVISION HISTORY:
!  27 Mar 2013 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER             :: I, J
      CHARACTER(LEN=255)  :: LOC

      ! ================================================================
      ! Add_ScalIDs begins here
      ! ================================================================

      ! Enter
      LOC = 'Add_ScalIDs (hco_config_mod.F90)'

      ! Loop over all entries in UseScalIDs
      DO J = 1, SclMax
 
         ! Skip this entry in UseScalIDs if invalid
         IF ( UseScalIDs(J) < 0 ) CYCLE

         ! Loop over all RgstScls
         DO I = 1, SIZE(RgstScls)
            
            ! Jump to next ID in UseScalIDs if this ID is already listed
            ! in the RgstScls vector.
            IF ( RgstScls(I) == UseScalIDs(J) ) EXIT
            
            ! The ScalIDs vector is filled starting with the first empty
            ! slot, so as soon as we encounter an empty slot (i.e. -1)
            ! we have passed all currently defined scale IDs and can
            ! thus add this ID to ScalIDs.
            IF ( RgstScls(I) < 0 ) THEN
               RgstScls(I) = UseScalIDs(J)
               EXIT
            ENDIF

            ! Make sure that ScalIDs is long enough 
            IF ( I == nnScals ) THEN
               CALL HCO_ERROR ( 'Increase nnScals!', LOC, RC )
               RETURN 
            ENDIF
         ENDDO !I
      ENDDO !J

      ! Leave w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE Add_ScalIDs 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ReadAndSplit_Line
!
! !DESCRIPTION: Subroutine ReadAndSplit\_Line reads a line from the HEMCO 
! config file and parses the specified columns into the passed integer
! and character variables.
! If the optional argument inLine is provided, this line will be parsed,
! otherwise a new line will be read from the config file. 
! If the optional argument outLine is provided, this variable will hold 
! the parsed line. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE ReadAndSplit_Line ( AIR,     IU_HCO,  char1, chr1cl, &
                                     char2,   chr2cl,  char3, chr3cl, &
                                     char4,   chr4cl,  char5, chr5cl, &
                                     char6,   chr6cl,  char7, chr7cl, &
                                     char8,   chr8cl,  int1,  int1cl, &
                                     int2,    int2cl,  int3,  int3cl, &
                                     STAT,    inLine,  outLine         )
!
! !USES:
!
      USE CHARPAK_MOD,        ONLY : STRREPL, STRSPLIT
!
! !ARGUMENTS:
!
      LOGICAL,            INTENT(IN)              :: AIR  
      INTEGER,            INTENT(IN)              :: IU_HCO     
      CHARACTER(LEN=*),   INTENT(  OUT)           :: char1
      INTEGER,            INTENT(IN   )           :: chr1cl
      CHARACTER(LEN=*),   INTENT(  OUT)           :: char2
      INTEGER,            INTENT(IN   )           :: chr2cl
      CHARACTER(LEN=*),   INTENT(  OUT)           :: char3
      INTEGER,            INTENT(IN   )           :: chr3cl
      CHARACTER(LEN=*),   INTENT(  OUT)           :: char4
      INTEGER,            INTENT(IN   )           :: chr4cl
      CHARACTER(LEN=*),   INTENT(  OUT)           :: char5
      INTEGER,            INTENT(IN   )           :: chr5cl
      CHARACTER(LEN=*),   INTENT(  OUT)           :: char6
      INTEGER,            INTENT(IN   )           :: chr6cl
      CHARACTER(LEN=*),   INTENT(  OUT)           :: char7
      INTEGER,            INTENT(IN   )           :: chr7cl
      CHARACTER(LEN=*),   INTENT(  OUT)           :: char8
      INTEGER,            INTENT(IN   )           :: chr8cl
      INTEGER,            INTENT(  OUT)           :: int1
      INTEGER,            INTENT(IN   )           :: int1cl
      INTEGER,            INTENT(  OUT)           :: int2
      INTEGER,            INTENT(IN   )           :: int2cl
      INTEGER,            INTENT(  OUT)           :: int3
      INTEGER,            INTENT(IN   )           :: int3cl
      INTEGER,            INTENT(INOUT)           :: STAT
      CHARACTER(LEN=255), INTENT(IN   ), OPTIONAL :: inLINE
      CHARACTER(LEN=255), INTENT(  OUT), OPTIONAL :: outLINE
!
! !REVISION HISTORY:
!  28 Aug 2013 - C. Keller: Initial version
!  11 Dec 2013 - C. Keller: Added optional arguments inLine and outLine 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
      INTEGER               :: N
      CHARACTER(LEN=255)    :: LINE
      CHARACTER(LEN=255)    :: LOC
      CHARACTER(LEN=255)    :: SUBSTR(255) 
      LOGICAL               :: EOF

      !======================================================================
      ! ReadAndSplit_Line begins here 
      !======================================================================

      ! Output status
      STAT = 0

      ! ---------------------------------------------------------------------
      ! Read line and split column 
      ! ---------------------------------------------------------------------
      IF ( PRESENT(inLINE) ) THEN
         LINE = inLINE
      ELSE
         LINE = ConfigRead_Line ( IU_HCO, EOF )
         IF ( EOF ) THEN
            STAT = -999
            RETURN
         ENDIF
      ENDIF

      ! Check for output line
      IF ( PRESENT(outLINE) ) outLINE = LINE

      ! Return here with flag = 10 if line starts with 'END SECTION'. 
      IF ( INDEX ( LINE, 'END SECTION' ) > 0 ) THEN
         STAT = 10
         RETURN
      ENDIF

      ! Return here with flag = 1 if line is commented 
      IF ( LINE(1:1) == COMMENT ) THEN
         STAT = 1
         RETURN
      ENDIF

      ! Split line into columns
      CALL STRREPL ( LINE, TAB,   SPACE     )
      CALL STRSPLIT( LINE, SPACE, SUBSTR, N ) 

      ! Also ignore empty lines 
      IF ( N <= 1 ) THEN
         STAT = 1
         RETURN
      ENDIF

      ! ---------------------------------------------------------------------
      ! Read characters as specified and write them into given variables 
      ! ---------------------------------------------------------------------

      ! Character 1 
      IF ( chr1cl > 0 ) THEN
         IF ( chr1cl > N ) THEN
            STAT = 100
            RETURN
         ELSE
            READ( SUBSTR(chr1cl), '(a)' ) char1
         ENDIF
      ENDIF 

      ! Character 2 
      IF ( chr2cl > 0 ) THEN
         IF ( chr2cl > N ) THEN
            STAT = 100
            RETURN
         ELSE
            READ( SUBSTR(chr2cl), '(a)' ) char2
         ENDIF
      ENDIF 

      ! Character 3 
      IF ( chr3cl > 0 ) THEN
         IF ( chr3cl > N ) THEN
            STAT = 100
            RETURN
         ELSE
            READ( SUBSTR(chr3cl), '(a)' ) char3
         ENDIF
      ENDIF 

      ! Character 4 
      IF ( chr4cl > 0 ) THEN
         IF ( chr4cl > N ) THEN
            STAT = 100
            RETURN
         ELSE
            READ( SUBSTR(chr4cl), '(a)' ) char4
         ENDIF
      ENDIF 

      ! Character 5 
      IF ( chr5cl > 0 ) THEN
         IF ( chr5cl > N ) THEN
            STAT = 100
            RETURN
         ELSE
            READ( SUBSTR(chr5cl), '(a)' ) char5
         ENDIF
      ENDIF 

      ! Character 6 
      IF ( chr6cl > 0 ) THEN
         IF ( chr6cl > N ) THEN
            STAT = 100
            RETURN
         ELSE
            READ( SUBSTR(chr6cl), '(a)' ) char6
         ENDIF
      ENDIF 

      ! Character 7 
      IF ( chr7cl > 0 ) THEN
         IF ( chr7cl > N ) THEN
            STAT = 100
            RETURN
         ELSE
            READ( SUBSTR(chr7cl), '(a)' ) char7
         ENDIF
      ENDIF 

      ! Character 8 
      IF ( chr8cl > 0 ) THEN
         IF ( chr8cl > N ) THEN
            STAT = 100
            RETURN
         ELSE
            READ( SUBSTR(chr8cl), '(a)' ) char8
         ENDIF
      ENDIF 

      ! ---------------------------------------------------------------------
      ! Read integers as specified and write them into given variables 
      ! ---------------------------------------------------------------------

      ! Integer 1 
      IF ( int1cl > 0 ) THEN
         IF ( int1cl > N ) THEN
            STAT = 100
            RETURN
         ELSE
            READ( SUBSTR(int1cl), * ) int1 
         ENDIF
      ENDIF 

      ! Integer 2 
      IF ( int2cl > 0 ) THEN
         IF ( int2cl > N ) THEN
            STAT = 100
            RETURN
         ELSE
            READ( SUBSTR(int2cl), * ) int2
         ENDIF
      ENDIF 

      ! Integer 3 
      IF ( int3cl > 0 ) THEN
         IF ( int3cl > N ) THEN
            STAT = 100
            RETURN
         ELSE
            READ( SUBSTR(int3cl), * ) int3
         ENDIF
      ENDIF

      END SUBROUTINE ReadAndSplit_Line
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CleanBuffer 
!
! !DESCRIPTION: Subroutine CleanBuffer deletes all HEMCO input file
! entries from buffer. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CleanBuffer 
!
! !INPUT PARAMETERS:
!
!
! !REVISION HISTORY: 
!  18 Sep 2013 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      TYPE(HcoCfgLin), POINTER  :: ThisLine => NULL()
      TYPE(HcoCfgLin), POINTER  :: NextLine => NULL()

      ! Loop over all lines
      ThisLine => HcoCfgFil
      DO WHILE ( ASSOCIATED ( ThisLine ) ) 

         ! Archive NextLine
         NextLine => ThisLine%NextLine

         ! Remove this line
         ThisLine%NextLine => NULL()
         DEALLOCATE ( ThisLine )

         ! Move to archived next line
         ThisLine => NextLine
      ENDDO

      ! Cleanup
      ThisLine => NULL()
      NextLine => NULL()

      END SUBROUTINE CleanBuffer 
!EOC
      END MODULE HCO_CONFIG_MOD
!EOF
