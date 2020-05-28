!BOP
!
! !MODULE: hco_config_mod.F90
!
! !DESCRIPTION: Module HCO\_Config\_Mod contains routines related
! to the HEMCO configuration file. It reads the content of the
! configuration file, checks which entires therein are actually used
! for this simulation run, and stores these information. This occurs
! in two calls: Config\_ReadFile and SetReadList. Config\_ReadFile
! writes the entire content of the configuration file into buffer
! except for the input data associated with a disabled extension.
! SetReadList does many more logical checks and adds all data used
! by HEMCO to ReadList. Scale factors not used by any of the base
! emissions and base emission fields (e.g. scale factors that won't
! be used ) are removed in this step.
!\\
!\\
! All data fields are saved in individual data containers, which are
! organized in the ConfigList. Hence, ConfigList is a collection of
! all HEMCO data containers, with every container representing an
! entry of the configuration file. Each data container has its unique
! container ID for identification. All HEMCO lists (ConfigList,
! ReadList, EmisList) access the same containers.
!\\
!\\
! The configuration file provides all source file information of the
! emission fields and scale factors to be used. It must be read at the
! beginning of a simulation run.
!\\
!\\
! As of HEMCO v2.0, the ConfigList linked list sits within the HEMCO
! configuration object (HcoConfig). HcoConfig must be passed to all
! routines. This allows the parallel usage of multiple invocations
! of HEMCO that use different input data. HcoConfig is initialized
! upon reading the HEMCO configuration file (within subroutine
! Config\_ReadFile).
!\\
!\\
! !INTERFACE:
!
MODULE HCO_Config_Mod
!
! !USES:
!
  USE HCO_ERROR_MOD
  USE HCO_DIAGN_MOD
  USE HCO_CHARTOOLS_MOD
  USE HCO_TYPES_MOD
  USE HCO_STATE_MOD,          ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: SetReadList
  PUBLIC  :: Config_ReadFile
  PUBLIC  :: Config_GetnSpecies
  PUBLIC  :: Config_GetSpecNames
  PUBLIC  :: ConfigInit
!
! !PRIVATE:
!
  PRIVATE :: ReadSettings
  PRIVATE :: ExtSwitch2Buffer
  PRIVATE :: ConfigList_AddCont
  PRIVATE :: Config_ReadCont
  PRIVATE :: RegisterPrepare
  PRIVATE :: Get_targetID
  PRIVATE :: Calc_Coverage
  PRIVATE :: Register_Base
  PRIVATE :: Register_Scal
  PRIVATE :: ReadAndSplit_Line
  PRIVATE :: Config_GetSpecAttr
  PRIVATE :: BracketCheck
  PRIVATE :: AddZeroScal
  PRIVATE :: AddShadowFields
  PRIVATE :: ParseEmisL
  PRIVATE :: CheckForDuplicateName
  PRIVATE :: Hco_GetTagInfo
!
! !REVISION HISTORY:
!  18 Jun 2013 - C. Keller   -  Initialization
!  08 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  08 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  15 Feb 2015 - C. Keller   - Added BracketCheck, AddZeroScal, AddShadowFields
!  15 Feb 2016 - C. Keller   - Update to v2.0: ConfigList now sits in HcoConfig
!  23 Oct 2018 - M. Sulprizio- Make routine ConfigInit public to allow for
!                              initialization of HcoConfig%ModelSpc from the
!                              external model. Also add routine Hco_GetTagInfo
!                              to get information for wildcard strings (e.g.
!                              ?ALL?) used in HEMCO_Config.rc.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE TYPES/ARGUMENTS:
!

  !----------------------------------------------------------------
  ! MODULE ROUTINES follow below
  !----------------------------------------------------------------

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Config_Readfile
!
! !DESCRIPTION: Subroutine CONFIG\_READFILE reads the HEMCO configuration file,
! archives all HEMCO options and settings (including traceback/error setup),
! and creates a data container for every (used) emission field in the config.
! file. All containers become linked through the ConfigList linked list.
! Note that lists EmisList and ReadList (created lateron)  will point to the
! same containers, but will order the containers in a manner that is most
! efficient for the respective purpose.
! Argument HcoConfig represents the HEMCO configuration object. It contains
! pointers to the HEMCO traceback and error information as well as a pointer
! to ConfigList. If undefined, HcoConfig becomes initialized as part of this
! routine.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_ReadFile( am_I_Root, HcoConfig, ConfigFile, Phase, &
                              RC,        IsNest,    IsDryRun )
!
! !USES:
!
    USE inquireMod,        ONLY : findFreeLUN
    USE CharPak_Mod,       ONLY : STRREPL
    USE HCO_EXTLIST_MOD,   ONLY : AddExt, CoreNr, ExtNrInUse
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN)              :: am_I_Root  ! root CPU?
    TYPE(ConfigObj),    POINTER                 :: HcoConfig  ! HEMCO config obj
    CHARACTER(LEN=*),   INTENT(IN)              :: ConfigFile ! Full file name
    INTEGER,            INTENT(IN)              :: Phase      ! 0: all
                                                              ! 1: Settings and switches only
                                                              ! 2: fields only
    LOGICAL,            INTENT(IN),    OPTIONAL :: IsNest     ! Nested call?
    LOGICAL,            INTENT(IN),    OPTIONAL :: IsDryRun   ! Dry-run?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(INOUT)           :: RC         ! Success?
!
! !REVISION HISTORY:
!  17 Sep 2012 - C. Keller   - Initialization
!  03 Jan 2014 - C. Keller   - Now use Config_ReadCont calls.
!  30 Sep 2014 - R. Yantosca - Now declare LINE w/ 2047 characters.  This lets
!                              us handle extra-long species lists
!  13 Feb 2015 - C. Keller   - Removed section extension data: these are now
!                              listed in section base emissions.
!  11 Dec 2015 - C. Keller   - Read settings and extension switches even for
!                              nested configuration files.
!  15 Feb 2016 - C. Keller   - Now pass HcoConfig argument.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL              :: EOF
    LOGICAL              :: EXISTS, NEST, DoDryRUn
    INTEGER              :: NN
    INTEGER              :: IU_HCO, IOS

    ! Strings
    CHARACTER(LEN=255)   :: MSG,    LOC,  FileMsg
    CHARACTER(LEN=2047)  :: CFDIR
    CHARACTER(LEN=2047)  :: LINE

    !======================================================================
    ! Config_ReadFile begins here
    !======================================================================

    ! Enter
    RC  = HCO_SUCCESS
    LOC = 'Config_ReadFile (hco_config_mod.F90)'

    ! Initialize config object if not already initialized
    IF ( .NOT. ASSOCIATED(HcoConfig) ) THEN
       CALL ConfigInit( HcoConfig, RC )
       HcoConfig%ConfigFileName =  TRIM(ConfigFile)
    ENDIF

    ! Initialize
    HcoConfig%amIRoot = am_I_Root

    ! Leave here if configuration file is already read
    IF ( HcoConfig%ConfigFileRead ) THEN
       RETURN
    ENDIF

    ! Nested call?
    IF ( PRESENT( IsNest ) ) THEN
       NEST = IsNest
    ELSE
       NEST = .FALSE.
    ENDIF

    ! Is this a dry-run simulation?
    IF ( PRESENT( IsDryRun) ) THEN
       DoDryRun= IsDryRun
    ELSE
       DoDryRun = .FALSE.
    ENDIF

    ! Prompt to standard output (only on the root core
    IF ( HcoConfig%amIRoot ) THEN

       IF ( DoDryRun ) THEN

          !-----------------------------------------------------------------
          ! For dry-run simulations: state if the configuration file
          ! is found on disk, or not.  Only write this message once.
          !-----------------------------------------------------------------

          ! Test if the file exists
          INQUIRE( FILE=TRIM( ConfigFile ), EXIST=Exists )

          ! Test if the file exists and define an output string
          IF ( Exists ) THEN
             FileMsg = 'HEMCO (INIT): Opening '
          ELSE
             FileMsg = 'HEMCO (INIT): REQUIRED FILE NOT FOUND '
          ENDIF

          ! Write message to stdout
          WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( ConfigFile )
 300      FORMAT( a, ' ./', a )

       ELSE

          !-----------------------------------------------------------------
          ! For regular simulations, write a message containing
          ! the configuration file as well as the Phase value.
          !-----------------------------------------------------------------
          WRITE(6,*) ' '
          IF ( Phase == 1 ) THEN
             WRITE( 6, 310 ) TRIM(ConfigFile)
 310         FORMAT( 'Reading settings & switches of HEMCO configuration file: ', a )

          ELSEIF ( Phase == 2 ) THEN
             WRITE( 6, 320 ) TRIM(ConfigFile)
 320         FORMAT( 'Reading fields of HEMCO configuration file: ', a )

          ELSE
             WRITE( 6, 330 ) TRIM(ConfigFile)
 330         FORMAT( 'Reading entire HEMCO configuration file: ', a )
          ENDIF
       ENDIF
    ENDIF

    ! Extract configuration file directory. This is the directory containing the
    ! configuration file. Any tokens $CFDIR in the given configuration file will
    ! be replaced with the configuration file directory
    CALL HCO_GetBase( ConfigFile, CFDIR, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Find free LUN
    IU_HCO = findFreeLUN()

    INQUIRE( FILE=TRIM(ConfigFile), EXIST=EXISTS )
    IF ( .NOT. EXISTS ) THEN
       IF ( HcoConfig%amIRoot ) THEN
          WRITE(*,*) 'Cannot read file - it does not exist: ', TRIM(ConfigFile)
       ENDIF
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Open configuration file
    OPEN ( IU_HCO, FILE=TRIM( ConfigFile ), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       IF ( HcoConfig%amIRoot ) THEN
          WRITE(*,*) 'Error reading ', TRIM(ConfigFile)
       ENDIF
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Register HEMCO core as extension Nr. CoreNr (default). The
    ! core module is used by all HEMCO simulations, and the overall
    ! HEMCO settings are stored as options of this extension.
    ! Note: cannot use HCO_GetOpt('Wildcard') for species here because
    ! this is linked to the core extension...
    IF ( .NOT. ExtNrInUse( HcoConfig%ExtList, CoreNr ) ) THEN
       CALL AddExt( HcoConfig, 'CORE', CoreNr, .TRUE., 'all', RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          WRITE(*,*) 'Error adding CORE extension'
          RC = HCO_FAIL
          RETURN
       ENDIF
    ENDIF

    ! NN counts how many sections have ben read already
    NN = 0

    ! Loop until EOF
    DO

       ! Read a line from the file, exit if EOF
       CALL HCO_ReadLine ( IU_HCO, LINE, EOF, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( EOF ) EXIT

       ! Replace tab characters in LINE (if any) w/ spaces
       CALL STRREPL( LINE, HCO_TAB, HCO_SPC )

       ! Read settings if this is beginning of settings section
       ! This reads all settings into buffer and initializes the
       ! HEMCO traceback/error options.
       IF ( INDEX ( LINE, 'BEGIN SECTION SETTINGS' ) > 0 ) THEN

          IF ( PHASE < 2 ) THEN
             CALL ReadSettings( HcoConfig, IU_HCO, EOF, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
             IF ( EOF ) EXIT

             ! Increase counter
             NN = NN + 1

             ! Can we leave here?
             IF ( PHASE == 1 .AND. NN == 2 ) EXIT
          ENDIF

       ! Read extension switches. This registers all enabled extensions.
       ! This must include the core extension.
       ELSEIF ( INDEX ( LINE, 'BEGIN SECTION EXTENSION SWITCHES' ) > 0 ) THEN

          IF ( PHASE < 2 ) THEN
             CALL ExtSwitch2Buffer( HcoConfig, IU_HCO, EOF, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
             IF ( EOF ) EXIT

             ! Increase counter
             NN = NN + 1

             ! Can we leave here?
             IF ( PHASE == 1 .AND. NN == 2 ) EXIT
          ENDIF

       ! Read base emissions. This creates a new data container for each
       ! base emission field.
       ELSEIF ( INDEX ( LINE, 'BEGIN SECTION BASE EMISSIONS' ) > 0 ) THEN

          ! Read data and write into container
          ! For dry-run simulations, print name of any
          ! nested configuration files (e.g. for standalone)
          IF ( PHASE == 0 .OR. PHASE == 2 ) THEN
             CALL Config_ReadCont( HcoConfig,        IU_HCO, CFDIR,          &
                                   HCO_DCTTYPE_BASE, EOF,    RC,             &
                                   IsDryRun=IsDryRun                        )
             IF ( RC /= HCO_SUCCESS ) RETURN
             IF ( EOF ) EXIT

             ! Increase counter
             NN = NN + 1
          ENDIF

       ! Read scale factors. This creates a new data container for each
       ! scale factor.
       ELSE IF ( INDEX ( LINE, 'BEGIN SECTION SCALE FACTORS' ) > 0 ) THEN

          CALL Config_ReadCont( HcoConfig, IU_HCO, CFDIR, &
                                HCO_DCTTYPE_SCAL, EOF, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( EOF ) EXIT

       ! Read masks. This creates a new data container for each mask
       ELSE IF ( INDEX ( LINE, 'BEGIN SECTION MASKS' ) > 0 ) THEN

          IF ( PHASE == 0 .OR. PHASE == 2 ) THEN
             CALL Config_ReadCont( HcoConfig, IU_HCO, CFDIR, &
                                   HCO_DCTTYPE_MASK, EOF, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
             IF ( EOF ) EXIT

             ! Increase counter
             NN = NN + 1
          ENDIF

       ENDIF
    ENDDO

    ! Check if we caught all sections. Do that only for phase 1.
    ! Sections SETTINGS and extension switches are needed.
    IF ( PHASE == 1 .AND. NN /= 2 .AND. .NOT. NEST ) THEN
       WRITE(*,*) 'Expected 2 sections, found/read ', NN
       WRITE(*,*) 'Should read SETTINGS and EXTENSION SWITCHES'
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Close file
    CLOSE( UNIT=IU_HCO, IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       WRITE(*,*) 'Error closing ' // TRIM(ConfigFile)
       RC = HCO_FAIL
       RETURN
    ENDIF

    ! Configuration file is now read
    IF ( .NOT. NEST ) THEN
       IF ( PHASE == 0 .OR. PHASE == 2 ) THEN
          HcoConfig%ConfigFileRead = .TRUE.
       ENDIF
    ENDIF

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Config_ReadFile
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetReadList
!
! !DESCRIPTION: Subroutine SetReadList writes data to the data reading
! lists (ReadList). This routine assumes that the configuration file has
! been read beforehand (via Config\_ReadFile).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SetReadList( HcoState, RC )
!
! !USES:
!
    USE HCO_DATACONT_Mod,    ONLY : cIDList_Create
    USE HCO_READLIST_Mod,    ONLY : ReadList_Init, ReadList_Print
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState    ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC          ! Error stat
!
! !REVISION HISTORY:
!  18 Jun 2013 - C. Keller: Initialization
!  17 Sep 2013 - C. Keller: Now get data from buffer
!EOP
!------------------------------------------------------------------------------
!BOC

    CHARACTER(LEN=255)  :: MSG

    !======================================================================
    ! SetReadList begins here
    !======================================================================

    ! Init
    CALL HCO_ENTER ( HcoState%Config%Err, 'SetReadList (hco_config_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       PRINT *,'Error in HCO_ENTER called from HEMCO SetReadList'
       RETURN
    ENDIF

    ! Return w/ error if configuration file hasn't been read yet!
    IF ( .NOT. ASSOCIATED(HcoState%Config) ) THEN
       PRINT *,'HEMCO configuration object in HEMCO state is empty! Error in HEMCO SetReadList.'
       RETURN
    ENDIF
    IF ( .NOT. HcoState%Config%ConfigFileRead ) THEN
       PRINT *,'HEMCO configuration file not read! Error in HEMCO SetReadList.'
       RETURN
    ENDIF

    ! Only if not yet done so...
    IF ( .NOT. HcoState%SetReadListCalled ) THEN

       ! Initialize ReadList
       CALL ReadList_Init ( HcoState%ReadLists, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *,'Error in HCO_ENTER called from HEMCO SetReadList'
          RETURN
       ENDIF

       ! Prepare data in buffer. This call identifies all base fields
       ! that have to be read by this CPU. It also kicks out base
       ! fields for emissions with an invalid species ID (if any) or
       ! if there are other base fields with higher priority.
       CALL RegisterPrepare ( HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *,'Error in RegisterPrepare called from HEMCO SetReadList'
          RETURN
       ENDIF

       ! Register base emissions. In this step, we also redefine the
       ! list UnqScalIDs to make sure that only those scale factors
       ! will be registered that are effectively used in the next step.
       CALL Register_Base( HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *,'Error in Register_Base called from HEMCO SetReadList'
          RETURN
       ENDIF

       ! Register scale factors based upon UnqScalIDs.
       CALL Register_Scal( HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *,'Error in Register_Scal called from HEMCO SetReadList'
          RETURN
       ENDIF

       ! Create cIDList which allows quick access to all data containers
       ! based on their container IDs cID
       CALL cIDList_Create ( HcoState, HcoState%Config%ConfigList, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *,'Error in cIDList_Create called from HEMCO SetReadlist'
          RETURN
       ENDIF

       ! Don't need internal lists anymore.
       CALL ScalID_Cleanup  ( HcoState%Config%ScalIDList   )
       CALL SpecName_Cleanup( HcoState%Config%SpecNameList )

    ENDIF ! SetReadListCalled

    ! SetReadList has now been called
    HcoState%SetReadListCalled = .TRUE.

    ! Debug
    IF ( HCO_IsVerb( HcoState%Config%Err, 1 ) ) THEN
       CALL ReadList_Print( HcoState, HcoState%ReadLists, 1 )
    ENDIF

    ! Leave w/ success
    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE SetReadList
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Config_ReadCont
!
! !DESCRIPTION: Subroutine CONFIG\_READCONT reads the given line into a
! list container. Depending on the specified data type, the line is
! assumed to hold base emissions, scale factors, or mask information.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_ReadCont( HcoConfig, IU_HCO,  CFDIR,          &
                              DctType,   EOF,       RC,      IsDryRun       )
!
! !USES:
!
    USE CHARPAK_MOD,      ONLY : StrSplit
    USE HCO_EXTLIST_MOD,  ONLY : ExtNrInUse, HCO_GetOpt
    USE HCO_TIDX_Mod,     ONLY : HCO_ExtractTime
    USE HCO_FILEDATA_Mod, ONLY : FileData_Init
    USE HCO_DATACONT_Mod, ONLY : CatMax, ZeroScalID
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER       :: HcoConfig ! Config object
    INTEGER,          INTENT(IN   ) :: IU_HCO    ! Logfile LUN
    CHARACTER(LEN=*), INTENT(IN   ) :: CFDIR     ! Configuration file directory
    INTEGER,          INTENT(IN   ) :: DctType   ! 1=base; 2=scale; 3=mask
    LOGICAL,          OPTIONAL      :: IsDryRun  ! Is this a HEMCO dry-run?
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(INOUT) :: EOF       ! end of file encountered?
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT) :: RC        ! error code
!
! !REVISION HISTORY:
!  03 Jan 2014 - C. Keller - Initial version
!  29 Dec 2014 - C. Keller - Added optional 11th element for scale factors. This
!                            value will be interpreted as mask field (applied to
!                            this scale factor only).
!  27 Feb 2015 - C. Keller - Added CycleFlag 'I' (interpolation)
!  13 Mar 2015 - C. Keller - Added include files (nested configuration files)
!                            and CFDIR argument.
!  23 Sep 2015 - C. Keller - Added cycle flags 'A' and 'RA' (for averaging).
!  06 Oct 2015 - C. Keller - Added cycle flags 'EF' and 'RF' (fields must be
!                            found).
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!  20 Jul 2018 - C. Keller   - Return error if duplicate container name
!  05 Oct 2018 - R. Yantosca - Cycle flag "E" now will read the target file
!                              only once (e.g. use for restart files).
!                              Cycle flag "EC" will now continously attempt
!                              to read/query from the target file.
!  23 Oct 2018 - M. Sulprizio- Add option to use wildcard (e.g. ?ALL?) in
!                              variable name to simplify reading all species
!                              concentration fields from the GEOS-Chem restart
!                              file, but may be expanded for other purposes
!  02 Nov 2018 - M. Sulprizio- Add cycle flag "CS" to skip fields not found
!  08 Mar 2019 - M. Sulprizio- Add "*Y" options to TmCycle to force always using
!                              simulation year (eg, instead of emissions year)
!  23 Oct 2019 - M. Sulprizio- Added cycle flag "ID" to denote when dataset is
!                              discontinous and needs to be interpolated
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                   :: I, N, nEdges
    INTEGER                   :: nScl
    INTEGER                   :: STAT
    INTEGER                   :: Int1
    INTEGER                   :: Int2
    INTEGER                   :: Int3
    INTEGER                   :: Int4
    INTEGER                   :: nCat
    INTEGER                   :: Cats(CatMax)
    INTEGER                   :: STRLEN
    INTEGER                   :: levScal1
    INTEGER                   :: levScal2
    INTEGER                   :: nTags
    LOGICAL                   :: SKIP
    LOGICAL                   :: Found
    CHARACTER(LEN= 63)        :: cName
    CHARACTER(LEN=255)        :: srcFile
    CHARACTER(LEN= 50)        :: srcVar
    CHARACTER(LEN= 31)        :: srcTime
    CHARACTER(LEN= 31)        :: TmCycle
    CHARACTER(LEN=  1)        :: WildCard
    CHARACTER(LEN=  1)        :: Separator
    CHARACTER(LEN= 31)        :: srcDim
    CHARACTER(LEN= 31)        :: srcUnit
    CHARACTER(LEN= 31)        :: SpcName
    CHARACTER(LEN=255)        :: Char1
    CHARACTER(LEN=255)        :: Char2
    CHARACTER(LEN=255)        :: LOC, MSG
    CHARACTER(LEN=255)        :: LINE
    CHARACTER(LEN=255)        :: tagId
    CHARACTER(LEN=255)        :: tagName
    CHARACTER(LEN=255)        :: tagcName
    CHARACTER(LEN=255)        :: ItemPrefix
    CHARACTER(LEN=255)        :: ErrMsg

    ! Arrays
    INTEGER                   :: SplitInts(255)
    CHARACTER(LEN=255)        :: SubStrs(255)

    ! Pointers
    TYPE(ListCont), POINTER   :: Lct
    TYPE(ListCont), POINTER   :: Tmp
    TYPE(FileData), POINTER   :: Dta

    !=================================================================
    ! Config_ReadCont begins here!
    !=================================================================

    ! Enter
    LOC = 'Config_ReadCont (hco_config_mod.F90)'

    ! Initialize
    SKIP           = .FALSE.
    nCat           = -1
    Lct            => NULL()
    Tmp            => NULL()
    Dta            => NULL()

    ! Get tokens
    WildCard  = HCO_GetOpt( HcoConfig%ExtList, 'Wildcard'  )
    Separator = HCO_GetOpt( HcoConfig%ExtList, 'Separator' )

    ! Repeat until end of the given section is found
    DO

       !==============================================================
       ! Read line and get desired character strings
       ! Since base emissions, scale factors and masks have different
       ! configuration file input parameter, need to use a different
       ! call for the three data types.
       !==============================================================
       IF ( DctType == HCO_DCTTYPE_BASE ) THEN
          CALL ReadAndSplit_Line ( HcoConfig, IU_HCO, cName,    2,  &
                                   srcFile,   3,      srcVar,   4,  &
                                   srcTime,   5,      TmCycle,  6,  &
                                   srcDim,    7,      srcUnit,  8,  &
                                   SpcName,   9,      Char1,   10,  &
                                   Char2,    11,                    &
                                   Int1,     -1,      Int2,    12,  &
                                   Int3,      1,      STAT,         &
                                   OutLine=LINE                      )

       ELSEIF ( DctType == HCO_DCTTYPE_SCAL ) THEN
          CALL ReadAndSplit_Line ( HcoConfig, IU_HCO, cName,    2,  &
                                   srcFile,   3,      srcVar,   4,  &
                                   srcTime,   5,      TmCycle,  6,  &
                                   srcDim,    7,      srcUnit,  8,  &
                                   SpcName,  -1,      Char1,   -1,  &
                                   Char2,    -1,                    &
                                   Int1,      1,      Int2,     9,  &
                                   Int3,     10,      STAT,         &
                                   optcl=10    ,      OutLine=LINE   )

       ELSEIF ( DctType == HCO_DCTTYPE_MASK ) THEN
          CALL ReadAndSplit_Line ( HcoConfig, IU_HCO, cName,    2,  &
                                   srcFile,   3,      srcVar,   4,  &
                                   srcTime,   5,      TmCycle,  6,  &
                                   srcDim,    7,      srcUnit,  8,  &
                                   SpcName,  -1,      Char1,   10,  &
                                   Char2,    11,                    &
                                   Int1,      1,      Int2,     9,  &
                                   Int3,     -1,      STAT,         &
                                   optcl=11,          OutLine=LINE   )
       ENDIF

       !--------------------------------------------------------------
       ! Error checks
       !--------------------------------------------------------------

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
          CALL HCO_ERROR ( HcoConfig%Err, 'STAT == 100', RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! -------------------------------------------------------------
       ! Check for emission shortcuts. Fields can be bracketed into
       ! 'collections'.
       ! -------------------------------------------------------------
       CALL BracketCheck( HcoConfig, STAT, LINE, SKIP, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Skip if needed
       IF ( SKIP ) CYCLE

       ! Can advance to next line if this was a bracket line: nothing
       ! else to do with this line.
       IF ( STAT == 5 .OR. STAT == 6 ) CYCLE

       ! Read include file. Configuration files can be 'nested', e.g.
       ! configuration files can be included into the 'main' configuration
       ! file. These files must be listed as '>>>include FileName', where
       ! FileName is the actual path of the file.
       IF ( STAT == 1000 ) THEN

          ! Call the parser. This is to make sure that any $ROOT statements
          ! will be evaluated properly. The configuration file must not
          ! contain any data tokens ($YR, $MM, etc.).
          CALL HCO_CharParse ( HcoConfig, LINE, 0, 0, 0, 0, 0, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          CALL Config_ReadFile( HcoConfig%amIRoot, HcoConfig, LINE, 0, RC,  &
                                IsNest=.TRUE., IsDryRun=IsDryRun            )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! All done with this line
          CYCLE
       ENDIF

       ! Output status should be 0 if none of the statuses above applies
       IF ( STAT /= 0 ) THEN
          CALL HCO_ERROR ( HcoConfig%Err, 'STAT /= 0', RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! For base fields, check if this extension number is indeed in
       ! use. Otherwise, we can ignore this line completely! The
       ! extension switches are read and evaluated prior to the
       ! extension data!
       IF ( DctType == HCO_DCTTYPE_BASE .AND. &
            .NOT. ExtNrInUse( HcoConfig%ExtList, Int3 ) ) CYCLE

       !==============================================================
       ! Create and fill list container and add to ConfigList
       !==============================================================

       ! Test if wildcard is present
       IF ( INDEX( srcVar, '?' ) > 0 ) THEN

          ! Split the name to get wildcard and string prior to wildcard
          CALL StrSplit( srcVar, '?', SubStrs, N )
          tagId = SubStrs(N-1)
          ItemPrefix = SubStrs(1)

          ! Get number of tags for this wildcard
          CALL Hco_GetTagInfo( tagId, HcoConfig, Found, RC, nTags=nTags )

          ! Add each tagged name as a separate item in the collection
          DO N = 1, nTags
             ! Construct the item name
             tagcName = ''

             ! Get tag, if any
             CALL Hco_GetTagInfo( tagId, HcoConfig, Found, RC, N=N, &
                                  tagName=tagName )
             IF ( RC /= HCO_SUCCESS ) THEN
                ErrMsg = 'Error retrieving tag name for' //            &
                         ' wildcard ' // TRIM(tagId)
                CALL HCO_Error( HcoConfig%Err, ErrMsg, RC )
                RETURN
             ENDIF

             ! Append the tag name to the output name
             srcVar   = TRIM( ItemPrefix ) // TRIM( tagName )
             tagcName = TRIM( cName      ) // TRIM( tagName )
             ! Do not overwrite species name for now. Use HEMCO wildcard
             ! to read all external model species.
             !spcName  = TRIM( tagName    )

             ! -------------------------------------------------------------
             ! Fill data container
             ! -------------------------------------------------------------

             ! Add blank list container to ConfigList list. The container
             ! is placed at the beginning of the list.
             CALL ConfigList_AddCont ( Lct, HcoConfig%ConfigList )

             ! Check if name exists already
             CALL CheckForDuplicateName( HcoConfig, tagcName, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN

             ! Attributes used by all data types: data type number and
             ! container name.
             Lct%Dct%DctType      = DctType
             Lct%Dct%cName        = ADJUSTL(tagcName)

             ! Set species name, extension number, emission category,
             ! hierarchy
             Lct%Dct%SpcName       = ADJUSTL(SpcName)
             Lct%Dct%Hier          = Int2
             Lct%Dct%ExtNr         = Int3

             ! Extract category from character 2. This can be up to
             ! CatMax integers, or empty.
             CALL HCO_CharSplit( Char2, Separator, Wildcard, Cats, nCat, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
             IF ( nCat == 0 ) THEN
                Lct%Dct%Cat = -999
             ELSE
                Lct%Dct%Cat = Cats(1)
             ENDIF

             ! Set scale factor IDs into Scal_cID. These values will be
             ! replaced lateron with the container IDs (in register_base)!
             CALL HCO_CharSplit( Char1, Separator, Wildcard, &
                                 SplitInts, nScl, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
             IF ( nScl > 0 ) THEN
                ALLOCATE ( Lct%Dct%Scal_cID(nScl) )
                Lct%Dct%Scal_cID(1:nScl) = SplitInts(1:nScl)
                Lct%Dct%nScalID          = nScl
             ENDIF

             ! Register species name. A list of all species names can be
             ! returned to the atmospheric model to match HEMCO species
             ! with model species (see Config\_GetSpecNames).
             CALL SpecName_Register ( HcoConfig, ADJUSTL(SpcName), RC )
             IF ( RC /= HCO_SUCCESS ) RETURN

             ! -------------------------------------------------------------
             ! Create and fill file data object. Use previous file data
             ! object if filename is undefined. Do not yet update the
             ! DoShare and DtaHome flags of the data file object and data
             ! container, respectively, since we still don't know which
             ! data containers will be effectively used for emission
             ! calculation (containers may be dropped lateron because the
             ! emission category / hierarchy is too low, species is not
             ! used, etc.). The DoShare and DtaHome flags will be set the
             ! first time that data is read (in hco_readlist_mod.F90).
             ! Here, we only set the DtaHome flag to -1000 instead of the
             ! default value of -999 to be able to identify data objects
             ! used by multiple containers.
             ! -------------------------------------------------------------
             IF ( TRIM(srcFile) == '-' ) THEN
                IF ( .NOT. ASSOCIATED(Dta) ) THEN
                   MSG = 'Cannot use previous data container: '//TRIM(tagcName)
                   CALL HCO_ERROR ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
                   RETURN
                ENDIF
                Lct%Dct%DtaHome = Lct%Dct%DtaHome - 1
             ELSE
                Dta => NULL()
                CALL FileData_Init ( Dta )

                ! Set source file name. Check if the read file name starts
                ! with the configuration file token '$CFDIR', in which case
                ! we replace this value with the passed CFDIR value.
                STRLEN = LEN(srcFile)
                IF ( STRLEN > 6 ) THEN
                   IF ( srcFile(1:6) == '$CFDIR' ) THEN
                      srcFile = TRIM(CFDIR) // TRIM(srcFile(7:STRLEN))
                   ENDIF
                ENDIF
                Dta%ncFile    = srcFile

                ! Set source variable and original data unit.
                Dta%ncPara    = ADJUSTL(srcVar)
                Dta%OrigUnit  = ADJUStL(srcUnit)

                ! If the parameter ncPara is not defined, attempt to read data
                ! directly from configuration file instead of netCDF.
                ! These data are always assumed to be in local time. Gridded
                ! data read from netCDF is always in UTC, except for weekdaily
                ! data that is treated in local time. The corresponding
                ! IsLocTime flag is updated when reading the data (see
                ! hcoio_dataread_mod.F90).
                IF ( TRIM(Dta%ncPara) == '-' ) THEN
                   Dta%ncRead    = .FALSE.
                   Dta%IsLocTime = .TRUE.
                ENDIF

                ! Extract information from time stamp character and pass values
                ! to the corresponding container variables. If no time string is
                ! defined, keep default values (-1 for all of them)
                IF ( TRIM(srcTime) /= '-' ) THEN
                   CALL HCO_ExtractTime( HcoConfig, srcTime, Dta, RC )
                   IF ( RC /= HCO_SUCCESS ) RETURN
                ENDIF

                ! In an ESMF environment, the source data will be imported
                ! through ExtData by name, hence need to set ncFile equal to
                ! container name!
#if defined(ESMF_)
                IF ( Dta%ncRead ) THEN
                   Dta%ncFile = ADJUSTL(tagcName)
                ENDIF
#endif

                ! Set time cycling behaviour. Possible values are:
                ! - "C"  : cycling <-- DEFAULT
                ! - "CS" : cycling, skip if not exist
                ! - "CY" : cycling, always use simulation year
                ! - "CYS": cycling, always use simulation yr, skip if not exist
                ! - "R"  : range
                ! - "RA" : range, average outside
                ! - "RF" : range, forced (error if not in range)
                ! - "RFY": range, forced, always use simulation year
                ! - "RFY3: range, forced, always use simulation year, 3-hourly
                ! - "RY" : range, always use simulation year
                ! - "E"  : exact (read file once)
                ! - "EF" : exact, forced (error if not exist, read/query once)
                ! - "EC" : exact (read/query continuously, e.g. for ESMF interface)
                ! - "ECF": exact, forced (error if not exist, read/query continuously)
                ! - "EY" : exact, always use simulation year
                ! - "A"  : average
                ! - "I"  : interpolate
                ! - "ID" : interpolate, discontinuous dataset
                Dta%MustFind  = .FALSE.
                Dta%UseSimYear= .FALSE.
                Dta%Discontinuous = .FALSE.
                IF ( TRIM(TmCycle) == "C" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_CYCLE
                   Dta%MustFind  = .TRUE.
                ELSEIF ( TRIM(TmCycle) == "CS" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_CYCLE
                   Dta%MustFind  = .FALSE.
                ELSEIF ( TRIM(TmCycle) == "CY" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_CYCLE
                   Dta%MustFind  = .TRUE.
                   Dta%UseSimYear= .TRUE.
                ELSEIF ( TRIM(TmCycle) == "CYS" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_CYCLE
                   Dta%MustFind  = .FALSE.
                   Dta%UseSimYear= .TRUE.
                ELSEIF ( TRIM(TmCycle) == "R" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_RANGE
                ELSEIF ( TRIM(TmCycle) == "RA" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_RANGEAVG
                ELSEIF ( TRIM(TmCycle) == "RF" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_RANGE
                   Dta%MustFind  = .TRUE.
                ELSEIF ( TRIM(TmCycle) == "RFY" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_RANGE
                   Dta%MustFind  = .TRUE.
                   Dta%UseSimYear= .TRUE.
                ELSEIF ( TRIM(TmCycle) == "RFY3" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_RANGE
                   Dta%MustFind  = .TRUE.
                   Dta%UseSimYear= .TRUE.
                   Dta%UpdtFlag  = HCO_UFLAG_3HR
                ELSEIF ( TRIM(TmCycle) == "RY" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_RANGE
                   Dta%UseSimYear= .TRUE.
                ELSEIF ( TRIM(TmCycle) == "E" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_EXACT
                   Dta%UpdtFlag  = HCO_UFLAG_ONCE
                ELSEIF ( TRIM(TmCycle) == "EF" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_EXACT
                   Dta%UpdtFlag  = HCO_UFLAG_ONCE
                   Dta%MustFind  = .TRUE.
                ELSEIF ( TRIM(TmCycle) == "EC" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_EXACT
                ELSEIF ( TRIM(TmCycle) == "ECF" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_EXACT
                   Dta%MustFind  = .TRUE.
                ELSEIF ( TRIM(TmCycle) == "EY" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_EXACT
                   Dta%UpdtFlag  = HCO_UFLAG_ONCE
                   Dta%UseSimYear=.TRUE.
                ELSEIF ( TRIM(TmCycle) == "A" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_AVERG
                ELSEIF ( TRIM(TmCycle) == "I" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_INTER
                ELSEIF ( TRIM(TmCycle) == "ID" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_INTER
                   Dta%Discontinuous = .TRUE.
                ELSEIF ( TRIM(TmCycle) == "-" ) THEN
                   Dta%CycleFlag = HCO_CFLAG_CYCLE
                ELSE
                   MSG = 'Invalid time cycling attribute: ' // &
                         TRIM(TmCycle) // ' - in ' // TRIM(tagcName)
                   CALL HCO_ERROR ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
                   RETURN
                ENDIF

                ! Set space dimension. This will determine the dimension of the
                ! data array vector, i.e. 3D or 2D. Different time slices will
                ! be stored as different vector elements.
                ! For 3D data, it is now possible to explicitly set the number
                ! of vertical levels to be used, as well as the 'reading
                ! direction' (up or down). These information is also extracted
                ! from srcDim and will be stored in variable Dta%Levels.
                ! (ckeller, 5/20/15)
                ! ExtractSrcDim now also returns possible scale factors for the
                ! injection level, which will be stored in container variable
                ! levScalID1 (bottom level) and levScalID2 (top level).
                CALL ExtractSrcDim( HcoConfig, srcDim, Dta, &
                                    levScal1,  levScal2,  RC )
                IF ( RC /= HCO_SUCCESS ) RETURN

                ! Set level scale factor index
                IF ( levScal1 > 0 ) Lct%Dct%levScalID1 = levScal1
                IF ( levScal2 > 0 ) Lct%Dct%levScalID2 = levScal2

                ! For scale factors: check if a mask is assigned to this scale
                ! factor. In this case, pass mask ID to first slot of Scal_cID
                ! vector. This value will be set to the container ID of the
                ! corresponding mask field lateron.
                IF ( DctType == HCO_DCTTYPE_SCAL .AND. Int3 > 0 ) THEN
                   ALLOCATE ( Lct%Dct%Scal_cID(1) )
                   Lct%Dct%Scal_cID(1) = Int3
                   Lct%Dct%nScalID     = 1
                ENDIF

                ! For masks: extract grid box edges. These will be used later
                ! on to determine if emissions have to be considered by this
                ! CPU.
                IF ( DctType == HCO_DCTTYPE_MASK ) THEN

                   ! Extract grid box edges. Need to be four values.
                   CALL HCO_CharSplit ( Char1, Separator, Wildcard, &
                                        SplitInts, nEdges, RC )
                   IF ( RC /= HCO_SUCCESS ) RETURN
                   IF ( nEdges /= 4 ) THEN
                      MSG = 'Cannot properly read mask coverage: ' // &
                           TRIM(Lct%Dct%cName)
                      CALL HCO_ERROR ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
                      RETURN
                   ENDIF

                   ! Save temporarily in year and month range. Will be
                   ! reset lateron.
                   Dta%ncYrs(1) = SplitInts(1)
                   Dta%ncYrs(2) = SplitInts(2)
                   Dta%ncMts(1) = SplitInts(3)
                   Dta%ncMts(2) = SplitInts(4)

                   ! Make sure that masks are always being read if specified so.
                   IF ( Char2(1:1) == 'y' .OR. Char2(1:1) == 'Y' ) THEN
                      CALL ScalID2List( HcoConfig%ScalIDList, Lct%Dct%ScalID, RC )
                      IF ( RC /= HCO_SUCCESS ) RETURN
                   ENDIF
                ENDIF
             ENDIF

             ! Connect file data object of this data container.
             Lct%Dct%Dta => Dta

             ! Free list container for next cycle
             Lct => NULL()

          ENDDO

       ! If no wildcard is present in variable name
       ELSE

          ! -------------------------------------------------------------
          ! Fill data container
          ! -------------------------------------------------------------

          ! Add blank list container to ConfigList list. The container
          ! is placed at the beginning of the list.
          CALL ConfigList_AddCont ( Lct, HcoConfig%ConfigList )

          ! Check if name exists already
          CALL CheckForDuplicateName( HcoConfig, cName, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Attributes used by all data types: data type number and
          ! container name.
          Lct%Dct%DctType      = DctType
          Lct%Dct%cName        = ADJUSTL(cName)

          ! Base container specific attributes
          IF ( DctType == HCO_DCTTYPE_BASE ) THEN

             ! Set species name, extension number, emission category,
             ! hierarchy
             Lct%Dct%SpcName       = ADJUSTL(SpcName)
             Lct%Dct%Hier          = Int2
             Lct%Dct%ExtNr         = Int3

             ! Extract category from character 2. This can be up to
             ! CatMax integers, or empty.
             CALL HCO_CharSplit( Char2, Separator, Wildcard, Cats, nCat, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
             IF ( nCat == 0 ) THEN
                Lct%Dct%Cat = -999
             ELSE
                Lct%Dct%Cat = Cats(1)
             ENDIF

             ! Set scale factor IDs into Scal_cID. These values will be
             ! replaced lateron with the container IDs (in register_base)!
             CALL HCO_CharSplit( Char1, Separator, Wildcard, &
                                 SplitInts, nScl, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
             IF ( nScl > 0 ) THEN
                ALLOCATE ( Lct%Dct%Scal_cID(nScl) )
                Lct%Dct%Scal_cID(1:nScl) = SplitInts(1:nScl)
                Lct%Dct%nScalID          = nScl
             ENDIF

             ! Register species name. A list of all species names can be
             ! returned to the atmospheric model to match HEMCO species
             ! with model species (see Config\_GetSpecNames).
             CALL SpecName_Register ( HcoConfig, ADJUSTL(SpcName), RC )
             IF ( RC /= HCO_SUCCESS ) RETURN

          ! Scale factor & mask specific attributes
          ELSE IF ( DctType == HCO_DCTTYPE_SCAL .OR. &
                    DctType == HCO_DCTTYPE_MASK       ) THEN

             ! Set scale factor ID and data operator
             Lct%Dct%ScalID = Int1
             Lct%Dct%Oper   = Int2

             ! Make sure that negative scale factors are always read
             IF ( Lct%Dct%ScalID < 0 ) THEN
                CALL ScalID2List( HcoConfig%ScalIDList, Lct%Dct%ScalID, RC )
                IF ( RC /= HCO_SUCCESS ) RETURN
             ENDIF

          ELSE
             CALL HCO_ERROR ( HcoConfig%Err, 'Invalid data type!', RC, &
                              THISLOC=LOC )
             RETURN
          ENDIF

          ! -------------------------------------------------------------
          ! Create and fill file data object. Use previous file data
          ! object if filename is undefined. Do not yet update the
          ! DoShare and DtaHome flags of the data file object and data
          ! container, respectively, since we still don't know which
          ! data containers will be effectively used for emission
          ! calculation (containers may be dropped lateron because the
          ! emission category / hierarchy is too low, species is not
          ! used, etc.). The DoShare and DtaHome flags will be set the
          ! first time that data is read (in hco_readlist_mod.F90).
          ! Here, we only set the DtaHome flag to -1000 instead of the
          ! default value of -999 to be able to identify data objects
          ! used by multiple containers.
          ! -------------------------------------------------------------
          IF ( TRIM(srcFile) == '-' ) THEN
             IF ( .NOT. ASSOCIATED(Dta) ) THEN
                MSG = 'Cannot use previous data container: '//TRIM(cName)
                CALL HCO_ERROR ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             Lct%Dct%DtaHome = Lct%Dct%DtaHome - 1
          ELSE
             Dta => NULL()
             CALL FileData_Init ( Dta )

             ! Set source file name. Check if the read file name starts
             ! with the configuration file token '$CFDIR', in which case
             ! we replace this value with the passed CFDIR value.
             STRLEN = LEN(srcFile)
             IF ( STRLEN > 6 ) THEN
                IF ( srcFile(1:6) == '$CFDIR' ) THEN
                   srcFile = TRIM(CFDIR) // TRIM(srcFile(7:STRLEN))
                ENDIF
             ENDIF
             Dta%ncFile    = srcFile

             ! Set source variable and original data unit.
             Dta%ncPara    = ADJUSTL(srcVar)
             Dta%OrigUnit  = ADJUStL(srcUnit)

             ! If the parameter ncPara is not defined, attempt to read data
             ! directly from configuration file instead of netCDF.
             ! These data are always assumed to be in local time. Gridded
             ! data read from netCDF is always in UTC, except for weekdaily
             ! data that is treated in local time. The corresponding
             ! IsLocTime flag is updated when reading the data (see
             ! hcoio_dataread_mod.F90).
             IF ( TRIM(Dta%ncPara) == '-' ) THEN
                Dta%ncRead    = .FALSE.
                Dta%IsLocTime = .TRUE.
             ENDIF

             ! Extract information from time stamp character and pass values
             ! to the corresponding container variables. If no time string is
             ! defined, keep default values (-1 for all of them)
             IF ( TRIM(srcTime) /= '-' ) THEN
                CALL HCO_ExtractTime( HcoConfig, srcTime, Dta, RC )
                IF ( RC /= HCO_SUCCESS ) RETURN
             ENDIF

             ! In an ESMF environment, the source data will be imported
             ! through ExtData by name, hence need to set ncFile equal to
             ! container name!
#if defined(ESMF_)
             IF ( Dta%ncRead ) THEN
                Dta%ncFile = ADJUSTL(cName)
             ENDIF
#endif

             ! Set time cycling behaviour. Possible values are:
             ! - "C"  : cycling <-- DEFAULT
             ! - "CS" : cycling, skip if not exist
             ! - "CY" : cycling, always use simulation year
             ! - "CYS": cycling, always use simulation yr, skip if not exist
             ! - "R"  : range
             ! - "RA" : range, average outside
             ! - "RF" : range, forced (error if not in range)
             ! - "RFY": range, forced, always use simulation year
             ! - "RFY3: range, forced, always use simulation year, 3-hourly
             ! - "RY" : range, always use simulation year
             ! - "E"  : exact (read file once)
             ! - "EF" : exact, forced (error if not exist, read/query once)
             ! - "EC" : exact (read/query continuously, e.g. for ESMF interface)
             ! - "ECF": exact, forced (error if not exist, read/query continuously)
             ! - "EY" : exact, always use simulation year
             ! - "A"  : average
             ! - "I"  : interpolate
             Dta%MustFind  = .FALSE.
             Dta%UseSimYear= .FALSE.
             IF ( TRIM(TmCycle) == "C" ) THEN
                Dta%CycleFlag = HCO_CFLAG_CYCLE
                Dta%MustFind  = .TRUE.
             ELSEIF ( TRIM(TmCycle) == "CS" ) THEN
                Dta%CycleFlag = HCO_CFLAG_CYCLE
                Dta%MustFind  = .FALSE.
             ELSEIF ( TRIM(TmCycle) == "CY" ) THEN
                Dta%CycleFlag = HCO_CFLAG_CYCLE
                Dta%MustFind  = .TRUE.
                Dta%UseSimYear= .TRUE.
             ELSEIF ( TRIM(TmCycle) == "CYS" ) THEN
                Dta%CycleFlag = HCO_CFLAG_CYCLE
                Dta%MustFind  = .FALSE.
                Dta%UseSimYear= .TRUE.
             ELSEIF ( TRIM(TmCycle) == "R" ) THEN
                Dta%CycleFlag = HCO_CFLAG_RANGE
             ELSEIF ( TRIM(TmCycle) == "RA" ) THEN
                Dta%CycleFlag = HCO_CFLAG_RANGEAVG
             ELSEIF ( TRIM(TmCycle) == "RF" ) THEN
                Dta%CycleFlag = HCO_CFLAG_RANGE
                Dta%MustFind  = .TRUE.
             ELSEIF ( TRIM(TmCycle) == "RFY" ) THEN
                Dta%CycleFlag = HCO_CFLAG_RANGE
                Dta%MustFind  = .TRUE.
                Dta%UseSimYear= .TRUE.
             ELSEIF ( TRIM(TmCycle) == "RFY3" ) THEN
                Dta%CycleFlag = HCO_CFLAG_RANGE
                Dta%MustFind  = .TRUE.
                Dta%UseSimYear= .TRUE.
                Dta%UpdtFlag  = HCO_UFLAG_3HR
             ELSEIF ( TRIM(TmCycle) == "RY" ) THEN
                Dta%CycleFlag = HCO_CFLAG_RANGE
                Dta%UseSimYear= .TRUE.
             ELSEIF ( TRIM(TmCycle) == "E" ) THEN
                Dta%CycleFlag = HCO_CFLAG_EXACT
                Dta%UpdtFlag  = HCO_UFLAG_ONCE
             ELSEIF ( TRIM(TmCycle) == "EF" ) THEN
                Dta%CycleFlag = HCO_CFLAG_EXACT
                Dta%UpdtFlag  = HCO_UFLAG_ONCE
                Dta%MustFind  = .TRUE.
             ELSEIF ( TRIM(TmCycle) == "EC" ) THEN
                Dta%CycleFlag = HCO_CFLAG_EXACT
             ELSEIF ( TRIM(TmCycle) == "ECF" ) THEN
                Dta%CycleFlag = HCO_CFLAG_EXACT
                Dta%MustFind  = .TRUE.
             ELSEIF ( TRIM(TmCycle) == "EY" ) THEN
                Dta%CycleFlag = HCO_CFLAG_EXACT
                Dta%UpdtFlag  = HCO_UFLAG_ONCE
                Dta%UseSimYear=.TRUE.
             ELSEIF ( TRIM(TmCycle) == "A" ) THEN
                Dta%CycleFlag = HCO_CFLAG_AVERG
             ELSEIF ( TRIM(TmCycle) == "I" ) THEN
                Dta%CycleFlag = HCO_CFLAG_INTER
             ELSEIF ( TRIM(TmCycle) == "-" ) THEN
                Dta%CycleFlag = HCO_CFLAG_CYCLE
             ELSE
                MSG = 'Invalid time cycling attribute: ' // &
                     TRIM(TmCycle) // ' - in ' // TRIM(tagcName)
                CALL HCO_ERROR ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF

             ! Set space dimension. This will determine the dimension of the
             ! data array vector, i.e. 3D or 2D. Different time slices will
             ! be stored as different vector elements.
             ! For 3D data, it is now possible to explicitly set the number
             ! of vertical levels to be used, as well as the 'reading
             ! direction' (up or down). These information is also extracted
             ! from srcDim and will be stored in variable Dta%Levels.
             ! (ckeller, 5/20/15)
             ! ExtractSrcDim now also returns possible scale factors for the
             ! injection level, which will be stored in container variable
             ! levScalID1 (bottom level) and levScalID2 (top level).
             CALL ExtractSrcDim( HcoConfig, srcDim, Dta, &
                                 levScal1,  levScal2,  RC )
             IF ( RC /= HCO_SUCCESS ) RETURN

             ! Set level scale factor index
             IF ( levScal1 > 0 ) Lct%Dct%levScalID1 = levScal1
             IF ( levScal2 > 0 ) Lct%Dct%levScalID2 = levScal2

             ! For scale factors: check if a mask is assigned to this scale
             ! factor. In this case, pass mask ID to first slot of Scal_cID
             ! vector. This value will be set to the container ID of the
             ! corresponding mask field lateron.
             IF ( DctType == HCO_DCTTYPE_SCAL .AND. Int3 > 0 ) THEN
                ALLOCATE ( Lct%Dct%Scal_cID(1) )
                Lct%Dct%Scal_cID(1) = Int3
                Lct%Dct%nScalID     = 1
             ENDIF

             ! For masks: extract grid box edges. These will be used later
             ! on to determine if emissions have to be considered by this
             ! CPU.
             IF ( DctType == HCO_DCTTYPE_MASK ) THEN

                ! Extract grid box edges. Need to be four values.
                CALL HCO_CharSplit ( Char1, Separator, Wildcard, &
                                     SplitInts, nEdges, RC )
                IF ( RC /= HCO_SUCCESS ) RETURN
                IF ( nEdges /= 4 ) THEN
                   MSG = 'Cannot properly read mask coverage: ' // &
                         TRIM(Lct%Dct%cName)
                   CALL HCO_ERROR ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
                   RETURN
                ENDIF

                ! Save temporarily in year and month range. Will be
                ! reset lateron.
                Dta%ncYrs(1) = SplitInts(1)
                Dta%ncYrs(2) = SplitInts(2)
                Dta%ncMts(1) = SplitInts(3)
                Dta%ncMts(2) = SplitInts(4)

                ! Make sure that masks are always being read if specified so.
                IF ( Char2(1:1) == 'y' .OR. Char2(1:1) == 'Y' ) THEN
                   CALL ScalID2List( HcoConfig%ScalIDList, Lct%Dct%ScalID, RC )
                   IF ( RC /= HCO_SUCCESS ) RETURN
                ENDIF
             ENDIF
          ENDIF

          ! Connect file data object of this data container.
          Lct%Dct%Dta => Dta

          ! If a base emission field covers multiple emission categories,
          ! create a 'shadow' container for each additional category.
          ! These shadow container have the same information as the main
          ! container except that a scale factor of zero will be applied in
          ! addition. This makes sure that the inventory cancels out other
          ! inventories with lower hierarchy for every specified category,
          ! while emission totals are not changed. All emissions of a base
          ! field with multiple categories is written into the category
          ! listed first.
          IF ( nCat > 1 ) THEN

             ! nCat cannot exceed CatMax
             IF ( nCat > CatMax ) THEN
                MSG = 'Category max exceeded'
                CALL HCO_ERROR ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF

             CALL AddShadowFields( HcoConfig, Lct, Cats, nCat, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN

             ! Reset nCat
             nCat = -1
          ENDIF

          ! Free list container for next cycle
          Lct => NULL()

       ENDIF

    ENDDO

    ! Leave w/ success
    Dta => NULL()
    RC  =  HCO_SUCCESS

  END SUBROUTINE Config_ReadCont
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: BracketCheck
!
! !DESCRIPTION: Subroutine BracketCheck checks if base emission data is within
! a bracket and if that field shall be ignored or not. Brackets can be used to
! lump entires of the HEMCO configuration file into collections that can be
! collectively enabled or disabled. The first entry of a collection is marked
! adding an 'opening bracket' to the HEMCO configuration file (on the line
! above the entry). Opening brackets must start with three opening brackets,
! e.g.: '(((TEST'. Similarly, the end of a collection is marked by placing a
! closing bracket after the last entry of the collection: '))))TEST'.
! Brackets can be enabled / disabled in the EXTENSION SWITCH section of the
! HEMCO configuration file:
! \# ExtNr ExtName           on/off  Species
! 0       Base              : on    *
!     --> TEST              :       true
!\\
!\\
! It is also possible to use 'opposite' brackets, e.g. to use a collection
! only if the given setting is *disabled*. This can be achieved by precede
! the collection word with '.not.', e.g. '(((.not.TEST' and '))).not.TEST'.
! Similarly, multiple collections can be combined to be evaluated together,
! e.g. NAME1.or.NAME2.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE BracketCheck( HcoConfig, STAT, LINE, SKIP, RC )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,  ONLY : GetExtOpt, GetExtNr
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN)    :: STAT        !
    CHARACTER(LEN=*), INTENT(IN)    :: LINE        !
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER       :: HcoConfig ! Config object
    LOGICAL,          INTENT(INOUT) :: SKIP        ! Skip
    INTEGER,          INTENT(INOUT) :: RC          ! Success/failure
!
! !REVISION HISTORY:
!  15 Feb 2015 - C. Keller   - Initial version.
!  12 Mar 2015 - C. Keller   - Added 'mirror' option.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Maximum number of nested brackets
    INTEGER, PARAMETER            :: MAXBRACKNEST = 5
    INTEGER                       :: IDX, STRLEN, ExtNr
    LOGICAL                       :: FOUND
    LOGICAL                       :: UseBracket, UseThis
    LOGICAL                       :: verb
    LOGICAL                       :: REV
    INTEGER, SAVE                 :: NEST      = 0
    INTEGER, SAVE                 :: SKIPLEVEL = 0
    CHARACTER(LEN=255), SAVE      :: AllBrackets(MAXBRACKNEST) = ''
    CHARACTER(LEN=255)            :: TmpBracket, CheckBracket, ThisBracket
    CHARACTER(LEN=255)            :: MSG

    CHARACTER(LEN=255), PARAMETER :: LOC = 'BracketCheck (hco_config_mod.F90)'

    !======================================================================
    ! BracketCheck begins here
    !======================================================================

    ! Init
    verb = HCO_IsVerb( HcoConfig%Err, 1 )

    ! Get name of this bracket
    IF ( STAT == 5 .OR. STAT == 6 ) THEN
       STRLEN     = LEN(LINE)
       IF ( STRLEN < 4 ) THEN
          CALL HCO_ERROR ( HcoConfig%Err, &
                          'Illegal bracket length: '//TRIM(LINE), &
                           RC, THISLOC=LOC )
          RETURN
       ELSE
          TmpBracket = TRIM(LINE(4:STRLEN))
       ENDIF
    ENDIF

    ! Open a bracket. Save out the bracket name in the list of all
    ! opened brackets. This is primarily to ensure that every opening
    ! brackets is properly closed. Only register it as skipping bracket
    ! if needed.
    IF ( STAT == 5 ) THEN

       ! Archive bracket name
       NEST = NEST + 1
       IF ( NEST > MAXBRACKNEST ) THEN
          MSG = 'Too many nested brackets'
          CALL HCO_ERROR( HcoConfig%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       AllBrackets(NEST) = TmpBracket

       ! Check if this bracket content shall be skipped. Always skip
       ! if this is a nested bracket in an already skipped bracket.
       IF ( .NOT. SKIP ) THEN

          ! Check for 'inverse' bracket. These start with '.not.'
          CheckBracket = TmpBracket
          REV          = .FALSE.
          IF ( STRLEN > 5 ) THEN
             IF ( TmpBracket(1:5) == '.not.' ) THEN
                STRLEN       = LEN(TmpBracket)
                CheckBracket = TmpBracket(6:STRLEN)
                REV          = .TRUE.
             ENDIF
          ENDIF

          ! Check if the evaluation of CheckBracket returns true, i.e.
          ! if any of the elements of CheckBracket is enabled. These
          ! can be multiple settings separated by '.or.'.
          ! By default, don't use the content of the bracket
          UseBracket = .FALSE.

          ! Make sure variable ThisBracket is initialized. Needed in the
          ! DO loop below
          ThisBracket = ''

          ! Pack the following into a DO loop to check for multiple
          ! flags separated by '.or.'.
          DO

             ! Leave do loop if ThisBracket is equal to CheckBracket.
             ! In this case, the entire bracket has already been
             ! evaluated.
             IF ( TRIM(CheckBracket) == TRIM(ThisBracket) ) EXIT

             ! Evaluate bracket for '.or.':
             IDX = INDEX(TRIM(CheckBracket),'.or.')

             ! If '.or.' is a substring of the whole bracket, get
             ! substring up to the first '.or.' and write it into variable
             ! ThisBracket, which will be evaluated below. The tail
             ! (everything after the first '.or.') is written into
             ! CheckBracket.
             IF ( IDX > 0 ) THEN
                ThisBracket  = CheckBracket(1:(IDX-1))
                STRLEN       = LEN(CheckBracket)
                CheckBracket = CheckBracket((IDX+4):STRLEN)

             ! If there is no '.or.' in the bracket, simply evaluate the
             ! whole bracket.
             ELSE
                ThisBracket = CheckBracket
             ENDIF

             ! Check if this bracket has been registered as being used.
             ! Scan all extensions, including the core one.
             CALL GetExtOpt( HcoConfig, -999, TRIM(ThisBracket), &
                OptValBool=UseThis, FOUND=FOUND, RC=RC )
             IF ( RC /= HCO_SUCCESS ) RETURN

             ! If bracket name was found in options, update the UseBracket
             ! variable accordingly.
             IF ( FOUND ) THEN
                UseBracket = UseThis

             ! If bracket name was not found, check if this is an extension
             ! name
             ELSE
                ExtNr = GetExtNr( HcoConfig%ExtList, TRIM(ThisBracket) )
                IF ( ExtNr > 0 ) THEN
                   UseBracket = .TRUE.
                ENDIF

             ENDIF

             ! As soon as UseBracket is true, we don't need to evaluate
             ! further
             IF ( UseBracket ) EXIT

          ENDDO

          ! We need to skip the content of this bracket?
          SKIP = .NOT. UseBracket

          ! Eventually reverse the skip flag
          IF ( REV ) THEN
             SKIP = .NOT. SKIP
          ENDIF

          ! If bracket is skipped, adjust skip level accordingly.
          ! This is so that we know when it's time to flip back to
          ! a bracket that is being used (if brackets are nested).
          IF ( SKIP ) THEN
             SKIPLEVEL = NEST
          ENDIF
       ENDIF

       ! Verbose mode
       IF ( verb ) THEN
          MSG = 'Opened shortcut bracket: '//TRIM(TmpBracket)
          CALL HCO_MSG( HcoConfig%Err, MSG )
          WRITE(MSG,*) ' - Skip content of this bracket: ', SKIP
          CALL HCO_MSG( HcoConfig%Err, MSG )
       ENDIF
    ENDIF

    ! Close a bracket
    IF ( STAT == 6 ) THEN

       ! This must be the latest opened bracket
       IF ( TRIM(TmpBracket) /= TRIM(AllBrackets(NEST)) ) THEN
          MSG = 'Closing bracket does not match opening bracket: '// &
             TRIM(TmpBracket)//', expected: '//TRIM(AllBrackets(NEST))
          CALL HCO_ERROR( HcoConfig%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! If that was the latest opened bracket that was disabled
       IF ( SKIPLEVEL == NEST ) THEN
          SKIP = .FALSE.
       ENDIF

       ! Update nesting level
       AllBrackets(NEST) = ''
       NEST              = NEST - 1

       ! Verbose mode
       IF ( verb ) THEN
          MSG = 'Closed shortcut bracket: '//TRIM(TmpBracket)
          CALL HCO_MSG( HcoConfig%Err, MSG )
          WRITE(MSG,*) ' - Skip following lines: ', SKIP
          CALL HCO_MSG( HcoConfig%Err, MSG )
       ENDIF
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE BracketCheck
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: AddShadowFields
!
! !DESCRIPTION: Subroutine AddShadowFields adds a shadow container for every
! additional category of a base emission field. These container contain the
! same container as the 'mother' container, but an additional scale factor
! of zero will be applied to them. This makes sure that no additional emissions
! are created by the virtue of the shadow container.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AddShadowFields( HcoConfig, Lct, Cats, nCat, RC )
!
! !USES:
!
    USE HCO_DATACONT_MOD,  ONLY : CatMax, ZeroScalID
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj), POINTER        :: HcoConfig    ! Config object
    TYPE(ListCont),  POINTER        :: Lct          ! List container of interest
    INTEGER,         INTENT(IN)     :: Cats(CatMax) ! Category numbers
    INTEGER,         INTENT(IN)     :: nCat         ! number of categories
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT)  :: RC           ! Success/failure
!
! !REVISION HISTORY:
!  15 Feb 2015 - C. Keller   - Initial version.
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL                       :: verb
    INTEGER                       :: I, N
    TYPE(ListCont), POINTER       :: Shd
    CHARACTER(LEN=255)            :: MSG
    CHARACTER(LEN=5)              :: C5

    CHARACTER(LEN=255), PARAMETER :: LOC   = 'AddShadowFields (hco_config_mod.F90)'

    !======================================================================
    ! AddShadowFields begins here
    !======================================================================

    ! Nothing to do if ncat is only 1
    IF ( nCat <= 1 ) THEN
       RC = HCO_SUCCESS
       RETURN
    ENDIF

    ! Init
    verb = HCO_IsVerb( HcoConfig%Err, 1 )
    Shd  => NULL()


!    ! Get number of currently used scale factors
!    N = 0
!    DO I = 1, SclMax
!       IF ( Lct%Dct%Scal_cID(I) < 0 ) EXIT
!       N = N + 1
!    ENDDO
!
!    ! There has to be space for scale factor zero.
!    IF ( N >= SclMax ) THEN
!       MSG = 'Cannot add shadow scale factor (zeros) - : ' // &
!             'All scale factors already used: ' // TRIM(Lct%Dct%cName)
!       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
!       RETURN
!    ENDIF

    ! Get number of scale factor IDs. Will add one more scale factor of
    ! zero to this list.
    N = Lct%Dct%nScalID

    ! Create 'shadow' container for every additional category.
    ! Add scale factor zero to it, so that emissions will all be zero.
    DO I = 2, nCat

       ! Create new data container
       CALL ConfigList_AddCont ( Shd, HcoConfig%ConfigList )

       ! Character of category
       Write(C5,'(I5.5)') Cats(I)

       ! Shadow variables. Append category name to name
       Shd%Dct%DctType       = Lct%Dct%DctType
       Shd%Dct%cName         = TRIM(Lct%Dct%cName) // '_Cat' // TRIM(C5)
       Shd%Dct%SpcName       = Lct%Dct%SpcName
       Shd%Dct%Hier          = Lct%Dct%Hier
       Shd%Dct%ExtNr         = Lct%Dct%ExtNr
       Shd%Dct%Cat           = Cats(I)

       ! Pass scale factors, add scale factor of zero to it
       ALLOCATE ( Shd%Dct%Scal_cID(N+1) )
       IF ( N > 0 ) THEN
          Shd%Dct%Scal_cID(1:N) = Lct%Dct%Scal_cID(1:N)
       ENDIF
       Shd%Dct%Scal_cID(N+1) = ZeroScalID
       Shd%Dct%nScalID       = N + 1

       ! Connect to data from main container. Make sure the new container
       ! is not identified as the home container (only points to the file
       ! data container of another data container.
       Shd%Dct%DtaHome =  Shd%Dct%DtaHome - 1
       Shd%Dct%Dta     => Lct%Dct%Dta

       ! verbose mode
       IF ( verb ) THEN
          MSG = 'Created shadow base emission field: ' // TRIM(Shd%Dct%cName)
          CALL HCO_MSG(HcoConfig%Err,MSG)
       ENDIF

       ! Cleanup
       Shd => NULL()
    ENDDO !I

    ! Add zero scale factor container
    CALL AddZeroScal( HcoConfig, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE AddShadowFields
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: AddZeroScal
!
! !DESCRIPTION: Subroutine AddZeroScal adds a scale factor of zero to the
! configuration container list. This scale factor is an internal scale factor
! used in combination with the 'shadow' containers. Its scale factor ID is
! defined in Hco\_DataCont\_Mod and must not be used otherwise, e.g. there
! must not be another scale factor in the HEMCO configuration file with the
! same scale factor ID. Otherwise, HEMCO will create an error lateron.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AddZeroScal( HcoConfig, RC )
!
! !USES:
!
    USE HCO_DATACONT_MOD,  ONLY : ZeroScalID
    USE HCO_DATACONT_MOD,  ONLY : ListCont_Find
    USE HCO_FILEDATA_MOD,  ONLY : FileData_Init
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj), POINTER       :: HcoConfig   ! Config object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT)         :: RC          ! Success/failure
!
! !REVISION HISTORY:
!  15 Feb 2015 - C. Keller   - Initial version.
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ListCont), POINTER       :: Lct
    TYPE(FileData), POINTER       :: Dta
    CHARACTER(LEN=255)            :: MSG

    LOGICAL                       :: FOUND
    CHARACTER(LEN=255), PARAMETER :: LOC   = 'AddZeroScal (hco_config_mod.F90)'

    !======================================================================
    ! AddZeroScal begins here
    !======================================================================

    ! Initialize
    Lct => NULL()
    Dta => NULL()

    ! Check if this container already exists
    CALL ListCont_Find ( HcoConfig%ConfigList, 'DUMMYSCALE_ZERO', FOUND )

    ! Only do on first call
    IF ( .NOT. FOUND ) THEN

       ! Add new container to configuration list and set data container
       ! attributes.
       CALL ConfigList_AddCont ( Lct, HcoConfig%ConfigList )
       Lct%Dct%DctType      = HCO_DCTTYPE_SCAL
       Lct%Dct%cName        = 'DUMMYSCALE_ZERO'
       Lct%Dct%ScalID       = ZeroScalID
       Lct%Dct%Oper         = 1

       ! Create new file data container and fill it with values.
       CALL FileData_Init ( Dta )
       Dta%ncFile    = '0.0'
       Dta%ncPara    = '-'
       Dta%OrigUnit  = 'unitless'
       Dta%CycleFlag = HCO_CFLAG_CYCLE
       Dta%SpaceDim  = 2
       Dta%ncRead    = .FALSE.
       Dta%IsLocTime = .TRUE.

       ! Connect data container
       Lct%Dct%Dta => Dta

       ! verbose mode
       IF ( HCO_IsVerb( HcoConfig%Err, 2 ) ) THEN
          MSG = 'Created a fake scale factor with zeros'
          CALL HCO_MSG(HcoConfig%Err,MSG)
          MSG = 'This field will be used to artificially expand ' // &
                'over multiple emission categories'
          CALL HCO_MSG(HcoConfig%Err,MSG)
       ENDIF

       ! Cleanup
       Lct => NULL()
       Dta => NULL()
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE AddZeroScal
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtSwitch2Buffer
!
! !DESCRIPTION: Subroutine ExtSwitch2Buffer reads the HEMCO extension
! switches and registers all enabled extensions.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtSwitch2Buffer( HcoConfig, IU_HCO, EOF, RC )
!
! !USES:
!
    USE CHARPAK_Mod,        ONLY : STRREPL, STRSPLIT, TRANLC
    USE HCO_EXTLIST_MOD,    ONLY : AddExt, AddExtOpt, HCO_GetOpt
    USE HCO_EXTLIST_MOD,    ONLY : GetExtNr
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj), POINTER       :: HcoConfig   ! Config object
    INTEGER,         INTENT(IN)    :: IU_HCO      ! HEMCO configfile LUN
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,         INTENT(INOUT) :: EOF         ! End of file?
    INTEGER,         INTENT(INOUT) :: RC          ! Success/failure
!
! !REVISION HISTORY:
!  17 Sep 2013 - C. Keller   - Initialization (update)
!  30 Sep 2014 - R. Yantosca - Declare SUBSTR and SPECS w/ 2047 characters,
!                              which lets us handle extra-long species lists
!  21 Apr 2015 - R. Yantosca - Bug fix: now look for END_SECTION before
!                              testing if the line is a comment.  This will
!                              allow for tags labeled "### END SECTION".
!  12 Dec 2015 - C. Keller   - Added argument IgnoreIfExist to AddExtOpt to
!                              make sure that nested configuration files do
!                              use the settings set at highest level.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: I, N, Idx, ExtNr
    LOGICAL               :: Enabled, NewExt
    CHARACTER(LEN=255)    :: LOC
    CHARACTER(LEN=1023)   :: OPTS
    CHARACTER(LEN=2047)   :: LINE
    CHARACTER(LEN=2047)   :: SUBSTR(255), SPECS(255)

    !======================================================================
    ! ExtSwitch2Buffer begins here
    !======================================================================

    ! Enter
    LOC   = 'ExtSwitch2Buffer (hco_config_mod.F90)'
    RC    = HCO_SUCCESS
    ExtNr = -1

    ! Do until exit
    DO

       ! Read line
       CALL HCO_ReadLine ( IU_HCO, LINE, EOF, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Return if EOF
       IF ( EOF ) RETURN

       ! Exit here if end of section encountered.  Place this before the
       ! test for comment to allow for "### END SECTION" tags (bmy, 4/21/15)
       IF ( INDEX ( LINE, 'END SECTION' ) > 0 ) RETURN

       ! Jump to next line if line is commented out
       IF ( LINE(1:1) == HCO_CMT ) CYCLE

       ! Check if these are options
       IF ( INDEX(LINE,'-->') > 0 ) THEN
          ! Only add if extension is defined!
          IF ( ExtNr >= 0 .AND. Enabled ) THEN
             CALL AddExtOpt( HcoConfig, TRIM(LINE), &
                             ExtNr, RC, IgnoreIfExist=.TRUE. )
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF
          CYCLE
       ENDIF

       ! ---------------------------------------------------------------------
       ! If the line is not an extension option, treat it as an extension
       ! definition (e.g. 108 MEGAN : on ISOP/ACET/PRPE/C2H4/ALD2)
       ! ---------------------------------------------------------------------

       ! Split character string
       CALL STRREPL ( LINE, HCO_TAB, HCO_TAB )
       CALL STRSPLIT( LINE, HCO_SPC, SUBSTR, N )

       ! Jump to next line if this line is empty
       IF ( N <= 1 ) CYCLE

       ! Check if extension already exists, e.g. if this is a nested HEMCO configuration
       ! file and the same extension has already been defined. In that case, use the
       ! on/off toggle that has already been defined.
       ExtNr = GetExtNr( HcoConfig%ExtList, TRIM(SUBSTR(2)) )

       ! Three possibilities:
       ! - ExtNr is -999              --> extension does not yet exist
       ! - ExtNr is a positive number --> extension exists and is enabled
       ! - ExtNr is -1                --> extension exists and is disabled
       IF ( ExtNr == -999 ) THEN
          NewExt = .TRUE.
       ELSEIF ( ExtNr >= 0 ) THEN
          NewExt  = .FALSE.
          Enabled = .TRUE.
       ELSE
          NewExt  = .FALSE.
          Enabled = .FALSE.
       ENDIF

       ! The following needs to be done for new extensions only
       IF ( NewExt ) THEN

          ! Check for on-switch. This is either the
          ! 3rd or the 4th substring, depending on the
          ! location of the colon sign!
          IF ( TRIM(SUBSTR(3)) /= ':' ) THEN
             idx = 3
          ELSE
             idx = 4
          ENDIF
          CALL TRANLC( TRIM(SUBSTR(idx)) )
          IF ( TRIM(SUBSTR(idx)) == 'on' ) THEN
             Enabled = .TRUE.
          ELSE
             Enabled = .FALSE.
          ENDIF

          ! Register extension name, number and species
          ! idx is the position of the species names
          idx = idx+1
          READ( SUBSTR(1), * ) ExtNr
          CALL AddExt ( HcoConfig, TRIM(SUBSTR(2)), &
                        ExtNr, Enabled, SUBSTR(idx), RC )
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Register species (specNames)
          IF ( Enabled ) THEN

             CALL STRSPLIT( SUBSTR(idx), &
                     HCO_GetOpt(HcoConfig%ExtList,'Separator'), SPECS, N )
             IF ( N < 1 ) THEN
                CALL HCO_ERROR ( HcoConfig%Err, 'No species defined', RC, THISLOC=LOC )
                RETURN
             ENDIF
             DO I = 1, N
                CALL SpecName_Register ( HcoConfig, SPECS(I), RC )
                IF ( RC /= HCO_SUCCESS ) RETURN
             ENDDO
          ENDIF
       ENDIF ! NextExt
    ENDDO

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ExtSwitch2Buffer
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ReadSettings
!
! !DESCRIPTION: Subroutine ReadSettings reads the HEMCO settings,
! stores them as HEMCO core extension options, and also evaluates
! some of the values (e.g. to initialize the HEMCO error module).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadSettings( HcoConfig, IU_HCO, EOF, RC )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,    ONLY : AddExtOpt, GetExtOpt, CoreNr
    USE HCO_EXTLIST_MOD,    ONLY : HCO_SetDefaultToken
    USE HCO_EXTLIST_MOD,    ONLY : HCO_GetOpt
    USE CHARPAK_MOD,        ONLY : STRREPL, STRSPLIT, TRANLC
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj), POINTER       :: HcoConfig   ! Config obj
    INTEGER,         INTENT(IN)    :: IU_HCO      ! HEMCO configfile LUN
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,         INTENT(INOUT) :: EOF         ! End of file?
    INTEGER,         INTENT(INOUT) :: RC          ! Success/failure
!
! !REVISION HISTORY:
!  17 Sep 2013 - C. Keller   - Initialization (update)
!  21 Apr 2015 - R. Yantosca - Bug fix: now look for END_SECTION before
!                              testing if the line is a comment.  This will
!                              allow for tags labeled "### END SECTION".
!  12 Dec 2015 - C. Keller   - Added argument IgnoreIfExist to AddExtOpt to
!                              make sure that nested configuration files do
!                              use the settings set at highest level.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL               :: FOUND
    INTEGER               :: I, N, POS
    INTEGER               :: verb
    INTEGER               :: warn

    ! Strings
    CHARACTER(LEN=255)    :: Line
    CHARACTER(LEN=255)    :: Loc
    CHARACTER(LEN=255)    :: Msg
    CHARACTER(LEN=255)    :: LogFile
    CHARACTER(LEN=255)    :: DiagnPrefix
    CHARACTER(LEN=255)    :: MetField
    CHARACTER(LEN=255)    :: GridRes

    !======================================================================
    ! ReadSettings begins here
    !======================================================================

    ! Enter
    Loc = 'ReadSettings (hco_config_mod.F90)'

    !-----------------------------------------------------------------------
    ! Read settings and add them as options to core extensions
    !-----------------------------------------------------------------------

    ! Do until exit
    DO

       ! Read line
       CALL HCO_ReadLine ( IU_HCO, LINE, EOF, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Return if EOF
       IF ( EOF ) EXIT

       ! Exit here if end of section encountered
       IF ( INDEX ( LINE, 'END SECTION' ) > 0 ) EXIT

       ! Jump to next line if line is commented out
       IF ( LINE(1:1) == HCO_CMT ) CYCLE

       ! Ignore empty lines
       IF ( TRIM(LINE) == '' ) CYCLE

       ! Add this option to HEMCO core
       CALL AddExtOpt ( HcoConfig, TRIM(LINE), &
                        CoreNr, RC, IgnoreIfExist=.TRUE. )
       IF ( RC /= HCO_SUCCESS ) RETURN

    ENDDO

    !-----------------------------------------------------------------------
    ! Initialize error object if needed.
    ! Extract values to initialize error module and set some further
    ! HEMCO variables. Only the first time the settings are read (settings
    ! can be read multiple times if nested HEMCO configuration files are
    ! used)
    !-----------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED(HcoConfig%Err) ) THEN

       ! Verbose mode?
       CALL GetExtOpt( HcoConfig, CoreNr, 'Verbose', &
                       OptValInt=verb, FOUND=FOUND, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( .NOT. FOUND ) THEN
          verb = 3
          WRITE(*,*) 'Setting `Verbose` not found in HEMCO logfile - use 3'
       ENDIF

       ! Logfile to write into
       CALL GetExtOpt( HcoConfig, CoreNr, 'Logfile', &
                       OptValChar=Logfile, FOUND=FOUND, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( .NOT. FOUND ) THEN
          LogFile = 'HEMCO.log'
          WRITE(*,*) 'Setting `Logfile` not found in HEMCO logfile - use `HEMCO.log`'
       ENDIF

       ! Prompt warnings to logfile?
       CALL GetExtOpt( HcoConfig, CoreNr, 'Warnings', &
                       OptValInt=warn, FOUND=FOUND, RC=RC  )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( .NOT. FOUND ) THEN
          warn = 3
          WRITE(*,*) 'Setting `Warnings` not found in HEMCO logfile - use 3'
       ENDIF

       ! Initialize (standard) HEMCO tokens
       CALL HCO_SetDefaultToken( HcoConfig, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! If LogFile is equal to wildcard character, set LogFile to asterik
       ! character. This will ensure that all output is written to standard
       ! output!
       IF ( TRIM(LogFile) == HCO_GetOpt(HcoConfig%ExtList,'Wildcard') ) LogFile = '*'

       ! We should now have everything to define the HEMCO error settings
       CALL HCO_ERROR_SET( HcoConfig%amIRoot, HcoConfig%Err, LogFile, &
                           verb, warn, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

    ENDIF

#ifndef MODEL_GEOS
#ifndef MODEL_WRF
#ifndef MODEL_CESM
#ifndef ESMF_
    !=======================================================================
    ! Look for met field and grid resolution.  When running the HEMCO
    ! standalone these will need to be read from the configuration file.
    ! Otherwise, HEMCO will inherit the met field and grid resolution
    ! of the parent model (GC-Classic, GCHP, etc.)
    !
    ! NOTE: Only do this check if not using GEOS-Chem in an external ESM!
    !=======================================================================

    ! Look for met field
    CALL GetExtOpt( HcoConfig, CoreNr, 'MET', &
                    OptValChar=MetField, FOUND=FOUND, RC=RC )
    IF ( FOUND ) THEN
       HcoConfig%MetField = TRIM( MetField )
    ENDIF

    ! Look for grid resolution
    ! Make sure resolution string is in the proper FlexGrid format
    CALL GetExtOpt( HcoConfig, CoreNr, 'RES', &
                    OptValChar=GridRes, FOUND=FOUND, RC=RC )
    IF ( FOUND ) THEN
       SELECT CASE( TRIM( GridRes ) )
          CASE( '4x5' )
             GridRes = '4.0x5.0'
          CASE( '2x25', '2x2.5' )
             GridRes = '2.0x2.5'
          CASE( '05x0625' )
             GridRes = '0.5x0.625'
          CASE( '025x03125' )
             GridRes = '0.25x0.3125'
          CASE DEFAULT
             Msg = 'Improperly formatted grid resolution: ' // TRIM( GridRes )
             CALL HCO_Error( HcoConfig%Err, Msg, RC, Loc )
             RETURN
       END SELECT
       HcoConfig%GridRes = TRIM( GridRes )
    ENDIF
#endif
#endif
#endif
#endif

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ReadSettings
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RegisterPrepare
!
! !DESCRIPTION: Subroutine RegisterPrepare extracts the spatial
! coverages of all mask fields as well as the HEMCO species IDs of
! all base emissions.
!\\
!\\
! The species IDs are determined by matching the species name read
! from the configuration file (in ConfigList) and the species names
! defined in the HEMCO state object HcoState.
!\\
!\\
! Mask coverages are defined based upon the passed horizontal grid
! extensions on this CPU (xrng and yrng).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE RegisterPrepare( HcoState, RC )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,  ONLY : ExtNrInUse
    USE HCO_STATE_Mod,    ONLY : HCO_GetHcoID
    USE HCO_DATACONT_MOD, ONLY : ListCont_NextCont
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState   ! HEMCO state obj.
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  18 Sep 2013 - C. Keller - Initial version (update)
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ListCont), POINTER      :: Lct
    INTEGER                      :: ThisCover, ThisHcoID, FLAG
    INTEGER                      :: lon1, lon2, lat1, lat2
    INTEGER                      :: cpux1, cpux2, cpuy1, cpuy2
    CHARACTER(LEN=255)           :: MSG

    !=================================================================
    ! RegisterPrepare begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER ( HcoState%Config%Err, 'RegisterPrepare', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Initialize
    Lct => NULL()

    ! Grid boundaries on this CPU. Will be needed to calculate
    ! coverages.
    ! NOTE: Use midpoints here because only those become defined in
    ! the ESMF environment (xedge and yedge are not used anywhere
    ! else in ESMF!).
    cpux1 = FLOOR(MINVAL(HcoState%Grid%XMID%Val))
    cpux2 = FLOOR(MAXVAL(HcoState%Grid%XMID%Val))
    cpuy1 = CEILING(MINVAL(HcoState%Grid%YMID%Val))
    cpuy2 = CEILING(MAXVAL(HcoState%Grid%YMID%Val))

    ! Make sure values are within -180.0 to 180.0
    IF ( cpux1 >= 180 ) cpux1 = cpux1 - 360
    IF ( cpux2 >= 180 ) cpux2 = cpux2 - 360

    ! verbose
    IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
       WRITE(MSG,*) 'Start to prepare fields for registering!'
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) 'This CPU x-range: ', cpux1, cpux2
       CALL HCO_MSG(HcoState%Config%Err,MSG)
       WRITE(MSG,*) 'This CPU y-range: ', cpuy1, cpuy2
       CALL HCO_MSG(HcoState%Config%Err,MSG)
    ENDIF

    ! Get next (first) line of ConfigList
    CALL ListCont_NextCont ( HcoState%Config%ConfigList, Lct, FLAG )

    ! Loop over all lines
    DO WHILE ( FLAG == HCO_SUCCESS )

       ! Check if data container defined
       IF ( .NOT. ASSOCIATED(Lct%Dct) ) THEN
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, Lct, FLAG )
          CYCLE
       ENDIF

       ! verbose
       IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
          WRITE(MSG,*) 'Prepare ', TRIM(Lct%Dct%cName)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       ! For base fields or data fields used in one of the HEMCO
       ! extensions:
       IF ( Lct%Dct%DctType == HCO_DCTTYPE_BASE ) THEN

          ! Only do for entries that will be used!
          IF ( ExtNrInUse( HcoState%Config%ExtList, Lct%Dct%ExtNr ) ) THEN

             ! Extract HEMCO species ID. This will return -1 for
             ! undefined species and 0 for wildcard character.
             ThisHcoID     = HCO_GetHcoID( Lct%Dct%SpcName, HcoState )
             Lct%Dct%HcoID = ThisHcoID

             ! verbose
             IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
                WRITE(MSG,*) 'Assigned HEMCO species ID: ', Lct%Dct%HcoID
                CALL HCO_MSG(HcoState%Config%Err,MSG)
             ENDIF

          ! Else: assign default value. These containers will be
          ! removed in the next step!
          ELSE
             Lct%Dct%HcoID = -999
          ENDIF

       ! Calculate coverage for masks
       ELSE IF ( Lct%Dct%DctType   == HCO_DCTTYPE_MASK .AND. &
                 Lct%Dct%Dta%Cover == -999                   ) THEN

          If (HcoState%Options%isESMF) Then
             ThisCover = -1
          Else
             ! Get mask edges
             lon1 = Lct%Dct%Dta%ncYrs(1)
             lat1 = Lct%Dct%Dta%ncYrs(2)
             lon2 = Lct%Dct%Dta%ncMts(1)
             lat2 = Lct%Dct%Dta%ncMts(2)

             ThisCover = CALC_COVERAGE( lon1,  lon2,  &
                                        lat1,  lat2,  &
                                        cpux1, cpux2, &
                                        cpuy1, cpuy2   )

             ! There appear to be some issues with full masks coverages
             ! when working in an MPI environment. Specifically, masks
             ! can be seen as fully covering a given CPU even though in
             ! reality it may only cover parts of it. Thus, in ESMF mode
             ! always set coverage to zero or partial (ckeller, 3/17/16).
             IF ( HcoState%Options%isESMF ) THEN
                IF ( ThisCover == 1 ) ThisCover = -1
             ENDIF
          ENDIF

          ! Update container information
          Lct%Dct%Dta%Cover    = ThisCover
          Lct%Dct%Dta%ncYrs(:) = -999
          Lct%Dct%Dta%ncMts(:) = -999

          IF ( HCO_IsVerb(HcoSTate%Config%Err,3) ) THEN
             WRITE(MSG,*) 'Coverage: ', Lct%Dct%Dta%Cover
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF
       ENDIF

       ! Advance to next line
       CALL ListCont_NextCont ( HcoState%Config%ConfigList, Lct, FLAG )
    ENDDO

    ! Cleanup
    Lct => NULL()

    ! Return w/ success
    CALL HCO_LEAVE ( HcoState%Config%Err, RC )

  END SUBROUTINE RegisterPrepare
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_Base
!
! !DESCRIPTION: Subroutine Register\_Base registers all base emission
! data and writes out all associated scale factor IDs.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_Base( HcoState, RC )
!
! !USES:
!
    USE HCO_READLIST_Mod,      ONLY : ReadList_Set
    USE HCO_DATACONT_Mod,      ONLY : DataCont_Cleanup
    USE HCO_DATACONT_MOD,      ONLY : ListCont_NextCont
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState    ! HEMCO state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  18 Jun 2013 - C. Keller: Initialization
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ListCont), POINTER   :: Lct

    ! Scalars
    INTEGER               :: N, cID, HcoID
    INTEGER               :: targetID, FLAG
    LOGICAL               :: Ignore, Add
    CHARACTER(LEN=255)    :: MSG

    !======================================================================
    ! Register_Base begins here
    !======================================================================

    ! Enter
    CALL HCO_ENTER ( HcoState%Config%Err, 'Register_Base (hco_config_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       PRINT *,'Error in HCO_ENTER called from Register_Base'
       RETURN
    ENDIF

    ! Initialize
    Lct => NULL()

    ! Point to next (first) line in ConfigList
    CALL ListCont_NextCont ( HcoState%Config%ConfigList, Lct, FLAG )

    ! Loop over temporary arrays
    DO WHILE ( FLAG == HCO_SUCCESS )

       ! Reset ignore flag
       Ignore = .FALSE.

       ! Skip entry if data container not defined
       IF ( .NOT. ASSOCIATED(Lct%Dct) ) THEN
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, Lct, FLAG )
          CYCLE
       ENDIF

       ! Skip entry if it's not a base field
       IF ( (Lct%Dct%DctType /= HCO_DCTTYPE_BASE) ) THEN
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, Lct, FLAG )
          CYCLE
       ENDIF

       ! If this base field is not used (either because it belongs to
       ! an extension that is not enabled or because its HEMCO or
       ! model species ID is undefined), we don't need this container
       ! any more. Hence remove it.
       ! Note: Routine RegisterPrepare assigns negative HcoID's to all
       ! base fields with invalid ExtNr's, so it is ok to check only
       ! for HcoID here. If data is used in one of the HEMCO extensions
       ! and has a species flag of '*' (= always read), its species ID
       ! becomes set to 0 in RegisterPrepare.
       HcoID = Lct%Dct%HcoID
       IF ( HcoID < 0 ) THEN
          Ignore = .TRUE.
       ELSE IF ( HcoID > 0 ) THEN
          IF ( HcoState%Spc(HcoID)%ModID < 0 ) Ignore = .TRUE.
       ENDIF

       IF ( Ignore ) THEN
          IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
             WRITE(MSG,*) &
                  'Register_Base: Ignore (and remove) base field ', &
                  TRIM(Lct%Dct%cName)
             CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1='-')
          ENDIF

          ! Remove data container from list.
          CALL DataCont_Cleanup ( Lct%Dct )
          Lct%Dct => NULL()
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, Lct, FLAG )
          CYCLE
       ENDIF

       ! Verbose mode
       IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
          WRITE(MSG,*) 'Register_Base: Checking ', TRIM(Lct%Dct%cName)
          CALL HCO_MSG(HcoState%Config%Err,MSG,SEP1='-')
       ENDIF

       ! -------------------------------------------------------------
       ! Extract vector of scale factor container IDs to be applied
       ! to this base field (vector Scal_cID). For now, this container
       ! contains the scale factor IDs, hence need to convert to
       ! container IDs. Beforehand, add scale factor IDs to internal
       ! list of used scale factors (UnqScalIDs).
       CALL ScalID_Register ( Lct%Dct, HcoState%Config, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *,'Error in ScaleID_Register called from Register_Base'
          RETURN
       ENDIF

       ! Get target ID of this container. The targetID corresponds
       ! to the container ID cID into which emissions data of the
       ! current container (Lct) will be added to. Typically,
       ! targetID is equal to cID, i.e. a container holds the data
       ! array corresponding to its source data information. In
       ! cases where multiple base emissions have same properties,
       ! however, we can merge these fields prior to emission
       ! calculation to save some calculation operations.
       ! Requirement is that these base emissions have same species
       ! ID, emission category and hierarchy, ext. number, scale
       ! factors, and update frequency.
       CALL Get_targetID( HcoState, Lct, targetID, RC)
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *,'Error in Get_targetID called from Register_Base'
          RETURN
       ENDIF

       ! verbose
       IF ( HCO_IsVerb(HcoState%Config%Err,3) ) THEN
          WRITE(MSG,*) 'Container ID     : ', Lct%Dct%cID
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(MSG,*) 'Assigned targetID: ', targetID
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       ! Negative targetID is assigned to base data that doesn't need
       ! to be considered either because it's a regional inventory
       ! with no spatial overlap with the region domain on this CPU,
       ! or because there exist another inventory with higher
       ! priority (and same category) that will overwrite these
       ! emissions data anyway!
       IF ( targetID <= 0 ) THEN
          CALL DataCont_Cleanup ( Lct%Dct )
          Lct%Dct => NULL()
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, Lct, FLAG )
          CYCLE
       ENDIF

       ! Pass targetID to container
       Lct%Dct%targetID = targetID

       ! Register container in ReadList. Containers will be listed
       ! in the reading lists sorted by cID.
       CALL ReadList_Set( HcoState, Lct%Dct, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *,'Error in ReadList_Set called from Register_Base'
          RETURN
       ENDIF

       ! Print some information if verbose mode is on
       IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
          WRITE(MSG,*) 'Base field registered: ', TRIM(Lct%Dct%cName)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       ! Advance to next line
       CALL ListCont_NextCont ( HcoState%Config%ConfigList, Lct, FLAG )
    ENDDO

    ! Cleanup
    Lct => NULL()

    ! Return w/ success
    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE Register_Base
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_Scal
!
! !DESCRIPTION: Subroutine Register\_Scal registers all scale factors.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_Scal( HcoState, RC )
!
! !USES:
!
    USE HCO_ReadList_Mod,      ONLY : ReadList_Set
    USE HCO_DATACONT_MOD,      ONLY : ListCont_NextCont

!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER       :: HcoState    ! HEMCO state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  18 Jun 2013 - C. Keller - Initialization
!  29 Dec 2014 - C. Keller - Now check for masks assigned to scale factors.
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ListCont),   POINTER :: Lct
    TYPE(ScalIDCont), POINTER :: TmpScalIDCont

    ! Scalars
    INTEGER                   :: cID, FLAG
    CHARACTER(LEN=255)        :: MSG
    CHARACTER(LEN=  5)        :: strID
    INTEGER                   :: ThisScalID

    !======================================================================
    ! Register_Scal begins here
    !======================================================================

    ! Enter
    CALL HCO_ENTER ( HcoState%Config%Err, 'Register_Scal (hco_config_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Loop over all scale factor ids
    Lct           => NULL()
    TmpScalIDCont => HcoState%Config%ScalIDList
    DO WHILE ( ASSOCIATED( TmpScalIDCont ) )

       ! Extract this scale factor ID
       ThisScalID = TmpScalIDCont%ScalID

       ! Make ThisLine point to first element of ConfigList
       Lct => NULL()
       CALL ListCont_NextCont ( HcoState%Config%ConfigList, Lct, FLAG )

       ! Loop over all lines in input file and find the one with the
       ! correct scale factor ID
       DO WHILE ( FLAG == HCO_SUCCESS )

          ! Leave if this is the wanted container.
          IF ( ASSOCIATED(Lct%Dct)) THEN
             IF ( Lct%Dct%ScalID == ThisScalID ) EXIT
          ENDIF

          ! Advance to next line otherwise
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, Lct, FLAG )
       ENDDO

       ! Return error if scale factor ID not found
       IF ( .NOT. ASSOCIATED(Lct) ) THEN
          WRITE ( strID, * ) ThisScalID
          MSG = 'Container ID not found: ' // strID
          CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC)
          RETURN
       ENDIF

       ! Return w/ error if container not scale factor or mask
       IF ( Lct%Dct%DctType == HCO_DCTTYPE_BASE ) THEN
          WRITE ( strID, * ) ThisScalID
          MSG = 'Container ID belongs to base field: ' // strID
          CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC)
          RETURN
       ENDIF

       ! Check if this scale factor has a mask field assigned to
       ! it, in which case we have to make sure that the mask field
       ! becomes registered in the scale factor list ScalIDList.
       ! We can do this while evaluating ScalIDList due to the dynamic
       ! structure of the linked list with new containers simply being
       ! added to the end of the list.
       IF ( Lct%Dct%nScalID > 0 ) THEN
          CALL ScalID_Register ( Lct%Dct, HcoState%Config, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       ! Register container in ReadList. Containers will be listed
       ! in the reading lists sorted by cID.
       CALL ReadList_Set( HcoState, Lct%Dct, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Print some information if verbose mode is on
       IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
          WRITE(MSG,*) 'Scale field registered: ', TRIM(Lct%Dct%cName)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF

       ! Advance
       TmpScalIDCont => TmpScalIDCont%NEXT

    ENDDO

    ! Cleanup
    Lct           => NULL()
    TmpScalIDCont => NULL()

    ! Return w/ success
    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE Register_Scal
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_targetID
!
! !DESCRIPTION: Subroutine Get\_targetID returns the target ID of a
! container. The target ID can point to the container ID (cID) of
! another base field if multiple emissions shall be added together
! prior to emission calculation, e.g. sectoral emissions data with
! same species ID, category, hierarchy, extension number, scale factors,
! etc.
!\\
!\\
! Target ID is set to -999 if there exists another inventory over
! the full spatial region covered by this CPU for this species but
! with higher hierarchy. In this case, we can ignore the current
! container from here onwards!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_targetID( HcoState, Lct, targetID, RC )
!
! !USES:
!
    USE HCO_DataCont_Mod, ONLY : ListCont_Find
    USE HCO_DataCont_Mod, ONLY : ListCont_NextCont
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState
    TYPE(ListCont),  POINTER       :: Lct
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(  OUT) :: targetID
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(INOUT) :: RC
!
! !NOTE: If data from multiple containers are added, the target ID
! is always set to the lowest cID of all involved containers, i.e.
! data are added to the container with the lowest cID. This makes
! sure that data is not accidentally overwritten, e.g. when updating
! container contents!
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller - Initialization
!  07 Dec 2015 - C. Keller - Make sure emissions with limited time range do
!                            never erase lower hierarchy base emissions.
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ListCont), POINTER   :: tmpLct
    TYPE(ListCont), POINTER   :: mskLct

    ! Scalars
    INTEGER                   :: HcoID, Cat, Hier, Scal, ExtNr, cID
    INTEGER                   :: tmpID
    INTEGER                   :: I, J, FLAG1, tmpCov
    LOGICAL                   :: found, sameCont
    CHARACTER(LEN=255)        :: MSG
    CHARACTER(LEN=  7)        :: strID

    !======================================================================
    ! Get_targetID begins here
    !======================================================================

    ! Enter
    CALL HCO_ENTER ( HcoState%Config%Err, 'Get_targetID (hco_config_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Initialize
    tmpLct => NULL()
    mskLct => NULL()

    ! Get Tracer ID, category and hierarchy of entry to be checked
    cID   = Lct%Dct%cID
    ExtNr = Lct%Dct%ExtNr
    Cat   = Lct%Dct%Cat
    Hier  = Lct%Dct%Hier
    HcoID = Lct%Dct%HcoID

    ! By default, set target ID to container ID
    targetID = cID

    ! If ExtNr is -999, always read this field. ExtNr becomes zero
    ! if the extension number entry in the configuration file is the
    ! wildcard character
    IF ( ExtNr == -999 ) THEN
       CALL HCO_LEAVE( HcoState%Config%Err, RC )
       RETURN
    ENDIF

    ! If species ID is zero, always read this field as is, i.e. don't
    ! skip it and don't add it to another field!
    ! Species ID become zero if the species ID entry in the
    ! configuration file is the wildcard character.
    IF ( HcoID == 0 ) THEN
       CALL HCO_LEAVE( HcoState%Config%Err, RC )
       RETURN
    ENDIF

    ! Check all scale factors of the current container to see if one
    ! of them is a mask that has no valid entries over the domain of
    ! this CPU. In this case we don't have to consider this field at
    ! all!
    IF ( Lct%Dct%nScalID > 0 ) THEN
       DO I = 1, Lct%Dct%nScalID

          ! Check if it's a valid scale factor
          IF ( Lct%Dct%Scal_cID(I) < 0 ) CYCLE

          ! Find container with this container ID
          ! Note: this should always look up the container ID, but make
          ! check for safety's sake.
          tmpID = Lct%Dct%Scal_cID(I)
          IF ( .NOT. Lct%Dct%Scal_cID_set ) THEN
             CALL ListCont_Find ( HcoState%Config%ConfigList, &
                                  tmpID, 1, FOUND, mskLct )
          ELSE
             CALL ListCont_Find ( HcoState%Config%ConfigList, &
                                  tmpID, 0, FOUND, mskLct )
          ENDIF

          ! Error if scale factor not found
          IF ( .NOT. FOUND ) THEN
             WRITE ( strID, * ) Lct%Dct%Scal_cID(I)
             MSG = 'No scale factor with cID: ' // TRIM(strID)
             CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC)
             RETURN
          ENDIF

          ! Check if this is a mask with zero coverage over this CPU, in
          ! which case we don't need to consider the base field at all!
          IF ( (mskLct%Dct%DctType  == HCO_DCTTYPE_MASK ) .AND. &
               (mskLct%Dct%Dta%Cover == 0 )        ) THEN
             targetID = -999
             IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
                WRITE(MSG,*) 'Data not defined over this CPU, skip ' // &
                     TRIM(Lct%Dct%cName)
                CALL HCO_MSG(HcoState%Config%Err,MSG)
             ENDIF

             ! Return
             CALL HCO_LEAVE( HcoState%Config%Err, RC )
             RETURN
          ENDIF
       ENDDO ! I
    ENDIF

    ! Now find out if there is another base field for the same species,
    ! emission category and extension number, but higher hierarchy.
    ! Such a field also needs to have full coverage over this CPU,
    ! then we can ignore the current container.

    ! Initialize looping pointer
    tmpLct => NULL()
    CALL ListCont_NextCont ( HcoState%Config%ConfigList, tmpLct, FLAG1 )

    ! Loop over containers
    DO WHILE ( FLAG1 == HCO_SUCCESS )

       ! Advance to next container if data container not defined
       IF ( .NOT. ASSOCIATED(tmpLct%Dct) ) THEN
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, tmpLct, FLAG1 )
          CYCLE
       ENDIF

       ! Advance to next container if this is the current container
       IF ( tmpLct%Dct%cID == cID ) THEN
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, tmpLct, FLAG1 )
          CYCLE
       ENDIF

       ! Advance to next container if this is not a base field
       IF ( tmpLct%Dct%DctType /= HCO_DCTTYPE_BASE ) THEN
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, tmpLct, FLAG1 )
          CYCLE
       ENDIF

       ! Advance to next container if not the same extension nr
       IF ( tmpLct%Dct%ExtNr /= ExtNr ) THEN
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, tmpLct, FLAG1 )
          CYCLE
       ENDIF

       ! Advance to next container if not the same species
       IF ( tmpLct%Dct%HcoID /= HcoID ) THEN
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, tmpLct, FLAG1 )
          CYCLE
       ENDIF

       ! Advance to next container if not the same category
       IF ( tmpLct%Dct%Cat /= Cat ) THEN
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, tmpLct, FLAG1 )
          CYCLE
       ENDIF

       ! Advance to next container if lower hierarchy
       IF ( tmpLct%Dct%Hier < Hier ) THEN
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, tmpLct, FLAG1 )
          CYCLE
       ENDIF

       ! Advance to next container if this container has limited time
       ! coverage. Emissions with limited time coverage may not be used
       ! during all of the simulation time, so it's important to keep the
       ! lower hierarchy emission fields in memory in case that those need
       ! to be used instead (e.g. if EDGAR shall only be used between years
       ! 2005 and 2013, we should keep GEIA in case that we are outside of
       ! that time window).
       IF ( ( tmpLct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGE    ) .OR. &
            ( tmpLct%Dct%Dta%CycleFlag == HCO_CFLAG_EXACT    ) .OR. &
            ( tmpLct%Dct%Dta%CycleFlag == HCO_CFLAG_RANGEAVG )      ) THEN
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, tmpLct, FLAG1 )
          CYCLE
       ENDIF

       ! Check for coverage of tmpLct. Default = full coverage (1)
       tmpCov = 1

       ! Check all scale factors of tmpLct to see if this base
       ! field has full coverage over this CPU domain or not.
       IF ( tmpLct%Dct%nScalID > 0 ) THEN
          DO I = 1, tmpLct%Dct%nScalID

             ! Check if it's a valid scale factor
             IF ( tmpLct%Dct%Scal_cID(I) < 0 ) CYCLE

             tmpID = tmpLct%Dct%Scal_cID(I)
             IF ( .NOT. tmpLct%Dct%Scal_cID_set ) THEN
                CALL ListCont_Find ( HcoState%Config%ConfigList, &
                                     tmpID, 1, FOUND, mskLct )
             ELSE
                CALL ListCont_Find ( HcoState%Config%ConfigList, &
                                     tmpID, 0, FOUND, mskLct )
             ENDIF

             ! Error if container not found
             IF ( .NOT. FOUND ) THEN
                WRITE(MSG,*) 'No scale factor with ID: ', tmpID
                CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC)
                RETURN
             ENDIF

             ! Write out coverage.
             ! Note: If one mask has only partial coverage, retain that
             ! value! If we encounter a mask with no coverage, set coverage
             ! to zero and leave immediately.
             IF ( (mskLct%Dct%DctType == HCO_DCTTYPE_MASK) ) THEN
                IF ( mskLct%Dct%Dta%Cover == -1 ) THEN
                   tmpCov = -1
                ELSEIF ( mskLct%Dct%Dta%Cover == 0 ) THEN
                   tmpCov = 0
                   EXIT
                ENDIF
             ENDIF
          ENDDO ! I
       ENDIF

       ! If tmpLct has no coverage, we can ignore this tmpLct as
       ! it will never overwrite data of currCont
       IF ( tmpCov == 0 ) THEN
          CALL ListCont_NextCont ( HcoState%Config%ConfigList, tmpLct, FLAG1 )
          CYCLE
       ENDIF

       ! If we made it up to here and tmpLct has full coverage, then
       ! tmpLct has the same species ID, category, ext. nr.,
       ! and a higher (or the same) hierarchy as Lct.

       ! If hierarchy of tmpLct is higher than Lct and this
         ! container has total coverage over this CPU, it will always
         ! replace all values of Lct. Hence, set targetID to -999
         ! (= ignore container) and return here.
       IF ( (tmpLct%Dct%Hier > Hier) .AND. (tmpCov==1) ) THEN
          IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
             WRITE(MSG,*) 'Skip container ', TRIM(Lct%Dct%cName), &
                          ' because of ', TRIM(tmpLct%Dct%cName)
             CALL HCO_MSG(HcoState%Config%Err,MSG)
          ENDIF

          ! Return
          targetID = -999
          CALL HCO_LEAVE( HcoState%Config%Err, RC )
          RETURN
       ENDIF

       ! If currCont and tmpLct have same hierarchy, scale factors,
       ! and update frequencies, we may add the two fields together
       ! in order to make emission calculation more efficient.
       ! Thus, set target ID to the lower of the two container IDs.
       ! This procedure will ensure that emission data are added
       ! when registering/updating containers in the emissions list
       ! (EmisList). Since containers are sorted in ReadList with
       ! increasing cID, pick the lowest cID to make sure that all
       ! fields are properly added.
       ! Note: this option is currently disabled for ESMF applications.
       IF ( tmpLct%Dct%Hier == Hier .AND. .NOT. HcoState%Options%isESMF ) THEN

          ! temporary flag
          sameCont = .TRUE.

          ! Check for same scale factors
          IF ( tmpLct%Dct%nScalID /= Lct%Dct%nScalID ) THEN
             sameCont = .FALSE.
          ELSE
             DO I = 1, tmpLct%Dct%nScalID
                IF (  tmpLct%Dct%Scal_cID(I) /= &
                     Lct%Dct%Scal_cID(I)     ) THEN
                   sameCont = .FALSE.
                   EXIT
                ENDIF
             ENDDO
          ENDIF

          ! Check for same update frequencies
          IF ( sameCont ) THEN
             IF    (tmpLct%Dct%Dta%ncYrs(1)/=Lct%Dct%Dta%ncYrs(1)) THEN
                sameCont = .FALSE.
             ELSEIF(tmpLct%Dct%Dta%ncYrs(2)/=Lct%Dct%Dta%ncYrs(2)) THEN
                sameCont = .FALSE.
             ELSEIF(tmpLct%Dct%Dta%ncMts(1)/=Lct%Dct%Dta%ncMts(1)) THEN
                sameCont = .FALSE.
             ELSEIF(tmpLct%Dct%Dta%ncMts(2)/=Lct%Dct%Dta%ncMts(2)) THEN
                sameCont = .FALSE.
             ELSEIF(tmpLct%Dct%Dta%ncDys(1)/=Lct%Dct%Dta%ncDys(1)) THEN
                sameCont = .FALSE.
             ELSEIF(tmpLct%Dct%Dta%ncDys(2)/=Lct%Dct%Dta%ncDys(2)) THEN
                sameCont = .FALSE.
             ELSEIF(tmpLct%Dct%Dta%ncHrs(1)/=Lct%Dct%Dta%ncHrs(1)) THEN
                sameCont = .FALSE.
             ELSEIF(tmpLct%Dct%Dta%ncHrs(2)/=Lct%Dct%Dta%ncHrs(2)) THEN
                sameCont = .FALSE.
             ENDIF
          ENDIF

          ! Check for same emitted level
          IF ( sameCont ) THEN
             IF ( ( tmpLct%Dct%Dta%SpaceDim /= Lct%Dct%Dta%SpaceDim ) .OR. &
                  ( tmpLct%Dct%Dta%Levels   /= Lct%Dct%Dta%Levels   ) .OR. &
                  ( tmpLct%Dct%Dta%EmisL1   /= Lct%Dct%Dta%EmisL1   ) .OR. &
                  ( tmpLct%Dct%Dta%EmisL2   /= Lct%Dct%Dta%EmisL2   ) ) THEN
                sameCont = .FALSE.
             ENDIF
          ENDIF

          ! Finally, check for "same" container names. This checks the
          ! container names ignoring the name 'tags'.
          IF ( sameCont ) THEN
             sameCont = Check_ContNames( tmpLct, Lct )
          ENDIF

          ! If "same" containers, set target ID to container ID of
          ! tmpLct if this value is lower than current target ID.
          IF ( sameCont ) THEN
             targetID = MIN( targetID, tmpLct%Dct%cID )
          ENDIF

       ENDIF

       ! Advance to next line
       ! Don't return here, because it is still possible that there is
       ! another inventory in the list coming up which overwrites this
       ! inventory (or another field emissions shall be added to which
       ! has lower container ID and hence needs to be the target
       ! container!).
       CALL ListCont_NextCont ( HcoState%Config%ConfigList, tmpLct, FLAG1 )

    ENDDO !Loop over all entries in ConfigList (tmpLct)

    ! Free pointers
    tmpLct => NULL()
    mskLct => NULL()

    ! Leave w/ success
    CALL HCO_LEAVE( HcoState%Config%Err, RC )

  END SUBROUTINE Get_targetID
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Coverage
!
! !DESCRIPTION: Function Calc\_Coverage calculates the coverage of
! the specified lon/lat box with the area covered by the inventory.
! Returns 0 if no overlap, 1 if complete overlap, and -1 for partial
! overlap.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Calc_Coverage( msk_x1, msk_x2, msk_y1, msk_y2,  &
                          cpu_x1, cpu_x2, cpu_y1, cpu_y2 ) RESULT ( COVERAGE )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: msk_x1
    INTEGER, INTENT(IN) :: msk_x2
    INTEGER, INTENT(IN) :: msk_y1
    INTEGER, INTENT(IN) :: msk_y2
    INTEGER, INTENT(IN) :: cpu_x1
    INTEGER, INTENT(IN) :: cpu_x2
    INTEGER, INTENT(IN) :: cpu_y1
    INTEGER, INTENT(IN) :: cpu_y2
!
! !RETURN VALUE:
!
    INTEGER             :: COVERAGE
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! CALC_COVERAGE begins here
    !======================================================================

    ! Check if specified area does not overlap with inventory
    COVERAGE = 1
    IF ( (msk_x1 > cpu_x2) .OR. (msk_x2 < cpu_x1) .OR. &
         (msk_y1 > cpu_y2) .OR. (msk_y2 < cpu_y1)        ) THEN
       COVERAGE = 0
       RETURN
    ENDIF

    ! Check for partial coverage
    IF ( (msk_x1 > cpu_x1) .OR. (msk_x2 < cpu_x2) .OR. &
         (msk_y1 > cpu_y1) .OR. (msk_y2 < cpu_y2)        ) THEN
       COVERAGE = -1
    ENDIF

  END FUNCTION Calc_Coverage
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadAndSplit_Line
!
! !DESCRIPTION: Subroutine ReadAndSplit\_Line reads a line from the HEMCO
! config file and parses the specified columns into the passed integer
! and character variables.  If the optional argument inLine is provided,
! this line will be parsed, otherwise a new line will be read from the config
! file. If the optional argument outLine is provided, this variable will hold
! the parsed line.
!\\
!\\
! This routine splits the input line (or the next line of an open file with
! ID IU\_HCO), using the HEMCO separator (default: space) as separator. The
! resulting elements are then passed to the specified output characters and
! integers. For example, to pass the 5th element of a line to variable int1,
! set int1cl to 5, etc. An error will be returned (STAT=100) if any of the
! output columns exceeds the number of line elements. The optional argument
! optcl can be used to denote an optional value, e.g. no error is returned
! if the value at position optcl cannot be read. Only one optional value can
! be specified.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadAndSplit_Line( HcoConfig, IU_HCO,  char1, chr1cl, &
                                char2,     chr2cl,  char3, chr3cl, &
                                char4,     chr4cl,  char5, chr5cl, &
                                char6,     chr6cl,  char7, chr7cl, &
                                char8,     chr8cl,  char9, chr9cl, &
                                char10,    chr10cl,                &
                                int1,      int1cl,  int2,  int2cl, &
                                int3,      int3cl,  STAT,  inLine, &
                                outLine,   optcl                    )
!
! !USES:
!
    USE CHARPAK_Mod,        ONLY : STRREPL, STRSPLIT
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),    POINTER                 :: HcoConfig
    INTEGER,            INTENT(IN   )           :: IU_HCO
    INTEGER,            INTENT(IN   )           :: chr1cl
    INTEGER,            INTENT(IN   )           :: chr2cl
    INTEGER,            INTENT(IN   )           :: chr3cl
    INTEGER,            INTENT(IN   )           :: chr4cl
    INTEGER,            INTENT(IN   )           :: chr5cl
    INTEGER,            INTENT(IN   )           :: chr6cl
    INTEGER,            INTENT(IN   )           :: chr7cl
    INTEGER,            INTENT(IN   )           :: chr8cl
    INTEGER,            INTENT(IN   )           :: chr9cl
    INTEGER,            INTENT(IN   )           :: chr10cl
    INTEGER,            INTENT(IN   )           :: int1cl
    INTEGER,            INTENT(IN   )           :: int2cl
    INTEGER,            INTENT(IN   )           :: int3cl
    CHARACTER(LEN=255), INTENT(IN   ), OPTIONAL :: inLINE
    INTEGER,            INTENT(IN   ), OPTIONAL :: optcl
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*),   INTENT(INOUT)           :: char1
    CHARACTER(LEN=*),   INTENT(INOUT)           :: char2
    CHARACTER(LEN=*),   INTENT(INOUT)           :: char3
    CHARACTER(LEN=*),   INTENT(INOUT)           :: char4
    CHARACTER(LEN=*),   INTENT(INOUT)           :: char5
    CHARACTER(LEN=*),   INTENT(INOUT)           :: char6
    CHARACTER(LEN=*),   INTENT(INOUT)           :: char7
    CHARACTER(LEN=*),   INTENT(INOUT)           :: char8
    CHARACTER(LEN=*),   INTENT(INOUT)           :: char9
    CHARACTER(LEN=*),   INTENT(INOUT)           :: char10
    INTEGER,            INTENT(INOUT)           :: int1
    INTEGER,            INTENT(INOUT)           :: int2
    INTEGER,            INTENT(INOUT)           :: int3
    CHARACTER(LEN=255), INTENT(  OUT), OPTIONAL :: outLINE
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(INOUT)           :: STAT
!
! !REVISION HISTORY:
!  28 Aug 2013 - C. Keller - Initial version
!  11 Dec 2013 - C. Keller - Added optional arguments inLine and outLine
!  29 Dec 2014 - C. Keller - Added optional argument optcl. Now use wrapper
!                            routines READCHAR and READINT.
!  13 Mar 2015 - C. Keller - Added check for include files.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: N, OPT, STRLEN, RC
    CHARACTER(LEN=255)    :: LINE
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
       ! Read line
       CALL HCO_READLINE( IU_HCO, LINE, EOF, RC )

       ! Return w/ error
       IF ( RC /= HCO_SUCCESS ) THEN
          STAT = 999
          RETURN
       ENDIF

       ! End of file
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
    IF ( LINE(1:1) == HCO_CMT ) THEN
       STAT = 1
       RETURN
    ENDIF

    ! Get string length
    STRLEN = LEN(TRIM(LINE))

    ! Return here with flag = 5 is line is opening a (shortcut) bracket.
    IF ( STRLEN > 3 ) THEN
       IF ( LINE(1:3) == '(((' ) THEN
          STAT = 5
          RETURN
       ENDIF

       ! Return here with flag = 6 is line is opening a (shortcut) bracket.
       IF ( LINE(1:3) == ')))' ) THEN
          STAT = 6
          RETURN
       ENDIF
    ENDIF

    ! Return with flag = 1000 if this is a link to an include file.
    IF ( STRLEN > 11 ) THEN
       IF ( LINE(1:10) == '>>>include' ) THEN
          IF ( PRESENT(outLINE) ) outLINE = outLINE(12:STRLEN)
          STAT = 1000
          RETURN
       ENDIF
    ENDIF

    ! Split line into columns
    CALL STRREPL ( LINE, HCO_TAB, HCO_SPC )
    CALL STRSPLIT( LINE, HCO_SPC, SUBSTR, N )

    ! Also ignore empty lines
    IF ( N <= 1 ) THEN
       STAT = 1
       RETURN
    ENDIF

    ! Are there any optional lines?
    IF ( PRESENT(optcl) ) THEN
       OPT = optcl
    ELSE
       OPT = -1
    ENDIF

    ! ---------------------------------------------------------------------
    ! Read characters as specified and write them into given variables
    ! ---------------------------------------------------------------------

    CALL READCHAR( LINE, SUBSTR, N, chr1cl,  char1,  OPT, STAT )
    IF ( STAT == 100 ) RETURN
    CALL READCHAR( LINE, SUBSTR, N, chr2cl,  char2,  OPT, STAT )
    IF ( STAT == 100 ) RETURN
    CALL READCHAR( LINE, SUBSTR, N, chr3cl,  char3,  OPT, STAT )
    IF ( STAT == 100 ) RETURN
    CALL READCHAR( LINE, SUBSTR, N, chr4cl,  char4,  OPT, STAT )
    IF ( STAT == 100 ) RETURN
    CALL READCHAR( LINE, SUBSTR, N, chr5cl,  char5,  OPT, STAT )
    IF ( STAT == 100 ) RETURN
    CALL READCHAR( LINE, SUBSTR, N, chr6cl,  char6,  OPT, STAT )
    IF ( STAT == 100 ) RETURN
    CALL READCHAR( LINE, SUBSTR, N, chr7cl,  char7,  OPT, STAT )
    IF ( STAT == 100 ) RETURN
    CALL READCHAR( LINE, SUBSTR, N, chr8cl,  char8,  OPT, STAT )
    IF ( STAT == 100 ) RETURN
    CALL READCHAR( LINE, SUBSTR, N, chr9cl,  char9,  OPT, STAT )
    IF ( STAT == 100 ) RETURN
    CALL READCHAR( LINE, SUBSTR, N, chr10cl, char10, OPT, STAT )
    IF ( STAT == 100 ) RETURN

    ! ---------------------------------------------------------------------
    ! Read integers as specified and write them into given variables.
    ! Value -999 is returned for wildcard characters.
    ! ---------------------------------------------------------------------

    CALL READINT( HcoConfig%ExtList, LINE, SUBSTR, N, int1cl, int1, OPT, STAT )
    IF ( STAT == 100 ) RETURN
    CALL READINT( HcoConfig%ExtList, LINE, SUBSTR, N, int2cl, int2, OPT, STAT )
    IF ( STAT == 100 ) RETURN
    CALL READINT( HcoConfig%ExtList, LINE, SUBSTR, N, int3cl, int3, OPT, STAT )
    IF ( STAT == 100 ) RETURN

  END SUBROUTINE ReadAndSplit_Line
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: READCHAR
!
! !DESCRIPTION: Subroutine READCHAR is a helper routine to read character
! values from the HEMCO configuration file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READCHAR ( LINE, SUBSTR, N, chrcl, charout, OPT, STAT )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=255), INTENT(IN   )    :: LINE
    CHARACTER(LEN=255), INTENT(IN   )    :: SUBSTR(255)
    INTEGER,            INTENT(IN   )    :: N
    INTEGER,            INTENT(IN   )    :: chrcl
    INTEGER,            INTENT(IN   )    :: OPT
!
! !INPUT/OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*),   INTENT(INOUT)    :: charout
    INTEGER,            INTENT(INOUT)    :: STAT
!
! !REVISION HISTORY:
!  29 Dec 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    IF ( chrcl > 0 ) THEN
       IF ( chrcl > N ) THEN
          IF ( chrcl /= OPT ) THEN
             WRITE(*,*) 'Not enough elements in: '//TRIM(LINE)
             STAT = 100
             RETURN
          ELSE
             charout = ''
          ENDIF
       ELSE
          READ( SUBSTR(chrcl), '(a)' ) charout
       ENDIF
    ENDIF
    charout = ADJUSTL(charout)

  END SUBROUTINE READCHAR
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: READINT
!
! !DESCRIPTION: Subroutine READINT is a helper routine to read integer
! values from the HEMCO configuration file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE READINT ( ExtList, LINE, SUBSTR, N, intcl, intout, OPT, STAT )
!
! !USES:
!
    USE HCO_EXTLIST_MOD, ONLY : HCO_GetOpt
!
! !INPUT PARAMETERS:
!
    TYPE(Ext),          POINTER          :: ExtList
    CHARACTER(LEN=255), INTENT(IN   )    :: LINE
    CHARACTER(LEN=255), INTENT(IN   )    :: SUBSTR(255)
    INTEGER,            INTENT(IN   )    :: N
    INTEGER,            INTENT(IN   )    :: intcl
    INTEGER,            INTENT(IN   )    :: OPT
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(INOUT)    :: intout
    INTEGER,            INTENT(INOUT)    :: STAT
!
! !REVISION HISTORY:
!  29 Dec 2014 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    IF ( intcl > 0 ) THEN
       IF ( intcl > N ) THEN
          IF ( intcl /= OPT ) THEN
             WRITE(*,*) 'Not enough elements in: '//TRIM(LINE)
             STAT = 100
             RETURN
          ELSE
             intout = -999
          ENDIF
       ELSE
          ! Check for wildcard
          IF ( SUBSTR(intcl) == TRIM(HCO_GetOpt(ExtList,'Wildcard')) ) THEN
             intout = -999
          ELSE
             READ( SUBSTR(intcl), * ) intout
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE READINT
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_cID
!
! !DESCRIPTION: Subroutine Get\_cID searches the whole ConfigList for an entry
! with the given ScalID and returns the corresponding container ID cID.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_cID( ScalID, HcoConfig, cID, RC )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN   )   :: scalID
    TYPE(ConfigObj), POINTER :: HcoConfig
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(  OUT)   :: cID
!
! !INPUT/OUTPUTP PARAMETERS:
!
    INTEGER, INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  18 Sep 2013 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ListCont), POINTER  :: Lct

    ! Scalars
    CHARACTER(LEN=255)       :: MSG, LOC
    CHARACTER(LEN= 31)       :: strID

    ! Enter
    LOC = 'Get_cID (hco_config_mod.F90)'
    cID = -999

    ! Loop over all containers
    Lct => HcoConfig%ConfigList
    DO WHILE ( ASSOCIATED ( Lct ) )

       ! Skip if data container not defined
       IF ( .NOT. ASSOCIATED(Lct%Dct) ) THEN
          Lct => Lct%NextCont
          CYCLE
       ENDIF

       ! Check if this container has desired scalID
       IF ( Lct%Dct%ScalID == ScalID ) THEN
          cID = Lct%Dct%cID
          EXIT
       ENDIF

       ! Move to archived next line
       Lct => Lct%NextCont
    ENDDO

    ! Free pointer
    Lct => NULL()

    ! cID must be positive!
    IF ( cID <= 0 ) THEN
       WRITE ( strID, * ) ScalID
       MSG = 'Cannot find ScalID' // TRIM(strID)
       PRINT *,'cID negative in HEMCO Get_cID'
       PRINT *, TRIM(MSG)
       RETURN
    ENDIF

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Get_cID
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ConfigList_AddCont
!
! !DESCRIPTION: Subroutine ConfigList\_AddCont adds a new (blank) container to
! the ConfigList list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConfigList_AddCont( Lct, List )
!
! !USES:
!
    USE HCO_DATACONT_Mod, ONLY : DataCont_Init
    USE HCO_DATACONT_Mod, ONLY : ListCont_Length
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont), POINTER       :: Lct
    TYPE(ListCont), POINTER       :: List
!
! !REVISION HISTORY:
!  17 Sep 2013 - C. Keller: Initialization (update)
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ListCont), POINTER :: NewLct
    INTEGER                 :: cID

    !======================================================================
    ! ConfigList_AddCont begins here
    !======================================================================

    ! Allocate container and create data structure.
    ! The DataCont_Init call creates a new data container (type DataCont)
    ! All HEMCO lists (ConfigList, ReadList, EmisList) point to this
    ! container!
    ALLOCATE ( NewLct )
    NewLct%Dct      => NULL()
    NewLct%NextCont => NULL()

    ! Get # of containers in list. Set new container ID (cID) to # of
    ! containers + 1.
    cID = ListCont_Length( List )
    cID = cID + 1
    CALL DataCont_Init ( NewLct%Dct, cID )

    ! Connect blank container with ConfigList list.
    NewLct%NextCont => List
    List            => NewLct

    ! Output pointer points to the new container
    Lct => NewLct

  END SUBROUTINE ConfigList_AddCont
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ScalID_Register
!
! !DESCRIPTION: Subroutine ScalID\_Register adds the scale factor IDs ScalIDs
! to the list of scale factor IDs.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ScalID_Register( Dct, HcoConfig, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(DataCont),  POINTER :: Dct
    TYPE(ConfigObj), POINTER :: HcoConfig
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!  29 Dec 2014 - C. Keller: Now add new container to end of list to allow
!                           list being updated while calling Register_Scal.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                    :: N, cID

    !======================================================================
    ! ScalID_Register begins here
    !======================================================================

    ! Check for every element of ScalIDs, if this scale factor ID is
    ! already a member of ScalIDList. If not, add it.
    DO N = 1, Dct%nScalID
       IF ( Dct%Scal_cID(N) < 0 ) CYCLE

       CALL ScalID2List( HcoConfig%ScalIDList, Dct%Scal_cID(N), RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *,'Error in ScaleID2List called from HEMCO ScalID_Register (1)'
          RETURN
       ENDIF
 
       ! Replace scale factor ID with container ID.
       CALL Get_cID ( Dct%Scal_cID(N), HcoConfig, cID, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *,'Error in Get_cID called from HEMCO ScalID_Register (1)'
          RETURN
       ENDIF
       Dct%Scal_cID(N) = cID

    ENDDO

    ! Also check for level scale factor IDs
    IF ( Dct%levScalID1 > 0 ) THEN
       CALL ScalID2List( HcoConfig%ScalIDList, Dct%levScalID1, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *,'Error in ScalID2List called from HEMCO ScalID_Register (2)'
          RETURN
       ENDIF
       CALL Get_cID ( Dct%levScalID1, HcoConfig, cID, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *,'Error in Get_cID called from HEMCO ScalID_Register (2)'
          RETURN
       ENDIF
       Dct%levScalID1 = cID
    ENDIF
    IF ( Dct%levScalID2 > 0 ) THEN
       CALL ScalID2List( HcoConfig%ScalIDList, Dct%levScalID2, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *,'Error in ScaleID2List called from HEMCO ScalID_Register (3)'
          RETURN
       ENDIF
       CALL Get_cID ( Dct%levScalID2, HcoConfig, cID, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          PRINT *,'Error in Get_cID called from HEMCO ScalID_Register (3)'
          RETURN
       ENDIF
       Dct%levScalID2 = cID
    ENDIF

    ! Vector Scal_cID of this container now points to cIDs
    Dct%Scal_cID_Set = .TRUE.

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ScalID_Register
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ScalID2List
!
! !DESCRIPTION: Subroutine ScalID2List adds the scale factor IDs ScalIDs
! to the list of scale factor IDs.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ScalID2List( ScalIDList, ID, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(ScalIDCont),   POINTER        :: ScalIDList
    INTEGER,            INTENT(IN   )  :: ID
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!  29 Dec 2014 - C. Keller: Now add new container to end of list to allow
!                           list being updated while calling Register_Scal.
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ScalIDCont), POINTER  :: NewScalIDCont
    TYPE(ScalIDCont), POINTER  :: TmpScalIDCont
    TYPE(ScalIDCont), POINTER  :: PrvScalIDCont

    ! Scalars
    LOGICAL                    :: IsInList

    !======================================================================
    ! ScalID2List begins here
    !======================================================================

    ! Initialize
    NewScalIDCont => NULL()
    TmpScalIDCont => NULL()
    PrvScalIDCont => NULL()

    ! Check for every element of ScalIDs, if this scale factor ID is
    ! already a member of ScalIDList. If not, add it.

    ! Check if already in list
    IsInList = .FALSE.
    TmpScalIDCont => ScalIDList
    PrvScalIDCont => TmpScalIDCont
    DO WHILE ( ASSOCIATED(TmpScalIDCont) )
       IF ( TmpScalIDCont%ScalID == ID ) THEN
          IsInList = .TRUE.
          EXIT
       ENDIF
       PrvScalIDCont => TmpScalIDCont
       TmpScalIDCont => TmpScalIDCont%NEXT
    ENDDO

    ! Add new container w/ this scal ID to (end of) list
    IF ( .NOT. IsInList ) THEN
       ALLOCATE ( NewScalIDCont )
       NewScalIDCont%ScalID =  ID
       NewScalIDCont%NEXT   => NULL()
       IF ( .NOT. ASSOCIATED(PrvScalIDCont) ) THEN
          ScalIDList => NewScalIDCont
       ELSE
          PrvScalIDCont%NEXT => NewScalIDCont
       ENDIF
!       NewScalIDCont%NEXT   => ScalIDList
!       ScalIDList           => NewScalIDCont
!       NewScalIDCont        => NULL()
    ENDIF

    ! Cleanup
    TmpScalIDCont => NULL()

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ScalID2List
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ScalID_Cleanup
!
! !DESCRIPTION: Subroutine ScalID\_Cleanup cleans up the internal ScalID
! list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ScalID_Cleanup( ScalIDList )
!
! !INPUT ARGUMENTS:
!
  TYPE(ScalIDCont),   POINTER  :: ScalIDList
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ScalIDCont), POINTER  :: TmpScalIDCont
    TYPE(ScalIDCont), POINTER  :: NxtScalIDCont

    !======================================================================
    ! ScalID_Cleanup begins here
    !======================================================================

    ! Walk through list and remove each element
    NxtScalIDCont => NULL()
    TmpScalIDCont => ScalIDList
    DO WHILE ( ASSOCIATED(TmpScalIDCont) )

       NxtScalIDCont      => TmpScalIDCont%NEXT
       TmpScalIDCont%NEXT => NULL()
       DEALLOCATE ( TmpScalIDCont )

       TmpScalIDCont => NxtScalIDCont
    ENDDO

    ! Exit
    TmpScalIDCont => NULL()
    NxtScalIDCont => NULL()
    ScalIDList    => NULL()

  END SUBROUTINE ScalID_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SpecName_Register
!
! !DESCRIPTION: Subroutine SpecName\_Register adds the species name SpecName
! to the list of species names.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SpecName_Register( HcoConfig, SpecName, RC )
!
! !USES:
!
    USE HCO_EXTLIST_MOD, ONLY : HCO_GetOpt
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),    POINTER       :: HcoConfig
    CHARACTER(LEN=*),   INTENT(IN   ) :: SpecName
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
    TYPE(SpecNameCont), POINTER  :: NewSpecNameCont
    TYPE(SpecNameCont), POINTER  :: TmpSpecNameCont
    LOGICAL                      :: IsInList

    !======================================================================
    ! SpecName_Register begins here
    !======================================================================

    ! Ignore if wildcard character. These fields will always be used!
    IF ( TRIM(SpecName) == TRIM(HCO_GetOpt(HcoConfig%ExtList,'Wildcard')) ) THEN
       RC = HCO_SUCCESS
       RETURN
    ENDIF

    ! Initialize
    NewSpecNameCont => NULL()
    TmpSpecNameCont => NULL()

    ! Check if already in list
    IsInList = .FALSE.
    TmpSpecNameCont => HcoConfig%SpecNameList
    DO WHILE ( ASSOCIATED(TmpSpecNameCont) )
       IF ( TRIM(TmpSpecNameCont%SpecName) == TRIM(SpecName) ) THEN
          IsInList = .TRUE.
          EXIT
       ENDIF
       TmpSpecNameCont => TmpSpecNameCont%NEXT
    ENDDO

    ! Add new container w/ this scal ID to (beginning) of list
    IF ( .NOT. IsInList ) THEN
       ALLOCATE ( NewSpecNameCont )
       NewSpecNameCont%SpecName =  SpecName
       NewSpecNameCont%NEXT     => HcoConfig%SpecNameList
       HcoConfig%SpecNameList   => NewSpecNameCont
       NewSpecNameCont          => NULL()
    ENDIF

    ! Cleanup
    TmpSpecNameCont => NULL()

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE SpecName_Register
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SpecName_Cleanup
!
! !DESCRIPTION: Subroutine SpecName\_Cleanup cleans up the internal SpecName
! list.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SpecName_Cleanup ( SpecNameList )
!
! !INPUT/OUTPUT ARGUMENT:
!
    TYPE(SpecNameCont), POINTER       :: SpecNameList
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
    ! Pointers
    TYPE(SpecNameCont), POINTER  :: TmpSpecNameCont
    TYPE(SpecNameCont), POINTER  :: NxtSpecNameCont

    !======================================================================
    ! SpecName_Cleanup begins here
    !======================================================================

    ! Initialize
    TmpSpecNameCont => NULL()
    NxtSpecNameCont => NULL()

    ! Walk through list and remove each element
    TmpSpecNameCont => SpecNameList
    DO WHILE ( ASSOCIATED(TmpSpecNameCont) )

       NxtSpecNameCont      => TmpSpecNameCont%NEXT
       TmpSpecNameCont%NEXT => NULL()
       DEALLOCATE ( TmpSpecNameCont )

       TmpSpecNameCont => NxtSpecNameCont
    ENDDO

    ! Exit
    TmpSpecNameCont => NULL()
    NxtSpecNameCont => NULL()
    SpecNameList    => NULL()

  END SUBROUTINE SpecName_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Config_GetnSpecies
!
! !DESCRIPTION: Function Config\_GetnSpecies is a wrapper function to
! get the number of (unique) species names in SpecNameList.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Config_GetnSpecies( HcoConfig ) RESULT( nSpecies )
!
! !INPUT ARGUMENT:
!
    TYPE(ConfigObj), POINTER       :: HcoConfig
!
! !RETURN VALUE:
!
    INTEGER :: nSpecies
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: THISRC

    !======================================================================
    ! Config_GetnSpecies begins here
    !======================================================================

    CALL Config_GetSpecAttr( HcoConfig, N=nSpecies, RC = THISRC )

  END FUNCTION Config_GetnSpecies
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Config_GetSpecNames
!
! !DESCRIPTION: Subroutine Config\_GetSpecNames is a wrapper routine to
! obtain the list of (unique) species names defined in SpecNameList.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_GetSpecNames( HcoConfig, SpecNames, nSpecies, RC )
!
! !INPUT ARGUMENT:
!
    TYPE(ConfigObj), POINTER          :: HcoConfig
!
! !OUTPUT PARAMTERS:
!
    CHARACTER(LEN=*), POINTER         :: SpecNames(:)
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)   :: nSpecies
    INTEGER,          INTENT(INOUT)   :: RC
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!EOP
!------------------------------------------------------------------------------
!BOC
    !======================================================================
    ! Config_GetSpecNames begins here
    !======================================================================

    CALL Config_GetSpecAttr( HcoConfig, N=nSpecies, SpecNames=SpecNames, RC=RC )

  END SUBROUTINE Config_GetSpecNames
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Config_getSpecAttr
!
! !DESCRIPTION: Subroutine Config\_GetSpecAttr returns the number of
! species names N and the vector of species names SpecNames.
! SpecNames must be of length nnSpecs, i.e. in order to obtain
! SpecNames, Config\_getSpecAttr has to be called twice:
! N = 0
! CALL Config\_getSpecAttr ( N=N, RC=RC )
! ALLOCATE(SpecNames(N))
! CALL Config\_getSpecAttr ( N=N, SpecNames=SpecNames, RC=RC )
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_GetSpecAttr( HcoConfig, N, SpecNames, RC )
!
! !INPUT ARGUMENT:
!
    TYPE(ConfigObj),    POINTER                    :: HcoConfig
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(INOUT)              :: N
    INTEGER,            INTENT(INOUT)              :: RC
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*),   POINTER,       OPTIONAL    :: SpecNames(:)
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(SpecNameCont), POINTER   :: TmpSpecNameCont
    INTEGER                       :: AS
    CHARACTER(LEN=255), PARAMETER :: &
         LOC = 'Config_GetSpecAttr (hco_config_mod.F90)'

    !======================================================================
    ! Config_GetSpecAttr begins here
    !======================================================================

    ! Initialize
    TmpSpecNameCont => NULL()

    ! Eventually allocate pointer
    IF ( PRESENT(SpecNames) ) THEN
       IF ( .NOT. ASSOCIATED(SpecNames) ) THEN
          IF ( N <= 0 ) THEN
             CALL HCO_ERROR ( HcoConfig%Err, &
                'Cannot allocate SpecNames - N is size 0 or smaller', RC, THISLOC=LOC )
             RETURN
          ENDIF
          ALLOCATE(SpecNames(N), STAT=AS )
          IF ( AS/= 0 ) THEN
             CALL HCO_ERROR ( HcoConfig%Err, &
                'SpecNames allocation error', RC, THISLOC=LOC )
             RETURN
          ENDIF
          SpecNames(:) = ''
       ELSEIF ( SIZE(SpecNames) /= N ) THEN
          CALL HCO_ERROR ( HcoConfig%Err, &
             'SpecNames size error', RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDIF

    ! Init
    N = 0

    ! Loop over entire list. Count number of containers and eventually
    ! write out the species names.
    TmpSpecNameCont => HcoConfig%SpecNameList
    DO WHILE ( ASSOCIATED(TmpSpecNameCont) )
       N = N + 1
       IF ( PRESENT(SpecNames) ) THEN
          SpecNames(N) = TRIM(TmpSpecNameCont%SpecName)
       ENDIF
       TmpSpecNameCont => TmpSpecNameCont%NEXT
    ENDDO

    ! Cleanup and return w/ success
    TmpSpecNameCont => NULL()

    RC = HCO_SUCCESS

  END SUBROUTINE Config_GetSpecAttr
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_ContNames
!
! !DESCRIPTION: Function Check\_Contnames compares the container names of
! two containers, ignoring the name 'tags', i.e. ignoring everything that
! follows double underscore (\_\_). For example, two containers with names
! "EDGAR\_NOX\_\_PNT" and "EDGAR\_NOX\_\_MOB" are considered equal, while
! "EDGAR\_NOX\_PNT" and "EDGAR\_NOX\_MOB" are not.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Check_ContNames( Lct1, Lct2 ) RESULT( SameName )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont), POINTER :: Lct1
    TYPE(ListCont), POINTER :: Lct2
!
! !RETURN VALUE:
!
    LOGICAL                 :: SameName
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=63)  :: name1, name2
    INTEGER            :: idx

    !======================================================================
    ! Check_ContNames begins here!
    !======================================================================

    SameName = .FALSE.
    name1 = 'a'
    name2 = 'b'

    idx = INDEX( TRIM(Lct1%Dct%cName), '__' )
    IF ( idx > 0 ) THEN
       name1 = Lct1%Dct%cName(1:idx)
    ELSE
       name1 = Lct1%Dct%cName
    ENDIF

    idx = INDEX( TRIM(Lct2%Dct%cName), '__' )
    IF ( idx > 0 ) THEN
       name2 = Lct2%Dct%cName(1:idx)
    ELSE
       name2 = Lct2%Dct%cName
    ENDIF

    IF ( TRIM(name1) == TRIM(name2) ) THEN
       SameName = .TRUE.
    ELSE
       SameName = .FALSE.
    ENDIF

  END FUNCTION Check_ContNames
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtractSrcDim
!
! !DESCRIPTION: Subroutine ExtractSrcDim extracts the source dimension
! attribute. Specifically, it checks if the field is expected to be 2D
! (xy) or 3D. Default 3D data is xyz, but it is also possible to explicitly
! define the number of vertical levels to be read, as well as the reading
! direction (up or down). For example, 'xy1' will be interpreted as reading
! only the first level, and 'xy27' will only read the first 27 levels. To
! reverse the vertical axis, use e.g. 'xy-1' to read only the top level,
! or 'xy-27' to read the top 27 levels, with the topmost level being put
! into the surface level.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtractSrcDim( HcoConfig, SrcDim, Dta, Lscal1, Lscal2, RC )
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER         :: HcoConfig
    CHARACTER(LEN=*), INTENT(IN   )   :: SrcDim
    TYPE(FileData),   POINTER         :: Dta
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)   :: Lscal1
    INTEGER,          INTENT(  OUT)   :: Lscal2
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  20 May 2015 - C. Keller   - Initial version
!  22 Jan 2016 - R. Yantosca - Bug fix, removed & in the middle of the line
!                              since the PGI compiler chokes on it.
!  26 Jan 2018 - C. Keller   - Add L1 & L2
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: i, idx, idx2
    INTEGER            :: strLen
    INTEGER            :: EmisUnit
    REAL(hp)           :: EmisL
    CHARACTER(LEN=255) :: str1, str2, tmpstr
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'ExtractSrcDim (hco_config_mod.F90)'

    !======================================================================
    ! ExtractSrcDim begins here
    !======================================================================

    MSG = 'Illegal source dimension ' // TRIM(srcDim) // &
          ' for file ' // TRIM(Dta%ncFile) // &
          '. Valid entries are e.g. xy or xyz.'

    ! Init output
    Lscal1 = -1
    Lscal2 = -1

    ! See if there is an arbitrary additional dimension. This must be added
    ! at the end of the string and be separated by a '+' sign
    idx = INDEX( TRIM(srcDim), '+' )
    IF ( idx > 0 ) THEN
       str1 = srcDim(1:(idx-1))
       str2 = srcDim((idx+1):LEN(srcDim))
    ELSE
       str1 = srcDim
       str2 = ''
    ENDIF

    ! 2D data:
    IF ( TRIM(str1) == 'xy' .OR. TRIM(str1) == '-' ) THEN
       Dta%SpaceDim = 2

    ! All other cases
    ELSE
       ! Character length
       strLen = LEN(TRIM(str1))

       ! There must be at least 3 characters (e.g. xyz)
       IF ( strLen < 3 ) THEN
          CALL HCO_ERROR ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! First two entries must be xy
       IF ( str1(1:2) /= 'xy' ) THEN
          CALL HCO_ERROR ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! If third entry is 'L', this means we have 2D data that shall be put
       ! into a particular level, e.g. xyL4 will cause the 2D data to be
       ! emitted into level 4.
       IF ( str1(3:3) == 'L' .OR. str1(3:3) == 'l' ) THEN
          IF ( strLen < 4 ) THEN
             CALL HCO_ERROR ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
          Dta%SpaceDim = 2
          ! Read levels to put emissions into:
          i=4
          IF ( str1(i:i) == '=' ) i = i + 1

          ! Reduce to data to be read
          tmpstr = str1(i:strLen)

          ! check if range of levels is provided, i.e. xyL=1:5
          idx = INDEX( TRIM(tmpstr), ':' )

          ! if multiple levels are provided (e.g. xyL=1:5)
          IF ( idx > 0 ) THEN

             ! Check for PBL flag. It is possible to emit stuff
             ! from the PBL up to e.g. level 30 (xyL=PBL:30)
             ! The call to ParseEmisL now returns three arguments: the emission
             ! level, the emission unit, and the emission scale factor. Ignore
             ! emission level and unit if scale factor is given.
             CALL ParseEmisL( tmpstr(1:(idx-1)), EmisL, EmisUnit, Lscal1 )
             Dta%EmisL1     = EmisL
             Dta%EmisL1Unit = EmisUnit
             CALL ParseEmisL( tmpstr((idx+1):LEN(tmpstr)), EmisL, EmisUnit, Lscal2 )
             Dta%EmisL2     = EmisL
             Dta%EmisL2Unit = EmisUnit

          ! if only one level is provided (e.g. xyL=5)
          ELSE
             CALL ParseEmisL( tmpstr, EmisL, EmisUnit, Lscal1 )
             Dta%EmisL1     = EmisL
             Dta%EmisL1Unit = EmisUnit
             Lscal2         = Lscal1
             Dta%EmisL2     = Dta%EmisL1
             Dta%EmisL2Unit = Dta%EmisL1Unit
          ENDIF
       ELSE

          ! If we get to here, it's 3D data
          Dta%SpaceDim = 3

          ! The third entry determines the vertical dimension.
          ! This can be 'z' (standard) or a number to explicitly define
          ! the vertical extension and direction.
          IF ( str1(3:3) /= 'z' ) THEN
             READ(str1(3:strLen),*) Dta%Levels
          ENDIF
       ENDIF
    ENDIF

    ! Eventually set additional dimension name and value
    IF ( TRIM(str2) /= '' ) THEN
       MSG = 'Cannot extract arbitrary dimension from ' &
           // TRIM(srcDim) // ' for file ' // TRIM(Dta%ncFile) &
           // ' - arbitrary dimensions must follow a `+` sign ' &
           // 'and contain the name/value pair, e.g. xyz+"ens"=3'
       idx = INDEX( TRIM(str2), '=' )
       IF ( idx <= 0 ) THEN
          CALL HCO_ERROR( HcoConfig%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Extract dimension name. Eventually remove '"' character at
       ! beginning
       IF ( str2(1:1) == '"' .OR. &
            str2(1:1) == '`'       ) THEN
          Dta%ArbDimName = str2(2:(idx-1))
       ELSE
          Dta%ArbDimName = str2(1:(idx-1))
       ENDIF

       ! Extract dimension value. Eventually remove trailing '"'
       ! character. The string value itself will be evaluated when
       ! reading the file (in hcoio_dataread_mod.F90).
       strlen = LEN(TRIM(str2))
       IF ( str2(strlen:strlen) == '"' .OR. &
            str2(strlen:strlen) == '`'       ) THEN
          Dta%ArbDimVal = str2((idx+1):(strlen-1))
       ELSE
          Dta%ArbDimVal = str2((idx+1):(strlen))
       ENDIF

       ! Verbose
       IF ( HcoConfig%amIRoot .AND. HCO_IsVerb(HcoConfig%Err,2) ) THEN
          WRITE(MSG,*) 'Will use additional dimension on file ', &
             TRIM(Dta%ncFile), ': ', TRIM(Dta%ArbDimName), ' = ', &
             TRIM(Dta%ArbDimVal)
          CALL HCO_MSG(HcoConfig%Err,MSG)
       ENDIF
    ENDIF

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ExtractSrcDim
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ConfigInit
!
! !DESCRIPTION: Subroutine ConfigInit is a wrapper routine to initialize the
!  HEMCO configuration object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ConfigInit ( HcoConfig, RC, nModelSpecies )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN), OPTIONAL  :: nModelSpecies  ! # model species
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ConfigObj), POINTER       :: HcoConfig
    INTEGER,         INTENT(INOUT) :: RC             ! Success/fail
!
! !REVISION HISTORY:
!  16 Feb 2016 - C. Keller: Initialization (update)
!  23 Oct 2018 - M. Sulprizio- Add nModelSpecies to represent all species from
!                              external model (i.e. advected+chemical species)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, AS

    !=====================================================================
    ! ConfigInit begins here!
    !=====================================================================

    ALLOCATE(HcoConfig)
    HcoConfig%ConfigFileName = ''
    HcoConfig%ROOT           = ''
    HcoConfig%ConfigFileRead = .FALSE.
    HcoConfig%ConfigList     => NULL()
    HcoConfig%ScalIDList     => NULL()
    HcoConfig%SpecNameList   => NULL()
    HcoConfig%ExtList        => NULL()
    HcoConfig%Err            => NULL()

    IF ( PRESENT( nModelSpecies ) ) THEN

       ! Initialize vector w/ species information
       HcoConfig%nModelSpc = nModelSpecies
       IF ( nModelSpecies > 0 ) THEN
          ALLOCATE ( HcoConfig%ModelSpc( nModelSpecies ), STAT=AS )
          IF ( AS /= 0 ) THEN
             CALL HCO_ERROR( HcoConfig%Err, 'ModelSpecies', RC )
             RETURN
          ENDIF

          ! Initalize species information. The effective values for species
          ! names, model IDs, etc. are set in the HEMCO-model interface
          ! routine.
          DO I = 1, nModelSpecies
             HcoConfig%ModelSpc(I)%HcoID      = I
             HcoConfig%ModelSpc(I)%ModID      = -1
             HcoConfig%ModelSpc(I)%SpcName    = ''
          ENDDO
       ENDIF

    ENDIF

  END SUBROUTINE ConfigInit
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ParseEmisL
!
! !DESCRIPTION: parses the emission level.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ParseEmisL ( str, EmisL, EmisUnit, ScalID )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),  INTENT(IN ) :: str
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(hp),          INTENT(OUT) :: EmisL
    INTEGER,           INTENT(OUT) :: EmisUnit
    INTEGER,           INTENT(OUT) :: ScalID
!
! !REVISION HISTORY:
!  09 May 2016 - C. Keller: Intial version.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: nchar, idx

    !======================================================================
    ! ParseEmisL begins here!
    !======================================================================

    ! Init
    EmisUnit = HCO_EMISL_LEV
    ScalID   = -1

    IF ( TRIM(str) == 'PBL' ) THEN
      EmisL    = 0.0_hp
      EmisUnit = HCO_EMISL_PBL
    ELSE
       ! extract scale factor if string starts with 'SCAL' or 'scal'
       nchar = LEN(str)
       IF ( nchar > 4 ) THEN
          IF ( str(1:4)=='SCAL' .OR. str(1:4)=='scal' ) THEN
             READ(str(5:nchar),*) ScalID
             EmisUnit = -1
             EmisL    = -1.0
          ENDIF
       ENDIF

       ! check for elevation unit flag (e.g. 1000m)
       IF ( ScalID < 0 ) THEN
          idx = INDEX(TRIM(str),'m')
          IF ( idx > 0 ) THEN
             READ(str(1:(idx-1)),*) EmisL
             EmisUnit = HCO_EMISL_M
          ELSE
             READ(str,*) EmisL
          ENDIF
       ENDIF
    ENDIF

  END SUBROUTINE ParseEmisL
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CheckForDuplicateName
!
! !DESCRIPTION: Subroutine CheckForDuplicateName checks if there is a
! container in the container linked list that has the same name as the
! name given as input argument.
!\\
!\\
! !INTERFACE:
!
  Subroutine CheckForDuplicateName( HcoConfig, cName, RC )
!
! !INPUT ARGUMENT:
!
    TYPE(ConfigObj) , POINTER    :: HcoConfig  ! HEMCO config obj
    CHARACTER(LEN=*), INTENT(IN) :: cName
!
! !OUTPUT ARGUMENT:
!
    INTEGER, INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  20 Jul 2018 - C. Keller: Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ListCont), POINTER :: ThisLct => NULL()
    LOGICAL                 :: Duplicate
    CHARACTER(LEN=255)      :: tmpName, MSG

    !======================================================================
    ! CheckForDuplicateName begins here!
    !======================================================================

    ! Init
    RC = HCO_SUCCESS
    Duplicate = .FALSE.

    ! Pass name to clear spaces
    tmpName = ADJUSTL(cName)

    ! Walk through list and check for duplicate. Exit if found
    ThisLct => HcoConfig%ConfigList
    DO WHILE ( ASSOCIATED ( ThisLct ) )

       ! Skip if data container not defined
       IF ( .NOT. ASSOCIATED(ThisLct%Dct) ) THEN
          ThisLct => ThisLct%NextCont
          CYCLE
       ENDIF

       ! Check if this container has desired scalID
       IF ( TRIM(ThisLct%Dct%cName) == TRIM(tmpName) ) THEN
          Duplicate = .TRUE.
          EXIT
       ENDIF

       ! Move to next container
       ThisLct => ThisLct%NextCont
    ENDDO

    IF ( Duplicate ) THEN
       MSG = 'Error: HEMCO field already exists:'//TRIM(cName)
       CALL HCO_ERROR ( HcoConfig%Err, MSG, RC )
       RETURN
    ENDIF

  END SUBROUTINE CheckForDuplicateName
  !EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Hco_GetTagInfo
!
! !DESCRIPTION: Subroutine HCO\_GETTAGINFO retrieves basic information about
! tags given a wildcard string.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Hco_GetTagInfo( tagID, HcoConfig, Found, &
                             RC,    N,         tagName,  nTags            )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*),   INTENT(IN)  :: tagID       ! ID of tag (e.g. wildcard)
    TYPE(ConfigObj),    POINTER     :: HcoConfig   ! HEMCO Config object
    INTEGER,            OPTIONAL    :: N           ! index (1 to # tags)
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,            INTENT(OUT) :: Found       ! Item found?
    INTEGER,            INTENT(OUT) :: RC          ! Return code
    CHARACTER(LEN=255), OPTIONAL    :: tagName     ! tag name for index N
    INTEGER,            OPTIONAL    :: nTags       ! # tags
!
! !REMARKS:
!
! !REVISION HISTORY:
!  23 Oct 2018 - M. Sulprizio- Initial version based on routine Get_TagInfo in
!                              GEOS-Chem's Headers/state_diag_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: D,         numTags
    LOGICAL            :: isNumTags, isTagName, isN

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg,    ThisLoc,   Nstr

    !=======================================================================
    ! Hco_GetTagInfo begins here
    !=======================================================================

    ! Initialize
    ErrMsg     = ''
    Found      = .TRUE.
    numTags    = 0

    ! Optional arguments present?
    isN        = PRESENT( N       )
    isTagName  = PRESENT( TagName )
    isNumTags  = PRESENT( nTags   )

    ! Exit with error if getting tag name but index not specified
    IF ( isTagName .AND. .NOT. isN ) THEN
       ErrMsg = 'Index must be specified if retrieving an individual tag name'
       CALL HCO_ERROR( HcoConfig%Err, ErrMsg, RC )
       RETURN
    ENDIF

    !=======================================================================
    ! Get number of tags
    !=======================================================================
    SELECT CASE( TRIM( tagId ) )
       CASE( 'ALL'     )
          numTags = HcoConfig%nModelSpc
       CASE( 'ADV'     )
          numTags = HcoConfig%nModelAdv
       CASE DEFAULT
          FOUND = .FALSE.
          ErrMsg = 'Handling of tagId ' // TRIM(tagId) // &
                   ' is not implemented for getting number of tags'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC )
          RETURN
    END SELECT

    !=======================================================================
    ! Sanity checks -- exit under certain conditions
    !=======================================================================

    ! If not getting tag name then set nTags and exit
    IF ( .NOT. isTagName ) THEN
       nTags = numTags
       RETURN
    ENDIF

    ! Exit with error if index exceeds number of tags for this wildcard
    IF ( isTagName .AND. .NOT. isN ) THEN
       ErrMsg = 'Index must be greater than total number of tags for wildcard' &
                // TRIM(tagId)
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC )
       RETURN
    ENDIF

    !=======================================================================
    ! Get mapping index
    !=======================================================================
    SELECT CASE( TRIM( tagID ) )
       CASE( 'ALL', 'ADV' )
          D = N
       CASE DEFAULT
          FOUND = .FALSE.
          ErrMsg = 'Handling of tagId ' // TRIM( tagId ) // &
                   ' is not implemented for getting tag name'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC )
          RETURN
    END SELECT

    !=======================================================================
    ! Return the tag name
    !=======================================================================
    tagName = HcoConfig%ModelSpc(D)%SpcName

  END SUBROUTINE Hco_GetTagInfo
!EOC
END MODULE HCO_Config_Mod
