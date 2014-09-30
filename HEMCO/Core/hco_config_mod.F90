!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
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
! except for the input data associated with an disabled extension. 
! SetReadList does many more logical checks and adds all data used
! by HEMCO to ReadList. Scale factors not used by any of the base 
! emissions and base emission fields that won't be used are skipped
! in this step.
!\\
!\\ 
! All data fields are saved in individual data containers, which are 
! organized in the ConfigList. Hence, ConfigList is a collection of 
! data containers, with every container representing an entry of the
! configuration file. Each data container has its unique container ID
! for identification. Note that all HEMCO lists (ConfigList, ReadList,
! EmisList) access the same containers.
!\\
!\\
! The configuration file provides all source file information of the
! emission fields and scale factors to be used. It is read at the
! beginning of a simulation run.
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
  USE HCO_FILEDATA_MOD,       ONLY : FileData
  USE HCO_DATACONT_MOD,       ONLY : DataCont, ListCont, SclMax 
  USE HCO_STATE_MOD,          ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: SetReadList 
  PUBLIC  :: Config_ReadFile
  PUBLIC  :: GetNextCont
  PUBLIC  :: Config_Cleanup
  PUBLIC  :: Config_ScalIDinUse
  PUBLIC  :: Config_GetnSpecies
  PUBLIC  :: Config_GetSpecNames
!
! !PRIVATE:
!
  PRIVATE :: ReadSettings
  PRIVATE :: ExtSwitch2Buffer
  PRIVATE :: ConfigList_AddCont
  PRIVATE :: Config_ReadLine
  PRIVATE :: Config_ReadCont
  PRIVATE :: RegisterPrepare
  PRIVATE :: Get_targetID 
  PRIVATE :: Calc_Coverage
  PRIVATE :: Register_Base
  PRIVATE :: Register_Scal
  PRIVATE :: ReadAndSplit_Line
  PRIVATE :: Config_GetSpecAttr
!
! !REVISION HISTORY:
!  18 Jun 2013 - C. Keller: Initialization
!  08 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  08 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !MODULE TYPES/ARGUMENTS:
!
  ! Variables to store (unique) scale factor IDs and species names
  TYPE ScalIDCont
     INTEGER                   :: ScalID
     TYPE(ScalIDCont), POINTER :: NEXT
  END TYPE
  TYPE SpecNameCont
     CHARACTER(LEN=31)           :: SpecName
     TYPE(SpecNameCont), POINTER :: NEXT
  END TYPE

  ! Store unique scale factor IDs and species names in these lists
  TYPE(ScalIDCont),   POINTER  :: ScalIDList   => NULL()
  TYPE(SpecNameCont), POINTER  :: SpecNameList => NULL()

  ! Linked list w/ all input file information. For every 
  ! line of the input file, a separate data container will 
  ! be added to ConfigList
  TYPE(ListCont), POINTER     :: ConfigList => NULL()

  ! Configuration file read or not? 
  LOGICAL              :: ConfigFileRead = .FALSE. 

  ! SetReadList called or not?
  LOGICAL              :: SetReadListCalled = .FALSE.

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
! !IROUTINE: Config_Readfile
!
! !DESCRIPTION: Subroutine CONFIG\_READFILE reads the HEMCO configuration file
! and creates a data container for every (used) emission field in the config. 
! file. All containers become linked through the ConfigList linked list. 
! Note that lists EmisList and ReadList (created lateron)  will point to the 
! same containers, but ordering the containers in a manner that is most 
! efficient for the respective purpose. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_ReadFile( am_I_Root, ConfigFile, RC )
!
! !USES:
!
    USE inquireMod,       ONLY : findFreeLUN
    USE CharPak_Mod,      ONLY : STRREPL
    USE HCO_EXTLIST_MOD,  ONLY : AddExt
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN)    :: am_I_Root   ! root CPU?
    CHARACTER(LEN=*),   INTENT(IN)    :: ConfigFile  ! Full file name
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(INOUT) :: RC          ! Success?
!
! !REVISION HISTORY:
!  17 Sep 2012 - C. Keller: Initialization
!  03 Jan 2014 - C. Keller: Now use Config_ReadCont calls.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: IU_HCO, IOS
    LOGICAL              :: AIR, EOF
    CHARACTER(LEN=255)   :: LINE, MSG, LOC

    !======================================================================
    ! Config_ReadFile begins here
    !======================================================================

    ! Enter
    RC  = HCO_SUCCESS
    LOC = 'Config_ReadFile (hco_config_mod.F90)'

    ! Leave here if configuration file is already read 
    IF ( ConfigFileRead ) THEN
       RETURN
    ENDIF

    ! For convenience only
    AIR = am_I_Root

    ! Find free LUN
    IU_HCO = findFreeLUN()

    ! Open configuration file
    OPEN ( IU_HCO, FILE=TRIM( ConfigFile ), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       WRITE(*,*) 'Error reading ' // TRIM(ConfigFile)
       RC = HCO_FAIL
       RETURN 
    ENDIF 

    ! Register HEMCO core as extension Nr. 0 (default). The core 
    ! module is used by all HEMCO simulations, and the overall
    ! HEMCO settings are stored as options of this extension.
    CALL AddExt ( 'CORE', 0, 'all' )

    ! Loop until EOF 
    DO

       ! Read a line from the file, exit if EOF
       LINE = Config_ReadLine ( IU_HCO, EOF )
       IF ( EOF ) EXIT

       ! Replace tab characters in LINE (if any) w/ spaces
       CALL STRREPL( LINE, HCO_TAB(), HCO_SPC() )

       ! Read settings if this is beginning of settings section 
       IF ( INDEX ( LINE, 'BEGIN SECTION SETTINGS' ) > 0 ) THEN

          CALL ReadSettings( AIR, IU_HCO, EOF, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( EOF ) EXIT

       ! Read base emissions. This creates a new data container for each 
       ! base emission field. 
       ELSEIF ( INDEX ( LINE, 'BEGIN SECTION BASE EMISSIONS' ) > 0 ) THEN

          ! Read data and write into container
          CALL Config_ReadCont( AIR, IU_HCO, 1, EOF, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( EOF ) EXIT

       ! Read scale factors. This creates a new data container for each 
       ! scale factor.
       ELSE IF ( INDEX ( LINE, 'BEGIN SECTION SCALE FACTORS' ) > 0 ) THEN

          CALL Config_ReadCont( AIR, IU_HCO, 2, EOF, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( EOF ) EXIT

       ! Read masks. This creates a new data container for each mask 
       ELSE IF ( INDEX ( LINE, 'BEGIN SECTION MASKS' ) > 0 ) THEN

          CALL Config_ReadCont( AIR, IU_HCO, 3, EOF, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( EOF ) EXIT

       ! Read extension switches. This registers all enabled extensions.
       ELSE IF ( INDEX ( LINE, &
                         'BEGIN SECTION EXTENSION SWITCHES' ) > 0 ) THEN

          CALL ExtSwitch2Buffer( AIR, IU_HCO, EOF, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( EOF ) EXIT

       ! Read extension data. This reads the extension data of all 
       ! activated extensions (as base data).
       ELSE IF ( INDEX ( LINE, &
                         'BEGIN SECTION EXTENSION DATA' ) > 0 ) THEN

          CALL Config_ReadCont( AIR, IU_HCO, 1, EOF, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          IF ( EOF ) EXIT
       ENDIF
    ENDDO

    ! Close file
    CLOSE( UNIT=IU_HCO, IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       WRITE(*,*) 'Error closing ' // TRIM(ConfigFile)
       RC = HCO_FAIL
       RETURN 
    ENDIF 

    ! Configuration file is now read
    ConfigFileRead = .TRUE.

    ! Leave w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE Config_ReadFile
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE SetReadList( am_I_Root, HcoState, RC )
!
! !USES:
!
    USE HCO_DATACONT_Mod,    ONLY : cIDList_Create
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)    :: am_I_Root   ! root CPU?
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
    CALL HCO_ENTER ( 'SetReadList (hco_config_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ error if configuration file hasn't been read yet! 
    IF ( .NOT. ConfigFileRead ) THEN
       MSG = 'HEMCO configuration file not read!'
       CALL HCO_ERROR ( MSG, RC ) 
       RETURN
    ENDIF

    ! Prepare data in buffer. This call identifies all base fields
    ! that have to be read by this CPU. It also kicks out base 
    ! fields for emissions with an invalid species ID (if any) or
    ! if there are other base fields with higher priority.
    CALL RegisterPrepare ( am_I_Root, HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Register base emissions. In this step, we also redefine the 
    ! list UnqScalIDs to make sure that only those scale factors
    ! will be registered that are effectively used in the next step.
    CALL Register_Base( am_I_Root, HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Register scale factors based upon UnqScalIDs. 
    CALL Register_Scal( am_I_Root, HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Create cIDList which allows quick access to all data containers
    ! based on their container IDs cID
    CALL cIDList_Create ( am_I_Root, HcoState, ConfigList, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Don't need internal lists anymore.
    CALL ScalID_Cleanup
    CALL SpecName_Cleanup

    ! SetReadList has now been called
    SetReadListCalled = .TRUE.

    ! Leave w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE SetReadList 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE Config_ReadCont( am_I_Root, IU_HCO, DctType, EOF, RC )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,  ONLY : ExtNrInUse
    USE HCO_TIDX_Mod,     ONLY : HCO_ExtractTime
    USE HCO_FILEDATA_Mod, ONLY : FileData_Init
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN   ) :: am_I_Root  ! Root CPU?
    INTEGER, INTENT(IN   ) :: IU_HCO     ! Logfile LUN
    INTEGER, INTENT(IN   ) :: DctType    ! 1=base; 2=scale; 3=mask
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(INOUT) :: EOF        ! end of file encountered?
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(  OUT) :: RC         ! error code
! 
! !REVISION HISTORY: 
!  03 Jan 2014 - C. Keller - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                   :: STAT
    CHARACTER(LEN= 31)        :: cName
    CHARACTER(LEN=255)        :: srcFile
    CHARACTER(LEN= 31)        :: srcVar
    CHARACTER(LEN= 31)        :: srcTime
    CHARACTER(LEN=  1)        :: TmCycle 
    CHARACTER(LEN= 31)        :: srcDim
    CHARACTER(LEN= 31)        :: srcUnit
    CHARACTER(LEN= 31)        :: SpcName 
    CHARACTER(LEN=255)        :: Char1
    INTEGER                   :: Int1
    INTEGER                   :: Int2
    INTEGER                   :: Int3
    INTEGER                   :: Int4
    INTEGER                   :: N
    CHARACTER(LEN=255)        :: LOC, MSG 

    ! Arrays
    INTEGER                   :: SplitInts(SclMax)

    ! Pointers
    TYPE(ListCont), POINTER   :: Lct => NULL()
    TYPE(FileData), POINTER   :: Dta => NULL() 

    !=================================================================
    ! Config_ReadCont begins here!
    !=================================================================

    ! Enter
    LOC = 'Config_ReadCont (hco_config_mod.F90)'

    ! Repeat until end of base emissions section is found 
    DO

       !==============================================================
       ! Read line and get desired character strings
       ! Since base emissions, scale factors and masks have different 
       ! configuration file input parameter, need to use a different
       ! call for the three data types.
       !==============================================================
       IF ( DctType == 1 ) THEN
          CALL ReadAndSplit_Line ( am_I_Root, IU_HCO, cName,    2,  &
                                   srcFile,   3,      srcVar,   4,  &
                                   srcTime,   5,      TmCycle,  6,  &
                                   srcDim,    7,      srcUnit,  8,  &
                                   SpcName,   9,      Char1,   10,  &
                                   Int1,     11,      Int2,    12,  &
                                   Int3,      1,      STAT           )

       ELSE IF ( DctType == 2 ) THEN
          CALL ReadAndSplit_Line ( am_I_Root, IU_HCO, cName,    2,  &
                                   srcFile,   3,      srcVar,   4,  &
                                   srcTime,   5,      TmCycle,  6,  &
                                   srcDim,    7,      srcUnit,  8,  &
                                   SpcName,  -1,      Char1,   -1,  &
                                   Int1,      1,      Int2,     9,  &
                                   Int3,     -1,      STAT           )
  
       ELSE IF ( DctType == 3 ) THEN
          CALL ReadAndSplit_Line ( am_I_Root, IU_HCO, cName,    2,  &
                                   srcFile,   3,      srcVar,   4,  &
                                   srcTime,   5,      TmCycle,  6,  &
                                   srcDim,    7,      srcUnit,  8,  &
                                   SpcName,  -1,      Char1,   10,  &
                                   Int1,      1,      Int2,     9,  &
                                   Int3,     -1,      STAT           )
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
          CALL HCO_ERROR ( 'STAT == 100', RC, THISLOC=LOC )
          RETURN 
       ENDIF

       ! Output status should be 0 if none of the statuses above applies 
       IF ( STAT /= 0 ) THEN
          CALL HCO_ERROR ( 'STAT /= 0', RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! For base fields, check if this extension number is indeed in 
       ! use. Otherwise, we can ignore this line completely! The 
       ! extension switches are read and evaluated prior to the 
       ! extension data!
       IF ( DctType == 1 .AND. .NOT. ExtNrInUse( Int3 ) ) CYCLE 

       !==============================================================
       ! Create and fill list container and add to ConfigList 
       !==============================================================

       ! Add blank list container to ConfigList list. The container 
       ! is placed at the beginning of the list.
       CALL ConfigList_AddCont ( Lct, ConfigList )

       ! -------------------------------------------------------------
       ! Fill data container. 
       ! -------------------------------------------------------------

       ! Attributes used by all data types: data type number and 
       ! container name.
       Lct%Dct%DctType      = DctType
       Lct%Dct%cName        = cName

       ! Base container specific attributes
       IF ( DctType == 1 ) THEN       
  
          ! Set species name, extension number, emission category, 
          ! hierarchy
          Lct%Dct%SpcName       = SpcName 
          Lct%Dct%Cat           = Int1
          Lct%Dct%Hier          = Int2
          Lct%Dct%ExtNr         = Int3

          ! Set scale factor IDs into Scal_cID. These values will be
          ! replaced lateron with the container IDs (in register_base)!
          ! Note: SplitInts is of same lenghth SclMax as Scal_cID.
          ALLOCATE ( Lct%Dct%Scal_cID(SclMax) )
          CALL HCO_CharSplit( Char1, HCO_SEP(), HCO_WCD(), &
                              SplitInts, N, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          Lct%Dct%Scal_cID(1:SclMax) = SplitInts(1:SclMax)

          ! Register species name. A list of all species names can be
          ! returned to the atmospheric model to match HEMCO species 
          ! with model species (see Config\_GetSpecNames). 
          CALL SpecName_Register ( SpcName, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
    
       ! Scale factor & mask specific attributes
       ELSE IF ( DctType == 2 .OR. DctType == 3 ) THEN
            
          ! Set scale factor ID and data operator
          Lct%Dct%ScalID = Int1 
          Lct%Dct%Oper   = Int2

       ELSE
          CALL HCO_ERROR ( 'Invalid data type!', RC, THISLOC=LOC )
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
       ! used by multiple containers (at the moment, this is only
       ! important for an ESMF environment).
       ! -------------------------------------------------------------
       IF ( TRIM(srcFile) == '-' ) THEN
          IF ( .NOT. ASSOCIATED(Dta) ) THEN
             MSG = 'Cannot use previous data container: '//TRIM(cName)
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
          Lct%Dct%DtaHome = Lct%Dct%DtaHome - 1
       ELSE
          Dta => NULL()
          CALL FileData_Init ( Dta )

          ! Set source file name, file variable, and original data unit.
          ! In an ESMF environment, the source data will be imported 
          ! through ExtData by name, hence need to set ncFile equal to
          ! container name!
#if defined(ESMF_)
          Dta%ncFile    = cName
#else 
          Dta%ncFile    = srcFile
#endif
          Dta%ncPara    = srcVar
          Dta%OrigUnit  = srcUnit

          ! Extract information from time stamp character and pass values 
          ! to the corresponding container variables. If no time string is
          ! defined, keep default values (-1 for all of them)
          IF ( TRIM(srcTime) /= '-' ) THEN
             CALL HCO_ExtractTime( srcTime, Dta, RC ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF

          ! If the parameter ncPara is not defined, attempt to read data
          ! directly from configuration file instead of netCDF.
          IF ( TRIM(Dta%ncPara) == '-' ) THEN
             Dta%ncRead = .FALSE. 
          ENDIF

          ! Set time cycling behaviour. Possible values are: 
          ! - "C": cycling (CycleFlag = 1) --> Default
          ! - "R": range   (CycleFlag = 2)
          ! - "E": exact   (CycleFlag = 3)
          IF ( TRIM(TmCycle) == "R" ) THEN
             Dta%CycleFlag = 2
          ELSEIF ( TRIM(TmCycle) == "E" ) THEN
             Dta%CycleFlag = 3
          ELSEIF ( TRIM(TmCycle) == "C" ) THEN
             Dta%CycleFlag = 1
          ELSEIF ( TRIM(TmCycle) == "-" ) THEN
             Dta%CycleFlag = 1
          ELSE
             MSG = 'Invalid time cycling attribute: ' // &
                  TRIM(TmCycle) // ' --> ' // TRIM(cName)
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF

          ! Set space dimension. This will determine the dimension of the
          ! data array vector, i.e. 3D or 2D. Different time slices will
          ! be stored as different vector elements.
          IF ( TRIM(srcDim) == 'xy' ) THEN
             Dta%SpaceDim = 2
          ELSEIF ( TRIM(srcDim) == 'xyz' ) THEN
             Dta%SpaceDim = 3
          ENDIF

          ! For masks only: extract grid box edges. These will be used
          ! lateron to determine if emissions have to be considered by
          ! this CPU.
          IF ( DctType == 3 ) THEN
               
             ! Extract grid box edges. Need to be four values.
             CALL HCO_CharSplit ( Char1, HCO_SEP(), HCO_WCD(), & 
                                  SplitInts, N, RC ) 
             IF ( RC /= HCO_SUCCESS ) RETURN
             IF ( N /= 4 ) THEN
                MSG = 'Cannot properly read mask coverage: ' // &
                     TRIM(Lct%Dct%cName)
                CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF

             ! Save temporarily in year and month range. Will be
             ! reset lateron.
             Dta%ncYrs(1) = SplitInts(1)
             Dta%ncYrs(2) = SplitInts(2)
             Dta%ncMts(1) = SplitInts(3)
             Dta%ncMts(2) = SplitInts(4)
          ENDIF
       ENDIF

       ! Connect file data object of this data container.
       Lct%Dct%Dta => Dta

       ! Free list container for next cycle
       Lct => NULL()
    ENDDO

    ! Leave w/ success
    Dta => NULL()
    RC  =  HCO_SUCCESS 

  END SUBROUTINE Config_ReadCont
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE ExtSwitch2Buffer( am_I_Root, IU_HCO, EOF, RC )
!
! !USES:
!
    USE CHARPAK_Mod,        ONLY : STRREPL, STRSPLIT, TRANLC
    USE HCO_EXTLIST_MOD,    ONLY : AddExt, AddExtOpt
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)    :: am_I_Root   ! root CPU?
    INTEGER, INTENT(IN)    :: IU_HCO      ! HEMCO configfile LUN
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(INOUT) :: EOF         ! End of file?
    INTEGER, INTENT(INOUT) :: RC          ! Success/failure
!
! !REVISION HISTORY:
!  17 Sep 2013 - C. Keller: Initialization (update)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: I, N, Idx, ExtNr
    LOGICAL               :: Enabled
    CHARACTER(LEN=255)    :: LINE, LOC 
    CHARACTER(LEN=255)    :: SUBSTR(255), SPECS(255) 
    CHARACTER(LEN=1023)   :: OPTS

    !======================================================================
    ! ExtSwitch2Buffer begins here
    !======================================================================

    ! Enter
    LOC   = 'ExtSwitch2Buffer (hco_config_mod.F90)'
    RC    = HCO_SUCCESS
    ExtNr = 0

    ! Do until exit 
    DO 

       ! Read line 
       LINE = Config_ReadLine ( IU_HCO, EOF )

       ! Return if EOF
       IF ( EOF ) RETURN 

       ! Jump to next line if line is commented out
       IF ( LINE(1:1) == HCO_CMT() ) CYCLE

       ! Exit here if end of section encountered 
       IF ( INDEX ( LINE, 'END SECTION' ) > 0 ) RETURN 

       ! Check if these are options
       IF ( INDEX(LINE,'-->') > 0 ) THEN
          ! Only add if extension is defined!
          IF ( ExtNr > 0 .AND. Enabled ) THEN 
             CALL AddExtOpt( TRIM(LINE), ExtNr, RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDIF
          CYCLE
       ENDIF

       ! Split character string
       CALL STRREPL ( LINE, HCO_TAB(), HCO_TAB() )
       CALL STRSPLIT( LINE, HCO_SPC(), SUBSTR, N ) 

       ! Jump to next line if this line is empty
       IF ( N <= 1 ) CYCLE
        
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

       ! Register extension if enabled
       IF ( Enabled ) THEN

          ! Register extension name, number and species
          ! idx is the position of the species names
          idx = idx+1 
          READ( SUBSTR(1), * ) ExtNr
          CALL AddExt ( TRIM(SUBSTR(2)), ExtNr, SUBSTR(idx) )

          ! Register species (specNames)
          CALL STRSPLIT( SUBSTR(idx), HCO_SEP(), SPECS, N ) 
          IF ( N < 1 ) THEN
             CALL HCO_ERROR ( 'No species defined', RC, THISLOC=LOC )
             RETURN
          ENDIF
          DO I = 1, N
             CALL SpecName_Register ( SPECS(I), RC )
             IF ( RC /= HCO_SUCCESS ) RETURN
          ENDDO
       ENDIF
    ENDDO

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ExtSwitch2Buffer
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE ReadSettings( am_I_Root, IU_HCO, EOF, RC )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,    ONLY : AddExtOpt, GetExtOpt
    USE CHARPAK_MOD,        ONLY : STRREPL, STRSPLIT, TRANLC
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)    :: am_I_Root   ! root CPU?
    INTEGER, INTENT(IN)    :: IU_HCO      ! HEMCO configfile LUN
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(INOUT) :: EOF         ! End of file?
    INTEGER, INTENT(INOUT) :: RC          ! Success/failure
!
! !REVISION HISTORY:
!  17 Sep 2013 - C. Keller: Initialization (update)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: I, N, POS
    LOGICAL               :: verb, warn, track
    CHARACTER(LEN=255)    :: LINE, LOC
    CHARACTER(LEN=255)    :: LogFile
    CHARACTER(LEN=255)    :: DiagnPrefix

    !======================================================================
    ! ReadSettings begins here
    !======================================================================

    ! Enter
    LOC = 'ReadSettings (hco_config_mod.F90)'

    ! Defaults
    LogFile   = 'HEMCO.log'
    verb      = .FALSE.
    track     = .FALSE.
    warn      = .TRUE.

    !-----------------------------------------------------------------------
    ! Read settings and add them as options to core extensions
    !-----------------------------------------------------------------------

    ! Do until exit 
    DO 

       ! Read line 
       LINE = Config_ReadLine ( IU_HCO, EOF )

       ! Return if EOF
       IF ( EOF ) EXIT 

       ! Jump to next line if line is commented out
       IF ( LINE(1:1) == HCO_CMT() ) CYCLE

       ! Exit here if end of section encountered 
       IF ( INDEX ( LINE, 'END SECTION' ) > 0 ) EXIT 

       ! Ignore empty lines
       IF ( TRIM(LINE) == '' ) CYCLE

       ! Add this option to HEMCO core (extension 0)
       CALL AddExtOpt ( TRIM(LINE), 0, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

    ENDDO

    !-----------------------------------------------------------------------
    ! Extract values to initialize error module and set some further
    ! HEMCO variables. 
    !-----------------------------------------------------------------------

    ! Verbose mode?
    CALL GetExtOpt( 0, 'Verbose', OptValBool=verb, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Logfile to write into
    CALL GetExtOpt( 0, 'Logfile', OptValChar=Logfile, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Prompt warnings to logfile? 
    CALL GetExtOpt( 0, 'Show warnings', OptValBool=warn, RC=RC  )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Track code? 
    CALL GetExtOpt( 0, 'Track', OptValBool=track, RC=RC  )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If LogFile is equal to wildcard character, set LogFile to asterik 
    ! character. This will ensure that all output is written to standard
    ! output!
    IF ( TRIM(LogFile) == HCO_WCD() ) LogFile = '*'

    ! We should now have everything to define the HEMCO error settings
    CALL HCO_ERROR_SET ( LogFile, verb, warn, track, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ReadSettings
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Config_ReadLine 
!
! !DESCRIPTION: Subroutine ConfigRead\_Line reads a line from the
! emissions configuration file. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION Config_ReadLine( IU_HCO, EOF ) RESULT( LINE )
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN)  :: IU_HCO   ! HEMCO configfile LUN
!
! !OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(OUT) :: EOF      ! End of file?
! 
! !REVISION HISTORY: 
!  18 Sep 2013 - C. Keller - Initial version (adapted from B. Yantosca's code) 
!  15 Jul 2014 - R. Yantosca - Remove dependency on routine IOERROR
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: IOS
    CHARACTER(LEN=255) :: LINE, MSG

    !=================================================================
    ! Config_ReadLine begins here!
    !=================================================================

    ! Initialize
    EOF = .FALSE.

    ! Read a line from the file
    READ( IU_HCO, '(a)', IOSTAT=IOS ) LINE

    ! IO Status < 0: EOF condition
    IF ( IOS < 0 ) THEN
       EOF = .TRUE.
       RETURN
    ENDIF

    ! IO Status > 0: true I/O error condition
    IF ( IOS > 0 ) THEN
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, 100   ) IOS
100    FORMAT( 'ERROR ', i5, ' in Config_Readline (hco_config_mod.F90)' )
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       STOP
    ENDIF

  END FUNCTION Config_ReadLine
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE RegisterPrepare( am_I_Root, HcoState, RC ) 
!
! !USES:
!
    USE HCO_EXTLIST_MOD,  ONLY : ExtNrInUse
    USE HCO_STATE_Mod,    ONLY : HCO_GetHcoID
!
! !INPUT PARAMETERS: 

    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! Root CPU
    TYPE(HCO_State),  POINTER       :: HcoState   ! HEMCO state obj.
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC
! 
! !REVISION HISTORY: 
!  18 Sep 2013 - C. Keller - Initial version (update) 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ListCont), POINTER      :: Lct => NULL()
    INTEGER                      :: ThisCover, ThisHcoID, FLAG
    INTEGER                      :: lon1, lon2, lat1, lat2
    INTEGER                      :: cpux1, cpux2, cpuy1, cpuy2
    CHARACTER(LEN=255)           :: MSG
    LOGICAL                      :: verb

    !=================================================================
    ! RegisterPrepare begins here!
    !=================================================================

    ! Enter
    CALL HCO_ENTER ( 'RegisterPrepare', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN 

    ! Check for verbose flag
    verb = HCO_VERBOSE_CHECK() .and. am_I_Root

    ! Grid boundaries on this CPU. Will be needed to calculate 
    ! coverages. 
    ! NOTE: Use midpoints here because only those become defined in
    ! the ESMF environment (xedge and yedge are not used anywhere
    ! else in ESMF!). 
    cpux1 = FLOOR(MINVAL(HcoState%Grid%XMID))
    cpux2 = FLOOR(MAXVAL(HcoState%Grid%XMID))
    cpuy1 = CEILING(MINVAL(HcoState%Grid%YMID))
    cpuy2 = CEILING(MAXVAL(HcoState%Grid%YMID))

    ! verbose
    IF ( verb ) THEN
       WRITE(MSG,*) 'Start to prepare fields for registering!'
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) 'This CPU x-range: ', cpux1, cpux2
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) 'This CPU y-range: ', cpuy1, cpuy2
       CALL HCO_MSG(MSG)
    ENDIF

    ! Get next (first) line of ConfigList 
    CALL GetNextCont ( Lct, FLAG ) 

    ! Loop over all lines
    DO WHILE ( FLAG == HCO_SUCCESS ) 

       ! Check if data container defined
       IF ( .NOT. ASSOCIATED(Lct%Dct) ) THEN
          CALL GetNextCont( Lct, FLAG ); CYCLE
       ENDIF

       ! verbose
       IF ( verb ) THEN
          WRITE(MSG,*) 'Prepare ', TRIM(Lct%Dct%cName) 
          CALL HCO_MSG(MSG)
       ENDIF
 
       ! For base fields or data fields used in one of the HEMCO
       ! extensions:
       IF ( Lct%Dct%DctType == 1 ) THEN

          ! Only do for entries that will be used! 
          IF ( ExtNrInUse( Lct%Dct%ExtNr ) ) THEN

             ! Extract HEMCO species ID. This will return -1 for 
             ! undefined species and 0 for wildcard character.
             ThisHcoID = HCO_GetHcoID( Lct%Dct%SpcName, HcoState)
             Lct%Dct%HcoID = ThisHcoID

             ! verbose
             IF ( verb ) THEN
                WRITE(MSG,*) 'Assigned HEMCO species ID: ', Lct%Dct%HcoID
                CALL HCO_MSG(MSG)
             ENDIF

          ! Else: assign default value. These containers will be 
          ! removed in the next step!
          ELSE
             Lct%Dct%HcoID = -999
          ENDIF

       ! Calculate coverage for masks 
       ELSE IF ( Lct%Dct%DctType  == 3    .AND. &
                 Lct%Dct%Dta%Cover == -999        ) THEN

          ! Get mask edges
          lon1 = Lct%Dct%Dta%ncYrs(1)
          lat1 = Lct%Dct%Dta%ncYrs(2)
          lon2 = Lct%Dct%Dta%ncMts(1)
          lat2 = Lct%Dct%Dta%ncMts(2) 

          ThisCover = CALC_COVERAGE( lon1,  lon2,  &
                                     lat1,  lat2,  &
                                     cpux1, cpux2, &
                                     cpuy1, cpuy2   )
        
          ! Update container information
          Lct%Dct%Dta%Cover    = ThisCover 
          Lct%Dct%Dta%ncYrs(:) = -999
          Lct%Dct%Dta%ncMts(:) = -999

          IF ( verb ) THEN
             WRITE(MSG,*) 'Coverage: ', Lct%Dct%Dta%Cover
             CALL HCO_MSG(MSG)
          ENDIF
       ENDIF

       ! Advance to next line
       CALL GetNextCont ( Lct, FLAG ) 
    ENDDO

    ! Cleanup
    Lct => NULL()

    ! Return w/ success
    CALL HCO_LEAVE ( RC )

  END SUBROUTINE RegisterPrepare
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE Register_Base ( am_I_Root, HcoState, RC )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,       ONLY : ExtNrInUse
    USE HCO_READLIST_Mod,      ONLY : ReadList_Set
    USE HCO_DATACONT_Mod,      ONLY : DataCont_Cleanup
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(HCO_State),  POINTER       :: HcoState    ! HEMCO state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  18 Jun 2013 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ListCont), POINTER   :: Lct => NULL()

    ! Scalars
    INTEGER               :: N, cID, HcoID
    INTEGER               :: targetID, FLAG
    LOGICAL               :: Ignore, Add, Verb
    CHARACTER(LEN=255)    :: MSG

    !======================================================================
    ! Register_Base begins here
    !======================================================================

      ! Enter
    CALL HCO_ENTER ( 'Register_Base (hco_config_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    Verb = HCO_VERBOSE_CHECK() .and. am_I_Root

    ! Point to next (first) line in ConfigList 
    CALL GetNextCont ( Lct, FLAG )

    ! Loop over temporary arrays
    DO WHILE ( FLAG == HCO_SUCCESS )  

       ! Reset ignore flag
       Ignore = .FALSE.

       ! Skip entry if data container not defined
       IF ( .NOT. ASSOCIATED(Lct%Dct) ) THEN
          CALL GetNextCont ( Lct, FLAG ); CYCLE
       ENDIF

       ! Skip entry if it's not a base field 
       IF ( (Lct%Dct%DctType /= 1) ) THEN
          CALL GetNextCont ( Lct, FLAG ); CYCLE
       ENDIF

       ! If this base field is not used (either because it belongs to 
       ! an extension that is not enabled or because its HEMCO or 
       ! model species ID is undefined), we don't need this container 
       ! any more. Hence remove it.
       ! Note: Routine RegisterPrepare assigns negative HcoID's to all
       ! base fields with invalid ExtNr's, so it is ok to check only
       ! for HcoID here. If data is used in one of the HEMCO extensions
       ! and has a species flag of '*' (= always read), its species ID
       ! becomes set to 0 in registerPrepare. 
       HcoID = Lct%Dct%HcoID
       IF ( HcoID < 0 ) THEN
          Ignore = .TRUE.
       ELSE IF ( HcoID > 0 ) THEN
          IF ( HcoState%Spc(HcoID)%ModID < 0 ) Ignore = .TRUE.
       ENDIF

       IF ( Ignore ) THEN
          IF ( Verb ) THEN
             WRITE(MSG,*) &
                  'Register_Base: Ignore (and remove) base field ', &
                  TRIM(Lct%Dct%cName)
             CALL HCO_MSG(MSG)
          ENDIF

          ! Remove data container from list.
          CALL DataCont_Cleanup ( Lct%Dct )
          Lct%Dct => NULL()
          CALL GetNextCont ( Lct, FLAG ); CYCLE
       ENDIF

       ! Verbose mode 
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Register_Base: Checking ', TRIM(Lct%Dct%cName)
          CALL HCO_MSG(MSG)
       ENDIF

!       ! -------------------------------------------------------------
!       ! Eventually read data from file 
!       IF ( .NOT. Lct%Dct%Dta%ncRead ) THEN
!          CALL ReadFromConfig ( am_I_Root, HcoState, Lct, RC )
!          IF ( RC /= HCO_SUCCESS ) RETURN
!       ENDIF       

       ! -------------------------------------------------------------
       ! Extract vector of scale factor container IDs to be applied 
       ! to this base field (vector Scal_cID). For now, this container
       ! contains the scale factor IDs, hence need to convert to 
       ! container IDs. Beforehand, add scale factor IDs to internal 
       ! list of used scale factors (UnqScalIDs).
       CALL ScalID_Register ( Lct%Dct%Scal_cID, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       DO N = 1,SclMax
          IF ( Lct%Dct%Scal_cID(N) <= 0 ) CYCLE 
          CALL Get_cID ( Lct%Dct%Scal_cID(N), cID, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
          Lct%Dct%Scal_cID(N) = cID

          ! Vector Scal_cID of this container now points to cIDs
          Lct%Dct%Scal_cID_Set = .TRUE.
       ENDDO

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
       CALL Get_targetID( am_I_Root, HcoState, Lct, targetID, RC)
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! verbose 
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Container ID     : ', Lct%Dct%cID
          CALL HCO_MSG(MSG)
          WRITE(MSG,*) 'Assigned targetID: ', targetID
          CALL HCO_MSG(MSG)
       ENDIF

       ! Negative targetID is assigned to base data that doesn't need 
       ! to be considered either because it's a regional inventory
       ! with no spatial overlap with the region domain on this CPU
       ! or because there exist another inventory with higher 
       ! priority (and same category) that will overwrite these 
       ! emissions data anyway!
       IF ( targetID <= 0 ) THEN
          CALL DataCont_Cleanup ( Lct%Dct )
          Lct%Dct => NULL()
          CALL GetNextCont ( Lct, FLAG ); CYCLE
       ENDIF
        
       ! Pass targetID to container
       Lct%Dct%targetID = targetID

       ! Register container in ReadList. Containers will be listed 
       ! in the reading lists sorted by cID.  
       CALL ReadList_Set( am_I_Root, HcoState, Lct%Dct, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Print some information if verbose mode is on 
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Base field registered: ', TRIM(Lct%Dct%cName)
          CALL HCO_MSG(MSG)
       ENDIF

       ! Advance to next line
       CALL GetNextCont ( Lct, FLAG )
    ENDDO

    ! Cleanup 
    Lct => NULL()

    ! Return w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE Register_Base
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE Register_Scal( am_I_Root, HcoState, RC )
!
! !USES:
!
    USE HCO_ReadList_Mod, ONLY : ReadList_Set
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root   ! Are we on root CPU? 
    TYPE(HCO_State),  POINTER       :: HcoState    ! HEMCO state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  18 Jun 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ListCont),   POINTER :: Lct => NULL()
    TYPE(ScalIDCont), POINTER :: TmpScalIDCont => NULL() 

    ! Scalars
    INTEGER                   :: FLAG
    CHARACTER(LEN=255)        :: MSG
    CHARACTER(LEN=  5)        :: strID
    INTEGER                   :: ThisScalID
    LOGICAL                   :: Verb
  
    !======================================================================
    ! Register_Scal begins here
    !======================================================================

    ! Enter
    CALL HCO_ENTER ( 'Register_Scal (hco_config_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    Verb = HCO_VERBOSE_CHECK() .and. am_I_Root

    ! Loop over all scale factor ids
    TmpScalIDCont => ScalIDList
    DO WHILE ( ASSOCIATED( TmpScalIDCont ) )  

       ! Extract this scale factor ID
       ThisScalID = TmpScalIDCont%ScalID

       ! Make ThisLine point to first element of ConfigList
       Lct => NULL()
       CALL GetNextCont ( Lct, FLAG )         

       ! Loop over all lines in Input file and find the one with the
       ! correct scale factor ID 
       DO WHILE ( FLAG == HCO_SUCCESS )  

          ! Leave if this is the wanted container. 
          IF ( ASSOCIATED(Lct%Dct)) THEN
             IF ( Lct%Dct%ScalID == ThisScalID ) EXIT
          ENDIF

          ! Advance to next line otherwise
          CALL GetNextCont ( Lct, FLAG )         
       ENDDO

       ! Return error if scale factor ID not found
       IF ( .NOT. ASSOCIATED(Lct) ) THEN
          WRITE ( strID, * ) ThisScalID 
          MSG = 'Container ID not found: ' // strID
          CALL HCO_ERROR ( MSG, RC)
          RETURN 
       ENDIF
 
       ! Return w/ error if container not scale factor or mask
       IF ( Lct%Dct%DctType == 1 ) THEN
          WRITE ( strID, * ) ThisScalID 
          MSG = 'Container ID belongs to base field: ' // strID
          CALL HCO_ERROR ( MSG, RC)
          RETURN 
       ENDIF
 
       ! Eventually read scalar data (scale fields only)
       ! Spatially homogenous scale factors ('scalars') can directly 
       ! be defined in the configuration file (column 'Scalar'). All 
       ! file information (sourceFile, sourceVar, and time stamp) will 
       ! be ignored in this case. More than one scalar can be 
       ! specified, in which case these are interpreted as temporal 
       ! variations (7=day of week, 12=monthly, 24=hourly). 
!       IF ( .NOT. Lct%Dct%Dta%ncRead ) THEN 
!          CALL ReadFromConfig ( am_I_Root, HcoState, Lct, RC )
!          IF ( RC /= HCO_SUCCESS ) RETURN
!       ENDIF

       ! Register container in ReadList. Containers will be listed 
       ! in the reading lists sorted by cID.  
       CALL ReadList_Set( am_I_Root, HcoState, Lct%Dct, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Print some information if verbose mode is on 
       IF ( Verb ) THEN
          WRITE(MSG,*) 'Scale field registered: ', TRIM(Lct%Dct%cName)
          CALL HCO_MSG(MSG)
       ENDIF

       ! Advance
       TmpScalIDCont => TmpScalIDCont%NEXT 

    ENDDO

    ! Cleanup
    Lct           => NULL()
    TmpScalIDCont => NULL()

    ! Return w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE Register_Scal
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE Get_targetID( am_I_Root, HcoState, Lct, targetID, RC ) 
!
! !USES:
!
    USE HCO_DataCont_Mod, ONLY : ListCont_Find
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN   ) :: am_I_Root
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
!  11 Apr 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ListCont), POINTER   :: tmpLct => NULL()
    TYPE(ListCont), POINTER   :: mskLct => NULL()

    ! Scalars
    INTEGER                   :: HcoID, Cat, Hier, Scal, ExtNr, cID
    INTEGER                   :: tmpID
    INTEGER                   :: I, J, FLAG1, tmpCov
    LOGICAL                   :: found, verb, sameCont
    CHARACTER(LEN=255)        :: MSG
    CHARACTER(LEN=  7)        :: strID

    !======================================================================
    ! Get_targetID begins here
    !======================================================================

    ! Enter
    CALL HCO_ENTER ( 'Get_targetID (hco_config_mod.F90)', RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get verbose flag from HCO_State
    verb = HCO_VERBOSE_CHECK() .AND. am_I_Root

    ! Get Tracer ID, category and hierarchy of entry to be checked
    cID   = Lct%Dct%cID
    ExtNr = Lct%Dct%ExtNr 
    Cat   = Lct%Dct%Cat
    Hier  = Lct%Dct%Hier
    HcoID = Lct%Dct%HcoID

    ! By default, set target ID to container ID
    targetID = cID

    ! If species ID is zero, always read this field as is, i.e. don't
    ! skip it and don't add it to another field!
    ! Species ID become zero if the species ID entry in the
    ! configuration file is the wildcard character. 
    IF ( HcoID == 0 ) THEN
       CALL HCO_LEAVE ( RC )
       RETURN
    ENDIF

    ! Check all scale factors of the current container to see if one 
    ! of them is a mask that has no valid entries over the domain of 
    ! this CPU. In this case we don't have to consider this field at 
    ! all!
    DO I = 1, SclMax

       ! Check if it's a valid scale factor
       IF ( Lct%Dct%Scal_cID(I) < 0 ) CYCLE

       ! Find container with this container ID
       ! Note: this should always look up the container ID, but make
       ! check for safety's sake. 
       tmpID = Lct%Dct%Scal_cID(I)
       IF ( .NOT. Lct%Dct%Scal_cID_set ) THEN
          CALL ListCont_Find ( ConfigList, tmpID, 1, FOUND, mskLct )
       ELSE
          CALL ListCont_Find ( ConfigList, tmpID, 0, FOUND, mskLct )
       ENDIF

       ! Error if scale factor not found
       IF ( .NOT. FOUND ) THEN
          WRITE ( strID, * ) Lct%Dct%Scal_cID(I) 
          MSG = 'No scale factor with cID: ' // TRIM(strID) 
          CALL HCO_ERROR ( MSG, RC)
          RETURN
       ENDIF

       ! Check if this is a mask with zero coverage over this CPU, in
       ! which case we don't need to consider the base field at all! 
       IF ( (mskLct%Dct%DctType  == 3 ) .AND. &
            (mskLct%Dct%Dta%Cover == 0 )        ) THEN 
          targetID = -999 
          IF ( verb ) THEN 
             WRITE(MSG,*) 'Data not defined over this CPU, skip ' // &
                  TRIM(Lct%Dct%cName)
             CALL HCO_MSG(MSG)
          ENDIF

          ! Return
          CALL HCO_LEAVE ( RC )
          RETURN 
       ENDIF
    ENDDO ! I 
      
    ! Now find out if there is another base field for the same species, 
    ! emission category and extension number, but higher hierarchy.
    ! Such a field also needs to have full coverage over this CPU, 
    ! then we can ignore the current container. 
    
    ! Initialize looping pointer
    tmpLct => NULL()
    CALL GetNextCont ( tmpLct, FLAG1 ) 

    ! Loop over containers
    DO WHILE ( FLAG1 == HCO_SUCCESS )

       ! Advance to next container if data container not defined 
       IF ( .NOT. ASSOCIATED(tmpLct%Dct) ) THEN
          CALL GetNextCont ( tmpLct, FLAG1 ); CYCLE
       ENDIF
 
       ! Advance to next container if this is the current container 
       IF ( tmpLct%Dct%cID == cID ) THEN
          CALL GetNextCont ( tmpLct, FLAG1 ); CYCLE
       ENDIF

       ! Advance to next container if this is not a base field
       IF ( tmpLct%Dct%DctType /= 1 ) THEN
          CALL GetNextCont ( tmpLct, FLAG1 ); CYCLE
       ENDIF

       ! Advance to next container if not the same extension nr
       IF ( tmpLct%Dct%ExtNr /= ExtNr ) THEN
          CALL GetNextCont ( tmpLct, FLAG1 ); CYCLE
       ENDIF

       ! Advance to next container if not the same species 
       IF ( tmpLct%Dct%HcoID /= HcoID ) THEN
          CALL GetNextCont ( tmpLct, FLAG1 ); CYCLE
       ENDIF

       ! Advance to next container if not the same category
       IF ( tmpLct%Dct%Cat /= Cat ) THEN
          CALL GetNextCont ( tmpLct, FLAG1 ); CYCLE
       ENDIF

       ! Advance to next container if lower hierarchy
       IF ( tmpLct%Dct%Hier < Hier ) THEN
          CALL GetNextCont ( tmpLct, FLAG1 ); CYCLE
       ENDIF

       ! Check for coverage of tmpLct. Default = full coverage (1) 
       tmpCov = 1 

       ! Check all scale factors of tmpLct to see if this base 
       ! field has full coverage over this CPU domain or not. 
       DO I = 1, SclMax

          ! Check if it's a valid scale factor
          IF ( tmpLct%Dct%Scal_cID(I) < 0 ) CYCLE

          tmpID = tmpLct%Dct%Scal_cID(I)
          IF ( .NOT. tmpLct%Dct%Scal_cID_set ) THEN
             CALL ListCont_Find ( ConfigList, tmpID, 1, FOUND, mskLct )
          ELSE
             CALL ListCont_Find ( ConfigList, tmpID, 0, FOUND, mskLct )
          ENDIF

          ! Error if container not found
          IF ( .NOT. FOUND ) THEN
             MSG = 'No container with cID: ' // TRIM(tmpLct%Dct%cName)
             CALL HCO_ERROR ( MSG, RC)
             RETURN
          ENDIF

          ! Write out coverage.
          ! Note: If one mask has only partial coverage, retain that
          ! value! If we encounter a mask with no coverage, set coverage
          ! to zero and leave immediately. 
          IF ( (mskLct%Dct%DctType == 3) ) THEN
             IF ( mskLct%Dct%Dta%Cover == -1 ) THEN
                tmpCov = -1
             ELSEIF ( mskLct%Dct%Dta%Cover == 0 ) THEN
                tmpCov = 0
                EXIT
             ENDIF
          ENDIF
       ENDDO ! I 

       ! If tmpLct has no coverage, we can ignore this tmpLct as 
       ! it will never overwrite data of currCont
       IF ( tmpCov == 0 ) THEN
          CALL GetNextCont ( tmpLct, FLAG1 ); CYCLE
       ENDIF

       ! If we made it up to here and tmpLct has full coverage, then
       ! tmpLct has the same species ID, category, ext. nr., 
       ! and a higher (or the same) hierarchy as Lct. 

       ! If hierarchy of tmpLct is higher than Lct and this 
         ! container has total coverage over this CPU, it will always 
         ! replace all values of Lct. Hence, set targetID to -999
         ! (= ignore container) and return here.
       IF ( (tmpLct%Dct%Hier > Hier) .AND. (tmpCov==1) ) THEN
          IF ( verb ) THEN
             WRITE(MSG,*) 'Skip container ', TRIM(Lct%Dct%cName), &
                          ' because of ', TRIM(tmpLct%Dct%cName)
             CALL HCO_MSG(MSG)
          ENDIF
            
          ! Return
          targetID = -999 
          CALL HCO_LEAVE ( RC )
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
       IF ( tmpLct%Dct%Hier == Hier ) THEN 

          ! temporary flag
          sameCont = .TRUE.

          ! Check for same scale factors
          DO I = 1, SclMax
             IF (  tmpLct%Dct%Scal_cID(I) /= &
                  Lct%Dct%Scal_cID(I)     ) THEN
                sameCont = .FALSE.
                EXIT
             ENDIF
          ENDDO

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

          ! Finally, check for "same" container names. This checks the 
          ! container names ignoring the name 'tags'. 
          sameCont = Check_ContNames( tmpLct, Lct )

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
       CALL GetNextCont ( tmpLct, FLAG1 )  

    ENDDO !Loop over all entries in ConfigList (tmpLct) 

    ! Free pointers 
    tmpLct => NULL()
    mskLct => NULL()

    ! Leave w/ success
    CALL HCO_LEAVE ( RC ) 

  END SUBROUTINE Get_targetID
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
! !INTERFACE:
!
  SUBROUTINE ReadAndSplit_Line( AIR,     IU_HCO,  char1, chr1cl, &
                                char2,   chr2cl,  char3, chr3cl, &
                                char4,   chr4cl,  char5, chr5cl, &
                                char6,   chr6cl,  char7, chr7cl, &
                                char8,   chr8cl,  char9, chr9cl, &
                                int1,    int1cl,  int2,  int2cl, &
                                int3,    int3cl,  STAT,  inLine, &
                                outLine                           )
!
! !USES:
!
    USE CHARPAK_Mod,        ONLY : STRREPL, STRSPLIT
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN   )           :: AIR  
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
    INTEGER,            INTENT(IN   )           :: int1cl
    INTEGER,            INTENT(IN   )           :: int2cl
    INTEGER,            INTENT(IN   )           :: int3cl
    CHARACTER(LEN=255), INTENT(IN   ), OPTIONAL :: inLINE
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*),   INTENT(  OUT)           :: char1
    CHARACTER(LEN=*),   INTENT(  OUT)           :: char2
    CHARACTER(LEN=*),   INTENT(  OUT)           :: char3
    CHARACTER(LEN=*),   INTENT(  OUT)           :: char4
    CHARACTER(LEN=*),   INTENT(  OUT)           :: char5
    CHARACTER(LEN=*),   INTENT(  OUT)           :: char6
    CHARACTER(LEN=*),   INTENT(  OUT)           :: char7
    CHARACTER(LEN=*),   INTENT(  OUT)           :: char8
    CHARACTER(LEN=*),   INTENT(  OUT)           :: char9
    INTEGER,            INTENT(  OUT)           :: int1
    INTEGER,            INTENT(  OUT)           :: int2
    INTEGER,            INTENT(  OUT)           :: int3
    CHARACTER(LEN=255), INTENT(  OUT), OPTIONAL :: outLINE
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(INOUT)           :: STAT
!
! !REVISION HISTORY:
!  28 Aug 2013 - C. Keller: Initial version
!  11 Dec 2013 - C. Keller: Added optional arguments inLine and outLine 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER               :: N
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
       LINE = Config_ReadLine ( IU_HCO, EOF )
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
    IF ( LINE(1:1) == HCO_CMT() ) THEN
       STAT = 1
       RETURN
    ENDIF

    ! Split line into columns
    CALL STRREPL ( LINE, HCO_TAB(), HCO_SPC() )
    CALL STRSPLIT( LINE, HCO_SPC(), SUBSTR, N ) 

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
          WRITE(*,*) TRIM(LINE)
          STAT = 100
          RETURN
       ELSE
          READ( SUBSTR(chr1cl), '(a)' ) char1
       ENDIF
    ENDIF 

    ! Character 2 
    IF ( chr2cl > 0 ) THEN
       IF ( chr2cl > N ) THEN
          WRITE(*,*) TRIM(LINE)
          STAT = 100
          RETURN
       ELSE
          READ( SUBSTR(chr2cl), '(a)' ) char2
       ENDIF
    ENDIF 

    ! Character 3 
    IF ( chr3cl > 0 ) THEN
       IF ( chr3cl > N ) THEN
          WRITE(*,*) TRIM(LINE)
          STAT = 100
          RETURN
       ELSE
          READ( SUBSTR(chr3cl), '(a)' ) char3
       ENDIF
    ENDIF 

    ! Character 4 
    IF ( chr4cl > 0 ) THEN
       IF ( chr4cl > N ) THEN
          WRITE(*,*) TRIM(LINE)
          STAT = 100
          RETURN
       ELSE
          READ( SUBSTR(chr4cl), '(a)' ) char4
       ENDIF
    ENDIF 

    ! Character 5 
    IF ( chr5cl > 0 ) THEN
       IF ( chr5cl > N ) THEN
          WRITE(*,*) TRIM(LINE)
          STAT = 100
          RETURN
       ELSE
          READ( SUBSTR(chr5cl), '(a)' ) char5
       ENDIF
    ENDIF 

    ! Character 6 
    IF ( chr6cl > 0 ) THEN
       IF ( chr6cl > N ) THEN
          WRITE(*,*) TRIM(LINE)
          STAT = 100
          RETURN
       ELSE
          READ( SUBSTR(chr6cl), '(a)' ) char6
       ENDIF
    ENDIF 

    ! Character 7 
    IF ( chr7cl > 0 ) THEN
       IF ( chr7cl > N ) THEN
          WRITE(*,*) TRIM(LINE)
          STAT = 100
          RETURN
       ELSE
          READ( SUBSTR(chr7cl), '(a)' ) char7
       ENDIF
    ENDIF 

    ! Character 8 
    IF ( chr8cl > 0 ) THEN
       IF ( chr8cl > N ) THEN
          WRITE(*,*) TRIM(LINE)
          STAT = 100
          RETURN
       ELSE
          READ( SUBSTR(chr8cl), '(a)' ) char8
       ENDIF
    ENDIF 

    ! Character 9 
    IF ( chr9cl > 0 ) THEN
       IF ( chr9cl > N ) THEN
          WRITE(*,*) TRIM(LINE)
          STAT = 100
          RETURN
       ELSE
          READ( SUBSTR(chr9cl), '(a)' ) char9
       ENDIF
    ENDIF 

    ! ---------------------------------------------------------------------
    ! Read integers as specified and write them into given variables 
    ! ---------------------------------------------------------------------

    ! Integer 1 
    IF ( int1cl > 0 ) THEN
       IF ( int1cl > N ) THEN
          WRITE(*,*) TRIM(LINE)
          STAT = 100
          RETURN
       ELSE
          READ( SUBSTR(int1cl), * ) int1 
       ENDIF
    ENDIF 

    ! Integer 2 
    IF ( int2cl > 0 ) THEN
       IF ( int2cl > N ) THEN
          WRITE(*,*) TRIM(LINE)
          STAT = 100
          RETURN
       ELSE
          READ( SUBSTR(int2cl), * ) int2
       ENDIF
    ENDIF 

    ! Integer 3 
    IF ( int3cl > 0 ) THEN
       IF ( int3cl > N ) THEN
          WRITE(*,*) TRIM(LINE)
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
! !IROUTINE: ListCont_Remove 
!
! !DESCRIPTION: Subroutine ListCont\_Remove removes the passed container
! Lct from the ConfigList (it also cleans up the corresponding data
! container Lct%Dct). Lct is reset to the container IN FRONT
! of Lct.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ListCont_Remove( Lct, RC )
!
! !USES:
!
    USE HCO_DATACONT_Mod, ONLY : DataCont_Cleanup
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont), POINTER       :: Lct
    INTEGER,        INTENT(INOUT) :: RC
!
! !REVISION HISTORY: 
!  18 Sep 2013 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Pointers
    TYPE(ListCont), POINTER   :: TmpLct  => NULL()
    TYPE(ListCont), POINTER   :: PrvLct => NULL()

    ! Scalars
    CHARACTER(LEN=255)        :: MSG, LOC

    ! Enter
    LOC = 'ListCont_Remove (hco_config_mod.F90)'

    ! Special case where Lct is first container in list:
    IF ( ConfigList%Dct%cID == Lct%Dct%cID ) THEN
       ConfigList => Lct%NextCont
       Lct%NextCont => NULL()
       CALL DataCont_Cleanup ( Lct%Dct )
       DEALLOCATE ( Lct )
       Lct => ConfigList
         
    ELSE

       ! TmpLct is the temporary container
       PrvLct => ConfigList
       TmpLct  => ConfigList%NextCont

       ! Go through list until TmpLct is Lct. Also keep previous
       ! container in memory (PrvLct) 
       DO WHILE ( ASSOCIATED ( TmpLct ) ) 
      
          IF ( TmpLct%Dct%cID == Lct%Dct%cID ) EXIT

          ! Advance
          PrvLct => TmpLct
          TmpLct  => TmpLct%NextCont
       ENDDO

       ! Return w/ error if container not found
       IF ( .NOT. ASSOCIATED ( TmpLct ) ) THEN
          MSG = 'Cannot remove container ' // TRIM(Lct%Dct%cName)
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Remove container:

       ! TmpLct is equal to Lct, so get rid of it. 
       TmpLct => NULL()

       ! Link previous container to following one in list to avoid 
       ! a break in the list chain.
       PrvLct%NextCont => Lct%NextCont

       ! Detach Lct from list, then remove
       Lct%NextCont => NULL()
       CALL DataCont_Cleanup ( Lct%Dct ) 
       DEALLOCATE ( Lct )      

       ! Reset Lct to previous container in list
       Lct => PrvLct
    ENDIF

    ! Free pointer
    TmpLct => NULL()
    PrvLct => NULL()

    ! Return w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE ListCont_Remove
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Config_Cleanup 
!
! !DESCRIPTION: Subroutine Config\_Cleanup removes ConfigList. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Config_Cleanup( RemoveDct ) 
!
! !USES:
!
    USE HCO_DATACONT_Mod, ONLY : ListCont_Cleanup
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN) :: RemoveDct
!
! !REVISION HISTORY: 
!  18 Sep 2013 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Clean up Config list
    CALL ListCont_Cleanup ( ConfigList, RemoveDct )
    ConfigList => NULL()

    ! Reset internal variables to default values
    ConfigFileRead = .FALSE. 

  END SUBROUTINE Config_Cleanup
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE Get_cID( ScalID, cID, RC ) 
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN   )  :: scalID
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(  OUT)  :: cID 
!
! !INPUT/OUTPUTP PARAMETERS:
!
    INTEGER, INTENT(INOUT)  :: RC 
!
! !REVISION HISTORY: 
!  18 Sep 2013 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ListCont), POINTER  :: Lct => NULL()

    ! Scalars
    CHARACTER(LEN=255)       :: MSG, LOC
    CHARACTER(LEN= 31)       :: strID

    ! Enter
    LOC = 'Get_cID (hco_config_mod.F90)'
    cID = -999

    ! Loop over all containers 
    Lct => ConfigList
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
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    
    ! Leave w/ success
    RC = HCO_SUCCESS 

  END SUBROUTINE Get_cID
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetNextCont
!
! !DESCRIPTION: Subroutine GetNextCont is a simple wrapper routine to get 
! the next/first pointer of the ConfigList. See ListCont\_NextCont for more
! details. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GetNextCont( Lct, FLAG ) 
!
! !USES:
!
    USE HCO_DATACONT_Mod, ONLY : ListCont_NextCont
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont), POINTER       :: Lct
    INTEGER,        INTENT(INOUT) :: FLAG 
!
! !REVISION HISTORY:
!  21 Nov 2013 - C. Keller: Initialization 
!EOP
!------------------------------------------------------------------------------
!BOC

    !======================================================================
    ! GetNextCont begins here
    !======================================================================

    CALL ListCont_NextCont ( ConfigList, Lct, FLAG )

  END SUBROUTINE GetNextCont
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ListCont), POINTER :: Lct
    TYPE(ListCont), POINTER :: List 
!
! !REVISION HISTORY:
!  17 Sep 2013 - C. Keller: Initialization (update)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(ListCont), POINTER :: NewLct => NULL()

    !======================================================================
    ! ConfigList_AddCont begins here
    !======================================================================

    ! Allocate container and create data structure.
    ! The DataCont_Init call creates a new data container (type DataCont)
    ! and assigns a unique container ID (cID) to it. All HEMCO lists 
    ! (ConfigList, ReadList, EmisList) point to this container!
    ALLOCATE ( NewLct ) 
    NewLct%Dct      => NULL()
    NewLct%NextCont => NULL()
    CALL DataCont_Init ( NewLct%Dct )

    ! Connect blank container with ConfigList list. 
    NewLct%NextCont => List 
    List            => NewLct

    ! Output pointer points to the new container 
    Lct => NewLct 

  END SUBROUTINE ConfigList_AddCont
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE ScalID_Register( ScalIDs, RC ) 
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN   ) :: ScalIDs(SclMax)
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ScalIDCont), POINTER  :: NewScalIDCont => NULL()
    TYPE(ScalIDCont), POINTER  :: TmpScalIDCont => NULL() 

    ! Scalars
    INTEGER                    :: N
    LOGICAL                    :: IsInList
      
    !======================================================================
    ! ScalID_Register begins here
    !======================================================================

    ! Check for every element of ScalIDs, if this scale factor ID is
    ! already a member of ScalIDList. If not, add to it. 
    DO N = 1, SclMax
       IF ( ScalIDs(N) < 0 ) CYCLE
     
       ! Check if already in list 
       IsInList = .FALSE.
       TmpScalIDCont => ScalIDList
       DO WHILE ( ASSOCIATED(TmpScalIDCont) ) 
          IF ( TmpScalIDCont%ScalID == ScalIDs(N) ) THEN
             IsInList = .TRUE.
             EXIT
          ENDIF
          TmpScalIDCont => TmpScalIDCont%NEXT
       ENDDO

       ! Add new container w/ this scal ID to (beginning) of list 
       IF ( .NOT. IsInList ) THEN
          ALLOCATE ( NewScalIDCont ) 
          NewScalIDCont%ScalID =  ScalIDs(N)
          NewScalIDCont%NEXT   => ScalIDList
          ScalIDList           => NewScalIDCont
          NewScalIDCont        => NULL()
       ENDIF
  
       TmpScalIDCont => NULL()
    ENDDO

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ScalID_Register
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE ScalID_Cleanup()
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ScalIDCont), POINTER  :: TmpScalIDCont => NULL()
    TYPE(ScalIDCont), POINTER  :: NxtScalIDCont => NULL() 

    !======================================================================
    ! ScalID_Cleanup begins here
    !======================================================================

    ! Walk through list and remove each element
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE SpecName_Register( SpecName, RC ) 
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   ) :: SpecName 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC 
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
    TYPE(SpecNameCont), POINTER  :: NewSpecNameCont => NULL()
    TYPE(SpecNameCont), POINTER  :: TmpSpecNameCont => NULL() 
    LOGICAL                      :: IsInList

    !======================================================================
    ! SpecName_Register begins here
    !======================================================================

    ! Ignore if wildcard character. These fields will always be used!
    IF ( TRIM(SpecName) == TRIM(HCO_WCD() ) ) THEN
       RC = HCO_SUCCESS
       RETURN
    ENDIF

    ! Check if already in list 
    IsInList = .FALSE.
    TmpSpecNameCont => SpecNameList
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
       NewSpecNameCont%NEXT     => SpecNameList
       SpecNameList             => NewSpecNameCont
       NewSpecNameCont          => NULL()
    ENDIF
 
    ! Cleanup 
    TmpSpecNameCont => NULL()

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE SpecName_Register
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE SpecName_Cleanup
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
    ! Pointers
    TYPE(SpecNameCont), POINTER  :: TmpSpecNameCont => NULL()
    TYPE(SpecNameCont), POINTER  :: NxtSpecNameCont => NULL() 

    !======================================================================
    ! SpecName_Cleanup begins here
    !======================================================================

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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  FUNCTION Config_GetnSpecies() RESULT( nSpecies )
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
   
    CALL Config_GetSpecAttr( N=nSpecies, RC = THISRC ) 
 
  END FUNCTION Config_GetnSpecies
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE Config_GetSpecNames( SpecNames, nSpecies, RC ) 
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
   
    CALL Config_GetSpecAttr( N=nSpecies, SpecNames=SpecNames, RC=RC )
 
  END SUBROUTINE Config_GetSpecNames
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
  SUBROUTINE Config_GetSpecAttr( N, SpecNames, RC ) 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(INOUT)           :: N
    INTEGER,           INTENT(INOUT)           :: RC 
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*),  POINTER,       OPTIONAL    :: SpecNames(:)
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(SpecNameCont), POINTER   :: TmpSpecNameCont => NULL() 
    INTEGER                       :: AS
    CHARACTER(LEN=255), PARAMETER :: &
       LOC = 'Config_GetSpecAttr (hco_config_mod.F90)'

    !======================================================================
    ! Config_GetSpecAttr begins here
    !======================================================================

    ! Eventually allocate pointer
    IF ( PRESENT(SpecNames) ) THEN
       IF ( .NOT. ASSOCIATED(SpecNames) ) THEN
          IF ( N <= 0 ) THEN
             CALL HCO_ERROR ( 'Cannot allocate SpecNames', RC, THISLOC=LOC )
             RETURN
          ENDIF
          ALLOCATE(SpecNames(N), STAT=AS )
          IF ( AS/= 0 ) THEN
             CALL HCO_ERROR ( 'SpecNames allocation error', RC, THISLOC=LOC )
             RETURN
          ENDIF
          SpecNames(:) = ''
       ELSEIF ( SIZE(SpecNames) /= N ) THEN
          CALL HCO_ERROR ( 'SpecNames size error', RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDIF
 
    ! Init
    N = 0
  
    ! Loop over entire list. Count number of containers and eventually
    ! write out the species names. 
    TmpSpecNameCont => SpecNameList
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Config_ScalIDinUse
!
! !DESCRIPTION: Function Config\_ScalIDinUse returns TRUE if scale factor
! ID scalID is registered as scale factor to be used by any of the base 
! data (FALSE otherwise). 
!\\
!\\
! !INTERFACE:
!
  FUNCTION Config_ScalIDinUse(ScalID) RESULT( inUse )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN   )   :: ScalID
!
! !RETURN VALUE:
!
    LOGICAL                  :: inUse 
!
! !REVISION HISTORY:
!  10 Jan 2014 - C. Keller: Initialization (update)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(ScalIDCont), POINTER  :: TmpScalIDCont => NULL() 

    ! For error handling
    INTEGER            :: RC
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'Config_ScalIDinUse (hco_config_mod.F90)'

    !======================================================================
    ! Config_ScalIDinUse begins here
    !======================================================================

    ! This routine should only be called after SetReadList! If not, assume
    ! that all scale factors are actually used!
    IF ( .NOT. SetReadListCalled ) THEN
       MSG = 'You are calling Config_ScalIDinUse before SetReadList - ' // &
             'this may lead to unexpected behavior!!'
       CALL HCO_WARNING ( MSG, RC, THISLOC=LOC )
       inUse = .TRUE.
       RETURN
    ENDIF
  
    ! Init 
    inUse = .FALSE.

    TmpScalIDCont => ScalIDList
    DO WHILE ( ASSOCIATED(TmpScalIDCont) )
       IF ( TmpScalIDCont%ScalID == ScalID ) THEN
          inUse = .TRUE. 
          EXIT
       ENDIF
       TmpScalIDCont => TmpScalIDCont%NEXT
    ENDDO

    ! Return
    TmpScalIDCont => NULL()

  END FUNCTION Config_ScalIDinUse
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
    CHARACTER(LEN=31)  :: name1, name2
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
END MODULE HCO_Config_Mod
!EOM
