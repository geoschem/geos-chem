!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_diagn_mod.F90
!
! !DESCRIPTION: Module HCO\_Diagn\_mod contains routines and 
! variables to handle the HEMCO diagnostics. The HEMCO diagnostics
! consist of a flexible suite of diagnostics container organized 
! in list DiagnList. Each diagnostics container contains information 
! about the diagnostics type (extension number, emission category / 
! hierarchy, species ID), data structure (Scalar, 2D, 3D), and output
! units (mass, area, time).
!\\
!\\
! The HEMCO diagnostics module can store multiple, independent 
! diagnostics `collections`, identifiable through the assigned 
! collection number. Each collection has an output frequency assigned
! to it, as well as an output file name (prefix). All containers of
! the same collection will have the same output frequency. Currently,
! the following output frequencies are defined: 'Hourly', 'Daily',
! 'Monthly', 'Annually', 'End', 'Manual'. 
!\\
!\\
! HEMCO has three default built-in diagnostic collections: default,
! manual, and restart. These three collections become automatically
! defined during initialization of HEMCO, and diagnostic containers
! can be added to them anytime afterwards.
! The output frequency of the default collection can be specified
! in the settings section of the HEMCO configuration file 
! (DiagnFreq), along with its output file prefix (DiagnPrefix). The
! restart collection always gets an output frequency of 'End', but
! writing its content to disk can be forced at any given time using
! routine HcoDiagn\_Write (see below). The manual diagnostics has
! an output frequency of 'Manual', which means that it's content is
! never written to disk. Instead, it's fields need to be fetched
! explicitly in other routine via routine Diagn\_Get. 
!\\
!\\ 
! The public module variables HcoDiagnIDDefault, HcoDiagnIDManual, 
! and HcoDiagnRestart can be used to refer to these collections. 
! The user can also define its own collections. It is recommended 
! to do this outside of this module, e.g. at the model - HEMCO 
! interface.
!\\
!\\
! Diagnostic collections are written to disk using the routines in
! module hcoio\_diagn\_mod.F90. Routine HcoDiagn\_Write will write
! out the three built-in HEMCO collections. Other collections need 
! be written out explicitly using routine HCOIO\_Diagn\_WriteOut.
! The HEMCO option 'HcoWritesDiagn' determines if the three HEMCO
! collections are automatically written out by the HEMCO driver 
! routines (hco\_driver\_mod.F90). If HcoWritesDiagn is set to 
! FALSE, the user can freely decide when to write out the 
! diagnostics. This is useful if the HEMCO diagnostics contain
! fields that are used/filled outside of HEMCO.
!\\
!\\
! Diagnostics container are created at the beginning of a simulation
! using subroutine Diagn\_Create. During the simulation, content is
! added to the individual containers via Diagn\_Update. Diagnostics
! data is fetched using Diagn\_Get. All emissions are stored in units
! of [kg/m2] and only converted to desired output units when returning 
! the data. The container variable IsOutFormat denotes whether data
! is currently stored in output units or internal units. Variable 
! nnGetCalls counts the number of times a diagnostics is called through 
! Diagn\_Get without updating its content. This is useful if you want to 
! make sure that data is only written once per time step.
!\\
!\\
! There are two types of emission diagnostics: automatic (`AutoFill`) 
! and manual diagnostics. AutoFill diagnostics become automatically 
! filled during execution of HEMCO. AutoFill diagnostics can be at 
! species level (level 1), ExtNr level (level 2), emission category level 
! (level 3), or hierarchy level (level 4). Level 1 diagnostics write out
! the collected emissions of the specified species, level 2 diagnostics
! write out emissions for the given ExtNr only (ignoring emissions from
! all other ExtNr's), etc. 
! Manual diagnostics can represent any content. They never become filled
! automatically and all update calls (Diagn\_Update) have to be set
! manually.
!\\
!\\!\\
! Individual diagnostics are identified by its name and/or container ID.
! Both are specified when creating the diagnostics (Diagn\_Create). 
!\\
!\\
! Before adding diagnostics to a collection, the collection needs to be 
! created using subroutine DiagnCollection\_Create. The collection number 
! argument (COL) should always be specified when creating, editing or 
! obtaining a diagnostics. If this argument is omitted, the default HEMCO
! collection (HcoDiagnIDDefault) is taken. 
!\\
!\\
! !INTERFACE: 
!
MODULE HCO_Diagn_Mod 
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Arr_Mod
  USE HCO_Clock_Mod  ! Contains all the reset flags

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HcoDiagn_AutoUpdate
  PUBLIC  :: HcoDiagn_Init
  PUBLIC  :: Diagn_Create
  PUBLIC  :: Diagn_Update 
  PUBLIC  :: Diagn_Get
  PUBLIC  :: Diagn_TotalGet
  PUBLIC  :: Diagn_AutoFillLevelDefined
  PUBLIC  :: Diagn_Print
  PUBLIC  :: Diagn_DefineFromConfig
  PUBLIC  :: DiagnCont_Find
  PUBLIC  :: DiagnCollection_Create
  PUBLIC  :: DiagnCollection_Cleanup
  PUBLIC  :: DiagnCollection_Get
  PUBLIC  :: DiagnFileOpen
  PUBLIC  :: DiagnFileGetNext
  PUBLIC  :: DiagnFileClose
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: DiagnList_Cleanup 
  PRIVATE :: DiagnCont_Init
  PRIVATE :: DiagnCont_PrepareOutput
  PRIVATE :: DiagnCont_Link_2D
  PRIVATE :: DiagnCont_Link_3D
  PRIVATE :: DiagnCont_Cleanup
  PRIVATE :: DiagnCollection_DefineID 
  PRIVATE :: DiagnCollection_Find
  PRIVATE :: Diagn_UpdateDriver

  INTERFACE Diagn_Update
     MODULE PROCEDURE Diagn_UpdateSP 
     MODULE PROCEDURE Diagn_UpdateDP 
  END INTERFACE
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller   - Initialization
!  08 Jul 2014 - R. Yantosca - Now use F90 free-format indentation  
!  08 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  01 Aug 2014 - C. Keller   - Added manual output frequency
!  12 Aug 2014 - C. Keller   - Added cumulative sum option
!  09 Jan 2015 - C. Keller   - Added diagnostics collections
!  03 Apr 2015 - C. Keller   - Now tie output frequency to collection instead
!                              of individual diagnostic containers.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  !------------------------------------------------------------------------
  ! DiagnCont: Diagnostics container derived type declaration
  !------------------------------------------------------------------------
  TYPE, PUBLIC :: DiagnCont
     CHARACTER(LEN= 31)          :: cName          ! Cont. name
     CHARACTER(LEN=255)          :: long_name      ! ncdf long_name attribute 
     INTEGER                     :: cID            ! Cont. ID
     INTEGER                     :: ExtNr          ! Extension #
     INTEGER                     :: Cat            ! Category 
     INTEGER                     :: Hier           ! Hierarchy
     INTEGER                     :: HcoID          ! HEMCO species ID
     INTEGER                     :: AutoFill       ! fill automatically? 
     INTEGER                     :: SpaceDim       ! Space dimension (1-3) 
     REAL(sp)                    :: Scalar         ! 1D scalar 
     TYPE(Arr2D_SP),     POINTER :: Arr2D          ! 2D array
     TYPE(Arr3D_SP),     POINTER :: Arr3D          ! 3D array
     REAL(sp)                    :: Total          ! Diagnostics total 
     LOGICAL                     :: DtaIsPtr       ! Is data just a pointer?
     INTEGER                     :: LevIdx         ! Level index to be used 
     CHARACTER(LEN= 31)          :: OutUnit        ! Output unit 
     INTEGER                     :: AreaFlag       ! 2=per area, 3=per volume, 0 otherwise 
     REAL(hp)                    :: AreaScal       ! Scale factor for area
     REAL(hp)                    :: MassScal       ! Scale factor for mass
     REAL(hp)                    :: ScaleFact      ! Uniform scale factor 
     INTEGER                     :: TimeAvg        ! Scale flag for time unit 
     INTEGER                     :: Counter        ! time steps since 
                                                   ! last output
     CHARACTER(LEN= 31)          :: AvgName        ! Output averaging operation 
     INTEGER                     :: AvgFlag        ! Averaging flag for 
                                                   !  non-standard units
     INTEGER                     :: LastUpdateID   ! Last update time
     INTEGER                     :: nnGetCalls     ! # of Diagn_Get calls w/o update 
     LOGICAL                     :: IsOutFormat    ! Data is in output format?
     INTEGER                     :: CollectionID   ! Collection diagnostics belongs to
     TYPE(DiagnCont),    POINTER :: NextCont       ! Ptr to next item in list
  END TYPE DiagnCont

  !------------------------------------------------------------------------
  ! Diagnostcs collection derived type.
  ! DiagnList      : Linked list with all diagnostics container of
  !                  this collection.
  ! nnDiag         : Number of diagnostics in this collection.
  ! AF_LevelDefined: Set to true if there is at least one autofill 
  !                  diagnostics at the given level (1-4).
  ! PREFIX         : Prefix to be used for diagnostics output file name.
  ! WriteFreq      : Output write frequency
  ! ResetFlag      : Reset flag of this collection. Will be determined 
  !                  based on WriteFreq.
  ! NX, NY, NZ     : Grid dimensions.
  ! TS             : Time step. This is only of relevance for emission 
  !                  diagnostics that are internally converted from
  !                  kg/m2/s to kg/m2.
  ! AREA_M2        : Surface grid box areas. May be required for unit 
  !                  conversions.
  !------------------------------------------------------------------------
  TYPE DiagnCollection
     TYPE(DiagnCont),       POINTER :: DiagnList          => NULL()
     INTEGER                        :: nnDiagn            =  0
     LOGICAL                        :: AF_LevelDefined(4) =  .FALSE.
     INTEGER                        :: CollectionID       = -1
     CHARACTER(LEN=255)             :: PREFIX             =  ''
     CHARACTER(LEN=31)              :: WriteFreq          = ''
     INTEGER                        :: ResetFlag          =  ResetFlagManually 
     INTEGER                        :: NX                 =  0
     INTEGER                        :: NY                 =  0
     INTEGER                        :: NZ                 =  0
     REAL(sp)                       :: TS                 =  0       ! Time step
     REAL(hp),              POINTER :: AREA_M2(:,:)       => NULL()
     TYPE(DiagnCollection), POINTER :: NextCollection     => NULL()
  END TYPE DiagnCollection

  ! Pointer to beginning of collections linked list 
  TYPE(DiagnCollection),  POINTER :: Collections => NULL()

  ! HEMCO diagnostic collection IDs. Used to identify collections. Its values
  ! will be set upon creation of the collection. 
  INTEGER, PUBLIC     :: HcoDiagnIDDefault = -999 
  INTEGER, PUBLIC     :: HcoDiagnIDRestart = -999
  INTEGER, PUBLIC     :: HcoDiagnIDManual  = -999

  ! Total number of collections in collection linked list
  INTEGER             :: nnCollections     = 0
!
! !DEFINED PARAMETERS:
!
  ! Parameter for averaging and summing non-standard data
  ! AvgFlagMean    : calculates the arithmetic mean
  ! AvgFlagSum     : calculates the sum, resets after every writeout
  ! AvgFlagCumulSum: calculates the cumulative sum, never resets.
  ! AvgFlagInst    : uses the instantaneous value, overwrites existing 
  INTEGER, PARAMETER             :: AvgFlagMean     = 1
  INTEGER, PARAMETER             :: AvgFlagSum      = 2
  INTEGER, PARAMETER             :: AvgFlagCumulSum = 3
  INTEGER, PARAMETER             :: AvgFlagInst     = 4

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoDiagn_autoupdate
!
! !DESCRIPTION: Subroutine HCODIAGN\_AUTOUPDATE updates the AutoFill
! diagnostics at species level. This routine should be called after
! running HEMCO core and all extensions. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoDiagn_AutoUpdate( am_I_Root, HcoState, RC ) 
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_GetHcoID
    USE HCO_STATE_MOD, ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  19 Dec 2013 - C. Keller   - Initial version 
!  11 Jun 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  11 Jun 2014 - R. Yantosca - Now use F90 freeform indentation
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)        :: MSG, LOC
    INTEGER                   :: I, tmpID
    REAL(hp), POINTER         :: Arr3D(:,:,:) => NULL()
    REAL(hp), POINTER         :: Arr2D(:,:)   => NULL()

    !=================================================================
    ! HCODIAGN_AUTOUPDATE begins here!
    !=================================================================
    
    ! Init 
    LOC = 'HCODIAGN_AUTOUPDATE (hco_diagn_mod.F90)'
    RC  = HCO_SUCCESS
    
    ! ================================================================
    ! AutoFill diagnostics: only write diagnostics at species level
    ! (level 1). Higher level diagnostics have been written in the
    ! respective subroutines (hco_calc & extension modules). 
    ! ================================================================
    DO I = 1, HcoState%nSpc
       IF ( ASSOCIATED(HcoState%Spc(I)%Emis) ) THEN
          IF ( ASSOCIATED(HcoState%Spc(I)%Emis%Val) ) THEN
             Arr3D => HcoState%Spc(I)%Emis%Val
             CALL Diagn_Update( am_I_Root,            &
                                ExtNr    = -1,        &
                                Cat      = -1,        &
                                Hier     = -1,        &
                                HcoID    = I,         &
                                AutoFill = 1,         &
                                Array3D  = Arr3D,     &
                                COL      = -1,        &
                                RC       = RC          ) 
             IF ( RC/= HCO_SUCCESS ) RETURN 
             Arr3D => NULL() 
          ENDIF
       ENDIF
    ENDDO
    
    ! Return
    RC = HCO_SUCCESS
    
  END SUBROUTINE HcoDiagn_AutoUpdate
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HcoDiagn_Init
!
! !DESCRIPTION: Subroutine HCODIAGN\_INIT initializes the three built-in 
! HEMCO diagnostic collections: default, restart, and manual. The
! identification ID of each collection is written into public variable
! HcoDiagnIDDefault, HcoDiagnIDRestart, and HcoDiagnIDManual, respectively.
! Those are used to easily refer to one of the diagnostics when adding
! fields ('containers') to a collection or fetching it's content. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HcoDiagn_Init( am_I_Root, HcoState, RC ) 
!
! !USES:
!
    USE HCO_STATE_MOD,   ONLY : HCO_GetHcoID
    USE HCO_STATE_MOD,   ONLY : HCO_State
    USE HCO_ExtList_Mod, ONLY : GetExtOpt
    USE HCO_ExtList_Mod, ONLY : CoreNr 
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER        :: HcoState   ! HEMCO state object 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  03 Apr 2015 - C. Keller   - Initial version 
!  10 Apr 2015 - C. Keller   - Now create diagnostics based on entries 
!                              in the HEMCO diagnostics definition file.
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: CollectionID
    LOGICAL             :: FOUND
    CHARACTER(LEN=15)   :: WriteFreq
    CHARACTER(LEN=255)  :: LOC, DiagnPrefix

    !=================================================================
    ! HCODIAGN_INIT begins here!
    !=================================================================
    
    ! Init 
    LOC = 'HCODIAGN_INIT (hco_diagn_mod.F90)'

    ! ------------------------------------------------------------------
    ! Default diagnostics
    ! ------------------------------------------------------------------

    ! Output frequency. Try to read from configuration file. 
    CALL GetExtOpt ( CoreNr, 'DiagnFreq', OptValChar=WriteFreq, &
                     FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       WriteFreq = 'Daily'
    ENDIF

    ! Force to 'Always' in ESMF environment to make sure that
    ! diagnostics are passed to MAPL HISTORY every time.
#if defined ( ESMF_ )
    WriteFreq = 'Always'
#endif

    ! Try to get prefix from configuration file
    CALL GetExtOpt ( CoreNr, 'DiagnPrefix', OptValChar=DiagnPrefix, &
                     FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       DiagnPrefix = 'HEMCO_Diagnostics_' // TRIM(WriteFreq)
    ENDIF

    CALL DiagnCollection_Create( am_I_Root,                             &
                                 NX        = HcoState%NX,               &
                                 NY        = HcoState%NY,               &
                                 NZ        = HcoState%NZ,               &
                                 TS        = HcoState%TS_EMIS,          &
                                 AM2       = HcoState%Grid%AREA_M2%Val, &
                                 COL       = CollectionID,              & 
                                 PREFIX    = TRIM(DiagnPrefix),         &
                                 WriteFreq = TRIM(WriteFreq),           & 
                                 RC        = RC                          )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Pass this collection ID to fixed variable for easy further 
    ! reference to this collection
    HcoDiagnIDDefault = CollectionID
 
    ! ------------------------------------------------------------------
    ! HEMCO restart 
    ! ------------------------------------------------------------------
#if defined ( ESMF_ )
    WriteFreq = 'Always'
#else
    WriteFreq = 'End'
#endif
    CALL DiagnCollection_Create( am_I_Root,                             &
                                 NX        = HcoState%NX,               &
                                 NY        = HcoState%NY,               &
                                 NZ        = HcoState%NZ,               &
                                 TS        = HcoState%TS_EMIS,          &
                                 AM2       = HcoState%Grid%AREA_M2%Val, &
                                 COL       = CollectionID,              & 
                                 PREFIX    = 'HEMCO_restart',           &
                                 WriteFreq = TRIM(WriteFreq),           & 
                                 RC        = RC                          )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Pass this collection ID to fixed variable for easy further 
    ! reference to this collection
    HcoDiagnIDRestart = CollectionID
 
    ! ------------------------------------------------------------------
    ! Manual diagnostics
    ! ------------------------------------------------------------------
#if defined ( ESMF_ )
    WriteFreq = 'Always'
#else
    WriteFreq = 'Manual'
#endif
    CALL DiagnCollection_Create( am_I_Root,                             &
                                 NX        = HcoState%NX,               &
                                 NY        = HcoState%NY,               &
                                 NZ        = HcoState%NZ,               &
                                 TS        = HcoState%TS_EMIS,          &
                                 AM2       = HcoState%Grid%AREA_M2%Val, &
                                 COL       = CollectionID,              & 
                                 PREFIX    = 'HEMCO_manual',            &
                                 WriteFreq = TRIM(WriteFreq),           & 
                                 RC        = RC                          )
    IF ( RC /= HCO_SUCCESS ) RETURN
 
    ! Pass this collection ID to fixed variable for easy further 
    ! reference to this collection
    HcoDiagnIDManual  = CollectionID

    ! ------------------------------------------------------------------
    ! Now that collections are defined, add diagnostics specified in the
    ! HEMCO diagnostics definition file. The latter can be specified in 
    ! the HEMCO configuration file. These diagnostics are all written
    ! into the default HEMCO collection.
    ! ------------------------------------------------------------------
    CALL Diagn_DefineFromConfig( am_I_Root, HcoState, RC ) 
    IF ( RC /= HCO_SUCCESS ) RETURN
 
    ! Return w/ success
    RC = HCO_SUCCESS
    
  END SUBROUTINE HcoDiagn_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_DefineFromConfig
!
! !DESCRIPTION: Subroutine Diagn\_DefineFromConfig defines HEMCO 
! diagnostic containers as specified in the diagnostics input file.
!\\
!\\
! This routine reads information from a HEMCO diagnostics definition 
! file (specified in the main HEMCO configuration file) and creates
! HEMCO diagnostic containers for each entry of the diagnostics
! definition file. Each line of the diagnostics definition file 
! represents a diagnostics container and is expected to consist of
! 7 entries: container name (character), HEMCO species (character), 
! extension number (integer), emission category (integer), emission 
! hierarchy (integer), space dimension (2 or 3), output unit 
! (character).
!\\
!\\
! The HEMCO setting 'DiagnFile' can be used to specify a diagnostics
! file. This setting should be placed in the settings section of the
! HEMCO configuration file.
!\\
!\\
! If argument `Add2MaplExp` is set to true, the diagnostics field
! defined in the diagnostics definition file are not added to the
! HEMCO diagnostics collection (yet), but rather added to the MAPL
! export state. This is useful in an ESMF environment to automate
! the coupling of HEMCO diagnostics, e.g. subroutine 
! Diagn\_DefineFromConfig can be called during SetServices to make
! sure that all diagnostic fields defined in DiagnFile have a 
! corresponding Export state object (and can thus be written out
! via the MAPL History component).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_DefineFromConfig( am_I_Root, HcoState, RC )
!
! !USES:
!
    USE HCO_CharTools_Mod
    USE CHARPAK_Mod,       ONLY : STRREPL, STRSPLIT
    USE inquireMod,        ONLY : findFreeLUN
    USE HCO_STATE_MOD,     ONLY : HCO_GetHcoID
    USE HCO_STATE_MOD,     ONLY : HCO_State
    USE HCO_ExtList_Mod,   ONLY : CoreNr, GetExtOpt
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root   ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState    ! HEMCO state object
    INTEGER,          INTENT(INOUT)           :: RC          ! Failure or success
!
! !REVISION HISTORY: 
!  10 Apr 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: N,      LUN
    LOGICAL             :: EOF
    CHARACTER(LEN=31)   :: cName,  SpcName, OutUnit
    INTEGER             :: HcoID,  ExtNr,   Cat, Hier, SpaceDim
    CHARACTER(LEN=255)  :: LOC,    MSG

    !=================================================================
    ! Diagn_DefineFromConfig begins here!
    !=================================================================
    
    ! Init 
    LOC = 'Diagn_DefineFromConfig (hco_diagn_mod.F90)'

    ! Load DiagnFile into buffer
    CALL DiagnFileOpen( am_I_Root, LUN, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If defined, sequentially get all entries 
    IF ( LUN > 0 ) THEN

       ! Do for every line
       DO

          ! Get next line
          CALL DiagnFileGetNext( am_I_Root, LUN, cName, &
             SpcName, ExtNr, Cat, Hier, SpaceDim, OutUnit, EOF, RC ) 
          IF ( RC /= HCO_SUCCESS ) RETURN

          ! Leave here if end of file
          IF ( EOF ) EXIT

          ! Get HEMCO species ID. Skip entry if HEMCO ID not
          ! defined for this species
          HcoID = HCO_GetHcoID( TRIM(SpcName), HcoState ) 
          IF ( HcoID <= 0 ) CYCLE
   
          ! ------------------------------------------------------------------
          ! Add it to the HEMCO diagnostics collection
          ! ------------------------------------------------------------------
          CALL Diagn_Create( am_I_Root,                     &
                             HcoState  = HcoState,          &
                             cName     = cName,             &
                             HcoID     = HcoID,             &  
                             ExtNr     = ExtNr,             &  
                             Cat       = Cat,               &  
                             Hier      = Hier,              &  
                             SpaceDim  = SpaceDim,          &  
                             OutUnit   = OutUnit,           &  
                             AutoFill  = 1,                 &  
                             COL       = HcoDiagnIDDefault, &
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN
  
       ENDDO

       ! Close file
       CALL DiagnFileClose ( LUN )

    ENDIF ! LUN > 0
 
    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Diagn_DefineFromConfig
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_Create
!
! !DESCRIPTION: Subroutine Diagn\_Create creates a new diagnostics. This
! routine takes the following input arguments:
!\begin{itemize} 
!\item am\_I\_Root: is this the root CPU?
!\item cName: distinct diagnostics (container) name. 
!\item long\_name: long\_name attribute used for netCDF output. 
!\item ExtNr: emissions extension number. 
!\item Cat: emissions category. 
!\item Hier: emissions  hierarchy.
!\item HcoID: HEMCO species ID of diagnostics species.
!\item SpaceDim: spatial dimension: 1 (scalar), 2 (lon-lat), 
!      or 3 (lon-lat-lev). 
!\item OutUnit: output unit. Emissions will be converted to this unit.
!      Conversion factors will be determined using the HEMCO unit
!      module (see HCO\_UNITS\_Mod.F90). No unit conversions will be
!      performed if the argument OutOper is set (see below). 
!\item HcoState: HEMCO state object. Used to determine the species
!      properties if any of arguments MW\_g, EmMW\_g or MolecRatio
!      is missing.
!\item OutOper: output operation for non-standard units. If this 
!      argument is used, the specified operation is performed and all
!      unit specifications are ignored. Can be one of 'Mean', 'Sum', 
!      'CumulSum', or 'Instantaneous'.
!\item AutoFill: containers with an AutoFill flag of 1 will be auto-
!      matically updated by the HEMCO standard diagnostics calls
!      (e.g. in hco\_calc\_mod.F90). If set to 0, the diagnostics 
!      updates have to be set manually.
!\item Trgt2D: 2D target array. If specified, the diagnostics array
!      will point to this data. This disables all time averaging, 
!      unit conversions, etc., and the data will be written to disk
!      as is.
!\item Trgt3D: as Trgt2D, but for 3D data. 
!\item MW\_g: species molecular weight. Used to determine unit
!      conversion factors. Not needed for target containers or if
!      argument OutOper is specified. Can be omitted if HcoState is
!      given.
!\item EmMW\_g: Molecular weight of emitted species. Used to determine
!      unit conversion factors. Not needed for target containers or if
!      argument OutOper is specified. Can be omitted if HcoState is
!      given.
!\item MolecRatio: Molecules of species per emitted molecule. Used to 
!      determine unit conversion factors. Not needed for target 
!      containers or if argument OutOper is specified. Can be omitted
!      if HcoState is given.
!\item ScaleFact: constant scale factor. If provided, the diagnostics
!      are scaled uniformly by this value before outputting. Will be
!      applied on top of any other unit conversions. Does not work on 
!      data pointers.
!\item cID: assigned container ID. Useful for later reference to this
!      diagnostics container.
!\item RC: HEMCO return code.
!\end{itemize} 
!
! !INTERFACE:
!
  SUBROUTINE Diagn_Create( am_I_Root, cName,      HcoState,   &
                           ExtNr,     Cat,        Hier,       &
                           HcoID,     SpaceDim,   OutUnit,    &  
                           OutOper,   LevIdx,     AutoFill,   &
                           Trgt2D,    Trgt3D,     MW_g,       &
                           EmMW_g,    MolecRatio, ScaleFact,  &
                           cID,       RC,    COL, OkIfExist,  &
                           long_name                           )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
    USE HCO_Unit_Mod,  ONLY : HCO_Unit_GetMassScal
    USE HCO_Unit_Mod,  ONLY : HCO_Unit_GetAreaScal
    USE HCO_Unit_Mod,  ONLY : HCO_Unit_GetTimeScal
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root     ! Root CPU?
    CHARACTER(LEN=*), INTENT(IN   )           :: cName         ! Diagnostics name 
    CHARACTER(LEN=*), INTENT(IN   )           :: OutUnit       ! Output units
    INTEGER,          INTENT(IN   ), OPTIONAL :: SpaceDim      ! Spatial dimension 
    INTEGER,          INTENT(IN   ), OPTIONAL :: ExtNr         ! Extension #    
    INTEGER,          INTENT(IN   ), OPTIONAL :: Cat           ! Category 
    INTEGER,          INTENT(IN   ), OPTIONAL :: Hier          ! Hierarchy 
    INTEGER,          INTENT(IN   ), OPTIONAL :: HcoID         ! HEMCO species ID 
    TYPE(HCO_State),  POINTER,       OPTIONAL :: HcoState      ! HEMCO state obj.
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: OutOper       ! Output operation 
    INTEGER,          INTENT(IN   ), OPTIONAL :: LevIdx        ! Level index to use 
    INTEGER,          INTENT(IN   ), OPTIONAL :: AutoFill      ! 1=fill auto.;0=don't
    REAL(sp),         INTENT(IN   ), OPTIONAL :: Trgt2D(:,:)   ! 2D target data
    REAL(sp),         INTENT(IN   ), OPTIONAL :: Trgt3D(:,:,:) ! 3D target data
    REAL(hp),         INTENT(IN   ), OPTIONAL :: MW_g          ! species MW (g/mol) 
    REAL(hp),         INTENT(IN   ), OPTIONAL :: EmMW_g        ! emission MW (g/mol)
    REAL(hp),         INTENT(IN   ), OPTIONAL :: MolecRatio    ! molec. emission ratio
    REAL(hp),         INTENT(IN   ), OPTIONAL :: ScaleFact     ! uniform scale factor 
    INTEGER,          INTENT(IN   ), OPTIONAL :: COL           ! Collection number 
    INTEGER,          INTENT(IN   ), OPTIONAL :: cID           ! Container ID 
    LOGICAL,          INTENT(IN   ), OPTIONAL :: OkIfExist     ! Is it ok if already exists? 
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: long_name     ! long name attribute 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC            ! Return code 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller - Initialization
!  05 Mar 2015 - C. Keller - container ID can now be set by the user
!  31 Mar 2015 - C. Keller - added argument OkIfExist
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(DiagnCont),       POINTER :: ThisDiagn => NULL()
    TYPE(DiagnCont),       POINTER :: TmpDiagn  => NULL()
    TYPE(DiagnCollection), POINTER :: ThisColl => NULL()

    ! Scalars
    CHARACTER(LEN=255)             :: LOC, MSG
    INTEGER                        :: PS, Flag
    REAL(hp)                       :: Scal
    REAL(hp)                       :: MWg, EmMWg, MolR 
    LOGICAL                        :: ForceMean, FOUND

    !======================================================================
    ! Diagn_Create begins here!
    !======================================================================

    ! Init
    LOC = 'Diagn_Create (hco_diagn_mod.F90)'
    CALL DiagnCollection_DefineID( PS, RC, COL=COL, &
                                   InUse=FOUND, ThisColl=ThisColl )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Error if collection does not exist
    IF ( .NOT. FOUND ) THEN
       WRITE(MSG,*) 'Cannot create diagnostics ', TRIM(cName), &
                    ' - collection does not exist: ', PS
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
    ENDIF

    !----------------------------------------------------------------------
    ! Check if diagnostics already exists 
    !----------------------------------------------------------------------
    IF ( PRESENT(OkIfExist) ) THEN
       IF ( OkIfExist ) THEN
          IF ( PRESENT(cID) ) THEN
             CALL DiagnCont_Find( cID, -1, -1, -1, -1, &
                           '', -1, FOUND, TmpDiagn, COL=PS )
          ELSE
             CALL DiagnCont_Find( -1, -1, -1, -1, -1, &
                           TRIM(cName), -1, FOUND, TmpDiagn, COL=PS )
          ENDIF
          TmpDiagn => NULL()

          ! Exit if found
          IF ( FOUND ) THEN
             IF ( HCO_IsVerb(2) ) THEN
                WRITE(MSG,*) 'Diagnostics already exists - ', &
                             'will not added again: ', TRIM(cName) 
                CALL HCO_MSG ( MSG )
             ENDIF
             RC = HCO_SUCCESS
             RETURN
          ENDIF
       ENDIF
    ENDIF

    !----------------------------------------------------------------------
    ! Initalize diagnostics container.
    !----------------------------------------------------------------------
    CALL DiagnCont_Init( ThisDiagn )

    !----------------------------------------------------------------------
    ! Pass input variables
    !----------------------------------------------------------------------
    ThisDiagn%cName   = cName
    ThisDiagn%OutUnit = TRIM(OutUnit)

    ! Optional arguments. If not provided, use default values set in
    ! DiagnCont_Init
    IF ( PRESENT(ExtNr)    ) ThisDiagn%ExtNr    = ExtNr
    IF ( PRESENT(Cat  )    ) ThisDiagn%Cat      = Cat
    IF ( PRESENT(Hier )    ) ThisDiagn%Hier     = Hier
    IF ( PRESENT(HcoID)    ) ThisDiagn%HcoID    = HcoID
    IF ( PRESENT(SpaceDim) ) ThisDiagn%SpaceDim = SpaceDim
    IF ( PRESENT(LevIdx)   ) ThisDiagn%LevIdx   = LevIdx
    IF ( PRESENT(AutoFill) ) ThisDiagn%AutoFill = AutoFill 

    ! long_name attribute. Defaults to container name
    IF ( PRESENT(long_name) ) THEN
       ThisDiagn%long_name = TRIM(long_name)
    ELSE
       ThisDiagn%long_name = TRIM(cName)
    ENDIF

    !----------------------------------------------------------------------
    ! Eventually link to data array. This will disable all time averaging,
    ! unit conversions, etc. (data will just be returned as is). 
    !----------------------------------------------------------------------
    IF ( PRESENT(Trgt2D) ) THEN
       CALL DiagnCont_Link_2D( am_I_Root, ThisDiagn, ThisColl, Trgt2D, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF
    IF ( PRESENT(Trgt3D) ) THEN
       CALL DiagnCont_Link_3D( am_I_Root, ThisDiagn, ThisColl, Trgt3D, RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    ! Update module variable AF_LevelDefined. For all AutoFill diagnostics,
    ! we store whether or not there is (at least one) diagnostics container
    ! defined at species level, ExtNr level, etc. 
    IF ( ThisDiagn%AutoFill == 1 ) THEN

       ! At species level: no ExtNr defined
       IF ( ThisDiagn%ExtNr < 0 ) THEN
          ThisColl%AF_LevelDefined(1) = .TRUE.
     
       ! At ExtNr level: no category defined
       ELSEIF ( ThisDiagn%Cat < 0 ) THEN
          ThisColl%AF_LevelDefined(2) = .TRUE.

       ! At category level: no hierarchy defined
       ELSEIF ( ThisDiagn%Hier < 0 ) THEN
          ThisColl%AF_LevelDefined(3) = .TRUE.

       ! At hierarchy level: all defined
       ELSE
          ThisColl%AF_LevelDefined(4) = .TRUE.
       ENDIF

    ENDIF

    !----------------------------------------------------------------------
    ! Determine scale factors to be applied to data. These values are
    ! determined from the specified output unit, assuming that the original
    ! HEMCO data is in kg/m2/s. If the optional argument OutOper is set,
    ! the output unit is ignored and the specified operation ('mean' or 
    ! 'sum') is performed. This is particular useful for data with 
    ! non-standard units, e.g. unitless data. 
    ! Don't need to be done for pointer data, which ignores all time 
    ! averaging, unit conversions, etc.
    !----------------------------------------------------------------------

    ! Uniform scale factor
    IF ( PRESENT(ScaleFact) ) THEN
       IF ( ThisDiagn%DtaIsPtr ) THEN
          MSG = 'Cannot use scale factor on diagnostics that '// &
                'are pointers to other data: '//TRIM(cName)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       ENDIF
       IF ( TRIM(OutOper) == 'CumulSum' ) THEN
          MSG = 'Cannot use scale factor on diagnostics that '// &
                'are cumulative sums: '//TRIM(cName)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       ENDIF
       ThisDiagn%ScaleFact = ScaleFact
    ENDIF

    ! Unit conversion factors don't need be defined for pointers
    IF ( ThisDiagn%DtaIsPtr ) THEN

       ! Pointer diagnostics are always instantaneous
       ThisDiagn%AvgName = 'Instantaneous' 

    ! Unit conversion factors for containers that are not pointers 
    ELSE

       ! Enforce specified output operator 
       IF ( PRESENT(OutOper) ) THEN

          ! Pass to diagnostics
          ThisDiagn%AvgName = TRIM(OutOper)

          ! Set flag accordingly
          IF ( TRIM(OutOper) == 'Mean' ) THEN
             ThisDiagn%AvgFlag = AvgFlagMean
          ELSEIF ( TRIM(OutOper) == 'Sum' ) THEN
             ThisDiagn%AvgFlag = AvgFlagSum
          ELSEIF ( TRIM(OutOper) == 'CumulSum' ) THEN
             ThisDiagn%AvgFlag = AvgFlagCumulSum
          ELSEIF ( TRIM(OutOper) == 'Instantaneous' ) THEN
             ThisDiagn%AvgFlag = AvgFlagInst
          ELSE
             MSG = 'Illegal output operator: ' // TRIM(OutOper)
             MSG = TRIM(MSG) // '. Allowed are `Mean`, `Sum`, '// &
                   '`CumulSum`, `Instantaneous`.'
             MSG = TRIM(MSG) // ' (' // TRIM(cName) // ')'
             CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          ENDIF
 
       ! If OutOper is not set, determine scale factors from output unit:
       ELSE
  
          ! Will calculate the mean 
          ThisDiagn%AvgName = 'mean'

          !----------------------------------------------------------------
          ! Scale factor for mass. This determines the scale factor from 
          ! HEMCO mass unit (kg) to the desired output unit. 
          ! HCO_UNIT_MassCal returns the mass scale factor from OutUnit to
          ! HEMCO unit, hence need to invert this value!
          !----------------------------------------------------------------
          IF ( .NOT. PRESENT(MW_g)       .OR. &
               .NOT. PRESENT(EmMW_g)     .OR. &
               .NOT. PRESENT(MolecRatio)       ) THEN
             IF ( .NOT. PRESENT(HcoState) ) THEN
                MSG = 'Not all species properties (MW, EmMW, MolecRatio) '// &
                      'defined for diagnostics '//TRIM(cName)//'. Please '// &
                      'specify OutOper, HcoState or species properties!.'
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
             ELSE
                MWg   = HcoState%Spc(HcoID)%MW_g
                EmMWg = HcoState%Spc(HcoID)%EmMW_g
                MolR  = HcoState%Spc(HcoID)%MolecRatio
             ENDIF
          ELSE
                MWg   = MW_g
                EmMWg = EmMW_g
                MolR  = MolecRatio
          ENDIF 

          CALL HCO_UNIT_GetMassScal(                              &
                    unt         = OutUnit,                        &
                    MW_IN       = MWg,                            &
                    MW_OUT      = EmMWg,                          &
                    MOLEC_RATIO = MolR,                           &
                    Scal        = Scal                             )
          IF ( Scal <= 0.0_hp ) THEN
             MSG = 'Cannot find mass scale factor for unit '//TRIM(OutUnit)
             CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
          ThisDiagn%MassScal = 1.0_hp / Scal
      
          !----------------------------------------------------------------
          ! Scale factor for area. This determines the scale factor from 
          ! HEMCO area unit (m2) to the desired output unit. 
          ! HCO_UNIT_AreaScal returns the area scale factor from OutUnit to
          ! HEMCO unit, hence need to invert this value!
          !----------------------------------------------------------------
          CALL HCO_UNIT_GetAreaScal( OutUnit, Scal, Flag )
          ThisDiagn%AreaFlag = Flag 
          IF ( Flag > 0 ) THEN
             ThisDiagn%AreaScal = 1.0_hp / Scal
          ENDIF
      
          !----------------------------------------------------------------
          ! Determine the normalization factors applied to the diagnostics 
          ! before they are written out. Diagnostics are always stored 
          ! internally in units of kg/m2, and the following flags make sure 
          ! that they are normalized by the desired time interval, e.g. to 
          ! get units of per second, per hour, etc.
          !----------------------------------------------------------------
      
          ! HCO_UNIT_GetTimeScal returns 1.0 for units of per second, 1/3600
          ! for per hour, etc. Returns -999.0 if no time unit could be found.
          CALL HCO_UNIT_GetTimeScal( OutUnit, 1, 2001, Scal, Flag )
          Scal = 1.0_dp / Scal
        
          ! No time unit found: don't enable any switch
          IF ( Scal < 0.0_dp ) THEN 
             ! Nothing to do
         
          ! Normalize by seconds
          ELSEIF ( Scal == 1.0_dp ) THEN
             ThisDiagn%TimeAvg = 1 
                  
          ! Normalize by hours 
          ELSEIF ( Scal == 3600.0_dp ) THEN
             ThisDiagn%TimeAvg = 2 
      
          ! Normalize by days 
          ELSEIF ( Scal == 86400.0_dp ) THEN
             ThisDiagn%TimeAvg = 3 
      
          ! Normalize by months 
          ELSEIF ( Scal == 2678400.0_dp ) THEN
             ThisDiagn%TimeAvg = 4 
      
          ! Normalize by years 
          ELSEIF ( Scal == 31536000.0_dp ) THEN
             ThisDiagn%TimeAvg = 5 
      
          ! Error otherwise
          ELSE
             MSG = 'Cannot determine time normalization: '//TRIM(OutUnit)
             CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF  
       ENDIF ! OutOper not set
    ENDIF ! .NOT. DtaIsPtr
   
    !-----------------------------------------------------------------------
    ! Make sure that there is no other diagnostics with this name 
    !-----------------------------------------------------------------------
    CALL DiagnCont_Find( -1, -1, -1, -1, -1, &
                        Trim(cName), -1, FOUND, TmpDiagn, COL=PS )
    IF ( FOUND ) THEN
       MSG = 'There is already a diagnostics with this name: ' // TRIM(cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Set container ID (if defined). There must not be two containers with
    ! the same container ID. 
    !-----------------------------------------------------------------------
    IF ( PRESENT(cID) ) THEN
       IF ( cID > 0 ) THEN

          ! Check if there is already a diagnostics with this container ID.
          TmpDiagn => NULL()
          CALL DiagnCont_Find( cID, -1, -1, -1, -1, '', -1, FOUND, TmpDiagn, COL=PS )
          IF ( FOUND ) THEN
             WRITE(MSG,*) 'Diagnostics ', TRIM(TmpDiagn%cName), ' already has ID ', &
                cID, ' - cannot create diagnostics ', TRIM(cName)
             
             CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF

          ! Set container ID
          ThisDiagn%cID = cID 
       ENDIF
    ENDIF

    !-----------------------------------------------------------------------
    ! Add to diagnostics list of this collection. 
    ! Insert at the beginning of the list.
    !-----------------------------------------------------------------------
    IF ( ThisColl%nnDiagn > 0 ) THEN
       ThisDiagn%NextCont => ThisColl%DiagnList
    ENDIF
    ThisColl%DiagnList => ThisDiagn

    ! This diagnostics is now part of this collection
    ThisDiagn%CollectionID = PS

    ! Increase diagnostics counter and set container ID accordingly.
    ThisColl%nnDiagn = ThisColl%nnDiagn + 1

    ! Verbose mode
    IF ( HCO_IsVerb( 1 ) ) THEN
       WRITE(MSG,*) 'Successfully added diagnostics to collection ' , PS
       CALL HCO_MSG ( MSG )
       CALL Diagn_Print( ThisDiagn, 3 )
    ENDIF

    ! Cleanup
    ThisDiagn => NULL()
    ThisColl  => NULL()

    ! Return
    RC  = HCO_SUCCESS

  END SUBROUTINE Diagn_Create
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_UpdateSp
!
! !DESCRIPTION: Subroutine Diagn\_UpdateSp is the wrapper routine to update 
! the diagnostics for single precision arrays. It invokes the main diagnostics
! update routine with the appropriate arguments. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_UpdateSP( am_I_Root, cID,      cName,                  &
                           ExtNr,     Cat,        Hier,       HcoID,      &
                           AutoFill,  Scalar,     Array2D,    Array3D,    &
                           Total,     PosOnly,    COL,        RC           )
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root         ! Root CPU?
    INTEGER,          INTENT(IN   ), OPTIONAL :: cID               ! Assigned 
                                                                   !  container ID
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: cName             ! Diagnostics 
                                                                   !  name
    INTEGER,          INTENT(IN   ), OPTIONAL :: ExtNr             ! Extension #
    INTEGER,          INTENT(IN   ), OPTIONAL :: Cat               ! Category 
    INTEGER,          INTENT(IN   ), OPTIONAL :: Hier              ! Hierarchy 
    INTEGER,          INTENT(IN   ), OPTIONAL :: HcoID             ! HEMCO species
                                                                   !  ID number 
    INTEGER,          INTENT(IN   ), OPTIONAL :: AutoFill          ! 1=yes; 0=no; 
                                                                   ! -1=either 
    REAL(sp),         INTENT(IN   ), OPTIONAL :: Scalar            ! 1D scalar 
    REAL(sp),         INTENT(IN   ), OPTIONAL :: Array2D   (:,:)   ! 2D array 
    REAL(sp),         INTENT(IN   ), OPTIONAL :: Array3D   (:,:,:) ! 3D array 
    REAL(sp),         INTENT(IN   ), OPTIONAL :: Total             ! Total 
    LOGICAL,          INTENT(IN   ), OPTIONAL :: PosOnly           ! Use only vals
                                                                   !  >= 0?
    INTEGER,          INTENT(IN   ), OPTIONAL :: COL               ! Collection Nr.
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC                ! Return code 
!
! !REVISION HISTORY:
!  20 Apr 2015 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Call down to driver routine
    CALL Diagn_UpdateDriver( am_I_Root, & 
                       cID = cID, & 
                       cName = cName, & 
                       ExtNr = ExtNr, & 
                       Cat = Cat, & 
                       Hier = Hier, & 
                       HcoID = HcoID, & 
                       AutoFill = AutoFill, & 
                       Scalar_SP = Scalar, & 
                       Array2D_SP = Array2D, & 
                       Array3D_SP = Array3D, & 
                       Total_SP   = Total,   & 
                       PosOnly = PosOnly, & 
                       COL = COL, & 
                       RC = RC )

  END SUBROUTINE Diagn_UpdateSp
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_UpdateDp
!
! !DESCRIPTION: Subroutine Diagn\_UpdateDp is the wrapper routine to update 
! the diagnostics for double precision arrays. It invokes the main diagnostics
! update routine with the appropriate arguments. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_UpdateDP( am_I_Root, cID,      cName,                  &
                           ExtNr,     Cat,        Hier,       HcoID,      &
                           AutoFill,  Scalar,     Array2D,    Array3D,    &
                           Total,     PosOnly,    COL,        RC           )
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root         ! Root CPU?
    INTEGER,          INTENT(IN   ), OPTIONAL :: cID               ! Assigned 
                                                                   !  container ID
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: cName             ! Diagnostics 
                                                                   !  name
    INTEGER,          INTENT(IN   ), OPTIONAL :: ExtNr             ! Extension #
    INTEGER,          INTENT(IN   ), OPTIONAL :: Cat               ! Category 
    INTEGER,          INTENT(IN   ), OPTIONAL :: Hier              ! Hierarchy 
    INTEGER,          INTENT(IN   ), OPTIONAL :: HcoID             ! HEMCO species
                                                                   !  ID number 
    INTEGER,          INTENT(IN   ), OPTIONAL :: AutoFill          ! 1=yes; 0=no; 
                                                                   ! -1=either 
    REAL(dp),         INTENT(IN   ), OPTIONAL :: Scalar            ! 1D scalar 
    REAL(dp),         INTENT(IN   ), OPTIONAL :: Array2D   (:,:)   ! 2D array 
    REAL(dp),         INTENT(IN   ), OPTIONAL :: Array3D   (:,:,:) ! 3D array 
    REAL(dp),         INTENT(IN   ), OPTIONAL :: Total             ! Total 
    LOGICAL,          INTENT(IN   ), OPTIONAL :: PosOnly           ! Use only vals
                                                                   !  >= 0?
    INTEGER,          INTENT(IN   ), OPTIONAL :: COL               ! Collection Nr.
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC                ! Return code 
!
! !REVISION HISTORY:
!  20 Apr 2015 - C. Keller - Initialization
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Call down to driver routine
    CALL Diagn_UpdateDriver( am_I_Root, & 
                       cID = cID, & 
                       cName = cName, & 
                       ExtNr = ExtNr, & 
                       Cat = Cat, & 
                       Hier = Hier, & 
                       HcoID = HcoID, & 
                       AutoFill = AutoFill, & 
                       Scalar = Scalar, & 
                       Array2D = Array2D, & 
                       Array3D = Array3D, & 
                       Total   = Total,   & 
                       PosOnly = PosOnly, & 
                       COL = COL, & 
                       RC = RC )

  END SUBROUTINE Diagn_UpdateDp
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_UpdateDriver
!
! !DESCRIPTION: Subroutine Diagn\_UpdateDriver updates the content of a 
! diagnostics container. The container to be updated is determined
! from the passed variables. If a valid (i.e. positive) container 
! ID is provided, this container is used. Otherwise, if a valid 
! HEMCO species ID (HcoID) is provided, the container with the same 
! combination of HcoID, extension number (ExtNr), emission category 
! (Cat) and hierarchy (Hier) is used. If no valid HcoID and no valid
! cID is given, the container name has to be provided. The passed 
! data array (Scalar, Array2D, or Array3D) needs to match the 
! spatial dimension of the given container. For 2D diagnostics, a 3D
! array can be passed, in which case the level index specified
! during initialization (`LevIdx`) is used. If LevIdx is set to -1,
! the column sum is used (default).
!\\
!\\
! If no matching container is found, the subroutine leaves with no 
! error. This allows automatic diagnostics generation, e.g. of 
! intermediate emission fields created in HCO\_CALC\_Mod.F90.
!\\
!\\
! Notes:
! - For a given time step, the same diagnostics container can be 
!   updated multiple times. The field average is always defined as
!   temporal average, e.g. multiple updates on the same time step
!   will not increase the averaging weight of that time step.
! - If the passed array is empty (i.e. not associated), it is 
!   treated as empty values (i.e. zeros).
! - The collection number can be set to -1 to scan trough all 
!   existing diagnostic collections.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_UpdateDriver( am_I_Root, cID,        cName,                  &
                                 ExtNr,     Cat,        Hier,       HcoID,      &
                                 AutoFill,  Scalar,     Array2D,    Array3D,    &
                                 Total,     Scalar_SP,  Array2D_SP, Array3D_SP, &
                                 Total_SP,  PosOnly,    COL,        RC           )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
    USE HCO_Arr_Mod,   ONLY : HCO_ArrAssert
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )                   :: am_I_Root         ! Root CPU?
    INTEGER,          INTENT(IN   ), OPTIONAL         :: cID               ! container ID 
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL         :: cName             ! Dgn name
    INTEGER,          INTENT(IN   ), OPTIONAL         :: ExtNr             ! Extension #
    INTEGER,          INTENT(IN   ), OPTIONAL         :: Cat               ! Category 
    INTEGER,          INTENT(IN   ), OPTIONAL         :: Hier              ! Hierarchy 
    INTEGER,          INTENT(IN   ), OPTIONAL         :: HcoID             ! HEMCO species ID
    INTEGER,          INTENT(IN   ), OPTIONAL         :: AutoFill          ! 1=yes; 0=no; 
                                                                           ! -1=either 
    REAL(dp),         INTENT(IN   ), OPTIONAL         :: Scalar            ! 1D scalar 
    REAL(dp),         INTENT(IN   ), OPTIONAL, TARGET :: Array2D   (:,:)   ! 2D array 
    REAL(dp),         INTENT(IN   ), OPTIONAL, TARGET :: Array3D   (:,:,:) ! 3D array 
    REAL(dp),         INTENT(IN   ), OPTIONAL         :: Total             ! Total 
    REAL(sp),         INTENT(IN   ), OPTIONAL         :: Scalar_SP         ! 1D scalar 
    REAL(sp),         INTENT(IN   ), OPTIONAL, TARGET :: Array2D_SP(:,:)   ! 2D array 
    REAL(sp),         INTENT(IN   ), OPTIONAL, TARGET :: Array3D_SP(:,:,:) ! 3D array 
    REAL(sp),         INTENT(IN   ), OPTIONAL         :: Total_SP          ! Total 
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: PosOnly           ! Use only vals
                                                                           ! >= 0?
    INTEGER,          INTENT(IN   ), OPTIONAL         :: COL               ! Collection Nr.
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)                   :: RC                ! Return code 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller - Initialization
!  25 Sep 2014 - C. Keller - Now allow updating multiple diagnostics
!  11 Mar 2015 - C. Keller - Now allow scanning of all diagnostic collections
!  13 Mar 2015 - C. Keller - Bug fix: only prompt warning if it's a new timestep
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(DiagnCollection), POINTER :: ThisColl      => NULL()
    TYPE(DiagnCont),       POINTER :: ThisDiagn     => NULL()
    REAL(sp),              POINTER :: Arr2D (:,:)   => NULL()
    REAL(sp),              POINTER :: Tmp2D (:,:)   => NULL()
    REAL(sp),              POINTER :: Arr3D (:,:,:) => NULL()
    REAL(sp)                       :: TmpScalar

    ! Scalars
    CHARACTER(LEN=255)             :: LOC, MSG
    REAL(hp)                       :: Fact
    REAL(hp)                       :: Tmp
    CHARACTER(LEN=31)              :: DgnName
    INTEGER                        :: I, J, L, PS, AS
    INTEGER                        :: DgncID,  DgnExtNr, DgnCat
    INTEGER                        :: DgnHier, DgnHcoID
    INTEGER                        :: MinResetFlag, ThisUpdateID
    INTEGER                        :: AutoFlag
    INTEGER                        :: CNT
    LOGICAL                        :: Found, OnlyPos, VertSum, IsAssoc, IsNewTS
    LOGICAL                        :: InUse, SearchAll

    !======================================================================
    ! Diagn_UpdateDriver begins here!
    !======================================================================

    ! Init
    LOC = 'Diagn_UpdateDriver (hco_diagn_mod.F90)'
    RC  = HCO_SUCCESS

    ! Get collection number. 
    CALL DiagnCollection_DefineID( PS, RC, COL=COL, DEF=-1, &
           OKIfAll=.TRUE., InUse=InUse, ThisColl=ThisColl )
    IF ( RC /= HCO_SUCCESS ) RETURN
 
    ! Check if we need to scan through all collections. This is only the
    ! case if PS is set to -1
    IF ( PS == -1 ) THEN
       SearchAll = .TRUE.
    ELSE
       SearchAll = .FALSE.
    ENDIF

    ! Nothing to do if this collection is empty
    IF ( .NOT. SearchAll .AND. .NOT. InUse ) THEN
       RETURN 
    ENDIF

    !----------------------------------------------------------------------
    ! Make sure all attributes are defined
    !----------------------------------------------------------------------
    DgnName  = ''
    DgncID   = -1
    DgnExtNr = -1
    DgnCat   = -1
    DgnHier  = -1
    DgnHcoID = -1
    OnlyPos  = .FALSE.
    AutoFlag = -1
    IF ( PRESENT(cName   ) ) DgnName  = cName
    IF ( PRESENT(cID     ) ) DgncID   = cID  
    IF ( PRESENT(ExtNr   ) ) DgnExtNr = ExtNr
    IF ( PRESENT(Cat     ) ) DgnCat   = Cat  
    IF ( PRESENT(Hier    ) ) DgnHier  = Hier 
    IF ( PRESENT(HcoID   ) ) DgnHcoID = HcoID
    IF ( PRESENT(PosOnly ) ) OnlyPos  = PosOnly
    IF ( PRESENT(AutoFill) ) AutoFlag = AutoFill

    ! Get current minimum reset flag as well as the update time ID.
    MinResetFlag = HcoClock_GetMinResetFlag()
    CALL HcoClock_Get( nSteps = ThisUpdateID, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Count # of containers that are updated
    CNT = 0 

    !-----------------------------------------------------------------
    ! Loop over collections
    !-----------------------------------------------------------------
    DO WHILE ( ASSOCIATED(ThisColl) ) 

       ! Reset Diagnostics
       ThisDiagn => NULL()

       !-----------------------------------------------------------------
       ! Do for every container in the diagnostics list that matches the 
       ! specified arguments (ID, ExtNr, etc.). This can be more than one
       ! container (ckeller, 09/25/2014).
       !-----------------------------------------------------------------
       DO
  
          ! Search for diagnostics that matches the given arguments.
          ! If ThisDiagn is empty (first call), the search will start
          ! at the first diagnostics container. Otherwise, the search
          ! will resume from this diagnostics container.
          CALL DiagnCont_Find( DgncID,    DgnExtNr, DgnCat,   DgnHier, &
                               DgnHcoID,  DgnName,  AutoFlag, Found,   &
                               ThisDiagn, RESUME=.TRUE., COL=ThisColl%CollectionID )
   
          ! Exit while loop if no diagnostics found
          IF ( .NOT. Found ) EXIT
   
          ! If container holds just a pointer to external data, don't do
          ! anything!
          IF ( ThisDiagn%DtaIsPtr ) THEN
             MSG = 'You try to update a container that holds a '  // &
                  'pointer to data - this should never happen! ' // &
                  TRIM(ThisDiagn%cName)
             CALL HCO_WARNING( MSG, RC, THISLOC=LOC )
             CYCLE 
          ENDIF
     
          ! Increase counter
          CNT = CNT + 1
    
          !----------------------------------------------------------------------
          ! Check if this is a new time step for this diagnostics. 
          !----------------------------------------------------------------------
          IsNewTS = .TRUE.
          IF ( ThisDiagn%LastUpdateID == ThisUpdateID ) IsNewTS = .FALSE. 
   
          !----------------------------------------------------------------------
          ! Sanity check: if data is beyond its averaging interval, it should
          ! be in output format. Otherwise, the content of this diagnostics has 
          ! never been passed to the output yet (via routine Diagn\_Get) and 
          ! will thus be lost!
          !----------------------------------------------------------------------
          !IF ( (ThisDiagn%ResetFlag >= MinResetFlag) &
          IF ( (ThisColl%ResetFlag >= MinResetFlag) &
               .AND. .NOT. ThisDiagn%IsOutFormat     &
               .AND.      (ThisDiagn%Counter > 0)    &
               .AND.       IsNewTS                    ) THEN
             MSG = 'Diagnostics is at end of its output interval '    // &
                   'but was not passed to output - data may be lost ' // &
                   'or diagnostics may be wrong: ' // TRIM(ThisDiagn%cName)
             CALL HCO_WARNING( MSG, RC, THISLOC=LOC )
          ENDIF
 
          !----------------------------------------------------------------------
          ! If data is in output format, set counter to zero. This will make
          ! sure that the new data is not added to the existing data.
          !----------------------------------------------------------------------
          IF ( ThisDiagn%IsOutFormat ) THEN
             ThisDiagn%Counter    = 0
          ENDIF
      
          !----------------------------------------------------------------------
          ! Determine scale factor to be applied to data. Diagnostics are
          ! stored in kg/m2, hence need to multiply HEMCO emissions, which
          ! are in kg/m2/s, by emission time step. Don't do anything for 
          ! data with non-standard units (e.g. unitless factors) or pointers.
          ! Note: conversion to final output units is done when writing out
          ! the diagnostics.
          !----------------------------------------------------------------------
          IF ( ThisDiagn%AvgFlag > 0 ) THEN
             Fact = 1.0_hp 
          ELSE
             Fact = ThisColl%TS
          ENDIF
          
          !----------------------------------------------------------------------
          ! Fill shadow arrays. Cast any input data to single precision, as this
          ! is the default diagnostics precision. 
          !----------------------------------------------------------------------
   
          ! Only need to do this on first container. Afterwards, arrays are
          ! already set!
          IF ( CNT == 1 ) THEN
   
             ! 3D array
             IF ( PRESENT(Array3D_SP) ) THEN
                Arr3D => Array3D_SP
             ELSEIF( PRESENT(Array3D) ) THEN
                ALLOCATE( Arr3D(ThisColl%NX,ThisColl%NY,ThisColl%NZ),STAT=AS)
                IF ( AS /= 0 ) THEN
                   CALL HCO_ERROR( 'Allocation error Arr3D', RC, THISLOC=LOC )
                   RETURN
                ENDIF
                Arr3D = Array3D
             ENDIF
     
             ! 2D array 
             IF ( PRESENT(Array2D_SP) ) THEN
                Arr2D => Array2D_SP
             ELSEIF( PRESENT(Array2D) ) THEN
                ALLOCATE( Arr2D(ThisColl%NX,ThisColl%NY),STAT=AS)
                IF ( AS /= 0 ) THEN
                   CALL HCO_ERROR( 'Allocation error Arr2D', RC, THISLOC=LOC )
                   RETURN
                ENDIF
                Arr2D = Array2D
             ENDIF
     
             ! Scalar
             IF ( PRESENT(Scalar_SP) ) THEN
                TmpScalar = Scalar_SP
             ELSEIF ( PRESENT(Scalar) ) THEN
                TmpScalar = Scalar
             ENDIF
   
          ENDIF ! Counter = 1
       
          !----------------------------------------------------------------------
          ! To add 3D array
          !----------------------------------------------------------------------
          IF ( ThisDiagn%SpaceDim == 3 ) THEN
      
             ! Make sure dimensions agree and diagnostics array is allocated
             IF ( PRESENT(Array3D_SP) .OR. PRESENT(Array3D) ) THEN
      
                ! By default, write into single precision array 
                CALL HCO_ArrAssert( ThisDiagn%Arr3D, ThisColl%NX,   &
                                    ThisColl%NY,     ThisColl%NZ, RC ) 
                IF ( RC /= HCO_SUCCESS ) RETURN 
                  
                ! Pass array to diagnostics: reset to zero if counter 
                ! is zero, add to it otherwise.
                ! Never reset containers with cumulative sums!
                IF ( ThisDiagn%Counter == 0 .AND. &
                     ThisDiagn%AvgFlag /= AvgFlagCumulSum ) ThisDiagn%Arr3D%Val = 0.0_sp
      
                ! Always reset containers with instantaneous values if it's a new
                ! time step.
                IF ( ThisDiagn%AvgFlag == AvgFlagInst .AND. IsNewTS ) ThisDiagn%Arr3D%Val = 0.0_sp
      
                ! Only if associated ...
                IF ( ASSOCIATED(Arr3D) ) THEN
                   IF ( OnlyPos ) THEN
                      WHERE ( Arr3D >= 0.0_sp )
                         ThisDiagn%Arr3D%Val = ThisDiagn%Arr3D%Val + ( Arr3D * Fact )
                      END WHERE
                   ELSE
                      ThisDiagn%Arr3D%Val = ThisDiagn%Arr3D%Val + ( Arr3D * Fact )
                   ENDIF
                ENDIF
             ENDIF    
 
          !----------------------------------------------------------------------
          ! To add 2D array
          !----------------------------------------------------------------------
          ELSEIF ( ThisDiagn%SpaceDim == 2 ) THEN
     
             IF ( PRESENT(Array3D_SP) .OR. PRESENT(Array3D) .OR. & 
                  PRESENT(Array2D_SP) .OR. PRESENT(Array2D)       ) THEN
 
                ! Make sure dimensions agree and diagnostics array is allocated
                CALL HCO_ArrAssert( ThisDiagn%Arr2D, ThisColl%NX, &
                                    ThisColl%NY,     RC            ) 
                IF ( RC /= HCO_SUCCESS ) RETURN 
                  
                ! Pass array to diagnostics: ignore existing data if counter 
                ! is zero, add to it otherwise.
                ! Never reset containers with cumulative sums!
                IF ( ThisDiagn%Counter == 0 .AND. &
                     ThisDiagn%AvgFlag /= AvgFlagCumulSum ) ThisDiagn%Arr2D%Val = 0.0_sp
         
                ! Always reset containers with instantaneous values if it's a new time step
                IF ( ThisDiagn%AvgFlag == AvgFlagInst .AND. IsNewTS ) ThisDiagn%Arr2D%Val = 0.0_sp
      
                ! Assume that we don't have to take the vertical sum
                VertSum = .FALSE.
         
                ! Assume data pointer is associated
                IsAssoc = .TRUE.
         
                ! Convert 3D array to 2D if necessary - only use first level!!
                IF ( PRESENT(Array2D) .OR. PRESENT(Array2D_SP) ) THEN
                   IF ( .NOT. ASSOCIATED(Arr2D) ) THEN
                      IsAssoc = .FALSE.
                   ELSE
                      Tmp2D => Arr2D
                   ENDIF
                ELSEIF ( PRESENT(Array3D) .OR. PRESENT(Array3D_SP) ) THEN
                   IF ( .NOT. ASSOCIATED(Arr3D) ) THEN
                      IsAssoc = .FALSE.
                   ELSE
                      IF ( ThisDiagn%LevIdx == -1 ) THEN
                         VertSum = .TRUE.
                      ELSE
                         Tmp2D => Arr3D(:,:,ThisDiagn%LevIdx)
                      ENDIF
                   ENDIF
                ELSE
                   MSG = 'No array passed for updating ' // TRIM(ThisDiagn%cName)
                   CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
                   RETURN
                ENDIF
         
                ! Do only if data pointer associated ...
                IF ( IsAssoc ) THEN
         
                   ! only positive values
                   IF ( OnlyPos ) THEN
            
                      ! need to do vertical summation
                      IF ( VertSum ) THEN
                         DO J=1,ThisColl%NY
                         DO I=1,ThisColl%NX
                            TMP = 0.0_hp
                            DO L=1,ThisColl%NZ
                               IF ( Arr3D(I,J,L) >= 0.0_sp ) &
                                  TMP = TMP + ( Arr3D(I,J,L) * Fact )
                            ENDDO
                            ThisDiagn%Arr2D%Val(I,J) = &
                               ThisDiagn%Arr2D%Val(I,J) + TMP
                         ENDDO
                         ENDDO
            
                      ! no vertical summation
                      ELSE
                         WHERE ( Tmp2D >= 0.0_sp )
                            ThisDiagn%Arr2D%Val = ThisDiagn%Arr2D%Val + ( Tmp2D * Fact )
                         END WHERE
                      ENDIF
            
                   ! all values
                   ELSE
             
                      ! need to do vertical summation
                      IF ( VertSum ) THEN
                         DO J=1,ThisColl%NY
                         DO I=1,ThisColl%NX
                            TMP = SUM(Arr3D(I,J,:)) * Fact
                            ThisDiagn%Arr2D%Val(I,J) = &
                               ThisDiagn%Arr2D%Val(I,J) + TMP
                         ENDDO
                         ENDDO
         
                      ! no vertical summation
                      ELSE
                         ThisDiagn%Arr2D%Val = ThisDiagn%Arr2D%Val + ( Tmp2D * Fact )
                      ENDIF
                   ENDIF
                ENDIF ! pointer is associated
             ENDIF ! Array present      

          !----------------------------------------------------------------------
          ! To add scalar (1D) 
          !----------------------------------------------------------------------
          ELSEIF ( ThisDiagn%SpaceDim == 1 ) THEN
      
             ! Make sure dimensions agree and diagnostics array is allocated
             IF ( PRESENT(Scalar_SP) .OR. PRESENT(Scalar) ) THEN
   
                ! Pass array to diagnostics: ignore existing data if counter 
                ! is zero, add to it otherwise.
                ! Never reset containers with cumulative sums!
                IF ( ThisDiagn%Counter == 0 .AND. &
                     ThisDiagn%AvgFlag /= AvgFlagCumulSum ) ThisDiagn%Scalar = 0.0_sp
      
                ! Always reset containers with instantaneous values if it's a new time step
                IF ( ThisDiagn%AvgFlag == AvgFlagInst .AND. IsNewTS ) ThisDiagn%Scalar = 0.0_sp
      
                ! Update scalar value
                IF ( OnlyPos ) THEN
                   IF ( TmpScalar >= 0.0_sp ) & 
                      ThisDiagn%Scalar = ThisDiagn%Scalar + ( TmpScalar * Fact )
                ELSE
                   ThisDiagn%Scalar = ThisDiagn%Scalar + ( TmpScalar * Fact )  
                ENDIF
             ENDIF      
          ENDIF

          !----------------------------------------------------------------------
          ! Eventually update total
          !----------------------------------------------------------------------
          IF ( PRESENT(Total) ) THEN
             ThisDiagn%Total = ThisDiagn%Total + Total
          ENDIF
          IF ( PRESENT(Total_SP) ) THEN
             ThisDiagn%Total = ThisDiagn%Total + Total_SP
          ENDIF
      
          !----------------------------------------------------------------------
          ! Update counter ==> Do only if last update time is not equal to 
          ! current one! This allows the same diagnostics to be updated
          ! multiple time on the same time step without increasing the
          ! time step counter.
          !----------------------------------------------------------------------
          IF ( IsNewTS ) THEN
             ThisDiagn%Counter      = ThisDiagn%Counter + 1
             ThisDiagn%LastUpdateID = ThisUpdateID
          ENDIF
      
          !----------------------------------------------------------------------
          ! Data is not in output format and hasn't been called yet by Diagn_Get.
          !----------------------------------------------------------------------
          ThisDiagn%IsOutFormat = .FALSE.
          ThisDiagn%nnGetCalls  = 0

          ! Verbose mode 
          IF ( HCO_IsVerb( 2 ) ) THEN
             WRITE(MSG,'(a,a,a,I3,a)') 'Successfully updated diagnostics: ', &
                TRIM(ThisDiagn%cName), ' (counter:', ThisDiagn%Counter, ')'
             CALL HCO_MSG ( MSG )
          ENDIF
       ENDDO ! loop over containers in collection
   
       ! Advance to next collection
       IF ( SearchAll ) THEN
          ThisColl => ThisColl%NextCollection
       ELSE
          ThisColl => NULL()
       ENDIF

    ENDDO ! loop over collections

    ! Cleanup
    IF (PRESENT(Array3D_SP) ) THEN
       Arr3D => NULL()
    ELSEIF (PRESENT(Array3D) ) THEN
       IF ( ASSOCIATED(Arr3D) ) DEALLOCATE(Arr3D)
    ENDIF
    IF (PRESENT(Array2D_SP) ) THEN
       Arr2D => NULL()
    ELSEIF (PRESENT(Array2D) ) THEN
       IF ( ASSOCIATED(Arr2D) ) DEALLOCATE(Arr2D)
    ENDIF

    ! Return
    Tmp2D     => NULL()
    ThisDiagn => NULL()
    ThisColl  => NULL()
    RC        =  HCO_SUCCESS

  END SUBROUTINE Diagn_UpdateDriver
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_Get
!
! !DESCRIPTION: Subroutine Diagn\_Get returns a diagnostics container from
! the diagnostics list, with the data converted to the output unit specified
! during initialization. Only diagnostics that contain data, i.e. with an
! update counter higher than zero, are returned. If EndOfIntvOnly is set to 
! TRUE, only containers at the end of their time averaging interval are 
! returned. The current HEMCO time will be used to determine which containers
! are at the end of their interval. The IsOutFormat flag of the container is 
! set to true, making sure that the currently saved data will be erased during 
! the next update (Diagn\_Update).
!\\
!\\
! If DgnCont is already associated, the search continues from the container 
! next to DgnCont. If DgnCont is empty (null), the search starts from the 
! first container of the diagnostics list ListDiagn. If the optional attribute 
! cName or cID is provided, this particular container is searched (through the 
! entire diagnostics list), but is only returned if it is at the end of it's 
! interval or if EndOfIntvOnly is disabled.
!\\
!\\
! The optional argument InclManual denotes whether or not containers with
! a manual update frequency shall be considered. This argument is only valid
! if EndOfIntvOnly is set to FALSE.
!\\
!\\ 
! The return flag FLAG is set to HCO\_SUCCESS if a container is found, and to 
! HCO\_FAIL otherwise.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Get( am_I_Root, EndOfIntvOnly, DgnCont, &
                        FLAG,      RC,            cName,   &
                        cID,       AutoFill,      COL,     &
                        SkipZeroCount                       )
!
! !INPUT PARAMETERS::
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root       ! Root CPU?
    LOGICAL,          INTENT(IN   )           :: EndOfIntvOnly   ! End of 
                                                                 ! interval 
                                                                 ! only? 
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: cName           ! container name
    INTEGER,          INTENT(IN   ), OPTIONAL :: cID             ! container ID
    INTEGER,          INTENT(IN   ), OPTIONAL :: AutoFill        ! 0=no; 1=yes; 
                                                                 ! -1=either
    INTEGER,          INTENT(IN   ), OPTIONAL :: COL             ! Collection Nr. 
    LOGICAL,          INTENT(IN   ), OPTIONAL :: SkipZeroCount   ! Skip if counter
                                                                 ! is zero 
!
! !OUTPUT PARAMETERS:
!

    TYPE(DiagnCont),  POINTER                 :: DgnCont         ! Return 
                                                                 ! container
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: FLAG            ! Return flag
    INTEGER,          INTENT(INOUT)           :: RC              ! Return code 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DiagnCollection), POINTER :: ThisColl => NULL()
    INTEGER                        :: MinResetFlag
    INTEGER                        :: PS, AF
    LOGICAL                        :: FOUND, CF
    LOGICAL                        :: SKIPZERO 

    !======================================================================
    ! Diagn_Get begins here!
    !======================================================================

    ! Init
    FLAG   = HCO_FAIL
    RC     = HCO_SUCCESS
    CF     = .FALSE.

    ! Get collection number
    CALL DiagnCollection_DefineID( PS, RC, COL=COL, ThisColl=ThisColl )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set AutoFill flag
    AF = -1
    IF ( PRESENT(AutoFill  ) ) AF = AutoFill

    ! Check if diagnostics with counter = 0 shall be skipped
    SKIPZERO = .FALSE.
    IF ( PRESENT(SkipZeroCount) ) SKIPZERO = SkipZeroCount

    ! Get minimum reset flag for current time. Set reset flag to -1 if
    ! EndOFIntvOnly flag is disabled. This will make sure that all 
    ! diagnostics are considered.
    IF ( .NOT. EndOfIntvOnly ) THEN
       MinResetFlag = ResetFlagManually
    ELSE
       MinResetFlag = HcoClock_GetMinResetFlag()
    ENDIF

    ! If current time stamp is not at the end of an interval - or if 
    ! there is no diagnostics container in the list with a reset flag 
    ! smaller or equal to MinResetFlag - there will be no matching 
    ! container whatsoever. Can leave right here.
    IF ( MinResetFlag > ThisColl%ResetFlag ) THEN
       DgnCont => NULL()
       RETURN
    ENDIF

    ! If container name is given, search for diagnostics with 
    ! the given name. 
    IF ( PRESENT( cName ) ) THEN
       CALL DiagnCont_Find( -1, -1, -1, -1, -1, cName, &
                            AF, FOUND, DgnCont, COL=PS )

       IF ( .NOT. FOUND ) THEN
          DgnCont => NULL()
       ELSE
          ! Don't consider container if counter is zero. 
          IF ( SKIPZERO .AND. DgnCont%Counter == 0 ) THEN
             DgnCont => NULL()
          ENDIF
       ENDIF
       CF = .TRUE.

    ENDIF
   
    ! If container id is given, search for diagnostics with 
    ! the given container ID.
    IF ( PRESENT( cID ) .AND. .NOT. CF ) THEN
       CALL DiagnCont_Find( cID, -1, -1, -1, -1, '', &
                            AF, FOUND, DgnCont, COL=PS )
       IF ( .NOT. FOUND ) THEN
          DgnCont => NULL()
       ELSE
          ! Don't consider container if counter is zero. 
          IF ( SKIPZERO .AND. DgnCont%Counter == 0 ) THEN
             DgnCont => NULL()
          ENDIF
       ENDIF
       CF = .TRUE.
    ENDIF

    ! If no container selected yet, point to next container in 
    ! list (or to head of list if DgnCont is not yet associated). 
    ! Number of updates since last output must be larger than zero!
    IF ( .NOT. CF ) THEN 
       IF ( .NOT. ASSOCIATED( DgnCont ) ) THEN
          DgnCont => ThisColl%DiagnList
       ELSE
          DgnCont => DgnCont%NextCont
       ENDIF
       DO WHILE ( ASSOCIATED ( DgnCont ) )
          ! Skip zero counters
          IF ( SKIPZERO .AND. DgnCont%Counter <= 0 ) THEN
             DgnCont => DgnCont%NextCont
             CYCLE
          ENDIF

          ! Exit if we reach this loop here
          EXIT 
       ENDDO
    ENDIF

    ! Before returning container, make sure its data is ready for output.
    IF ( ASSOCIATED (DgnCont ) ) THEN
       CALL DiagnCont_PrepareOutput ( DgnCont, RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
       FLAG = HCO_SUCCESS

       ! Increase number of times this container has been called by
       ! Diagn_Get
       DgnCont%nnGetCalls = DgnCont%nnGetCalls + 1
    ENDIF

    ! Cleanup
    ThisColl => NULL()

  END SUBROUTINE Diagn_Get
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_TotalGet
!
! !DESCRIPTION: Subroutine Diagn\_TotalGet returns the total of a given
! diagnostics container.
!\\
!\\ 
! !INTERFACE:
!
  SUBROUTINE Diagn_TotalGet( am_I_Root, cName, cID,   COL, &
                             FOUND,     Total, Reset, RC    ) 
!
! !INPUT PARAMETERS::
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root      ! Root CPU?
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: cName          ! container name
    INTEGER,          INTENT(IN   ), OPTIONAL :: cID            ! container ID
    INTEGER,          INTENT(IN   ), OPTIONAL :: COL            ! Collection Nr. 
    LOGICAL,          INTENT(IN   ), OPTIONAL :: Reset          ! Reset total? 
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(  OUT), OPTIONAL :: FOUND          ! Container found 
    REAL(sp),         INTENT(  OUT)           :: Total          ! Container total 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC             ! Return code 
!
! !REVISION HISTORY:
!  15 Mar 2015 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DiagnCont),  POINTER  :: DgnCont => NULL() 
    INTEGER                    :: PS
    LOGICAL                    :: FND

    !======================================================================
    ! Diagn_TotalGet begins here!
    !======================================================================

    ! Init
    RC    = HCO_FAIL
    Total = 0.0_sp
    FND   = .FALSE.
    IF ( PRESENT(FOUND) ) THEN
       FOUND = .FALSE.
    ENDIF

    ! Get collection number
    CALL DiagnCollection_DefineID( PS, RC, COL=COL )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If container name is given, search for diagnostics with 
    ! the given name. 
    IF ( PRESENT( cName ) ) THEN
       CALL DiagnCont_Find( -1, -1, -1, -1, -1, cName, &
                            -1, FND, DgnCont, COL=PS )
    ENDIF
   
    ! If container id is given, search for diagnostics with 
    ! the given container ID.
    IF ( PRESENT( cID ) .AND. .NOT. FND ) THEN
       CALL DiagnCont_Find( cID, -1, -1, -1, -1, '', &
                            -1, FND, DgnCont, COL=PS )
    ENDIF

    ! Pass total to output
    IF ( FND .AND. ASSOCIATED ( DgnCont ) ) THEN
       Total = DgnCont%Total

       ! Eventually reset
       IF ( PRESENT(Reset) ) THEN
          IF ( Reset ) THEN
             DgnCont%Total = 0.0_sp
          ENDIF
       ENDIF

       ! Eventually update FOUND argument
       IF ( PRESENT(FOUND) ) THEN
          FOUND = .TRUE.
       ENDIF
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Diagn_TotalGet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DiagnList_Cleanup
!
! !DESCRIPTION: Subroutine DiagnList\_Cleanup cleans up all the diagnostics
! containers of the given diagnostics list.
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnList_Cleanup ( DiagnList )
! !INPUT ARGUMENTS:
!
    TYPE(DiagnCont), POINTER  :: DiagnList   ! List to be removed 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(DiagnCont), POINTER  :: TmpCont => NULL()
    TYPE(DiagnCont), POINTER  :: NxtCont => NULL()

    !======================================================================
    ! Diagn_Cleanup begins here!
    !======================================================================

    ! Walk through entire list and remove all containers
    TmpCont => DiagnList
    DO WHILE ( ASSOCIATED( TmpCont ) ) 

       ! Detach from list
       NxtCont => TmpCont%NextCont 

       ! Clean up this container 
       CALL DiagnCont_Cleanup( TmpCont )
       DEALLOCATE ( TmpCont )

       ! Advance
       TmpCont => NxtCont
    ENDDO

    ! Nullify DiagnList pointer
    DiagnList => NULL()

  END SUBROUTINE DiagnList_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_AutoFillLevelDefined
!
! !DESCRIPTION: Function Diagn\_AutoFillLevelDefined returns .TRUE. if there
! is at least one AutoFill diagnostics container defined for the given level
! (1=Species level, 2=ExtNr level, 3=Category level, 4=Hierarchy level).
!\\
!\\
! !INTERFACE:
!
  FUNCTION Diagn_AutoFillLevelDefined( Level, COL ) RESULT ( IsDefined )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN)           :: Level     ! Level of interest
    INTEGER, INTENT(IN), OPTIONAL :: COL       ! Collection Nr.
!
! !RETURN VALUE:
! 
    LOGICAL                       :: IsDefined ! Return argument 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DiagnCollection), POINTER :: ThisColl => NULL()
    INTEGER                        :: I, RC, PS
    LOGICAL                        :: InUse

    !======================================================================
    ! Diagn_AutoFillLevelDefined begins here!
    !======================================================================

    ! Initialize
    IsDefined = .FALSE.

    ! Get collection number
    CALL DiagnCollection_DefineID( PS, RC, COL=COL, DEF=-1, &
            OKIfAll=.TRUE., InUse=InUse, ThisColl=ThisColl )
    IF ( RC /= HCO_SUCCESS ) RETURN 

    ! Nothing to do if collection is not in use
    IF ( .NOT. InUse ) RETURN

    ! Do for every collection to be searched
    DO WHILE ( ASSOCIATED(ThisColl) )

       ! Check if autofill level is defined for this collection.
       ! If so, can leave here 
       IsDefined = ThisColl%AF_LevelDefined( Level )
       IF ( IsDefined ) EXIT
      
       ! Eventually go to next collection
       IF ( PS == -1 ) THEN
          ThisColl => ThisColl%NextCollection
       ELSE
          ThisColl => NULL()
       ENDIF
    ENDDO

    ! Cleanup
    ThisColl => NULL()
 
  END FUNCTION Diagn_AutoFillLevelDefined
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DiagnCollection_Get
!
! !DESCRIPTION: Subroutine DiagnCollection\_Get returns variables assigned to
! a given diagnostics collection. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCollection_Get( COL,       InUse,   Prefix, WriteFreq, &
                                  ResetFlag, nnDiagn, RC                  )
!
! !INPUT ARGUMENTS:
!
    INTEGER,          INTENT(IN), OPTIONAL :: COL       ! Collection Nr.
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(OUT), OPTIONAL :: InUse 
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: Prefix
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: WriteFreq 
    INTEGER,          INTENT(OUT), OPTIONAL :: ResetFlag 
    INTEGER,          INTENT(OUT), OPTIONAL :: nnDiagn 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    TYPE(DiagnCollection), POINTER :: ThisColl => NULL()
    INTEGER                        :: PS
    LOGICAL                        :: FOUND

    !======================================================================
    ! DiagnCollection_Get begins here!
    !======================================================================

    ! Init
    IF ( PRESENT(Prefix   ) ) Prefix    = ''
    IF ( PRESENT(WriteFreq) ) WriteFreq = ''
    IF ( PRESENT(ResetFlag) ) ResetFlag = -999
    IF ( PRESENT(InUse    ) ) InUse     = .FALSE. 
    IF ( PRESENT(nnDiagn  ) ) nnDiagn   = 0 

    ! Get collection number
    CALL DiagnCollection_DefineID( PS, RC, COL=COL, InUse=FOUND, ThisColl=ThisColl )
    IF ( RC /= HCO_SUCCESS ) RETURN

    IF ( PRESENT(InUse) ) THEN
       InUse = FOUND
    ENDIF   

    ! Get variables from collection 
    IF ( FOUND ) THEN
       IF ( PRESENT(Prefix) ) THEN
          Prefix = ThisColl%PREFIX
       ENDIF
   
       IF ( PRESENT(WriteFreq) ) THEN
          WriteFreq = ThisColl%WriteFreq
       ENDIF
   
       IF ( PRESENT(ResetFlag) ) THEN
          ResetFlag = ThisColl%ResetFlag
       ENDIF

       IF ( PRESENT(nnDiagn) ) THEN
          nnDiagn = ThisColl%nnDiagn
       ENDIF
    ENDIF
   
    ! Cleanup
    ThisColl => NULL()
    
    ! Return w/ success
    RC = HCO_SUCCESS
 
  END SUBROUTINE DiagnCollection_Get
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DiagnCont_Init
!
! !DESCRIPTION: Subroutine DiagnCont\_Init initializes a new (blank) 
! diagnostics container DgnCont. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCont_Init( OutCont )
!
! !OUTPUT PARAMETERS:
!
    TYPE(DiagnCont),       POINTER      :: OutCont    ! Created container
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(DiagnCont), POINTER :: DgnCont => NULL()

    !======================================================================
    ! DiagnCont_Init begins here!
    !======================================================================

    ! Allocate the new container
    ALLOCATE( DgnCont )

    ! Initialize ponters and scalar value 
    DgnCont%NextCont => NULL()
    DgnCont%Arr2D    => NULL()
    DgnCont%Arr3D    => NULL()
    DgnCont%DtaIsPtr = .FALSE.
    DgnCont%Scalar   =  0.0_sp
    DgnCont%Total    =  0.0_sp
    DgnCont%LevIdx   = -1
    DgnCont%AutoFill =  0
    DgnCont%SpaceDim =  2

    ! Default values for unit conversion factors 
    DgnCont%MassScal  = 1.0_hp
    DgnCont%AreaScal  = 1.0_hp
    DgnCont%ScaleFact = 1.0_hp
    DgnCont%AreaFlag  = 2
    DgnCont%Counter   = 0
    DgnCont%TimeAvg   = -1   
    DgnCont%AvgFlag   = -1
    DgnCont%AvgName   = 'mean'

    ! Set last update time to -1 to start with
    DgnCont%LastUpdateID = -1

    ! By default, data is not in output format
    DgnCont%IsOutFormat = .FALSE.
    DgnCont%nnGetCalls  = 0

    ! Default container ID and collection 
    DgnCont%cID          = -1
    DgnCont%CollectionID = -1

    ! Initialize other varaibles
    DgnCont%HcoID        = -1
    DgnCont%ExtNr        = -1
    DgnCont%Cat          = -1
    DgnCont%Hier         = -1

    ! Pass to output container
    OutCont => DgnCont

  END SUBROUTINE DiagnCont_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DiagnCont_Cleanup
!
! !DESCRIPTION: Subroutine DiagnCont\_Cleanup cleans up diagnostics 
! container DgnCont.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCont_Cleanup( DgnCont )
!
! !USES:
!
    USE HCO_ARR_Mod, ONLY : HCO_ArrCleanup 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DiagnCont), POINTER :: DgnCont  ! Container to be cleaned
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC

    LOGICAL :: DeepClean

    !======================================================================
    ! DiagnCont_Cleanup begins here!
    !======================================================================

    ! Only if associated... 
    IF ( ASSOCIATED( DgnCont ) ) THEN
       IF ( DgnCont%DtaIsPtr ) THEN
          DeepClean = .FALSE.
       ELSE
          DeepClean = .TRUE.
       ENDIF
       CALL HCO_ArrCleanup( DgnCont%Arr2D, DeepClean )
       CALL HCO_ArrCleanup( DgnCont%Arr3D, DeepClean )
       DgnCont%NextCont => NULL()
       DEALLOCATE ( DgnCont )
    ENDIF

  END SUBROUTINE DiagnCont_Cleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DiagnCont_PrepareOutput 
!
! !DESCRIPTION: Subroutine DiagnCont\_PrepareOutput converts the data of
! the given diagnostics container to proper output units.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCont_PrepareOutput ( DgnCont, RC )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
    USE HCO_Arr_Mod,   ONLY : HCO_ArrAssert
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DiagnCont),   POINTER       :: DgnCont  ! diagnostics container 
    INTEGER,           INTENT(INOUT) :: RC       ! Return code 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DiagnCollection), POINTER :: ThisColl => NULL()
    LOGICAL                        :: FOUND
    INTEGER                        :: I, J, YYYY, MM
    REAL(hp)                       :: norm1, mult1, DPY, totscal
    CHARACTER(LEN=255)             :: MSG, LOC
    INTEGER                        :: DPM(12) = (/ 31, 28, 31, 30, 31, 30, &
                                                   31, 31, 30, 31, 30, 31   /)

    !======================================================================
    ! DiagnCont_PrepareOutput begins here!
    !======================================================================

    ! Init
    RC  = HCO_SUCCESS
    LOC = 'DiagnCont_PrepareOutput (hco_diagn_mod.F90) '
    
    !-----------------------------------------------------------------------
    ! Don't do anything for pointer data and/or if data is already in 
    ! output format
    !-----------------------------------------------------------------------
    IF ( DgnCont%IsOutFormat ) RETURN
    IF ( DgnCont%DtaIsPtr    ) RETURN

    !-----------------------------------------------------------------------
    ! Get pointer to this collection
    !-----------------------------------------------------------------------
    CALL DiagnCollection_Find( DgnCont%CollectionID, FOUND, RC, ThisColl=ThisColl )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! This should never happen
    IF ( .NOT. FOUND .OR. .NOT. ASSOCIATED(ThisColl) ) THEN
       WRITE(MSG,*) 'Diagnostics ', TRIM(DgnCont%cName), ' has invalid ', &
                    'collection ID of ', DgnCont%CollectionID
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Return zero array if counter is still zero 
    !-----------------------------------------------------------------------
    IF ( DgnCont%Counter == 0 ) THEN
      
       ! Make sure array is defined and zero
       IF ( DgnCont%SpaceDim == 2 ) THEN
          CALL HCO_ArrAssert( DgnCont%Arr2D, ThisColl%NX, &
                              ThisColl%NY,   RC            ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 

          ! Make sure it's zero
          DgnCont%Arr2D%Val = 0.0_sp

       ELSEIF ( DgnCont%SpaceDim == 3 ) THEN
          CALL HCO_ArrAssert( DgnCont%Arr3D, ThisColl%NX, &
                              ThisColl%NY,   ThisColl%NZ, RC ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 

          ! Make sure it's zero
          DgnCont%Arr3D%Val = 0.0_sp
       ENDIF

       ! Prompt warning
       MSG = 'Diagnostics counter is zero - return empty array: ' // &
             TRIM(DgnCont%cName)
       CALL HCO_WARNING( MSG, RC, THISLOC=LOC, WARNLEV=1 )    
       RETURN
    ENDIF 

    !-----------------------------------------------------------------------
    ! Output data is calculated as: 
    ! out = saved / norm1 * mult1 * MassScal * AreaScal
    ! The normalization factor norm1 and multiplication factor mult1
    ! are used to average to the desired time interval, i.e. per
    ! second, per day, etc. 
    ! Since all diagnostics are internally stored in units of [kg/m2] 
    ! first convert to [kg/m2/s] and then multiply by the desired time
    ! averaging interval (e.g. seconds/hour to get kg/m2/s). Factors
    ! MassScal and AreaScal convert mass and area to desired units, as
    ! determined during initialization of the diagnostics.
    !-----------------------------------------------------------------------

    ! If the averaging is forced to the sum: 
    IF ( DgnCont%AvgFlag == AvgFlagSum .OR. DgnCont%AvgFlag == AvgFlagCumulSum ) THEN
       norm1 = 1.0_hp
       mult1 = 1.0_dp

    ! If the averaging is forced to the arithmetic mean: 
    ELSEIF ( DgnCont%AvgFlag == AvgFlagMean ) THEN
       norm1 = REAL(DgnCont%Counter,kind=hp)
       mult1 = 1.0_dp

    ! If there is no time averaging interval defined 
    ELSEIF ( DgnCont%TimeAvg < 0 ) THEN
       norm1 = 1.0_hp
       mult1 = 1.0_dp

    ! For other, time averaging intervals 
    ELSE

       ! Get current month and year
       CALL HcoClock_Get( cYYYY = YYYY, cMM = MM, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Days per year
       IF ( (MOD(YYYY,4) == 0) .AND. (MOD(YYYY,400) /= 0) ) THEN
          DPY    = 366.0_hp
          DPM(2) = 29
       ELSE
          DPY = 365.0_hp
       ENDIF

       ! Seconds since last reset
       norm1 = REAL(DgnCont%Counter,kind=hp) * ThisColl%TS

       ! Factors depends on averaging time
       IF ( DgnCont%TimeAvg == 1 ) THEN
          mult1 = 1.0_hp                     ! seconds / second

       ELSEIF ( DgnCont%TimeAvg == 2 ) THEN
          mult1 = 3600.0_hp                  ! seconds / hour

       ELSEIF ( DgnCont%TimeAvg == 3 ) THEN
          mult1 = 86400.0_hp                 ! seconds / day 

       ELSEIF ( DgnCont%TimeAvg == 4 ) THEN
          mult1 = 86400.0_hp * DPM(MM)       ! seconds / month 

       ELSEIF ( DgnCont%TimeAvg == 5 ) THEN
          mult1 = 86400.0_hp * DPY           ! seconds / year

       ! We shouldn't get here!
       ELSE
          WRITE(MSG,*) 'Illegal time averaging of ', DgnCont%TimeAvg, &
                       ' for diagnostics ', TRIM(DgnCont%cName)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

    ENDIF

    ! Error trap
    IF ( norm1 <= 0.0_hp ) THEN
       MSG = 'Illegal normalization factor: ' // TRIM(DgnCont%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! totscal is the combined scale factor
    totscal = mult1             &
            / norm1             &
            * DgnCont%MassScal  &
            * DgnCont%AreaScal  &
            * DgnCont%ScaleFact

    ! For 3D:
    IF ( DgnCont%SpaceDim == 3 ) THEN
       DO J = 1, ThisColl%NY
       DO I = 1, ThisColl%NX

          ! Multiply by area if output unit is not per area 
          IF ( DgnCont%AreaFlag == 0 ) THEN
             IF( ASSOCIATED(DgnCont%Arr3D) ) THEN
                DgnCont%Arr3D%Val(I,J,:) = DgnCont%Arr3D%Val(I,J,:)  & 
                                         * ThisColl%AREA_M2(I,J)
             ENDIF 
          ENDIF

          ! Apply scale factors
          IF ( ASSOCIATED(DgnCont%Arr3D) ) THEN
             DgnCont%Arr3D%Val(I,J,:) = DgnCont%Arr3D%Val(I,J,:) & 
                                      * totscal
          ENDIF
       ENDDO !I
       ENDDO !J

    ! For 2D:
    ELSEIF ( DgnCont%SpaceDim == 2 ) THEN
       DO J = 1, ThisColl%NY
       DO I = 1, ThisColl%NX

          ! Multiply by area if output unit is not per area 
          IF ( DgnCont%AreaFlag == 0 ) THEN
             IF ( ASSOCIATED(DgnCont%Arr2D) ) THEN
                DgnCont%Arr2D%Val(I,J) = DgnCont%Arr2D%Val(I,J) &
                                       * ThisColl%AREA_M2(I,J) 
             ENDIF
          ENDIF

          ! Apply scale factors
          IF ( ASSOCIATED(DgnCont%Arr2D) ) THEN
             DgnCont%Arr2D%Val(I,J) = DgnCont%Arr2D%Val(I,J) &
                                    * totscal
          ENDIF

       ENDDO !I
       ENDDO !J

    ! For 1D:
    ELSE
       DgnCont%Scalar = DgnCont%Scalar * totscal

    ENDIF

    ! Data is now in output format
    DgnCont%IsOutFormat = .TRUE.

    ! Cleanup
    ThisColl => NULL()

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE DiagnCont_PrepareOutput
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DiagnCont_Find 
!
! !DESCRIPTION: Subroutine DiagnCont\_Find searches for a diagnostics
! container in ListDiagn. If a valid container ID (>0) is given, the 
! container with this ID is searched. Otherwise, if a valid HEMCO
! specied ID (>0) is given, the container with the same combination
! of HcoID, extension number (ExtNr), and emission category (Cat) and
! hierarchy (Hier), is searched. If no valid HcoID and no valid cID is
! provided, the container with the given container name is searched.
!\\
!\\
! If the optional resume flag is set to TRUE, search will resume after 
! OutCnt. If OutCnt is not associated or resume flag is FALSE, search
! starts at the beginning of the diagnostics list.
!\\
!\\
! This subroutine does return the diagnostics as is, i.e. in the internal
! units. It should NOT be used to access the content of a diagnostics but is
! rather intended to be used in the background, e.g. to check if a 
! diagnostics exists at all. To get the values of a diagnostics, use routine
! Diagn\_Get.
!
! !INTERFACE:
!
  SUBROUTINE DiagnCont_Find ( cID,   ExtNr,    Cat,   Hier,   HcoID, &
                              cName, AutoFill, FOUND, OutCnt, Resume, COL )
!
! !INPUT PARAMETERS:
!
    INTEGER,           INTENT(IN)   :: cID      ! wanted cont. ID
    INTEGER,           INTENT(IN)   :: ExtNr    ! wanted ExtNr
    INTEGER,           INTENT(IN)   :: Cat      ! wanted category
    INTEGER,           INTENT(IN)   :: Hier     ! wanted hierarchy
    INTEGER,           INTENT(IN)   :: HcoID    ! wanted spec. ID
    CHARACTER(LEN=*),  INTENT(IN)   :: cName    ! wanted name 
    INTEGER,           INTENT(IN)   :: AutoFill ! 0=no; 1=yes; -1=either 
    LOGICAL, OPTIONAL, INTENT(IN)   :: Resume   ! Resume at OutCnt?
    INTEGER, OPTIONAL, INTENT(IN)   :: COL      ! Collection number 
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,           INTENT(OUT)  :: FOUND    ! container found?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DiagnCont),   POINTER      :: OutCnt   ! data container 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!  25 Sep 2014 - C. Keller: Added Resume flag
!  09 Apr 2015 - C. Keller: Can now search all collections
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                            :: RC, PS
    TYPE(DiagnCont),       POINTER     :: CurrCnt  => NULL() 
    TYPE(DiagnCollection), POINTER     :: ThisColl => NULL() 
    LOGICAL                            :: IsMatch, InUse, Rsm
 
    !======================================================================
    ! DiagnCont_Find begins here!
    !======================================================================

    ! Initialize
    FOUND  = .FALSE.

    ! Get collection number
    CALL DiagnCollection_DefineID( PS, RC, COL=COL, Def=-1, InUse=InUse, &
                                   OkIfAll=.TRUE., ThisColl=ThisColl )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave if collection not in use
    IF ( .NOT. InUse ) RETURN

    ! Resume from OutCnt
    IF ( PRESENT(Resume) ) THEN
       RSM = Resume
    ELSE
       RSM = .FALSE.
    ENDIF

    ! Make CurrCnt point to first element of the diagnostics list or to
    ! the container after OutCnt if resume flag is activated.
    IF ( RSM .AND. ASSOCIATED(OutCnt) ) THEN
       CurrCnt => OutCnt%NextCont
    ELSE 
       CurrCnt => ThisColl%DiagnList 
    ENDIF

    ! Error trap
    IF ( .NOT. ASSOCIATED(CurrCnt) ) THEN
       OutCnt => NULL()
       RETURN
    ENDIF

    ! Loop over all collections
    DO

       ! Now reset OutCnt. Will be defined again when diagnostics is found.
       OutCnt => NULL()

       ! Loop over list until container found
       DO WHILE ( ASSOCIATED ( CurrCnt ) )
   
          ! Check if this is the container of interest. If a valid 
          ! container ID is given, use this attribute. Otherwise, check
          ! for correct match of ExtNr, HcoID, Cat, and Hier attributes
          ! if a valid HcoID is given. Otherwise, use the container name.
          IsMatch = .FALSE.
   
          ! Check AutoFill flag.
          IF ( CurrCnt%AutoFill /= AutoFill .AND. AutoFill >= 0 ) THEN
             CurrCnt => CurrCnt%NextCont
             CYCLE
          ENDIF
   
          ! For valid container ID:
          IF ( cID > 0 ) THEN
             IF ( CurrCnt%cID == cID ) IsMatch = .TRUE.
        
          ! For valid HcoID, check for correct match of HcoID, ExtNr, 
          ! category, and hierarchy.
          ELSEIF ( HcoID > 0 ) THEN
             IF ( CurrCnt%HcoID == HcoID .AND. &
                  CurrCnt%ExtNr == ExtNr .AND. &
                  CurrCnt%Hier  == Hier  .AND. &
                  CurrCnt%Cat   == Cat          ) IsMatch = .TRUE.
   
          ! Use container name otherwise:
          ELSE
             IF ( TRIM(CurrCnt%cName) == TRIM(cName) ) IsMatch = .TRUE.
          ENDIF
   
          IF ( IsMatch ) THEN
             OutCnt  => CurrCnt
             FOUND   = .TRUE.
             EXIT 
          ENDIF
   
          ! Advance to next field otherwise
          CurrCnt => CurrCnt%NextCont
       ENDDO
  
       ! Leave loop over all collections if container was found
       IF ( FOUND ) EXIT

       ! Advance to next collection
       IF ( PS == -1 ) THEN
       
          ! Point to next collection
          ThisColl => ThisColl%NextCollection

          ! Leave if collection is empty
          IF ( .NOT. ASSOCIATED(ThisColl) ) EXIT

          ! Make working pointer point to first container in this collection,
          ! then resume list search
          CurrCnt => ThisColl%DiagnList
          CYCLE
       ENDIF

       ! Leave collection loop if we get here 
       EXIT
    ENDDO ! Loop over all collections
 
    ! Cleanup
    CurrCnt  => NULL()
    ThisColl => NULL()

  END SUBROUTINE DiagnCont_Find
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DiagnCont_Link_2D 
!
! !DESCRIPTION: Subroutine DiagnCont\_Link\_2D links the data of container
! DgnCont to the 2D array Tgt2D. This will disable all time averaging, 
! unit conversion, etc., i.e. the data will be returned as is.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCont_Link_2D( am_I_Root, DgnCont, ThisColl, Trgt2D, RC )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !INPUT ARGUMENTS:
!
    LOGICAL,               INTENT(IN   )          :: am_I_Root   ! Root CPU?
    REAL(sp),              INTENT(IN   ), TARGET  :: Trgt2D(:,:) ! 2D target data 
    TYPE(DiagnCollection), POINTER                :: ThisColl    ! Collection
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DiagnCont),       POINTER                :: DgnCont     ! diagnostics container 
    INTEGER,               INTENT(INOUT)          :: RC          ! Return code 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)             :: MSG, LOC 

    !======================================================================
    ! DiagnCont_Link_2D begins here!
    !======================================================================

    ! Init
    LOC = 'DiagnCont_Link_2D (hco_diagn_mod.F90)'

    ! Check if dimensions match. Also, containers with pointers must not
    ! be set to AutoFill
    IF ( DgnCont%AutoFill == 1 ) THEN
       MSG = 'Target diagnostics has AutoFill flag of 1 - reset to 0: ' & 
           // TRIM(DgnCont%cName)
       CALL HCO_WARNING( MSG, RC, THISLOC=LOC, WARNLEV=2 )
       DgnCont%AutoFill = 0
    ENDIF
    IF ( DgnCont%SpaceDim /= 2 ) THEN
       MSG = 'Diagnostics is not 2D: ' // TRIM(DgnCont%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    IF ( SIZE(Trgt2D,1) /= ThisColl%NX .OR. &
         SIZE(Trgt2D,2) /= ThisColl%NY       ) THEN
       MSG = 'Incorrect target array size: ' // TRIM(DgnCont%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Define 2D array pointer
    CALL HCO_ArrInit( DgnCont%Arr2D, 0, 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Point to data
    DgnCont%Arr2D%Val => Trgt2D

    ! Update pointer switch. This will make sure that data is not modified. 
    ! Also set counter to non-zero to make sure that diagnostics will be 
    ! correctly written. 
    DgnCont%DtaIsPtr = .TRUE. 
    DgnCont%Counter  = 1

    ! Return
    RC = HCO_SUCCESS

  END SUBROUTINE DiagnCont_Link_2D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DiagnCont_Link_3D 
!
! !DESCRIPTION: Subroutine DiagnCont\_Link\_3D links the data of container
! DgnCont to the 3D array Tgt3D. This will disable all time averaging, 
! unit conversion, etc., i.e. the data will be returned as is.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCont_Link_3D( am_I_Root, DgnCont, ThisColl, Trgt3D, RC )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !INPUT PARAEMTERS:
!
    LOGICAL,               INTENT(IN   )         :: am_I_Root     ! Root CPU?
    REAL(sp),              INTENT(IN   ), TARGET :: Trgt3D(:,:,:) ! 3D target data 
    TYPE(DiagnCollection), POINTER               :: ThisColl      ! Collection
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DiagnCont),       POINTER               :: DgnCont       ! diagnostics 
                                                                  !  container
    INTEGER,               INTENT(INOUT)         :: RC            ! Return code 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
    CHARACTER(LEN=255)             :: MSG, LOC 

    !======================================================================
    ! DiagnCont_Link_3D begins here!
    !======================================================================

    ! Init
    LOC = 'DiagnCont_Link_3D (hco_diagn_mod.F90)'

    ! Check if dimensions match. Also, containers with pointers must not
    ! be set to AutoFill
    IF ( DgnCont%AutoFill == 1 ) THEN
       MSG = 'Cannot link AutoFill container: ' // TRIM(DgnCont%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( DgnCont%AutoFill == 1 ) THEN
       MSG = 'Target diagnostics has autofill flag of 1 - reset to 0: ' & 
           // TRIM(DgnCont%cName)
       CALL HCO_WARNING( MSG, RC, THISLOC=LOC, WARNLEV=2 )
       DgnCont%AutoFill = 0
    ENDIF
    IF ( DgnCont%SpaceDim /= 3 ) THEN
       MSG = 'Diagnostics is not 3D: ' // TRIM(DgnCont%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Check array size
    IF ( SIZE(Trgt3D,1) /= ThisColl%NX .OR. &
         SIZE(Trgt3D,2) /= ThisColl%NY .OR. &
         SIZE(Trgt3D,3) /= ThisColl%NZ       ) THEN
       MSG = 'Incorrect target array size: ' // TRIM(DgnCont%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
   
    ! Define 3D array pointer
    CALL HCO_ArrInit( DgnCont%Arr3D, 0, 0, 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
   
    ! Point to data
    DgnCont%Arr3D%Val => Trgt3D

    ! Update pointer switch. This will make sure that data is not modified.
    ! Also set counter to non-zero to make sure that diagnostics will be 
    ! correctly written. 
    DgnCont%DtaIsPtr    = .TRUE. 
    DgnCont%Counter     = 1
    DgnCont%IsOutFormat = .TRUE.

    ! Return
    RC = HCO_SUCCESS

  END SUBROUTINE DiagnCont_Link_3D
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_Print
!
! !DESCRIPTION: Subroutine Diagn\_Print displays the content of the
! passed diagnostics container.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Print ( Dgn, VerbNr )
!
! !INPUT ARGUMENTS:
!
    TYPE(DiagnCont),       POINTER    :: Dgn
    INTEGER,               INTENT(IN) :: VerbNr  
!
! !REVISION HISTORY:
!  01 Aug 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ARGUMENTS:
!
    TYPE(DiagnCollection), POINTER    :: ThisColl => NULL()
    CHARACTER(LEN=255)                :: MSG 
    CHARACTER(LEN= 31)                :: WriteFreq 
    INTEGER                           :: RC, PS, nx, ny, nz
    REAL(sp)                          :: sm

    ! ================================================================
    ! Diagn_Print begins here
    ! ================================================================

    ! Get collection number
    CALL DiagnCollection_DefineID( PS, RC, &
       COL=Dgn%CollectionID, ThisColl=ThisColl )
    IF ( RC /= HCO_SUCCESS ) RETURN

    sm = 0.0_sp
    nx = 0 
    ny = 0
    nz = 0
    IF ( Dgn%SpaceDim<=2 ) THEN
       IF ( ASSOCIATED(Dgn%Arr2D) ) THEN
          nx = SIZE(Dgn%Arr2D%Val,1)
          ny = SIZE(Dgn%Arr2D%Val,2)
          sm = SUM(Dgn%Arr2D%Val)
       ENDIF
    ELSE
       IF ( ASSOCIATED(Dgn%Arr3D) ) THEN
          nx = SIZE(Dgn%Arr3D%Val,1)
          ny = SIZE(Dgn%Arr3D%Val,2)
          nz = SIZE(Dgn%Arr3D%Val,3)
          sm = SUM(Dgn%Arr3D%Val)
       ENDIF
    ENDIF

    ! Always print name 
    MSG = 'Container ' // TRIM(Dgn%cName)
    CALL HCO_MSG(MSG)

    ! Eventually add details
    IF ( HCO_IsVerb( VerbNr ) ) THEN

       ! Write frequency
       WriteFreq = ThisColl%WriteFreq

       ! General information
       WRITE(MSG,*) '   --> Collection         : ', Dgn%CollectionID
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) '   --> Diagn ID           : ', Dgn%cID
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) '   --> Extension Nr       : ', Dgn%ExtNr
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) '   --> Category           : ', Dgn%Cat
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) '   --> Hierarchy          : ', Dgn%Hier
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) '   --> HEMCO species ID   : ', Dgn%HcoID
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) '   --> Autofill?            ', Dgn%AutoFill
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) '   --> Space dimension    : ', Dgn%SpaceDim
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) '   --> Used level index   : ', Dgn%LevIdx 
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) '   --> Output unit        : ', TRIM(Dgn%OutUnit)
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) '   --> Uniform scaling    : ', Dgn%ScaleFact
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) '   --> Write frequency    : ', TRIM(WriteFreq)
       CALL HCO_MSG(MSG)
       WRITE(MSG,*) '   --> Current array sum  : ', sm
       CALL HCO_MSG(MSG)
    ENDIF

    ! Cleanup
    ThisColl => NULL()

  END SUBROUTINE Diagn_Print
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DiagnCollection_Create
!
! !DESCRIPTION: Subroutine DiagnCollection\_Create creates a new diagnostics 
! collection at position COL. The class arguments are set as specified by the 
! input arguments.
!\\
!\\
! If the given position is already occupied, the routine returns an error if 
! the input argument do not match with the corresponding arguments of the 
! diagnostics class at that position. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCollection_Create ( am_I_Root, NX,     NY, NZ,    &
                                      TS,   AM2, PREFIX, WriteFreq, &
                                      RC,   COL                      )
!
! !USES:
!
!
! !INPUT ARGUMENTS:
!
    LOGICAL,            INTENT(IN)           :: am_I_Root  ! Root CPU?
    INTEGER,            INTENT(IN)           :: NX         ! # of lons
    INTEGER,            INTENT(IN)           :: NY         ! # of lats
    INTEGER,            INTENT(IN)           :: NZ         ! # of levels
    REAL(sp),           INTENT(IN)           :: TS         ! timestep [s] 
    REAL(hp),           POINTER              :: AM2(:,:)   ! grid box areas [m2]
    CHARACTER(LEN=*),   INTENT(IN)           :: PREFIX     ! Output prefix
    CHARACTER(LEN=*),   INTENT(IN)           :: WriteFreq  ! Output frequency 
!
! !OUTPUT ARGUMENTS:
!
    INTEGER,            INTENT(  OUT)        :: COL        ! Collection Nr. 
!
! !INPUT/OUTPUT ARGUMENTS:
!
    INTEGER,            INTENT(INOUT)        :: RC         ! Return code
!
! !REVISION HISTORY:
!  08 Jan 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ARGUMENTS:
!
    TYPE(DiagnCollection), POINTER :: NewCollection => NULL()
    INTEGER                        :: PS
    CHARACTER(LEN=255)             :: MSG
    CHARACTER(LEN=255)             :: LOC = 'DiagnCollection_Create (hco_diagn_mod.F90)' 

    ! ================================================================
    ! DiagnCollection_Create begins here
    ! ================================================================

    ! Allocate new collection
    ALLOCATE(NewCollection)

    ! Pass arguments
    NewCollection%NX      =  NX 
    NewCollection%NY      =  NY 
    NewCollection%NZ      =  NZ 
    NewCollection%TS      =  TS 
    NewCollection%AREA_M2 => AM2 

    ! Set prefix
    NewCollection%PREFIX = TRIM(PREFIX)

    ! Determine output frequency
    IF ( TRIM(WriteFreq) == 'Annually' ) THEN
       NewCollection%ResetFlag = ResetFlagAnnually 

    ! Write out every month
    ELSEIF ( TRIM(WriteFreq) == 'Monthly' ) THEN
       NewCollection%ResetFlag = ResetFlagMonthly
         
    ! Write out every day
    ELSEIF ( TRIM(WriteFreq) == 'Daily' ) THEN
       NewCollection%ResetFlag = ResetFlagDaily

    ! Write out every hour
    ELSEIF ( TRIM(WriteFreq) == 'Hourly' ) THEN
       NewCollection%ResetFlag = ResetFlagHourly

    ! Write out all the time 
    ELSEIF ( TRIM(WriteFreq) == 'Always' ) THEN
       NewCollection%ResetFlag = ResetFlagAlways

    ! Write out only at end of simulation
    ELSEIF ( TRIM(WriteFreq) == 'End' ) THEN
       NewCollection%ResetFlag = ResetFlagEnd

    ! Manually write out.
    ELSEIF ( TRIM(WriteFreq) == 'Manual' ) THEN
       NewCollection%ResetFlag = ResetFlagManually

    ! Error otherwise
    ELSE
       MSG = 'Illegal averaging interval: ' // TRIM(WriteFreq) // &
             ' - cannot create diagnostics ' // TRIM(NewCollection%PREFIX)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    NewCollection%WriteFreq = TRIM(WriteFreq)

    ! Add to collections list. Put at the beginning
    NewCollection%NextCollection => Collections
    Collections                  => NewCollection

    ! Define this collection ID
    nnCollections              = nnCollections + 1
    NewCollection%CollectionID = nnCollections
    COL                        = NewCollection%CollectionID 

    ! verbose
    IF ( HCO_IsVerb( 1 ) ) THEN
       MSG = 'Created diagnostics collection: '
       CALL HCO_MSG(MSG)
       WRITE(MSG,'(a21,i2)') ' - Collection ID  : ', COL 
       CALL HCO_MSG(MSG)
       WRITE(MSG,'(a21,a)' ) ' - PREFIX         : ', TRIM(NewCollection%PREFIX)
       CALL HCO_MSG(MSG)
       WRITE(MSG,'(a21,a)' ) ' - Output interval: ', TRIM(NewCollection%WriteFreq)
       CALL HCO_MSG(MSG)
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE DiagnCollection_Create
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DiagnCollection_Cleanup
!
! !DESCRIPTION: Subroutine DiagnCollection\_Cleanup cleans up the diagnostics
! collection for collection COL. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCollection_Cleanup
!
! !INPUT ARGUMENTS:
!
!
! !REVISION HISTORY:
!  08 Jan 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ARGUMENTS:
!
    TYPE(DiagnCollection), POINTER ::  ThisColl => NULL()
    TYPE(DiagnCollection), POINTER ::  NextColl => NULL()

    ! ================================================================
    ! DiagnCollection_Cleanup begins here
    ! ================================================================

    ! Do for every collection in list
    ThisColl => Collections

    DO WHILE ( ASSOCIATED(ThisColl) ) 

       ! Cleanup
       CALL DiagnList_Cleanup( ThisColl%DiagnList )
       ThisColl%nnDiagn = 0
       ThisColl%AREA_M2 => NULL()

       ! Advance
       NextColl                => ThisColl%NextCollection
       ThisColl%NextCollection => NULL()
       ThisColl                => NextColl 
    ENDDO

  END SUBROUTINE DiagnCollection_Cleanup
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DiagnCollection_DefineID 
!
! !DESCRIPTION: Subroutine DiagnCollection\_DefineID is a helper routine to 
! return the collection ID. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCollection_DefineID ( PS, RC, COL, DEF, OkIfAll, InUse, ThisColl ) 
!
! !INPUT ARGUMENTS:
!
    INTEGER,               INTENT(IN   ), OPTIONAL :: COL        ! desired collection number 
    INTEGER,               INTENT(IN   ), OPTIONAL :: DEF        ! default collection number 
    LOGICAL,               INTENT(IN   ), OPTIONAL :: OkIfAll    ! Ok if all (PS=-1) 
!
! !INPUT/OUTPUT ARGUMENTS:
!
    INTEGER,               INTENT(INOUT)           :: PS         ! Assigned collection number 
    INTEGER,               INTENT(INOUT)           :: RC         ! Return code 
    LOGICAL,               INTENT(  OUT), OPTIONAL :: InUse      ! Is this in use? 
!
! !OUTPUT ARGUMENTS:
!
    TYPE(DiagnCollection), POINTER,      OPTIONAL :: ThisColl    ! Pointer to collection 
!
! !REVISION HISTORY:
!  01 Apr 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ARGUMENTS:
!
    LOGICAL                        :: AllOk, FOUND
    CHARACTER(LEN=255)             :: MSG
    CHARACTER(LEN=255), PARAMETER  :: LOC = 'DiagnCollection_DefineID (hco_diagn_mod.F90)'

    ! ================================================================
    ! DiagnCollection_DefineID begins here
    ! ================================================================

    ! Init
    IF ( PRESENT(ThisColl) ) ThisColl => NULL()
    IF ( PRESENT(InUse   ) ) InUse    = .FALSE. 

    ! Check if it's negative
    AllOk = .FALSE.
    IF ( PRESENT(OkIfAll) ) AllOK = OkIfAll 

    ! Get collection position
    IF ( PRESENT(DEF) ) THEN
       PS = DEF
    ELSE
       PS = HcoDiagnIDDefault
    ENDIF 
    IF ( PRESENT(COL) ) PS = COL
 
    ! Check if all collections are selected (-1) 
    IF ( PS == -1 ) THEN
       IF ( AllOK ) THEN
          IF ( PRESENT(InUse)    ) InUse    =  .TRUE.
          IF ( PRESENT(ThisColl) ) ThisColl => Collections
          RC = HCO_SUCCESS
          RETURN
       ELSE
          WRITE(MSG,*) 'Not allowed to select all collections ', PS
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

    ! If individual collection is selected 
    ELSE

       ! Try to find collection
       CALL DiagnCollection_Find( PS, FOUND, RC, ThisColl=ThisColl )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Eventually fill argumnet
       IF ( PRESENT(InUse) ) THEN
          InUse = FOUND 

       ELSEIF ( .NOT. FOUND ) THEN
          WRITE(MSG,*) 'Diagnostics collection not defined: ', PS
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE DiagnCollection_DefineID 
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: DiagnCollection_Find
!
! !DESCRIPTION: Subroutine DiagnCollection\_Find searches the collection
! linked list for the collection with the given collection ID.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCollection_Find ( PS, FOUND, RC, ThisColl ) 
!
! !INPUT ARGUMENTS:
!
    INTEGER,               INTENT(IN   )          :: PS       ! desired collection number 
!
! !INPUT/OUTPUT ARGUMENTS:
!
    LOGICAL,               INTENT(  OUT)          :: FOUND    ! Collection exists?
    INTEGER,               INTENT(INOUT)          :: RC       ! Return code 
    TYPE(DiagnCollection), POINTER,      OPTIONAL :: ThisColl ! Pointer to collection 
!
! !REVISION HISTORY:
!  01 Apr 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ARGUMENTS:
!
    TYPE(DiagnCollection), POINTER :: TmpColl => NULL()
    CHARACTER(LEN=255), PARAMETER  :: LOC = 'DiagnCollection_Find (hco_diagn_mod.F90)'

    ! ================================================================
    ! DiagnCollection_Find begins here
    ! ================================================================

    ! Check if it's negative
    FOUND = .FALSE.

    ! Loop over all collections
    TmpColl => Collections
    DO WHILE ( ASSOCIATED(TmpColl) ) 

       ! Check if this is the collection of insterest
       IF ( TmpColl%CollectionID == PS ) THEN
          FOUND = .TRUE.
          EXIT
       ENDIF

       ! Advance in list
       TmpColl => TmpColl%NextCollection
    ENDDO

    ! Eventually pass to output argument
    IF ( PRESENT(ThisColl) ) ThisColl => TmpColl

    ! Cleanup
    TmpColl => NULL()

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE DiagnCollection_Find
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DiagnFileOpen
!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnFileOpen( am_I_Root, LUN, RC )
!
! !USES:
!
    USE inquireMod,        ONLY : findFreeLUN
    USE HCO_ExtList_Mod,   ONLY : CoreNr, GetExtOpt
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root   ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC          ! Failure or success
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)           :: LUN         ! File LUN 
!
! !REVISION HISTORY: 
!  10 Apr 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: IOS
    LOGICAL             :: EXISTS, FOUND
    CHARACTER(LEN=255)  :: MSG, DiagnFile
    CHARACTER(LEN=255)  :: LOC = 'DiagnFileOpen (hco_diagn_mod.F90)' 

    !=================================================================
    ! DiagnFileOpen begins here!
    !=================================================================
    
    ! Init
    LUN = -1    

    ! Try to get name of diagnostics file
    CALL GetExtOpt ( CoreNr, 'DiagnFile', OptValChar=DiagnFile, &
                     FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Read file and define diagnostics for each entry
    IF ( FOUND ) THEN

       ! Find free LUN
       LUN = findFreeLUN()

       INQUIRE( FILE=TRIM(DiagnFile), EXIST=EXISTS )
       IF ( .NOT. EXISTS ) THEN
          MSG = 'Cannot read file - it does not exist: ' // TRIM(DiagnFile)
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Open configuration file
       OPEN ( LUN, FILE=TRIM( DiagnFile ), STATUS='OLD', IOSTAT=IOS )
       IF ( IOS /= 0 ) THEN
          MSG = 'Error opening ' // TRIM(DiagnFile)
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF 

    ENDIF !FOUND

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE DiagnFileOpen
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DiagnFileGetNext
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnFileGetNext( am_I_Root, LUN,      cName,   &
                               SpcName,   ExtNr,    Cat,     &
                               Hier,      SpaceDim, OutUnit, & 
                               EOF,       RC ) 
!
! !USES:
!
    USE HCO_CharTools_Mod
    USE CHARPAK_Mod,       ONLY : STRREPL, STRSPLIT
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root   ! root CPU?
    INTEGER,          INTENT(IN   )           :: LUN         ! file LUN 
!
! !INPUT/OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(INOUT)           :: EOF
    INTEGER,          INTENT(INOUT)           :: RC          ! Failure or success
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(  OUT)           :: cName  
    CHARACTER(LEN=*), INTENT(  OUT)           :: SpcName
    INTEGER,          INTENT(  OUT)           :: ExtNr 
    INTEGER,          INTENT(  OUT)           :: Cat 
    INTEGER,          INTENT(  OUT)           :: Hier
    INTEGER,          INTENT(  OUT)           :: SpaceDim 
    CHARACTER(LEN=*), INTENT(  OUT)           :: OutUnit
!
! !REVISION HISTORY: 
!  10 Apr 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: N
    CHARACTER(LEN=255)  :: LINE
    CHARACTER(LEN=255)  :: MSG
    CHARACTER(LEN=255)  :: SUBSTR(255) 
    CHARACTER(LEN=255)  :: LOC = 'DiagnFileGetNext (hco_diagn_mod.F90)'

    !=================================================================
    ! DiagnFileGetNext begins here!
    !=================================================================

    ! Init
    cName    = ''
    SpcName  = '' 
    OutUnit  = ''
    ExtNr    = -1
    Cat      = -1
    Hier     = -1
    SpaceDim = -1

    ! Get next line
    CALL GetNextLine ( am_I_Root, LUN, LINE, EOF, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave here if end of file
    IF ( .NOT. EOF ) THEN 

       ! Parse diagnostics information from line
       CALL STRREPL( LINE, HCO_TAB(), HCO_SPC() )

       ! Split into substrings
       CALL STRSPLIT( LINE, HCO_SPC(), SUBSTR, N ) 

       ! There must be at least 7 entries
       IF ( N < 7 ) THEN 
          MSG = 'Diagnostics entries must have 7 elements: '// TRIM(LINE)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Extract diagnostics properties
       cName   = TRIM(SUBSTR(1))
       SpcName = TRIM(SUBSTR(2))

       ! Extension number, category, hierarchy, space dimension
       READ( SUBSTR(3), * ) ExtNr
       READ( SUBSTR(4), * ) Cat 
       READ( SUBSTR(5), * ) Hier
       READ( SUBSTR(6), * ) SpaceDim

       ! Read output unit
       OutUnit = TRIM(SUBSTR(7))

    ENDIF !EOF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE DiagnFileGetNext
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DiagnFileClose
!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnFileClose ( LUN )
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: LUN         ! File LUN 
!
! !REVISION HISTORY: 
!  10 Apr 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
    CLOSE ( LUN ) 

  END SUBROUTINE DiagnFileClose
!EOC
END MODULE HCO_Diagn_Mod
