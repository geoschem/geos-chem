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
! hierarchy, species ID), data structure (Scalar, 2D, 3D), output
! units (mass, area, time), and output frequency (every hour / day / 
! month / year).
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
!\\
! The HEMCO diagnostics module can store multiple, independent diagnostics 
! `collections`, identifiable through the assigned collection number. 
! Before adding diagnostics to a collection, the collection needs to be 
! created using subroutine DiagnCollection\_Create. The collection number 
! argument (COL) should always be specified when creating, editing or 
! obtaining a diagnostics. If this argument is omitted, the collection
! number is set to 1.
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
  PUBLIC  :: HCO_Diagn_AutoUpdate
  PUBLIC  :: Diagn_Create
  PUBLIC  :: Diagn_Update 
  PUBLIC  :: Diagn_Get
  PUBLIC  :: Diagn_AutoFillLevelDefined
  PUBLIC  :: Diagn_GetMaxResetFlag
  PUBLIC  :: Diagn_GetDiagnPrefix
  PUBLIC  :: Diagn_Print
  PUBLIC  :: DiagnCont_Find
  PUBLIC  :: DiagnCollection_Create
  PUBLIC  :: DiagnCollection_Cleanup
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: DiagnList_Cleanup 
  PRIVATE :: DiagnCont_Init
  PRIVATE :: DiagnCont_PrepareOutput
  PRIVATE :: DiagnCont_Link_2D
  PRIVATE :: DiagnCont_Link_3D
  PRIVATE :: DiagnCont_Cleanup
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller   - Initialization
!  08 Jul 2014 - R. Yantosca - Now use F90 free-format indentation  
!  08 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  01 Aug 2014 - C. Keller   - Added manual output frequency
!  12 Aug 2014 - C. Keller   - Added cumulative sum option
!  09 Jan 2015 - C. Keller   - Added diagnostics collections
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
     INTEGER                     :: cID            ! Cont. ID
     INTEGER                     :: ExtNr          ! Extension #
     INTEGER                     :: Cat            ! Category 
     INTEGER                     :: Hier           ! Hierarchy
     INTEGER                     :: HcoID          ! HEMCO species ID
     INTEGER                     :: AutoFill       ! fill automatically? 
     INTEGER                     :: SpaceDim       ! Space dimension (1-3) 
     REAL(hp)                    :: Scalar         ! 1D scalar 
     TYPE(Arr2D_HP),     POINTER :: Arr2D          ! 2D array
     TYPE(Arr3D_HP),     POINTER :: Arr3D          ! 3D array
     LOGICAL                     :: DtaIsPtr       ! Is data just a pointer?
     INTEGER                     :: LevIdx         ! Level index to be used 
     CHARACTER(LEN= 31)          :: OutUnit        ! Output unit 
     INTEGER                     :: AreaFlag       ! 2=per area, 3=per volume, 0 otherwise 
     REAL(hp)                    :: AreaScal       ! Scale factor for area
     REAL(hp)                    :: MassScal       ! Scale factor for mass
     INTEGER                     :: TimeAvg        ! Scale flag for time unit 
     INTEGER                     :: Counter        ! time steps since 
                                                   !  last output
     INTEGER                     :: AvgFlag        ! Averaging flag for 
                                                   !  non-standard units
     INTEGER                     :: ResetFlag      ! Diagn. output frequency
     INTEGER                     :: LastUpdateID   ! Last update time
     INTEGER                     :: nnGetCalls     ! # of Diagn_Get calls w/o update 
     LOGICAL                     :: IsOutFormat    ! Data is in output format?
     TYPE(DiagnCont),    POINTER :: NextCont       ! Ptr to next item in list
  END TYPE DiagnCont

  !------------------------------------------------------------------------
  ! Diagnostcs collection derived type.
  ! DiagnList      : Linked list with all diagnostics container of
  !                  this collection.
  ! nnDiag         : Number of diagnostics in this collection.
  ! MaxResetFlag   : Highest reset flag found in this collection.
  ! AF_LevelDefined: Set to true if there is at least one autofill 
  !                  diagnostics at the given level (1-4).
  ! InUse          : Is this collection in use?
  ! PREFIX         : Prefix to be used for diagnostics output file name.
  ! NX, NY, NZ     : Grid dimensions.
  ! TS             : Time step. This is only of relevance for emission 
  !                  diagnostics that are internally converted from
  !                  kg/m2/s to kg/m2.
  ! AREA_M2        : Surface grid box areas. May be required for unit 
  !                  conversions.
  !------------------------------------------------------------------------
  TYPE DiagnCollection
     TYPE(DiagnCont),    POINTER :: DiagnList          => NULL()
     INTEGER                     :: nnDiagn            =  0
     INTEGER                     :: MaxResetFlag       =  ResetFlagManually
     LOGICAL                     :: AF_LevelDefined(4) =  .FALSE.
     LOGICAL                     :: InUse              =  .FALSE.
     CHARACTER(LEN=255)          :: PREFIX             =  ''
     INTEGER                     :: NX                 =  0
     INTEGER                     :: NY                 =  0
     INTEGER                     :: NZ                 =  0
     REAL(sp)                    :: TS                 =  0       ! Time step
     REAL(hp),           POINTER :: AREA_M2(:,:)       => NULL()
  END TYPE DiagnCollection

  ! Vector of diagnostic collections. The maximum number of collections
  ! to be used is defined below. 
  TYPE(DiagnCollection),  POINTER :: Collections(:) => NULL()
!
! !DEFINED PARAMETERS:
!
  ! Parameter for averaging and summing non-standard data
  ! AvgFlagMean    : calculates the arithmetic mean
  ! AvgFlagSum     : calculates the sum, resets after every writeout
  ! AvgFlagCumSum  : calculates the cumulative sum, never resets.
  ! AvgFlagInst    : uses the instantaneous value, overwrites existing 
  INTEGER, PARAMETER             :: AvgFlagMean    = 1
  INTEGER, PARAMETER             :: AvgFlagSum     = 2
  INTEGER, PARAMETER             :: AvgFlagCumsum  = 3
  INTEGER, PARAMETER             :: AvgFlagInst    = 4

  ! Maximum number of diagnostics collections
  INTEGER, PARAMETER             :: MaxCollections = 3

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: hco_diagn_autoupdate
!
! !DESCRIPTION: Subroutine HCO\_DIAGN\_AUTOUPDATE updates the AutoFill
! diagnostics at species level. This routine should be called after
! running HEMCO core and all extensions. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_Diagn_AutoUpdate( am_I_Root, HcoState, RC ) 
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_GetHcoID
    USE HCO_STATE_MOD, ONLY : HCO_State

    ! temp only
    USE HCO_ARR_MOD,   ONLY : HCO_ArrAssert
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
    ! HCO_DIAGN_AUTOUPDATE begins here!
    !=================================================================
    
    ! Init 
    LOC = 'HCO_DIAGN_AUTOUPDATE (hco_diagn_mod.F90)'
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
                                RC       = RC          ) 
             IF ( RC/= HCO_SUCCESS ) RETURN 
             Arr3D => NULL() 
          ENDIF
       ENDIF
    ENDDO
    
    ! Return
    RC = HCO_SUCCESS
    
  END SUBROUTINE HCO_Diagn_AutoUpdate
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
!\item WriteFreq: output frequency. Can be one of 'Hourly', 'Daily',
!      'Monthly', 'Annualy', 'End', 'Manual'.
!      Manual diagnostics are expected to be manually received and 
!      written out. These diagnostics may or may not be written out 
!      at the end of the simulation run, depending on the corresponding
!      attribute set in Diagn\_Get.
!\item HcoState: HEMCO state object. Used to determine the species
!      properties if any of arguments MW\_g, EmMW\_g or MolecRatio
!      is missing.
!\item OutOper: output operation for non-standard units. If this 
!      argument is used, the specified operation is performed and all
!      unit specifications are ignored. Can be one of 'Mean', 'Sum', 
!      'Cumsum', or 'Instantaneous'.
!\item AutoFill: containers with an AutoFill flag of 1 will be auto-
!      matically updated by the HEMCO standard diagnostics calls
!      (e.g. in hco\_calc\_mod.F90). If set to 0, the diagnostics 
!      updates have to be set manually.
!\item Trgt2D: 2D target array. If specified, the diagnostics array
!      will point to this data. This disables all time averaging, 
!      unit conversions, etc., and the data will be written to disk
!      as is.
!\item Trgt3D: as Trgt2D, but for 3D data. 
!\item MW_g: species molecular weight. Used to determine unit
!      conversion factors. Not needed for target containers or if
!      argument OutOper is specified. Can be omitted if HcoState is
!      given.
!\item EmMW_g: Molecular weight of emitted species. Used to determine
!      unit conversion factors. Not needed for target containers or if
!      argument OutOper is specified. Can be omitted if HcoState is
!      given.
!\item MolecRatio: Molecules of species per emitted molecule. Used to 
!      determine unit conversion factors. Not needed for target 
!      containers or if argument OutOper is specified. Can be omitted
!      if HcoState is given.
!\item cID: assigned container ID. Useful for later reference to this
!      diagnostics container.
!\item RC: HEMCO return code.
!\end{itemize} 
!
! !INTERFACE:
!
  SUBROUTINE Diagn_Create( am_I_Root, cName,     HcoState,   &
                           ExtNr,     Cat,       Hier,       &
                           HcoID,     SpaceDim,  OutUnit,    &  
                           WriteFreq, OutOper,   LevIdx,     &
                           AutoFill,  Trgt2D,    Trgt3D,     &
                           MW_g,      EmMW_g,    MolecRatio, &
                           cID,       RC,        COL          )
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
    CHARACTER(LEN=*), INTENT(IN   )           :: cName         ! Diagnostics 
                                                               !  name
    INTEGER,          INTENT(IN   )           :: ExtNr         ! Extension #    
    INTEGER,          INTENT(IN   )           :: Cat           ! Category 
    INTEGER,          INTENT(IN   )           :: Hier          ! Hierarchy 
    INTEGER,          INTENT(IN   )           :: HcoID         ! HEMCO species 
                                                               !  ID # 
    INTEGER,          INTENT(IN   )           :: SpaceDim      ! Spatial 
                                                               !  dimension 
    CHARACTER(LEN=*), INTENT(IN   )           :: OutUnit       ! Output units
    CHARACTER(LEN=*), INTENT(IN   )           :: WriteFreq     ! Write out 
                                                               !  frequency
    TYPE(HCO_State),  POINTER,       OPTIONAL :: HcoState      ! HEMCO state obj.
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: OutOper       ! Output 
                                                               !  operation 
    INTEGER,          INTENT(IN   ), OPTIONAL :: LevIdx        ! Level index 
                                                               !  to use
    INTEGER,          INTENT(IN   ), OPTIONAL :: AutoFill      ! 1=auto fill
                                                               ! 0=don't 
    REAL(hp),         INTENT(IN   ), OPTIONAL :: Trgt2D(:,:)   ! 2D target data
    REAL(hp),         INTENT(IN   ), OPTIONAL :: Trgt3D(:,:,:) ! 3D target data
    REAL(hp),         INTENT(IN   ), OPTIONAL :: MW_g          ! species MW (g/mol) 
    REAL(hp),         INTENT(IN   ), OPTIONAL :: EmMW_g        ! emission MW (g/mol)
    REAL(hp),         INTENT(IN   ), OPTIONAL :: MolecRatio    ! molec. emission ratio
    INTEGER,          INTENT(IN   ), OPTIONAL :: COL           ! Collection number 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)           :: cID           ! Assigned 
                                                               !  container ID
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC            ! Return code 
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
    TYPE(DiagnCont), POINTER :: ThisDiagn => NULL()
    TYPE(DiagnCont), POINTER :: TmpDiagn  => NULL()

    ! Scalars
    CHARACTER(LEN=255)       :: LOC, MSG
    INTEGER                  :: PS, Flag
    REAL(hp)                 :: Scal
    REAL(hp)                 :: MWg, EmMWg, MolR 
    LOGICAL                  :: ForceMean, FOUND

    !======================================================================
    ! Diagn_Create begins here!
    !======================================================================

    ! Init
    LOC = 'Diagn_Create (hco_diagn_mod.F90)'

    ! DEBUGGING - ewl, 2/6/15
    PRINT *, " "
    PRINT *, "In ", TRIM( LOC )
    ! END DEBUGGING

    ! Get collection position
    PS = 1
    IF ( PRESENT(COL) ) PS = COL

    IF ( .NOT. Collections(PS)%InUse ) THEN
       WRITE(MSG,*) 'Diagnostics collection not defined, cannot add ',  &
                    'diagnostcs ', TRIM(cName), ' to collection ', PS,  &
                    '. Please call DiagnCollection_Create before ',     &
                    'adding diagnostics to that collection.'
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    !----------------------------------------------------------------------
    ! Initalize diagnostics container. This will automatically add the
    ! container to the diagnostics list.
    !----------------------------------------------------------------------
    CALL DiagnCont_Init( ThisDiagn )

    !----------------------------------------------------------------------
    ! Pass input variables
    !----------------------------------------------------------------------
    ThisDiagn%cName    = cName
    ThisDiagn%ExtNr    = ExtNr
    ThisDiagn%Cat      = Cat
    ThisDiagn%Hier     = Hier
    ThisDiagn%HcoID    = HcoID
    ThisDiagn%SpaceDim = SpaceDim
    ThisDiagn%OutUnit  = TRIM(OutUnit)
    IF ( PRESENT(LevIdx)   ) ThisDiagn%LevIdx   = LevIdx
    IF ( PRESENT(AutoFill) ) ThisDiagn%AutoFill = AutoFill 

    !----------------------------------------------------------------------
    ! Eventually link to data array. This will disable all time averaging,
    ! unit conversions, etc. (data will just be returned as is). 
    !----------------------------------------------------------------------
    IF ( PRESENT(Trgt2D) ) THEN
       CALL DiagnCont_Link_2D( am_I_Root, ThisDiagn, Trgt2D, PS, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF
    IF ( PRESENT(Trgt3D) ) THEN
       CALL DiagnCont_Link_3D( am_I_Root, ThisDiagn, Trgt3D, PS, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    !----------------------------------------------------------------------
    ! Determine output frequency. This is the frequency with which the
    ! diagnostics will be written into output. 
    !----------------------------------------------------------------------

    ! Write out every Year
    IF ( TRIM(WriteFreq) == 'Annually' ) THEN
       ThisDiagn%ResetFlag = ResetFlagAnnually 

    ! Write out every month
    ELSEIF ( TRIM(WriteFreq) == 'Monthly' ) THEN
       ThisDiagn%ResetFlag = ResetFlagMonthly
         
    ! Write out every day
    ELSEIF ( TRIM(WriteFreq) == 'Daily' ) THEN
       ThisDiagn%ResetFlag = ResetFlagDaily

    ! Write out every hour
    ELSEIF ( TRIM(WriteFreq) == 'Hourly' ) THEN
       ThisDiagn%ResetFlag = ResetFlagHourly

    ! Write out only at end of simulation
    ELSEIF ( TRIM(WriteFreq) == 'End' ) THEN
       ThisDiagn%ResetFlag = ResetFlagEnd

    ! Manually write out.
    ELSEIF ( TRIM(WriteFreq) == 'Manual' ) THEN
       ThisDiagn%ResetFlag = ResetFlagManually

    ! Error otherwise
    ELSE
       MSG = 'Illegal averaging interval: ' // TRIM(WriteFreq)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Update module variable MaxResetFlag. This variable defines the 
    ! highest reset flag used by any of the diagnostics container.
    Collections(PS)%MaxResetFlag = &
       MAX( Collections(PS)%MaxResetFlag, ThisDiagn%ResetFlag )

    ! Update module variable AF_LevelDefined. For all AutoFill diagnostics,
    ! we store whether or not there is (at least one) diagnostics container
    ! defined at species level, ExtNr level, etc. 
    IF ( ThisDiagn%AutoFill == 1 ) THEN

       ! At species level: no ExtNr defined
       IF ( ThisDiagn%ExtNr < 0 ) THEN
          Collections(PS)%AF_LevelDefined(1) = .TRUE.
     
       ! At ExtNr level: no category defined
       ELSEIF ( ThisDiagn%Cat < 0 ) THEN
          Collections(PS)%AF_LevelDefined(2) = .TRUE.

       ! At category level: no hierarchy defined
       ELSEIF ( ThisDiagn%Hier < 0 ) THEN
          Collections(PS)%AF_LevelDefined(3) = .TRUE.

       ! At hierarchy level: all defined
       ELSE
          Collections(PS)%AF_LevelDefined(4) = .TRUE.
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
    IF ( .NOT. ThisDiagn%DtaIsPtr ) THEN
 
       ! Enforce specified output operator 
       IF ( PRESENT(OutOper) ) THEN
          IF ( TRIM(OutOper) == 'Mean' ) THEN
             ThisDiagn%AvgFlag = AvgFlagMean
          ELSEIF ( TRIM(OutOper) == 'Sum' ) THEN
             ThisDiagn%AvgFlag = AvgFlagSum
          ELSEIF ( TRIM(OutOper) == 'Cumsum' ) THEN
             ThisDiagn%AvgFlag = AvgFlagCumsum
          ELSEIF ( TRIM(OutOper) == 'Instantaneous' ) THEN
             ThisDiagn%AvgFlag = AvgFlagInst
          ELSE
             MSG = 'Illegal output operator: ' // TRIM(OutOper)
             MSG = TRIM(MSG) // '. Allowed are `Mean`, `Sum`, '// &
                   '`Cumsum`, `Instantaneous`.'
             MSG = TRIM(MSG) // ' (' // TRIM(cName) // ')'
             CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          ENDIF

          ! In an ESMF/MAPL environment, treat all data as instantaneous.
          ! The HEMCO diagnostics are expected to be passed to the MAPL
          ! history component on every time step and we should only pass
          ! the instantaneous values so that the history component can 
          ! properly perform its own data operations (based on the settings
          ! in HISTORY.rc). 
#if defined(ESMF_)
          ThisDiagn%AvgFlag = AvgFlagInst
#endif
 
       ! If OutOper is not set, determine scale factors from output unit:
       ELSE
   
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
    ! Make sure that there is no duplicate entry (for AutoFill only)
    ! Now allow multiple autofill diagnostics with same ExtNr, Cat, Hier, 
    ! and species ID (ckeller, 09/25/2014).
    !-----------------------------------------------------------------------
    CALL DiagnCont_Find( -1, -1, -1, -1, -1, &
                         Trim(cName), 1, FOUND, TmpDiagn, COL=PS )
    IF ( FOUND ) THEN
       MSG = 'There is already a diagnostics with this name: ' // TRIM(cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Add to diagnostics list of this collection. 
    ! Insert at the beginning of the list.
    !-----------------------------------------------------------------------
    
    ! DEBUGGING - ewl, 2/2/15
    PRINT *, "   Creating container for diag ", TRIM( cName )
    ! END DEBUGGING

    IF ( Collections(PS)%nnDiagn > 0 ) THEN
       ThisDiagn%NextCont => Collections(PS)%DiagnList
       
       ! DEBUGGING - ewl, 2/2/15
       PRINT *, "   NextCont now points to previous diag"
       ! END DEBUGGING
    ENDIF
    Collections(PS)%DiagnList => ThisDiagn

    ! Increase diagnostics counter and set container ID accordingly.
    Collections(PS)%nnDiagn = Collections(PS)%nnDiagn + 1
    ThisDiagn%cID           = Collections(PS)%nnDiagn

    ! Verbose mode
    IF ( am_I_Root .AND. HCO_VERBOSE_CHECK() ) THEN
       WRITE(MSG,*) 'Successfully added diagnostics to collection ' , PS
       CALL HCO_MSG ( MSG )
       CALL Diagn_Print( ThisDiagn, .TRUE. )
    ENDIF

    ! DEBUGGING - ewl, 2/2/15
    PRINT *, "   Collection: ", PS 
    PRINT *, "   Index in collection: ", ThisDiagn%cID
    PRINT *, "   Counter: ", ThisDiagn%Counter
    PRINT *, "Exiting ", TRIM( LOC ), " for ", TRIM( ThisDiagn%cName )
    ! END DEBUGGING

    ! Return
    cID = ThisDiagn%cID 
    RC  = HCO_SUCCESS
    ThisDiagn => NULL()

  END SUBROUTINE Diagn_Create
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_Update
!
! !DESCRIPTION: Subroutine Diagn\_Update updates the content of a 
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
! Note that for a given time step, the same diagnostics container can 
! be updated multiple times.
!\\
!\\
! If the passed array is empty (i.e. not associated), it is treated as
! empty values (i.e. zeros).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Update( am_I_Root, cID,      cName,            &
                           ExtNr,     Cat,      Hier,    HcoID,   &
                           AutoFill,  Scalar,   Array2D, Array3D, &
                           PosOnly,   COL,      RC                 ) 
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
    USE HCO_Arr_Mod,   ONLY : HCO_ArrAssert
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root      ! Root CPU?
    INTEGER,          INTENT(IN   ), OPTIONAL :: cID            ! Assigned 
                                                                !  container ID
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: cName          ! Diagnostics 
                                                                !  name
    INTEGER,          INTENT(IN   ), OPTIONAL :: ExtNr          ! Extension #
    INTEGER,          INTENT(IN   ), OPTIONAL :: Cat            ! Category 
    INTEGER,          INTENT(IN   ), OPTIONAL :: Hier           ! Hierarchy 
    INTEGER,          INTENT(IN   ), OPTIONAL :: HcoID          ! HEMCO species
                                                                !  ID number 
    INTEGER,          INTENT(IN   ), OPTIONAL :: AutoFill       ! 1=yes; 0=no; 
                                                                ! -1=either 
    REAL(hp),         INTENT(IN   ), OPTIONAL :: Scalar         ! 1D scalar 
    REAL(hp),         POINTER,       OPTIONAL :: Array2D(:,:)   ! 2D array 
    REAL(hp),         POINTER,       OPTIONAL :: Array3D(:,:,:) ! 3D array 
    LOGICAL,          INTENT(IN   ), OPTIONAL :: PosOnly        ! Use only vals
                                                                !  >= 0?
    INTEGER,          INTENT(IN   ), OPTIONAL :: COL            ! Collection Nr.
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC             ! Return code 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!  25 Sep 2014 - C. Keller: Now allow updating multiple diagnostics
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(DiagnCont), POINTER :: ThisDiagn  => NULL()
    REAL(hp),        POINTER :: Arr2D(:,:) => NULL()

    ! Scalars
    CHARACTER(LEN=255)       :: LOC, MSG
    REAL(hp)                 :: Fact
    REAL(hp)                 :: Tmp
    CHARACTER(LEN=31)        :: DgnName
    INTEGER                  :: I, J, L
    INTEGER                  :: DgncID,  DgnExtNr, DgnCat
    INTEGER                  :: DgnHier, DgnHcoID
    INTEGER                  :: MinResetFlag, ThisUpdateID
    INTEGER                  :: AutoFlag
    INTEGER                  :: PS
    LOGICAL                  :: Found, OnlyPos, VertSum, IsAssoc, IsNewTS

    !======================================================================
    ! Diagn_Update begins here!
    !======================================================================

    ! Init
    LOC = 'Diagn_Update (hco_diagn_mod.F90)'
    RC  = HCO_SUCCESS

    ! DEBUGGING, ewl - 2/2/15
    PRINT *, " "
    PRINT *, "Now in " // TRIM( LOC )
    ! END DEBUGGING

    ! Get collection number
    PS = 1
    IF ( PRESENT(COL) ) PS = COL
    IF ( PS > MaxCollections ) THEN
       WRITE(MSG,*) 'Illegal diagnostics collection number:', PS
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Nothing to do if this collection is empty
    IF ( .NOT. Collections(PS)%InUse ) RETURN 

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
  
    !-----------------------------------------------------------------
    ! Do for every container in the diagnostics list that matches the 
    ! specified arguments (ID, ExtNr, etc.). This can be more than one
    ! container (ckeller, 09/25/2014).
    !-----------------------------------------------------------------

    ! DEBUGGING - ewl, 2/2/15
    IF ( PRESENT( cName ) ) THEN
       PRINT *, "   Calling DiagnCont_Find for diagnostic " // TRIM( DgnName )
    ELSE
       PRINT *, "   Calling DiagnCont_Find for unnamed diagnostic."
    ENDIF
    PRINT *, "      DgncID: ", DgncID
    PRINT *, "      DgnExtNr: ", DgnExtNr
    PRINT *, "      DgnCat: ",DgnCat
    PRINT *, "      DgnHier: ", DgnHier
    PRINT *, "      DgnHcoID: ", DgnHcoID
    PRINT *, "      OnlyPos: ", OnlyPos
    PRINT *, "      AutoFlag: ", AutoFlag
    ! END DEBUGGING

    DO

       ! Search for diagnostics that matches the given arguments.
       ! If ThisDiagn is empty (first call), the search will start
       ! at the first diagnostics container. Otherwise, the search
       ! will resume from this diagnostics container.


       CALL DiagnCont_Find( DgncID,    DgnExtNr, DgnCat,   DgnHier, &
                            DgnHcoID,  DgnName,  AutoFlag, Found,   &
                            ThisDiagn, RESUME=.TRUE., COL=PS         )

       ! Exit while loop if no diagnostics found
       IF ( .NOT. Found ) THEN

          ! DEBUGGING - ewl, 2/2/15
          PRINT *, "   No other diagnostic found. Exiting loop over containers."
          ! END DEBUGGING

          EXIT

       ENDIF

       ! DEBUGGING - ewl, 2/5/15
       PRINT *, "   Diagnostic found. Counter = ", ThisDiagn%Counter
       ! END DEBUGGING

       ! If container holds just a pointer to external data, don't do
       ! anything!
       IF ( ThisDiagn%DtaIsPtr ) THEN
          MSG = 'You try to update a container that holds a '  // &
               'pointer to data - this should never happen! ' // &
               TRIM(ThisDiagn%cName)
          CALL HCO_WARNING( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
   
       !----------------------------------------------------------------------
       ! Sanity check: if data is beyond its averaging interval, it should
       ! be in output format. Otherwise, the content of this diagnostics has 
       ! never been passed to the output yet (via routine Diagn\_Get) and 
       ! will thus be lost!
       !----------------------------------------------------------------------
       IF ( ThisDiagn%ResetFlag >= MinResetFlag .AND. &
            .NOT. ThisDiagn%IsOutFormat ) THEN
          MSG = 'Diagnostics is at end of its output interval '      // &
                'but was not passed to output - data will be lost: ' // &
                TRIM(ThisDiagn%cName)
          CALL HCO_WARNING( MSG, RC, THISLOC=LOC )
       ENDIF
   
       !----------------------------------------------------------------------
       ! If data is in output format, set counter to zero. This will make
       ! sure that the new data is not added to the existing data.
       !----------------------------------------------------------------------
       IF ( ThisDiagn%IsOutFormat ) THEN
          ThisDiagn%Counter    = 0

          ! DEBUGGING - ewl, 2/2/15
          PRINT *, "   ThisDiagn%IsOutFormat is true. Counter reset to 0."
          ! END DEBUGGING
          
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
          Fact = Collections(PS)%TS
       ENDIF
       
       !----------------------------------------------------------------------
       ! Check if this is a new time step for this diagnostics. 
       !----------------------------------------------------------------------
       IsNewTS = .TRUE.
       IF ( ThisDiagn%LastUpdateID == ThisUpdateID ) THEN

          ! DEBUGGING - ewl, 2/2/15
          PRINT *, "   Not a new time step. Setting IsNewTS to FALSE."
          ! END DEBUGGING
          
          IsNewTS = .FALSE. 
          
       ENDIF

       !----------------------------------------------------------------------
       ! To add 3D array
       !----------------------------------------------------------------------
       IF ( ThisDiagn%SpaceDim == 3 ) THEN
   
          ! Make sure dimensions agree and diagnostics array is allocated
          IF ( .NOT. PRESENT(Array3D) ) THEN
             MSG = 'Diagnostics must be 3D: ' // TRIM(ThisDiagn%cName)
             CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
          CALL HCO_ArrAssert( ThisDiagn%Arr3D,    Collections(PS)%NX,   &
                              Collections(PS)%NY, Collections(PS)%NZ, RC ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 
            
          ! Pass array to diagnostics: reset to zero if counter 
          ! is zero, add to it otherwise.
          ! Never reset containers with cumulative sums!
          IF ( ThisDiagn%Counter == 0 .AND. &
               ThisDiagn%AvgFlag /= AvgFlagCumsum ) ThisDiagn%Arr3D%Val = 0.0_hp

          ! Always reset containers with instantaneous values if it's a new
          ! time step.
          IF ( ThisDiagn%AvgFlag == AvgFlagInst .AND. IsNewTS ) ThisDiagn%Arr3D%Val = 0.0_hp

          ! Only if associated ...
          IF ( ASSOCIATED(Array3D) ) THEN
             IF ( OnlyPos ) THEN
                WHERE ( Array3D >= 0.0_hp )
                   ThisDiagn%Arr3D%Val = ThisDiagn%Arr3D%Val + ( Array3D * Fact )
                END WHERE
             ELSE
                ThisDiagn%Arr3D%Val = ThisDiagn%Arr3D%Val + ( Array3D * Fact )
             ENDIF
          ENDIF
   
       !----------------------------------------------------------------------
       ! To add 2D array
       !----------------------------------------------------------------------
       ELSEIF ( ThisDiagn%SpaceDim == 2 ) THEN
   
          ! Make sure dimensions agree and diagnostics array is allocated
          CALL HCO_ArrAssert( ThisDiagn%Arr2D,    Collections(PS)%NX, &
                              Collections(PS)%NY, RC                   ) 
          IF ( RC /= HCO_SUCCESS ) RETURN 
            
          ! Pass array to diagnostics: ignore existing data if counter 
          ! is zero, add to it otherwise.
          ! Never reset containers with cumulative sums!
          IF ( ThisDiagn%Counter == 0 .AND. &
               ThisDiagn%AvgFlag /= AvgFlagCumsum ) ThisDiagn%Arr2D%Val = 0.0_hp
   
          ! Always reset containers with instantaneous values if it's a new time step
          IF ( ThisDiagn%AvgFlag == AvgFlagInst .AND. IsNewTS ) ThisDiagn%Arr2D%Val = 0.0_hp

          ! Assume that we don't have to take the vertical sum
          VertSum = .FALSE.
   
          ! Assume data pointer is associated
          IsAssoc = .TRUE.
   
          ! Convert 3D array to 2D if necessary - only use first level!!
          IF ( PRESENT(Array2D) ) THEN
             IF ( .NOT. ASSOCIATED(Array2D) ) THEN
                IsAssoc = .FALSE.
             ELSE
                Arr2D => Array2D
             ENDIF
          ELSEIF ( PRESENT(Array3D) ) THEN
             IF ( .NOT. ASSOCIATED(Array3D) ) THEN
                IsAssoc = .FALSE.
             ELSE
                IF ( ThisDiagn%LevIdx == -1 ) THEN
                   VertSum = .TRUE.
                ELSE
                   Arr2D => Array3D(:,:,ThisDiagn%LevIdx)
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
                   DO J=1,Collections(PS)%NY
                   DO I=1,Collections(PS)%NX
                      TMP = 0.0_hp
                      DO L=1,Collections(PS)%NZ
                         IF ( Array3D(I,J,L) >= 0.0_hp ) &
                            TMP = TMP + ( Array3D(I,J,L) * Fact )
                      ENDDO
                      ThisDiagn%Arr2D%Val(I,J) = &
                         ThisDiagn%Arr2D%Val(I,J) + TMP
                   ENDDO
                   ENDDO
      
                ! no vertical summation
                ELSE
                   WHERE ( Arr2D >= 0.0_hp )
                      ThisDiagn%Arr2D%Val = ThisDiagn%Arr2D%Val + ( Arr2D * Fact )
                   END WHERE
                ENDIF
      
             ! all values
             ELSE
       
                ! need to do vertical summation
                IF ( VertSum ) THEN
                   DO J=1,Collections(PS)%NY
                   DO I=1,Collections(PS)%NX
                      TMP = SUM(Array3D(I,J,:)) * Fact
                      ThisDiagn%Arr2D%Val(I,J) = &
                         ThisDiagn%Arr2D%Val(I,J) + TMP
                   ENDDO
                   ENDDO
      
                ! no vertical summation
                ELSE
                   ThisDiagn%Arr2D%Val = ThisDiagn%Arr2D%Val + ( Arr2D * Fact )
                ENDIF
             ENDIF
          ENDIF ! pointer is associated
   
       !----------------------------------------------------------------------
       ! To add scalar (1D) 
       !----------------------------------------------------------------------
       ELSEIF ( ThisDiagn%SpaceDim == 1 ) THEN
   
          ! Make sure dimensions agree and diagnostics array is allocated
          IF ( .NOT. PRESENT(Scalar) ) THEN
             MSG = 'Diagnostics must be scalar: '// TRIM(ThisDiagn%cName)
             CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF

          ! Pass array to diagnostics: ignore existing data if counter 
          ! is zero, add to it otherwise.
          ! Never reset containers with cumulative sums!
          IF ( ThisDiagn%Counter == 0 .AND. &
               ThisDiagn%AvgFlag /= AvgFlagCumsum ) ThisDiagn%Scalar = 0.0_hp

          ! Always reset containers with instantaneous values if it's a new time step
          IF ( ThisDiagn%AvgFlag == AvgFlagInst .AND. IsNewTS ) ThisDiagn%Scalar = 0.0_hp

          ! Update scalar value
          IF ( OnlyPos ) THEN
             IF ( Scalar >= 0.0_hp ) & 
                ThisDiagn%Scalar = ThisDiagn%Scalar + ( Scalar * Fact )
          ELSE
             ThisDiagn%Scalar = ThisDiagn%Scalar + ( Scalar * Fact )  
          ENDIF
   
       !----------------------------------------------------------------------
       ! Diangostics space dimension must be 1-3. 
       !----------------------------------------------------------------------
       ELSE
          WRITE(MSG,*) 'Space dimension must be 1-3: ',TRIM(ThisDiagn%cName), &
                       '--> dimension is ', ThisDiagn%SpaceDim
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
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

          ! DEBUGGING - ewl, 2/5/15
          PRINT *, "   IsNewTs is TRUE. New counter: ", ThisDiagn%Counter
          ! END DEBUGGING

       ENDIF
   
       !----------------------------------------------------------------------
       ! Data is not in output format and hasn't been called yet by Diagn_Get.
       !----------------------------------------------------------------------
       ThisDiagn%IsOutFormat = .FALSE.
       ThisDiagn%nnGetCalls  = 0

       ! Verbose mode 
       IF ( am_I_Root .AND. HCO_VERBOSE_CHECK() ) THEN
          WRITE(MSG,'(a,a,a,I3,a)') 'Successfully updated diagnostics: ', &
             TRIM(ThisDiagn%cName), ' (counter:', ThisDiagn%Counter, ')'
!          MSG = 'Successfully updated diagnostics: ' // TRIM(ThisDiagn%cName)
          CALL HCO_MSG ( MSG )
       ENDIF

    ENDDO

    ! DEBUGGING - ewl, 2/2/15
    IF ( PRESENT ( cName ) ) THEN
       PRINT *, "Exiting Diagn_Update for diagnostic: " // TRIM( cName )
    ELSE
       PRINT *, "Exiting Diagn_Update for unnamed diagnostic."
    ENDIF
    ! END DEBUGGING

    ! Return
    ThisDiagn => NULL()
    Arr2D     => NULL()
    RC = HCO_SUCCESS

  END SUBROUTINE Diagn_Update
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
  SUBROUTINE Diagn_Get( am_I_Root, EndOfIntvOnly, DgnCont,       &
                        FLAG,      RC,            cName,         &
                        cID,       AutoFill,      InclManual, COL )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !INPUT PARAMETERS::
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root      ! Root CPU?
    LOGICAL,          INTENT(IN   )           :: EndOfIntvOnly  ! End of 
                                                                ! interval 
                                                                ! only? 
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: cName          ! container name
    INTEGER,          INTENT(IN   ), OPTIONAL :: cID            ! container ID
    INTEGER,          INTENT(IN   ), OPTIONAL :: AutoFill       ! 0=no; 1=yes; 
                                                                ! -1=either
    LOGICAL,          INTENT(IN   ), OPTIONAL :: InclManual     ! Include manual cont.? 
    INTEGER,          INTENT(IN   ), OPTIONAL :: COL            ! Collection Nr. 
!
! !OUTPUT PARAMETERS:
!

    TYPE(DiagnCont),  POINTER                 :: DgnCont        ! Return 
                                                                ! container
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: FLAG           ! Return flag
    INTEGER,          INTENT(INOUT)           :: RC             ! Return code 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER  :: MinResetFlag
    INTEGER  :: PS, AF
    LOGICAL  :: FOUND, CF, Manual 

    !======================================================================
    ! Diagn_Get begins here!
    !======================================================================

    ! Init
    FLAG   = HCO_FAIL
    RC     = HCO_SUCCESS
    CF     = .FALSE.
    Manual = .FALSE.

    ! Get collection number
    PS = 1
    IF ( PRESENT(COL) ) PS = COL
    IF ( PS > MaxCollections ) THEN
       DgnCont => NULL()
       RETURN
    ENDIF 
    IF ( .NOT. Collections(PS)%InUse ) THEN
       DgnCont => NULL()
       RETURN
    ENDIF

    ! DEBUGGING - ewl, 2/2/15
    PRINT *, "   In subroutine Diagn_Get in hco_diagn_mod.F90"
    IF ( PRESENT (cName) ) THEN
       PRINT *, "      Diag name passed as argument: " // TRIM( cName )
    ELSE
       PRINT *, "      Name not passed to subroutine."
    ENDIF
    ! END DEBUGGING

    ! Set AutoFill flag
    AF = -1
    IF ( PRESENT(AutoFill  ) ) AF     = AutoFill
    IF ( PRESENT(InclManual) ) Manual = InclManual

    ! Get minimum reset flag for current time. Set reset flag to -1 if
    ! EndOFIntvOnly flag is disabled. This will make sure that all 
    ! diagnostics are considered.
    IF ( .NOT. EndOfIntvOnly ) THEN
       IF ( Manual ) THEN
          MinResetFlag = ResetFlagManually
       ELSE
          MinResetFlag = ResetFlagEnd
       ENDIF
    ELSE
       MinResetFlag = HcoClock_GetMinResetFlag()
    ENDIF

    ! If current time stamp is not at the end of an interval - or if 
    ! there is no diagnostics container in the list with a reset flag 
    ! smaller or equal to MinResetFlag - there will be no matching 
    ! container whatsoever. Can leave right here.
    IF ( MinResetFlag > Collections(PS)%MaxResetFlag ) THEN
       DgnCont => NULL()
       RETURN
    ENDIF

    ! If container name is given, search for diagnostics with 
    ! the given name. 
    IF ( PRESENT( cName ) ) THEN

       CALL DiagnCont_Find( -1, -1, -1, -1, -1, cName, &
                            AF, FOUND, DgnCont, COL=PS )

       IF ( .NOT. FOUND ) THEN

          ! DEBUGGING - ewl, 2/2/15
          PRINT *, "      Diagnostic with name " // TRIM(cName) // " not found."
          ! END DEBUGGING
          
          DgnCont => NULL()
       ELSE

          ! DEBUGGING - ewl, 2/2/15
          PRINT *, "      Diagnostic found."
          ! END DEBUGGING

          ! Don't consider container if not at the desired
          ! time interval or if counter is zero.
          IF ( DgnCont%ResetFlag <  MinResetFlag .OR. &
               DgnCont%Counter   == 0                  ) THEN

             ! DEBUGGING - ewl, 2/2/15
             PRINT *, "      Diagnostic at wrong time interval or counter is 0."
             ! END DEBUGGING

             DgnCont => NULL()
          ENDIF
       ENDIF
       CF = .TRUE.

    ENDIF
   
    ! If container id is given, search for diagnostics with 
    ! the given container ID.
    IF ( PRESENT( cID ) ) THEN

       ! DEBUGGING - ewl, 2/2/15
       PRINT *, "      Searching for diagnostic using id ", cID
       ! END DEBUGGING

       CALL DiagnCont_Find( cID, -1, -1, -1, -1, '', &
                            AF, FOUND, DgnCont, COL=PS )
       IF ( .NOT. FOUND ) THEN

          ! DEBUGGING - ewl, 2/2/15
          PRINT *, "      Diagnostic with that id not found."
          ! END DEBUGGING

          DgnCont => NULL()
       ELSE

          ! DEBUGGING - ewl, 2/2/15
          PRINT *, "      Diagnostic with that id found."
          ! END DEBUGGING

          ! Don't consider container if not at the desired
          ! time interval or if counter is zero.
          IF ( DgnCont%ResetFlag <  MinResetFlag .OR. &
               DgnCont%Counter   == 0                  ) THEN
             DgnCont => NULL()

             ! DEBUGGING - ewl, 2/2/15
             PRINT *, "      Diagnostic at wrong time interval or counter is 0."
             ! END DEBUGGING

          ENDIF
       ENDIF
       CF = .TRUE.
    ENDIF

    ! If no container selected yet, point to next container in 
    ! list (or to head of list if DgnCont is not yet associated). 
    ! Number of updates since last output must be larger than zero!
    IF ( .NOT. CF ) THEN 

       ! DEBUGGING - ewl, 2/2/15
       PRINT *, "      No container selected yet."
       ! END DEBUGGING

       IF ( .NOT. ASSOCIATED( DgnCont ) ) THEN
          DgnCont => Collections(PS)%DiagnList

          ! DEBUGGING - ewl, 2/2/15
          PRINT *, "      DgnCont not associated so pointing to head of list."
          ! END DEBUGGING

       ELSE
          DgnCont => DgnCont%NextCont

          ! DEBUGGING - ewl, 2/2/15
          PRINT *, "      DiagCont associated so pointing to next container."
          ! END DEBUGGING

       ENDIF

       DO WHILE ( ASSOCIATED ( DgnCont ) ) 

          ! DEBUGGING - ewl, 2/2/15
          PRINT *, "         Next container is associated" 
          PRINT *, "            DgnCont%cName = ", DgnCont%Cname
          PRINT *, "            Counter = ", DgnCont%Counter
          ! END DEBUGGING
          
          IF ( DgnCont%Counter > 0 ) THEN
             
             ! DEBUGGING - ewl, 2/2/15
             PRINT *,"            Exiting do loop since counter > 0"
             ! END DEBUGGING

             EXIT
          ENDIF
          DgnCont => DgnCont%NextCont

       ENDDO
  
       ! If EndOfIntvOnly flag is enabled, make sure that the
       ! selected container is at the end of its interval.
       IF ( EndOfIntvOnly ) THEN

          ! DEBUGGING - ewl, 2/2/15
          PRINT *, "      EndOfIntvOnly is TRUE"
          ! END DEBUGGING

          ! If MinResetFlag > 0, search for first container with a
          ! ResetFlag equal or larger than MinResetFlag and where
          ! updates since last output (counter) is not zero.
          DO WHILE ( ASSOCIATED ( DgnCont ) ) 
             IF ( DgnCont%ResetFlag >= MinResetFlag .AND. &
                  DgnCont%Counter   >  0                   ) EXIT
             DgnCont => DgnCont%NextCont
          ENDDO
   
       ENDIF ! EndOfIntvOnly
    ENDIF

    ! Before returning container, make sure its data is ready for output.
    IF ( ASSOCIATED (DgnCont ) ) THEN
       CALL DiagnCont_PrepareOutput ( DgnCont, PS, RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN
       FLAG = HCO_SUCCESS

       ! Increase number of times this container has been called by
       ! Diagn_Get
       DgnCont%nnGetCalls = DgnCont%nnGetCalls + 1

       ! DEBUGGING - ewl, 2/2/15
       PRINT *, "      Diagnostic prepared for output and Diagn_Get is done."
       ! END DEBUGGING

    ENDIF

    ! DEBUGGING - ewl, 2/2/15
    PRINT *, "   Exiting Diagn_Get."
    ! END DEBUGGING

  END SUBROUTINE Diagn_Get
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

    ! Reset all internal variables to default initial values
!    nnDiagn            = 0
!    MaxResetFlag       = -1 
!    AF_LevelDefined(:) = .FALSE.
!    DiagnPrefix        = 'HEMCO_Diagn'

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
    INTEGER :: PS

    ! Initialize
    IsDefined = .FALSE.

    ! Get collection number
    PS = 1
    IF ( PRESENT(COL) ) PS = COL
    IF ( PS > MaxCollections         ) RETURN
    IF ( .NOT. Collections(PS)%InUse ) RETURN 

    IsDefined = Collections(PS)%AF_LevelDefined( Level )

  END FUNCTION Diagn_AutoFillLevelDefined
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagn_GetMaxResetFlag
!
! !DESCRIPTION: Function Diagn\_GetMaxResetFlag returns the highest reset
! flag used by any of the containers in the diagnostics list. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION Diagn_GetMaxResetFlag ( COL ) RESULT ( MaxRF )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN), OPTIONAL :: COL       ! Collection Nr.
!
! !RETURN VALUE:
!
    INTEGER                       :: MaxRF !Maximum reset flag 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER :: PS

    ! Begins here

    ! Initialize
    MaxRF = ResetFlagManually 

    ! Get collection number
    PS = 1
    IF ( PRESENT(COL) ) PS = COL
    IF ( PS > MaxCollections         ) RETURN
    IF ( .NOT. Collections(PS)%InUse ) RETURN 

    MaxRF = Collections(PS)%MaxResetFlag 
 
  END FUNCTION Diagn_GetMaxResetFlag
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_GetDiagnPrefix
!
! !DESCRIPTION: Subroutine Diagn\_GetDiagnPrefix returns the HEMCO diagnostics
! file prefix as set in the HEMCO configuration file. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_GetDiagnPrefix( Prefix, RC, COL )
!
! !INPUT ARGUMENTS:
!
    INTEGER,          INTENT(IN), OPTIONAL :: COL       ! Collection Nr.
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(OUT)          :: Prefix
    INTEGER,          INTENT(OUT)          :: RC
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER            :: PS
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'Diagn_GetDiagnPrefix (hco_diagn_mod.F90)'

    !======================================================================
    ! Diagn_GetDiagnPrefix begins here!
    !======================================================================

    ! Get collection number
    PS = 1
    IF ( PRESENT(COL) ) PS = COL
    IF ( PS > MaxCollections ) THEN
       CALL HCO_ERROR ( 'PS > MaxCollections!', RC, THISLOC=LOC )
       RETURN
    ENDIF 
    IF ( .NOT. Collections(PS)%InUse ) THEN
       WRITE(MSG,*) 'Collection not defined: ', PS
       CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Get prefix from collection
    Prefix = Collections(PS)%PREFIX

    ! Return w/ success
    RC = HCO_SUCCESS
 
  END SUBROUTINE Diagn_GetDiagnPrefix
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
    DgnCont%Scalar   =  0.0_hp
    DgnCont%LevIdx   = -1
    DgnCont%AutoFill =  1

    ! Default values for unit conversion factors 
    DgnCont%MassScal  = 1.0_dp
    DgnCont%AreaScal  = 1.0_dp
    DgnCont%AreaFlag  = 2
    DgnCont%Counter   = 0
    DgnCont%TimeAvg   = -1   
    DgnCont%AvgFlag   = -1

    ! Set last update time to -1 to start with
    DgnCont%LastUpdateID = -1

    ! By default, data is not in output format
    DgnCont%IsOutFormat = .FALSE.
    DgnCont%nnGetCalls  = 0

    ! Default container ID 
    DgnCont%cID         = 0 

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
  SUBROUTINE DiagnCont_PrepareOutput ( DgnCont, COL, RC )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !INPUT PARAMETERS::
!
    INTEGER,           INTENT(IN   ) :: COL      ! Collection number
!
! !OUTPUT PARAMETERS:
!
    TYPE(DiagnCont),   POINTER       :: DgnCont  ! diagnostics container 
!
! !INPUT/OUTPUT PARAMETERS:
!
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
    INTEGER             :: I, J, YYYY, MM, PS
    REAL(hp)            :: norm1, mult1, DPY, totscal
    CHARACTER(LEN=255)  :: MSG, LOC
    INTEGER             :: DPM(12) = (/ 31, 28, 31, 30, 31, 30, &
                                        31, 31, 30, 31, 30, 31   /)

    !======================================================================
    ! DiagnCont_PrepareOutput begins here!
    !======================================================================

    ! Init
    RC  = HCO_SUCCESS
    LOC = 'DiagnCont_PrepareOutput (hco_diagn_mod.F90) '
   
    ! Get collection number
    PS = COL
 
    !-----------------------------------------------------------------------
    ! Don't do anything for pointer data and/or if data is already in 
    ! output format
    !-----------------------------------------------------------------------
    IF ( DgnCont%IsOutFormat ) RETURN
    IF ( DgnCont%DtaIsPtr    ) RETURN

    !-----------------------------------------------------------------------
    ! Return w/ error if counter is still zero. This should not happen!
    !-----------------------------------------------------------------------
    IF ( DgnCont%Counter == 0 ) THEN
       MSG = 'Counter is zero : ' // TRIM(DgnCont%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
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
    IF ( DgnCont%AvgFlag == AvgFlagSum .OR. DgnCont%AvgFlag == AvgFlagCumsum ) THEN
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
       norm1 = REAL(DgnCont%Counter,kind=hp) * Collections(PS)%TS

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
    totscal = mult1 / norm1 * DgnCont%MassScal * DgnCont%AreaScal 

    ! For 3D:
    IF ( DgnCont%SpaceDim == 3 ) THEN
       DO J = 1, Collections(PS)%NY
       DO I = 1, Collections(PS)%NX

          ! Multiply by area if output unit is not per area 
          IF ( DgnCont%AreaFlag == 0 ) THEN
             DgnCont%Arr3D%Val(I,J,:) = DgnCont%Arr3D%Val(I,J,:)  & 
                                      * Collections(PS)%AREA_M2(I,J) 
          ENDIF

          ! Apply scale factors
          DgnCont%Arr3D%Val(I,J,:) = DgnCont%Arr3D%Val(I,J,:) & 
                                   * totscal
       ENDDO !I
       ENDDO !J

    ! For 2D:
    ELSEIF ( DgnCont%SpaceDim == 2 ) THEN
       DO J = 1, Collections(PS)%NY
       DO I = 1, Collections(PS)%NX

          ! Multiply by area if output unit is not per area 
          IF ( DgnCont%AreaFlag == 0 ) THEN
             DgnCont%Arr2D%Val(I,J) = DgnCont%Arr2D%Val(I,J) &
                                    * Collections(PS)%AREA_M2(I,J) 
          ENDIF

          ! Apply scale factors
          DgnCont%Arr2D%Val(I,J) = DgnCont%Arr2D%Val(I,J) &
                                 * totscal

       ENDDO !I
       ENDDO !J

    ! For 1D:
    ELSE
       DgnCont%Scalar = DgnCont%Scalar * totscal

    ENDIF

    ! Data is now in output format
    DgnCont%IsOutFormat = .TRUE.

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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                            :: PS
    TYPE(DiagnCont),   POINTER         :: CurrCnt => NULL() 
    LOGICAL                            :: IsMatch, Rsm
 
    !======================================================================
    ! DiagnCont_Find begins here!
    !======================================================================

    ! Initialize
    FOUND  = .FALSE.

    ! Collection position
    PS = 1
    IF ( PRESENT(COL) ) PS = COL 
    IF ( PS > MaxCollections ) RETURN
    IF ( .NOT. Collections(PS)%InUse ) RETURN

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
       CurrCnt => Collections(PS)%DiagnList 
    ENDIF

    ! Error trap
    IF ( .NOT. ASSOCIATED(CurrCnt) ) THEN
       OutCnt => NULL()
       RETURN
    ENDIF

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

    ! Cleanup
    CurrCnt => NULL()

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
  SUBROUTINE DiagnCont_Link_2D( am_I_Root, DgnCont, Tgt2D, COL, RC )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !ARGUMENTS:
!
    LOGICAL,               INTENT(IN   )         :: am_I_Root  ! Root CPU?
    INTEGER,               INTENT(IN   )         :: COL        ! Collection Nr. 
    REAL(hp),              INTENT(IN   ), TARGET :: Tgt2D(:,:) ! 2D target data 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DiagnCont),       POINTER               :: DgnCont    ! diagnostics 
                                                           !  container
    INTEGER,               INTENT(INOUT)         :: RC         ! Return code 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)   :: MSG, LOC 

    !======================================================================
    ! DiagnCont_Link_2D begins here!
    !======================================================================

    ! Init
    LOC = 'DiagnCont_Link_2D (hco_diagn_mod.F90)'

    ! Check if dimensions match. Also, containers with pointers must not
    ! be set to AutoFill
    IF ( DgnCont%AutoFill == 1 ) THEN
       MSG = 'Cannot link AutoFill container: ' // TRIM(DgnCont%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( DgnCont%SpaceDim /= 2 ) THEN
       MSG = 'Diagnostics is not 2D: ' // TRIM(DgnCont%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( SIZE(Tgt2D,1) /= Collections(COL)%NX .OR. &
         SIZE(Tgt2D,2) /= Collections(COL)%NY       ) THEN
       MSG = 'Incorrect target array size: ' // TRIM(DgnCont%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Define 2D array pointer
    CALL HCO_ArrInit( DgnCont%Arr2D, 0, 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Point to data
    DgnCont%Arr2D%Val => Tgt2D

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
  SUBROUTINE DiagnCont_Link_3D( am_I_Root, DgnCont, Tgt3D, COL, RC )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !INPUT PARAEMTERS:
!
    LOGICAL,               INTENT(IN   )         :: am_I_Root    ! Root CPU?
    INTEGER,               INTENT(IN   )         :: COL          ! Collection Nr. 
    REAL(hp),              INTENT(IN   ), TARGET :: Tgt3D(:,:,:) ! 3D target data 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DiagnCont),       POINTER               :: DgnCont      ! diagnostics 
                                                                 !  container
    INTEGER,               INTENT(INOUT)         :: RC           ! Return code 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL ARGUMENTS:
!
    CHARACTER(LEN=255)   :: MSG, LOC 

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
    IF ( DgnCont%SpaceDim /= 3 ) THEN
       MSG = 'Diagnostics is not 3D: ' // TRIM(DgnCont%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( SIZE(Tgt3D,1) /= Collections(COL)%NX .OR. &
         SIZE(Tgt3D,2) /= Collections(COL)%NY .OR. &
         SIZE(Tgt3D,3) /= Collections(COL)%NZ       ) THEN
       MSG = 'Incorrect target array size: ' // TRIM(DgnCont%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Define 3D array pointer
    CALL HCO_ArrInit( DgnCont%Arr3D, 0, 0, 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Point to data
    DgnCont%Arr3D%Val => Tgt3D

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
  SUBROUTINE Diagn_Print ( Dgn, Verbose )
!
! !INPUT ARGUMENTS:
!
    TYPE(DiagnCont), POINTER    :: Dgn
    LOGICAL,         INTENT(IN) :: Verbose
!
! !REVISION HISTORY:
!  01 Aug 2014 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ARGUMENTS:
!
    CHARACTER(LEN=255) :: MSG 
    CHARACTER(LEN= 31) :: WriteFreq 
    INTEGER            :: nx, ny, nz
    REAL(dp)           :: sm

    ! ================================================================
    ! Diagn_Print begins here
    ! ================================================================

    sm = 0.0_dp
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
    IF ( verbose ) THEN

       ! Write frequency
       SELECT CASE ( Dgn%ResetFlag )
          CASE ( ResetFlagAnnually )
             WriteFreq = 'Annually' 
          CASE ( ResetFlagMonthly  )
             WriteFreq = 'Monthly'
          CASE ( ResetFlagDaily    )
             WriteFreq = 'Daily'
          CASE ( ResetFlagHourly   )
             WriteFreq = 'Hourly'
          CASE ( ResetFlagEnd      )
             WriteFreq = 'End'
          CASE ( ResetFlagManually )
             WriteFreq = 'Manual'
       END SELECT

       ! General information
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
       WRITE(MSG,*) '   --> Write frequency    : ', TRIM(WriteFreq)
       CALL HCO_MSG(MSG)
    ENDIF

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
  SUBROUTINE DiagnCollection_Create ( am_I_Root, NX,     NY, NZ,  &
                                      TS,   AM2, PREFIX, RC, COL, & 
                                      Overwrite )
!
! !USES:
!
      USE HCO_EXTLIST_MOD,        ONLY : GetExtOpt
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
    INTEGER,            INTENT(IN), OPTIONAL :: COL        ! Collection Nr. 
    LOGICAL,            INTENT(IN), OPTIONAL :: OVERWRITE  ! OverWrite existing? 
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
    INTEGER            :: PS
    LOGICAL            :: SAME
    LOGICAL            :: OW
    CHARACTER(LEN=255) :: MSG,  MyPrefix 
    CHARACTER(LEN=255) :: LOC = 'DiagnCollection_Create (hco_diagn_mod.F90)' 

    ! ================================================================
    ! DiagnCollection_Create begins here
    ! ================================================================

    ! Set POSITION
    PS = 1
    IF ( PRESENT(COL) ) PS = COL

    ! OverWrite existing?
    OW = .FALSE.
    IF ( PRESENT(OVERWRITE) ) OW = OVERWRITE

    ! Position must not exceed max. number of collections
    IF ( PS > MaxCollections ) THEN
       WRITE(MSG,*) 'Collection position too high. Please increase ', &
                    'parameter MaxCollections in ', TRIM(LOC)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Eventually initialize the Collections vector
    IF ( .NOT. ASSOCIATED(Collections) ) THEN
       ALLOCATE( Collections(MaxCollections) )
    ENDIF 

    ! Check if collection is already defined at this position. If so,
    ! all input arguments must exactly match the existing parameter
    IF ( Collections(PS)%InUse .AND. .NOT. OW ) THEN
       SAME = .TRUE.
       IF (      Collections(PS)%NX      /= NX           ) SAME = .FALSE. 
       IF (      Collections(PS)%NY      /= NY           ) SAME = .FALSE. 
       IF (      Collections(PS)%NZ      /= NZ           ) SAME = .FALSE. 
       IF (      Collections(PS)%TS      /= TS           ) SAME = .FALSE. 
       IF (  ANY(Collections(PS)%AREA_M2 /= AM2 )        ) SAME = .FALSE.
       IF ( TRIM(Collections(PS)%PREFIX) /= TRIM(PREFIX) ) SAME = .FALSE.

       IF ( .NOT. SAME ) THEN
          WRITE(MSG,*) 'Collection at position ', PS, ' already exists'
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ELSE
          RC = HCO_SUCCESS
          RETURN
       ENDIF
    ENDIF

    ! Pass arguments
    Collections(PS)%NX      =  NX 
    Collections(PS)%NY      =  NY 
    Collections(PS)%NZ      =  NZ 
    Collections(PS)%TS      =  TS 
    Collections(PS)%AREA_M2 => AM2 

    ! For emissions diagnostics collections and if the prefix is empty,
    ! try to get prefix from the HEMCO configuration file. 
    IF ( PS == 1 .AND. TRIM(PREFIX) == '' ) THEN
       CALL GetExtOpt ( 0, 'DiagnPrefix', OptValChar=MyPrefix, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       Collections(PS)%PREFIX = TRIM(MyPrefix)
    ELSE
       Collections(PS)%PREFIX = TRIM(PREFIX)
    ENDIF

    ! This collection is now in use
    Collections(PS)%InUse = .TRUE.

    ! verbose
    IF ( am_I_Root .AND. HCO_VERBOSE_CHECK() ) THEN
       MSG = 'Created/updated diagnostics collection: '
       CALL HCO_MSG(MSG)
       WRITE(MSG,'(a13,i2)') ' - POSITION: ', PS
       CALL HCO_MSG(MSG)
       WRITE(MSG,'(a13,a)' ) ' - PREFIX  : ', TRIM(Collections(PS)%PREFIX)
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
  SUBROUTINE DiagnCollection_Cleanup ( COL )
!
! !INPUT ARGUMENTS:
!
    INTEGER,            INTENT(IN), OPTIONAL :: COL  ! Collection number 
!
! !REVISION HISTORY:
!  08 Jan 2015 - C. Keller - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ARGUMENTS:
!
    INTEGER            :: PS

    ! ================================================================
    ! DiagnCollection_Cleanup begins here
    ! ================================================================

    ! Set POSITION
    PS = 1
    IF ( PRESENT(COL) ) PS = COL
    IF ( PS > MaxCollections ) RETURN

    ! Cleanup if in use
    IF ( Collections(PS)%InUse ) THEN
       CALL DiagnList_Cleanup( Collections(PS)%DiagnList )
       Collections(PS)%nnDiagn = 0
       Collections(PS)%AREA_M2 => NULL()
       Collections(PS)%InUse   = .FALSE.
    ENDIF 

  END SUBROUTINE DiagnCollection_Cleanup
!EOC
END MODULE HCO_Diagn_Mod
