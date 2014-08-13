!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_diagn_mod.F90
!
! !DESCRIPTION: Module HCO\_Diagn\_mod contains routines and 
! variables to handle the HEMCO diagnostics. The HEMCO diagnostics
! consist of a collection of diagnostics container organized through
! list DiagnList. Each diagnostics container contains information 
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
! of [kg/m2] and only converted to desired output unit when returning 
! the data. The container variable IsOutFormat denotes whether data
! is in output units or internal units. Variable nnGetCalls counts the
! number of times a diagnostics is called through Diagn\_Get without
! updating its content. This is useful if you want to make sure that
! data is only written once per time step.
!\\
!\\
! There are two types of diagnostics: automatic (`AutoFill`) and 
! manual diagnostics. AutoFill diagnostics become automatically filled 
! during HEMCO execution. AutoFill diagnostics can be at species level
! (level 1), ExtNr level (level 2), emission category level (level 3),
! or hierarchy level (level 4). Level 1 diagnostics write out the
! collected emissions of the specified species, level 2 diagnostics
! write out emissions for the given ExtNr only (ignoring emissions from
! all other ExtNr's), etc.
!\\
!\\
! Manual diagnostics can represent any content. They never become filled
! automatically and all update calls (Diagn\_Update) have to be set
! manually.
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
  PUBLIC  :: Diagn_Cleanup 
  PUBLIC  :: Diagn_Create
  PUBLIC  :: Diagn_Update 
  PUBLIC  :: Diagn_Get
  PUBLIC  :: Diagn_AutoFillLevelDefined
  PUBLIC  :: Diagn_GetMaxResetFlag
  PUBLIC  :: Diagn_GetDiagnPrefix
  PUBLIC  :: Diagn_Print
  PUBLIC  :: DiagnCont_Find
!
! !PRIVATE MEMBER FUNCTIONS:
!
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
     INTEGER                     :: SpaceDim       ! 1, 2, or 3
     REAL(hp)                    :: Scalar         ! 1D scalar 
     TYPE(Arr2D_HP),     POINTER :: Arr2D          ! 2D array
     TYPE(Arr3D_HP),     POINTER :: Arr3D          ! 3D array
     LOGICAL                     :: DtaIsPtr       ! Data is just a pointer?
     INTEGER                     :: LevIdx         ! Level index to be used 
     CHARACTER(LEN= 31)          :: OutUnit        ! Output unit 
     LOGICAL                     :: PerArea        ! Is output unit per area?
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
  ! Other diagnostic variables
  !------------------------------------------------------------------------

  ! Diagnostics list 
  TYPE(DiagnCont),       POINTER :: DiagnList => NULL()

  ! Number of diagnostics container
  INTEGER                        :: nnDiagn = 0

  ! Highest reset flag by any of the containers 
  INTEGER                        :: MaxResetFlag = ResetFlagManually

  ! For AutoFill level flags 
  LOGICAL                        :: AF_LevelDefined(4) = .FALSE.

  ! File prefix of HEMCO diagnostics. Will be expanded by date
  ! & time (YYYYMMDDhm).
  CHARACTER(LEN=255)             :: DiagnPrefix = 'HEMCO_Diagn'
!
! !DEFINED PARAMETERS:
!
  ! Parameter for averaging and summing non-standard data
  ! AvgFlagMean    : calculates the arithmetic mean
  ! AvgFlagSum     : calculates the sum, resets after every writeout
  ! AvgFlagCumSum  : calculates the cumulative sum, never resets.
  INTEGER, PARAMETER             :: AvgFlagMean   = 1
  INTEGER, PARAMETER             :: AvgFlagSum    = 2
  INTEGER, PARAMETER             :: AvgFlagCumsum = 2

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
!  12 Sep 2013 - C. Keller   - Initial version 
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
             CALL Diagn_Update( am_I_Root,  HcoState, ExtNr=-1, &
                                Cat=-1,     Hier=-1,  HcoID=I,  &
                                AutoFill=1, Array3D=Arr3D, RC=RC ) 
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
!\item HcoState: HEMCO state object. 
!\item cName: distinct diagnostics (container) name. 
!\item ExtNr: emissions extension number. 
!\item Cat: emissions category. 
!\item Hier: emissions  hierarchy.
!\item HcoID: HEMCO species ID of diagnostics species.
!\item SpaceDim: spatial dimension: 1 (scalar), 2 (lon-lat), 
!      or 3 (lon-lat-lev). 
!\item OutUnit: output unit. Emissions will be converted to this unit.
!      Conversion factors will be determined using the HEMCO unit
!      module (see HCO\_UNITS\_Mod.F90).
!\item WriteFreq: output frequency. Can be one of 'Hourly', 'Daily',
!      'Monthly', 'Annualy', 'End', 'Manual'.
!      Manual diagnostics are expected to be manually received and 
!      written out. These diagnostics may or may not be written out 
!      at the end of the simulation run, depending on the corresponding
!      attribute set in Diagn\_Get.
!\item OutOper: output operation for non-standard units. If this 
!      argument is used, the specified operation is performed and all
!      unit specifications are ignored. Can be one of 'Mean', 'Sum', 
!      or 'Cumsum'.
!\item AutoFill: containers with an AutoFill flag of 1 will be auto-
!      matically updated by the HEMCO standard diagnostics calls
!      (e.g. in hco\_calc\_mod.F90). If set to 0, the diagnostics 
!      updates have to be set manually.
!\item Trgt2D: 2D target array. If specified, the diagnostics array
!      will point to this data. This disables all time averaging, 
!      unit conversions, etc., and the data will be written to disk
!      as is.
!\item Trgt3D: as Trgt2D, but for 3D data. 
!\item cID: assigned container ID. Useful for later reference to this
!      diagnostics container.
!\item RC: HEMCO return code.
!\end{itemize} 
!
! !INTERFACE:
!
  SUBROUTINE Diagn_Create( am_I_Root, HcoState,  cName,      &
                           ExtNr,     Cat,       Hier,       &
                           HcoID,     SpaceDim,  OutUnit,    &  
                           WriteFreq, OutOper,   LevIdx,     &
                           AutoFill,  Trgt2D,    Trgt3D,     &
                           cID,       RC                      )
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
    TYPE(HCO_State),  POINTER                 :: HcoState      ! HEMCO state object
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
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: OutOper       ! Output 
                                                               !  operation 
    INTEGER,          INTENT(IN   ), OPTIONAL :: LevIdx        ! Level index 
                                                               !  to use
    INTEGER,          INTENT(IN   ), OPTIONAL :: AutoFill      ! 1=auto fill
                                                               ! 0=don't 
    REAL(hp),         INTENT(IN   ), OPTIONAL :: Trgt2D(:,:)   ! 2D target data
    REAL(hp),         INTENT(IN   ), OPTIONAL :: Trgt3D(:,:,:) ! 3D target data
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
!
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
    INTEGER                  :: ThiscID
    REAL(hp)                 :: Scal
    LOGICAL                  :: ForceMean, FOUND

    !======================================================================
    ! Diagn_Create begins here!
    !======================================================================

    ! Init
    LOC = 'Diagn_Create (hco_diagn_mod.F90)'

    !----------------------------------------------------------------------
    ! Initalize diagnostics container. This will automatically add the
    ! container to the diagnostics list.
    !----------------------------------------------------------------------
    CALL DiagnCont_Init( ThisDiagn, ThiscID )

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
       CALL DiagnCont_Link_2D( am_I_Root, HcoState, ThisDiagn, Trgt2D, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF
    IF ( PRESENT(Trgt3D) ) THEN
       CALL DiagnCont_Link_3D( am_I_Root, HcoState, ThisDiagn, Trgt3D, RC )
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
    MaxResetFlag = MAX( MaxResetFlag, ThisDiagn%ResetFlag )

    ! Update module variable AF_LevelDefined. For all AutoFill diagnostics,
    ! we store whether or not there is (at least one) diagnostics container
    ! defined at species level, ExtNr level, etc. 
    IF ( ThisDiagn%AutoFill == 1 ) THEN

       ! At species level: no ExtNr defined
       IF ( ThisDiagn%ExtNr < 0 ) THEN
          AF_LevelDefined(1) = .TRUE.
     
       ! At ExtNr level: no category defined
       ELSEIF ( ThisDiagn%Cat < 0 ) THEN
          AF_LevelDefined(2) = .TRUE.

       ! At category level: no hierarchy defined
       ELSEIF ( ThisDiagn%Hier < 0 ) THEN
          AF_LevelDefined(3) = .TRUE.

       ! At hierarchy level: all defined
       ELSE
          AF_LevelDefined(4) = .TRUE.
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
          ELSE
             MSG = 'Illegal output operator: ' // TRIM(OutOper)
             MSG = TRIM(MSG) // '. Allowed are `Mean`, `Sum`, `Cumsum`.'
             MSG = TRIM(MSG) // ' (' // TRIM(cName) // ')'
             CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          ENDIF
   
       ! If OutOper is not set, determine scale factors from output unit:
       ELSE
   
          !----------------------------------------------------------------
          ! Scale factor for mass. This determines the scale factor from 
          ! HEMCO mass unit (kg) to the desired output unit. 
          ! HCO_UNIT_MassCal returns the mass scale factor from OutUnit to
          ! HEMCO unit, hence need to invert this value!
          !----------------------------------------------------------------
          Scal = HCO_UNIT_GetMassScal(                          &
                    unt         = OutUnit,                      &
                    MW_IN       = HcoState%Spc(HcoID)%MW_g,     &
                    MW_OUT      = HcoState%Spc(HcoID)%EmMW_g,   &
                    MOLEC_RATIO = HcoState%Spc(HcoID)%MolecRatio )
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
          ! Returns a negative scale factor if no valid area unit is found 
          ! in OutUnit. In this case, set PerArea to False so that the final 
          ! diagnostics will be multiplied by the surface grid box area!
          !----------------------------------------------------------------
          Scal =  HCO_UNIT_GetAreaScal( OutUnit )
          IF ( Scal < 0.0_hp ) THEN
             ThisDiagn%PerArea  = .FALSE.
          ELSE
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
          Scal = HCO_UNIT_GetTimeScal( OutUnit, 1, 2001 )
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
    !-----------------------------------------------------------------------
    CALL DiagnCont_Find( -1, ThisDiagn%ExtNr, ThisDiagn%Cat,         &
                             ThisDiagn%Hier,  ThisDiagn%HcoID, '', 1, &
                             FOUND, TmpDiagn )
    IF ( FOUND .AND. TmpDiagn%AutoFill==1 .AND. ThisDiagn%AutoFill==1 ) THEN
       MSG = 'These two diagnostics seem to be the same:'
       CALL HCO_MSG(MSG)
       CALL Diagn_Print( ThisDiagn, .TRUE. ) 
       CALL Diagn_Print( TmpDiagn,  .TRUE. )
       MSG = 'Cannot add diagnostics - duplicate entry!'
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Add to diagnostics list. Insert at the beginning of the list.
    !-----------------------------------------------------------------------
    IF ( nnDiagn > 1 ) THEN
       ThisDiagn%NextCont => DiagnList
    ENDIF
    DiagnList => ThisDiagn

    ! Verbose mode
    IF ( am_I_Root .AND. HCO_VERBOSE_CHECK() ) THEN
       MSG = 'Successfully added diagnostics: ' 
       CALL HCO_MSG ( MSG )
       CALL Diagn_Print( ThisDiagn, .TRUE. )
    ENDIF

    ! Return
    cID = ThiscID
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
  SUBROUTINE Diagn_Update( am_I_Root, HcoState, cID,     cName,   &
                           ExtNr,     Cat,      Hier,    HcoID,   &
                           AutoFill,  Scalar,   Array2D, Array3D, &
                           PosOnly,   RC ) 
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
    USE HCO_Arr_Mod,   ONLY : HCO_ArrAssert
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root      ! Root CPU?
    TYPE(HCO_State),  POINTER                 :: HcoState       ! HEMCO state 
                                                                !  object
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
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC             ! Return code 
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(DiagnCont), POINTER :: ThisDiagn => NULL()
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
    LOGICAL                  :: Found, OnlyPos, VertSum, IsAssoc

    !======================================================================
    ! Diagn_Update begins here!
    !======================================================================

    ! Init
    LOC = 'Diagn_Update (hco_diagn_mod.F90)'
    RC  = HCO_SUCCESS

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

    !----------------------------------------------------------------------
    ! Find the wanted container. Return if container does not exist.
    !----------------------------------------------------------------------
    CALL DiagnCont_Find( DgncID,   DgnExtNr, DgnCat,   DgnHier, &
                         DgnHcoID, DgnName,  AutoFlag, Found,   &
                         ThisDiagn                               ) 
    IF ( .NOT. Found ) THEN
       RETURN
    ENDIF

    ! If container holds just a pointer to external data, don't do
    ! anything!
    IF ( ThisDiagn%DtaIsPtr ) THEN
       MSG = 'You try to update a container that holds a '  // &
            'pointer to data - this should never happen! ' // &
            TRIM(ThisDiagn%cName)
       CALL HCO_WARNING( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Get current minimum reset flag as well as the update time ID.
    MinResetFlag = HcoClock_GetMinResetFlag()
    CALL HcoClock_Get( nSteps = ThisUpdateID, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

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
       Fact = HcoState%TS_EMIS
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
       CALL HCO_ArrAssert( ThisDiagn%Arr3D, HcoState%NX, &
                           HcoState%NY,     HcoState%NZ, RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
         
       ! Pass array to diagnostics: reset to zero if counter 
       ! is zero, add to it otherwise.
       ! Never reset containers with cumulative sums!
       IF ( ThisDiagn%Counter == 0 .AND. &
            ThisDiagn%AvgFlag /= AvgFlagCumsum ) ThisDiagn%Arr3D%Val = 0.0_hp

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
       CALL HCO_ArrAssert( ThisDiagn%Arr2D, HcoState%NX, &
                           HcoState%NY,     RC            ) 
       IF ( RC /= HCO_SUCCESS ) RETURN 
         
       ! Pass array to diagnostics: ignore existing data if counter 
       ! is zero, add to it otherwise.
       ! Never reset containers with cumulative sums!
       IF ( ThisDiagn%Counter == 0 .AND. &
            ThisDiagn%AvgFlag /= AvgFlagCumsum ) ThisDiagn%Arr2D%Val = 0.0_hp

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
                DO J=1,HcoState%NY
                DO I=1,HcoState%NX
                   TMP = 0.0_hp
                   DO L=1,HcoState%NZ
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
                DO J=1,HcoState%NY
                DO I=1,HcoState%NX
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
       IF ( OnlyPos ) THEN
          IF ( Scalar >= 0.0_hp ) & 
             ThisDiagn%Scalar = ThisDiagn%Scalar + ( Scalar * Fact )
       ELSE
          ThisDiagn%Scalar = ThisDiagn%Scalar + ( Scalar * Fact )  
       ENDIF

    ELSE
       MSG = 'Invalid space dimension: ' // TRIM(ThisDiagn%cName)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    !----------------------------------------------------------------------
    ! Update counter ==> Do only if last update time is not equal to 
    ! current one! This allows the same diagnostics to be updated
    ! multiple time on the same time step without increasing the
    ! time step counter.
    !----------------------------------------------------------------------
    IF ( ThisDiagn%LastUpdateID /= ThisUpdateID ) THEN
       ThisDiagn%Counter      = ThisDiagn%Counter + 1
       ThisDiagn%LastUpdateID = ThisUpdateID
    ENDIF

    !----------------------------------------------------------------------
    ! Data is not in output format and hasn't been called yet by Diagn_Get.
    !----------------------------------------------------------------------
    ThisDiagn%IsOutFormat = .FALSE.
    ThisDiagn%nnGetCalls  = 0

    ! Verbose mode 
    IF ( am_I_Root .AND. HCO_VERBOSE_CHECK() ) THEN
       MSG = 'Successfully updated diagnostics: ' // TRIM(ThisDiagn%cName)
       CALL HCO_MSG ( MSG )
    ENDIF

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
  SUBROUTINE Diagn_Get( am_I_Root, HcoState, EndOfIntvOnly, &
                        DgnCont,   FLAG,     RC,     cName, &
                        cID,       AutoFill, InclManual      )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !INPUT PARAMETERS::
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root      ! Root CPU?
    TYPE(HCO_State),  POINTER                 :: HcoState       ! HEMCO state
                                                                ! object
    LOGICAL,          INTENT(IN   )           :: EndOfIntvOnly  ! End of 
                                                                ! interval 
                                                                ! only? 
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL :: cName          ! container name
    INTEGER,          INTENT(IN   ), OPTIONAL :: cID            ! container ID
    INTEGER,          INTENT(IN   ), OPTIONAL :: AutoFill       ! 0=no; 1=yes; 
                                                                ! -1=either
    LOGICAL,          INTENT(IN   ), OPTIONAL :: InclManual     ! Include manual cont.? 
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
    INTEGER  :: AF
    LOGICAL  :: FOUND, CF, Manual 

    !======================================================================
    ! Diagn_Get begins here!
    !======================================================================

    ! Init
    FLAG   = HCO_FAIL
    RC     = HCO_SUCCESS
    CF     = .FALSE.
    Manual = .FALSE.

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
    IF ( MinResetFlag > MaxResetFlag ) THEN
       DgnCont => NULL()
       RETURN
    ENDIF

    ! If container name is given, search for diagnostics with 
    ! the given name. 
    IF ( PRESENT( cName ) ) THEN
       CALL DiagnCont_Find( -1, -1, -1, -1, -1, cName, AF, FOUND, DgnCont)

       IF ( .NOT. FOUND ) THEN
          DgnCont => NULL()
       ELSE
          ! Don't consider container if not at the desired
          ! time interval or if counter is zero.
          IF ( DgnCont%ResetFlag <  MinResetFlag .OR. &
               DgnCont%Counter   == 0                  ) THEN
             DgnCont => NULL()
          ENDIF
       ENDIF
       CF = .TRUE.
    ENDIF
   
    ! If container id is given, search for diagnostics with 
    ! the given container ID.
    IF ( PRESENT( cID ) ) THEN
       CALL DiagnCont_Find( cID, -1, -1, -1, -1, '', AF, FOUND, DgnCont)
       IF ( .NOT. FOUND ) THEN
          DgnCont => NULL()
       ELSE
          ! Don't consider container if not at the desired
          ! time interval or if counter is zero.
          IF ( DgnCont%ResetFlag <  MinResetFlag .OR. &
               DgnCont%Counter   == 0                  ) THEN
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
          DgnCont => DiagnList
       ELSE
          DgnCont => DgnCont%NextCont
       ENDIF
       DO WHILE ( ASSOCIATED ( DgnCont ) ) 
          IF ( DgnCont%Counter > 0 ) EXIT 
          DgnCont => DgnCont%NextCont
       ENDDO
  
       ! If EndOfIntvOnly flag is enabled, make sure that the
       ! selected container is at the end of its interval.
       IF ( EndOfIntvOnly ) THEN
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
       CALL DiagnCont_PrepareOutput ( HcoState, DgnCont, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       FLAG = HCO_SUCCESS

       ! Increase number of times this container has been called by
       ! Diagn_Get
       DgnCont%nnGetCalls = DgnCont%nnGetCalls + 1
    ENDIF

  END SUBROUTINE Diagn_Get
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_Cleanup
!
! !DESCRIPTION: Subroutine Diagn\_Cleanup cleans up all the diagnostics
! containers.
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_Cleanup 
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
    nnDiagn            = 0
    MaxResetFlag       = -1 
    AF_LevelDefined(:) = .FALSE.
    DiagnPrefix        = 'HEMCO_Diagn'

  END SUBROUTINE Diagn_Cleanup
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
  FUNCTION Diagn_AutoFillLevelDefined( Level ) RESULT ( IsDefined )
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN) :: Level     ! Level of interest
!
! !RETURN VALUE:
! 
    LOGICAL             :: IsDefined ! Return argument 
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC

    IsDefined = AF_LevelDefined( Level )
 
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
  FUNCTION Diagn_GetMaxResetFlag RESULT ( MaxRF )
!
! !RETURN VALUE:
!
    INTEGER :: MaxRF !Maximum reset flag 
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC

    MaxRF = MaxResetFlag 
 
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
  SUBROUTINE Diagn_GetDiagnPrefix( Prefix, RC )
!
! !USES:
!
      USE HCO_EXTLIST_MOD,        ONLY : GetExtOpt
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(OUT)     :: Prefix
    INTEGER,          INTENT(OUT)     :: RC
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    LOGICAL, SAVE      :: FIRST = .TRUE.
    LOGICAL            :: FOUND
    CHARACTER(LEN=255) :: MyPrefix

    !======================================================================
    ! Diagn_GetDiagnPrefix begins here!
    !======================================================================

    ! On first call, try to get DiagnPrefix from settings. If setting
    ! not found, keep default value.
    IF ( FIRST ) THEN
       CALL GetExtOpt ( 0, 'DiagnPrefix', OptValChar=MyPrefix, &
                       Found=FOUND, RC=RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( FOUND ) DiagnPrefix = MyPrefix

       FIRST = .FALSE.
    ENDIF

      Prefix = DiagnPrefix 

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
  SUBROUTINE DiagnCont_Init( OutCont, cID )
!
! !OUTPUT PARAMETERS:
!
    TYPE(DiagnCont), POINTER      :: OutCont  ! Created container
    INTEGER,         INTENT(OUT)  :: cID      ! Assigned container ID
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
    DgnCont%PerArea   = .TRUE.
    DgnCont%Counter   = 0
    DgnCont%TimeAvg   = -1   
    DgnCont%AvgFlag   = -1

    ! Set last update time to -1 to start with
    DgnCont%LastUpdateID = -1

    ! By default, data is not in output format
    DgnCont%IsOutFormat = .FALSE.
    DgnCont%nnGetCalls  = 0

    ! Assign container ID.
    ! Set default target ID to cont. ID.
    nnDiagn              = nnDiagn + 1
    DgnCont%cID          = nnDiagn

    ! Pass to output container
    OutCont => DgnCont
    cID     =  DgnCont%cID

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
  SUBROUTINE DiagnCont_PrepareOutput ( HcoState, DgnCont, RC )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !INPUT PARAMETERS::
!
    TYPE(HCO_State),   POINTER       :: HcoState ! HEMCO state object
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
!  04 Dec 2012 - C. Keller: Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I, J, YYYY, MM
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
       norm1 = REAL(DgnCont%Counter,kind=hp) * HcoState%TS_EMIS

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
       DO J = 1, HcoState%NY
       DO I = 1, HcoState%NX

          ! Multiply by area if output unit is not per area 
          IF ( .NOT. DgnCont%PerArea ) THEN
             DgnCont%Arr3D%Val(I,J,:) = DgnCont%Arr3D%Val(I,J,:)  & 
                                      * HcoState%Grid%AREA_M2(I,J) 
          ENDIF

          ! Apply scale factors
          DgnCont%Arr3D%Val(I,J,:) = DgnCont%Arr3D%Val(I,J,:) & 
                                   * totscal
       ENDDO !I
       ENDDO !J

    ! For 2D:
    ELSEIF ( DgnCont%SpaceDim == 2 ) THEN
       DO J = 1, HcoState%NY
       DO I = 1, HcoState%NX

          ! Multiply by area if output unit is not per area 
          IF ( .NOT. DgnCont%PerArea ) THEN
             DgnCont%Arr2D%Val(I,J) = DgnCont%Arr2D%Val(I,J) &
                                    * HcoState%Grid%AREA_M2(I,J) 
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
! !INTERFACE:
!
  SUBROUTINE DiagnCont_Find ( cID,   ExtNr, Cat,      Hier,        &
                              HcoID, cName, AutoFill, FOUND, OutCnt )
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
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,           INTENT(OUT)  :: FOUND    ! container found?
    TYPE(DiagnCont),   POINTER      :: OutCnt   ! matched container 
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DiagnCont),   POINTER         :: CurrCnt => NULL() 
    LOGICAL                            :: IsMatch
 
    !======================================================================
    ! DiagnCont_Find begins here!
    !======================================================================

    ! Initialize
    FOUND  = .FALSE.
    OutCnt => NULL()

    ! Error trap
    IF ( .NOT. ASSOCIATED(DiagnList) ) RETURN

    ! Make CurrCnt point to first element of the diagnostics list 
    CurrCnt => DiagnList 

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
! unit conversion, etc., i.e. the data will returned as is.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCont_Link_2D( am_I_Root, HcoState, DgnCont, Tgt2D, RC )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !ARGUMENTS:
!
    LOGICAL,           INTENT(IN   )         :: am_I_Root  ! Root CPU?
    TYPE(HCO_State),   POINTER               :: HcoState   ! HEMCO state object
    REAL(hp),          INTENT(IN   ), TARGET :: Tgt2D(:,:) ! 2D target data 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DiagnCont),   POINTER               :: DgnCont    ! diagnostics 
                                                           !  container
    INTEGER,           INTENT(INOUT)         :: RC         ! Return code 
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
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
    IF ( SIZE(Tgt2D,1) /= HcoState%NX .OR. SIZE(Tgt2D,2) /= HcoState%NY ) THEN
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
! unit conversion, etc., i.e. the data will returned as is.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCont_Link_3D( am_I_Root, HcoState, DgnCont, Tgt3D, RC )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !INPUT PARAEMTERS:
!
    LOGICAL,           INTENT(IN   )         :: am_I_Root    ! Root CPU?
    TYPE(HCO_State),   POINTER               :: HcoState     ! HEMCO state object
    REAL(hp),          INTENT(IN   ), TARGET :: Tgt3D(:,:,:) ! 3D target data 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DiagnCont),   POINTER               :: DgnCont      ! diagnostics 
                                                             !  container
    INTEGER,           INTENT(INOUT)         :: RC           ! Return code 
!
! !REVISION HISTORY:
!  04 Dec 2012 - C. Keller: Initialization
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
    IF ( SIZE(Tgt3D,1) /= HcoState%NX .OR. &
         SIZE(Tgt3D,2) /= HcoState%NY .OR. &
         SIZE(Tgt3D,3) /= HcoState%NZ       ) THEN
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
END MODULE HCO_Diagn_Mod
