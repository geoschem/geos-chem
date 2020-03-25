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
! in the HEMCO configuration file through argument 'DiagnFreq'.
! This can be a character indicating the output frequency (valid
! entries are 'Always', 'Hourly', 'Daily', 'Monthly', 'Annually',
! 'Manual', and 'End') or by two integer strings of format
! '00000000 000000' denoting the year-month-day and hour-minute-
! second output interval, respectively. For example, setting
! DiagnFreq to '00000001 000000' would be equivalent to setting
! it to 'Daily'. A value of '00000000 030000' indicates that the
! diagnostics shall be written out every 3 hours.
!\\
!\\
! The restart collection always gets an output frequency of 'End',
! but writing its content to disk can be forced at any given time
! using routine HcoDiagn\_Write (see below). The manual diagnostics
! has an output frequency of 'Manual', which means that its content
! is never written to disk. Instead, its fields need to be fetched
! explicitly from other routines via routine Diagn\_Get.
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
  USE HCO_Types_Mod
  USE HCO_Arr_Mod
  USE HCO_Clock_Mod
  USE HCO_State_Mod, ONLY : HCO_State

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
  PUBLIC  :: DiagnCollection_Set
  PUBLIC  :: DiagnCollection_GetDefaultDelta
  PUBLIC  :: DiagnCollection_IsTimeToWrite
  PUBLIC  :: DiagnCollection_LastTimesSet
  PUBLIC  :: DiagnFileOpen
  PUBLIC  :: DiagnFileGetNext
  PUBLIC  :: DiagnFileClose
  PUBLIC  :: DiagnBundle_Cleanup
  PUBLIC  :: DiagnBundle_Init
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
  PRIVATE :: Diagn_UpdateSp0d
  PRIVATE :: Diagn_UpdateSp2d
  PRIVATE :: Diagn_UpdateSp3d
  PRIVATE :: Diagn_UpdateDp0d
  PRIVATE :: Diagn_UpdateDp2d
  PRIVATE :: Diagn_UpdateDp3d

  INTERFACE Diagn_Update
     MODULE PROCEDURE Diagn_UpdateSp0d
     MODULE PROCEDURE Diagn_UpdateSp2d
     MODULE PROCEDURE Diagn_UpdateSp3d
     MODULE PROCEDURE Diagn_UpdateDp0d
     MODULE PROCEDURE Diagn_UpdateDp2d
     MODULE PROCEDURE Diagn_UpdateDp3d
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
!  06 Nov 2015 - C. Keller   - Added argument OutTimeStamp to collection to
!                              control the file output time stamp (beginning,
!                              middle, end of diagnostics interval).
!  25 Jan 2016 - R. Yantosca - Added bug fixes for pgfortran compiler
!  19 Sep 2016 - R. Yantosca - Add extra overloaded functions to the
!                              Diagn_Update interface to avoid Gfortran errors
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
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

  ! Parameter for output time stamp. This is the time stamp that will be used
  ! on the output file. End means the simulation date at output time is used,
  ! 'Mid' uses the midpoint of the diagnostics windows, 'Start' uses the
  ! beginning of the window.
  INTEGER, PARAMETER, PUBLIC     :: HcoDiagnStart   = 1
  INTEGER, PARAMETER, PUBLIC     :: HcoDiagnMid     = 2
  INTEGER, PARAMETER, PUBLIC     :: HcoDiagnEnd     = 3

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
  SUBROUTINE HcoDiagn_AutoUpdate( HcoState, RC )
!
! !USES:
!
    USE HCO_STATE_MOD, ONLY : HCO_GetHcoID
    USE HCO_STATE_MOD, ONLY : HCO_State
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
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)        :: MSG, LOC
    INTEGER                   :: I, tmpID
    REAL(hp), POINTER         :: Arr3D(:,:,:)
    REAL(hp), POINTER         :: Arr2D(:,:)

    !=================================================================
    ! HCODIAGN_AUTOUPDATE begins here!
    !=================================================================

    ! Init
    LOC   = 'HCODIAGN_AUTOUPDATE (hco_diagn_mod.F90)'
    RC    =  HCO_SUCCESS
    Arr3D => NULL()
    Arr2D => NULL()

    ! ================================================================
    ! AutoFill diagnostics: only write diagnostics at species level
    ! (level 1). Higher level diagnostics have been written in the
    ! respective subroutines (hco_calc & extension modules).
    ! ================================================================
    DO I = 1, HcoState%nSpc
       IF ( ASSOCIATED(HcoState%Spc(I)%Emis) ) THEN
          IF ( ASSOCIATED(HcoState%Spc(I)%Emis%Val) ) THEN
             Arr3D => HcoState%Spc(I)%Emis%Val
             CALL Diagn_Update( HcoState,               &
                                ExtNr      = -1,        &
                                Cat        = -1,        &
                                Hier       = -1,        &
                                HcoID      = I,         &
                                AutoFill   = 1,         &
                                Array3D    = Arr3D,     &
                                COL        = -1,        &
                                RC         = RC          )
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
  SUBROUTINE HcoDiagn_Init( HcoState, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,   ONLY : HCO_GetHcoID
    USE HCO_STATE_MOD,   ONLY : HCO_State
    USE HCO_ExtList_Mod, ONLY : GetExtOpt
    USE HCO_ExtList_Mod, ONLY : CoreNr
    USE CHARPAK_MOD,     ONLY : TRANLC
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
!  06 Nov 2015 - C. Keller   - Added OutTimeStamp.
!  01 Nov 2017 - E. Lundgren - Change default OutTimeStamp from end to start
!                              for diagnostics collection
!  29 Dec 2017 - C. Keller   - Added datetime tokens to file prefixes.
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: CollectionID
    INTEGER             :: deltaYMD, deltaHMS
    INTEGER             :: OutTimeStamp
    LOGICAL             :: FOUND
    CHARACTER(LEN=255)  :: MSG, LOC, DiagnPrefix, OutTimeStampChar

    !=================================================================
    ! HCODIAGN_INIT begins here!
    !=================================================================

    ! Init
    LOC = 'HCODIAGN_INIT (hco_diagn_mod.F90)'

    ! Initialize diagnostics bundle
    CALL DiagnBundle_Init ( HcoState%Diagn )

    ! ------------------------------------------------------------------
    ! Default diagnostics
    ! ------------------------------------------------------------------
    CALL DiagnCollection_GetDefaultDelta ( HcoState, deltaYMD, deltaHMS, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Try to get prefix from configuration file
    CALL GetExtOpt ( HcoState%Config, CoreNr, 'DiagnPrefix', &
                     OptValChar=DiagnPrefix, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
#if defined( MODEL_GEOS )
       DiagnPrefix = 'HEMCO_Diagnostics.$YYYY$MM$DD$HH$MN.nc'
#else
       DiagnPrefix = 'HEMCO_diagnostics'
#endif
    ENDIF

    ! Output time stamp location
    CALL GetExtOpt ( HcoState%Config, CoreNr, 'DiagnTimeStamp', &
                     OptValChar=OutTimeStampChar, FOUND=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( .NOT. FOUND ) THEN
       OutTimeStamp = HcoDiagnStart
    ELSE
       CALL TRANLC( OutTimeStampChar )
       IF (     TRIM(OutTimeStampChar) == 'start' ) THEN
          OutTimeStamp = HcoDiagnStart

       ELSEIF ( TRIM(OutTimeStampChar) == 'mid'   ) THEN
          OutTimeStamp = HcoDiagnMid

       ELSEIF ( TRIM(OutTimeStampChar) == 'end'   ) THEN
          OutTimeStamp = HcoDiagnEnd

       ELSE
          WRITE(MSG,*) 'Unrecognized output time stamp location: ', &
             TRIM(OutTimeStampChar), ' - will use default (start)'
          CALL HCO_WARNING(HcoState%Config%Err,MSG,RC,THISLOC=LOC,WARNLEV=1)
          OutTimeStamp = HcoDiagnStart
       ENDIF
    ENDIF

    CALL DiagnCollection_Create( HcoState%Diagn,                           &
                                 NX           = HcoState%NX,               &
                                 NY           = HcoState%NY,               &
                                 NZ           = HcoState%NZ,               &
                                 TS           = HcoState%TS_EMIS,          &
                                 AM2          = HcoState%Grid%AREA_M2%Val, &
                                 COL          = CollectionID,              &
                                 PREFIX       = TRIM(DiagnPrefix),         &
                                 deltaYMD     = deltaYMD,                  &
                                 deltaHMS     = deltaHMS,                  &
                                 OutTimeStamp = OutTimeStamp,              &
                                 RC           = RC                          )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Pass this collection ID to fixed variable for easy further
    ! reference to this collection
    HcoState%Diagn%HcoDiagnIDDefault = CollectionID

    ! ------------------------------------------------------------------
    ! HEMCO restart
    ! ------------------------------------------------------------------
#if defined ( ESMF_ )
    deltaYMD = 0
    deltaHMS = 1
#else
    deltaYMD = 99999999
    deltaHMS = 999999
#endif
#if defined( MODEL_GEOS )
    DiagnPrefix = 'HEMCO_restart.$YYYY$MM$DD$HH$MN.nc'
#else
    DiagnPrefix = 'HEMCO_restart'
#endif
    CALL DiagnCollection_Create( HcoState%Diagn,                           &
                                 NX           = HcoState%NX,               &
                                 NY           = HcoState%NY,               &
                                 NZ           = HcoState%NZ,               &
                                 TS           = HcoState%TS_EMIS,          &
                                 AM2          = HcoState%Grid%AREA_M2%Val, &
                                 COL          = CollectionID,              &
                                 PREFIX       = TRIM(DiagnPrefix),         &
                                 deltaYMD     = deltaYMD,                  &
                                 deltaHMS     = deltaHMS,                  &
                                 OutTimeStamp = HcoDiagnEnd,               &
                                 RC           = RC                          )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Pass this collection ID to fixed variable for easy further
    ! reference to this collection
    HcoState%Diagn%HcoDiagnIDRestart = CollectionID

    ! ------------------------------------------------------------------
    ! Manual diagnostics
    ! ------------------------------------------------------------------
#if defined ( ESMF_ )
    deltaYMD = 0
    deltaHMS = 1
#else
    deltaYMD = -1
    deltaHMS = -1
#endif
#if defined( MODEL_GEOS )
    DiagnPrefix = 'HEMCO_manual.$YYYY$MM$DD$HH$MN.nc'
#else
    DiagnPrefix = 'HEMCO_manual'
#endif
    CALL DiagnCollection_Create( HcoState%Diagn,                        &
                                 NX        = HcoState%NX,               &
                                 NY        = HcoState%NY,               &
                                 NZ        = HcoState%NZ,               &
                                 TS        = HcoState%TS_EMIS,          &
                                 AM2       = HcoState%Grid%AREA_M2%Val, &
                                 COL       = CollectionID,              &
                                 PREFIX    = TRIM(DiagnPrefix),         &
                                 deltaYMD  = deltaYMD,                  &
                                 deltaHMS  = deltaHMS,                  &
                                 RC        = RC                          )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Pass this collection ID to fixed variable for easy further
    ! reference to this collection
    HcoState%Diagn%HcoDiagnIDManual  = CollectionID

    ! ------------------------------------------------------------------
    ! Now that collections are defined, add diagnostics specified in the
    ! HEMCO diagnostics definition file. The latter can be specified in
    ! the HEMCO configuration file. These diagnostics are all written
    ! into the default HEMCO collection.
    ! ------------------------------------------------------------------
    CALL Diagn_DefineFromConfig( HcoState, RC )
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
  SUBROUTINE Diagn_DefineFromConfig( HcoState, RC )
!
! !USES:
!
    USE HCO_CharTools_Mod
    USE CHARPAK_Mod,       ONLY : STRREPL, STRSPLIT
    USE inquireMod,        ONLY : findFreeLUN
    USE HCO_STATE_MOD,     ONLY : HCO_GetHcoID
    USE HCO_STATE_MOD,     ONLY : HCO_State
    USE HCO_EXTLIST_MOD,   ONLY : GetExtOpt
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState    ! HEMCO state object
    INTEGER,          INTENT(INOUT)           :: RC          ! Failure or success
!
! !REVISION HISTORY:
!  10 Apr 2015 - C. Keller   - Initial version
!  21 Feb 2016 - C. Keller   - Added default diagnostics (optional)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                    :: N,      LUN
    LOGICAL                    :: EOF,    FOUND,   DefaultSet
    CHARACTER(LEN=31)          :: SpcName, OutUnit
    CHARACTER(LEN=63)          :: cName,  dLname, dSname
    INTEGER                    :: HcoID,  ExtNr,   Cat, Hier, SpaceDim
    CHARACTER(LEN=255)         :: lName,  LOC,     MSG

    !=================================================================
    ! Diagn_DefineFromConfig begins here!
    !=================================================================

    ! Init
    LOC = 'Diagn_DefineFromConfig (hco_diagn_mod.F90)'

    ! Load DiagnFile into buffer
    CALL DiagnFileOpen( HcoState%Config, LUN, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If defined, sequentially get all entries
    IF ( LUN > 0 ) THEN

       ! Do for every line
       DO

          ! Get next line
          CALL DiagnFileGetNext( HcoState%Config, LUN,     cName,      &
                                 SpcName,         ExtNr,   Cat, Hier,  &
                                 SpaceDim,        OutUnit, EOF, RC,    &
                                 lName=lName )
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
          CALL Diagn_Create( HcoState,                      &
                             cName     = cName,             &
                             long_name = lName,             &
                             HcoID     = HcoID,             &
                             ExtNr     = ExtNr,             &
                             Cat       = Cat,               &
                             Hier      = Hier,              &
                             SpaceDim  = SpaceDim,          &
                             OutUnit   = OutUnit,           &
                             AutoFill  = 1,                 &
                             COL       = HcoState%Diagn%HcoDiagnIDDefault, &
                             RC        = RC                  )
          IF ( RC /= HCO_SUCCESS ) RETURN

       ENDDO

       ! Close file
       CALL DiagnFileClose ( LUN )

    ENDIF ! LUN > 0

    ! ---------------------------------------------------------------------
    ! Eventually prepare a diagnostics for every HEMCO species.
    ! This is optional and controlled by HEMCO setting DefaultDiagnSet.
    ! ---------------------------------------------------------------------
    CALL GetExtOpt( HcoState%Config, -999, 'DefaultDiagnOn', &
                    OptValBool=DefaultSet, FOUND=FOUND, RC=RC )
    IF ( .NOT. FOUND ) DefaultSet = .FALSE.
    IF ( DefaultSet ) THEN

       ! Search for default diagnostics variable prefix
       CALL GetExtOpt( HcoState%Config, -999, 'DefaultDiagnSname', &
                       OptValChar=dSname, FOUND=FOUND, RC=RC )
       IF ( .NOT. FOUND ) dSname = 'HEMCO_EMIS_'

       CALL GetExtOpt( HcoState%Config, -999, 'DefaultDiagnLname', &
                       OptValChar=dLname, FOUND=FOUND, RC=RC )
       IF ( .NOT. FOUND ) dLname = 'HEMCO_emissions_of_species_'

       ! Search for default diagnostics dimension
       CALL GetExtOpt( HcoState%Config, -999, 'DefaultDiagnDim', &
                       OptValInt=SpaceDim, FOUND=FOUND, RC=RC )
       IF ( .NOT. FOUND ) SpaceDim = 3
       SpaceDim = MAX(MIN(SpaceDim,3),2)

       ! Get units
       CALL GetExtOpt( HcoState%Config, -999, 'DefaultDiagnUnit', &
                       OptValChar=OutUnit, FOUND=FOUND, RC=RC )
       IF ( .NOT. FOUND ) OutUnit = 'kg m-2 s-1'

       ! Loop over all species and create diagnostics
       DO N = 1, HcoState%nSpc
          cName = TRIM(dSname)//TRIM(HcoState%Spc(N)%SpcName)
          lName = TRIM(dLname)//TRIM(HcoState%Spc(N)%SpcName)

          CALL Diagn_Create( HcoState,                                     &
                             cName     = cName,                            &
                             long_name = lName,                            &
                             HcoID     = HcoState%Spc(N)%HcoID,            &
                             ExtNr     = -1,                               &
                             Cat       = -1,                               &
                             Hier      = -1,                               &
                             SpaceDim  = SpaceDim,                         &
                             OutUnit   = OutUnit,                          &
                             AutoFill  = 1,                                &
                             COL       = HcoState%Diagn%HcoDiagnIDDefault, &
                             RC        = RC                                 )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDDO
    ENDIF

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
  SUBROUTINE Diagn_Create( HcoState,  cName,                  &
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
    TYPE(HCO_State),  POINTER                 :: HcoState      ! HEMCO state obj.
    CHARACTER(LEN=*), INTENT(IN   )           :: cName         ! Diagnostics name
    CHARACTER(LEN=*), INTENT(IN   )           :: OutUnit       ! Output units
    INTEGER,          INTENT(IN   ), OPTIONAL :: SpaceDim      ! Spatial dimension
    INTEGER,          INTENT(IN   ), OPTIONAL :: ExtNr         ! Extension #
    INTEGER,          INTENT(IN   ), OPTIONAL :: Cat           ! Category
    INTEGER,          INTENT(IN   ), OPTIONAL :: Hier          ! Hierarchy
    INTEGER,          INTENT(IN   ), OPTIONAL :: HcoID         ! HEMCO species ID
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
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(DiagnCont),       POINTER :: ThisDiagn
    TYPE(DiagnCont),       POINTER :: TmpDiagn
    TYPE(DiagnCollection), POINTER :: ThisColl

    ! Scalars
    CHARACTER(LEN=255)             :: LOC, MSG
    INTEGER                        :: PS, Flag
    REAL(hp)                       :: Scal
    REAL(hp)                       :: MWg, EmMWg, MolR
    LOGICAL                        :: ForceMean, FOUND

    !======================================================================
    ! Diagn_Create begins here!
    !======================================================================

    ! Nullify
    ThisDiagn => NULL()
    TmpDiagn  => NULL()
    ThisColl  => NULL()

    ! Init
    LOC = 'Diagn_Create (hco_diagn_mod.F90)'
    CALL DiagnCollection_DefineID( HcoState%Diagn, PS, RC, COL=COL, &
                                   HcoState=HcoState, InUse=FOUND, ThisColl=ThisColl )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Error if collection does not exist
    IF ( .NOT. FOUND ) THEN
       WRITE(MSG,*) 'Cannot create diagnostics ', TRIM(cName), &
                    ' - collection does not exist: ', PS
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    !----------------------------------------------------------------------
    ! Check if diagnostics already exists
    !----------------------------------------------------------------------
    IF ( PRESENT(OkIfExist) ) THEN
       IF ( OkIfExist ) THEN
          IF ( PRESENT(cID) ) THEN
             CALL DiagnCont_Find( HcoState%Diagn, cID, -1, -1, -1, -1, &
                           '', -1, FOUND, TmpDiagn, COL=PS )
          ELSE
             CALL DiagnCont_Find( HcoState%Diagn, -1, -1, -1, -1, -1, &
                           TRIM(ADJUSTL(cName)), -1, FOUND, TmpDiagn, COL=PS )
          ENDIF
          TmpDiagn => NULL()

          ! Exit if found
          IF ( FOUND ) THEN
             IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
                WRITE(MSG,*) 'Diagnostics already exists - ', &
                             'will not be added again: ', TRIM(cName)
                CALL HCO_MSG ( HcoState%Config%Err, MSG )
             ENDIF
             RC = HCO_SUCCESS
             RETURN
          ENDIF
       ENDIF
! End conflict (ewl, 1/8/15)
    ENDIF

    !----------------------------------------------------------------------
    ! Initalize diagnostics container.
    !----------------------------------------------------------------------
    CALL DiagnCont_Init( ThisDiagn )

    !----------------------------------------------------------------------
    ! Pass input variables
    !----------------------------------------------------------------------
    ThisDiagn%cName   = ADJUSTL(cName)
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
       ThisDiagn%long_name = TRIM(ADJUSTL(cName))
    ENDIF

    !----------------------------------------------------------------------
    ! Eventually link to data array. This will disable all time averaging,
    ! unit conversions, etc. (data will just be returned as is).
    !----------------------------------------------------------------------
    IF ( PRESENT(Trgt2D) ) THEN
       CALL DiagnCont_Link_2D( ThisDiagn, ThisColl, Trgt2D, RC, &
                               HcoState=HcoState )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF
    IF ( PRESENT(Trgt3D) ) THEN
       CALL DiagnCont_Link_3D( ThisDiagn, ThisColl, Trgt3D, RC, &
                               HcoState=HcoState )
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
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       IF ( TRIM(OutOper) == 'CumulSum' ) THEN
          MSG = 'Cannot use scale factor on diagnostics that '// &
                'are cumulative sums: '//TRIM(cName)
          CALL HCO_ERROR( HcoState%config%Err, MSG, RC, THISLOC=LOC )
          RETURN
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
             CALL HCO_ERROR( HcoState%config%Err, MSG, RC, THISLOC=LOC )
             RETURN
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
             MWg   = HcoState%Spc(HcoID)%MW_g
             EmMWg = HcoState%Spc(HcoID)%EmMW_g
             MolR  = HcoState%Spc(HcoID)%MolecRatio
          ELSE
                MWg   = MW_g
                EmMWg = EmMW_g
                MolR  = MolecRatio
          ENDIF

          CALL HCO_UNIT_GetMassScal( HcoState%Config,             &
                    unt         = OutUnit,                        &
                    MW_IN       = MWg,                            &
                    MW_OUT      = EmMWg,                          &
                    MOLEC_RATIO = MolR,                           &
                    Scal        = Scal                             )
          IF ( Scal <= 0.0_hp ) THEN
             MSG = 'Cannot find mass scale factor for unit '//TRIM(OutUnit)
             CALL HCO_ERROR( HcoState%config%Err, MSG, RC, THISLOC=LOC )
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
             CALL HCO_ERROR( HcoState%config%Err, MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF
       ENDIF ! OutOper not set
    ENDIF ! .NOT. DtaIsPtr

    !-----------------------------------------------------------------------
    ! Make sure that there is no other diagnostics with this name
    !-----------------------------------------------------------------------
    CALL DiagnCont_Find( HcoState%Diagn, -1, -1, -1, -1, -1, &
                        Trim(ADJUSTL(cName)), -1, FOUND, TmpDiagn, COL=PS )
    IF ( FOUND ) THEN
!       MSG = 'There is already a diagnostics with this name: ' // TRIM(cName)
!       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
!       RETURN
       ThisDiagn%cName = trim(cName) // '_a'
       MSG = 'Changed Diagn name to ' // trim(ThisDiagn%cName)
       CALL HCO_MSG ( HcoState%Config%Err, MSG )
    ENDIF

    !-----------------------------------------------------------------------
    ! Set container ID (if defined). There must not be two containers with
    ! the same container ID.
    !-----------------------------------------------------------------------
    IF ( PRESENT(cID) ) THEN
       IF ( cID > 0 ) THEN

          ! Check if there is already a diagnostics with this container ID.
          TmpDiagn => NULL()
          CALL DiagnCont_Find( HcoState%Diagn, cID, -1, -1, -1, -1, &
                               '', -1, FOUND, TmpDiagn, COL=PS )
          IF ( FOUND ) THEN
             WRITE(MSG,*) 'Diagnostics ', TRIM(TmpDiagn%cName), ' already has ID ', &
                cID, ' - cannot create diagnostics ', TRIM(cName)

             CALL HCO_ERROR( HcoState%config%Err, MSG, RC, THISLOC=LOC )
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
    IF ( HCO_IsVerb(HcoState%Config%Err,1) ) THEN
       WRITE(MSG,'(a, i4)') 'Successfully added diagnostic '// &
                             TRIM(ThisDiagn%cName) // ' to collection ', PS
       CALL HCO_MSG ( HcoState%Config%Err, MSG )
    ENDIF

    ! Cleanup
    ThisDiagn => NULL()
    ThisColl  => NULL()

    ! Return
    RC  = HCO_SUCCESS

  END SUBROUTINE Diagn_Create
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_UpdateSp0d
!
! !DESCRIPTION: Subroutine Diagn\_UpdateSp0d is the wrapper routine to update
! the diagnostics for single precision scalar values.  It invokes the main
! diagnostics update routine with the appropriate arguments.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_UpdateSp0d( HcoState,    cID,    cName,   ExtNr, &
                               Cat,         Hier,   HcoID,   AutoFill,       &
                               Scalar,      Total,  PosOnly, COL,            &
                               MinDiagnLev, RC                         )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState       ! HEMCO state obj
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
    REAL(sp),         INTENT(IN   )           :: Scalar         ! 0D scalar
    REAL(sp),         INTENT(IN   ), OPTIONAL :: Total          ! Total
    LOGICAL,          INTENT(IN   ), OPTIONAL :: PosOnly        ! Use only vals
                                                                !  >= 0?
    INTEGER,          INTENT(IN   ), OPTIONAL :: COL            ! Collection Nr.
    INTEGER,          INTENT(IN   ), OPTIONAL :: MinDiagnLev    ! minimum diagn level
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC             ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2015 - C. Keller - Initialization
!  19 Sep 2016 - R. Yantosca - Rewritten for Gfortran: remove Scalar and
!                              Array3d (put those in other overloaded methods)
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Call down to driver routine
    CALL Diagn_UpdateDriver( HcoState,                  &
                             cID         = cID,         &
                             cName       = cName,       &
                             ExtNr       = ExtNr,       &
                             Cat         = Cat,         &
                             Hier        = Hier,        &
                             HcoID       = HcoID,       &
                             AutoFill    = AutoFill,    &
                             Scalar_SP   = Scalar,      &
                             Total_SP    = Total,       &
                             PosOnly     = PosOnly,     &
                             COL         = COL,         &
                             MinDiagnLev = MinDiagnLev, &
                             RC          = RC )

  END SUBROUTINE Diagn_UpdateSp0d
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_UpdateSp2d
!
! !DESCRIPTION: Subroutine Diagn\_UpdateSp2d is the wrapper routine to update
! the diagnostics for single precision 2-D arrays.  It invokes the main
! diagnostics update routine with the appropriate arguments.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_UpdateSp2d( HcoState,    cID,   cName,   ExtNr,     &
                               Cat,         Hier,  HcoID,   AutoFill,  &
                               Array2D,     Total, PosOnly, COL,       &
                               MinDiagnLev, RC                        )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState       ! HEMCO state obj
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
    REAL(sp),         INTENT(IN   )           :: Array2D(:,:)   ! 2D array
    REAL(sp),         INTENT(IN   ), OPTIONAL :: Total          ! Total
    LOGICAL,          INTENT(IN   ), OPTIONAL :: PosOnly        ! Use only vals
                                                                !  >= 0?
    INTEGER,          INTENT(IN   ), OPTIONAL :: COL            ! Collection Nr.
    INTEGER,          INTENT(IN   ), OPTIONAL :: MinDiagnLev    ! minimum diagn level
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC             ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2015 - C. Keller - Initialization
!  19 Sep 2016 - R. Yantosca - Rewritten for Gfortran: remove Scalar and
!                              Array3d (put those in other overloaded methods)
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Call down to driver routine
    CALL Diagn_UpdateDriver( HcoState,                  &
                             cID         = cID,         &
                             cName       = cName,       &
                             ExtNr       = ExtNr,       &
                             Cat         = Cat,         &
                             Hier        = Hier,        &
                             HcoID       = HcoID,       &
                             AutoFill    = AutoFill,    &
                             Array2D_SP  = Array2D,     &
                             Total_SP    = Total,       &
                             PosOnly     = PosOnly,     &
                             COL         = COL,         &
                             MinDiagnLev = MinDiagnLev, &
                             RC          = RC )

  END SUBROUTINE Diagn_UpdateSp2d
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_UpdateSp3d
!
! !DESCRIPTION: Subroutine Diagn\_UpdateSp is the wrapper routine to update
! the diagnostics for single precision 3-D arrays. It invokes the main
! diagnostics update routine with the appropriate arguments.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_UpdateSp3d( HcoState,    cID,   cName,   ExtNr,     &
                               Cat,         Hier,  HcoID,   AutoFill,  &
                               Array3D,     Total, PosOnly, COL,       &
                               MinDiagnLev, RC                        )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState       ! HEMCO state obj
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
    REAL(sp),         INTENT(IN   )           :: Array3D(:,:,:) ! 3D array
    REAL(sp),         INTENT(IN   ), OPTIONAL :: Total          ! Total
    LOGICAL,          INTENT(IN   ), OPTIONAL :: PosOnly        ! Use only vals
                                                                !  >= 0?
    INTEGER,          INTENT(IN   ), OPTIONAL :: COL            ! Collection Nr.
    INTEGER,          INTENT(IN   ), OPTIONAL :: MinDiagnLev    ! minimum diagn level
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC             ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2015 - C. Keller - Initialization
!  19 Sep 2016 - R. Yantosca - Rewritten for Gfortran: remove Scalar and
!                              Array2d (put those in other overloaded methods)
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Call down to driver routine
    CALL Diagn_UpdateDriver( HcoState,                  &
                             cID         = cID,         &
                             cName       = cName,       &
                             ExtNr       = ExtNr,       &
                             Cat         = Cat,         &
                             Hier        = Hier,        &
                             HcoID       = HcoID,       &
                             AutoFill    = AutoFill,    &
                             Array3D_SP  = Array3D,     &
                             Total_SP    = Total,       &
                             PosOnly     = PosOnly,     &
                             COL         = COL,         &
                             MinDiagnLev = MinDiagnLev, &
                             RC          = RC )

  END SUBROUTINE Diagn_UpdateSp3d
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_UpdateDp0d
!
! !DESCRIPTION: Subroutine Diagn\_UpdateSp0d is the wrapper routine to update
! the diagnostics for double-precision scalar values.  It invokes the main
! diagnostics update routine with the appropriate arguments.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_UpdateDp0d( HcoState,    cID,   cName,   ExtNr,     &
                               Cat,         Hier,  HcoID,   AutoFill,  &
                               Scalar,      Total, PosOnly, COL,       &
                               MinDiagnLev, RC                        )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState       ! HEMCO state obj
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
    REAL(dp),         INTENT(IN   )           :: Scalar         ! 1D scalar
    REAL(dp),         INTENT(IN   ), OPTIONAL :: Total          ! Total
    LOGICAL,          INTENT(IN   ), OPTIONAL :: PosOnly        ! Use only vals
                                                                !  >= 0?
    INTEGER,          INTENT(IN   ), OPTIONAL :: COL            ! Collection Nr.
    INTEGER,          INTENT(IN   ), OPTIONAL :: MinDiagnLev    ! minimum diagn level
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC             ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2015 - C. Keller   - Initialization
!  19 Sep 2016 - R. Yantosca - Rewritten for Gfortran: remove Array2d and
!                              Array3d (put those in other overloaded methods)
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Call down to driver routine
    CALL Diagn_UpdateDriver( HcoState,                  &
                             cID         = cID,         &
                             cName       = cName,       &
                             ExtNr       = ExtNr,       &
                             Cat         = Cat,         &
                             Hier        = Hier,        &
                             HcoID       = HcoID,       &
                             AutoFill    = AutoFill,    &
                             Scalar      = Scalar,      &
                             Total       = Total,       &
                             PosOnly     = PosOnly,     &
                             COL         = COL,         &
                             MinDiagnLev = MinDiagnLev, &
                             RC          = RC )

  END SUBROUTINE Diagn_UpdateDp0d
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_UpdateDp2d
!
! !DESCRIPTION: Subroutine Diagn\_UpdateSp2d is the wrapper routine to update
! the diagnostics for single precision 2D arrays.  It invokes the main
! diagnostics update routine with the appropriate arguments.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_UpdateDp2d( HcoState,    cID,   cName,   ExtNr,     &
                               Cat,         Hier,  HcoID,   AutoFill,  &
                               Array2D,     Total, PosOnly, COL,       &
                               MinDiagnLev, RC                        )
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState       ! HEMCO state obj
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
    REAL(dp),         INTENT(IN   )           :: Array2D(:,:)   ! 2D array
    REAL(dp),         INTENT(IN   ), OPTIONAL :: Total          ! Total
    LOGICAL,          INTENT(IN   ), OPTIONAL :: PosOnly        ! Use only vals
                                                                !  >= 0?
    INTEGER,          INTENT(IN   ), OPTIONAL :: COL            ! Collection Nr.
    INTEGER,          INTENT(IN   ), OPTIONAL :: MinDiagnLev    ! minimum diagn level
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC             ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2015 - C. Keller   - Initialization
!  19 Sep 2016 - R. Yantosca - Rewritten for Gfortran: remove Scalar and
!                              Array3d (put those in other overloaded methods)
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Call down to driver routine
    CALL Diagn_UpdateDriver( HcoState,                  &
                             cID         = cID,         &
                             cName       = cName,       &
                             ExtNr       = ExtNr,       &
                             Cat         = Cat,         &
                             Hier        = Hier,        &
                             HcoID       = HcoID,       &
                             AutoFill    = AutoFill,    &
                             Array2D     = Array2D,     &
                             Total       = Total,       &
                             PosOnly     = PosOnly,     &
                             COL         = COL,         &
                             MinDiagnLev = MinDiagnLev, &
                             RC          = RC )

  END SUBROUTINE Diagn_UpdateDp2d
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diagn_UpdateDp3d
!
! !DESCRIPTION: Subroutine Diagn\_UpdateSp3d is the wrapper routine to update
! the diagnostics for single precision arrays. It invokes the main diagnostics
! update routine with the appropriate arguments.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagn_UpdateDp3d( HcoState,    cID,   cName,   ExtNr,           &
                               Cat,         Hier,  HcoID,   AutoFill,        &
                               Array3D,     Total, PosOnly, COL,             &
                               MinDiagnLev, RC                        )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState       ! HEMCO state obj
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
    REAL(dp),         INTENT(IN   )           :: Array3D(:,:,:) ! 3D array
    REAL(dp),         INTENT(IN   ), OPTIONAL :: Total          ! Total
    LOGICAL,          INTENT(IN   ), OPTIONAL :: PosOnly        ! Use only vals
                                                                !  >= 0?
    INTEGER,          INTENT(IN   ), OPTIONAL :: COL            ! Collection Nr.
    INTEGER,          INTENT(IN   ), OPTIONAL :: MinDiagnLev    ! minimum diagn level
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC             ! Return code
!
! !REVISION HISTORY:
!  20 Apr 2015 - C. Keller   - Initialization
!  19 Sep 2016 - R. Yantosca - Rewritten for Gfortran: remove Scalar and
!                              Array2d (put those in other overloaded methods)
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Call down to driver routine
    CALL Diagn_UpdateDriver( HcoState,                  &
                             cID         = cID,         &
                             cName       = cName,       &
                             ExtNr       = ExtNr,       &
                             Cat         = Cat,         &
                             Hier        = Hier,        &
                             HcoID       = HcoID,       &
                             AutoFill    = AutoFill,    &
                             Array3D     = Array3D,     &
                             Total       = Total,       &
                             PosOnly     = PosOnly,     &
                             COL         = COL,         &
                             MinDiagnLev = MinDiagnLev, &
                             RC          = RC )

  END SUBROUTINE Diagn_UpdateDp3d
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
! HEMCO species ID (HcoID) is provided, all containers with the same
! combination of HcoID, extension number (ExtNr), emission category
! (Cat) and hierarchy (Hier) are updated. If no valid HcoID and no
! valid cID is given, the container name has to be provided. The passed
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
! The optional input argument `MinDiagnLev` determines how `deep`
! this routine will search for diagnostics with matching HcoID, ExtNr,
! etc. For example, if a HcoID, an ExtNr, and a category is provided,
! HEMCO by default will only update diagnostics containers with exactly
! the same HcoID, ExtNr, and category - but not diagnostics of `lower
! level`, e.g. with the same HcoID and ExtNr but no assigned category.
! This behavior can be changed by explicitly setting MinDiagnLev to the
! minimum diagnostics level. In the given example, setting MinDiagnLev
! to 1 would also update level 1 and level 2 diagnostics of the same
! HcoID (e.g. diagnostics with the same HcoID and no assigned ExtNr and
! category; as well as diagnostics with the same HcoID and ExtNr and no
! assigned category).
!\\
!\\
! Notes:
! \begin{itemize}
! \item For a given time step, the same diagnostics container can be
!       updated multiple times. The field average is always defined as
!       temporal average, e.g. multiple updates on the same time step
!       will not increase the averaging weight of that time step.
!
! \item If the passed array is empty (i.e. not associated), it is
!       treated as empty values (i.e. zeros).
!
! \item The collection number can be set to -1 to scan trough all
!       existing diagnostic collections.
! \end{itemize}
!
! !INTERFACE:
!
  SUBROUTINE Diagn_UpdateDriver( HcoState,  cID,        cName,                   &
                                 ExtNr,     Cat,        Hier,       HcoID,       &
                                 AutoFill,  Scalar,     Array2D,    Array3D,     &
                                 Total,     Scalar_SP,  Array2D_SP, Array3D_SP,  &
                                 Total_SP,  Scalar_HP,  Array2D_HP, Array3D_HP, &
                                 Total_HP,  PosOnly,    COL,        MinDiagnLev, &
                                 RC                                               )
!
! !USES:
!
    USE HCO_Arr_Mod,   ONLY : HCO_ArrAssert
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                         :: HcoState          ! HEMCO state obj
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
    REAL(hp),         INTENT(IN   ), OPTIONAL         :: Scalar_HP         ! 1D scalar
    REAL(hp),         INTENT(IN   ), OPTIONAL, TARGET :: Array2D_HP(:,:)   ! 2D array
    REAL(hp),         INTENT(IN   ), OPTIONAL, TARGET :: Array3D_HP(:,:,:) ! 3D array
    REAL(hp),         INTENT(IN   ), OPTIONAL         :: Total_HP          ! Total
    LOGICAL,          INTENT(IN   ), OPTIONAL         :: PosOnly           ! Use only vals
                                                                           ! >= 0?
    INTEGER,          INTENT(IN   ), OPTIONAL         :: COL               ! Collection Nr.
    INTEGER,          INTENT(IN   ), OPTIONAL         :: MinDiagnLev       ! Collection Nr.
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
!  17 Jun 2015 - C. Keller - Added argument MinDiagnLev
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(DiagnCollection), POINTER :: ThisColl
    TYPE(DiagnCont),       POINTER :: ThisDiagn
    REAL(sp),              POINTER :: Arr2D (:,:)
    REAL(sp),              POINTER :: Tmp2D (:,:)
    REAL(sp),              POINTER :: Arr3D (:,:,:)
    REAL(sp)                       :: TmpScalar

    ! Scalars
    CHARACTER(LEN=255)             :: LOC, MSG
    REAL(hp)                       :: Fact
    REAL(hp)                       :: Tmp
    CHARACTER(LEN=63)              :: DgnName
    INTEGER                        :: I, J, L, PS, AS
    INTEGER                        :: DgncID,  DgnExtNr, DgnCat
    INTEGER                        :: iHier,   iExt,     iCat
    INTEGER                        :: DgnHier, DgnHcoID
    INTEGER                        :: ThisUpdateID
    INTEGER                        :: AutoFlag
    INTEGER                        :: CNT
    INTEGER                        :: MnDgnLev, OrigDgnLev, ThisDgnLev
    LOGICAL                        :: Found, OnlyPos, VertSum, IsAssoc, IsNewTS
    LOGICAL                        :: InUse, SearchAll

    !======================================================================
    ! Diagn_UpdateDriver begins here!
    !======================================================================

    ! Init
    LOC       =  'Diagn_UpdateDriver (hco_diagn_mod.F90)'
    RC        =  HCO_SUCCESS
    ThisColl  => NULL()
    ThisDiagn => NULL()
    Arr2D     => NULL()
    Tmp2D     => NULL()
    Arr3D     => NULL()

    ! Get collection number.
    CALL DiagnCollection_DefineID( HcoState%Diagn, PS, RC, COL=COL, DEF=-1, &
           OKIfAll=.TRUE., InUse=InUse, ThisColl=ThisColl, HcoState=HcoState )
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
    IF ( PRESENT(cName   ) ) DgnName  = ADJUSTL(cName)
    IF ( PRESENT(cID     ) ) DgncID   = cID
    IF ( PRESENT(ExtNr   ) ) DgnExtNr = ExtNr
    IF ( PRESENT(Cat     ) ) DgnCat   = Cat
    IF ( PRESENT(Hier    ) ) DgnHier  = Hier
    IF ( PRESENT(HcoID   ) ) DgnHcoID = HcoID
    IF ( PRESENT(PosOnly ) ) OnlyPos  = PosOnly
    IF ( PRESENT(AutoFill) ) AutoFlag = AutoFill

    ! Get the update time ID.
    CALL HcoClock_Get( HcoState%Clock, nSteps=ThisUpdateID, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Count # of containers that are updated
    CNT = 0

    !-----------------------------------------------------------------
    ! Diagnostics levels to be used. By default, use only diagnostics
    ! at the provided level. For instance, if a hierarchy number is
    ! given do not update diagnostics with the same species and
    ! extension number but a hierarchy number of -1. If a diagnostics
    ! level is given, update all diagnostics up to this diagnostics
    ! level.
    !-----------------------------------------------------------------

    ! Get original diagnostics level
    OrigDgnLev = 999
    IF ( DgnHcoID > -1 ) THEN
       OrigDgnLev = 1
       IF ( DgnExtNr > -1 ) THEN
          OrigDgnLev = 2
          IF ( DgnCat > -1 ) THEN
             OrigDgnLev = 3
             IF ( DgnHier > -1 ) THEN
                OrigDgnLev = 4
             ENDIF
          ENDIF
       ENDIF
    ENDIF

    ! Set diagnostics level
    IF ( PRESENT(MinDiagnLev) ) THEN
       MnDgnLev = MinDiagnLev
    ELSE
       MnDgnLev = OrigDgnLev
    ENDIF


    !-----------------------------------------------------------------
    ! Loop over collections
    !-----------------------------------------------------------------
    DO WHILE ( ASSOCIATED(ThisColl) )

       ! Reset Diagnostics
       ThisDiagn => NULL()

       ! Reset diagnostics level to use
       ThisDgnLev = OrigDgnLev

       !-----------------------------------------------------------------
       ! Do for every container in the diagnostics list that matches the
       ! specified arguments (ID, ExtNr, etc.). This can be more than one
       ! container (ckeller, 09/25/2014).
       !-----------------------------------------------------------------
       DO

          ! Set ExtNr, Cat, Hier based on current diagnostics level.
          iExt  = -1
          iCat  = -1
          iHier = -1
          IF ( ThisDgnLev > 1 ) iExt  = DgnExtNr
          IF ( ThisDgnLev > 2 ) iCat  = DgnCat
          IF ( ThisDgnLev > 3 ) iHier = DgnHier

          ! Search for diagnostics that matches the given arguments.
          ! If ThisDiagn is empty (first call), the search will start
          ! at the first diagnostics container. Otherwise, the search
          ! will resume from this diagnostics container.
          CALL DiagnCont_Find( HcoState%Diagn, &
                               DgncID,    iExt,     iCat,     iHier,   &
                               DgnHcoID,  DgnName,  AutoFlag, Found,   &
                               ThisDiagn, RESUME=.TRUE., COL=ThisColl%CollectionID )

          ! Exit while loop if no diagnostics found
          !IF ( .NOT. Found ) EXIT
          ! Now also check lower level diagnostics is specified so
          IF ( .NOT. Found ) THEN
             IF ( ThisDgnLev > MnDgnLev ) THEN
                ThisDgnLev =  ThisDgnLev - 1
                ThisDiagn  => NULL()
                CYCLE
             ELSE
                EXIT
             ENDIF
          ENDIF

          ! If container holds just a pointer to external data, don't do
          ! anything!
          IF ( ThisDiagn%DtaIsPtr ) THEN
             MSG = 'You try to update a container that holds a '  // &
                  'pointer to data - this should never happen! ' // &
                  TRIM(ThisDiagn%cName)
             CALL HCO_WARNING( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
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
             ELSEIF ( PRESENT(Array3D_HP) ) THEN
                ALLOCATE( Arr3D(ThisColl%NX,ThisColl%NY,ThisColl%NZ),STAT=AS)
                IF ( AS /= 0 ) THEN
                   CALL HCO_ERROR( HcoState%Config%Err,&
                                   'Allocation error Arr3D', RC, THISLOC=LOC )
                   RETURN
                ENDIF
                Arr3D = Array3D_HP
             ELSEIF( PRESENT(Array3D) ) THEN
                ALLOCATE( Arr3D(ThisColl%NX,ThisColl%NY,ThisColl%NZ),STAT=AS)
                IF ( AS /= 0 ) THEN
                   CALL HCO_ERROR( HcoState%Config%Err,&
                                   'Allocation error Arr3D', RC, THISLOC=LOC )
                   RETURN
                ENDIF
                Arr3D = Array3D
             ENDIF

             ! 2D array
             IF ( PRESENT(Array2D_SP) ) THEN
                Arr2D => Array2D_SP
             ELSEIF ( PRESENT(Array2D_HP) ) THEN
                ALLOCATE( Arr2D(ThisColl%NX,ThisColl%NY),STAT=AS)
                IF ( AS /= 0 ) THEN
                   CALL HCO_ERROR( HcoState%Config%Err,&
                                   'Allocation error Arr2D', RC, THISLOC=LOC )
                   RETURN
                ENDIF
                Arr2D = Array2D_HP
             ELSEIF( PRESENT(Array2D) ) THEN
                ALLOCATE( Arr2D(ThisColl%NX,ThisColl%NY),STAT=AS)
                IF ( AS /= 0 ) THEN
                   CALL HCO_ERROR( HcoState%Config%Err,&
                                   'Allocation error Arr2D', RC, THISLOC=LOC )
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
             IF ( PRESENT(Array3D_SP) .OR. PRESENT(Array3D) .OR. PRESENT(Array3D_HP) ) THEN

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

             IF ( PRESENT(Array3D_SP) .OR. PRESENT(Array3D) .OR. PRESENT(Array3D_HP) .OR. &
                  PRESENT(Array2D_SP) .OR. PRESENT(Array2D) .OR. PRESENT(Array2D_HP)      ) THEN

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
                IF ( PRESENT(Array2D) .OR. PRESENT(Array2D_SP) .OR. PRESENT(Array2D_HP) ) THEN
                   IF ( .NOT. ASSOCIATED(Arr2D) ) THEN
                      IsAssoc = .FALSE.
                   ELSE
                      Tmp2D => Arr2D
                   ENDIF
                ELSEIF ( PRESENT(Array3D) .OR. PRESENT(Array3D_SP) .OR. PRESENT(Array3D_HP) ) THEN
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
                   CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
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
             IF ( PRESENT(Scalar_SP) .OR. PRESENT(Scalar) .OR. PRESENT(Scalar_HP) ) THEN

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
          IF ( PRESENT(Total_HP) ) THEN
             ThisDiagn%Total = ThisDiagn%Total + Total_HP
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
          IF ( HCO_IsVerb(HcoState%Config%Err,2) ) THEN
             WRITE(MSG,'(a,a,a,I3,a)') 'Successfully updated diagnostics: ', &
                TRIM(ThisDiagn%cName), ' (counter:', ThisDiagn%Counter, ')'
             CALL HCO_MSG ( HcoState%Config%Err, MSG )
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
    ELSEIF (PRESENT(Array3D_HP) ) THEN
       IF ( ASSOCIATED(Arr3D) ) DEALLOCATE(Arr3D)
    ELSEIF (PRESENT(Array3D) ) THEN
       IF ( ASSOCIATED(Arr3D) ) DEALLOCATE(Arr3D)
    ENDIF
    IF (PRESENT(Array2D_SP) ) THEN
       Arr2D => NULL()
    ELSEIF (PRESENT(Array2D_HP) ) THEN
       IF ( ASSOCIATED(Arr2D) ) DEALLOCATE(Arr2D)
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
  SUBROUTINE Diagn_Get( HcoState,                          &
                        EndOfIntvOnly,            DgnCont, &
                        FLAG,      RC,            cName,   &
                        cID,       AutoFill,      COL,     &
                        SkipZeroCount                       )
!
! !USES:
!
    USE HCO_STATE_MOD,   ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),  POINTER                 :: HcoState        ! HEMCO state obj
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
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DiagnCollection), POINTER :: ThisColl
    INTEGER                        :: PS, AF
    LOGICAL                        :: TimeToWrite
    LOGICAL                        :: FOUND, CF
    LOGICAL                        :: SKIPZERO

    !======================================================================
    ! Diagn_Get begins here!
    !======================================================================

    ! Init
    FLAG     =  HCO_FAIL
    RC       =  HCO_SUCCESS
    CF       = .FALSE.
    ThisColl => NULL()

    ! Get collection number
    CALL DiagnCollection_DefineID( HcoState%Diagn, PS, RC, COL=COL, &
                                   ThisColl=ThisColl, HcoState=HcoState )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set AutoFill flag
    AF = -1
    IF ( PRESENT(AutoFill  ) ) AF = AutoFill

    ! Check if diagnostics with counter = 0 shall be skipped
    SKIPZERO = .FALSE.
    IF ( PRESENT(SkipZeroCount) ) SKIPZERO = SkipZeroCount

    IF ( EndOfIntvOnly ) THEN
       TimeToWrite =  DiagnCollection_IsTimeToWrite( HcoState, PS )
       IF ( .NOT. TimeToWrite ) THEN
          DgnCont => NULL()
          RETURN
       ENDIF
    ENDIF

    ! If container name is given, search for diagnostics with
    ! the given name.
    IF ( PRESENT( cName ) ) THEN
       CALL DiagnCont_Find( HcoState%Diagn, &
                            -1, -1, -1, -1, -1, cName, &
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
       CALL DiagnCont_Find( HcoState%Diagn, &
                            cID, -1, -1, -1, -1, '', &
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
       CALL DiagnCont_PrepareOutput ( HcoState, DgnCont, RC )
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
  SUBROUTINE Diagn_TotalGet( HcoState, Diagn, cName, cID,   COL, &
                             FOUND,     Total, Reset, RC    )
!
! !INPUT PARAMETERS::
!
    TYPE(HCO_State),  POINTER                 :: HcoState       ! HEMCO state obj
    TYPE(DiagnBundle),POINTER                 :: Diagn          ! Diagn bundle obj
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
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DiagnCont),  POINTER  :: DgnCont
    INTEGER                    :: PS
    LOGICAL                    :: FND

    !======================================================================
    ! Diagn_TotalGet begins here!
    !======================================================================

    ! Init
    RC      = HCO_FAIL
    Total   = 0.0_sp
    DgnCont => NULL()
    FND     = .FALSE.
    IF ( PRESENT(FOUND) ) THEN
       FOUND = .FALSE.
    ENDIF

    ! Get collection number
    CALL DiagnCollection_DefineID( Diagn, PS, RC, COL=COL )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If container name is given, search for diagnostics with
    ! the given name.
    IF ( PRESENT( cName ) ) THEN
       CALL DiagnCont_Find( Diagn, -1, -1, -1, -1, -1, cName, &
                            -1, FND, DgnCont, COL=PS )
    ENDIF

    ! If container id is given, search for diagnostics with
    ! the given container ID.
    IF ( PRESENT( cID ) .AND. .NOT. FND ) THEN
       CALL DiagnCont_Find( Diagn, cID, -1, -1, -1, -1, '', &
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
!
! !INPUT PARAMETERS:
!
    TYPE(DiagnCont), POINTER  :: DiagnList   ! List to be removed
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller   - Initialization
!  25 Jan 2016 - R. Yantosca - Bug fix for pgfortran compiler: Test if the
!                              TMPCONT object is associated before deallocating
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Pointers
    TYPE(DiagnCont), POINTER  :: TmpCont
    TYPE(DiagnCont), POINTER  :: NxtCont

    !======================================================================
    ! Diagn_Cleanup begins here!
    !======================================================================

    ! Walk through entire list and remove all containers
    NxtCont => NULL()
    TmpCont => DiagnList
    DO WHILE ( ASSOCIATED( TmpCont ) )

       ! Detach from list
       NxtCont => TmpCont%NextCont

       ! Clean up this container
       CALL DiagnCont_Cleanup( TmpCont )
       IF ( ASSOCIATED( TmpCont ) ) DEALLOCATE ( TmpCont )

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
  FUNCTION Diagn_AutoFillLevelDefined( Diagn, Level, COL ) RESULT ( IsDefined )
!
! !INPUT PARAMETERS:
!
    TYPE(DiagnBundle),POINTER     :: Diagn     ! Diagn bundle obj
    INTEGER, INTENT(IN)           :: Level     ! Level of interest
    INTEGER, INTENT(IN), OPTIONAL :: COL       ! Collection Nr.
!
! !RETURN VALUE:
!
    LOGICAL                       :: IsDefined ! Return argument
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DiagnCollection), POINTER :: ThisColl
    INTEGER                        :: I, RC, PS
    LOGICAL                        :: InUse

    !======================================================================
    ! Diagn_AutoFillLevelDefined begins here!
    !======================================================================

    ! Initialize
    IsDefined = .FALSE.
    ThisColl  => NULL()

    ! Get collection number
    CALL DiagnCollection_DefineID( Diagn, PS, RC, COL=COL, DEF=-1, &
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
  SUBROUTINE DiagnCollection_Get( Diagn,        COL,               &
                                  InUse,        Prefix,            &
                                  nnDiagn,      DeltaYMD,          &
                                  LastYMD,      DeltaHMS, LastHMS, &
                                  OutTimeStamp, RC                  )
!
! !INPUT ARGUMENTS:
!
    INTEGER,          INTENT(IN), OPTIONAL :: COL       ! Collection Nr.
!
! !OUTPUT PARAMETERS:
!
    TYPE(DiagnBundle),POINTER               :: Diagn
    LOGICAL,          INTENT(OUT), OPTIONAL :: InUse
    CHARACTER(LEN=*), INTENT(OUT), OPTIONAL :: Prefix
    INTEGER,          INTENT(OUT), OPTIONAL :: nnDiagn
    INTEGER,          INTENT(OUT), OPTIONAL :: DeltaYMD
    INTEGER,          INTENT(OUT), OPTIONAL :: LastYMD
    INTEGER,          INTENT(OUT), OPTIONAL :: DeltaHMS
    INTEGER,          INTENT(OUT), OPTIONAL :: LastHMS
    INTEGER,          INTENT(OUT), OPTIONAL :: OutTimeStamp
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    TYPE(DiagnCollection), POINTER :: ThisColl
    INTEGER                        :: PS
    LOGICAL                        :: FOUND

    !======================================================================
    ! DiagnCollection_Get begins here!
    !======================================================================

    ! Init
    ThisColl => NULL()
    IF ( PRESENT(Prefix      ) ) Prefix       = ''
    IF ( PRESENT(InUse       ) ) InUse        = .FALSE.
    IF ( PRESENT(nnDiagn     ) ) nnDiagn      = 0
    IF ( PRESENT(DeltaYMD    ) ) DeltaYMD     = 0
    IF ( PRESENT(LastYMD     ) ) LastYMD      = -1
    IF ( PRESENT(DeltaHMS    ) ) DeltaHMS     = 0
    IF ( PRESENT(LastHMS     ) ) LastHMS      = -1
    IF ( PRESENT(OutTimeStamp) ) OutTimeStamp = -1

    ! Get collection number
    CALL DiagnCollection_DefineID( Diagn, PS, RC, COL=COL, InUse=FOUND, ThisColl=ThisColl )
    IF ( RC /= HCO_SUCCESS ) RETURN

    IF ( PRESENT(InUse) ) THEN
       InUse = FOUND
    ENDIF

    ! Get variables from collection
    IF ( FOUND ) THEN
       IF ( PRESENT(Prefix       ) ) Prefix       = ThisColl%PREFIX
       IF ( PRESENT(nnDiagn      ) ) nnDiagn      = ThisColl%nnDiagn
       IF ( PRESENT(DeltaYMD     ) ) DeltaYMD     = ThisColl%DeltaYMD
       IF ( PRESENT(LastYMD      ) ) LastYMD      = ThisColl%LastYMD
       IF ( PRESENT(DeltaHMS     ) ) DeltaHMS     = ThisColl%DeltaHMS
       IF ( PRESENT(LastHMS      ) ) LastHMS      = ThisColl%LastHMS
       IF ( PRESENT(OutTimeStamp ) ) OutTimeStamp = ThisColl%OutTimestamp
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
! !ROUTINE: DiagnCollection_Set
!
! !DESCRIPTION: Subroutine DiagnCollection\_Set sets variables assigned to
! a given diagnostics collection.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCollection_Set( Diagn, COL, InUse, LastYMD, LastHMS, RC )
!
! !INPUT ARGUMENTS:
!
    TYPE(DiagnBundle),POINTER              :: Diagn     ! Diagn bundle
    INTEGER,          INTENT(IN), OPTIONAL :: COL       ! Collection Nr.
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,          INTENT(OUT), OPTIONAL :: InUse
    INTEGER,          INTENT(IN ), OPTIONAL :: LastYMD
    INTEGER,          INTENT(IN ), OPTIONAL :: LastHMS
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    TYPE(DiagnCollection), POINTER :: ThisColl
    INTEGER                        :: PS
    LOGICAL                        :: FOUND

    !======================================================================
    ! DiagnCollection_Set begins here!
    !======================================================================

    ! Init
    ThisColl => NULL()
    IF ( PRESENT(InUse    ) ) InUse     = .FALSE.

    ! Get collection number
    CALL DiagnCollection_DefineID( Diagn, PS, RC, COL=COL, InUse=FOUND, ThisColl=ThisColl )
    IF ( RC /= HCO_SUCCESS ) RETURN

    IF ( PRESENT(InUse) ) THEN
       InUse = FOUND
    ENDIF

    ! Get variables from collection
    IF ( FOUND ) THEN
       IF ( PRESENT(LastYMD ) ) ThisColl%LastYMD = LastYMD
       IF ( PRESENT(LastHMS ) ) ThisColl%LastHMS = LastHMS
    ENDIF

    ! Cleanup
    ThisColl => NULL()

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE DiagnCollection_Set
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
  SUBROUTINE DiagnCont_PrepareOutput( HcoState, DgnCont, RC )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
    USE HCO_Arr_Mod,   ONLY : HCO_ArrAssert
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State),   POINTER       :: HcoState  ! HEMCO state obj
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DiagnCont),   POINTER       :: DgnCont   ! diagnostics container
    INTEGER,           INTENT(INOUT) :: RC        ! Return code
!
! !REVISION HISTORY:
!  19 Dec 2013 - C. Keller: Initialization
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DiagnCollection), POINTER :: ThisColl
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
    ThisColl => NULL()

    !-----------------------------------------------------------------------
    ! Don't do anything for pointer data and/or if data is already in
    ! output format
    !-----------------------------------------------------------------------
    IF ( DgnCont%IsOutFormat ) RETURN
    IF ( DgnCont%DtaIsPtr    ) RETURN

    !-----------------------------------------------------------------------
    ! Get pointer to this collection
    !-----------------------------------------------------------------------
    CALL DiagnCollection_Find( HcoState%Diagn, DgnCont%CollectionID, &
                               FOUND, RC, ThisColl=ThisColl )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! This should never happen
    IF ( .NOT. FOUND .OR. .NOT. ASSOCIATED(ThisColl) ) THEN
       WRITE(MSG,*) 'Diagnostics ', TRIM(DgnCont%cName), ' has invalid ', &
                    'collection ID of ', DgnCont%CollectionID
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
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
       CALL HCO_WARNING( HcoState%Config%Err, MSG, RC, THISLOC=LOC, WARNLEV=1 )
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
       CALL HcoClock_Get( HcoState%Clock, cYYYY=YYYY, cMM=MM, RC=RC )
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
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

    ENDIF

    ! Error trap
    IF ( norm1 <= 0.0_hp ) THEN
       MSG = 'Illegal normalization factor: ' // TRIM(DgnCont%cName)
       CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
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
  SUBROUTINE DiagnCont_Find ( Diagn, cID,      ExtNr, Cat,    Hier,   HcoID, &
                              cName, AutoFill, FOUND, OutCnt, Resume, COL )
!
! !INPUT PARAMETERS:
!
    TYPE(DiagnBundle), POINTEr      :: Diagn    ! diagn bundle
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
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                            :: RC, PS
    TYPE(DiagnCont),       POINTER     :: CurrCnt
    TYPE(DiagnCollection), POINTER     :: ThisColl
    LOGICAL                            :: IsMatch, InUse, Rsm

    !======================================================================
    ! DiagnCont_Find begins here!
    !======================================================================

    ! Initialize
    FOUND    = .FALSE.
    CurrCnt  => NULL()
    ThisColl => NULL()

    ! Get collection number
    CALL DiagnCollection_DefineID( Diagn, PS, RC, COL=COL, Def=-1, &
                      InUse=InUse, OkIfAll=.TRUE., ThisColl=ThisColl )

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
  SUBROUTINE DiagnCont_Link_2D( DgnCont, ThisColl, Trgt2D, RC, HcoState )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !INPUT ARGUMENTS:
!
    TYPE(HCO_STATE),       POINTER, OPTIONAL      :: HcoState    ! HEMCO state obj
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
       IF ( PRESENT(HcoState) ) THEN
          CALL HCO_WARNING( HcoState%Config%Err, MSG, RC, THISLOC=LOC, WARNLEV=2 )
       ELSE
          WRITE(*,*) 'HEMCO WARNING: ', TRIM(MSG)
       ENDIF
       DgnCont%AutoFill = 0
    ENDIF
    IF ( DgnCont%SpaceDim /= 2 ) THEN
       MSG = 'Diagnostics is not 2D: ' // TRIM(DgnCont%cName)
       IF ( PRESENT(HcoState) ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       ELSE
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       ENDIF
       RETURN
    ENDIF

    IF ( SIZE(Trgt2D,1) /= ThisColl%NX .OR. &
         SIZE(Trgt2D,2) /= ThisColl%NY       ) THEN
       MSG = 'Incorrect target array size: ' // TRIM(DgnCont%cName)
       IF ( PRESENT(HcoState) ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       ELSE
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       ENDIF
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
  SUBROUTINE DiagnCont_Link_3D( DgnCont, ThisColl, Trgt3D, RC, HcoState )
!
! !USES:
!
    USE HCO_State_Mod, ONLY : HCO_State
!
! !INPUT PARAEMTERS:
!
    TYPE(HCO_STATE),       POINTER,  OPTIONAL    :: HcoState      ! HEMCO state obj
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
       IF ( PRESENT(HcoState) ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       ELSE
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       ENDIF
       RETURN
    ENDIF
    IF ( DgnCont%AutoFill == 1 ) THEN
       MSG = 'Target diagnostics has autofill flag of 1 - reset to 0: ' &
           // TRIM(DgnCont%cName)
       IF ( PRESENT(HcoState) ) THEN
          CALL HCO_WARNING( HcoState%Config%Err, MSG, RC, THISLOC=LOC, WARNLEV=2 )
       ELSE
          WRITE(*,*) 'HEMCO WARNING: ', TRIM(MSG)
       ENDIF
       DgnCont%AutoFill = 0
    ENDIF
    IF ( DgnCont%SpaceDim /= 3 ) THEN
       MSG = 'Diagnostics is not 3D: ' // TRIM(DgnCont%cName)
       IF ( PRESENT(HcoState) ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       ELSE
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       ENDIF
       RETURN
    ENDIF

    ! Check array size
    IF ( SIZE(Trgt3D,1) /= ThisColl%NX .OR. &
         SIZE(Trgt3D,2) /= ThisColl%NY .OR. &
         SIZE(Trgt3D,3) /= ThisColl%NZ       ) THEN
       MSG = 'Incorrect target array size: ' // TRIM(DgnCont%cName)
       IF ( PRESENT(HcoState) ) THEN
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
       ELSE
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       ENDIF
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
  SUBROUTINE Diagn_Print ( HcoState, Dgn, VerbNr )
!
! !USES:
!
    USE HCO_STATE_MOD,  ONLY : HCO_STATE
!
! !INPUT ARGUMENTS:
!
    TYPE(HCO_STATE),       POINTER    :: HcoState
    TYPE(DiagnCont),       POINTER    :: Dgn
    INTEGER,               INTENT(IN) :: VerbNr
!
! !REVISION HISTORY:
!  01 Aug 2014 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ARGUMENTS:
!
    TYPE(DiagnCollection), POINTER    :: ThisColl
    CHARACTER(LEN=255)                :: MSG
    INTEGER                           :: RC, PS, nx, ny, nz
    REAL(sp)                          :: sm

    ! ================================================================
    ! Diagn_Print begins here
    ! ================================================================

    ! Initialize
    ThisColl => NULL()

    ! Get collection number
    CALL DiagnCollection_DefineID( HcoState%Diagn, PS, RC, &
       COL=Dgn%CollectionID, ThisColl=ThisColl, HcoState=HcoState )
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
    CALL HCO_MSG(HcoState%Config%Err,MSG)

    ! Eventually add details
    IF ( HCO_IsVerb( HcoState%Config%Err, VerbNr ) ) THEN

       ! General information
       WRITE(MSG,*) '   --> Collection         : ', Dgn%CollectionID
       CALL HCO_MSG( HcoState%Config%Err, MSG)
       WRITE(MSG,*) '   --> Diagn ID           : ', Dgn%cID
       CALL HCO_MSG( HcoState%Config%Err, MSG)
       WRITE(MSG,*) '   --> Extension Nr       : ', Dgn%ExtNr
       CALL HCO_MSG( HcoState%Config%Err, MSG)
       WRITE(MSG,*) '   --> Category           : ', Dgn%Cat
       CALL HCO_MSG( HcoState%Config%Err, MSG)
       WRITE(MSG,*) '   --> Hierarchy          : ', Dgn%Hier
       CALL HCO_MSG( HcoState%Config%Err, MSG)
       WRITE(MSG,*) '   --> HEMCO species ID   : ', Dgn%HcoID
       CALL HCO_MSG( HcoState%Config%Err, MSG)
       WRITE(MSG,*) '   --> Autofill?            ', Dgn%AutoFill
       CALL HCO_MSG( HcoState%Config%Err, MSG)
       WRITE(MSG,*) '   --> Space dimension    : ', Dgn%SpaceDim
       CALL HCO_MSG( HcoState%Config%Err, MSG)
       WRITE(MSG,*) '   --> Used level index   : ', Dgn%LevIdx
       CALL HCO_MSG( HcoState%Config%Err, MSG)
       WRITE(MSG,*) '   --> Output unit        : ', TRIM(Dgn%OutUnit)
       CALL HCO_MSG( HcoState%Config%Err, MSG)
       WRITE(MSG,*) '   --> Uniform scaling    : ', Dgn%ScaleFact
       CALL HCO_MSG( HcoState%Config%Err, MSG)
       WRITE(MSG,*) '   --> Current array sum  : ', sm
       CALL HCO_MSG( HcoState%Config%Err, MSG)
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
  SUBROUTINE DiagnCollection_Create ( Diagn,    NX, NY, NZ,    &
                                      TS,       AM2, PREFIX,   &
                                      deltaYMD, deltaHMS, OutTimeStamp, &
                                      RC,       COL, HcoState )
!
! !USES:
!
    USE HCO_STATE_MOD,  ONLY : HCO_STATE
!
! !INPUT ARGUMENTS:
!
    TYPE(DiagnBundle),  POINTER              :: Diagn        ! Diagn bundle
    INTEGER,            INTENT(IN)           :: NX           ! # of lons
    INTEGER,            INTENT(IN)           :: NY           ! # of lats
    INTEGER,            INTENT(IN)           :: NZ           ! # of levels
    REAL(sp),           INTENT(IN)           :: TS           ! timestep [s]
    REAL(hp),           POINTER              :: AM2(:,:)     ! grid box areas [m2]
    CHARACTER(LEN=*),   INTENT(IN)           :: PREFIX       ! Output prefix
    INTEGER,            INTENT(IN), OPTIONAL :: deltaYMD     ! Output frequency
    INTEGER,            INTENT(IN), OPTIONAL :: deltaHMS     ! Output frequency
    INTEGER,            INTENT(IN), OPTIONAL :: OutTimeStamp ! Output time stamp
    TYPE(HCO_State),    POINTER,    OPTIONAL :: HcoState     ! HEMCO state obj
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
!  08 Jan 2015 - C. Keller   - Initial version
!  06 Nov 2015 - C. Keller   - Added OutTimeStamp.
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

    ! Add to collections list. Put at the beginning
    NewCollection%NextCollection => Diagn%Collections
    Diagn%Collections            => NewCollection

    ! Define this collection ID
    Diagn%nnCollections        = Diagn%nnCollections + 1
    NewCollection%CollectionID = Diagn%nnCollections
    COL                        = NewCollection%CollectionID

    ! New output frequency
    IF ( PRESENT(DeltaYMD ) ) NewCollection%DeltaYMD = DeltaYMD
    IF ( PRESENT(DeltaHMS ) ) NewCollection%DeltaHMS = DeltaHMS

    ! Determine output time stamp
    IF ( PRESENT(OutTimeStamp) ) THEN
       ! Make sure it's one of the valid values
       IF ( OutTimeStamp == HcoDiagnStart .OR. &
            OutTimeStamp == HcoDiagnMid   .OR. &
            OutTimeStamp == HcoDiagnEnd         ) THEN
          NewCollection%OutTimeStamp = OutTimeStamp
       ELSE
          WRITE(MSG,*) 'Error when creating diagnostics collection ',    &
             TRIM(NewCollection%PREFIX), ' the specified output time ',  &
             'stamp of ', OutTimeStamp, ' is invalid, must be one of: ', &
             HcoDiagnStart, HcoDiagnMid, HcoDiagnEnd
          IF ( PRESENT(HcoState) ) THEN
             CALL HCO_ERROR(HcoState%Config%Err,MSG,RC,THISLOC=LOC)
          ELSE
             CALL HCO_ERROR(MSG,RC,THISLOC=LOC)
          ENDIF
          RETURN
       ENDIF
    ELSE
       NewCollection%OutTimeStamp = HcoDiagnEnd
    ENDIF

    ! verbose
    IF ( PRESENT(HcoState) ) THEN
       IF ( HCO_IsVerb( HcoState%Config%Err, 1 ) ) THEN
          MSG = 'Created diagnostics collection: '
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(MSG,'(a21,i2)') ' - Collection ID  : ', COL
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(MSG,'(a21,a)' ) ' - PREFIX         : ', TRIM(NewCollection%PREFIX)
          CALL HCO_MSG(HcoState%Config%Err,MSG)
          WRITE(MSG,'(a21,i8,a1,i6)' ) ' - Output interval: ', NewCollection%DeltaYMD, &
                                       ' ',                    NewCollection%DeltaHMS
          CALL HCO_MSG(HcoState%Config%Err,MSG)
       ENDIF
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
  SUBROUTINE DiagnCollection_Cleanup ( Diagn )
!
! !INPUT ARGUMENTS:
!
    TYPE(DiagnBundle), POINTER       :: Diagn
!
! !REVISION HISTORY:
!  08 Jan 2015 - C. Keller - Initial version
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DiagnCollection), POINTER ::  ThisColl => NULL()
    TYPE(DiagnCollection), POINTER ::  NextColl => NULL()

    ! ================================================================
    ! DiagnCollection_Cleanup begins here
    ! ================================================================

    ! Do for every collection in list
    ThisColl => Diagn%Collections
    NextColl => NULL()

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
  SUBROUTINE DiagnCollection_DefineID ( Diagn, PS, RC, COL, DEF, OkIfAll, &
                                        InUse, ThisColl, HcoState )
!
! !USES:
!
    USE HCO_STATE_MOD,     ONLY : HCO_STATE
!
! !INPUT ARGUMENTS:
!
    INTEGER,               INTENT(IN   ), OPTIONAL :: COL        ! desired collection number
    INTEGER,               INTENT(IN   ), OPTIONAL :: DEF        ! default collection number
    LOGICAL,               INTENT(IN   ), OPTIONAL :: OkIfAll    ! Ok if all (PS=-1)
!
! !INPUT/OUTPUT ARGUMENTS:
!
    TYPE(DiagnBundle),     POINTER                 :: Diagn      ! Diagn bundle obj
    TYPE(HCO_STATE),       POINTER,       OPTIONAL :: HcoState   ! HEMCO state obj
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
       PS = Diagn%HcoDiagnIDDefault
    ENDIF
    IF ( PRESENT(COL) ) PS = COL

    ! Check if all collections are selected (-1)
    IF ( PS == -1 ) THEN
       IF ( AllOK ) THEN
          IF ( PRESENT(InUse)    ) InUse    =  .TRUE.
          IF ( PRESENT(ThisColl) ) ThisColl => Diagn%Collections
          RC = HCO_SUCCESS
          RETURN
       ELSE
          WRITE(MSG,*) 'Not allowed to select all collections ', PS
          IF ( PRESENT(HcoState) ) THEN
             CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          ELSE
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          ENDIF
          RETURN
       ENDIF

    ! If individual collection is selected
    ELSE

       ! Try to find collection
       CALL DiagnCollection_Find( Diagn, PS, FOUND, RC, ThisColl=ThisColl )
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Eventually fill argumnet
       IF ( PRESENT(InUse) ) THEN
          InUse = FOUND

       ELSEIF ( .NOT. FOUND ) THEN
          WRITE(MSG,*) 'Diagnostics collection not defined: ', PS
          IF ( PRESENT(HcoState) ) THEN
             CALL HCO_ERROR ( HcoState%Config%Err, MSG, RC, THISLOC=LOC )
          ELSE
             CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          ENDIF
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
  SUBROUTINE DiagnCollection_Find ( Diagn, PS, FOUND, RC, ThisColl )
!
! !INPUT PARAMETERS:
!
    INTEGER,               INTENT(IN   )          :: PS       ! desired collection number
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DiagnBundle),     POINTER                :: Diagn    ! Diagn bundle obj
    LOGICAL,               INTENT(  OUT)          :: FOUND    ! Collection exists?
    INTEGER,               INTENT(INOUT)          :: RC       ! Return code
    TYPE(DiagnCollection), POINTER,      OPTIONAL :: ThisColl ! Pointer to collection
!
! !REVISION HISTORY:
!  01 Apr 2015 - C. Keller   - Initial version
!  10 Jul 2015 - R. Yantosca - Fixed minor issues in ProTeX header
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !ARGUMENTS:
!
    TYPE(DiagnCollection), POINTER :: TmpColl
    CHARACTER(LEN=255), PARAMETER  :: LOC = 'DiagnCollection_Find (hco_diagn_mod.F90)'

    ! ================================================================
    ! DiagnCollection_Find begins here
    ! ================================================================

    ! Init
    TmpColl => NULL()

    ! Check if it's negative
    FOUND = .FALSE.

    ! Loop over all collections
    TmpColl => Diagn%Collections
    DO WHILE ( ASSOCIATED(TmpColl) )

       ! Check if this is the collection of interest
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
! !IROUTINE: DiagnCollection_GetDefaultDelta returns the default diagnostics
!  output intervals based on the 'DiagnFreq' entry of the HEMCO configuration
!  file. This can be one of the following character values: 'Hourly', 'Daily',
!  'Monthly', 'Annually', 'Always', or 'End'; or two integer explicitly denoting
!  the year-month-day and hour-minute-second interval, respectively (format
!  00000000 000000). For example, setting DiagnFreq to '00000000 010000' would
!  be the same as setting it to 'Hourly'.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnCollection_GetDefaultDelta ( HcoState, &
                                               deltaYMD,  deltaHMS, RC )
!
! !USES:
!
    USE HCO_STATE_MOD,   ONLY : HCO_State
    USE HCO_ExtList_Mod, ONLY : GetExtOpt
    USE HCO_ExtList_Mod, ONLY : CoreNr
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_STATE),  POINTER                 :: HcoState  ! HEMCO state obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(  OUT)           :: deltaYMD  ! delta YYYYMMDD
    INTEGER,          INTENT(  OUT)           :: deltaHMS  ! delta HHMMSS
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)           :: RC        ! Return code
!
! !REVISION HISTORY:
!  06 Aug 2015 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL             :: FOUND, SET
    CHARACTER(LEN=255)  :: MSG,   WriteFreq
    CHARACTER(LEN=255)  :: LOC = 'DiagnCollection_GetDefaultDelta (hco_diagn_mod.F90)'

    !=================================================================
    ! DiagnCollection_GetDefaultDelta begins here!
    !=================================================================

    ! Try to get name of diagnostics file

    ! Output frequency. Try to read from configuration file.
    CALL GetExtOpt ( HcoState%Config, CoreNr, 'DiagnFreq', &
                     OptValChar=WriteFreq, FOUND=FOUND, RC=RC )

    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Determine output frequency from given output frequency
    IF ( FOUND ) THEN

       ! Frequency set?
       SET = .FALSE.

       IF ( TRIM(WriteFreq) == 'Hourly' ) THEN
          DeltaYMD = 0
          DeltaHMS = 10000
          SET = .TRUE.
       ELSEIF ( TRIM(WriteFreq) == 'Daily' ) THEN
          DeltaYMD = 1
          DeltaHMS = 0
          SET = .TRUE.
       ELSEIF ( TRIM(WriteFreq) == 'Monthly' ) THEN
          DeltaYMD = 100
          DeltaHMS = 0
          SET = .TRUE.
       ELSEIF ( TRIM(WriteFreq) == 'Annually' ) THEN
          DeltaYMD = 10000
          DeltaHMS = 0
          SET = .TRUE.
       ELSEIF ( TRIM(WriteFreq) == 'Always' ) THEN
          DeltaYMD = 0
          DeltaHMS = 1
          SET = .TRUE.
       ELSEIF ( TRIM(WriteFreq) == 'End' ) THEN
          DeltaYMD = 99999999
          DeltaHMS = 999999
          SET = .TRUE.

       ! If none of the above works, assume that string explicitly gives integer
       ! intervals (YYYYMMDD HHMMSS)
       ELSE
          IF ( LEN(TRIM(WriteFreq)) == 15 ) THEN
             READ(WriteFreq(1 :8 ), * ) DeltaYMD
             READ(WriteFreq(10:15), * ) DeltaHMS
             IF ( DeltaYMD >= 0 .AND. DeltaHMS >= 0 ) SET=.TRUE.
          ENDIF
       ENDIF

       ! Error check
       IF ( .NOT. SET ) THEN
          MSG = 'Cannot define output frequency from string ' // &
                TRIM(WriteFreq) // '. The output frequency must be one of '  // &
                '`Hourly`, `Daily`, `Monthly`, `Annually`, `Always`, `End`,' // &
                ' or the explicit YYYYMMDD HHMMSS interval (15 characters).'
          CALL HCO_ERROR( HcoState%Config%Err, MSG, RC, THISLOC=LOC)
          RETURN
       ENDIF

    ! If output frequency is not explicitly given, set to a default of 1 day.
    ELSE
       DeltaYMD = 1 ! = 00000001 ==> 1 day
       DeltaHMS = 0 ! = 000000
    ENDIF

    ! Force to 'Always' in ESMF environment to make sure that
    ! diagnostics are passed to MAPL HISTORY every time.
#if defined ( ESMF_ )
    DeltaYMD = 0 ! = 00000000
    DeltaHMS = 1 ! = 000001   ==> 1 second!
#endif

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE DiagnCollection_GetDefaultDelta
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Function DiagnCollection_IsTimeToWrite returns true if it is time
!  to write the provided diagnostics collection (identified by the collection
!  number) to output. Whether it is time to write the diagnostics is based upon
!  the current simulation time, the diagnostics output frequency, and the time
!  span since the last output datetime.
!\\
!\\
! !INTERFACE:
!
  FUNCTION DiagnCollection_IsTimeToWrite( HcoState, PS ) &
     RESULT ( TimeToWrite )
!
! !USES:
!
    USE HCO_STATE_MOD,   ONLY : HCO_State
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN   )     :: PS          ! Diagnostics collection
    TYPE(HCO_State), POINTER   :: HcoState    ! HEMCO state obj
!
! !OUTPUT PARAMETERS:
!
    LOGICAL                    :: TimeToWrite ! Is it time to write?
!
! !REVISION HISTORY:
!  06 Aug 2015 - C. Keller   - Initial version
!  30 Sep 2015 - C. Keller   - Bug fix: now set current hour from 0 to 24 to
!                              make sure that it will be greater than previous
!                              hour.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: YYYY, MM,   DD,   h, m, s
    INTEGER             :: lh, lm, ls
    INTEGER             :: delta
    INTEGER             :: dymd, lymd, dhms, lhms
    INTEGER             :: RC
    LOGICAL             :: IsLast
    CHARACTER(LEN=255)  :: LOC = 'DiagnCollection_IsTimeToWrite (hco_diagn_mod.F90)'

    !=================================================================
    ! DiagnCollection_IsTimeToWrite begins here!
    !=================================================================

    ! Init
    TimeToWrite = .FALSE.

    ! Get collection time interval
    CALL DiagnCollection_Get( HcoState%Diagn, PS, DeltaYMD=dymd, &
                   LastYMD=lymd, DeltaHMS=dhms, LastHMS=lhms, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Get current simulation date
    CALL HcoClock_Get( HcoState%Clock, IsLast=IsLast, &
                       sYYYY=YYYY,sMM=MM,sDD=DD,sH=h,sM=m,sS=s,RC=RC)
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Check for last time step
    IF ( IsLast .AND. dymd == 99999999 .AND. dhms == 999999 ) THEN
       TimeToWrite = .TRUE.
    ENDIF

    ! Check if we need to write this collection now
    IF ( .NOT. TimeToWrite .AND. dhms > 0 .AND. lhms >= 0 ) THEN
       ! lh is the last hour of writeout
       lh = FLOOR( MOD(lhms*1.d0, 1000000.0d0 ) / 1.0d4 )
       IF ( h == 0 .AND. lh > 0 ) h = 24
       delta = ( h * 10000 + m * 100 + s ) - lhms
       IF ( delta >= dhms ) TimeToWrite = .TRUE.
    ENDIF

    IF ( .NOT. TimeToWrite .AND. dymd > 0 .AND. lymd >= 0 ) THEN
       delta = ( YYYY * 10000 + MM * 100 + DD ) - lymd
       IF ( delta >= dymd ) TimeToWrite = .TRUE.
    ENDIF

  END FUNCTION DiagnCollection_IsTimeToWrite
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Function DiagnCollection_LastTimesSet returns true if there
!  exists a valid entry for the last datetime that collection PS has been
!  written to disk. This is primarily important to check if the last
!  output date needs be initialized (to non-default values).
!\\
!\\
! !INTERFACE:
!
  FUNCTION DiagnCollection_LastTimesSet( Diagn, PS ) Result ( LastTimesSet )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    TYPE(DiagnBundle), POINTER :: Diagn
    INTEGER, INTENT(IN   )     :: PS   ! Diagnostics collection
!
! !OUTPUT PARAMETERS:
!
    LOGICAL                    :: LastTimesSet ! Are last times defined or not?
!
! !REVISION HISTORY:
!  09 Sep 2015 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER             :: lymd, lhms
    INTEGER             :: RC
    CHARACTER(LEN=255)  :: LOC = 'DiagnCollection_LastTimesSet (hco_diagn_mod.F90)'

    !=================================================================
    ! DiagnCollection_LastTimesSet begins here!
    !=================================================================

    ! Init
    LastTimesSet = .FALSE.

    CALL DiagnCollection_Get( Diagn, PS, LastYMD=lymd, LastHMS=lhms, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Last time stamp is defined if either of the values is greater equal zero.
    IF ( lymd >= 0 .OR. lhms >= 0 ) LastTimesSet = .TRUE.

  END FUNCTION DiagnCollection_LastTimesSet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !DESCRIPTION: Opens a diagnostic configuration file.  This is where you
!  tell HEMCO which diagnostics you would like to send directly to netCDF
!  output.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnFileOpen( HcoConfig, LUN, RC, IsDryRun )
!
! !USES:
!
    USE inquireMod,        ONLY : findFreeLUN
    USE HCO_ExtList_Mod,   ONLY : CoreNr, GetExtOpt
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ), OPTIONAL :: IsDryRun    ! Is it a dry run?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER                 :: HcoConfig   ! HEMCO config obj
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
    ! Scalars
    INTEGER             :: IOS
    LOGICAL             :: EXISTS, FOUND, DoDryRun

    ! Strings
    CHARACTER(LEN=255)  :: MSG, DiagnFile, FileMsg
    CHARACTER(LEN=255)  :: LOC = 'DiagnFileOpen (hco_diagn_mod.F90)'

    !=======================================================================
    ! DiagnFileOpen begins here!
    !=======================================================================

    ! Initialize
    RC  = HCO_SUCCESS
    LUN = -1

    ! Determine if we need to do a dry-run simulation
    IF ( PRESENT( IsDryRun ) ) THEN
       DoDryRun = IsDryRun
    ELSE
       DoDryRun = .FALSE.
    ENDIF

    ! Try to get name of HEMCO diagnostics file
    CALL GetExtOpt( HcoConfig,            CoreNr,     'DiagnFile',           &
                    OptValChar=DiagnFile, FOUND=FOUND, RC=RC                )

    ! Trap potential errors
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Could not find "DiagnFile" in configuration file!'
       CALL HCO_Error( HcoConfig%Err, MSG, RC, LOC )
       RETURN
    ENDIF

    ! If a "DiagnFile" entry is found in the configuration file ...
    IF ( FOUND ) THEN

       ! Test if the diagnostics file exists
       INQUIRE( FILE=TRIM(DiagnFile), EXIST=EXISTS )

       !====================================================================
       ! For dry-runs, print file status and then return
       !====================================================================
       IF ( DoDryRun ) THEN

          ! Test if the file exists and define an output string
          IF ( Exists ) THEN
             FileMsg = 'HEMCO (INIT): Opening'
          ELSE
             FileMsg = 'HEMCO (INIT): REQUIRED FILE NOT FOUND'
          ENDIF

          ! Write message to stdout and then return
          IF ( HcoConfig%amIRoot ) THEN
             WRITE( 6, 300 ) TRIM( FileMsg ), TRIM( DiagnFile )
 300         FORMAT( a, ' ./', a )
          ENDIF
          RETURN
       ENDIF

       !====================================================================
       ! For regular simulations, continue to open the diagnostic file.
       !====================================================================

       ! Find free LUN
       LUN = findFreeLUN()

       ! If the diagnostics file doesn't exist, then exit
       IF ( .NOT. EXISTS ) THEN
          MSG = 'Cannot read file - it does not exist: ' // TRIM(DiagnFile)
          CALL HCO_ERROR( HcoConfig%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Open configuration file
       OPEN ( LUN, FILE=TRIM( DiagnFile ), STATUS='OLD', IOSTAT=IOS )
       IF ( IOS /= 0 ) THEN
          MSG = 'Error opening ' // TRIM(DiagnFile)
          CALL HCO_ERROR( HcoConfig%Err, MSG, RC, THISLOC=LOC )
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
! !IROUTINE: DiagnFileGetNext returns the diagnostics entries of the next
!  line of the diagnostics list file.
!
! !DESCRIPTION: Gets information from the next line of the diagnostic
!  configuration file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnFileGetNext( HcoConfig, LUN,     cName,       &
                               SpcName,   ExtNr,   Cat,   Hier, &
                               SpaceDim,  OutUnit, EOF,   RC,   &
                               lName,     UnitName              )
!
! !USES:
!
    USE HCO_CharTools_Mod
    USE CHARPAK_Mod,       ONLY : STRREPL, STRSPLIT
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   )           :: LUN         ! file LUN
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER                 :: HcoConfig
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
    CHARACTER(LEN=*), INTENT(  OUT), OPTIONAL :: lName
    CHARACTER(LEN=*), INTENT(  OUT), OPTIONAL :: UnitName
!
! !REVISION HISTORY:
!  10 Apr 2015 - C. Keller   - Initial version
!  23 Feb 2016 - C. Keller   - Added lName and UnitName arguments
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
    CALL GetNextLine( LUN, LINE, EOF, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Leave here if end of file
    IF ( .NOT. EOF ) THEN

       ! Parse diagnostics information from line
       CALL STRREPL( LINE, HCO_TAB, HCO_SPC )

       ! Split into substrings
       CALL STRSPLIT( LINE, HCO_SPC, SUBSTR, N )

       ! There must be at least 7 entries
       IF ( N < 7 ) THEN
          MSG = 'Diagnostics entries must have 7 elements: '// TRIM(LINE)
          CALL HCO_ERROR( HcoConfig%Err, MSG, RC, THISLOC=LOC )
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

       ! Eventually get long name
       IF ( PRESENT(lname) ) THEN
          IF ( N > 7 ) THEN
             lName = TRIM(SUBSTR(8))
          ELSE
             lName = TRIM(cName)
          ENDIF
       ENDIF

       ! Eventually get unit name
       IF ( PRESENT(UnitName) ) THEN
          IF ( N > 8 ) THEN
             UnitName = TRIM(SUBSTR(9))
          ELSE
             UnitName = TRIM(OutUnit)
          ENDIF
       ENDIF

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
! !DESCRIPTION: Closes the diagnostic configuration file.
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
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DiagnBundle_Init
!
! !DESCRIPTION: Creates an empty diagnostics bundle
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnBundle_Init ( Diagn )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DiagnBundle), POINTER                :: Diagn
!
! !REVISION HISTORY:
!  17 Feb 2016 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( .NOT. ASSOCIATED(Diagn) ) THEN
       ALLOCATE(Diagn)
       Diagn%Collections => NULL()
       Diagn%HcoDiagnIDDefault = -999
       Diagn%HcoDiagnIDRestart = -999
       Diagn%HcoDiagnIDManual  = -999
       Diagn%nnCollections     = 0
    ENDIF

  END SUBROUTINE DiagnBundle_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DiagnBundle_Cleanup
!
! !DESCRIPTION: Cleans up a diagnostics bundle
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnBundle_Cleanup ( Diagn )
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DiagnBundle), POINTER                :: Diagn
!
! !REVISION HISTORY:
!  17 Feb 2016 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( ASSOCIATED(Diagn) ) THEN
       CALL DiagnCollection_Cleanup ( Diagn )
       DEALLOCATE(Diagn)
       Diagn => NULL()
    ENDIF

  END SUBROUTINE DiagnBundle_Cleanup
!EOC
END MODULE HCO_Diagn_Mod
