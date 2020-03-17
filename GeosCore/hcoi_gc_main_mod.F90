!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoi_gc_main_mod.F90
!
! !DESCRIPTION: Module hcoi\_gc\_main\_mod.F90 is the HEMCO-to-GEOS-Chem
! interface module, providing the link between GEOS-Chem and HEMCO. It
! contains wrapper routines to initialize, execute and finalize HEMCO from
! within GEOS-Chem. These routines are called from emissions\_mod.F90.
!\\
!\\
! Notes:
! \begin{itemize}
! \item HEMCO is used to calculate all emission fields. The emission tendencies
!  are passed to GEOS-Chem in module mixing\_mod.F90.
! \item Most meteorological fields needed by the HEMCO extensions are provided
!  through the GEOS-Chem meteorological state object Met\_State. Few fields
!  such as the pressure edges or J-values are defined and updated explicitly
!  within this module.
! \End{itemize}
! !INTERFACE:
!
MODULE HCOI_GC_Main_Mod
!
! !USES:
!
  USE Precision_Mod
  USE HCO_Error_Mod
  USE HCO_Interface_Mod
  USE HCOX_State_Mod, ONLY : Ext_State
  USE HCO_State_Mod,  ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOI_GC_Init
  PUBLIC  :: HCOI_GC_Run
  PUBLIC  :: HCOI_GC_Final
  PUBLIC  :: HCOI_GC_WriteDiagn
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Calc_SumCosZa
  PRIVATE :: ExtState_InitTargets
  PRIVATE :: ExtState_SetFields
  PRIVATE :: ExtState_UpdateFields
  PRIVATE :: Get_SzaFact
  PRIVATE :: GridEdge_Set
  PRIVATE :: CheckSettings
  PRIVATE :: SetHcoGrid
  PRIVATE :: SetHcoSpecies
#if !defined(ESMF_) && !defined( MODEL_WRF )
  !=========================================================================
  ! These are only needed for GEOS-Chem "Classic"
  !=========================================================================
  PRIVATE :: Get_GC_Restart
  PRIVATE :: Get_Met_Fields
  PRIVATE :: Get_Boundary_Conditions
#endif
!
! !REMARKS:
!  This module is ignored if you are using HEMCO in an ESMF environment.
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller   - Initial version.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  !------------------------------------
  ! %%% Species ID's %%%
  !------------------------------------
  INTEGER             :: id_HNO3
  INTEGER             :: id_LIMO
  INTEGER             :: id_NO
  INTEGER             :: id_NO2
  INTEGER             :: id_O3
  INTEGER             :: id_POPG

  !------------------------------------
  ! %%% Arrays, Pointers, Targets %%%
  !------------------------------------

  ! Internal met fields (will be used by some extensions)
  INTEGER,  TARGET    :: HCO_PBL_MAX                      ! level
  REAL(hp), POINTER   :: HCO_FRAC_OF_PBL(:,:,:)
  REAL(hp), POINTER   :: HCO_SZAFACT(:,:)

  ! Arrays to store J-values (used by Paranox extension)
  REAL(hp), POINTER   :: JNO2(:,:)
  REAL(hp), POINTER   :: JOH(:,:)

  ! Sigma coordinate (temporary)
  REAL(hp), POINTER   :: ZSIGMA(:,:,:)

  ! Sum of cosine of the solar zenith angle. Used to impose a
  ! diurnal variability on OH concentrations
  REAL(fp), POINTER   :: SUMCOSZA(:,:)
!
! !DEFINED PARAMETERS:
!
  ! Temporary toggle for diagnostics
  LOGICAL,  PARAMETER :: DoDiagn = .TRUE.

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_GC_Init
!
! !DESCRIPTION: Subroutine HCOI\_GC\_INIT initializes the HEMCO derived
! types and arrays. The HEMCO configuration is read from the HEMCO
! configuration file (as listed in Input\_Opt%HcoConfigFile) and stored in
! the HEMCO configuration object. The entire HEMCO setup is based upon the
! entries in the HEMCO configuration object. It is possible to explicitly
! provide a (previously read) HEMCO configuration object via input argument
! `HcoConfig`. In this case the HEMCO configuration file will not be read
! any more.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_GC_Init( Input_Opt, State_Chm, State_Grid, &
                           State_Met, RC,        HcoConfig )
!
! !USES:
!
    USE CMN_SIZE_Mod,       ONLY : NDSTBIN
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_EMIS, GET_TS_DYN
    USE TIME_MOD,           ONLY : GET_TS_CHEM
#ifdef TOMAS
    USE TOMAS_MOD,          ONLY : IBINS
    USE TOMAS_MOD,          ONLY : Xk
#endif

    ! HEMCO routines
    USE HCO_Types_Mod,      ONLY : ConfigObj
    USE HCO_Config_Mod,     ONLY : Config_ReadFile, ConfigInit
    USE HCO_State_Mod,      ONLY : HcoState_Init
    USE HCO_Diagn_Mod,      ONLY : DiagnFileOpen
    USE HCO_Driver_Mod,     ONLY : HCO_Init
    USE HCOI_GC_Diagn_Mod,  ONLY : HCOI_GC_Diagn_Init
    USE HCOX_Driver_Mod,    ONLY : HCOX_Init
    USE HCOX_State_Mod,     ONLY : ExtStateInit
!
! !INPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(IN   )          :: State_Chm  ! Chemistry state
    TYPE(GrdState),   INTENT(IN   )          :: State_Grid ! Grid state
    TYPE(MetState),   INTENT(IN   )          :: State_Met  ! Met state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)          :: Input_Opt  ! Input opts
    TYPE(ConfigObj),  POINTER,      OPTIONAL :: HcoConfig  ! HEMCO config object
    INTEGER,          INTENT(INOUT)          :: RC         ! Failure or success
!
! !REVISION HISTORY:
!  12 Sep 2013 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                   :: LSTRAT,  FOUND
    INTEGER                   :: nHcoSpc, HMRC
    INTEGER                   :: N

    ! Strings
    CHARACTER(LEN=255)        :: OptName, ThisLoc, Instr
    CHARACTER(LEN=512)        :: ErrMSg

    ! Pointers
    TYPE(ConfigObj), POINTER  :: iHcoConfig => NULL()
    TYPE(Species),   POINTER  :: SpcInfo

    !=======================================================================
    ! HCOI_GC_INIT begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at HCOI_GC_Init (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    ! Define all species ID's here, for use in module routines below
    id_HNO3  = Ind_('HNO3')
    id_LIMO  = Ind_('LIMO')
    id_NO    = Ind_('NO'  )
    id_NO2   = Ind_('NO2' )
    id_O3    = Ind_('O3'  )
    id_POPG  = Ind_('POPG')

    ! Create a splash page
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(a)' ) REPEAT( '%', 79 )
       WRITE( 6, 100   ) 'HEMCO: Harvard-NASA Emissions Component'
       WRITE( 6, 101   ) 'You are using HEMCO version ', ADJUSTL(HCO_VERSION)
       WRITE( 6, '(a)' ) REPEAT( '%', 79 )
 100   FORMAT( '%%%%%', 15x, a,      15x, '%%%%%' )
 101   FORMAT( '%%%%%', 15x, a, a12, 14x  '%%%%%' )
    ENDIF

    !=======================================================================
    ! Read HEMCO configuration file and save into buffer. This also
    ! sets the HEMCO error properties (verbose mode? log file name,
    ! etc.) based upon the specifications in the configuration file.
    ! The log file is now read in two phases: phase 1 reads only the
    ! settings and extensions; phase 2 reads all data fields. This
    ! way, settings and extension options can be updated before
    ! reading all the associated fields. For instance, if the LEMIS
    ! toggle is set to false (=no emissions), all extensions can be
    ! deactivated. Similarly, certain brackets can be set explicitly
    ! to make sure that these data is only read by HEMCO if the
    ! corresponding GEOS-Chem switches are turned on.
    ! (ckeller, 2/13/15).
    !=======================================================================

    ! If HcoConfig is provided
    IF ( PRESENT( HcoConfig ) ) iHcoConfig => HcoConfig

    !---------------------------------------
    ! Initialize HEMCO config object
    !---------------------------------------
    CALL ConfigInit ( iHcoConfig, HMRC, State_Chm%nSpecies )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ConfigInit"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
    ENDIF

    ! Is this the root CPU?
    iHcoConfig%amIRoot = Input_Opt%amIRoot

    ! Met and grid parameters
    iHcoConfig%MetField = Input_Opt%MetField
    iHcoConfig%GridRes  = State_Grid%GridRes

    ! Pass GEOS-Chem species information to HEMCO config object to
    ! facilitate reading GEOS-Chem restart file via HEMCO
    iHcoConfig%nModelSpc = State_Chm%nSpecies
    iHcoConfig%nModelAdv = State_Chm%nAdvect
    DO N = 1, State_Chm%nSpecies
       ! Get info for this species from the species database
       SpcInfo => State_Chm%SpcData(N)%Info

       ! Model ID and species name
       iHcoConfig%ModelSpc(N)%ModID      = SpcInfo%ModelID
       iHcoConfig%ModelSpc(N)%SpcName    = TRIM( SpcInfo%Name )
    ENDDO

    !---------------------------------------
    ! Phase 1: read settings and switches
    !---------------------------------------
    CALL Config_ReadFile( Input_Opt%amIRoot,        &
                          iHcoConfig,               &
                          Input_Opt%HcoConfigFile,  &
                          1,                        &
                          HMRC,                     &
                          IsDryRun=Input_Opt%DryRun )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "Config_Readfile" (Phase 1)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
    ENDIF

    ! Check settings
    CALL CheckSettings( iHcoConfig, Input_Opt, State_Met, State_Chm,  HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "CheckSettings"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    !---------------------------------------
    ! Open logfile
    !---------------------------------------
    IF ( Input_Opt%amIRoot ) THEN
       CALL HCO_LOGFILE_OPEN( iHcoConfig%Err, RC=HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCO_LogFile_Open"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
    ENDIF

    !---------------------------------------
    ! Phase 2: read fields
    !---------------------------------------
    CALL Config_ReadFile( Input_Opt%amIRoot,        &
                          iHcoConfig,               &
                          Input_Opt%HcoConfigFile,  &
                          2,                        &
                          HMRC,                     &
                          IsDryRun=Input_Opt%DryRun )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "Config_Readfile" (Phase 2)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    !=======================================================================
    ! Initialize HEMCO state object and populate it
    !=======================================================================

    !-----------------------------------------------------------------------
    ! Extract species to use in HEMCO. nHcoSpc denotes the number of
    ! species that shall be used in HEMCO. The species properties are
    ! defined in the Register_Species call below.
    ! Typically, nHcoSpc is just the number of species defined in both
    ! the HEMCO configuration file and GEOS-Chem. However, additional
    ! species can be defined, e.g. those not transported in GEOS-Chem
    ! (e.g. SESQ) or tagged species (e.g. specialty simulations).
    !-----------------------------------------------------------------------
    CALL SetHcoSpecies ( Input_Opt, State_Chm, HcoState,  nHcoSpc, 1, HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = &
         'Error encountered in ""SetHcoSpecies" (first call, to get species)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Now that number of HEMCO species are known, initialize HEMCO
    ! state object.  Links the HEMCO configuration file object
    ! iHcoConfig to HcoState%Config.
    !-----------------------------------------------------------------------
    CALL HcoState_Init( HcoState, iHcoConfig, nHcoSpc, HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCOState_Init"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Register species. This will define all species properties
    ! (names, molecular weights, etc.) of the HEMCO species.
    !-----------------------------------------------------------------------
    CALL SetHcoSpecies ( Input_Opt, State_Chm, HcoState, nHcoSpc, 2, HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = &
     'Error encountered in "SetHcoSpecies" (second call, to register species)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Set the HEMCO grid
    !-----------------------------------------------------------------------
    CALL SetHcoGrid( State_Grid, State_Met, HcoState, RC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "Set_Grid"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !=======================================================================
    ! Set misc. parameter
    !=======================================================================

    ! Emission, chemistry and dynamics timestep in seconds
    HcoState%TS_EMIS = GET_TS_EMIS()
    HcoState%TS_CHEM = GET_TS_CHEM()
    HcoState%TS_DYN  = GET_TS_DYN()

    ! Is this an ESMF simulation or not?
    ! The ESMF flag must be set before calling HCO_Init because the
    ! source file name is set differently in an ESMF environment
    ! compared to a stand-alone version: in ESMF, the source file name
    ! is set to the container name since this is the identifying name
    ! used by ExtData.
#ifdef ESMF_
    HcoState%Options%isESMF = .TRUE.
#else
    HcoState%Options%isESMF = .FALSE.
#endif

    ! Set deposition length scale. This determines if dry deposition
    ! frequencies are calculated over the entire PBL or the first
    ! model layer only.
    HcoState%Options%PBL_DRYDEP = Input_Opt%PBL_DRYDEP

    !----------------------------------------------------------------------
    ! Are we running HEMCO in a dry-run mode?
    ! This is dictated by the GEOS-Chem environment. If GEOS-Chem is in a
    ! dry-run mode, no compute is performed and files are only "checked".
    ! Simulations will NOT stop on missing files. This is intended to be a
    ! quick sanity check to make sure that GEOS-Chem IO are all correctly
    ! set up, which is why most of the runs fail to complete successfully.
    ! (hplin, 11/2/19)
    !
    ! Dry run simulations will send output to the stdout (which usually
    ! gets piped to the GEOS-Chem log file).
    !
    ! NOTE: The dry-run option is only invoked in GEOS-Chem "Classic",
    ! so these values will remain at their defaults (.FALSE. and -1,
    ! respectively) when we use HEMCO in external ESMs (bmy, 11/13/19)
    !----------------------------------------------------------------------
    HcoState%Options%IsDryRun = Input_Opt%DryRun

    !=======================================================================
    ! Initialize HEMCO internal lists and variables. All data
    ! information is written into internal lists (ReadList) and
    ! the HEMCO configuration file is removed from buffer in this
    ! step. This also initializes the HEMCO clock as well as the
    ! HEMCO emissions diagnostics collection.
    !=======================================================================
    CALL HCO_Init( HcoState, HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCO_Init"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    ! Save # of defined dust species in HcoState
    HcoState%nDust                     =  NDSTBIN

    ! Use marine organic aerosols?
    HcoState%MarinePOA                 =  Input_Opt%LMPOA

#ifdef TOMAS

    ! Save # of TOMAS size bins in HcoState
    HcoState%MicroPhys%nBins           =  IBINS

    ! Point to TOMAS bin boundaries array (Xk) in HcoState
    HcoState%MicroPhys%BinBound        => Xk

    ! Save # of TOMAS active mode bins in HcoState
#if defined( TOMAS40 )
    HcoState%MicroPhys%nActiveModeBins =  10
#elif defined( TOMAS15 )
    HcoState%MicroPhys%nActiveModeBins =  3
#else
    HcoState%MicroPhys%nActiveModeBins =  0
#endif
#endif

    !=======================================================================
    ! Initialize all HEMCO extensions. This also selects the required
    ! met fields used by each extension.
    !=======================================================================
    CALL HCOX_Init( HcoState, ExtState, HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCOX_Init"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !=======================================================================
    ! For dry-runs only: Print the status of the HEMCO diagnostic
    ! configuration file (e.g. HEMCO_Diagn.rc), and then exit
    !=======================================================================
    IF ( Input_Opt%DryRun ) THEN
       CALL DiagnFileOpen( HcoState%Config, N, RC, IsDryRun=.TRUE. )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Update and check logical switches in Input_Opt
    !-----------------------------------------------------------------------

    ! Soil NOx
    Input_Opt%LSOILNOX      = ( ExtState%SoilNOx > 0 )

    ! Ginoux dust emissions
    IF ( ExtState%DustGinoux > 0 ) THEN
       IF ( .not. Input_Opt%LDUST ) THEN
          ErrMsg = 'DustGinoux is on in HEMCO but LDUST=F in input.geos!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
       Input_Opt%LDEAD = .FALSE.
    ENDIF

    ! DEAD dust emissions
    IF ( ExtState%DustDead > 0 ) THEN
       IF ( .not. Input_Opt%LDUST ) THEN
          ErrMsg = 'DustDead is on in HEMCO but LDUST=F in input.geos!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
       Input_Opt%LDEAD = .TRUE.
    ENDIF

    ! Dust alkalinity
    IF ( ExtState%DustAlk > 0 ) THEN
       IF ( .not. Input_Opt%LDSTUP ) THEN
          ErrMsg = 'DustAlk is on in HEMCO but LDSTUP=F in input.geos'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------------
    ! Set constants for POPs simulation
    !-----------------------------------------------------------------------
    IF ( ExtState%GC_POPs > 0 ) THEN
       ExtState%POP_DEL_H   = Input_Opt%POP_DEL_H
       ExtState%POP_KOA     = Input_Opt%POP_KOA
       ExtState%POP_KBC     = Input_Opt%POP_KBC
       ExtState%POP_DEL_Hw  = Input_Opt%POP_DEL_Hw
       ExtState%POP_XMW     = Input_Opt%POP_XMW
       ExtState%POP_HSTAR   = Input_Opt%POP_HSTAR
    ENDIF

    !-----------------------------------------------------------------------
    ! Initialize ExtState target arrays.
    ! Extensions typically depend on environmental dependent met.
    ! variables such as wind speed, surface temp., etc. Pointers
    ! to these (2D or 3D) fields are defined in the extension object.
    ! Here, we need to make sure that these pointers are properly
    ! connected.
    !-----------------------------------------------------------------------
    CALL ExtState_InitTargets( HcoState, ExtState, State_Grid, HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtState_InitTargets"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Define diagnostics
    !-----------------------------------------------------------------------
    IF ( DoDiagn ) THEN

       ! Set up traditional GEOS-Chem NDxx diagnostics for emissions
       CALL HCOI_GC_Diagn_Init( Input_Opt, HcoState, ExtState, HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCOI_GC_Diagn_Init"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    !=================================================================
    ! Cleanup and quit
    !=================================================================

    ! Eventually remove pointer
    IF ( PRESENT( HcoConfig ) ) iHcoConfig => NULL()

    ! Leave w/ success
    RC = GC_SUCCESS

    END SUBROUTINE HCOI_GC_INIT
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_GC_Run
!
! !DESCRIPTION: Subroutine HCOI\_GC\_Run runs HEMCO from GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_GC_Run( Input_Opt, State_Chm, State_Grid,  &
                          State_Met, EmisTime,  Phase,  RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Get_Ndep_Mod,    ONLY : Reset_Dep_N   ! For soilnox
    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Met_Mod,   ONLY : MetState
    USE Time_Mod

    ! HEMCO routines
    USE HCO_Clock_Mod,   ONLY : HcoClock_Get
    USE HCO_Clock_Mod,   ONLY : HcoClock_EmissionsDone
    USE HCO_Diagn_Mod,   ONLY : HcoDiagn_AutoUpdate
    USE HCO_FluxArr_Mod, ONLY : HCO_FluxarrReset
    USE HCO_Driver_Mod,  ONLY : HCO_Run
    USE HCOX_Driver_Mod, ONLY : HCOX_Run
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: EmisTime   ! Is this an emission time step?
    INTEGER,          INTENT(IN   )  :: Phase      ! Run phase (see remarks)
    TYPE(GrdState),   INTENT(IN   )  :: State_Grid ! Grid state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input options
    TYPE(MetState),   INTENT(INOUT)  :: State_Met  ! Meteo state
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm  ! Chemistry state
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Phase -1 : Used for GCHP
!  Phase  0 : Simplified Phase 1 for reading initial met fields and restart file
!  Phase  1 : Update HEMCO clock and HEMCO data list and get met fields
!  Phase  2 : Perform emissions calculation
!
! !REVISION HISTORY:
!  12 Sep 2013 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd variables
    REAL(f8), SAVE     :: PrevTAU = -999999.0_f8

    ! Scalars
    LOGICAL            :: notDryRun
    INTEGER            :: HMRC
    LOGICAL            :: IsEmisTime
    LOGICAL            :: IsEndStep

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc, Instr
    CHARACTER(LEN=512) :: ErrMsg

    ! Arrays
    INTEGER            :: D(2)               ! Variable for date and time

    !=======================================================================
    ! HCOI_GC_RUN begins here!
    !=======================================================================

    ! Initialize
    RC        = GC_SUCCESS
    HMRC      = HCO_SUCCESS
    ErrMsg    = ''
    ThisLoc   = &
       ' -> at HCOI_GC_Run (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr     = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '      // &
                'HEMCO log file for additional error messages!'
    notDryRun = ( .not. Input_Opt%DryRun )

    !=======================================================================
    ! Make sure HEMCO time is in sync with simulation time
    ! This is now done in main.F90
    !=======================================================================
    CALL SetHcoTime( EmisTime, HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "SetHcoTime"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !=======================================================================
    ! See if it's time for emissions. Don't just use the EmisTime flag in
    ! case that we call this routine multiple times. IsEmisTime will only
    ! be true if this is an emission time step AND emissions have not yet
    ! been calculated for that time step.
    !=======================================================================
    CALL HcoClock_Get( HcoState%Clock, IsEmisTime=IsEmisTime, RC=HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HcoClock_Get"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    ! Check if this is the last GEOS-Chem timestep
    IF ((( HcoState%Clock%ThisYear * 10000 + HcoState%Clock%ThisMonth * 100 + &
           HcoState%Clock%ThisDay ) == Input_Opt%NYMDe ) .AND. &
        (( HcoState%Clock%ThisHour * 10000 + HcoState%Clock%ThisMin   * 100 + &
           HcoState%Clock%ThisSec ) == Input_Opt%NHMSe )) THEN
       IsEndStep = .TRUE.
    ELSE
       IsEndStep = .FALSE.
    ENDIF

    !=======================================================================
    ! Reset all emission and deposition values. Do this only if it is time
    ! for emissions, i.e. if those values will be refilled.
    !=======================================================================
    IF ( IsEmisTime .AND. Phase == 2 .and. notDryRun ) THEN
       CALL HCO_FluxArrReset( HcoState, HMRC                                )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCO_FluxArrReset"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Define pressure edges [Pa] on HEMCO grid.
    ! At Phase 0, the pressure field is not known yet.
    !=======================================================================
    IF ( Phase /= 0 .and. notDryRun ) THEN
       CALL GridEdge_Set( State_Grid, State_Met, HcoState, HMRC  )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "GridEdge_Set"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

#if !defined( ESMF_ )
    ! Check if HEMCO has already been called for this timestep
    IF ( ( Phase == 1 ) .and. ( GET_TAU() == PrevTAU ) .and. Input_Opt%amIRoot ) THEN
       Print*, 'HEMCO already called for this timestep. Returning.'
       RETURN
    ENDIF
#endif

    !=======================================================================
    ! Set HCO options
    !=======================================================================

    ! Range of species and emission categories.
    ! Set Extension number ExtNr to 0, indicating that the core
    ! module shall be executed.
    HcoState%Options%SpcMin = 1
    HcoState%Options%SpcMax = -1
    HcoState%Options%CatMin = 1
    HcoState%Options%CatMax = -1
    HcoState%Options%ExtNr  = 0

    ! Use temporary array?
    HcoState%Options%FillBuffer = .FALSE.

    !=======================================================================
    ! Run HCO core module
    ! Pass phase as argument. Phase 1 will update the emissions list,
    ! phase 2 will calculate the emissions. Emissions will be written into
    ! the corresponding flux arrays in HcoState.
    !=======================================================================
    CALL HCO_Run( HcoState, Phase, HMRC, IsEndStep )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCO_Run"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

#if !defined(ESMF_) && !defined( MODEL_WRF )
    !=======================================================================
    ! Get met fields from HEMCO (GEOS-Chem "Classic" only)
    !=======================================================================
    IF ( ( Phase == 0 .or. PHASE == 1 ) .and. notDryRun ) THEN
       CALL Get_Met_Fields( Input_Opt, State_Chm, State_Grid, State_Met, &
                            Phase,     RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Get_Met_Fields"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Get fields from GEOS-Chem restart file (GEOS-Chem "Classic" only)
    !=======================================================================
    IF ( Phase == 0 .and. notDryRun ) THEN
       CALL Get_GC_Restart( Input_Opt, State_Chm, State_Grid, State_Met, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Get_GC_Restart"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Get boundary conditions from HEMCO (GEOS-Chem "Classic" only)
    !=======================================================================

    ! Assume BCs are 3-hourly and only get from HEMCO when needed
    IF ( PHASE == 0 ) THEN
       D = GET_FIRST_BC_TIME()
    ELSE
       D = GET_BC_TIME()
    ENDIF
    IF ( State_Grid%NestedGrid .and. notDryRun .and. &
       ( Phase == 0 .or. ( PHASE == 1 .and. ITS_TIME_FOR_BC() ) ) ) THEN
       IF ( Input_Opt%LTRAN ) THEN
          CALL Get_Boundary_Conditions( Input_Opt, State_Chm, State_Grid, &
                                        State_Met, D(1), D(2), RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Get_Boundary_Conditions"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF
    ENDIF

#endif

    !=======================================================================
    ! Do the following only if it's time to calculate emissions
    !=======================================================================
    IF ( Phase == 2 .AND. IsEmisTime ) THEN

       !--------------------------------------------------------------------
       ! Set / update ExtState fields.
       ! Extensions typically depend on environmental dependent met.
       ! variables such as wind speed, surface temp., etc. Pointers
       ! to these (2D or 3D) fields are defined in the extension object.
       ! Here, we need to make sure that these pointers are properly
       ! connected.
       !--------------------------------------------------------------------
       IF ( notDryRun ) THEN
          CALL ExtState_SetFields( State_Chm, State_Met, &
                                   HcoState,  ExtState,  HMRC )

          ! Trap potential errors
          IF ( HMRC /= HCO_SUCCESS ) THEN
             RC     = HMRC
             ErrMsg = 'Error encountered in "ExtState_SetFields"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
             CALL Flush( HcoState%Config%Err%Lun )
             RETURN
          ENDIF

          CALL ExtState_UpdateFields( Input_Opt, State_Chm,            &
                                      State_Grid, State_Met, HcoState, &
                                      ExtState,   HMRC )

          ! Trap potential errors
          IF ( HMRC /= HCO_SUCCESS ) THEN
             RC     = HMRC
             ErrMsg = 'Error encountered in "ExtState_UpdateFields"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
             CALL Flush( HcoState%Config%Err%Lun )
             RETURN
          ENDIF
       ENDIF

       !====================================================================
       ! Run HCO extensions. Emissions will be added to corresponding
       ! flux arrays in HcoState.
       !====================================================================
       CALL HCOX_Run( HcoState, ExtState, HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCOX_Run"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF

       !====================================================================
       ! Update all 'AutoFill' diagnostics. This makes sure that all
       ! diagnostics fields with the 'AutoFill' flag are up-to-date. The
       ! AutoFill flag is specified when creating a diagnostics container
       ! (Diagn_Create).
       !====================================================================
       IF ( DoDiagn .and. notDryRun ) THEN
          CALL HcoDiagn_AutoUpdate( HcoState, HMRC )

          ! Trap potential errors
          IF ( HMRC /= HCO_SUCCESS ) THEN
             RC     = HMRC
             ErrMsg = 'Error encountered in "HcoDiagn_AutoUpdate"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
             CALL Flush( HcoState%Config%Err%Lun )
             RETURN
          ENDIF
       ENDIF

       !====================================================================
       ! Reset the accumulated nitrogen dry and wet deposition to zero.
       ! Will be re-filled in drydep and wetdep.
       !====================================================================
       IF ( ( Input_Opt%ITS_A_FULLCHEM_SIM   .or.                           &
              Input_Opt%ITS_AN_AEROSOL_SIM ) .and. notDryRun ) THEN
          CALL RESET_DEP_N( State_Chm )
       ENDIF

       !====================================================================
       ! Emissions are now done for this time step
       !====================================================================
       CALL HcoClock_EmissionsDone( HcoState%Clock, HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HcoClock_EmissionsDone"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF

    ENDIF

    ! Save TAU
    PrevTAU = GET_TAU()

    ! We are now back in GEOS-Chem environment, hence set
    ! return flag accordingly!
    RC = GC_SUCCESS

  END SUBROUTINE HCOI_GC_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_GC_Final
!
! !DESCRIPTION: Subroutine HCOI\_GC\_Final cleans up HEMCO. This routine
! should be called before the finalize routines of State\_Chm in order to
! make sure that the emissions flux pointers are properly removed!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_GC_Final( Error, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Driver_Mod,   ONLY : HCO_Final
    USE HCO_Diagn_Mod,    ONLY : DiagnBundle_Cleanup
    USE HCO_State_Mod,    ONLY : HcoState_Final
    USE HCOX_Driver_Mod,  ONLY : HCOX_Final
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)  :: Error       ! Cleanup after exit?
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  12 Sep 2013 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: HMRC

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc, Instr
    CHARACTER(LEN=512) :: ErrMsg

    !=======================================================================
    ! HCOI_GC_Final begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at HCOI_GC_Final (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! Cleanup HEMCO core
    !-----------------------------------------------------------------------
    CALL HCO_Final( HcoState, Error, HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCO_Final"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Cleanup extensions and ExtState object
    ! This will also nullify all pointer to the met fields.
    !-----------------------------------------------------------------------
    CALL HCOX_Final( HcoState, ExtState, HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCOX_Final"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Cleanup diagnostics
    !-----------------------------------------------------------------------
    CALL DiagnBundle_Cleanup( HcoState%Diagn )

    !-----------------------------------------------------------------------
    ! Cleanup HcoState object
    !-----------------------------------------------------------------------
    CALL HcoState_Final( HcoState )

    !-----------------------------------------------------------------------
    ! Deallocate module variables
    !-----------------------------------------------------------------------
    IF ( ASSOCIATED( ZSIGMA ) ) THEN
       DEALLOCATE( ZSIGMA, STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:ZSIGMA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( HCO_FRAC_OF_PBL ) ) THEN
       DEALLOCATE( HCO_FRAC_OF_PBL, STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:HCO_FRAC_OF_PBL', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( HCO_SZAFACT ) ) THEN
       DEALLOCATE( HCO_SZAFACT, STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:HCO_SZAFACT', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( JNO2 ) ) THEN
       DEALLOCATE( JNO2, STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:JNO2', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( JOH ) ) THEN
       DEALLOCATE( JOH, STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:JOH', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( SUMCOSZA ) ) THEN
       DEALLOCATE( SUMCOSZA, STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:SUMCOSZA', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE HCOI_GC_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_GC_WriteDiagn
!
! !DESCRIPTION: Subroutine HCOI\_GC\_WriteDiagn is the wrapper routine to
! write the HEMCO diagnostics. This will only write the diagnostics of
! diagnostics collection 1.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_GC_WriteDiagn( Input_Opt, Restart, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCOIO_Diagn_Mod, ONLY : HcoDiagn_Write
    USE Input_Opt_Mod,   ONLY : OptInput
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN   ) :: Input_Opt    ! Input options
    LOGICAL,        INTENT(IN   ) :: Restart      ! write restart (enforced)?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  01 Apr 2015 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: HMRC

    ! Strings
    CHARACTER(LEN=255)  :: ThisLoc, Instr
    CHARACTER(LEN=512)  :: ErrMsg

    !=======================================================================
    ! HCOI_GC_WriteDiagn begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at HCOI_GC_WriteDiagn (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! Make sure HEMCO time is in sync
    !-----------------------------------------------------------------------
    CALL SetHcoTime ( .FALSE., HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "SetHcoTime"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Write diagnostics
    !-----------------------------------------------------------------------
    CALL HcoDiagn_Write( HcoState, RESTART, HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HcoDiagn_Write"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

  END SUBROUTINE HCOI_GC_WriteDiagn
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtState_InitTargets
!
! !DESCRIPTION: SUBROUTINE ExtState\_InitTargets allocates some local arrays
! that act as targets for the ExtState object.
!\\
! Note that for now, this explicitly assumes that the HEMCO emissions grid is
! the same as the GEOS-Chem simulation grid. To support other grids, the met
! data has to be regridded explicitly at every time step!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtState_InitTargets( HcoState, ExtState, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Arr_Mod,    ONLY : HCO_ArrAssert
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState),   INTENT(IN   )  :: State_Grid ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_STATE),  POINTER        :: HcoState
    TYPE(EXT_STATE),  POINTER        :: ExtState
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller    - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: HMRC

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc, Instr
    CHARACTER(LEN=512) :: ErrMsg

    !=======================================================================
    ! ExtState_InitTargets begins here
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at ExtState_InitTargets (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! HCO_SZAFACT is not defined in Met_State.  Hence need to
    ! define here so that we can point to them.
    !
    ! Now include HCO_FRAC_OF_PBL and HCO_PBL_MAX for POPs specialty
    ! simulation (mps, 8/20/14)
    ! ----------------------------------------------------------------------
    IF ( ExtState%SZAFACT%DoUse ) THEN

       ALLOCATE( SUMCOSZA( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:SUMCOSZA', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       SUMCOSZA = 0.0_fp

       ALLOCATE( HCO_SZAFACT( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:HCO_SZAFACT', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       HCO_SZAFACT = 0e0_hp

    ENDIF

    IF ( ExtState%FRAC_OF_PBL%DoUse ) THEN
       ALLOCATE( HCO_FRAC_OF_PBL( State_Grid%NX, State_Grid%NY, &
                                  State_Grid%NZ ), STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:HCO_FRAC_OF_PBL', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       HCO_FRAC_OF_PBL = 0.0_hp
    ENDIF

    ! Initialize max. PBL
    HCO_PBL_MAX = 0

    ! ----------------------------------------------------------------------
    ! The J-Values for NO2 and O3 are not defined in Met_State. We
    ! need to compute them separately.
    ! ----------------------------------------------------------------------
    IF ( ExtState%JNO2%DoUse ) THEN
       ALLOCATE( JNO2( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:JNO2', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       JNO2 = 0.0e0_hp
    ENDIF

    IF ( ExtState%JOH%DoUse ) THEN
       ALLOCATE( JOH( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:JOH', 0, RC )
       JOH = 0.0e0_hp
    ENDIF

    ! ----------------------------------------------------------------------
    ! Arrays to be copied physically because HEMCO units are not the
    ! same as in GEOS-Chem
    ! ----------------------------------------------------------------------

    ! TROPP: GEOS-Chem TROPP is in hPa, while HEMCO uses Pa.
    IF ( ExtState%TROPP%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%TROPP%Arr, State_Grid%NX, State_Grid%NY, &
                           HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCO_ArrAssert( TROPP )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    ! SPHU: GEOS-Chem SPHU is in g/kg, while HEMCO uses kg/kg.
    ! NOTE: HEMCO only uses SPHU surface values.
    IF ( ExtState%SPHU%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%SPHU%Arr, State_Grid%NX, State_Grid%NY, &
                           1, HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCO_ArrAssert( SPHU )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    ! FLASH_DENS
    IF ( ExtState%FLASH_DENS%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%FLASH_DENS%Arr, State_Grid%NX, &
                           State_Grid%NY, HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCO_ArrAssert( FLASH_DENS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    ! CONV_DEPTH
    IF ( ExtState%CONV_DEPTH%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%CONV_DEPTH%Arr, State_Grid%NX, &
                           State_Grid%NY, HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCO_ArrAssert( CONV_DEPTH )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    ! SUNCOS: HEMCO now calculates SUNCOS values based on its own
    ! subroutine
    IF ( ExtState%SUNCOS%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%SUNCOS%Arr, State_Grid%NX, State_Grid%NY, &
                           HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCO_ArrAssert( SUNCOS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    ! Leave with success
    RC = GC_SUCCESS

  END SUBROUTINE ExtState_InitTargets
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtState_SetFields
!
! !DESCRIPTION: SUBROUTINE ExtState\_SetFields connects the ExtState fields
! of the HEMCO ExtState object to its target data. This can be a field in
! State\_Met, State\_Chm, or any other 2D/3D field defined within GEOS-Chem or
! even explicitly calculated in this module. All these fields are expected
! to be of the same type as the corresponding ExtState object, and a pointer
! link is established between the two fields on the first call.
!\\
!\\
! The ExtState object fields can also be linked to data fields read through
! the HEMCO configuration file. In this case, the data fields will be copied
! from the HEMCO data list into the ExtState object every time this routine
! is called. The field name of the HEMCO field must match the field name
! passed to ExtState\_Set.
!\\
!\\
! Fields from the HEMCO data list are given priority over the target fields from
! State\_Met, State\_Chm, etc. For example, if the HEMCO data list contains
! a field named 'U10M', this field will be used in ExtState%U10M in lieu of
! State\_Met%U10M.
!\\
!\\
! Note that for now, this explicitly assumes that the HEMCO emissions grid is
! the same as the GEOS-Chem simulation grid. To support other grids, the met
! data has to be regridded explicitly at every time step!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtState_SetFields( State_Chm, State_Met, HcoState, ExtState, RC )
!
! !USES:
!
    USE Hcox_State_Mod, ONLY : ExtDat_Set
    USE ErrCode_Mod
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE Drydep_Mod,     ONLY : DryCoeff
#ifdef ESMF_
    USE HCOI_Esmf_Mod,  ONLY : HCO_SetExtState_ESMF
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState),   INTENT(INOUT)  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm  ! Chemistry state
    TYPE(HCO_STATE),  POINTER        :: HcoState   ! HEMCO state
    TYPE(EXT_STATE),  POINTER        :: ExtState   ! HEMCO ext. state
    INTEGER,          INTENT(INOUT)  :: RC         ! Return code
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller    - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Pointers
    REAL(hp), POINTER  :: Trgt3D(:,:,:) => NULL()

    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Scalars
    INTEGER            :: HMRC

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc, Instr
    CHARACTER(LEN=512) :: ErrMsg

    !=======================================================================
    ! ExtState_SetFields begins here
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at ExtState_SetFields (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! Pointers to local module arrays
    !-----------------------------------------------------------------------

    ! SZAFACT
    CALL ExtDat_Set( HcoState, ExtState%SZAFACT, 'SZAFACT_FOR_EMIS', &
                     HMRC,     FIRST,            HCO_SZAFACT )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( SZAFACT_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! JNO2
    CALL ExtDat_Set( HcoState, ExtState%JNO2, 'JNO2_FOR_EMIS', &
                     HMRC,     FIRST,         JNO2 )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( JNO2_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! JOH
    CALL ExtDat_Set( HcoState, ExtState%JOH, 'JOH_FOR_EMIS', &
                     HMRC,     FIRST,        JOH )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( JOH_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Frac of PBL
    CALL ExtDat_Set( HcoState, ExtState%FRAC_OF_PBL, 'FRAC_OF_PBL_FOR_EMIS', &
                     HMRC,     FIRST,                HCO_FRAC_OF_PBL )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRAC_OF_PBL_FOR_EMIS)"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! ----------------------------------------------------------------------
    ! 2D fields
    ! ----------------------------------------------------------------------

    ! U10M
    CALL ExtDat_Set( HcoState, ExtState%U10M, 'U10M_FOR_EMIS', &
                     HMRC,     FIRST,         State_Met%U10M )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( U10M_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! V10M
    CALL ExtDat_Set( HcoState, ExtState%V10M, 'V10M_FOR_EMIS', &
                     HMRC,     FIRST,         State_Met%V10M )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( V10M_FOR_EMIS)"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    ! ALBD
    CALL ExtDat_Set( HcoState, ExtState%ALBD, 'ALBD_FOR_EMIS', &
                     HMRC,     FIRST,         State_Met%ALBD )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( ALBD_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! WLI
    CALL ExtDat_Set( HcoState, ExtState%WLI, 'WLI_FOR_EMIS', &
                     HMRC,     FIRST,        State_Met%LWI )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( WLI_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    CALL ExtDat_Set( HcoState, ExtState%T2M, 'T2M_FOR_EMIS', &
                     HMRC,     FIRST,        State_Met%TS )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( T2M_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! TSKIN
    CALL ExtDat_Set( HcoState, ExtState%TSKIN, 'TSKIN_FOR_EMIS', &
                     HMRC,     FIRST,          State_Met%TSKIN )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( TSKIN_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! GWETROOT
    CALL ExtDat_Set( HcoState, ExtState%GWETROOT, 'GWETROOT_FOR_EMIS', &
                     HMRC,     FIRST,             State_Met%GWETROOT )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( GWETROOT_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! GWETTOP
    CALL ExtDat_Set( HcoState, ExtState%GWETTOP, 'GWETTOP_FOR_EMIS', &
                     HMRC,     FIRST,            State_Met%GWETTOP )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( GWETTOP_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    CALL ExtDat_Set( HcoState, ExtState%USTAR, 'USTAR_FOR_EMIS', &
                     HMRC,     FIRST,          State_Met%USTAR )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( USTAR_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Z0
    CALL ExtDat_Set( HcoState, ExtState%Z0, 'Z0_FOR_EMIS', &
                     HMRC,     FIRST,       State_Met%Z0 )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( Z0_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    CALL ExtDat_Set( HcoState, ExtState%PARDR, 'PARDR_FOR_EMIS', &
                     HMRC,     FIRST,          State_Met%PARDR )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( PARDR_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! PARDF
    CALL ExtDat_Set( HcoState, ExtState%PARDF, 'PARDF_FOR_EMIS', &
                     HMRC, FIRST,              State_Met%PARDF )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( PARDF_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! PSC2_WET
    CALL ExtDat_Set( HcoState, ExtState%PSC2_WET, 'PSC2_WET_FOR_EMIS', &
                     HMRC,     FIRST,             State_Met%PSC2_WET )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( PSC2_WET_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! RADSWG
    CALL ExtDat_Set( HcoState, ExtState%RADSWG, 'RADSWG_FOR_EMIS', &
                     HMRC,     FIRST,           State_Met%SWGDN )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( RADSWG_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! FRCLND
    CALL ExtDat_Set( HcoState, ExtState%FRCLND, 'FRCLND_FOR_EMIS', &
                     HMRC,     FIRST,           State_Met%FRCLND )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRCLND_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! CLDFRC
    CALL ExtDat_Set( HcoState, ExtState%CLDFRC, 'CLDFRC_FOR_EMIS', &
                     HMRC,     FIRST,           State_Met%CLDFRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( CLDFRC_FOR_EMIS)"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! SNOWHGT is is mm H2O, which is the same as kg H2O/m2.
    ! This is the unit of SNOMAS.
    CALL ExtDat_Set( HcoState, ExtState%SNOWHGT, 'SNOWHGT_FOR_EMIS', &
                     HMRC,     FIRST,            State_Met%SNOMAS )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( SNOWHGT_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! SNOWDP is in m
    CALL ExtDat_Set( HcoState, ExtState%SNODP, 'SNODP_FOR_EMIS', &
                     HMRC,    FIRST,           State_Met%SNODP )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( SNOWDP_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! FRLAND
    CALL ExtDat_Set( HcoState, ExtState%FRLAND, 'FRLAND_FOR_EMIS', &
                     HMRC,     FIRST,           State_Met%FRLAND )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRLAND_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! FROCEAN
    CALL ExtDat_Set( HcoState, ExtState%FROCEAN, 'FROCEAN_FOR_EMIS', &
                     HMRC,     FIRST,            State_Met%FROCEAN )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FROCEAN_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! FRLAKE
    CALL ExtDat_Set( HcoState, ExtState%FRLAKE, 'FRLAKE_FOR_EMIS', &
                     HMRC,     FIRST,           State_Met%FRLAKE )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRLAKE_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! FRLANDIC
    CALL ExtDat_Set( HcoState, ExtState%FRLANDIC, 'FRLANDIC_FOR_EMIS', &
                     HMRC,     FIRST,             State_Met%FRLANDIC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRLANDIC_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! LAI
    CALL ExtDat_Set( HcoState, ExtState%LAI, 'LAI_FOR_EMIS', &
                     HMRC,     FIRST,        State_Met%MODISLAI )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( LAI_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Convective fractions
    CALL ExtDat_Set( HcoState, ExtState%CNV_FRC,  'CNV_FRC_FOR_EMIS', &
                     HMRC,     FIRST,             State_Met%CNV_FRC,  &
                     NotFillOk=.TRUE. )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( CNV_FRC_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! 3D fields
    !-----------------------------------------------------------------------

    ! CNV_MFC
    CALL ExtDat_Set( HcoState,  ExtState%CNV_MFC, 'CNV_MFC_FOR_EMIS', &
                     HMRC,      FIRST,            State_Met%CMFMC,    &
                     OnLevEdge=.TRUE. )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( CNV_MFC_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    CALL ExtDat_Set( HcoState, ExtState%TK, 'TK_FOR_EMIS', &
                     HMRC,     FIRST,       State_Met%T )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( TK_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Air mass [kg/grid box]
    CALL ExtDat_Set( HcoState, ExtState%AIR, 'AIR_FOR_EMIS', &
                     HMRC,     FIRST,        State_Met%AD )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( AIR_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! AIRVOL_FOR_EMIS
    CALL ExtDat_Set( HcoState, ExtState%AIRVOL, 'AIRVOL_FOR_EMIS', &
                     HMRC,     FIRST,           State_Met%AIRVOL )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( AIRVOL_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Dry air density [kg/m3]
    CALL ExtDat_Set( HcoState, ExtState%AIRDEN, 'AIRDEN', &
                     HMRC,     FIRST,           State_Met%AIRDEN )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( AIRDEN )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Tropopause level
    CALL ExtDat_Set( HcoState, ExtState%TropLev, 'TropLev', &
                     HMRC,     FIRST,            State_Met%TropLev )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( TropLev )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! ----------------------------------------------------------------
    ! Species concentrations
    ! ----------------------------------------------------------------
    IF ( id_O3 > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_O3)
       CALL ExtDat_Set( HcoState, ExtState%O3, 'HEMCO_O3_FOR_EMIS', &
                        HMRC,     FIRST,       Trgt3D )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_O3_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

       Trgt3D => NULL()
    ENDIF

    IF ( id_NO2 > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_NO2)
       CALL ExtDat_Set( HcoState, ExtState%NO2, 'HEMCO_NO2_FOR_EMIS', &
                        HMRC,     FIRST,        Trgt3D )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_NO2_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

       Trgt3D => NULL()
    ENDIF

    IF ( id_NO > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_NO)
       CALL ExtDat_Set( HcoState, ExtState%NO, 'HEMCO_NO_FOR_EMIS', &
                        HMRC,     FIRST,       Trgt3D )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_NO_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

       Trgt3D => NULL()
    ENDIF

    IF ( id_HNO3 > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_HNO3)
       CALL ExtDat_Set( HcoState, ExtState%HNO3, 'HEMCO_HNO3_FOR_EMIS', &
                        HMRC,     FIRST,         Trgt3D )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_HNO3_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

       Trgt3D => NULL()
    ENDIF

    IF ( id_POPG > 0 ) THEN
       Trgt3D => State_Chm%Species(:,:,:,id_POPG)
       CALL ExtDat_Set( HcoState, ExtState%POPG, 'HEMCO_POPG_FOR_EMIS', &
                        HMRC,     FIRST,         Trgt3D )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_POPG_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

       Trgt3D => NULL()
    ENDIF

    ! ----------------------------------------------------------------------
    ! Deposition parameter
    ! ----------------------------------------------------------------------

    ! DRY_TOTN
    CALL ExtDat_Set( HcoState, ExtState%DRY_TOTN, 'DRY_TOTN_FOR_EMIS', &
                     HMRC,     FIRST,             State_Chm%DryDepNitrogen )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( DRY_TOTN )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! WET_TOTN
    CALL ExtDat_Set( HcoState, ExtState%WET_TOTN, 'WET_TOTN_FOR_EMIS', &
                     HMRC,     FIRST,             State_Chm%WetDepNitrogen )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( WET_TOTN )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! ----------------------------------------------------------------
    ! Other pointers to be set on first call
    ! ----------------------------------------------------------------
    IF ( FIRST ) THEN
       IF ( ExtState%WET_TOTN%DoUse .OR. ExtState%DRY_TOTN%DoUse ) THEN
          ExtState%DRYCOEFF => DRYCOEFF
       ENDIF

       ExtState%PBL_MAX => HCO_PBL_MAX
    ENDIF

    ! ----------------------------------------------------------------
    ! ESMF environment: add some additional variables to ExtState.
    ! These values must be defined here and not in the initialization
    ! because it seems like the IMPORT state is not yet properly
    ! defined during initialization.
    ! ckeller, 06/02/17: now call this on every time step. Routine
    ! HCO_SetExtState_ESMF copies the fields to ExtState.
    ! ----------------------------------------------------------------
#ifdef ESMF_
    ! IF ( FIRST ) THEN
    CALL HCO_SetExtState_ESMF ( HcoState, ExtState, RC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCO_SetExtState_ESMF"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF
    !ENDIF
#endif

    ! Not first call any more
    FIRST = .FALSE.

    ! Leave with success
    RC = GC_SUCCESS

  END SUBROUTINE ExtState_SetFields
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtState_UpdateFields
!
! !DESCRIPTION: SUBROUTINE ExtState\_UpdateFields updates the extension
! object data pointers. Updates are only required for the shadow arrays
! defined in this module, such as J-values, SZAFACT, etc.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtState_UpdateFields( Input_Opt,  State_Chm,             &
                                    State_Grid, State_Met, HcoState,   &
                                    ExtState,   RC )
!
! !USES:
!
    USE CMN_FJX_MOD,      ONLY : ZPJ
    USE ErrCode_Mod
    USE FAST_JX_MOD,      ONLY : RXN_NO2, RXN_O3_1, RXN_O3_2a
    USE HCO_GeoTools_Mod, ONLY : HCO_GetSUNCOS
    USE Input_Opt_Mod,    ONLY : OptInput
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Grid_Mod,   ONLY : GrdState
    USE State_Met_Mod,    ONLY : MetState
#ifdef ESMF_
    USE HCOI_ESMF_MOD,    ONLY : HCO_SetExtState_ESMF
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input options
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state
    TYPE(GrdState),   INTENT(IN   )  :: State_Grid ! Grid state
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Meteorology state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_STATE),  POINTER        :: HcoState   ! HEMCO state
    TYPE(EXT_STATE),  POINTER        :: ExtState   ! HEMCO ext. state
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  23 Oct 2012 - C. Keller   - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, J, L
    INTEGER            :: HMRC

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc, Instr
    CHARACTER(LEN=512) :: ErrMsg

    !=======================================================================
    ! ExtState_UpdateFields begins here
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at ExtState_UpdateFields (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !=======================================================================
    ! Update fields in the HEMCO Extension state
    !=======================================================================

    ! TROPP: convert from hPa to Pa
    IF ( ExtState%TROPP%DoUse ) THEN
       ExtState%TROPP%Arr%Val = State_Met%TROPP * 100.0_hp
    ENDIF

    ! SPHU: convert from g/kg to kg/kg. Only need surface value.
    IF ( ExtState%SPHU%DoUse ) THEN
       ExtState%SPHU%Arr%Val(:,:,1) = State_Met%SPHU(:,:,1) / 1000.0_hp
    ENDIF

    ! FLASH_DENS: flash density [#/km2/s]
    IF ( ExtState%FLASH_DENS%DoUse ) THEN
       ExtState%FLASH_DENS%Arr%Val = State_Met%FLASH_DENS
    ENDIF

    ! CONV_DEPTH: convective cloud depth [m]
    IF ( ExtState%CONV_DEPTH%DoUse ) THEN
       ExtState%CONV_DEPTH%Arr%Val = State_Met%CONV_DEPTH
    ENDIF

    ! SUNCOS
    IF ( ExtState%SUNCOS%DoUse ) THEN
       CALL HCO_GetSUNCOS( HcoState, ExtState%SUNCOS%Arr%Val, 0, HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "HCO_GetSuncos"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
    ENDIF

    ! If we need to use the SZAFACT scale factor (i.e. to put a diurnal
    ! variation on monthly mean OH concentrations), then call CALC_SUMCOSZA
    ! here.  CALC_SUMCOSZA computes the sum of cosine of the solar zenith
    ! angle over a 24 hour day, as well as the total length of daylight.
    ! This information is required by GET_SZAFACT. (bmy, 3/11/15)
    IF ( ExtState%SZAFACT%DoUse ) THEN
       CALL Calc_SumCosZa( State_Grid )
    ENDIF

!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J, L )
    ! Loop over all grid boxes
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Current SZA divided by total daily SZA (2D field only)
       ! (This is mostly needed for offline simulations where a diurnal
       ! scale factor has to be imposed on monthly mean OH concentrations.)
       IF ( ExtState%SZAFACT%DoUse .AND. L==1 ) THEN
          HCO_SZAFACT(I,J) = GET_SZAFACT(I,J,State_Met)
       ENDIF

       ! Fraction of PBL for each box [unitless]
       IF ( ExtState%FRAC_OF_PBL%DoUse ) THEN
          HCO_FRAC_OF_PBL(I,J,L) = State_Met%F_OF_PBL(I,J,L)
       ENDIF

       ! Maximum extent of the PBL [model level]
       HCO_PBL_MAX = State_Met%PBL_MAX_L

       ! J-values for NO2 and O3 (2D field only)
       ! This code was moved from hcox_paranox_mod.F90 to break
       ! dependencies to GC specific code (ckeller, 07/28/14).
       IF ( L==1 .AND.                                    &
            (ExtState%JNO2%DoUse .OR. ExtState%JOH%DoUse) ) THEN

          ! Check if sun is up
          IF ( State_Met%SUNCOS(I,J) == 0d0 ) THEN
             IF ( ExtState%JNO2%DoUse ) JNO2 = 0.0_hp
             IF ( ExtState%JOH%DoUse  ) JOH  = 0.0_hp
          ELSE
             IF ( ExtState%JNO2%DoUse ) THEN
                ! RXN_NO2: NO2 + hv --> NO  + O
                JNO2(I,J) = ZPJ(L,RXN_NO2,I,J)
             ENDIF
             IF ( ExtState%JOH%DoUse ) THEN
                IF ( Input_Opt%LUCX ) THEN
                   ! RXN_O3_1: O3  + hv --> O2  + O
                   JOH(I,J) = ZPJ(L,RXN_O3_1,I,J)
                ELSE
                   ! RXN_O3_2a: O3 + hv --> 2OH
                   JOH(I,J) = ZPJ(L,RXN_O3_2a,I,J)
                ENDIF
             ENDIF
          ENDIF

       ENDIF
    ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

  END SUBROUTINE ExtState_UpdateFields
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GridEdge_Set
!
! !DESCRIPTION: SUBROUTINE GridEdge\_Set sets the grid edge pressure values
! on the HEMCO grid.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GridEdge_Set( State_Grid, State_Met, HcoState, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCOX_STATE_MOD,   ONLY : ExtDat_Set
    USE HCO_GeoTools_Mod, ONLY : HCO_CalcVertGrid
    USE HCO_GeoTools_Mod, ONLY : HCO_SetPBLm
    USE State_Grid_Mod,   ONLY : GrdState
    USE State_Met_Mod,    ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState),   INTENT(IN   )  :: State_Grid ! Grid state
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(HCO_STATE),  POINTER        :: HcoState   ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC          ! Success or failure?
!
! !REMARKS:
!  GridEdge_Set defines the HEMCO vertical grid used in GEOS-Chem "classic"
!  simulations.  (GCHP uses its own interface to HEMCO.)
!
! !REVISION HISTORY:
!  08 Oct 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, J
    INTEGER            :: HMRC

    ! Strings
    CHARACTER(LEN=255) :: ThisLoc, Instr
    CHARACTER(LEN=512) :: ErrMsg

    ! Pointers
    REAL(hp), POINTER  :: PBLM    (:,:  )    ! PBL height           [m ]
    REAL(hp), POINTER  :: BXHEIGHT(:,:,:)    ! Grid box height      [m ]
    REAL(hp), POINTER  :: PEDGE   (:,:,:)    ! Pressure @ lvl edges [Pa]
    REAL(hp), POINTER  :: PSFC    (:,:  )    ! Surface pressure     [Pa]
    REAL(hp), POINTER  :: TK      (:,:,:)    ! Temperature          [K ]
    REAL(hp), POINTER  :: ZSFC    (:,:  )    ! Surface geopotential [m ]

    !=======================================================================
    ! GridEdge_Set begins here
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at GridEdge_Set (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! Allocate the PEDGE array, which holds level edge pressures [Pa]
    ! NOTE: Hco_CalcVertGrid expects pointer-based arguments, so we must
    ! make PEDGE be a pointer and allocate/deallocate it on each call.
    !-----------------------------------------------------------------------
    ALLOCATE( PEDGE( State_Grid%NX, State_Grid%NY, State_Grid%NZ+1 ), STAT=RC )
    CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:GridEdge_Set:PEDGE', 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Edge and surface pressures [Pa]
    PEDGE    =  State_Met%PEDGE * 100.0_hp  ! Convert hPa -> Pa
    PSFC     => PEDGE(:,:,1)

    ! Point to other fields of State_Met
    ZSFC     => State_Met%PHIS
    BXHEIGHT => State_Met%BXHEIGHT
    TK       => State_Met%T

    !-----------------------------------------------------------------------
    ! Calculate vertical grid properties
    !-----------------------------------------------------------------------
    CALL HCO_CalcVertGrid( HcoState, PSFC,  ZSFC, TK, BXHEIGHT, PEDGE, HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCO_CalcVertGrid"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Set PBL heights
    !-----------------------------------------------------------------------
    ALLOCATE( PBLM( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:GridEdge_Set:PBLM', 0, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J )
    DO J=1,State_Grid%NY
    DO I=1,State_Grid%NX
       PBLM(I,J) = State_Met%PBL_TOP_m(I,J)
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

    ! Use the met field PBL field to initialize HEMCO
    CALL HCO_SetPBLm( HcoState, FldName='PBL_HEIGHT', &
                      PBLM=PBLM, DefVal=1000.0_hp, RC=HMRC )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCO_SetPblM"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Cleanup and quit
    !-----------------------------------------------------------------------

    ! Deallocate and PEDGE
    IF ( ASSOCIATED( PEDGE ) ) THEN
       DEALLOCATE( PEDGE, STAT=RC )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:GridEdge_Set:PEDGE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! Deallocate PBLM
    IF ( ASSOCIATED( PBLM ) ) THEN
       DEALLOCATE( PBLM  )
       CALL GC_CheckVar( 'hcoi_gc_main_mod.F90:GridEdge_Set:PBLM', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! Free pointers
    ZSFC     => NULL()
    BXHEIGHT => NULL()
    TK       => NULL()
    PSFC     => NULL()
    PEDGE    => NULL()
    PBLM     => NULL()

  END SUBROUTINE GridEdge_Set
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetHcoSpecies
!
! !DESCRIPTION: Subroutine SetHcoSpecies defines the HEMCO species. These
! are typically just the GEOS-Chem species. Some additional species may be
! manually added, e.g. SESQ (which is not a species) or individual CO2 species
! per emission source (for CO2 specialty sim).
!\\
!\\
! This routine has two phases: phase 1 simply returns the number of species
! to be used by HEMCO. This is useful as this number needs to be passed to
! the HEMCO initialization call.
! Phase 2 sets the HEMCO species information in the HEMCO state object. This
! needs to be done after initialization of the HEMCO state object.
! !INTERFACE:
!
  SUBROUTINE SetHcoSpecies( Input_Opt, State_Chm, HcoState, nSpec, Phase, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_LogFile_Mod, ONLY : HCO_SPEC2LOG
    USE Input_Opt_Mod,   ONLY : OptInput
    USE Species_Mod,     ONLY : Species
    USE HCO_Types_Mod,   ONLY : ConfigObj
    USE State_Chm_Mod,   ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   )   :: Phase      ! 1=Init, 2=Run
    TYPE(ChmState),   INTENT(IN   )   :: State_Chm  ! Chemistry State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)   :: Input_Opt  ! Input Options object
    TYPE(Hco_State),  POINTER         :: HcoState   ! HEMCO state
    INTEGER,          INTENT(INOUT)   :: nSpec      ! # of species for HEMCO
    INTEGER,          INTENT(INOUT)   :: RC         ! Success or failure?
!
! !REMARKS:
!  (1) We now get physical parameters for species from the species database,
!       which is part of the State_Chm object.
!  (2) In the future, it will be easier to specify non-advected species
!       like SESQ and the CO2 regional species from the species database.
!       The species database flags if a species is advected or not.
!
! !REVISION HISTORY:
!  06 Mar 2015 - C. Keller   - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: nSpc, HMRC
    INTEGER                :: N,    L,    M
    REAL(dp)               :: K0,   CR,   pKa

    ! Strings
    CHARACTER(LEN= 31)     :: ThisName
    CHARACTER(LEN=255)     :: ThisLoc, Instr
    CHARACTER(LEN=512)     :: ErrMsg,  Msg

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! SetHcoSpecies begins here
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at SetHcoSpecies (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------
    ! For most simulations (e.g. full-chem simulation, most of the
    ! specialty sims), just use the GEOS-Chem species definitions.
    !-----------------------------------------------------------------
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM   .or. &
         Input_Opt%ITS_AN_AEROSOL_SIM   .or. &
         Input_OPt%ITS_A_CH4_SIM        .or. &
         Input_Opt%ITS_A_MERCURY_SIM    .or. &
         Input_Opt%ITS_A_POPS_SIM       .or. &
         Input_Opt%ITS_A_RnPbBe_SIM     .or. &
         Input_Opt%ITS_A_TAGO3_SIM      .or. &
         Input_Opt%ITS_A_TAGCO_SIM      .or. &
         Input_Opt%ITS_A_CO2_SIM              ) THEN

       ! Get number of model species
       nSpc = State_Chm%nAdvect

       !%%%%% FOR SOA SIMULATIONS %%%%%
       ! Check for SESQ: SESQ is not transported due to its short lifetime,
       ! but emissions are still calculated (in MEGAN). SESQ is only used
       ! in the SOA simulation, i.e. if LIMO is defined. Thus, add one more
       ! species here if LIMO is a model species and calculate SESQ emissions
       ! along with LIMO!
       IF ( id_LIMO > 0 ) THEN
          nSpc = nSpc + 1
       ENDIF

       !%%%%% FOR THE TAGGED CO SIMULATION %%%%%
       ! Add 5 extra species (ISOP, ACET, MTPA, LIMO, MTPO) for tagged CO
       IF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN
          nSpc = nSpc + 5
       ENDIF

       ! Assign species variables
       IF ( PHASE == 2 ) THEN

          ! Verbose
          IF ( Input_Opt%amIRoot ) THEN
             Msg = 'Registering HEMCO species:'
             CALL HCO_MSG( HcoState%Config%Err, Msg, SEP1='-' )
          ENDIF

          ! Sanity check: number of input species should agree with nSpc
          IF ( nSpec /= nSpc ) THEN
             WRITE(ErrMsg,*) 'Input species /= expected species: ', nSpec, nSpc
             CALL HCO_ERROR( HcoState%Config%Err, ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          DO N = 1, State_Chm%nAdvect

             ! Get info for this species from the species database
             SpcInfo => State_Chm%SpcData(N)%Info

             ! Model ID and species name
             HcoState%Spc(N)%ModID      = SpcInfo%ModelID
             HcoState%Spc(N)%SpcName    = TRIM( SpcInfo%Name )

             ! Actual molecular weight of species [g/mol]
             HcoState%Spc(N)%MW_g       = SpcInfo%MW_g

             ! Emitted molecular weight of species [g/mol].
             ! Some hydrocarbon species (like ISOP) are emitted and
             ! transported as a number of equivalent carbon atoms.
             ! For these species, the emitted molecular weight will
             ! be 12.0 (the weight of 1 carbon atom).
             HcoState%Spc(N)%EmMW_g     = SpcInfo%EmMw_g

             ! Emitted molecules per molecules of species [1].
             ! For most species, this will be 1.0.  For hydrocarbon
             ! species (like ISOP) that are emitted and transported
             ! as equivalent carbon atoms, this will be be the number
             ! of moles carbon per mole species.
             HcoState%Spc(N)%MolecRatio = SpcInfo%MolecRatio

             ! Set Henry's law coefficients
             HcoState%Spc(N)%HenryK0    = SpcInfo%Henry_K0   ! [M/atm]
             HcoState%Spc(N)%HenryCR    = SpcInfo%Henry_CR   ! [K    ]
             HcoState%Spc(N)%HenryPKA   = SpcInfo%Henry_pKa  ! [1    ]

             ! Write to logfile
             IF ( Input_Opt%amIRoot ) CALL HCO_SPEC2LOG( HcoState, N )

             ! Free pointer memory
             SpcInfo => NULL()
          ENDDO

          !------------------------------------------------------------------
          ! %%%%% FOR SOA SIMULATIONS %%%%%
          !
          ! Add the non-advected species SESQ in the last species slot
          !------------------------------------------------------------------
          IF ( id_LIMO > 0 ) THEN
             N                           = nSpec
             HcoState%Spc(N)%ModID       = N
             HcoState%Spc(N)%SpcName     = 'SESQ'
             HcoState%Spc(N)%MW_g        = 150.0_hp
             HcoState%Spc(N)%EmMW_g      = 150.0_hp
             HcoState%Spc(N)%MolecRatio  = 1.0_hp
             HcoState%Spc(N)%HenryK0     = 0.0_hp
             HcoState%Spc(N)%HenryCR     = 0.0_hp
             HcoState%Spc(N)%HenryPKa    = 0.0_hp

             ! Write to logfile
             IF ( Input_Opt%amIRoot ) CALL HCO_SPEC2LOG(  HcoState, N )
          ENDIF

          !------------------------------------------------------------------
          ! %%%%% FOR THE TAGGED CO SIMULATION %%%%%
          !
          ! Add the non-advected species ISOP, ACET, MTPA, LIMO, MTPO
          ! in the last 5 species slots (bmy, ckeller, 6/1/16)
          !------------------------------------------------------------------
          IF ( Input_Opt%ITS_A_TAGCO_SIM ) THEN

             ! Add 5 additional species
             DO L = 1, 5

                ! ISOP, ACET, MONX follow the regular tagged CO species
                M = State_Chm%nAdvect + L

                ! Get the species name
                SELECT CASE( L )
                   CASE( 1 )
                      ThisName = 'ISOP'
                   CASE( 2 )
                      ThisName = 'ACET'
                   CASE( 3 )
                      ThisName = 'MTPA'
                   CASE( 4 )
                      ThisName = 'LIMO'
                   CASE( 5 )
                      ThisName = 'MTPO'
                END SELECT

                ! Add physical properties to the HEMCO state
                HcoState%Spc(M)%ModID      = M
                HcoState%Spc(M)%SpcName    = TRIM( ThisName )
                HcoState%Spc(M)%MW_g       = 12.0_hp
                HcoState%Spc(M)%EmMW_g     = 12.0_hp
                HcoState%Spc(M)%MolecRatio = 1.0_hp
                HcoState%Spc(M)%HenryK0    = 0.0_hp
                HcoState%Spc(M)%HenryCR    = 0.0_hp
                HcoState%Spc(M)%HenryPKa   = 0.0_hp

                ! Write to log file
                IF ( Input_Opt%amIRoot ) CALL HCO_SPEC2LOG( HcoState, M )
             ENDDO
          ENDIF

          ! Add line to log-file
          IF ( Input_Opt%amIRoot ) CALL HCO_MSG( HcoState%Config%Err, SEP1='-' )
       ENDIF ! Phase = 2

    !-----------------------------------------------------------------
    ! DEFAULT (RETURN W/ ERROR)
    !-----------------------------------------------------------------
    ELSE
       ErrMsg = 'Invalid simulation type - cannot define model species'
       CALL HCO_ERROR( HcoState%Config%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! For phase 1, pass species to output
    nSpec = nSpc

    ! Return w/ success
    RC = HCO_SUCCESS

    END SUBROUTINE SetHcoSpecies
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetHcoGrid
!
! !DESCRIPTION: Subroutine SetHcoGrid tells HEMCO about the grid that is being
!  used by the GEOS-Chem simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SetHcoGrid( State_Grid, State_Met, HcoState, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_ARR_MOD,        ONLY : HCO_ArrInit
    USE HCO_VERTGRID_MOD,   ONLY : HCO_VertGrid_Define
    USE PRESSURE_MOD,       ONLY : GET_AP, GET_BP
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT ARGUMENTS:
!
    TYPE(GrdState),   INTENT(IN   )  :: State_Grid ! Grid state
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
!
! !INPUT/OUTPUT ARGUMENTS:
!
    TYPE(Hco_State),  POINTER        :: HcoState   ! HEMCO state
    INTEGER,          INTENT(INOUT)  :: RC
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller   - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    INTEGER               :: L
    INTEGER               :: HMRC

    ! Arrays
    REAL(hp), ALLOCATABLE :: Ap(:),   Bp(:)

    ! Strings
    CHARACTER(LEN=255)    :: ThisLoc, Instr
    CHARACTER(LEN=512)    :: ErrMsg

    !=======================================================================
    ! SetHcoGrid begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at SetHcoGrid (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !=======================================================================
    ! NOTE: for now, just copy GEOS-Chem grid, i.e. HEMCO calculations
    ! are performed on the GEOS-Chem simulation grid.
    ! It is possible to define a different emissions grid below.
    ! In this case, all arrays have to be regridded when passing
    ! them between HEMCO and GEOS-Chem (this is also true for the
    ! met-fields used by the extensions)!
    !=======================================================================

    ! Grid dimensions
    HcoState%NX = State_Grid%NX
    HcoState%NY = State_Grid%NY
    HcoState%NZ = State_Grid%NZ

    ! Allocate Ap array
    ALLOCATE( Ap( State_Grid%NZ+1 ), STAT=RC )
    CALL GC_CheckVar( 'hcoi_gc_main_mod:SetHcoGrid:Ap', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Allocate Bp array
    ALLOCATE( Bp( State_Grid%NZ+1 ), STAT=RC )
    CALL GC_CheckVar( 'hcoi_gc_main_mod:SetHcoGrid:Bp', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Get Ap and Bp values from GEOS-Chem
    DO L = 1, State_Grid%NZ+1
       Ap(L) = GET_AP(L) * 100_hp ! hPa to Pa
       Bp(L) = GET_BP(L)          ! unitless
    ENDDO

    ! Define the vertical grid
    CALL HCO_VertGrid_Define( HcoState%Config,                               &
                              zGrid      = HcoState%Grid%zGrid,              &
                              nz         = State_Grid%NZ,                    &
                              Ap         = Ap,                               &
                              Bp         = Bp,                               &
                              RC         = HMRC                             )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "HCO_VertGrid_Define"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Set pointers to grid variables
    HcoState%Grid%XMID%Val       => State_Grid%XMid   (:,:)
    HcoState%Grid%YMID%Val       => State_Grid%YMid   (:,:)
    HcoState%Grid%XEDGE%Val      => State_Grid%XEdge  (:,:)
    HcoState%Grid%YEDGE%Val      => State_Grid%YEdge  (:,:)
    HcoState%Grid%YSIN%Val       => State_Grid%YSIN   (:,:)
    HcoState%Grid%AREA_M2%Val    => State_Grid%Area_M2(:,:)
!    HcoState%Grid%ZSFC%Val       => State_Met%PHIS      ! Surface geopotential height
!    HcoState%Grid%BXHEIGHT_M%Val => State_Met%BXHEIGHT  ! Grid box heights

!    ! Allocate PEDGE. Will be updated every time step!
!    CALL HCO_ArrInit( HcoState%Grid%PEDGE, HcoState%NX, HcoState%NY, HcoState%NZ+1, RC )
!    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Return w/ success
    RC = HCO_SUCCESS

    END SUBROUTINE SetHcoGrid
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CheckSettings
!
! !DESCRIPTION: Subroutine CheckSettings performs some sanity checks of the
! switches provided in the HEMCO configuration file (in combination with the
! settings specified in input.geos).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CheckSettings( HcoConfig, Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Types_Mod,      ONLY : ConfigObj
    USE HCO_ExtList_Mod,    ONLY : GetExtNr,  SetExtNr
    USE HCO_ExtList_Mod,    ONLY : GetExtOpt, AddExtOpt
    USE HCO_ExtList_Mod,    ONLY : CoreNr
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(ConfigObj),  POINTER        :: HcoConfig  ! HEMCO config obj
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!  Gfortran will choke unless we use the .eqv. operator to compare LOGICAL
!  variables for equality (or .neqv. for inequality).

! !REVISION HISTORY:
!  18 Feb 2015 - C. Keller   - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: ExtNr
    INTEGER            :: HMRC
    LOGICAL            :: LTMP
    LOGICAL            :: FOUND

    ! Strings
    CHARACTER(LEN=31 ) :: OptName
    CHARACTER(LEN=255) :: ThisLoc, Instr
    CHARACTER(LEN=512) :: ErrMsg

    !=======================================================================
    ! CheckSettings begins here
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    HMRC     = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at CheckSettings (in module GeosCore/hcoi_gc_main_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! If emissions are turned off, do not read emissions data
    !-----------------------------------------------------------------------
    IF ( .NOT. Input_Opt%LEMIS ) THEN

       IF ( Input_Opt%amIRoot ) THEN
          Print*, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
          Print*, '% Emissions are set to false in input.geos so emissions %'
          Print*, '% data will not be read by HEMCO (hcoi_gc_main_mod.F90) %'
          Print*, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
       END IF

       OptName = 'EMISSIONS : false'
       CALL AddExtOpt( HcoConfig, TRIM(OptName), CoreNr, RC=HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "AddExtOpt( EMISSIONS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

       ! Reset all extension numbers to -999.
       ! This will make sure that none of the extensions will be initialized
       ! and none of the input data related to any of the extensions will be
       ! used.
       CALL SetExtNr( HcoConfig, -999, RC=HMRC                   )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "SetExtNr"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

    ENDIF

    !-----------------------------------------------------------------------
    ! If chemistry is turned off, do not read chemistry input data
    !-----------------------------------------------------------------------
    IF ( .NOT. Input_Opt%LCHEM ) THEN

       IF ( Input_Opt%amIRoot ) THEN
          Print*, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
          Print*, '% Chemistry is set to false in input.geos so chemistry  %'
          Print*, '% data will not be read by HEMCO (hcoi_gc_main_mod.F90) %'
          Print*, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
       ENDIF

       OptName = 'CHEMISTRY_INPUT : false'
       CALL AddExtOpt( HcoConfig, TRIM(OptName), CoreNr, RC=HMRC )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "AddExtOpt( CHEMISTRY_INPUT )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF

    ENDIF

    !-----------------------------------------------------------------------
    ! UV Albedo
    !
    ! UV albedoes are needed for photolysis.  Photolysis is only used in
    ! fullchem and aerosol-only simulations that have chemistry switched on.
    !-----------------------------------------------------------------------
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       CALL GetExtOpt( HcoConfig, -999, 'UVALBEDO',           &
                       OptValBool=LTMP, FOUND=FOUND,  RC=HMRC )

       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "GetExtOpt( UVALBEDO )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
       IF ( .not. FOUND ) THEN
          ErrMsg = 'UVALBEDO not found in HEMCO_Config.rc file!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( .not. LTMP ) THEN
          ErrMsg = 'UVALBEDO is set to false in HEMCO_Config.rc ' // &
                   'but should be set to true for this simulation.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    !-----------------------------------------------------------------------
    ! PSC STATE (for UCX)
    !-----------------------------------------------------------------------
#if !defined( ESMF_ )
    IF ( Input_Opt%LUCX ) THEN

       CALL GetExtOpt( HcoConfig, -999, 'STATE_PSC',          &
                       OptValBool=LTMP, FOUND=FOUND,  RC=HMRC )

       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "GetExtOpt( STATE_PSC )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
       IF ( .not. FOUND ) THEN
          ErrMsg = 'STATE_PSC not found in HEMCO_Config.rc file!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( .not. LTMP ) THEN
          ErrMsg = 'STATE_PSC is set to false in HEMCO_Config.rc ' // &
                   'but should be set to true for this simulation.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF
#endif

#ifdef ESMF_
    !-----------------------------------------------------------------------
    ! Also check that HEMCO_RESTART is not set in ESMF
    !-----------------------------------------------------------------------
    CALL GetExtOpt( HcoConfig,       -999,        'HEMCO_RESTART',           &
                    OptValBool=LTMP, FOUND=FOUND,  RC=HMRC                  )

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "GetExtOpt( HEMCO_RESTART in ESMF)"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF
    If ( FOUND .and. LTMP ) Then
       ErrMsg = 'Error encountered in "ESMF HEMCO_RESTART"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    End If
#endif

    !-----------------------------------------------------------------------
    ! TOMS/SBUV overhead O3 columns
    !
    ! Do not read in the TOMS/SBUV O3 columns unless running a mercury
    ! simulation. We will instead use the O3 columns from the GEOS-FP
    ! or MERRA-2 met fields.
    !-----------------------------------------------------------------------
    CALL GetExtOpt( HcoConfig, -999, 'TOMS_SBUV_O3',       &
                    OptValBool=LTMP, FOUND=FOUND,  RC=HMRC )
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "GetExtOpt( TOMS_SBUV_O3 )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    IF ( Input_Opt%ITS_A_MERCURY_SIM .and. Input_Opt%LKRedUV ) THEN

       IF ( .not. FOUND ) THEN
          ErrMsg = 'TOMS_SBUV_O3 not found in HEMCO_Config.rc file ' // &
                   'but is required for this simulation.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( .not. LTMP ) THEN
          ErrMsg = 'TOMS_SBUV_O3 is set to false in HEMCO_Config.rc ' // &
                   'but should be set to true for this simulation.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ELSE

       IF ( LTMP ) THEN
          ErrMsg = 'TOMS_SBUV_O3 is set to true in HEMCO_Config.rc ' // &
                   'but should be set to false for this simulation. '// &
                   'O3 columns are obtained from met fields instead.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    !-----------------------------------------------------------------------
    ! Ocean Hg input data (for Hg sims only)
    !
    ! If we have turned on the Ocean Mercury simulation in the
    ! input.geos file, then we will also toggle the OCEAN_Hg
    ! collection so that HEMCO reads the appropriate data.
    !-----------------------------------------------------------------------
    IF ( Input_Opt%ITS_A_MERCURY_SIM .and. Input_Opt%LDYNOCEAN ) THEN

       CALL GetExtOpt( HcoConfig, -999, 'OCEAN_Hg',           &
                       OptValBool=LTMP, FOUND=FOUND,  RC=HMRC )

       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "GetExtOpt( OCEAN_Hg )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
       IF ( .not. FOUND ) THEN
          ErrMsg = 'OCEAN_Hg not found in HEMCO_Config.rc file!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( .not. LTMP ) THEN
          ErrMsg = 'OCEAN_Hg is set to false in HEMCO_Config.rc ' // &
                   'but LDYNOCEAN is true in input.geos.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    !-----------------------------------------------------------------------
    ! RRTMG input data
    !
    ! If we have turned on the RRTMG simulation in the
    ! input.geos file, then we will also toggle the RRTMG
    ! collection so that HEMCO reads the appropriate data.
    !-----------------------------------------------------------------------
    IF ( Input_Opt%LRAD .and. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       CALL GetExtOpt( HcoConfig, -999, 'RRTMG',              &
                       OptValBool=LTMP, FOUND=FOUND,  RC=HMRC )

       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "GetExtOpt( RRTMG )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
       IF ( .not. FOUND ) THEN
          ErrMsg = 'RRTMG not found in HEMCO_Config.rc file!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       IF ( .not. LTMP ) THEN
          ErrMsg = 'RRTMG is set to false in HEMCO_Config.rc ' // &
                   'but should be set to false for this simulation.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE CheckSettings
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_szafact
!
! !DESCRIPTION:
!  Subroutine GET\_SZAFACT returns diurnal scale factors from dividing
!  the sza by the sum of the total sza per day. These factors are mainly
!  imposed to the monthly OH climatology.
!  However, the same scale factors are dimensionless and can hence be
!  applied to other compounds too (e.g. O3).
!\\
! !INTERFACE:
!
  FUNCTION Get_SzaFact( I, J, State_Met ) RESULT( FACT )
!
! !USES:
!
    USE State_Met_Mod,      ONLY : MetState
    USE Time_Mod,           ONLY : Get_TS_Chem
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: I, J
    TYPE(MetState), INTENT(IN) :: State_Met
!
! !RETURN VALUE:
!
    REAL(fp)                   :: FACT
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=======================================================================
    ! GET_SZAFACT begins here!
    !=======================================================================

    ! Test for sunlight...
    IF ( State_Met%SUNCOS(I,J) > 0e+0_fp  .AND. &
         SUMCOSZA(I,J)         > 0e+0_fp ) THEN

       ! Impose a diurnal variation on OH during the day
       FACT = ( State_Met%SUNCOS(I,J) / SUMCOSZA(I,J) ) &
              * ( 86400e+0_fp / GET_TS_CHEM() )

    ELSE

       ! At night, OH goes to zero
       FACT = 0e+0_fp

    ENDIF

  END FUNCTION Get_SzaFact
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_sumcosza
!
! !DESCRIPTION:
!  Subroutine CALC\_SUMCOSZA computes the sum of cosine of the solar zenith
!  angle over a 24 hour day, as well as the total length of daylight.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_SumCosZa( State_Grid )
!
! !USES:
!
    USE PhysConstants
    USE State_Grid_Mod, ONLY : GrdState
    USE TIME_MOD,       ONLY : GET_NHMSb,   GET_ELAPSED_SEC
    USE TIME_MOD,       ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER, SAVE :: SAVEDOY = -1
    INTEGER       :: I, J, L, N, NT, NDYSTEP
    REAL(fp)      :: A0, A1, A2, A3, B1, B2, B3
    REAL(fp)      :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
    REAL(fp)      :: SUNTMP(State_Grid%NX,State_Grid%NY)

    !=======================================================================
    ! CALC_SUMCOSZA begins here!
    !=======================================================================

    !  Solar declination angle (low precision formula, good enough for us):
    A0 = 0.006918
    A1 = 0.399912
    A2 = 0.006758
    A3 = 0.002697
    B1 = 0.070257
    B2 = 0.000907
    B3 = 0.000148
    R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.

    DEC = A0 - A1*cos(  R) + B1*sin(  R) &
             - A2*cos(2*R) + B2*sin(2*R) &
             - A3*cos(3*R) + B3*sin(3*R)

    LHR0 = int(float( GET_NHMSb() )/10000.)

    ! Only do the following at the start of a new day
    IF ( SAVEDOY /= GET_DAY_OF_YEAR() ) THEN

       ! Zero arrays
       SUMCOSZA(:,:) = 0e+0_fp

       ! NDYSTEP is # of chemistry time steps in this day
       NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 3600 / GET_TS_CHEM()

       ! NT is the elapsed time [s] since the beginning of the run
       NT = GET_ELAPSED_SEC()

       ! Loop forward through NDYSTEP "fake" timesteps for this day
       DO N = 1, NDYSTEP

          ! Zero SUNTMP array
          SUNTMP = 0e+0_fp

          ! Loop over surface grid boxes
!!$OMP PARALLEL DO
!!$OMP DEFAULT( SHARED )
!!$OMP PRIVATE( I, J, YMID_R, TIMLOC, AHR )
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Grid box latitude center [radians]
             YMID_R = State_Grid%YMid_R(I,J)

             TIMLOC = real(LHR0) + real(NT)/3600.0 + &
                      State_Grid%XMid(I,J) / 15.0

             DO WHILE (TIMLOC .lt. 0)
                TIMLOC = TIMLOC + 24.0
             ENDDO

             DO WHILE (TIMLOC .gt. 24.0)
                TIMLOC = TIMLOC - 24.0
             ENDDO

             AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

             !===========================================================
             ! The cosine of the solar zenith angle (SZA) is given by:
             !
             !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)
             !
             ! where LAT = the latitude angle,
             !       DEC = the solar declination angle,
             !       AHR = the hour angle, all in radians.
             !
             ! If SUNCOS < 0, then the sun is below the horizon, and
             ! therefore does not contribute to any solar heating.
             !===========================================================

             ! Compute Cos(SZA)
             SUNTMP(I,J) = sin(YMID_R) * sin(DEC) +          &
                           cos(YMID_R) * cos(DEC) * cos(AHR)

             ! SUMCOSZA is the sum of SUNTMP at location (I,J)
             ! Do not include negative values of SUNTMP
             SUMCOSZA(I,J) = SUMCOSZA(I,J) +             &
                             MAX(SUNTMP(I,J),0e+0_fp)

         ENDDO
         ENDDO
!!$OMP END PARALLEL DO

         ! Increment elapsed time [sec]
         NT = NT + GET_TS_CHEM()
      ENDDO

      ! Set saved day of year to current day of year
      SAVEDOY = GET_DAY_OF_YEAR()

   ENDIF

   ! Return to calling program
 END SUBROUTINE Calc_SumCosZa
!EOC
#if !defined(ESMF_) && !defined( MODEL_WRF )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_met_fields
!
! !DESCRIPTION: Subroutine GET\_MET\_FIELDS calls the various routines to get
! met fields from HEMCO.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE Get_Met_Fields( Input_Opt, State_Chm, State_Grid, &
                            State_Met, Phase, RC )
!
! ! USES:
!
   USE Calc_Met_Mod
   USE ErrCode_Mod
   USE FlexGrid_Read_Mod
   USE HCO_INTERFACE_MOD,      ONLY : HcoState
   USE HCO_EMISLIST_MOD,       ONLY : HCO_GetPtr
   USE Input_Opt_Mod,          ONLY : OptInput
   USE Pressure_Mod,           ONLY : Set_Floating_Pressures
   USE State_Chm_Mod,          ONLY : ChmState
   USE State_Grid_Mod,         ONLY : GrdState
   USE State_Met_Mod,          ONLY : MetState
   USE Time_Mod
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput),   INTENT(IN   )          :: Input_Opt  ! Input options
   TYPE(GrdState),   INTENT(IN   )          :: State_Grid ! Grid State
   INTEGER,          INTENT(IN   )          :: Phase      ! Run phase
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(MetState),   INTENT(INOUT)          :: State_Met  ! Meteorology State
   TYPE(ChmState),   INTENT(INOUT)          :: State_Chm  ! Chemistry State
   INTEGER,          INTENT(INOUT)          :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY:
!  07 Feb 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER              :: N_DYN              ! Dynamic timestep in seconds
   INTEGER              :: D(2)               ! Variable for date and time
   LOGICAL              :: FOUND              ! Found in restart file?
   LOGICAL              :: Update_MR          ! Update species mixing ratio?
   CHARACTER(LEN=255)   :: v_name             ! Variable name

   ! Pointers
   REAL*4,  POINTER     :: Ptr2D(:,:)
   REAL*4,  POINTER     :: Ptr3D(:,:,:)

   !=================================================================
   !    *****  R E A D   M E T   F I E L D S    *****
   !    *****  At the start of the GEOS-Chem simulation  *****
   !=================================================================

   ! Assume success
   RC        = GC_SUCCESS

   ! Initialize pointers
   Ptr2D       => NULL()
   Ptr3D       => NULL()

   !----------------------------------
   ! Read time-invariant data (Phase 0 only)
   !----------------------------------
   IF ( PHASE == 0 ) THEN
      CALL FlexGrid_Read_CN( Input_Opt, State_Grid, State_Met )
   ENDIF

   !----------------------------------
   ! Read 1-hr time-averaged data
   !----------------------------------
   IF ( PHASE == 0 ) THEN
      D = GET_FIRST_A1_TIME()
   ELSE
      D = GET_A1_TIME()
   ENDIF
   IF ( PHASE == 0 .or. ITS_TIME_FOR_A1() .and. &
        .not. ITS_TIME_FOR_EXIT() ) THEN
      CALL FlexGrid_Read_A1( D(1), D(2), Input_Opt, State_Grid, State_Met )
   ENDIF

   !----------------------------------
   ! Read 3-hr time averaged data
   !----------------------------------
   IF ( PHASE == 0 ) THEN
      D = GET_FIRST_A3_TIME()
   ELSE
      D = GET_A3_TIME()
   ENDIF
   IF ( PHASE == 0 .or. ITS_TIME_FOR_A3() .and. &
        .not. ITS_TIME_FOR_EXIT() ) THEN
      CALL FlexGrid_Read_A3( D(1), D(2), Input_Opt, State_Grid, State_Met )
   ENDIF

   !----------------------------------
   ! Read 3-hr instantanous data
   !----------------------------------
   IF ( PHASE == 0 ) THEN
      D = GET_FIRST_I3_TIME()
      CALL FlexGrid_Read_I3_1( D(1), D(2), Input_Opt, State_Grid, State_Met )

      ! On first call, attempt to get instantaneous met fields for prior
      ! timestep from the GEOS-Chem restart file. Otherwise, initialize
      ! to met fields for this timestep.

      !-------------
      ! TMPU
      !-------------

      ! Define variable name
      v_name = 'TMPU1'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Met%TMPU1 = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize TMPU1    from restart file'
         ENDIF
      ELSE
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'TMPU1    not found in restart, keep as value at t=0'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

      !-------------
      ! SPHU
      !-------------

      ! Define variable name
      v_name = 'SPHU1'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Met%SPHU1 = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize SPHU1    from restart file'
         ENDIF
      ELSE
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'SPHU1    not found in restart, keep as value at t=0'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

      !-------------
      ! PS1_WET
      !-------------

      ! Define variable name
      v_name = 'PS1WET'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr2D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Met%PS1_WET = Ptr2D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize PS1_WET  from restart file'
         ENDIF
      ELSE
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'PS1_WET  not found in restart, keep as value at t=0'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr2D => NULL()

      !-------------
      ! PS1_DRY
      !-------------

      ! Define variable name
      v_name = 'PS1DRY'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr2D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Met%PS1_DRY = Ptr2D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize PS1_DRY  from restart file'
         ENDIF
      ELSE
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'PS1_DRY  not found in restart, keep as value at t=0'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr2D => NULL()

      !-------------
      ! DELP_DRY
      !-------------

      ! Define variable name
      v_name = 'DELPDRY'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Met%DELP_DRY = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize DELP_DRY from restart file'
         ENDIF
      ELSE
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'DELP_DRY not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

      ! Set dry surface pressure (PS1_DRY) from State_Met%PS1_WET
      ! and compute avg dry pressure near polar caps
      CALL Set_Dry_Surface_Pressure( State_Grid, State_Met, 1 )
      CALL AvgPole( State_Grid, State_Met%PS1_DRY )

      ! Compute avg moist pressure near polar caps
      CALL AvgPole( State_Grid, State_Met%PS1_WET )

      ! Initialize surface pressures prior to interpolation
      ! to allow initialization of floating pressures
      State_Met%PSC2_WET = State_Met%PS1_WET
      State_Met%PSC2_DRY = State_Met%PS1_DRY
      CALL Set_Floating_Pressures( State_Grid, State_Met, RC )

      ! Call AIRQNT to compute initial air mass quantities
      ! Do not update initial tracer concentrations since not read
      ! from restart file yet (ewl, 10/28/15)
      CALL AirQnt( Input_Opt, State_Chm, State_Grid, State_Met, &
                   RC, update_mixing_ratio=.FALSE. )

   ENDIF

   ! Read in I3 fields at t+3hours for this timestep
   IF ( ITS_TIME_FOR_I3() .and. .not. ITS_TIME_FOR_EXIT() ) THEN

      D = GET_I3_TIME()
      CALL FlexGrid_Read_I3_2( D(1), D(2), Input_Opt, State_Grid, State_Met )
      
      ! Set dry surface pressure (PS2_DRY) from State_Met%PS2_WET
      ! and compute avg dry pressure near polar caps
      CALL Set_Dry_Surface_Pressure( State_Grid, State_Met, 2 )
      CALL AvgPole( State_Grid, State_Met%PS2_DRY )

      ! Compute avg moist pressure near polar caps
      CALL AvgPole( State_Grid, State_Met%PS2_WET )

   ENDIF

 END SUBROUTINE Get_Met_Fields
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_gc_restart
!
! !DESCRIPTION: Subroutine GET\_GC\_RESTART reads species concentrations
!  [mol/mol] from the GEOS-Chem restart file and uses them to initialize
!  species concentrations in [kg/kg dry]. If species data are missing from
!  the restart file, pre-configured background values are used. If using the
!  mercury simulation, additional restart data are read from file.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE Get_GC_Restart( Input_Opt, State_Chm, State_Grid, &
                            State_Met, RC )
!
! !USES:
!
   USE ErrCode_Mod
   USE Error_Mod
   USE HCO_INTERFACE_MOD,  ONLY : HcoState
   USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
   USE PHYSCONSTANTS,      ONLY : AIRMW
   USE Input_Opt_Mod,      ONLY : OptInput
   USE Species_Mod,        ONLY : Species
   USE State_Chm_Mod,      ONLY : ChmState
   USE State_Grid_Mod,     ONLY : GrdState
   USE State_Met_Mod,      ONLY : MetState
   USE TIME_MOD,           ONLY : EXPAND_DATE
   USE UnitConv_Mod,       ONLY : Convert_Spc_Units
#ifdef APM
   USE APM_INIT_MOD,       ONLY : APMIDS
#endif
   USE OCEAN_MERCURY_MOD,  ONLY : CHECK_OCEAN_MERCURY
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
   TYPE(GrdState), INTENT(IN)    :: State_Grid ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
   TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!
!  09 Feb 2016 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER              :: I, J, L, M, N      ! lon, lat, lev, indexes
   LOGICAL              :: FOUND              ! Found in restart file?
   CHARACTER(LEN=60)    :: Prefix             ! utility string
   CHARACTER(LEN=255)   :: LOC                ! routine location
   CHARACTER(LEN=255)   :: MSG                ! message
   CHARACTER(LEN=255)   :: v_name             ! variable name
   REAL(fp)             :: MW_g               ! species molecular weight
   REAL(fp)             :: SMALL_NUM          ! small number threshold
   CHARACTER(LEN=63)    :: OrigUnit

   ! Temporary arrays and pointers
   REAL*4,  TARGET           :: Temp2D(State_Grid%NX,State_Grid%NY)
   REAL*4,  TARGET           :: Temp3D(State_Grid%NX,State_Grid%NY, &
                                       State_Grid%NZ)
   REAL*4,  POINTER          :: Ptr2D(:,:  )
   REAL*4,  POINTER          :: Ptr3D(:,:,:)

   ! For Hg simulation
   INTEGER                   :: Num_Hg_Categories
   INTEGER                   :: Total_Hg_Id
   CHARACTER(LEN=60)         :: HgSpc
   CHARACTER(LEN=4), POINTER :: Hg_Cat_Name(:)

   ! Default background concentration
   REAL(fp)                  :: Background_VV

   ! Objects
   TYPE(Species),    POINTER :: SpcInfo

   !=================================================================
   ! READ_GC_RESTART begins here!
   !=================================================================

   ! Assume success
   RC        = GC_SUCCESS

   ! Initialize pointers
   Ptr2D       => NULL()
   Ptr3D       => NULL()
   SpcInfo     => NULL()
   Hg_Cat_Name => NULL()

   ! Name of this routine
   LOC = ' -> at Get_GC_Restart (in GeosCore/hcoi_gc_main_mod.F90)'

   ! Set minimum value threshold for [mol/mol]
   SMALL_NUM = 1.0e-30_fp

   !=================================================================
   ! If running Hg simulation, set Hg-specific local variables
   !=================================================================
   IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

      ! Set the # of tagHg categories from State_Chm
      Num_Hg_Categories   =  State_Chm%N_Hg_CATS

      ! Set variable storing names for each of the Hg categories
      Hg_Cat_Name => State_Chm%Hg_Cat_Name

      ! Set Hg species index corresponding to a given Hg category number;
      ! total is always the first category
      Total_Hg_Id   =  State_Chm%Hg0_Id_List(1)

   ENDIF

   !=================================================================
   ! Open GEOS-Chem restart file
   !=================================================================

   ! Write read message to log
   WRITE( 6, '(a)'   ) REPEAT( '=', 79 )
   WRITE( 6, '(a,/)' ) 'R E S T A R T   F I L E   I N P U T'

   !=================================================================
   ! Read species concentrations from NetCDF or use default
   ! background [mol/mol]; store in State_Chm%Species in [kg/kg dry]
   !=================================================================

   ! IMPORTANT NOTE: the unit conversion from mol/mol to kg/kg uses
   ! the molecular weight stored in the species database which is
   ! a meaningful value for advected species but is a bad value (-1)
   ! for all others. Non-advected species should NOT be used when
   ! State_Chm%Species units are in mass mixing ratio. Current
   ! units can be determined at any point by looking at
   ! State_Chm%Spc_Units. (ewl, 8/11/16)

   ! Print header for min/max concentration to log
   WRITE( 6, 110 )
110 FORMAT( 'Min and Max of each species in restart file [mol/mol]:' )

   ! Initialize species to all zeroes
   State_Chm%Species = 0.e+0_fp

   ! Loop over species
   DO N = 1, State_Chm%nSpecies

      ! Get info about this species from the species database
      SpcInfo => State_Chm%SpcData(N)%Info
      MW_g    =  SpcInfo%emMW_g

      ! Define variable name
      v_name = 'SPC_' // TRIM( SpcInfo%Name )

      ! Initialize temporary array for this species and point to it
      Temp3D = 0.0_fp
      Ptr3D => Temp3D

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), &
                       Ptr3D,    RC,       FOUND=FOUND )

      ! Check if species data is in file
      IF ( FOUND ) THEN
         SpcInfo%Is_InRestart = .TRUE.
      ELSE
         SpcInfo%Is_InRestart = .FALSE.
      ENDIF

      ! If data is in file, read in as [mol/mol] and convert to
      ! [kg/kg dry]. Otherwise, set to background value [mol/mol]
      ! either stored in species database (advected species all levels and
      ! non-advected species levels in the chemistry grid) or a small number
      ! (non-advected species levels above the chemistry grid) converted to
      ! [kg/kg dry]
      IF ( SpcInfo%Is_InRestart ) THEN

         ! Print the min & max of each species as it is read from
         ! the restart file in mol/mol
         IF ( Input_Opt%amIRoot ) THEN
            WRITE( 6, 120 ) N, TRIM( SpcInfo%Name ), &
                            MINVAL( Ptr3D ), MAXVAL( Ptr3D )
120         FORMAT( 'Species ', i3, ', ', a8, ': Min = ', es15.9, &
                    '  Max = ',es15.9)
         ENDIF

         ! Convert file value [mol/mol] to [kg/kg dry] for storage
!$OMP PARALLEL DO                                                       &
!$OMP DEFAULT( SHARED )                                                 &
!$OMP PRIVATE( I, J, L )
         DO L = 1, State_Grid%NZ
         DO J = 1, State_Grid%NY
         DO I = 1, State_Grid%NX
            State_Chm%Species(I,J,L,N) = Ptr3D(I,J,L) * MW_g / AIRMW
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ELSE

         ! Set species to the background value converted to [kg/kg dry]
!$OMP PARALLEL DO                                                       &
!$OMP DEFAULT( SHARED )                                                 &
!$OMP PRIVATE( I, J, L )
         ! Loop over all grid boxes
         DO L = 1, State_Grid%NZ
         DO J = 1, State_Grid%NY
         DO I = 1, State_Grid%NX

            ! For non-advected species at levels above chemistry grid,
            ! use a small number for background
            IF ( L > State_Grid%MaxChemLev .and. &
                     .NOT. SpcInfo%Is_Advected ) THEN

               State_Chm%Species(I,J,L,N) = SMALL_NUM * MW_g / AIRMW

            ! For all other cases, use the background value
            ! stored in the species database
            ELSE

               State_Chm%Species(I,J,L,N) = SpcInfo%BackgroundVV &
                                            * MW_g / AIRMW

               ! Print to log if debugging is on
               IF ( Input_Opt%amIRoot .AND. &
                    I == 1 .AND. J == 1 .AND. L == 1 ) THEN
                  WRITE( 6, 140 ) N, TRIM( SpcInfo%Name ), SpcInfo%BackgroundVV
140               FORMAT('Species ', i3, ', ', a9, &
                         ': Use background = ', es15.9)
               ENDIF


            ENDIF

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

#ifdef APM
         !================================================================
         ! APM MICROPHYSICS
         !================================================================
         WRITE(*,*)'APM run does not find '// TRIM( SpcInfo%Name ),N
         IF(SpcInfo%Name(1:9)=='APMSPBIN2')THEN
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED  ) &
!$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very
               ! low
               State_Chm%Species(I,J,L,N) = &
               State_Chm%Species(I,J,L,APMIDS%id_SO4)/20.D0
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:9)=='APMSPBIN3')THEN
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED  ) &
!$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very
               ! low
               State_Chm%Species(I,J,L,N) = &
               State_Chm%Species(I,J,L,APMIDS%id_SO4)/20.D0
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
           ENDIF
!GanLuotest
           IF(SpcInfo%Name(1:10)=='APMSEABIN0')THEN
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED  ) &
!$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very
               ! low
               State_Chm%Species(I,J,L,N) = &
               State_Chm%Species(I,J,L,APMIDS%id_SALA)/9.D0
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:10)=='APMSEABIN1')THEN
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED  ) &
!$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very
               ! low
               State_Chm%Species(I,J,L,N) = &
               State_Chm%Species(I,J,L,APMIDS%id_SALC)/10.D0
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:10)=='APMDSTBIN1')THEN
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED  ) &
!$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very low
               State_Chm%Species(I,J,L,N) = &
               SUM(State_Chm%Species(I,J,L,APMIDS%id_DST1:APMIDS%id_DST4))/6.D0
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:8)=='APMBCBIN')THEN
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED  ) &
!$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very low
               State_Chm%Species(I,J,L,N) = 1.D-30
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
           ENDIF
           IF(SpcInfo%Name(1:8)=='APMOCBIN')THEN
!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED  ) &
!$OMP PRIVATE( I, J, L )
            DO L = 1, State_Grid%NZ
            DO J = 1, State_Grid%NY
            DO I = 1, State_Grid%NX
               ! Apply minimum value threshold where input conc is very low
               State_Chm%Species(I,J,L,N) = 1.D-30
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
           ENDIF
#endif
      ENDIF

      ! Free pointer
      SpcInfo => NULL()

   ENDDO

   ! Set species units
   State_Chm%Spc_Units = 'kg/kg dry'

   ! If in debug mode, print out species min and max in [molec/cm3]
   IF ( Input_Opt%LPRT .and. Input_Opt%amIRoot ) THEN

      ! Convert units
      PRINT *, " "
      PRINT *, "Species min and max in molec/cm3"
      CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                              'molec/cm3', RC, OrigUnit=OrigUnit )

      ! Trap error
      IF ( RC /= GC_SUCCESS ) THEN
         Msg = 'Error returned from Convert_Spc_Units, call #1!'
         CALL GC_Error( Msg, RC, Loc )
         RETURN
      ENDIF

      ! Print values
      DO N = 1, State_Chm%nSpecies
         SpcInfo => State_Chm%SpcData(N)%Info
         WRITE(6,150) N, TRIM( SpcInfo%Name ),                 &
                         MINVAL( State_Chm%Species(:,:,:,N) ), &
                         MAXVAL( State_Chm%Species(:,:,:,N) )
150      FORMAT( 'Species ', i3, ', ', a9,                     &
                 ': Min = ', es15.9, ', Max = ', es15.9 )
         SpcInfo => NULL()
      ENDDO

      ! Convert units back
      CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                              OrigUnit,  RC )

      ! Trap error
      IF ( RC /= GC_SUCCESS ) THEN
         Msg = 'Error returned from Convert_Spc_Units, call #2!'
         CALL GC_Error( Msg, RC, Loc )
         RETURN
      ENDIF

   ENDIF

   !=================================================================
   ! Get variables for FlexChem
   !=================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM .and. Input_Opt%LCHEM ) THEN

      ! Define variable name
      v_name = 'KPP_HVALUE'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%KPPHvalue = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize KPP H-value from restart file'
            WRITE(6,160) MINVAL( State_Chm%KPPHvalue(:,:,:) ), &
                         MAXVAL( State_Chm%KPPHvalue(:,:,:) )
160         FORMAT( 'KPP_HVALUE: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%KPPHvalue = 0e+0_fp
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'KPP_HVALUE     not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

   ENDIF

   !=================================================================
   ! Get variables for Soil NOx emissions
   !=================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

      ! Define variable name
      v_name = 'WETDEP_N'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr2D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%WetDepNitrogen = Ptr2D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize wet deposited nitrogen from restart file'
            WRITE(6,170) MINVAL( State_Chm%WetDepNitrogen(:,:) ), &
                         MAXVAL( State_Chm%WetDepNitrogen(:,:) )
170         FORMAT( 12x, '  WETDEP_N: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%WetDepNitrogen = 0e+0_fp
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'WETDEP_N       not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr2D => NULL()

      ! Define variable name
      v_name = 'DRYDEP_N'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr2D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%DryDepNitrogen = Ptr2D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize dry deposited nitrogen from restart file'
            WRITE(6,180) MINVAL( State_Chm%DryDepNitrogen(:,:) ), &
                         MAXVAL( State_Chm%DryDepNitrogen(:,:) )
180         FORMAT( 12x, '  DRYDEP_N: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%DryDepNitrogen = 0e+0_fp
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'DRYDEP_N       not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr2D => NULL()

   ENDIF

   !=================================================================
   ! Read variables for sulfate chemistry
   !=================================================================
   IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. &
        Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

      ! Define variable name
      v_name = 'H2O2_AFTERCHEM'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%H2O2AfterChem = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize H2O2 from restart file'
            WRITE(6,190) MINVAL( State_Chm%H2O2AfterChem(:,:,:) ), &
                         MAXVAL( State_Chm%H2O2AfterChem(:,:,:) )
190         FORMAT( 12x, 'H2O2_AChem: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%H2O2AfterChem = 0e+0_fp
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'H2O2_AFTERCHEM not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

      ! Define variable name
      v_name = 'SO2_AFTERCHEM'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%SO2AfterChem = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize dry deposited nitrogen from restart file'
            WRITE(6,200) MINVAL( State_Chm%SO2AfterChem(:,:,:) ), &
                         MAXVAL( State_Chm%SO2AfterChem(:,:,:) )
200         FORMAT( 12x, ' SO2_AChem: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         State_Chm%SO2AfterChem = 0e+0_fp
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'SO2_AFTERCHEM  not found in restart, set to zero'
         ENDIF
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

   ENDIF

   !=================================================================
   ! Read variables for UCX
   !=================================================================
   IF ( Input_Opt%LUCX ) THEN

      ! Define variable name
      v_name = 'STATE_PSC'

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )

      ! Check if variable is in file
      IF ( FOUND ) THEN
         State_Chm%STATE_PSC = Ptr3D
         IF ( Input_Opt%amIRoot ) THEN
            WRITE(6,*) 'Initialize PSC from restart for UCX'
            WRITE(6,210) MINVAL( State_Chm%STATE_PSC(:,:,:) ), &
                         MAXVAL( State_Chm%STATE_PSC(:,:,:) )
210         FORMAT( 12x, ' STATE_PSC: Min = ', es15.9, ', Max = ', es15.9 )
         ENDIF
      ELSE
         IF ( Input_OPt%amIRoot ) THEN
#ifdef ESMF_
            ! ExtData and HEMCO behave ambiguously - if the file was found
            ! but was full of zeros throughout the domain of interest, it
            ! will result in the same output from ExtData as if the field
            ! was missing from the file. As such, HEMCO cannot distinguish
            ! between a missing file and a field of zeros
            WRITE(6,*) 'PSC restart either all zeros in the '
            WRITE(6,*) 'root domain, or the restart file did '
            WRITE(6,*) 'not contain STATE_PSC. Root domain '
            WRITE(6,*) 'will be initialized PSC-free'
         ENDIF
#else
            WRITE(6,*) 'STATE_PSC      not found in restart, initialize PSC-free'
         ENDIF
#endif
      ENDIF

      ! Nullify pointer
      Ptr3D => NULL()

   ENDIF

   !=================================================================
   ! Read ocean mercury variables
   !=================================================================
   IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

      ! Print total mass to log
      WRITE( 6, 220 )
220   FORMAT(/, 'Total mass of each ocean and snow Hg species:')

      !--------------------------------------------------------------
      ! Total Hg in ocean
      !--------------------------------------------------------------
      DO M = 1, 3

         ! Define variable name
         SELECT CASE( M )
         CASE ( 1 )
            HgSpc    = 'Hg0'
         CASE ( 2 )
            HgSpc    = 'Hg2'
         CASE ( 3 )
            HgSpc    = 'HgP'
         END SELECT
         v_name = 'OCEAN_' // TRIM( HgSpc )

         ! Get variable from HEMCO and store in local array
         CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr2D, RC, FOUND=FOUND )

         ! Check if variable is in file
         IF ( FOUND ) THEN

            ! Check for negative concentrations (jaf, 7/6/11)
            DO I = 1, State_Grid%NX
            DO J = 1, State_Grid%NY
               IF ( Ptr2D(I,J) < 0.0d4 ) THEN
                  Ptr2D(I,J) = 0.0d4
               ENDIF
            ENDDO
            ENDDO

            ! Assign ocean mercury data and write total mass to log file
            SELECT CASE( M )
            CASE ( 1 )
               State_Chm%OceanHg0(:,:,Total_Hg_Id) = Ptr2D
               WRITE( 6, 240 ) TRIM( v_name ), &
                            SUM( State_Chm%OceanHg0(:,:,Total_Hg_Id) ), 'kg'
            CASE ( 2 )
               State_Chm%OceanHg2(:,:,Total_Hg_Id) = Ptr2D
               WRITE( 6, 240 ) TRIM( v_name ),  &
                            SUM( State_Chm%OceanHg2(:,:,Total_Hg_Id) ), 'kg'
            CASE ( 3 )
               State_Chm%OceanHgP(:,:,Total_Hg_Id) = Ptr2D
               WRITE( 6, 240 ) TRIM( v_name ),  &
                            SUM( State_Chm%OceanHgP(:,:,Total_Hg_Id) ), 'kg'
            END SELECT

         ELSE
            WRITE( 6, 230 ) TRIM( v_name )
         ENDIF

         ! Nullify pointer
         Ptr2D => NULL()

      ENDDO

      !--------------------------------------------------------------
      ! Additional tagged ocean Hg species
      !--------------------------------------------------------------
      IF ( Input_Opt%LSPLIT ) THEN
         DO M = 1, 3
            DO N = 2, Num_Hg_Categories

               ! Define variable name. Include appended region.
               SELECT CASE( M )
               CASE ( 1 )
                  HgSpc = 'Hg0'
               CASE ( 2 )
                  HgSpc = 'Hg2'
               CASE ( 3 )
                  HgSpc = 'HgP'
               END SELECT
               v_name = 'OCEAN_' // TRIM( HgSpc ) //  &
                        '_'      // TRIM( Hg_Cat_Name(N) )

               ! Get variable from HEMCO and store in local array
               CALL HCO_GetPtr( HcoState, TRIM(v_name), &
                                Ptr2D, RC, FOUND=FOUND )

               ! Check if variable is in file
               IF ( FOUND ) THEN

                  ! Assign ocean mercury data and write total mass to log
                  SELECT CASE( M )
                  CASE ( 1 )
                     State_Chm%OceanHg0(:,:,N) = Ptr2D
                     WRITE( 6, 240 ) TRIM( v_name ),  &
                                     SUM( State_Chm%OceanHg0(:,:,N) ), 'kg'
                  CASE ( 2 )
                     State_Chm%OceanHg2(:,:,N) = Ptr2D
                     WRITE( 6, 240 ) TRIM( v_name ),  &
                                     SUM( State_Chm%OceanHg2(:,:,N) ), 'kg'
                  CASE ( 3 )
                     State_Chm%OceanHgP(:,:,N) = Ptr2D
                     WRITE( 6, 240 ) TRIM( v_name ),  &
                                     SUM( State_Chm%OceanHgP(:,:,N) ), 'kg'
                  END SELECT

               ELSE
                  WRITE( 6, 230 ) TRIM( v_name )
               ENDIF

               ! Nullify pointer
               Ptr2D => NULL()

            ENDDO
         ENDDO

         ! Make sure tagged & total species sum up
         IF ( Input_Opt%USE_CHECKS ) THEN
            CALL CHECK_OCEAN_MERCURY( State_Chm, State_Grid, &
                                      'end of READ_GC_RESTART' )
         ENDIF
      ENDIF

      !--------------------------------------------------------------
      ! Hg snowpack on land and ocean
      !--------------------------------------------------------------
      DO M = 1, 4
         DO N = 1, Num_Hg_Categories

            ! Define variable name prefix
            SELECT CASE( M )
            CASE ( 1 )
               Prefix = 'SNOW_HG_OCEAN'        ! Reducible on ocean
            CASE ( 2 )
               Prefix = 'SNOW_HG_OCEAN_STORED' ! Non-reducible on ocean
            CASE ( 3 )
               Prefix = 'SNOW_HG_LAND'         ! Reducible on land
            CASE ( 4 )
               Prefix = 'SNOW_HG_LAND_STORED'  ! Non-reducible on land
            END SELECT

            IF ( N == 1 ) THEN
               v_name = TRIM( Prefix )
            ELSE
               ! Append category name if tagged
               v_name = TRIM( Prefix         ) // '_' // &
                        TRIM( Hg_Cat_Name(N) )
            ENDIF

            ! Get variable from HEMCO and store in local array
            CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr2D, RC, FOUND=FOUND )

            ! Check if variable is in file
            IF ( FOUND ) THEN

               ! Assign ocean mercury data and write total mass to file
               SELECT CASE( M )
               CASE ( 1 )
                  State_Chm%SnowHgOcean(:,:,N) = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name ),  &
                                  SUM( State_Chm%SnowHgOcean(:,:,N) ), 'kg'
               CASE ( 2 )
                  State_Chm%SnowHgOceanStored(:,:,N) = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name ),  &
                                  SUM( State_Chm%SnowHgOceanStored(:,:,N) ),'kg'
               CASE ( 3 )
                  State_Chm%SnowHgLand(:,:,N) = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name ),  &
                                  SUM( State_Chm%SnowHgLand(:,:,N) ), 'kg'
               CASE ( 4 )
                  State_Chm%SnowHgLandStored(:,:,N) = Ptr2D
                  WRITE( 6, 240 ) TRIM( v_name ),  &
                                  SUM( State_Chm%SnowHgLandStored(:,:,N) ), 'kg'
               END SELECT

            ELSE
               WRITE( 6, 230 ) TRIM( v_name )
            ENDIF

            ! Nullify pointer
            Ptr2D => NULL()

         ENDDO
      ENDDO

      ! Format strings
230   FORMAT( a24, ' not found in restart file, set to zero')
240   FORMAT( a24, ':   ', es15.9, 1x, a4)

      ! Print note that variables are initialized to zero if not
      ! found (currently only happens in tagged Hg simulation)
      IF ( Input_Opt%LSPLIT ) THEN
         WRITE( 6, 250 )
250      FORMAT( /, 'NOTE: all variables not found in restart ', &
                    'are initialized to zero')
      ENDIF

      ! Free pointers for Hg indexing
      Hg_Cat_Name => NULL()

   ENDIF

   !=================================================================
   ! Clean up
   !=================================================================

   ! Mark end of section in log
   IF ( Input_Opt%LPRT .AND. Input_Opt%amIRoot ) THEN
      CALL DEBUG_MSG('### DONE GET_GC_RESTART')
   ENDIF
   WRITE( 6, '(a)' ) REPEAT( '=', 79 )

 END SUBROUTINE Get_GC_Restart
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_boundary_conditions
!
! !DESCRIPTION: Subroutine GET\_BOUNDARY\_CONDITIONS calls the various routines
! to get boundary conditions from HEMCO for nested grid simulations.
!\\
!\\
! !INTERFACE:
!
 SUBROUTINE Get_Boundary_Conditions( Input_Opt, State_Chm, State_Grid, &
                                     State_Met, YYYYMMDD,  HHMMSS, RC )
!
! ! USES:
!
   USE ErrCode_Mod
   USE HCO_INTERFACE_MOD,      ONLY : HcoState
   USE HCO_EMISLIST_MOD,       ONLY : HCO_GetPtr
   USE Input_Opt_Mod,          ONLY : OptInput
   USE PHYSCONSTANTS,          ONLY : AIRMW
   USE Species_Mod,            ONLY : Species
   USE State_Chm_Mod,          ONLY : ChmState
   USE State_Grid_Mod,         ONLY : GrdState
   USE State_Met_Mod,          ONLY : MetState
   USE Time_Mod
!
! !INPUT PARAMETERS:
!
   TYPE(OptInput),   INTENT(IN   )          :: Input_Opt  ! Input options
   TYPE(GrdState),   INTENT(IN   )          :: State_Grid ! Grid State
   INTEGER,          INTENT(IN   )          :: YYYYMMDD   ! GMT date
   INTEGER,          INTENT(IN   )          :: HHMMSS     ! GMT time
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(MetState),   INTENT(INOUT)          :: State_Met  ! Meteorology State
   TYPE(ChmState),   INTENT(INOUT)          :: State_Chm  ! Chemistry State
   INTEGER,          INTENT(INOUT)          :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Apr 2019 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER              :: I, J, L, N, NA     ! lon, lat, lev, spc indexes
   INTEGER              :: t_index            ! Time index
   LOGICAL              :: FOUND              ! Found in restart file?
   CHARACTER(LEN=60)    :: Prefix             ! utility string
   CHARACTER(LEN=255)   :: LOC                ! routine location
   CHARACTER(LEN=255)   :: MSG                ! message
   CHARACTER(LEN=255)   :: v_name             ! variable name
   REAL(fp)             :: MW_g               ! species molecular weight
   CHARACTER(LEN=16)    :: STAMP

   ! Temporary arrays and pointers
   REAL*4,  TARGET      :: Temp3D(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
   REAL*4,  POINTER     :: Ptr3D(:,:,:)
   REAL(fp), POINTER    :: Spc(:,:,:,:)

   ! Objects
   TYPE(Species), POINTER :: SpcInfo

   !=================================================================
   ! GET_BOUNDARY_CONDITIONS begins here!
   !=================================================================

   ! Assume success
   RC        = GC_SUCCESS

   ! Initialize pointers
   Ptr3D     => NULL()
   SpcInfo   => NULL()

   ! Point to species array [kg/kg]
   Spc       => State_Chm%Species

   ! Name of this routine
   LOC = ' -> at Get_Boundary_Conditions (in GeosCore/hcoi_gc_main_mod.F90)'

   ! Find the proper time-slice to read from disk
   t_index = ( HHMMSS / 030000 ) + 1

   ! Stop w/ error if the time index is invalid
   IF ( t_index < 1 .or. t_index > 8 ) THEN
      WRITE( MSG, 100 ) t_index
100   FORMAT( 'Time_index value ', i5, ' must be in the range 1 to 8!' )
      CALL GC_Error( MSG, RC, LOC)
      RETURN
   ENDIF

   !=================================================================
   ! Read species concentrations from NetCDF [mol/mol] and
   ! store in State_Chm%BoundaryCond in [kg/kg dry]
   !=================================================================

   ! IMPORTANT NOTE: the unit conversion from mol/mol to kg/kg uses
   ! the molecular weight stored in the species database which is
   ! a meaningful value for advected species but is a bad value (-1)
   ! for all others. Non-advected species should NOT be used when
   ! State_Chm%Species units are in mass mixing ratio. Current
   ! units can be determined at any point by looking at
   ! State_Chm%Spc_Units. (ewl, 8/11/16)

   ! Initialize BCs to all zeroes
   State_Chm%BoundaryCond = 0.e+0_fp

   ! Loop over advected species
   DO NA = 1, State_Chm%nAdvect

      ! Get the species ID from the advected species ID
      N = State_Chm%Map_Advect(NA)

      ! Get info about this species from the species database
      SpcInfo => State_Chm%SpcData(N)%Info
      MW_g    =  SpcInfo%emMW_g

      ! Define variable name
      v_name = 'BC_' // TRIM( SpcInfo%Name )

      ! Initialize temporary array for this species and point to it
      Temp3D = 0.0_fp
      Ptr3D => Temp3D

      ! Get variable from HEMCO and store in local array
      CALL HCO_GetPtr( HcoState, TRIM(v_name), Ptr3D, RC, &
                       TIDX=t_index, FOUND=FOUND )

      ! Check if BCs are found
      IF ( FOUND ) THEN

         ! Copy data from file to State_Chm%BoundaryCond
         ! and convert from [mol/mol] to [kg/kg dry]
         State_Chm%BoundaryCond(:,:,:,N) = Ptr3D(:,:,:) * MW_g / AIRMW

         ! Debug
!         Print*, 'BCs found for ', TRIM( SpcInfo%Name ), &
!                 MINVAL(State_Chm%BoundaryCond(:,:,:,N)), &
!                 MAXVAL(State_Chm%BoundaryCond(:,:,:,N))

         ! Loop over grid boxes and apply BCs to the specified buffer zone
!$OMP PARALLEL DO                                                       &
!$OMP DEFAULT( SHARED )                                                 &
!$OMP PRIVATE( I, J, L )
         DO L = 1, State_Grid%NZ

            ! First loop over all latitudes of the nested domain
            DO J = 1, State_Grid%NY

               ! West BC
               DO I = 1, State_Grid%WestBuffer
                  State_Chm%Species(I,J,L,N) = State_Chm%BoundaryCond(I,J,L,N)
               ENDDO

               ! East BC
               DO I = (State_Grid%NX-State_Grid%EastBuffer)+1, State_Grid%NX
                  State_Chm%Species(I,J,L,N) = State_Chm%BoundaryCond(I,J,L,N)
               ENDDO

            ENDDO

            ! Then loop over the longitudes of the nested domain
            DO I = 1+State_Grid%WestBuffer,(State_Grid%NX-State_Grid%EastBuffer)

               ! South BC
               DO J = 1, State_Grid%SouthBuffer
                  Spc(I,J,L,N) = State_Chm%BoundaryCond(I,J,L,N)
               ENDDO

               ! North BC
               DO J = (State_Grid%NY-State_Grid%NorthBuffer)+1, State_Grid%NY
                  Spc(I,J,L,N) = State_Chm%BoundaryCond(I,J,L,N)
               ENDDO
            ENDDO

         ENDDO
!OMP END PARALLEL DO

      ELSE

         MSG = 'No boundary condition found for '// TRIM( SpcInfo%Name )
         CALL GC_Error( MSG, RC, LOC)
         RETURN

      ENDIF

      ! Free pointer
      SpcInfo => NULL()

   ENDDO

   ! Echo output
   STAMP = TIMESTAMP_STRING()
   WRITE( 6, 110 ) STAMP
110 FORMAT( 'GET_BOUNDARY_CONDITIONS: Found All BCs at ', a )


 END SUBROUTINE Get_Boundary_Conditions
!EOC
#endif
END MODULE Hcoi_GC_Main_Mod
