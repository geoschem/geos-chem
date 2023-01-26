!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!                    Harmonized Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_interface_gc_mod.F90
!
! !DESCRIPTION: Module hco\_interface\_gc\_mod.F90 contains routines and
! variables to interface GEOS-Chem and HEMCO. It contains the HEMCO
! state object (HcoState) as well as init-run-finalize driver routines
! to run HEMCO within GEOS-Chem.
!\\
!\\
! The HEMCO driver is now present in this file as HEMCO is restructured to provide
! a unified point-of-entry for coupling with other models.
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
MODULE HCO_Interface_GC_Mod
!
! !USES:
!
  USE Precision_Mod

  USE HCO_Error_Mod
  USE HCO_Interface_Common

  ! Import the HEMCO states and their types from the state container
  USE HCOX_State_Mod, ONLY : Ext_State
  USE HCO_State_Mod,  ONLY : HCO_State
  USE HCO_State_GC_Mod, ONLY : HcoState, ExtState

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOI_GC_Init
  PUBLIC  :: HCOI_GC_Run
  PUBLIC  :: HCOI_GC_Final
  PUBLIC  :: HCOI_GC_WriteDiagn

  PUBLIC  :: Compute_Sflx_For_Vdiff
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: ExtState_InitTargets
  PRIVATE :: ExtState_SetFields
  PRIVATE :: ExtState_UpdateFields
  PRIVATE :: Get_SzaFact
  PRIVATE :: GridEdge_Set
  PRIVATE :: CheckSettings
  PRIVATE :: SetHcoGrid
  PRIVATE :: SetHcoSpecies

#if defined( MODEL_CLASSIC )
  PRIVATE :: Get_Met_Fields
#endif
!
! !REMARKS:
!  Formerly HCOI\_GC\_Main\_Mod.
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

#if defined( MODEL_CLASSIC )

  !------------------------------------------------
  ! %%% Internal HEMCO intermediate resolution %%%
  ! %%%        meteorological fields           %%%
  !------------------------------------------------

  ! These will ONLY be allocated and used if HEMCO-intermediate grid
  ! is enabled (resolutions differ)
  !
  ! Note: These could be migrated to inside a derived-type object
  ! by storing directly in ExtState through array assertions. However,
  ! this would involve a lot of code path branching as we do not want
  ! to assert arrays for other models when not running IMGrid and incurring
  ! waste of memory. For the sake of maintaining old code implementations
  ! "as-is" and making the intermediate an as-seamless-possible patch over
  ! the original implementation, we list all these module variables here.
  ! I am aware that there is a cleaner implementation. (hplin, 6/9/20)

  ! 2-D fields
  REAL(hp), POINTER :: H_SZAFACT   (:,:)
  REAL(hp), POINTER :: H_JNO2      (:,:)
  REAL(hp), POINTER :: H_JOH       (:,:)
  REAL(hp), POINTER :: H_PSC2_WET  (:,:)
  REAL(hp), POINTER :: H_FRCLND    (:,:)
  REAL(hp), POINTER :: H_MODISLAI  (:,:)
  REAL(hp), POINTER :: H_CNV_FRC   (:,:)
  INTEGER, POINTER  :: H_TropLev   (:,:)

  REAL(hp), POINTER :: H_TROPP     (:,:)
  REAL(hp), POINTER :: H_FLASH_DENS(:,:)
  REAL(hp), POINTER :: H_CONV_DEPTH(:,:)
  REAL(hp), POINTER :: H_SUNCOS    (:,:)

  REAL(hp), POINTER :: H_DRY_TOTN  (:,:)
  REAL(hp), POINTER :: H_WET_TOTN  (:,:)

  ! 3-D fields
  REAL(hp), POINTER :: H_T         (:,:,:)
  REAL(hp), POINTER :: H_AD        (:,:,:)
  REAL(hp), POINTER :: H_AIRVOL    (:,:,:)
  REAL(hp), POINTER :: H_AIRDEN    (:,:,:)
  REAL(hp), POINTER :: H_F_OF_PBL  (:,:,:)

  REAL(hp), POINTER :: H_SPHU      (:,:,:)                ! Note: Only need ZBND = 1 sfc value

  REAL(hp), POINTER :: H_SpcO3     (:,:,:)
  REAL(hp), POINTER :: H_SpcNO2    (:,:,:)
  REAL(hp), POINTER :: H_SpcNO     (:,:,:)
  REAL(hp), POINTER :: H_SpcHNO3   (:,:,:)
  REAL(hp), POINTER :: H_SpcPOPG   (:,:,:)

  ! Intermediate temporaries for regridding
  REAL(hp), POINTER :: REGR_3DI    (:,:,:)                ! Regridding temporary pointer (in)
  REAL(hp), POINTER :: REGR_3DO    (:,:,:)                ! (out)

#endif

!
! !DEFINED PARAMETERS:
!
  ! Temporary toggle for diagnostics
  LOGICAL,  PARAMETER :: DoDiagn = .TRUE.

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_GC_Init
!
! !DESCRIPTION: Subroutine HCOI\_GC\_INIT initializes the HEMCO derived
! types and arrays. The HEMCO configuration is read from the HEMCO
! configuration file (HEMCO_Config.rc) and stored in the HEMCO configuration
! object. The entire HEMCO setup is based upon the entries in the HEMCO
! configuration object. It is possible to explicitly provide a (previously
! read) HEMCO configuration object via input argument `HcoConfig`. In this
! case the HEMCO configuration file will not be read any more.
!\\
!\\
! It is also possible to specify an optional, State_Grid_HCO argument.
! If this is specified, the HEMCO 'intermediate grid' implementation will be
! enabled and HEMCO will operate on a distinct grid from the GEOS-Chem simulation,
! and met fields and emissions will be regridded on-demand in memory. This permits
! the use of higher resolution masks, scaling factors and HEMCO extensions.
! (hplin, 6/2/20)
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

#ifdef MODEL_CLASSIC
    ! HEMCO Intermediate Grid specification
    USE HCO_State_GC_Mod,   ONLY : State_Grid_HCO
    USE GC_Grid_Mod,        ONLY : Compute_Scaled_Grid
    USE HCO_Utilities_GC_Mod, ONLY : Init_IMGrid
#endif

!
! !INPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(IN   )           :: State_Chm  ! Chemistry state
    TYPE(GrdState),   INTENT(IN   )           :: State_Grid ! Grid state
    TYPE(MetState),   INTENT(IN   )           :: State_Met  ! Met state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)           :: Input_Opt  ! Input opts
    TYPE(ConfigObj),  POINTER,      OPTIONAL  :: HcoConfig  ! HEMCO config object
    INTEGER,          INTENT(INOUT)           :: RC         ! Failure or success
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
    CHARACTER(LEN=255)        :: HcoConfigFile
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
    ThisLoc  = ' -> at HCOI_GC_Init (in module GeosCore/hco_interface_gc_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    ! Name of HEMCO configuration file
    HcoConfigFile = 'HEMCO_Config.rc'

    ! Define all species ID's here, for use in module routines below
    id_HNO3  = Ind_('HNO3')
    id_LIMO  = Ind_('LIMO')
    id_NO    = Ind_('NO'  )
    id_NO2   = Ind_('NO2' )
    id_O3    = Ind_('O3'  )
    id_POPG  = Ind_('POPG')
    IF ( id_POPG < 0 ) THEN
       id_POPG  = Ind_('POPG_BaP')
    ENDIF
    IF ( id_POPG < 0 ) THEN
       id_POPG  = Ind_('POPG_PHE')
    ENDIF
    IF ( id_POPG < 0 ) THEN
       id_POPG  = Ind_('POPG_PYR')
    ENDIF

#ifdef MODEL_CLASSIC
    ! Initialize the intermediate grid descriptor.
    ! To disable the HEMCO intermediate grid feature, simply set this DY, DX to
    ! equal to State_Grid%DY, State_Grid%DX (e.g. 2x2.5, 4x5, ...)
    !
    ! TODO: Read in the grid parameters via geoschem_config.yml. For now,
    ! hardcode the scale factor.
    Input_Opt%IMGRID_XSCALE = 1
    Input_Opt%IMGRID_YSCALE = 1

    ! To test GC-Classic WITHOUT intermediate grid
    ! Input_Opt%IMGRID_XSCALE = 1
    ! Input_Opt%IMGRID_YSCALE = 1

    ! Are we using HEMCO intermediate grid implementation?
    ! Note: The mere presence of State_Grid_HCO does not mean
    ! that the intermediate grid is necessarily different.
    ! The code path is decided if the intermediate is actually a different grid.
    IF ( Input_Opt%IMGRID_XSCALE .ne. 1 .or. Input_Opt%IMGRID_XSCALE .ne. 1 ) THEN
      ! Force .or. .true. to waste CPU cycles in regridding and debug Map_A2A above
      Input_Opt%LIMGRID = .true.

      ! In intermediate grid implementation
      WRITE(6, *) "HEMCO INTERMEDIATE GRID:"
      CALL Compute_Scaled_Grid ( Input_Opt, State_Grid, State_Grid_HCO, Input_Opt%IMGRID_XSCALE, Input_Opt%IMGRID_YSCALE, RC )
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in "Compute_Scaled_Grid"!'
         CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
      ENDIF

      write(State_Grid_HCO%GridRes, '(f10.4,a,f10.4)') State_Grid_HCO%DX, "x", State_Grid_HCO%DY

      ! Initialize module temporaries for regridding
      ALLOCATE( REGR_3DI( State_Grid%NX, State_Grid%NY, State_Grid%NZ+1 ), STAT=RC )
      CALL GC_CheckVar( 'hco_interface_gc_mod.F90:HCOI_GC_Init:REGR_3DI', 0, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ALLOCATE( REGR_3DO( State_Grid_HCO%NX, State_Grid_HCO%NY, State_Grid_HCO%NZ+1 ), STAT=RC )
      CALL GC_CheckVar( 'hco_interface_gc_mod.F90:HCOI_GC_Init:REGR_3DO', 0, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! And within the utilities module
      CALL Init_IMGrid( Input_Opt, State_Grid, State_Grid_HCO )
    ELSE
      ! Intermediate grid is same as model grid. Maintain current implementation
      ! all computations about State_Grid_HCO can be skipped to save memory.
    ENDIF
#endif

    ! Create a splash page
    IF ( Input_Opt%amIRoot ) THEN
       WRITE( 6, '(a)' ) REPEAT( '%', 79 )
       WRITE( 6, 100   ) 'HEMCO: Harmonized Emissions Component'
#ifdef MODEL_CLASSIC
       IF ( Input_Opt%LIMGRID ) THEN
         WRITE( 6, '(a)' ) 'HEMCO is running on a different grid than the model', State_Grid_HCO%GridRes
       ENDIF
#endif
       WRITE( 6, 101   ) 'You are using HEMCO version ', ADJUSTL(HCO_VERSION)
       WRITE( 6, '(a)' ) REPEAT( '%', 79 )
 100   FORMAT( '%%%%%', 15x, a,      17x, '%%%%%' )
 101   FORMAT( '%%%%%', 15x, a, a12, 14x, '%%%%%' )
    ENDIF

    !=======================================================================
    ! Read HEMCO configuration file and save into buffer. This also
    ! sets the HEMCO error properties (verbose mode? log file name,
    ! etc.) based upon the specifications in the configuration file.
    ! The log file is now read in two phases: phase 1 reads only the
    ! settings and extensions; phase 2 reads all data fields. This
    ! way, settings and extension options can be updated before
    ! reading all the associated fields.
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

    ! In any way, even if the intermediate grid is enabled, the grid resolution
    ! tokens ($RES$) should be replaced with the token that is most appropriate
    ! with the model native resolution.
    !
    ! This might not be the case in the future, but it involves other
    ! science implications that should be discussed with the scientific GCSC!
    ! (hplin, 6/4/20)
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

       ! Free pointer
       SpcInfo => NULL()
    ENDDO

    !---------------------------------------
    ! Phase 1: read settings and switches
    !---------------------------------------
    CALL Config_ReadFile( Input_Opt%amIRoot,        &
                          iHcoConfig,               &
                          HcoConfigFile,            &
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
                          HcoConfigFile,            &
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
#ifdef MODEL_CLASSIC
    IF ( Input_Opt%LIMGRID ) THEN
      CALL SetHcoGrid( State_Grid_HCO, State_Met, HcoState, RC )
    ELSE
#endif
      CALL SetHcoGrid( State_Grid, State_Met, HcoState, RC )
#ifdef MODEL_CLASSIC
    ENDIF
#endif

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

#ifdef ADJOINT
    if ( Input_Opt%amIRoot ) WRITE(*,*) 'Setting isAdjoint to ', Input_Opt%is_adjoint
    HcoState%isAdjoint = Input_opt%is_adjoint
#endif
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
          ErrMsg = 'DustGinoux is on in HEMCO but activate dust is false ' // &
                   ' in geoschem_config.yml!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
       Input_Opt%LDEAD = .FALSE.
    ENDIF

    ! DEAD dust emissions
    IF ( ExtState%DustDead > 0 ) THEN
       IF ( .not. Input_Opt%LDUST ) THEN
          ErrMsg = 'DustDead is on in HEMCO but activate dust is false ' // &
                   'in geoschem_config.yml!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
       Input_Opt%LDEAD = .TRUE.
    ENDIF

    ! Dust alkalinity
    IF ( ExtState%DustAlk > 0 ) THEN
       IF ( .not. Input_Opt%LDSTUP ) THEN
          ErrMsg = 'DustAlk is on in HEMCO but acid_uptake_on_dust is ' // &
                   'false in geoschem_config.yml'
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
#ifdef MODEL_CLASSIC
    IF ( Input_Opt%LIMGRID ) THEN
      ! Use the HEMCO "intermediate" grid -- allocate arrays accordingly
      CALL ExtState_InitTargets( HcoState, ExtState, Input_Opt, State_Grid, HMRC, State_Grid_HCO )
    ELSE
      CALL ExtState_InitTargets( HcoState, ExtState, Input_Opt, State_Grid, HMRC )
    ENDIF
#else
    CALL ExtState_InitTargets( HcoState, ExtState, Input_Opt, State_Grid, HMRC )
#endif


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
!                  GEOS-Chem Global Chemical Transport Model                  !
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

#if defined( MODEL_CLASSIC )
    ! HEMCO utility routines for GEOS-Chem
    USE HCO_Utilities_GC_Mod, ONLY : Get_GC_Restart
    USE HCO_Utilities_GC_Mod, ONLY : Get_Boundary_Conditions
#endif
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

    INTEGER            :: year, month, day, dayOfYr, hour, minute, second

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
       ' -> at HCOI_GC_Run (in module GeosCore/hco_interface_gc_mod.F90)'
    Instr     = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '      // &
                'HEMCO log file for additional error messages!'
    notDryRun = ( .not. Input_Opt%DryRun )

    !=======================================================================
    ! Make sure HEMCO time is in sync with simulation time
    ! Now done through a universal function in HCO_Interface_Common
    !=======================================================================
    year      = GET_YEAR()
    month     = GET_MONTH()
    day       = GET_DAY()
    dayOfYr   = GET_DAY_OF_YEAR()
    hour      = GET_HOUR()
    minute    = GET_MINUTE()
    second    = GET_SECOND()

    CALL SetHcoTime( HcoState, ExtState, year,   month,     day, dayOfYr, &
                     hour,     minute,   second, EmisTime,  HMRC         )

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
       CALL GridEdge_Set( Input_Opt, State_Grid, State_Met, HcoState, HMRC  )

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "GridEdge_Set"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

#if !defined( ESMF_ ) && !defined( MODEL_WRF )
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

#if defined( MODEL_CLASSIC )
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
          ! Set fields either as pointer targets or else.
          ! Regrid those for ImGrid.
          CALL ExtState_SetFields( Input_Opt, State_Chm, State_Grid, State_Met, &
                                   HcoState,  ExtState,  HMRC )

          ! Trap potential errors
          IF ( HMRC /= HCO_SUCCESS ) THEN
             RC     = HMRC
             ErrMsg = 'Error encountered in "ExtState_SetFields"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
             CALL Flush( HcoState%Config%Err%Lun )
             RETURN
          ENDIF

          ! Update fields directly from State_Met.
          ! Regrid those for ImGrid.
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
    ThisLoc  = ' -> at HCOI_GC_Final (in module GeosCore/hco_interface_gc_mod.F90)'
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

#if defined( MODEL_CLASSIC )
    IF ( ASSOCIATED( REGR_3DI ) ) THEN
       DEALLOCATE( REGR_3DI  )
       CALL GC_CheckVar( 'hco_interface_gc_mod.F90:HCOI_GC_Final:REGR_3DI', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ASSOCIATED( REGR_3DO ) ) THEN
       DEALLOCATE( REGR_3DO  )
       CALL GC_CheckVar( 'hco_interface_gc_mod.F90:HCOI_GC_Final:REGR_3DO', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF
#endif

  END SUBROUTINE HCOI_GC_Final
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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

    USE Time_Mod,        ONLY : Get_Year, Get_Month, Get_Day, GET_DAY_OF_YEAR
    USE Time_Mod,        ONLY : GET_HOUR, GET_MINUTE, GET_SECOND
#if defined( ADJOINT )
    USE MAPL_CommsMod,   ONLY : MAPL_AM_I_ROOT
#endif
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
    INTEGER             :: year, month, day, dayOfYr, hour, minute, second

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
       ' -> at HCOI_GC_WriteDiagn (in module GeosCore/hco_interface_gc_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! Make sure HEMCO time is in sync
    !-----------------------------------------------------------------------

    ! Now done through a universal function in HCO_Interface_Common
    ! (hplin, 3/12/20)
    year      = GET_YEAR()
    month     = GET_MONTH()
    day       = GET_DAY()
    dayOfYr   = GET_DAY_OF_YEAR()
    hour      = GET_HOUR()
    minute    = GET_MINUTE()
    second    = GET_SECOND()

    CALL SetHcoTime( HcoState, ExtState, year,   month,   day, dayOfYr, &
                     hour,     minute,   second, .FALSE., HMRC         )

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
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtState_InitTargets
!
! !DESCRIPTION: SUBROUTINE ExtState\_InitTargets allocates some local arrays
! that act as targets for the ExtState object.
!\\
! If State_Grid_HCO is not given or Input_Opt%LIMGRID is set to false (from Init),
! this explicitly assumes that the HEMCO emissions grid is the same as the GEOS-Chem
! simulation grid.
!
! Otherwise, the met data has to be regridded explicitly at every time step!
! This is performed by creating a set of in-memory array temporaries within this module,
! that can be pointed to by ExtState. Note that a SET of temporaries MUST be always
! present -- you cannot cheat and swap them in memory.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtState_InitTargets( HcoState, ExtState, Input_Opt, State_Grid, RC, &
                                   State_Grid_HCO )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Arr_Mod,    ONLY : HCO_ArrAssert
    USE State_Grid_Mod, ONLY : GrdState
    USE Input_Opt_Mod,  ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN   )           :: Input_Opt      ! Input options
    TYPE(GrdState),   INTENT(IN   )           :: State_Grid     ! Grid State object
    TYPE(GrdState),   INTENT(IN   ), OPTIONAL :: State_Grid_HCO ! HEMCO ImGrid
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

    INTEGER            :: IM, JM, LM
    INTEGER            :: IMh, JMh              ! HEMCO grid sizes

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
       ' -> at ExtState_InitTargets (in module GeosCore/hco_interface_gc_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    ! Shorthands
    IM     = State_Grid%NX
    JM     = State_Grid%NY
    LM     = State_Grid%NZ

    ! By default, HEMCO grid sizes are the same as model. These values will be
    ! overwritten later
    IMh    = IM
    JMh    = JM

#if defined( MODEL_CLASSIC )
    IF ( Input_Opt%LIMGRID ) THEN
      IMh  = State_Grid_HCO%NX
      JMh  = State_Grid_HCO%NY
    ENDIF

    !-----------------------------------------------------------------------
    ! Notes for HEMCO Intermediate Grid:
    ! Due to the extensions requiring met fields and some of them are derived
    ! data, we need to initialize targets here.
    !
    ! Some variables are stored on the MODEL grid, i.e., SZAFACT.
    ! Because they need to be computed here, before a regrid.
    !
    ! Shadow targets for all derived met fields (not read from HEMCO) are
    ! initialized if Input_Opt%LIMGRID (Intermediate enabled) and kept updated
    ! at every ExtState_UpdateFields.
    !-----------------------------------------------------------------------
    IF ( Input_Opt%LIMGRID ) THEN

      ! 2-D Fields
      IF ( ExtState%SZAFACT%DoUse ) THEN
         ALLOCATE( H_SZAFACT( IMh, JMh ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_SZAFACT', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_SZAFACT = 0.0e0_hp
      ENDIF

      IF ( ExtState%JNO2%DoUse ) THEN
         ALLOCATE( H_JNO2( IMh, JMh ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_JNO2', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_JNO2 = 0.0e0_hp
      ENDIF

      IF ( ExtState%JOH%DoUse ) THEN
         ALLOCATE( H_JOH( IMh, JMh ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_JOH', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_JOH = 0.0e0_hp
      ENDIF

      IF ( ExtState%PSC2_WET%DoUse ) THEN
         ALLOCATE( H_PSC2_WET( IMh, JMh ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_PSC2_WET', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_PSC2_WET = 0.0e0_hp
      ENDIF

      IF ( ExtState%FRCLND%DoUse ) THEN
         ALLOCATE( H_FRCLND( IMh, JMh ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_FRCLND', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_FRCLND = 0.0e0_hp
      ENDIF

      IF ( ExtState%LAI%DoUse ) THEN
         ALLOCATE( H_MODISLAI( IMh, JMh ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_MODISLAI', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_MODISLAI = 0.0e0_hp
      ENDIF

      IF ( ExtState%CNV_FRC%DoUse ) THEN
         ALLOCATE( H_CNV_FRC( IMh, JMh ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_CNV_FRC', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_CNV_FRC = 0.0e0_hp
      ENDIF

      IF ( ExtState%TropLev%DoUse ) THEN
         ALLOCATE( H_TropLev( IMh, JMh ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_TropLev', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_TropLev = 0.0e0_hp
      ENDIF

      IF ( ExtState%TROPP%DoUse ) THEN
         ALLOCATE( H_TROPP( IMh, JMh ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_TROPP', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_TROPP = 0.0e0_hp
      ENDIF

      IF ( ExtState%FLASH_DENS%DoUse ) THEN
         ALLOCATE( H_FLASH_DENS( IMh, JMh ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_FLASH_DENS', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_FLASH_DENS = 0.0e0_hp
      ENDIF

      IF ( ExtState%CONV_DEPTH%DoUse ) THEN
         ALLOCATE( H_CONV_DEPTH( IMh, JMh ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_CONV_DEPTH', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_CONV_DEPTH = 0.0e0_hp
      ENDIF

      IF ( ExtState%DRY_TOTN%DoUse ) THEN
         ALLOCATE( H_DRY_TOTN( IMh, JMh ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_DRY_TOTN', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_DRY_TOTN = 0.0e0_hp
      ENDIF

      IF ( ExtState%WET_TOTN%DoUse ) THEN
         ALLOCATE( H_WET_TOTN( IMh, JMh ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_WET_TOTN', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_WET_TOTN = 0.0e0_hp
      ENDIF

      ! Almost 3-D Fields
      IF ( ExtState%SPHU%DoUse ) THEN
         ALLOCATE( H_SPHU( IMh, JMh, 1 ), STAT=RC )         ! Only need SURFACE SPHU
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_SPHU', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_SPHU = 0.0e0_hp
      ENDIF

      ! 3-D Fields
      IF ( ExtState%TK%DoUse ) THEN
         ALLOCATE( H_T( IMh, JMh, LM ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_T', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_T = 0.0e0_hp
      ENDIF

      IF ( ExtState%AIR%DoUse ) THEN
         ALLOCATE( H_AD( IMh, JMh, LM ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_AD', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_AD = 0.0e0_hp
      ENDIF

      IF ( ExtState%AIRVOL%DoUse ) THEN
         ALLOCATE( H_AIRVOL( IMh, JMh, LM ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_AIRVOL', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_AIRVOL = 0.0e0_hp
      ENDIF

      IF ( ExtState%FRAC_OF_PBL%DoUse ) THEN
         ALLOCATE( H_F_OF_PBL( IMh, JMh, LM ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_F_OF_PBL', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_F_OF_PBL = 0.0e0_hp
      ENDIF

      ! 3-D State_Chm fields
      IF ( id_O3 > 0 ) THEN
        ALLOCATE( H_SpcO3( IMh, JMh, LM ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_SpcO3', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_SpcO3 = 0.0e0_hp
      ENDIF

      IF ( id_NO2 > 0 ) THEN
        ALLOCATE( H_SpcNO2( IMh, JMh, LM ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_SpcNO2', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_SpcNO2 = 0.0e0_hp
      ENDIF

      IF ( id_NO > 0 ) THEN
        ALLOCATE( H_SpcNO( IMh, JMh, LM ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_SpcNO', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_SpcNO = 0.0e0_hp
      ENDIF

      IF ( id_HNO3 > 0 ) THEN
        ALLOCATE( H_SpcHNO3( IMh, JMh, LM ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_SpcHNO3', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_SpcHNO3 = 0.0e0_hp
      ENDIF

      IF ( id_POPG > 0 ) THEN
        ALLOCATE( H_SpcPOPG( IMh, JMh, LM ), STAT=RC )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:H_SpcPOPG', 0, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
         H_SpcPOPG = 0.0e0_hp
      ENDIF

    ENDIF
#endif

    ! Initialize max. PBL
    HCO_PBL_MAX = 0

    ! ----------------------------------------------------------------------
    ! Arrays to be copied physically because HEMCO units are not the
    ! same as in GEOS-Chem
    !
    ! In ImGrid, they also need to be regridded on the fly.
    ! ----------------------------------------------------------------------

    ! TROPP: GEOS-Chem TROPP is in hPa, while HEMCO uses Pa.
    IF ( ExtState%TROPP%DoUse ) THEN
       CALL HCO_ArrAssert( ExtState%TROPP%Arr, IMh, JMh, HMRC )

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
       CALL HCO_ArrAssert( ExtState%SPHU%Arr, IMh, JMh, 1, HMRC )

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
       CALL HCO_ArrAssert( ExtState%FLASH_DENS%Arr, IMh, JMh, HMRC )

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
       CALL HCO_ArrAssert( ExtState%CONV_DEPTH%Arr, IMh, JMh, HMRC )

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
       CALL HCO_ArrAssert( ExtState%SUNCOS%Arr, IMh, JMh, HMRC )

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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  SUBROUTINE ExtState_SetFields( Input_Opt, State_Chm, State_Grid, State_Met, HcoState, ExtState, RC )
!
! !USES:
!
    USE Hcox_State_Mod, ONLY : ExtDat_Set
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE Drydep_Mod,     ONLY : DryCoeff
#ifdef ESMF_
    USE HCOI_Esmf_Mod,  ONLY : HCO_SetExtState_ESMF
#endif
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input options
    TYPE(MetState),   INTENT(INOUT)  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm  ! Chemistry state
    TYPE(GrdState),   INTENT(IN   )  :: State_Grid ! Grid state
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
    REAL(hp), POINTER  :: Trgt3D(:,:,:)

    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.
#if defined ( MODEL_WRF )
    LOGICAL, DIMENSION(1:8), SAVE      :: FIRST_PERID = .TRUE.
#endif

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
    Trgt3d   => NULL()
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at ExtState_SetFields (in module GeosCore/hco_interface_gc_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    ! If using intermediate grid (MODEL_CLASSIC and LIMGRID), then load data
    ! from the shadow H_* arrays which have been regridded to HEMCO sizes.
    ! The data is regridded at every call to ExtState_UpdateFields.
    !
    ! Otherwise, simply load as before, either from module variables or
    ! directly from State_Met.
    ! To avoid code duplication, in GEOS-Chem classic data is also loaded
    ! directly from HEMCO met pointers, when possible.

#if defined( MODEL_WRF )
    ! For WRF-GC, the FIRST call needs to be domain-ID specific.
    FIRST = FIRST_PERID(State_Grid%ID)
#endif

    !-----------------------------------------------------------------------
    ! Pointers to local module arrays
    !-----------------------------------------------------------------------

    ! SZAFACT
#if defined( MODEL_CLASSIC )
    IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
      CALL ExtDat_Set( HcoState, ExtState%SZAFACT, 'SZAFACT_FOR_EMIS', &
                       HMRC,     FIRST,            State_Met%SZAFACT )
#if defined( MODEL_CLASSIC )
    ELSE
      CALL ExtDat_Set( HcoState, ExtState%SZAFACT, 'SZAFACT_FOR_EMIS', &
                       HMRC,     FIRST,            H_SZAFACT )
    ENDIF
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( SZAFACT_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! JNO2
#if defined( MODEL_CLASSIC )
    IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
      CALL ExtDat_Set( HcoState, ExtState%JNO2, 'JNO2_FOR_EMIS', &
                       HMRC,     FIRST,         State_Chm%JNO2 )
#if defined( MODEL_CLASSIC )
    ELSE
      CALL ExtDat_Set( HcoState, ExtState%JNO2, 'JNO2_FOR_EMIS', &
                       HMRC,     FIRST,         H_JNO2 )
    ENDIF
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( JNO2_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! JOH
#if defined( MODEL_CLASSIC )
    IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
      CALL ExtDat_Set( HcoState, ExtState%JOH, 'JOH_FOR_EMIS', &
                       HMRC,     FIRST,        State_Chm%JOH )
#if defined( MODEL_CLASSIC )
    ELSE
      CALL ExtDat_Set( HcoState, ExtState%JOH, 'JOH_FOR_EMIS', &
                       HMRC,     FIRST,        H_JOH )
    ENDIF
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( JOH_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! ----------------------------------------------------------------------
    ! 2D fields requiring interpolation (for GC-Classic)
    ! ----------------------------------------------------------------------

    ! PSC2_WET
    ! Computed in calc_met_mod
    IF ( ExtState%PSC2_WET%DoUse ) THEN
#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        CALL ExtDat_Set( HcoState, ExtState%PSC2_WET, 'PSC2_WET_FOR_EMIS', &
                         HMRC,     FIRST,             State_Met%PSC2_WET )
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%PSC2_WET, 'PSC2_WET_FOR_EMIS', &
                         HMRC,     FIRST,             H_PSC2_WET )
      ENDIF
#endif
    ENDIF

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( PSC2_WET_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! FRCLND
    ! Computed in olson_landmap_mod
    IF ( ExtState%FRCLND%DoUse ) THEN
#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        CALL ExtDat_Set( HcoState, ExtState%FRCLND, 'FRCLND_FOR_EMIS', &
                         HMRC,     FIRST,           State_Met%FRCLND )
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%FRCLND, 'FRCLND_FOR_EMIS', &
                         HMRC,     FIRST,           H_FRCLND )
      ENDIF
#endif
    ENDIF

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRCLND_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! LAI
    ! Calculated in modis_lai_mod from XLAI_NATIVE
    IF ( ExtState%LAI%DoUse ) THEN
#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        CALL ExtDat_Set( HcoState, ExtState%LAI, 'LAI_FOR_EMIS', &
                         HMRC,     FIRST,        State_Met%MODISLAI )
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%LAI, 'LAI_FOR_EMIS', &
                         HMRC,     FIRST,        H_MODISLAI )
      ENDIF
#endif
    ENDIF

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( LAI_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! CNV_FRC Convective fractions
    ! Apparently this is not used anywhere right now? maybe it is a GCHP thing
    IF ( ExtState%CNV_FRC%DoUse ) THEN
#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        CALL ExtDat_Set( HcoState, ExtState%CNV_FRC,  'CNV_FRC_FOR_EMIS', &
                         HMRC,     FIRST,             State_Met%CNV_FRC,  &
                         NotFillOk=.TRUE. )
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%CNV_FRC,  'CNV_FRC_FOR_EMIS', &
                         HMRC,     FIRST,             H_CNV_FRC,  &
                         NotFillOk=.TRUE. )
      ENDIF
#endif
    ENDIF

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( CNV_FRC_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Tropopause level
    IF ( ExtState%TropLev%DoUse ) THEN
#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        CALL ExtDat_Set( HcoState, ExtState%TropLev, 'TropLev', &
                         HMRC,     FIRST,            State_Met%TropLev )
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%TropLev, 'TropLev', &
                         HMRC,     FIRST,            H_TropLev )
      ENDIF
#endif
    ENDIF

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( TropLev )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! ----------------------------------------------------------------------
    ! 3D fields requiring interpolation (for GC-Classic)
    ! ----------------------------------------------------------------------

    ! TK
    IF ( ExtState%TK%DoUse ) THEN
#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        CALL ExtDat_Set( HcoState, ExtState%TK, 'TK_FOR_EMIS', &
                         HMRC,     FIRST,       State_Met%T )
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%TK, 'TK_FOR_EMIS', &
                         HMRC,     FIRST,       H_T )
      ENDIF
#endif
    ENDIF

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( TK_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Air mass [kg/grid box]
    IF ( ExtState%AIR%DoUse ) THEN
#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        CALL ExtDat_Set( HcoState, ExtState%AIR, 'AIR_FOR_EMIS', &
                         HMRC,     FIRST,        State_Met%AD )
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%AIR, 'AIR_FOR_EMIS', &
                         HMRC,     FIRST,        H_AD )
      ENDIF
#endif
    ENDIF

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( AIR_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! AIRVOL_FOR_EMIS
    IF ( ExtState%AIRVOL%DoUse ) THEN
#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        CALL ExtDat_Set( HcoState, ExtState%AIRVOL, 'AIRVOL_FOR_EMIS', &
                         HMRC,     FIRST,           State_Met%AIRVOL )
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%AIRVOL, 'AIRVOL_FOR_EMIS', &
                         HMRC,     FIRST,           H_AIRVOL )
      ENDIF
#endif
    ENDIF

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( AIRVOL_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Dry air density [kg/m3]
    IF ( ExtState%AIRDEN%DoUse ) THEN
#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        CALL ExtDat_Set( HcoState, ExtState%AIRDEN, 'AIRDEN', &
                         HMRC,     FIRST,           State_Met%AIRDEN )
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%AIRDEN, 'AIRDEN', &
                         HMRC,     FIRST,           H_AIRDEN )
      ENDIF
#endif
    ENDIF

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( AIRDEN )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Frac of PBL
    IF ( ExtState%FRAC_OF_PBL%DoUse ) THEN
#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        CALL ExtDat_Set( HcoState, ExtState%FRAC_OF_PBL, 'FRAC_OF_PBL_FOR_EMIS', &
                         HMRC,     FIRST,                State_Met%F_OF_PBL )
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%FRAC_OF_PBL, 'FRAC_OF_PBL_FOR_EMIS', &
                         HMRC,     FIRST,                H_F_OF_PBL )
      ENDIF
#endif
    ENDIF

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRAC_OF_PBL_FOR_EMIS)"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! ----------------------------------------------------------------------
    ! 2D fields directly readable from HEMCO (for GC-Classic)
    ! 2D fields (for other models)
    ! ----------------------------------------------------------------------

    ! For MODEL_CLASSIC with the optional HEMCO intermediate grid option,
    ! the meteorological field pointers point to the HEMCO pointers directly.
    ! This avoids an extra regridding in the GC -> HEMCO direction.
    !
    ! For this, simply specify the met field name directly in the FldName arg
    ! call to ExtDat_Set. This will prompt HEMCO to call HCO_EvalFld directly
    ! (hplin, 6/2/20)

    ! U10M
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%U10M, 'U10M', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%U10M, 'U10M_FOR_EMIS', &
                     HMRC,     FIRST,         State_Met%U10M )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( U10M_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! V10M
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%V10M, 'V10M', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%V10M, 'V10M_FOR_EMIS', &
                     HMRC,     FIRST,         State_Met%V10M )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( V10M_FOR_EMIS)"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    ! ALBD
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%ALBD, 'ALBEDO', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%ALBD, 'ALBD_FOR_EMIS', &
                     HMRC,     FIRST,         State_Met%ALBD )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( ALBD_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! T2M
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%T2M, 'T2M', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%T2M, 'T2M_FOR_EMIS', &
                     HMRC,     FIRST,        State_Met%TS )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( T2M_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! TSKIN
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%TSKIN, 'TS', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%TSKIN, 'TSKIN_FOR_EMIS', &
                     HMRC,     FIRST,          State_Met%TSKIN )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( TSKIN_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! GWETROOT
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%GWETROOT, 'GWETROOT', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%GWETROOT, 'GWETROOT_FOR_EMIS', &
                     HMRC,     FIRST,             State_Met%GWETROOT )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( GWETROOT_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! GWETTOP
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%GWETTOP, 'GWETTOP', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%GWETTOP, 'GWETTOP_FOR_EMIS', &
                     HMRC,     FIRST,            State_Met%GWETTOP )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( GWETTOP_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! USTAR
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%USTAR, 'USTAR', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%USTAR, 'USTAR_FOR_EMIS', &
                     HMRC,     FIRST,          State_Met%USTAR )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( USTAR_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! Z0
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%Z0, 'Z0M', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%Z0, 'Z0_FOR_EMIS', &
                     HMRC,     FIRST,       State_Met%Z0 )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( Z0_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! PARDR
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%PARDR, 'PARDR', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%PARDR, 'PARDR_FOR_EMIS', &
                     HMRC,     FIRST,          State_Met%PARDR )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( PARDR_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! PARDF
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%PARDF, 'PARDF', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%PARDF, 'PARDF_FOR_EMIS', &
                     HMRC, FIRST,              State_Met%PARDF )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( PARDF_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! RADSWG
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%RADSWG, 'SWGDN', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%RADSWG, 'RADSWG_FOR_EMIS', &
                     HMRC,     FIRST,           State_Met%SWGDN )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( RADSWG_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! CLDFRC
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%CLDFRC, 'CLDTOT', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%CLDFRC, 'CLDFRC_FOR_EMIS', &
                     HMRC,     FIRST,           State_Met%CLDFRC )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( CLDFRC_FOR_EMIS)"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! SNOWHGT is is mm H2O, which is the same as kg H2O/m2.
    ! This is the unit of SNOMAS.
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%SNOWHGT, 'SNOMAS', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%SNOWHGT, 'SNOWHGT_FOR_EMIS', &
                     HMRC,     FIRST,            State_Met%SNOMAS )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( SNOWHGT_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! SNOWDP [m]
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%SNODP, 'SNODP', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%SNODP, 'SNODP_FOR_EMIS', &
                     HMRC,     FIRST,          State_Met%SNODP )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( SNOWDP_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! FRLAND
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%FRLAND, 'FRLAND', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%FRLAND, 'FRLAND_FOR_EMIS', &
                     HMRC,     FIRST,           State_Met%FRLAND )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRLAND_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! FROCEAN
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%FROCEAN, 'FROCEAN', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%FROCEAN, 'FROCEAN_FOR_EMIS', &
                     HMRC,     FIRST,            State_Met%FROCEAN )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FROCEAN_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

! start blowing snow
    ! FRSEAICE (for blowing snow, huang & jaegle 04/12/20)
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%FRSEAICE, 'FRSEAICE', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%FRSEAICE, 'FRSEAICE_FOR_EMIS', &
                     HMRC,     FIRST,             State_Met%FRSEAICE )
#endif
    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRSEAICE_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! QV2M (for blowing snow, huang & jaegle 04/12/20)
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%QV2M, 'QV2M', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%QV2M, 'QV2M_FOR_EMIS', &
                     HMRC,     FIRST,          State_Met%QV2M )
#endif
    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( QV2M_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF
! end blowing snow

    ! FRLAKE
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%FRLAKE, 'FRLAKE', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%FRLAKE, 'FRLAKE_FOR_EMIS', &
                     HMRC,     FIRST,           State_Met%FRLAKE )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRLAKE_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! FRLANDIC
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState, ExtState%FRLANDIC, 'FRLANDIC', &
                     HMRC,     FIRST=FIRST )
#else
    CALL ExtDat_Set( HcoState, ExtState%FRLANDIC, 'FRLANDIC_FOR_EMIS', &
                     HMRC,     FIRST,             State_Met%FRLANDIC )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( FRLANDIC_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! 3D fields
    !-----------------------------------------------------------------------

    ! CNV_MFC
#if defined( MODEL_CLASSIC )
    CALL ExtDat_Set( HcoState,  ExtState%CNV_MFC, 'CMFMC', &
                     HMRC,      FIRST,            OnLevEdge=.TRUE. )
#else
    CALL ExtDat_Set( HcoState,  ExtState%CNV_MFC, 'CNV_MFC_FOR_EMIS', &
                     HMRC,      FIRST,            State_Met%CMFMC,    &
                     OnLevEdge=.TRUE. )
#endif

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( CNV_MFC_FOR_EMIS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! ----------------------------------------------------------------
    ! Species concentrations
    ! All of these require interpolation on-demand.
    ! ----------------------------------------------------------------

    ! Note: ExtDat_Set points DIRECTLY to the assigned target when received
    ! through ExtDat%Arr%Val => Trgt. This means that the array temporary
    ! must be maintained through time in memory. It may be particularly taxing
    ! for a GC-classic run with intermediate grid option, as all of these must
    ! be maintained in regridded HIGH-RESOLUTION (HEMCO resolution) array
    ! temporaries now! (hplin, 6/2/20)

    IF ( id_O3 > 0 ) THEN

#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        Trgt3D => State_Chm%Species(id_O3)%Conc
        CALL ExtDat_Set( HcoState, ExtState%O3, 'HEMCO_O3_FOR_EMIS', &
                         HMRC,     FIRST,       Trgt3D )

        Trgt3D => NULL()
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%O3, 'HEMCO_O3_FOR_EMIS', &
                         HMRC,     FIRST,       H_SpcO3 )
      ENDIF
#endif

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_O3_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
    ENDIF

    IF ( id_NO2 > 0 ) THEN

#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        Trgt3D => State_Chm%Species(id_NO2)%Conc
        CALL ExtDat_Set( HcoState, ExtState%NO2, 'HEMCO_NO2_FOR_EMIS', &
                         HMRC,     FIRST,        Trgt3D )

        Trgt3D => NULL()
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%NO2, 'HEMCO_NO2_FOR_EMIS', &
                         HMRC,     FIRST,        H_SpcNO2 )
      ENDIF
#endif

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_NO2_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
    ENDIF

    IF ( id_NO > 0 ) THEN

#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        Trgt3D => State_Chm%Species(id_NO)%Conc
        CALL ExtDat_Set( HcoState, ExtState%NO, 'HEMCO_NO_FOR_EMIS', &
                         HMRC,     FIRST,       Trgt3D )

        Trgt3D => NULL()
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%NO, 'HEMCO_NO_FOR_EMIS', &
                         HMRC,     FIRST,       H_SpcNO )
      ENDIF
#endif

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_NO_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
    ENDIF

    IF ( id_HNO3 > 0 ) THEN

#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        Trgt3D => State_Chm%Species(id_HNO3)%Conc
        CALL ExtDat_Set( HcoState, ExtState%HNO3, 'HEMCO_HNO3_FOR_EMIS', &
                         HMRC,     FIRST,         Trgt3D )

        Trgt3D => NULL()
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%HNO3, 'HEMCO_HNO3_FOR_EMIS', &
                         HMRC,     FIRST,         H_SpcHNO3 )
      ENDIF
#endif

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_HNO3_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
    ENDIF

    IF ( id_POPG > 0 ) THEN

#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        Trgt3D => State_Chm%Species(id_POPG)%Conc
        CALL ExtDat_Set( HcoState, ExtState%POPG, 'HEMCO_POPG_FOR_EMIS', &
                         HMRC,     FIRST,         Trgt3D )

        Trgt3D => NULL()
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%POPG, 'HEMCO_POPG_FOR_EMIS', &
                         HMRC,     FIRST,         H_SpcPOPG )
      ENDIF
#endif

       ! Trap potential errors
       IF ( HMRC /= HCO_SUCCESS ) THEN
          RC     = HMRC
          ErrMsg = 'Error encountered in "ExtDat_Set( HEMCO_POPG_FOR_EMIS )"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
          RETURN
       ENDIF
    ENDIF

    ! ----------------------------------------------------------------------
    ! Deposition parameter
    ! ----------------------------------------------------------------------

    ! DRY_TOTN
    IF ( ExtState%DRY_TOTN%DoUse ) THEN
#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        CALL ExtDat_Set( HcoState, ExtState%DRY_TOTN, 'DRY_TOTN_FOR_EMIS', &
                         HMRC,     FIRST,             State_Chm%DryDepNitrogen )
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%DRY_TOTN, 'DRY_TOTN_FOR_EMIS', &
                         HMRC,     FIRST,             H_DRY_TOTN )
      ENDIF
#endif
    ENDIF

    ! Trap potential errors
    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "ExtDat_Set( DRY_TOTN )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF

    ! WET_TOTN
    IF ( ExtState%WET_TOTN%DoUse ) THEN
#if defined( MODEL_CLASSIC )
      IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
        CALL ExtDat_Set( HcoState, ExtState%WET_TOTN, 'WET_TOTN_FOR_EMIS', &
                         HMRC,     FIRST,             State_Chm%WetDepNitrogen )
#if defined( MODEL_CLASSIC )
      ELSE
        CALL ExtDat_Set( HcoState, ExtState%WET_TOTN, 'WET_TOTN_FOR_EMIS', &
                         HMRC,     FIRST,             H_WET_TOTN )
      ENDIF
#endif
    ENDIF


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
          ! Polynomial coefficients for dry deposition, array target in drydep_mod
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
    FIRST  = .FALSE.
#if defined( MODEL_WRF )
    FIRST_PERID(State_Grid%ID) = .FALSE.
#endif
    Trgt3D => NULL()

    ! Leave with success
    RC = GC_SUCCESS

  END SUBROUTINE ExtState_SetFields
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
    USE CMN_FJX_MOD,          ONLY : ZPJ
    USE ErrCode_Mod
    USE FAST_JX_MOD,          ONLY : RXN_NO2, RXN_O3_1
    USE HCO_GeoTools_Mod,     ONLY : HCO_GetSUNCOS
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
#ifdef ESMF_
    USE HCOI_ESMF_MOD,        ONLY : HCO_SetExtState_ESMF
#endif
#if defined( MODEL_CLASSIC )
    USE HCO_State_GC_Mod,     ONLY : State_Grid_HCO    ! HEMCO intermediate grid
    USE HCO_Utilities_GC_Mod, ONLY : Regrid_MDL2HCO
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
       ' -> at ExtState_UpdateFields (in module GeosCore/hco_interface_gc_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !=======================================================================
    ! Update fields in the HEMCO Extension state
    ! Directly from State_Met if not intermediate
    !=======================================================================

    IF ( .not. Input_Opt%LIMGRID ) THEN

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

    ENDIF

    ! SUNCOS
    ! Calculated by HEMCO on its own grid - simply need to allocate
    ! at correct sizes
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

    ! Compute SZAFACT on MODEL GRID
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
          State_Met%SZAFACT(I,J) = GET_SZAFACT(I,J,State_Met)
       ENDIF

       ! Maximum extent of the PBL [model level]
       HCO_PBL_MAX = State_Met%PBL_MAX_L

    ENDDO
    ENDDO
    ENDDO
!$OMP END PARALLEL DO

#if defined ( MODEL_CLASSIC )
    IF ( Input_Opt%LIMGRID ) THEN
    !=======================================================================
    ! GEOS-Chem Classic HEMCO "Intermediate" Grid: Regrid appropriate met fields
    !=======================================================================

    ! Template for 3D Edge
    ! REGR_3DI(:,:,:) = State_Met%PEDGE * 100.0_hp
    ! CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
    !                      REGR_3DI,  REGR_3DO,   ZBND=State_Grid_HCO%NZ+1, &
    !                      ResetRegrName=.true. )
    ! PEDGE(:,:,:)    = REGR_3DO(:,:,:)

    ! Template for 2D
    ! REGR_3DI(:,:,1) = State_Met%PHIS(:,:)
    ! CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
    !                      REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
    !                      ResetRegrName=.true. )
    ! ZSFC(:,:)       = REGR_3DO(:,:,1)

    !-----------------------------------------------------------------------
    ! Local module arrays
    !-----------------------------------------------------------------------

    ! SZAFACT
    IF ( ExtState%SZAFACT%DoUse ) THEN
      REGR_3DI(:,:,1) = State_Met%SZAFACT(:,:)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      H_SZAFACT(:,:)  = REGR_3DO(:,:,1)
    ENDIF

    ! JNO2
    IF ( ExtState%JNO2%DoUse ) THEN
      REGR_3DI(:,:,1) = State_Chm%JNO2(:,:)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      H_JNO2(:,:)     = REGR_3DO(:,:,1)
    ENDIF

    ! JOH
    IF ( ExtState%JOH%DoUse ) THEN
      REGR_3DI(:,:,1) = State_Chm%JOH(:,:)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      H_JOH(:,:)     = REGR_3DO(:,:,1)
    ENDIF

    !-----------------------------------------------------------------------
    ! 2-D State_Met
    !-----------------------------------------------------------------------

    ! PSC2_WET
    IF ( ExtState%PSC2_WET%DoUse ) THEN
      REGR_3DI(:,:,1) = State_Met%PSC2_WET(:,:)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      H_PSC2_WET(:,:) = REGR_3DO(:,:,1)
    ENDIF

    ! FRCLND
    IF ( ExtState%FRCLND%DoUse ) THEN
      REGR_3DI(:,:,1) = State_Met%FRCLND(:,:)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      H_FRCLND(:,:) = REGR_3DO(:,:,1)
    ENDIF

    ! MODISLAI
    IF ( ExtState%LAI%DoUse ) THEN
      REGR_3DI(:,:,1) = State_Met%MODISLAI(:,:)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      H_MODISLAI(:,:) = REGR_3DO(:,:,1)
    ENDIF

    ! CNV_FRC
    IF ( ExtState%CNV_FRC%DoUse ) THEN
      REGR_3DI(:,:,1) = State_Met%FRCLND(:,:)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      H_FRCLND(:,:) = REGR_3DO(:,:,1)
    ENDIF

    ! TropLev
    IF ( ExtState%TropLev%DoUse ) THEN
      REGR_3DI(:,:,1) = State_Met%TropLev(:,:)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      H_TropLev(:,:) = REGR_3DO(:,:,1)
    ENDIF

    ! TROPP - requires conversion
    IF ( ExtState%TROPP%DoUse ) THEN
      REGR_3DI(:,:,1) = State_Met%TROPP(:,:) * 100.0_hp    ! hPa -> Pa
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      H_TROPP(:,:) = REGR_3DO(:,:,1)
    ENDIF

    ! FLASH_DENS
    IF ( ExtState%FLASH_DENS%DoUse ) THEN
      REGR_3DI(:,:,1) = State_Met%FLASH_DENS(:,:)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      H_FLASH_DENS(:,:) = REGR_3DO(:,:,1)
    ENDIF

    ! CONV_DEPTH
    IF ( ExtState%CONV_DEPTH%DoUse ) THEN
      REGR_3DI(:,:,1) = State_Met%CONV_DEPTH(:,:)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      H_CONV_DEPTH(:,:) = REGR_3DO(:,:,1)
    ENDIF

    !-----------------------------------------------------------------------
    ! 3-D State_Met
    !-----------------------------------------------------------------------

    ! SPHU - 3-D up to level 1
    IF ( ExtState%SPHU%DoUse ) THEN
      REGR_3DI(:,:,1) = State_Met%SPHU(:,:,1) / 1000.0_hp
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      H_SPHU(:,:,:) = 0.0_hp
      H_SPHU(:,:,1) = REGR_3DO(:,:,1)
    ENDIF

    ! TK/T
    IF ( ExtState%TK%DoUse ) THEN
      REGR_3DI(:,:,1:State_Grid%NZ) = State_Met%T(:,:,1:State_Grid%NZ)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=State_Grid%NZ,       & ! 3D data
                           ResetRegrName=.true. )
      H_T(:,:,1:State_Grid%NZ) = REGR_3DO(:,:,1:State_Grid%NZ)
    ENDIF

    ! AIR/AD
    IF ( ExtState%AIR%DoUse ) THEN
      REGR_3DI(:,:,1:State_Grid%NZ) = State_Met%AD(:,:,1:State_Grid%NZ)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=State_Grid%NZ,       & ! 3D data
                           ResetRegrName=.true. )
      H_AD(:,:,1:State_Grid%NZ) = REGR_3DO(:,:,1:State_Grid%NZ)
    ENDIF

    ! AIRVOL
    IF ( ExtState%AIRVOL%DoUse ) THEN
      REGR_3DI(:,:,1:State_Grid%NZ) = State_Met%AIRVOL(:,:,1:State_Grid%NZ)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=State_Grid%NZ,       & ! 3D data
                           ResetRegrName=.true. )
      H_AIRVOL(:,:,1:State_Grid%NZ) = REGR_3DO(:,:,1:State_Grid%NZ)
    ENDIF

    ! AIRDEN
    IF ( ExtState%AIRDEN%DoUse ) THEN
      REGR_3DI(:,:,1:State_Grid%NZ) = State_Met%AIRDEN(:,:,1:State_Grid%NZ)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=State_Grid%NZ,       & ! 3D data
                           ResetRegrName=.true. )
      H_AIRDEN(:,:,1:State_Grid%NZ) = REGR_3DO(:,:,1:State_Grid%NZ)
    ENDIF

    ! FRAC_OF_PBL
    IF ( ExtState%FRAC_OF_PBL%DoUse ) THEN
      REGR_3DI(:,:,1:State_Grid%NZ) = State_Met%F_OF_PBL(:,:,1:State_Grid%NZ)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=State_Grid%NZ,       & ! 3D data
                           ResetRegrName=.true. )
      H_F_OF_PBL(:,:,1:State_Grid%NZ) = REGR_3DO(:,:,1:State_Grid%NZ)
    ENDIF

    !-----------------------------------------------------------------------
    ! 3-D State_Chm
    !-----------------------------------------------------------------------

    ! O3
    IF ( id_O3 > 0 ) THEN
      REGR_3DI(:,:,1:State_Grid%NZ) = &
                 State_Chm%Species(id_O3)%Conc(:,:,1:State_Grid%NZ)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=State_Grid%NZ,       & ! 3D data
                           ResetRegrName=.true. )
      H_SpcO3(:,:,1:State_Grid%NZ) = REGR_3DO(:,:,1:State_Grid%NZ)
    ENDIF

    ! NO2
    IF ( id_NO2 > 0 ) THEN
      REGR_3DI(:,:,1:State_Grid%NZ) = &
                 State_Chm%Species(id_NO2)%Conc(:,:,1:State_Grid%NZ)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=State_Grid%NZ,       & ! 3D data
                           ResetRegrName=.true. )
      H_SpcNO2(:,:,1:State_Grid%NZ) = REGR_3DO(:,:,1:State_Grid%NZ)
    ENDIF

    ! NO
    IF ( id_NO > 0 ) THEN
      REGR_3DI(:,:,1:State_Grid%NZ) = &
                 State_Chm%Species(id_NO)%Conc(:,:,1:State_Grid%NZ)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=State_Grid%NZ,       & ! 3D data
                           ResetRegrName=.true. )
      H_SpcNO(:,:,1:State_Grid%NZ) = REGR_3DO(:,:,1:State_Grid%NZ)
    ENDIF

    ! HNO3
    IF ( id_HNO3 > 0 ) THEN
      REGR_3DI(:,:,1:State_Grid%NZ) = &
                 State_Chm%Species(id_HNO3)%Conc(:,:,1:State_Grid%NZ)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=State_Grid%NZ,       & ! 3D data
                           ResetRegrName=.true. )
      H_SpcHNO3(:,:,1:State_Grid%NZ) = REGR_3DO(:,:,1:State_Grid%NZ)
    ENDIF

    ! POPG
    IF ( id_POPG > 0 ) THEN
      REGR_3DI(:,:,1:State_Grid%NZ) = &
                 State_Chm%Species(id_POPg)%Conc(:,:,1:State_Grid%NZ)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=State_Grid%NZ,       & ! 3D data
                           ResetRegrName=.true. )
      H_SpcPOPG(:,:,1:State_Grid%NZ) = REGR_3DO(:,:,1:State_Grid%NZ)
    ENDIF

    ! DRY_TOTN
    IF ( ExtState%DRY_TOTN%DoUse ) THEN
      REGR_3DI(:,:,1) = State_Chm%DryDepNitrogen(:,:)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      H_DRY_TOTN(:,:) = REGR_3DO(:,:,1)
    ENDIF

    IF ( ExtState%WET_TOTN%DoUse ) THEN
      REGR_3DI(:,:,1) = State_Chm%WetDepNitrogen(:,:)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      H_WET_TOTN(:,:) = REGR_3DO(:,:,1)
    ENDIF


    ENDIF
#endif

  END SUBROUTINE ExtState_UpdateFields
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  SUBROUTINE GridEdge_Set( Input_Opt, State_Grid, State_Met, HcoState, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCOX_STATE_MOD,   ONLY : ExtDat_Set
    USE HCO_GeoTools_Mod, ONLY : HCO_CalcVertGrid
    USE HCO_GeoTools_Mod, ONLY : HCO_SetPBLm
    USE State_Grid_Mod,   ONLY : GrdState
    USE State_Met_Mod,    ONLY : MetState
    USE Input_Opt_Mod,    ONLY : OptInput
#if defined( MODEL_CLASSIC )
    USE HCO_State_GC_Mod, ONLY : State_Grid_HCO    ! HEMCO intermediate grid
    USE HCO_Utilities_GC_Mod, ONLY : Regrid_MDL2HCO
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input options
    TYPE(GrdState),   INTENT(IN   )  :: State_Grid ! Grid state
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(HCO_STATE),  POINTER        :: HcoState   ! HEMCO state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Success or failure?
!
! !REMARKS:
!  GridEdge_Set defines the HEMCO vertical grid used in GEOS-Chem "classic"
!  and WRF-GC simulations.  (GCHP and HEMCO-CESM-GC uses their own interface to HEMCO.)
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
       ' -> at GridEdge_Set (in module GeosCore/hco_interface_gc_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! Allocate all arrays.
    !-----------------------------------------------------------------------
#if defined( MODEL_CLASSIC )
    IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
      ! NOTE: Hco_CalcVertGrid expects pointer-based arguments, so we must
      ! make PEDGE be a pointer and allocate/deallocate it on each call.
      ALLOCATE( PEDGE( State_Grid%NX, State_Grid%NY, State_Grid%NZ+1 ), STAT=RC )
      CALL GC_CheckVar( 'hco_interface_gc_mod.F90:GridEdge_Set:PEDGE', 0, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Edge and surface pressures [Pa]
      PEDGE    =  State_Met%PEDGE * 100.0_hp  ! Convert hPa -> Pa
      PSFC     => PEDGE(:,:,1)

      ! Point to other fields of State_Met
      ZSFC     => State_Met%PHIS
      BXHEIGHT => State_Met%BXHEIGHT
      TK       => State_Met%T

      ALLOCATE( PBLM( State_Grid%NX, State_Grid%NY ), STAT=RC )
      CALL GC_CheckVar( 'hco_interface_gc_mod.F90:GridEdge_Set:PBLM', 0, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN
#if defined( MODEL_CLASSIC )
    ELSE
      ! If intermediate grid, allocate and manually regrid meteorological fields.
      ALLOCATE( PEDGE( State_Grid_HCO%NX, State_Grid_HCO%NY, State_Grid_HCO%NZ+1 ), STAT=RC )
      CALL GC_CheckVar( 'hco_interface_gc_mod.F90:GridEdge_Set:IMG_PEDGE', 0, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ALLOCATE( ZSFC( State_Grid_HCO%NX, State_Grid_HCO%NY ), STAT=RC )
      CALL GC_CheckVar( 'hco_interface_gc_mod.F90:GridEdge_Set:IMG_PHIS', 0, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ALLOCATE( BXHEIGHT( State_Grid_HCO%NX, State_Grid_HCO%NY, State_Grid_HCO%NZ ), STAT=RC )
      CALL GC_CheckVar( 'hco_interface_gc_mod.F90:GridEdge_Set:IMG_BXHEIGHT', 0, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ALLOCATE( TK( State_Grid_HCO%NX, State_Grid_HCO%NY, State_Grid_HCO%NZ ), STAT=RC )
      CALL GC_CheckVar( 'hco_interface_gc_mod.F90:GridEdge_Set:IMG_TK', 0, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ALLOCATE( PBLM( State_Grid_HCO%NX, State_Grid_HCO%NY ), STAT=RC )
      CALL GC_CheckVar( 'hco_interface_gc_mod.F90:GridEdge_Set:IMG_PBLM', 0, RC )
      IF ( RC /= HCO_SUCCESS ) RETURN

      ! Fill and respectively regrid to targets.
      ! IO, SG, SGH, PtrIn, PtrOut, ZBND, ResetRegrName=.true.

      ! Edge pressures [Pa] (hPa->Pa convert)
      REGR_3DI(:,:,:) = State_Met%PEDGE * 100.0_hp
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=State_Grid_HCO%NZ+1, &
                           ResetRegrName=.true. )
      PEDGE(:,:,:)    = REGR_3DO(:,:,:)

      ! Surface pressure
      PSFC            => PEDGE(:,:,1)

      ! ZSFC
      REGR_3DI(:,:,1) = State_Met%PHIS(:,:)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      ZSFC(:,:)       = REGR_3DO(:,:,1)

      ! BXHEIGHT
      REGR_3DI(:,:,1:State_Grid%NZ) = State_Met%BXHEIGHT(:,:,1:State_Grid%NZ)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=State_Grid_HCO%NZ,   & ! 2D data
                           ResetRegrName=.true. )
      BXHEIGHT(:,:,1:State_Grid%NZ) = REGR_3DO(:,:,1:State_Grid%NZ)

      ! TK
      REGR_3DI(:,:,1:State_Grid%NZ) = State_Met%T(:,:,1:State_Grid%NZ)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=State_Grid_HCO%NZ,   & ! 2D data
                           ResetRegrName=.true. )
      TK(:,:,1:State_Grid%NZ)       = REGR_3DO(:,:,1:State_Grid%NZ)
    ENDIF
#endif


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
#if defined( MODEL_CLASSIC )
    IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
!$OMP PARALLEL DO                                                 &
!$OMP DEFAULT( SHARED )                                           &
!$OMP PRIVATE( I, J )
      DO J=1,State_Grid%NY
      DO I=1,State_Grid%NX
         PBLM(I,J) = State_Met%PBL_TOP_m(I,J)
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
#if defined( MODEL_CLASSIC )
    ELSE
      ! Intermediate grid
      ! PBLM
      REGR_3DI(:,:,1) = State_Met%PBL_TOP_m(:,:)
      CALL Regrid_MDL2HCO( Input_Opt, State_Grid, State_Grid_HCO,           &
                           REGR_3DI,  REGR_3DO,   ZBND=1,                   & ! 2D data
                           ResetRegrName=.true. )
      PBLM(:,:)       = REGR_3DO(:,:,1)
    ENDIF
#endif

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
       CALL GC_CheckVar( 'hco_interface_gc_mod.F90:GridEdge_Set:PEDGE', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    ! Deallocate PBLM
    IF ( ASSOCIATED( PBLM ) ) THEN
       DEALLOCATE( PBLM  )
       CALL GC_CheckVar( 'hco_interface_gc_mod.F90:GridEdge_Set:PBLM', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

#if defined( MODEL_CLASSIC )
    IF ( .not. Input_Opt%LIMGRID ) THEN
#endif
      ! Free pointers
      ZSFC     => NULL()
      BXHEIGHT => NULL()
      TK       => NULL()
      PSFC     => NULL()
      PEDGE    => NULL()
      PBLM     => NULL()
#if defined( MODEL_CLASSIC )
    ELSE
      ! Deallocate arrays
      IF ( ASSOCIATED( ZSFC ) ) THEN
         DEALLOCATE( ZSFC  )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:GridEdge_Set:ZSFC', 2, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
      ENDIF

      IF ( ASSOCIATED( BXHEIGHT ) ) THEN
         DEALLOCATE( BXHEIGHT  )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:GridEdge_Set:BXHEIGHT', 2, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
      ENDIF

      IF ( ASSOCIATED( TK ) ) THEN
         DEALLOCATE( TK  )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:GridEdge_Set:TK', 2, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
      ENDIF

      PSFC     => NULL()

      IF ( ASSOCIATED( PEDGE ) ) THEN
         DEALLOCATE( PEDGE  )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:GridEdge_Set:PEDGE', 2, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
      ENDIF

      IF ( ASSOCIATED( PBLM ) ) THEN
         DEALLOCATE( PBLM  )
         CALL GC_CheckVar( 'hco_interface_gc_mod.F90:GridEdge_Set:PBLM', 2, RC )
         IF ( RC /= GC_SUCCESS ) RETURN
      ENDIF
    ENDIF
#endif

  END SUBROUTINE GridEdge_Set
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
       ' -> at SetHcoSpecies (in module GeosCore/hco_interface_gc_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------
    ! For most simulations (e.g. full-chem simulation, most of the
    ! specialty sims), just use the GEOS-Chem species definitions.
    !-----------------------------------------------------------------
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM     .or.                               &
         Input_Opt%ITS_AN_AEROSOL_SIM     .or.                               &
         Input_Opt%ITS_A_CO2_SIM          .or.                               &
         Input_Opt%ITS_A_CH4_SIM          .or.                               &
         Input_Opt%ITS_A_MERCURY_SIM      .or.                               &
         Input_Opt%ITS_A_POPS_SIM         .or.                               &
         Input_Opt%ITS_A_RnPbBe_SIM       .or.                               &
         Input_Opt%ITS_A_TAGO3_SIM        .or.                               &
         Input_Opt%ITS_A_TAGCO_SIM        .or.                               &
         Input_Opt%ITS_A_CARBON_SIM       .or.                               &
         Input_Opt%ITS_A_TRACEMETAL_SIM ) THEN


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

       !%%%%% FOR THE CARBON OR TAGGED CO SIMULATIONS %%%%%
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
       ' -> at SetHcoGrid (in module GeosCore/hco_interface_gc_mod.F90)'
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
    CALL GC_CheckVar( 'hco_interface_gc_mod:SetHcoGrid:Ap', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Allocate Bp array
    ALLOCATE( Bp( State_Grid%NZ+1 ), STAT=RC )
    CALL GC_CheckVar( 'hco_interface_gc_mod:SetHcoGrid:Bp', 0, RC )
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
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CheckSettings
!
! !DESCRIPTION: Subroutine CheckSettings performs some sanity checks of the
! switches provided in the HEMCO configuration file in combination with the
! settings specified in geoschem_config.yml.
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
       ' -> at CheckSettings (in module GeosCore/hco_interface_gc_mod.F90)'
    Instr    = 'THIS ERROR ORIGINATED IN HEMCO!  Please check the '       // &
               'HEMCO log file for additional error messages!'

    !-----------------------------------------------------------------------
    ! If chemistry is turned off, do not read chemistry input data
    !-----------------------------------------------------------------------
    IF ( .NOT. Input_Opt%LCHEM ) THEN

       IF ( Input_Opt%amIRoot ) THEN
          Print*, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
          Print*, '% WARNING: Activate chemistry is set to false in        %'
          Print*, '% geoschem_config.yml so chemistry data will not be     %'
          Print*, '% read by HEMCO(hco_interface_gc_mod.F90)               %'
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
    ! EMISSIONS switch in HEMCO_Config.rc
    !
    ! Create a shadow field (Input_Opt%DoEmissions) to determine if
    ! emissions fluxes should be applied in mixing_mod.F90
    !-----------------------------------------------------------------------
    CALL GetExtOpt( HcoConfig, -999, 'EMISSIONS',           &
                    OptValBool=LTMP, FOUND=FOUND,  RC=HMRC )

    IF ( HMRC /= HCO_SUCCESS ) THEN
       RC     = HMRC
       ErrMsg = 'Error encountered in "GetExtOpt( EMISSIONS )"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
       RETURN
    ENDIF
    IF ( .not. FOUND ) THEN
       ErrMsg = 'EMISSIONS not found in HEMCO_Config.rc file!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    Input_Opt%DoEmissions = LTMP

    !-----------------------------------------------------------------------
    ! Lightning NOx extension
    !
    ! The lightning NOx extension is only used in fullchem simulations. We
    ! will create a shadow field (Input_Opt%DoLightningNOx) to determine if
    ! the FLASH_DENS and CONV_DEPTH fields are needed in flexgrid_read_mod.F90
    !-----------------------------------------------------------------------
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
       ExtNr = GetExtNr( HcoConfig%ExtList, 'LightNOx' )
       IF ( ExtNr <= 0 ) THEN
          Input_Opt%DoLightNOx = .FALSE.
       ELSE
          Input_Opt%DoLightNOx = .TRUE.
       ENDIF
    ELSE
       Input_Opt%DoLightNOx = .FALSE.
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
    ! geoschem_config.yml file, then we will also toggle the OCEAN_Hg
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
                   'but use_dynamic_ocean_Hg is true in geoschem_config.yml.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    !-----------------------------------------------------------------------
    ! Input data for CH4 simulations only
    !
    ! If we have turned on CH4 options in geoschem_config.yml, then we
    ! also need to toggle switches so that HEMCO reads the appropriate data.
    !-----------------------------------------------------------------------
    IF ( Input_Opt%ITS_A_CH4_SIM ) THEN

       IF ( Input_Opt%AnalyticalInv ) THEN
          CALL GetExtOpt( HcoConfig, -999, 'AnalyticalInv', &
                          OptValBool=LTMP, FOUND=FOUND,  RC=HMRC )

          IF ( HMRC /= HCO_SUCCESS ) THEN
             RC     = HMRC
             ErrMsg = 'Error encountered in "GetExtOpt( AnalyticalInv )"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
             RETURN
          ENDIF
          IF ( .not. FOUND ) THEN
             ErrMsg = 'AnalyticalInv not found in HEMCO_Config.rc file!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
          IF ( .not. LTMP ) THEN
             ErrMsg = 'AnalyticalInv is set to false in HEMCO_Config.rc ' // &
                  'but should be set to true for this simulation.'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       IF ( Input_Opt%UseEmisSF ) THEN
          CALL GetExtOpt( HcoConfig, -999, 'Emis_ScaleFactor', &
                          OptValBool=LTMP, FOUND=FOUND,  RC=HMRC )

          IF ( HMRC /= HCO_SUCCESS ) THEN
             RC     = HMRC
             ErrMsg = 'Error encountered in "GetExtOpt( Emis_ScaleFactor )"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
             RETURN
          ENDIF
          IF ( .not. FOUND ) THEN
             ErrMsg = 'Emis_ScaleFactor not found in HEMCO_Config.rc file!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
          IF ( .not. LTMP ) THEN
             ErrMsg = 'Emis_ScaleFactor is set to false in HEMCO_Config.rc '// &
                  'but should be set to true for this simulation.'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       IF ( Input_Opt%UseOHSF ) THEN
          CALL GetExtOpt( HcoConfig, -999, 'OH_ScaleFactor',           &
                          OptValBool=LTMP, FOUND=FOUND,  RC=HMRC )

          IF ( HMRC /= HCO_SUCCESS ) THEN
             RC     = HMRC
             ErrMsg = 'Error encountered in "GetExtOpt( OH_ScaleFactor )"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc, Instr )
             RETURN
          ENDIF
          IF ( .not. FOUND ) THEN
             ErrMsg = 'OH_ScaleFactor not found in HEMCO_Config.rc file!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
          IF ( .not. LTMP ) THEN
             ErrMsg = 'OH_ScaleFactor is set to false in HEMCO_Config.rc ' // &
                  'but should be set to true for this simulation.'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

    ENDIF

    !-----------------------------------------------------------------------
    ! RRTMG input data
    !
    ! If we have turned on the RRTMG simulation in the
    ! geoschem_config.yml file, then we will also toggle the RRTMG
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
                   'but should be set to true for this simulation.'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    ! Print value of shadow fields
    IF ( Input_Opt%amIRoot ) THEN
       Print*, ''
       Print*, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
       Print*, 'Switches read from HEMCO_Config.rc:'
       Print*, '  EMISSIONS : ', Input_Opt%DoEmissions
       Print*, '  LightNOx  : ', Input_Opt%DoLightNOx
       Print*, '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
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
         State_Met%SUNCOSsum(I,J) > 0e+0_fp ) THEN

       ! Impose a diurnal variation on OH during the day
       FACT = ( State_Met%SUNCOS(I,J) / State_Met%SUNCOSsum(I,J) ) &
              * ( 86400e+0_fp / GET_TS_CHEM() )

    ELSE

       ! At night, OH goes to zero
       FACT = 0e+0_fp

    ENDIF

  END FUNCTION Get_SzaFact
!EOC
#if defined ( MODEL_CLASSIC )
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
   USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_GetPtr
   USE Input_Opt_Mod,        ONLY : OptInput
   USE Pressure_Mod,         ONLY : Set_Floating_Pressures
   USE State_Chm_Mod,        ONLY : ChmState
   USE State_Grid_Mod,       ONLY : GrdState
   USE State_Met_Mod,        ONLY : MetState
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

      ! Get delta pressure per grid box stored in restart file to allow
      ! mass conservation across consecutive runs.
      v_name = 'DELPDRY'
      CALL HCO_GC_GetPtr( Input_Opt, State_Grid, TRIM(v_name), Ptr3D, RC, FOUND=FOUND )
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
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute_Sflx_for_Vdiff
!
! !DESCRIPTION: Computes the surface flux (\= emissions - drydep) for the
!  non-local PBL mixing.  This code was removed from within the non-local
!  PBL mixing driver routine VDIFFDR.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Sflx_for_Vdiff( Input_Opt,  State_Chm, State_Diag,      &
                                     State_Grid, State_Met, RC              )
!
! ! USES:
!
    USE Depo_Mercury_Mod,     ONLY : Add_Hg2_DD
    USE Depo_Mercury_Mod,     ONLY : Add_HgP_DD
    USE Depo_Mercury_Mod,     ONLY : Add_Hg2_SnowPack
    USE ErrCode_Mod
    USE Get_Ndep_Mod,         ONLY : Soil_Drydep
    USE HCO_Utilities_GC_Mod, ONLY : GetHcoValEmis, GetHcoValDep, InquireHco
    USE HCO_Utilities_GC_Mod, ONLY : LoadHcoValEmis, LoadHcoValDep
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_GetDiagn
    USE HCO_State_GC_Mod,     ONLY : ExtState
    USE HCO_State_GC_Mod,     ONLY : HcoState
    USE Input_Opt_Mod,        ONLY : OptInput
#if !defined( MODEL_CESM )
    USE Mercury_Mod,          ONLY : Hg_Emis
#endif
    USE PhysConstants
    USE Species_Mod,          ONLY : Species,  SpcConc
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Chm_Mod,        ONLY : Ind_
    USE State_Diag_Mod,       ONLY : DgnState, DgnMap
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE Time_Mod,             ONLY : Get_Ts_Conv
    USE Time_Mod,             ONLY : Get_Ts_Emis
    USE UnitConv_Mod,         ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(GrdState),   INTENT(IN)    :: State_Grid  ! Grid State
    TYPE(MetState),   INTENT(IN)    :: State_Met   ! Meteorology State
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State
    TYPE(DgnState),   INTENT(INOUT) :: State_Diag  ! Diagnostics State
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  The loop order of this routine was changed from I,J,N to N,I,J. This allows
!  us to keep one species in memory at a time, to avoid regridding over and over
!  again when the HEMCO grid is not the same as the model grid.
!  The on-demand regridder caches the LAST SPECIES information, so the outermost
!  loop should be based on the species ID.
!
! !REVISION HISTORY:
!  18 May 2020 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    INTEGER, SAVE           :: id_O3    = -1
    INTEGER, SAVE           :: id_HNO3  = -1

    ! Scalars
    LOGICAL                 :: found,   zeroHg0Dep
    INTEGER                 :: I,       J
    INTEGER                 :: L,       NA
    INTEGER                 :: ND,      N
    INTEGER                 :: Hg_Cat,  topMix
    INTEGER                 :: S
    REAL(fp)                :: dep,     emis
    REAL(fp)                :: MW_kg,   fracNoHg0Dep
    REAL(fp)                :: tmpFlx

    LOGICAL                 :: EmisSpec, DepSpec

    ! Strings
    CHARACTER(LEN=63)       :: origUnit
    CHARACTER(LEN=255)      :: errMsg,  thisLoc

    ! Arrays
    REAL(fp), TARGET        :: eflx(State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Chm%nAdvect                       )
    REAL(fp), TARGET        :: colEflx(State_Grid%NX,                        &
                                       State_Grid%NY,                        &
                                       State_Chm%nAdvect                    )
    REAL(fp), TARGET        :: dflx(State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Chm%nAdvect                       )

    ! Pointers and Objects
    REAL(f4),       POINTER :: Ptr2D(:,:) => NULL()

    REAL(f4),       POINTER :: PNOxLoss_O3(:,:)
    REAL(f4),       POINTER :: PNOxLoss_HNO3(:,:)

    TYPE(Species),  POINTER :: ThisSpc
    TYPE(DgnMap),   POINTER :: mapData

    !=======================================================================
    ! Compute_Sflx_For_Vdiff begins here!
    !=======================================================================

    ! Initialize
    RC      =  GC_SUCCESS
    dflx    =  0.0_fp
    eflx    =  0.0_fp
    colEflx =  0.0_fp
    ThisSpc => NULL()
    errMsg  = ''
    thisLoc = &
    ' -> at Compute_Sflx_for_Vdiff (in module GeosCore/hco_interface_gc_mod.F90)'

    PNOXLoss_HNO3 => NULL()
    PNOxLoss_O3   => NULL()

    ! Reset DryDepMix diagnostic so as not to accumulate from prior timesteps
    IF ( State_Diag%Archive_DryDepMix .or. State_Diag%Archive_DryDep ) THEN
       State_Diag%DryDepMix = 0.0_f4
    ENDIF

    !=======================================================================
    ! Convert units to v/v dry
    !=======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met,     &
                            'v/v dry', RC,        OrigUnit=OrigUnit         )

    ! Trap potential error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountred in "Convert_Spc_Units" (to v/v dry)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Get pointers to the PARANOX loss fluxes.
    ! These are stored in diagnostics 'PARANOX_O3_DEPOSITION_FLUX' and
    ! 'PARANOX_HNO3_DEPOSITION_FLUX'. The call below links pointers
    ! PNOXLOSS_O3 and PNOXLOSS_HNO3 to the data values stored in the
    ! respective diagnostics. The pointers will remain unassociated if
    ! the diagnostics do not exist.
    !
    ! The arrays are now allocated and copied to to accommodate for the
    ! HEMCO intermediate grid feature, as we want the regridded data to
    ! be kept for the rest of this subroutine call.
    !=======================================================================

    ! Get species IDs
    id_O3   = Ind_('O3'  )
    id_HNO3 = Ind_('HNO3')

#if !defined( MODEL_CESM )
    IF ( id_O3 > 0 ) THEN
       CALL HCO_GC_GetDiagn(                                              &
            Input_Opt,  State_Grid,                                       &
            DiagnName      = 'PARANOX_O3_DEPOSITION_FLUX',                &
            StopIfNotFound = .FALSE.,                                     &
            Ptr2D          = Ptr2D,                                       &
            RC             = RC                                          )
    ENDIF
    IF( ASSOCIATED( Ptr2D )) THEN
       ALLOCATE ( PNOxLoss_O3( State_Grid%NX, State_Grid%NY ), STAT=RC )
       PNOxLoss_O3(:,:) = Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()

    IF ( id_HNO3 > 0 ) THEN
       CALL HCO_GC_GetDiagn(                                              &
            Input_Opt,  State_Grid,                                       &
            DiagnName      = 'PARANOX_HNO3_DEPOSITION_FLUX',              &
            StopIfNotFound = .FALSE.,                                     &
            Ptr2D          = Ptr2D,                                       &
            RC             = RC                                          )
    ENDIF
    IF( ASSOCIATED( Ptr2D )) THEN
       ALLOCATE ( PNOxLoss_HNO3( State_Grid%NX, State_Grid%NY ), STAT=RC )
       PNOxLoss_HNO3(:,:) = Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()
#endif

    !=======================================================================
    ! Add emissions & deposition values calculated in HEMCO.
    ! Here we only consider emissions below the PBL top.
    !
    ! The loop has been separated to go over N, J, I in order to optimize
    ! for retrieving regridded data from HEMCO. There is only one regrid
    ! buffer per field (emis/dep, N) so these must not be intertwined,
    ! or there will be a large performance penalty.
    !=======================================================================

    ! Advected species loop
    DO NA = 1, State_Chm%nAdvect

      ! Get the modelId
      N = State_Chm%Map_Advect(NA)

      ! Point to the corresponding entry in the species database
      ThisSpc => State_Chm%SpcData(N)%Info

      ! Check if there is emissions or deposition for this species
#if !defined( MODEL_CESM )
      CALL InquireHco ( N, Emis=EmisSpec, Dep=DepSpec )
#else
      ! Do not apply for MODEL_CESM as its handled by HEMCO-CESM independently
      EmisSpec = .False.
      DepSpec  = .False.
#endif

      ! If there is emissions for this species, it must be loaded into
      ! memory first.   This is achieved by attempting to retrieve a
      ! grid box while NOT in a parallel loop. Failure to load this will
      ! result in severe performance issues!! (hplin, 9/27/20)
      IF ( EmisSpec ) THEN
         CALL LoadHcoValEmis ( Input_Opt, State_Grid, NA )
      ENDIF

      IF ( DepSpec ) THEN
         CALL LoadHcoValDep ( Input_Opt, State_Grid, NA )
      ENDIF

      !$OMP PARALLEL DO                                                      &
      !$OMP DEFAULT( SHARED )                                                &
      !$OMP PRIVATE( I,       J,       topMix                               )&
      !$OMP PRIVATE( tmpFlx,  found,   emis,      dep                       )
      DO J = 1, State_Grid%NY
      DO I = 1, State_Grid%NX

      ! Below emissions. Do not apply for MODEL_CESM as its handled by
      ! HEMCO-CESM independently
#ifndef MODEL_CESM

        ! PBL top level [integral model levels]
        topMix = MAX( 1, FLOOR( State_Met%PBL_TOP_L(I,J) ) )

        !------------------------------------------------------------------
        ! Add total emissions in the PBL to the EFLX array
        ! which tracks emission fluxes.  Units are [kg/m2/s].
        !------------------------------------------------------------------
        IF ( Input_Opt%ITS_A_CH4_SIM ) THEN

           ! CH4 emissions become stored in state_chm_mod.F90.
           ! We use CH4_EMIS here instead of the HEMCO internal emissions
           ! only to make sure that total CH4 emissions are properly defined
           ! in a multi-tracer CH4 simulation. For a single-tracer simulation
           ! and/or all other source types, we could use the HEMCO internal
           ! values set above and would not need the code below.
           ! Units are already in kg/m2/s. (ckeller, 10/21/2014)
           !
           !%%% NOTE: MAYBE THIS CAN BE REMOVED SOON (bmy, 5/18/19)%%%
           eflx(I,J,NA) = State_Chm%CH4_EMIS(I,J,NA)

        ELSE IF ( EmisSpec ) THEN  ! Are there emissions for these species?

           ! Compute emissions for all other simulation
           tmpFlx = 0.0_fp
           DO L = 1, topMix
              CALL GetHcoValEmis( Input_Opt, State_Grid, NA,    I,           &
                                  J,         L,          found, emis        )
              IF ( .NOT. found ) EXIT
              tmpFlx = tmpFlx + emis
           ENDDO
           eflx(I,J,NA) = eflx(I,J,NA) + tmpFlx

           ! Compute column emission fluxes for satellite diagnostics
           IF ( State_Diag%Archive_SatDiagnColEmis ) THEN
              tmpFlx = 0.0_fp
              DO L = 1, State_Grid%NZ
                 CALL GetHcoValEmis( Input_Opt, State_Grid, NA,    I,        &
                                     J,         L,          found, emis     )
                 IF ( .NOT. found ) EXIT
                 tmpFlx = tmpFlx + emis
              ENDDO
              colEflx(I,J,NA) = colEflx(I,J,NA) + tmpFlx
           ENDIF

        ENDIF

        ! For Hg simulations, also add Hg emissions not handled by HEMCO
        IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN
           eflx(I,J,NA) = eflx(I,J,NA) + Hg_EMIS(I,J,NA)
        ENDIF
#endif

        !------------------------------------------------------------------
        ! Also add drydep frequencies calculated by HEMCO (e.g. from the
        ! air-sea exchange module) to DFLX.  These values are stored
        ! in 1/s.  They are added in the same manner as the DEPSAV values
        ! from drydep_mod.F90.  DFLX will be converted to kg/m2/s later.
        ! (ckeller, 04/01/2014)
        !------------------------------------------------------------------
        IF ( DepSpec ) THEN
          CALL GetHcoValDep( Input_Opt, State_Grid, NA, I, J, L, found, dep )
          IF ( found ) THEN
             dflx(I,J,NA) = dflx(I,J,NA) + ( dep                   &
                            * State_Chm%Species(NA)%Conc(I,J,1) &
                            / (AIRMW / ThisSpc%MW_g)  )
          ENDIF
        ENDIF
      ENDDO ! I
      ENDDO ! J
      !$OMP END PARALLEL DO

      ! Free pointers
      ThisSpc => NULL()
    ENDDO   ! NA

    !=======================================================================
    ! Add emissions & deposition values calculated in HEMCO.
    ! Here we only consider emissions below the PBL top.
    !
    ! For the full-chemistry simulations, emissions above the PBL
    ! top will be applied in routine SETEMIS, which occurs just
    ! before the SMVGEAR/KPP solvers are invoked.
    !
    ! For the specialty simulations, emissions above the PBL top
    ! will be applied in the chemistry routines for each
    ! specialty simulation.
    !
    ! For more information, please see this wiki page:
    ! http://wiki.geos-chem.org/Distributing_emissions_in_the_PBL
    !========================================================================
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED )                                                  &
    !$OMP PRIVATE( I,       J,            N                                 )&
    !$OMP PRIVATE( thisSpc, dep,          S                                 )&
    !$OMP PRIVATE( ND,      fracNoHg0Dep, zeroHg0Dep                        )&
    !$OMP COLLAPSE( 2                                                       )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !=====================================================================
       ! Apply dry deposition frequencies
       ! These are the frequencies calculated in drydep_mod.F90
       ! The HEMCO drydep frequencies (from air-sea exchange and
       ! PARANOX) were already added above.
       !
       ! NOTES:
       ! (1) Loops over only the drydep species
       ! (2) If drydep is turned off, nDryDep=0 and the loop won't execute
       ! (3) Tagged species are included in this loop. via species database
       !=====================================================================
       DO ND = 1, State_Chm%nDryDep

          ! Get the species ID from the drydep ID
          N = State_Chm%Map_DryDep(ND)

          IF ( N <= 0 ) CYCLE

          ! Point to the corresponding Species Database entry
          ThisSpc => State_Chm%SpcData(N)%Info

          ! only use the lowest model layer for calculating drydep fluxes
          ! given that spc is in v/v
          dflx(I,J,N) = dflx(I,J,N) + State_Chm%DryDepFreq(I,J,ND) &
                        * State_Chm%Species(N)%Conc(I,J,1)      &
                        /  ( AIRMW / ThisSpc%MW_g )


          IF ( Input_Opt%ITS_A_MERCURY_SIM .and. ThisSpc%Is_Hg0 ) THEN

             ! Hg(0) exchange with the ocean is handled by ocean_mercury_mod
             ! so disable deposition over water here.
             ! Turn off Hg(0) deposition to snow and ice because we haven't yet
             ! included emission from these surfaces and most field studies
             ! suggest Hg(0) emissions exceed deposition during sunlit hours.
             fracNoHg0Dep = MIN( State_Met%FROCEAN(I,J) + &
                                 State_Met%FRSNO(I,J)   + &
                                 State_Met%FRLANDIC(I,J), 1e+0_fp)
             zeroHg0Dep   = ( fracNoHg0Dep > 0e+0_fp )

             IF ( zeroHg0Dep ) THEN
                dflx(I,J,N) = dflx(I,J,N) * MAX( 1.0_fp-fracNoHg0Dep, 0.0_fp )
             ENDIF
          ENDIF

          ! Free species database pointer
          ThisSpc => NULL()
       ENDDO

       !=====================================================================
       ! Convert DFLX from 1/s to kg/m2/s
       !
       ! If applicable, add PARANOX loss to this term. The PARANOX
       ! loss term is already in kg/m2/s. PARANOX loss (deposition) is
       ! calculated for O3 and HNO3 by the PARANOX module, and data is
       ! exchanged via the HEMCO diagnostics.  The data pointers PNOXLOSS_O3
       ! and PNOXLOSS_HNO3 have been linked to these diagnostics at the
       ! beginning of this routine (ckeller, 4/10/15).
       !=====================================================================
       dflx(I,J,:) = dflx(I,J,:) * State_Met%AD(I,J,1)                        &
                                 / State_Grid%Area_M2(I,J)

       IF ( ASSOCIATED( PNOxLoss_O3 ) .AND. id_O3 > 0 ) THEN
          dflx(I,J,id_O3) = dflx(I,J,id_O3) + PNOxLoss_O3(I,J)
       ENDIF

       IF ( ASSOCIATED( PNOXLOSS_HNO3 ) .AND. id_HNO3 > 0 ) THEN
          dflx(I,J,id_HNO3) = dflx(I,J,id_HNO3) + PNOxLOss_HNO3(I,J)
       ENDIF

       !=====================================================================
       ! Surface flux (SFLX) = emissions (EFLX) - dry deposition (DFLX)
       !
       ! SFLX is what we need to pass into routine VDIFF
       !=====================================================================
       State_Chm%SurfaceFlux(I,J,:) = eflx(I,J,:) - dflx(I,J,:) ! kg/m2/s

       !=====================================================================
       ! Defining Satellite Diagnostics
       !=====================================================================

       ! Define emission satellite diagnostics
       IF ( State_Diag%Archive_SatDiagnColEmis ) THEN
          DO S = 1, State_Diag%Map_SatDiagnColEmis%nSlots
             N = State_Diag%Map_SatDiagnColEmis%slot2id(S)
             State_Diag%SatDiagnColEmis(:,:,S) = colEflx(:,:,N)
          ENDDO
       ENDIF

       ! N.B. SatDiagnSurfFlux contains within it the underlying eflx
       ! variable as opposed to colEflx.
       ! Thus, taking SatDiagnSurfFlux - SatDiagnColEmis will not
       ! necessarily equal the dry deposition flux (dflx)
       IF ( State_Diag%Archive_SatDiagnSurfFlux ) THEN
          DO S = 1, State_Diag%Map_SatDiagnSurfFlux%nSlots
             N = State_Diag%Map_SatDiagnSurfFlux%slot2id(S)
             State_Diag%SatDiagnSurfFlux(:,:,S) = State_Chm%SurfaceFlux(:,:,N)
          ENDDO
       ENDIF

       !=====================================================================
       ! Archive Hg deposition for surface reservoirs (cdh, 08/28/09)
       !=====================================================================
       IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

          ! Loop over only the drydep species
          ! If drydep is turned off, nDryDep=0 and the loop won't execute
          DO ND = 1, State_Chm%nDryDep

             ! Get the species ID from the drydep ID
             N = State_Chm%Map_DryDep(ND)

             ! Point to the Species Database entry for tracer N
             ThisSpc => State_Chm%SpcData(N)%Info

             ! Deposition mass, kg
             dep = dflx(I,J,N) * State_Grid%Area_M2(I,J) * GET_TS_CONV()

             IF ( ThisSpc%Is_Hg2 ) THEN

                ! Archive dry-deposited Hg2
                CALL ADD_Hg2_DD      ( I, J, dep                            )
                CALL ADD_Hg2_SNOWPACK( I, J, dep,                    &
                                       State_Met, State_Chm, State_Diag     )

             ELSE IF ( ThisSpc%Is_HgP ) THEN

                ! Archive dry-deposited HgP
                CALL ADD_HgP_DD      ( I, J, dep                            )
                CALL ADD_Hg2_SNOWPACK( I, J, dep,                            &
                                       State_Met, State_Chm, State_Diag     )

             ENDIF

             ! Free pointer
             ThisSpc => NULL()
          ENDDO
       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !### Uncomment for debug output
    !WRITE( 6, '(a)' ) 'eflx and dflx values HEMCO [kg/m2/s]'
    !DO NA = 1, State_Chm%nAdvect
    !   WRITE(6,*) 'eflx TRACER ', NA, ': ', SUM(eflx(:,:,NA))
    !   WRITE(6,*) 'dflx TRACER ', NA, ': ', SUM(dflx(:,:,NA))
    !   WRITE(6,*) 'sflx TRACER ', NA, ': ', SUM(State_Chm%SurfaceFlux(:,:,NA))
    !ENDDO

    !=======================================================================
    ! DIAGNOSTICS: Compute drydep flux loss due to mixing [molec/cm2/s]
    !
    ! NOTE: Dry deposition of "tagged" species (e.g. in tagO3, tagCO, tagHg
    ! specialty simulations) are accounted for in species 1..nDrydep,
    ! so we don't need to do any further special handling.
    !=======================================================================
    IF ( Input_Opt%LGTMM              .or. Input_Opt%LSOILNOX          .or.  &
         State_Diag%Archive_DryDepMix .or. State_Diag%Archive_DryDep ) THEN

       ! Loop over only the drydep species
       ! If drydep is turned off, nDryDep=0 and the loop won't execute
       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( ND, N, ThisSpc, MW_kg, tmpFlx, S                      )
       DO ND = 1, State_Chm%nDryDep

          ! Get the species ID from the drydep ID
          N = State_Chm%Map_DryDep(ND)

          ! Skip if not a valid species
          IF ( N <= 0 ) CYCLE

          ! Point to the Species Database entry for this tracer
          ! NOTE: Assumes a 1:1 tracer index to species index mapping
          ThisSpc => State_Chm%SpcData(N)%Info

          ! Get the molecular weight of the species in kg
          MW_kg = ThisSpc%MW_g * 1.e-3_fp

          !-----------------------------------------------------------------
          ! HISTORY: Update dry deposition flux loss [molec/cm2/s]
          !
          ! DFLX is in kg/m2/s.  We convert to molec/cm2/s by:
          !
          ! (1) multiplying by 1e-4 cm2/m2        => kg/cm2/s
          ! (2) multiplying by ( AVO / MW_KG )    => molec/cm2/s
          !
          ! The term AVO/MW_kg = (molec/mol) / (kg/mol) = molec/kg
          !
          ! NOTE: we don't need to multiply by the ratio of TS_CONV /
          ! TS_CHEM, as the updating frequency for HISTORY is determined
          ! by the "frequency" setting in the "HISTORY.rc"input file.
          ! The old bpch diagnostics archived the drydep due to chemistry
          ! every chemistry timestep = 2X the dynamic timestep.  So in
          ! order to avoid double-counting the drydep flux from mixing,
          ! you had to multiply by TS_CONV / TS_CHEM.
          !
          ! ALSO NOTE: When comparing History output to bpch output,
          ! you must use an updating frequency equal to the dynamic
          ! timestep so that the drydep fluxes due to mixing will
          ! be equivalent w/ the bpch output.  It is also recommended to
          ! turn off chemistry so as to be able to compare the drydep
          ! fluxes due to mixing in bpch vs. History as an "apples-to-
          ! apples" comparison.
          !
          !    -- Bob Yantosca (yantosca@seas.harvard.edu)
          !-----------------------------------------------------------------
          IF ( State_Diag%Archive_DryDepMix   .or.                           &
               State_Diag%Archive_DryDep    ) THEN
             S = State_Diag%Map_DryDepMix%id2slot(ND)
             IF ( S > 0 ) THEN
                State_Diag%DryDepMix(:,:,S) = Dflx(:,:,N)                    &
                                            * 1.0e-4_fp                      &
                                            * ( AVO / MW_kg  )
             ENDIF
          ENDIF

          !-----------------------------------------------------------------
          ! If Soil NOx is turned on, then call SOIL_DRYDEP to
          ! archive dry deposition fluxes for nitrogen species
          ! (SOIL_DRYDEP will exit if it can't find a match.
          !-----------------------------------------------------------------
          IF ( Input_Opt%LSOILNOX ) THEN
             tmpFlx = 0.0_fp
             DO J = 1, State_Grid%NY
             DO I = 1, State_Grid%NX
                tmpFlx = dflx(I,J,N)                                         &
                       / MW_kg                                               &
                       * AVO           * 1.e-4_fp                            &
                       * GET_TS_CONV() / GET_TS_EMIS()

                CALL Soil_DryDep( I, J, 1, N, tmpFlx, State_Chm )
             ENDDO
             ENDDO
          ENDIF

          ! Free species database pointer
          ThisSpc => NULL()
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !=======================================================================
    ! Unit conversion #2: Convert back to the original units
    !=======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid,                &
                            State_Met, OrigUnit,  RC                        )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountred in "Convert_Spc_Units" (from v/v dry)!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Cleanup
    IF ( ASSOCIATED( PNOxLoss_O3 ) )   DEALLOCATE( PNOxLoss_O3 )
    IF ( ASSOCIATED( PNOxLoss_HNO3 ) ) DEALLOCATE( PNOxLoss_HNO3 )

  END SUBROUTINE Compute_Sflx_For_Vdiff
!EOC
END MODULE Hco_Interface_GC_Mod
