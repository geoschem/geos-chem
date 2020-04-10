!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hcoi_standalone_mod.F90
!
! !DESCRIPTION: Module HCOI\_StandAlone\_Mod contains all wrapper routines
! to run HEMCO in standalone mode, i.e. without any external model connected
! to it. All HEMCO input variables (grid, species, times) are read from disk.
! All meteorological variables needed by the (enabled) HEMCO extensions must
! be provided through the HEMCO configuration file (see ExtOpt\_SetPointers).
!\\
! Subroutine HCOI\_StandAlone\_Run will execute the standalone version of
! HEMCO. The following input files are needed for a standalone run:
!
! \begin{itemize}
!
! \item HEMCO\_sa\_Config: the HEMCO configuration file. Must be passed
!  as argument to HCO\_StandAlone\_Run.
!
! \item HEMCO\_sa\_Spec: contains the HEMCO species definitions. The first row
!  must contain the total number of species. For each species, the following
!  parameter need to be specified (separated by at least one space character):
!  species ID, name, molecular weight [g/mol], emitted molecular weight
!  [g/mol], the molecule emission ratio, the liq. over gas Henry constant
!  [M/atm], the temperature dependency of the Henry constant (K0, in [K]), and
!  the pKa (for correction of the Henry constant).
!
! \item HEMCO\_sa\_Grid: contains the definition of the emission grid. Must
!  contain the grid dimensions (NX, NY, NZ) in the first three rows (e.g.
!  NX: 72), followed by the horizontal grid spaces (DX and DY) in rows four
!  and five, respectively. DX and DY can be only one value (applied to all grid
!  boxes), or a vector of length NX or NY, respectively. For now, no vertical
!  regridding is supported, e.g. all emissions input file need to be either
!  2D fields or already remapped onto the correct model levels.
!
! \item HEMCO\_sa\_Time: contains the time definitions. The first two rows must
!  contain the start and end date of the simulation, in format
!  Start/End: YYYY-MM-DD HH:MM:SS (e.g. 2013-07-01 00:00:00).
!  The third row must contain the emission time step (e.g. TS\_EMIS: 3600.0).
!
! \end{itemize}
!
! The file names of the species, grid, and time input files can be provided
! in the settings section of the HEMCO configuration file. For instance, to
! set the species file to 'mySpecFile', add the following line to the con-
! figuration file: 'SpecFile: mySpecFile'. The same applies to grid and time
! definitions (GridFile and TimeFile, respectively). If no file names are
! provided in the configuration file, the default file names (HEMCO\_sa\_Spec,
! HEMCO\_sa\_Grid, HEMCO\_sa\_Time) will be used.
!
! !INTERFACE:
!
MODULE HCOI_StandAlone_Mod
!
! !USES:
!
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
  USE HCO_CharTools_Mod
  USE HCO_Types_Mod
  USE HCOX_State_Mod,      ONLY : Ext_State
  USE HCO_State_Mod,       ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOI_StandAlone_Run
  PUBLIC  :: HCOI_SA_Init
  PUBLIC  :: HCOI_SA_Run
  PUBLIC  :: HCOI_SA_Final
  PUBLIC  :: HCOI_SA_InitCleanup
  PUBLIC  :: Get_nnMatch
  PUBLIC  :: Register_Species
  PUBLIC  :: Define_Diagnostics
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Model_GetSpecies
  PRIVATE :: Set_Grid
  PRIVATE :: Read_Time
  PRIVATE :: ExtState_SetFields
  PRIVATE :: ExtState_UpdateFields
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller   - Initial version.
!  14 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  14 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
!  09 Apr 2015 - C. Keller   - Now accept comments and empty lines in
!                              all input files.
!  15 Feb 2015 - C. Keller   - Update to v2.0
!  18 Jan 2019 - R. Yantosca - Improve error trapping.  Also now made
!                              compatible w/ met field names for FlexGrid
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Default values for HEMCO input files: contain definitions of
  ! species, grid, and time settings, etc.
  CHARACTER(LEN=255)             :: GridFile  = 'HEMCO_sa_Grid'
  CHARACTER(LEN=255)             :: SpecFile  = 'HEMCO_sa_Spec'
  CHARACTER(LEN=255)             :: TimeFile  = 'HEMCO_sa_Time'

  ! HEMCO state
  TYPE(HCO_State),       POINTER :: HcoState  => NULL()

  ! HEMCO extensions state
  TYPE(Ext_State),       POINTER :: ExtState  => NULL()

  ! HEMCO config object
  TYPE(ConfigObj),       POINTER :: HcoConfig => NULL()

  ! Pointers used during initialization (for species matching)
  INTEGER                        :: nHcoSpec
  CHARACTER(LEN= 31),    POINTER :: HcoSpecNames       (:) => NULL()
  INTEGER                        :: nModelSpec
  CHARACTER(LEN= 31),    POINTER :: ModelSpecNames     (:) => NULL()
  INTEGER,               POINTER :: ModelSpecIDs       (:) => NULL()
  REAL(hp),              POINTER :: ModelSpecMW        (:) => NULL()
  REAL(hp),              POINTER :: ModelSpecEmMW      (:) => NULL()
  REAL(hp),              POINTER :: ModelSpecMolecRatio(:) => NULL()
  REAL(hp),              POINTER :: ModelSpecK0        (:) => NULL()
  REAL(hp),              POINTER :: ModelSpecCR        (:) => NULL()
  REAL(hp),              POINTER :: ModelSpecPKA       (:) => NULL()
  INTEGER,               POINTER :: matchidx           (:) => NULL()

  ! Start and end time of simulation
  INTEGER                        :: YRS(2), MTS(2), DYS(2)
  INTEGER                        :: HRS(2), MNS(2), SCS(2)

  ! Grid
  REAL(hp), ALLOCATABLE, TARGET  :: XMID   (:,:,:)
  REAL(hp), ALLOCATABLE, TARGET  :: YMID   (:,:,:)
  REAL(hp), ALLOCATABLE, TARGET  :: XEDGE  (:,:,:)
  REAL(hp), ALLOCATABLE, TARGET  :: YEDGE  (:,:,:)
  REAL(hp), ALLOCATABLE, TARGET  :: YSIN   (:,:,:)
  REAL(hp), ALLOCATABLE, TARGET  :: AREA_M2(:,:,:)

  ! MAXIT is the maximum number of run calls allowed
  INTEGER, PARAMETER             :: MAXIT = 100000
!
! !MODULE INTERFACES:
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_StandAlone_Run
!
! !DESCRIPTION: Subroutine HCOI\_StandAlone\_Run runs the standalone version
! of HEMCO. All input variables are taken from input files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_StandAlone_Run( ConfigFile, IsDryRun, RC )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)  :: ConfigFile   ! HEMCO configuration file
    LOGICAL,          INTENT(IN)  :: IsDryRun     ! Is it a dry-run?
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT) :: RC           ! Success or failure?
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
    LOGICAL            :: am_I_Root

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! HCOI_STANDALONE_RUN begins here!
    !=================================================================

    ! Initialize
    RC        = HCO_SUCCESS
    am_I_Root = .TRUE.         ! Treat this as if we are on the root core!
    ErrMsg    = ''
    ThisLoc   = ' -> at HCOI_StandAlone_Run '                             // &
                '(in module HEMCO/Interfaces/hcoi_standalone_mod.F90)'

    ! Initialize the HEMCO standalone
    CALL HCOI_Sa_Init( am_I_Root, ConfigFile, IsDryRun, RC                  )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "HCO_Sa_Init"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc                   )
       RETURN
    ENDIF

    ! Run the HEMCO standalone
    CALL HCOI_Sa_Run( RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "HCOI_Sa_Run"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc                   )
       RETURN
    ENDIF

    ! Finalize the HEMCO standalone
    CALL HCOI_Sa_Final( )

  END SUBROUTINE HCOI_StandAlone_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_SA_Init
!
! !DESCRIPTION: Subroutine HCOI\_SA\_Init initializes the HEMCO derived
! types and arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_SA_Init( am_I_Root, ConfigFile, IsDryRun, RC )
!
! !USES:
!
    USE HCO_Config_Mod,    ONLY : Config_ReadFile
    USE HCO_State_Mod,     ONLY : HcoState_Init
    USE HCO_Driver_Mod,    ONLY : HCO_Init
    USE HCOX_Driver_Mod,   ONLY : HCOX_Init
    USE HCO_EXTLIST_Mod,   ONLY : GetExtOpt, CoreNr
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root   ! Are we on the root core?
    CHARACTER(LEN=*), INTENT(IN)    :: ConfigFile  ! Configuration file
    LOGICAL,          INTENT(IN)    :: IsDryRun    ! Is it a dry-run?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT) :: RC          ! Failure or success
!
! !REVISION HISTORY:
!  12 Sep 2013 - C. Keller   - Initial version
!  18 Jan 2019 - R. Yantosca - Improve error trapping
!  29 Jan 2019 - R. Yantosca - Now flush errmsgs to logfile before exiting
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: nnMatch, LUN
    LOGICAL            :: Dum,     Found

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg,  ThisLoc

    !=======================================================================
    ! HCOI_SA_INIT begins here!
    !=======================================================================

    ! Initialize
    RC      = HCO_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
     'HCOI_SA_Init (in module HEMCO/Interfaces/hcoi_standalone_mod.F90)'

    !=======================================================================
    ! Read HEMCO configuration file and save into buffer. This also
    ! sets the HEMCO error properties (verbose mode? log file name,
    ! etc.) based upon the specifications in the configuration file.
    !=======================================================================
    CALL Config_ReadFile( am_I_Root, HcoConfig, ConfigFile,       &
                          0,         RC,        IsDryRun=IsDryRun )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Config_Readfile!"'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Open logfile
    !======================================================================
    IF ( am_I_Root ) THEN
       CALL HCO_LogFile_Open( HcoConfig%Err, RC=RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Error encountered in routine "HCO_Logfile_Open_Readfile!"'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Initialize HEMCO state object and populate it
    !=======================================================================

    !-----------------------------------------------------------------------
    ! Extract species to use in HEMCO
    CALL Get_nnMatch( HcoConfig, nnMatch, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Get_nnMatch"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Initialize HCO state. Use only species that are used
    ! in HEMCO_sa_Spec.rc and are also found in the HEMCO config. file.
    CALL HcoState_Init( HcoState, HcoConfig, nnMatch, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "HcoState_Init"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Set grid
    CALL Set_Grid ( HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Set_Grid"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Register species
    CALL Register_Species( HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Register_Species"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Read time information, incl. timesteps and simulation time(s)
    CALL Read_Time( HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Read_Time"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !=======================================================================
    ! Set misc. parameter
    !=======================================================================

    ! Set ESMF flag
    HcoState%Options%isESMF = .FALSE.

    ! Let HEMCO schedule the diagnostics output
    HcoState%Options%HcoWritesDiagn = .TRUE.

    ! If not explicitly set, make sure that option Field2Diagn is true
    CALL GetExtOpt ( HcoState%Config, CoreNr, &
                    'ConfigField to diagnostics', &
                     OptValBool=Dum, Found=Found, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "GetExtOpt"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF
    IF ( .NOT. Found ) HcoState%Options%Field2Diagn = .TRUE.

    !=======================================================================
    ! Are we running the HEMCO standalone in a dry-run mode?
    ! This is dictated by the HEMCO environment. If HEMCO is in a
    ! dry-run mode, no compute is performed and files are only "checked".
    ! Simulations will NOT stop on missing files. This is intended to be a
    ! quick sanity check to make sure that GEOS-Chem IO are all correctly
    ! set up, which is why most of the runs fail to complete successfully.
    ! (hplin, 11/2/19)
    !
    ! Dry-run simulations now send output to a log file that is separate
    ! from the HEMCO log files. (bmy, 11/11/19)
    !
    ! NOTE: The dry-run option is not invoked when we use HEMCO
    ! in external ESMs. (bmy, 11/13/19)
    !=======================================================================
    CALL Init_Dry_Run( IsDryRun, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Init_Dry_Run"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !=======================================================================
    ! Initialize HEMCO internal lists and variables. All data
    ! information is written into internal lists (ReadList) and
    ! the HEMCO configuration file is removed from buffer in this
    ! step. Also initializes the HEMCO clock
    !=======================================================================
    CALL HCO_Init( HcoState, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "HCO_Init"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !=======================================================================
    ! Initialize extensions.
    ! This initializes all (enabled) extensions and selects all met.
    ! fields needed by them.
    !=======================================================================
    CALL HCOX_Init( HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "HCOX_Init"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       CALL Flush( HcoState%Config%Err%Lun )
       RETURN
    ENDIF

    !=======================================================================
    ! Define diagnostics
    !=======================================================================
    IF ( HcoState%Options%IsDryRun ) THEN

       !--------------------------------------------------------------------
       ! For dry-run simulations, print the status of the HEMCO
       ! diagnostic configurations file (but do not read from it)
       !--------------------------------------------------------------------
       CALL DiagnFileOpen( HcoConfig, LUN, RC, IsDryRun=.TRUE. )

    ELSE

       !--------------------------------------------------------------------
       ! For regular simulations, read diagnostics configuration file
       ! and define diagnostic variables for output
       !--------------------------------------------------------------------
       CALL Define_Diagnostics( HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Error encountered in routine "Define_Diagnostics"!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL Flush( HcoState%Config%Err%Lun )
          RETURN
       ENDIF
    ENDIF

    !=======================================================================
    ! Leave
    !=======================================================================
    CALL HCOI_SA_InitCleanup( RC )

    ! Leave w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOI_SA_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_SA_Run
!
! !DESCRIPTION: Subroutine HCOI\_SA\_RUN runs HCO from GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_SA_Run( RC )
!
! !USES:
!
    USE HCO_FluxArr_Mod, ONLY : HCO_FluxarrReset
    USE HCO_Clock_Mod,   ONLY : HcoClock_Set
    USE HCO_Clock_Mod,   ONLY : HcoClock_Get
    USE HCO_Clock_Mod,   ONLY : HcoClock_Increase
    USE HCO_Driver_Mod,  ONLY : HCO_RUN
    USE HCOX_Driver_Mod, ONLY : HCOX_RUN
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT) :: RC         ! Failure or success
!
! !REVISION HISTORY:
!  12 Sep 2013 - C. Keller   - Initial version
!  18 Jan 2019 - R. Yantosca - Improve error trapping
!  29 Jan 2019 - R. Yantosca - Bug fix: Call HCO_RUN twice, once with phase=1
!                              and again with phase=2.  This is necessary
!                              for emissions to be computed.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: notDryRun
    INTEGER            :: CNT
    INTEGER            :: YR, MT, DY, HR, MN, SC

    ! Strings
    CHARACTER(LEN=255) :: Msg, ErrMsg, ThisLoc

    !=======================================================================
    ! HCOI_SA_RUN begins here!
    !=======================================================================

    ! Initialize
    RC        = HCO_SUCCESS
    notDryRun = ( .not. HcoState%Options%IsDryRun )
    ErrMsg    = ''
    ThisLoc   = &
     ' -> at HCOI_SA_Run (in module HEMCO/Standalone/hcoi_standalone_mod.F90)'

    ! Time step counter
    CNT = 0

    ! Do until end of simulation
    DO

       ! Increase counter by one
       CNT = CNT + 1

       ! Set iteration limit to avoid infinite runs
       IF ( CNT > MAXIT ) THEN
          WRITE(*,*) 'Counter limit reached - Increase MAXIT if you don`t like that!'
          EXIT
       ENDIF

       !====================================================================
       ! Set HcoClock. On first call, use specified start date.
       ! Increase clock by one emission time step otherwise.
       !====================================================================
       IF ( CNT == 1 ) THEN
          CALL HcoClock_Set ( HcoState,  YRS(1), MTS(1), &
                              DYS(1),    HRS(1), MNS(1), SCS(1), &
                              IsEmisTime=.TRUE., RC=RC)
          IF ( RC /= HCO_SUCCESS) RETURN
       ELSE
          CALL HcoClock_Increase ( HcoState, HcoState%TS_EMIS, .TRUE., RC=RC )
          IF ( RC /= HCO_SUCCESS) RETURN
       ENDIF

       ! Get current time
       CALL HcoClock_Get ( HcoState%Clock, cYYYY=YR, &
                           cMM=MT, cDD=DY, cH=HR, cM=MN, cS=SC, RC=RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Error encountered in routine "HcoClock_Get"!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Leave loop if this is the end of the simulation
       IF ( IsEndOfSimulation(YR,MT,DY,HR,MN,SC) ) EXIT

       ! Write to logfile and standard output (skip for dry-run)
       IF ( notDryRun ) THEN
          WRITE( Msg, 100 ) YR, MT, DY, HR, MN, SC
100       FORMAT( 'Calculate emissions at ', i4,  '-', i2.2 ,'-', i2.2,' ',  &
                                             i2.2,':', i2.2, ':', i2.2      )
          CALL HCO_MSG(HcoState%Config%Err,Msg)
          WRITE(*,*) TRIM( MSG )
       ENDIF

       ! ================================================================
       ! Reset all emission and deposition values
       ! ================================================================
       IF ( notDryRun ) THEN
          CALL HCO_FluxArrReset( HcoState, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Error encountered in routine "HCO_FluxArrReset"!'
             CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       ! ================================================================
       ! Set HCO options and define all arrays needed by core module
       ! and the extensions
       ! ================================================================

       ! Range of tracers and emission categories.
       ! Set Extension number ExtNr to 0, indicating that the core
       ! module shall be executed.
       HcoState%Options%SpcMin = 1
       HcoState%Options%SpcMax = nModelSpec
       HcoState%Options%CatMin = 1
       HcoState%Options%CatMax = -1
       HcoState%Options%ExtNr  = 0

       ! Use temporary array?
       HcoState%Options%FillBuffer = .FALSE.

       ! ================================================================
       ! Run HCO core module
       ! Emissions will be written into the corresponding flux arrays
       ! in HcoState.
       !
       ! NOTE: Call HCO_Run explicitly twice, once for phase 1 and
       ! once for phase 2.  This will ensure emissions get computed.
       ! (bmy, 1/29/18)
       ! ================================================================

       ! Phase 1: Update reading data fields etc.
       CALL HCO_Run( HcoState, 1, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Error encountered in routine "Hco_Run", phase 1!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Phase 2: Compute emissions (skip for dry-run)
       IF ( notDryRun ) THEN
          CALL HCO_Run( HcoState, 2, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Error encountered in routine "Hco_Run", phase 2!'
             CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       ! ================================================================
       ! Run HCO extensions
       ! ================================================================
       IF ( notDryRun ) THEN

          ! Set ExtState fields (skip for dry-run)
          CALL ExtState_SetFields ( HcoState, ExtState, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Error encountered in routine "ExtState_SetFields"!'
             CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Update ExtState fields (skip for dry-run)
          CALL ExtState_UpdateFields( HcoState, ExtState, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Error encountered in routine "ExtState_Update_Fields"!'
             CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       ! Execute all enabled emission extensions. Emissions will be
       ! added to corresponding flux arrays in HcoState.
       CALL HCOX_Run ( HcoState, ExtState, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Error encountered in routine "HCOX_Run"!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       !=================================================================
       ! Update all autofill diagnostics (skip for dry-run)
       !=================================================================
       IF ( notDryRun ) THEN
          CALL HcoDiagn_AutoUpdate ( HcoState, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Error encountered in routine "HCOX_AutoUpdate"!'
             CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF
    ENDDO

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOI_SA_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_SA_Final
!
! !DESCRIPTION: Subroutine HCOI\_SA\_Final cleans up HEMCO.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_SA_Final( )
!
! !USES:
!
    USE HCO_Driver_Mod,  ONLY : HCO_Final
    USE HCOX_Driver_Mod, ONLY : HCOX_Final
    USE HCO_State_Mod,   ONLY : HcoState_Final
!
! !REVISION HISTORY:
!  12 Sep 2013 - C. Keller   - Initial version
!  18 Jan 2019 - R. Yantosca - Improve error trapping
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I,      RC
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    CHARACTER(LEN= 31) :: RST='HEMCO_restart'

    !=================================================================
    ! HCOI_SA_FINAL begins here!
    !=================================================================

    ! Initialize
    RC      = HCO_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
    'HCOI_SA_FINAL (in module HEMCO/Interfaces/hcoi_standalone_mod.F90)'

    ! Cleanup the dry-run
    CALL Cleanup_Dry_Run( RC )

    ! Cleanup HCO core
    CALL HCO_FINAL( HcoState, .FALSE., RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "HCO_Final"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Cleanup extensions and ExtState object
    ! This will also nullify all pointer to the met fields.
    CALL HCOX_FINAL( HcoState, ExtState, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "HCOX_Final"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Cleanup diagnostics (skip if dry-run)
    IF ( .not. HcoState%Options%IsDryRun ) THEN
       CALL DiagnBundle_Cleanup( HcoState%Diagn )
    ENDIF

    ! Deallocate module arrays/pointers
    IF ( ALLOCATED( XMID    ) ) DEALLOCATE ( XMID    )
    IF ( ALLOCATED( YMID    ) ) DEALLOCATE ( YMID    )
    IF ( ALLOCATED( XEDGE   ) ) DEALLOCATE ( XEDGE   )
    IF ( ALLOCATED( YEDGE   ) ) DEALLOCATE ( YEDGE   )
    IF ( ALLOCATED( YSIN    ) ) DEALLOCATE ( YSIN    )
    IF ( ALLOCATED( AREA_M2 ) ) DEALLOCATE ( AREA_M2 )

    ! Cleanup HcoState object
    CALL HcoState_Final( HcoState )

  END SUBROUTINE HCOI_SA_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Model_GetSpecies
!
! !DESCRIPTION: SUBROUTINE Model\_GetSpecies returns 'model' species
! information from the HEMCO standalone input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Model_GetSpecies( HcoConfig,                          &
                               nModelSpec,     ModelSpecNames,     &
                               ModelSpecIDs,   ModelSpecMW,        &
                               ModelSpecEmMW,  ModelSpecMolecRatio,&
                               ModelSpecK0,    ModelSpecCR,        &
                               ModelSpecPKA,   RC                   )
!
! !USES:
!
    USE inquireMod,      ONLY : findfreeLUN
    USE HCO_EXTLIST_Mod, ONLY : GetExtOpt, CoreNr
!
! !OUTPUT PARAMETERS:
!
    TYPE(ConfigObj),    POINTER     :: HcoConfig
    INTEGER,            INTENT(OUT) :: nModelSpec
    CHARACTER(LEN= 31), POINTER     :: ModelSpecNames     (:)
    INTEGER,            POINTER     :: ModelSpecIDs       (:)
    REAL(hp),           POINTER     :: ModelSpecMW        (:)
    REAL(hp),           POINTER     :: ModelSpecEmMW      (:)
    REAL(hp),           POINTER     :: ModelSpecMolecRatio(:)
    REAL(hp),           POINTER     :: ModelSpecK0        (:)
    REAL(hp),           POINTER     :: ModelSpecCR        (:)
    REAL(hp),           POINTER     :: ModelSpecPKA       (:)
    INTEGER,            INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER             :: I, N, LNG, LOW, UPP
    INTEGER             :: IU_FILE, IOS
    LOGICAL             :: FOUND,   EOF
    CHARACTER(LEN=255)  :: MSG, LOC
    CHARACTER(LEN=255)  :: MySpecFile
    CHARACTER(LEN=2047) :: DUM

    !=================================================================
    ! Model_GetSpecies begins here
    !=================================================================

    ! For error handling
    LOC = 'Model_GetSpecies (hcoi_standalone_mod.F90)'

    ! Try to get SpecFile from configuration file (in settings)
    CALL GetExtOpt ( HcoConfig, CoreNr, 'SpecFile', &
                     OptValChar=MySpecFile,   Found=FOUND, RC=RC )
    !IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( FOUND ) THEN
       SpecFile = MySpecFile
    ELSE
       MSG = 'Please provide filename with species definitions ' // &
             'in the configuration file settings, e.g. ' // &
             'SpecFile: MySpecies.rc'
       CALL HCO_Error ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Open spec file
    OPEN( IU_FILE, FILE=TRIM(SpecFile), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       MSG = 'Error 1 reading ' // TRIM(SpecFile)
       CALL HCO_Error( HcoConfig%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Get number of species
    nModelSpec = 0
    DO
       CALL GetNextLine( IU_FILE, DUM, EOF, RC )
       IF ( EOF               ) EXIT
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Error encountered in reading SpecFile!.  Please ' // &
                'doublecheck that all species information has '    // &
                'been correctly entered.'
          CALL HCO_Error ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
       ENDIF
       nModelSpec = nModelSpec + 1
    ENDDO

    ! Make sure we have one species
    IF ( nModelSpec == 0 ) THEN
       MSG = 'Species file ' // TRIM(SpecFile)      // &
             ' does not seem to have any content. ' // &
             'You must define at least one species.'
       CALL HCO_Error( HcoConfig%Err, MSG, RC, THISLOC=LOC )
    ENDIF

    ! Go back to line one
    REWIND( IU_FILE )

    ! Get next valid line
!    CALL GetNextLine( IU_FILE, DUM, EOF, RC )
!    IF ( RC /= HCO_SUCCESS .OR. EOF ) THEN
!       MSG = 'Error 2 reading ' // TRIM(SpecFile)
!       CALL HCO_Error( MSG, RC, THISLOC=LOC )
!       RETURN
!    ENDIF
!
!    LNG = LEN(TRIM(DUM))
!    LOW = NextCharPos ( TRIM(DUM), HCO_COL(), 1 )
!    IF ( LOW < 0 .OR. LOW == LNG ) THEN
!       MSG = 'Cannot extract index after colon: ' // TRIM(DUM)
!       CALL HCO_Error( MSG, RC, THISLOC=LOC )
!       RETURN
!    ENDIF
!    LOW = LOW + 1
!    READ ( DUM(LOW:LNG), * ) nModelSpec

    ! Allocate species arrays
    ALLOCATE(ModelSpecNames     (nModelSpec))
    ALLOCATE(ModelSpecIDs       (nModelSpec))
    ALLOCATE(ModelSpecMW        (nModelSpec))
    ALLOCATE(ModelSpecEmMW      (nModelSpec))
    ALLOCATE(ModelSpecMolecRatio(nModelSpec))
    ALLOCATE(ModelSpecK0        (nModelSpec))
    ALLOCATE(ModelSpecCR        (nModelSpec))
    ALLOCATE(ModelSpecPKA       (nModelSpec))

    ! Assign variables to each species
    DO N = 1, nModelSpec

       CALL GetNextLine( IU_FILE, DUM, EOF, RC )
       IF ( RC /= HCO_SUCCESS .OR. EOF ) THEN
          WRITE(MSG,100) N, TRIM(SpecFile)
          CALL HCO_Error( HcoConfig%Err, MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Start reading line from beginning
       LNG = LEN(TRIM(DUM))
       LOW = 0

       ! Read species ID, name, molecular weight, emitted molecular weight,
       ! molecular coefficient, and Henry coefficients K0, CR, pKa (in this
       ! order).
       DO I = 1, 8

          ! Get lower and upper index of species ID (first entry in row).
          ! Skip all leading spaces.
          UPP = LOW

          DO WHILE( UPP == LOW .AND. LOW /= LNG )
             LOW = LOW + 1
             IF ( LOW > LNG ) THEN
                WRITE(MSG,101) I, TRIM(DUM)
                CALL HCO_Error( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             UPP = NextCharPos( TRIM(DUM), HCO_SPC, LOW )
             IF ( UPP < 0 ) UPP = LNG
          ENDDO

          IF ( I < 8 ) THEN
             UPP = UPP - 1 ! Don't read space
          ENDIF

          ! Error check
          IF ( UPP > LNG ) THEN
             WRITE(MSG,*) 'Error reading species property ', I, &
                          ' on line ', TRIM(DUM), '. Each ', &
                          'species definition line is expected ', &
                          'to have 8 entries (ID, Name, MW, MWemis, ', &
                          'MOLECRATIO, K0, CR, PKA, e.g.: ', &
                          '1 CO   28.0 28.0 1.0 0.0 0.0 0.0'
             CALL HCO_Error ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
             RETURN
          ENDIF

          ! Read into vector
          SELECT CASE ( I )
             CASE ( 1 )
                READ( DUM(LOW:UPP), * ) ModelSpecIDs(N)
             CASE ( 2 )
                READ( DUM(LOW:UPP), * ) ModelSpecNames(N)
             CASE ( 3 )
                READ( DUM(LOW:UPP), * ) ModelSpecMW(N)
             CASE ( 4 )
                READ( DUM(LOW:UPP), * ) ModelSpecEmMW(N)
             CASE ( 5 )
                READ( DUM(LOW:UPP), * ) ModelSpecMolecRatio(N)
             CASE ( 6 )
                READ( DUM(LOW:UPP), * ) ModelSpecK0(N)
             CASE ( 7 )
                READ( DUM(LOW:UPP), * ) ModelSpecCR(N)
             CASE ( 8 )
                READ( DUM(LOW:UPP), * ) ModelSpecPKA(N)
          END SELECT

          ! Continue from upper position (+1 to skip space). The
          ! while loop at the beginning of the do-loop will advance
          ! low by another one position, so the next character
          ! search will start at position UPP + 2, which is exactly
          ! what we want (UPP is the position BEFORE the space!).
          LOW = UPP + 1

       ENDDO !I
    ENDDO !N

    ! Close file
    CLOSE( IU_FILE )

    ! Make sure that the species indexing starts at 1
    IF ( MINVAL( ModelSpecIDs ) /= 1 ) THEN
       MSG = 'Error encountered in reading SpecFile!.  The species '      // &
             'ID numbers do not start at 1!  Please check SpecFile '      // &
             'for typos.'
       CALL HCO_Error ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Make sure that the ID of the last species is the same as nModelSpec
    IF ( MAXVAL( ModelSpecIDs ) /= nModelSpec ) THEN
       MSG = 'Error encountered in reading SpecFile!.  The ID number '    // &
             'of the last species does not match the number of species '  // &
             'that were read from SpecFile!  Please check SpecFile for '  //&
             'typos.'
       CALL HCO_Error ( HcoConfig%Err, MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

100 FORMAT( 'Error reading species ', i3, ' in ', a )
101 FORMAT( 'Cannot extract element ', i1, ' in ', a )

  END SUBROUTINE Model_GetSpecies
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Grid
!
! !DESCRIPTION: SUBROUTINE SET\_GRID reads the grid information from the
!  HEMCO standalone grid file and sets all HEMCO grid arrays accordingly.
!  The grid file is expected to contain information on the grid edge lon/lat
!  range, as well as the number of grid cells in longitude and latitude
!  direction.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_Grid( HcoState, RC )
!
! !USES:
!
    USE inquireMod,       ONLY  : findFreeLUN
    USE HCO_ExtList_Mod,  ONLY  : HCO_GetOpt, GetExtOpt, CoreNr
    USE HCO_VertGrid_Mod, ONLY  : HCO_VertGrid_Define
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_STATE), POINTER       :: HcoState
    INTEGER,         INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller - Initial Version
!  11 May 2015 - C. Keller - Now provide lon/lat edges instead of assuming
!                            global grid.
!  10 Sep 2015 - C. Keller - Allow to provide mid-points instead of edges.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    INTEGER               :: NX, NY, NZ
    INTEGER               :: I, J, N, LNG, LOW, UPP
    INTEGER               :: IU_FILE, IOS, STRT
    REAL(hp)              :: XMIN, XMAX
    REAL(hp)              :: YMIN, YMAX
    REAL(hp)              :: DVAL
    REAL(hp)              :: DLON, DLAT
    REAL(hp)              :: PI_180, YDGR, YSN, SIN_DELTA, AM2
    LOGICAL               :: FOUND,   EOF

    ! Arrays
    INTEGER               :: SZ(3)
    REAL(hp)              :: RG(4)
    REAL(hp), ALLOCATABLE :: Ap(:), Bp(:)

    ! Strings
    CHARACTER(LEN=255)    :: LOC
    CHARACTER(LEN=  1)    :: COL
    CHARACTER(LEN=255)    :: MyGridFile, ThisLoc
    CHARACTER(LEN=4095)   :: DUM,        ErrMsg,  Msg

    !=================================================================
    ! SET_GRID begins here
    !=================================================================

    ! Initialize
    RC      = HCO_SUCCESS
    Msg     = ''
    ErrMsg  = ''
    ThisLoc = &
     'SET_GRID (in module HEMCO/Interfaces/hcoi_standalone_mod.F90)'

    ! Set PI_180
    PI_180 = HcoState%Phys%PI / 180.0_hp

    ! Try to get GridFile from configuration file (in settings)
    CALL GetExtOpt ( HcoState%Config, CoreNr, 'GridFile', &
                     OptValChar=MyGridFile,   Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "GetExtOpt"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( FOUND ) GridFile = MyGridFile

    ! Write colon character to local variable
    COL = HCO_GetOpt( HcoState%Config%ExtList, 'Colon' )

    ! ------------------------------------------------------------------
    ! Open grid file
    ! ------------------------------------------------------------------

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Open grid file
    OPEN( IU_FILE, FILE=TRIM(GridFile), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       ErrMsg = 'Error 1 reading ' // TRIM(GridFile)
       CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
       RETURN
    ENDIF

    ! ------------------------------------------------------------------
    ! Extract grid range
    ! The lon/lat grid ranges are expected to be provided first, with
    ! each range provided in a separate line:
    ! XMIN: -180.0
    ! XMAX:  180.0
    ! YMIN:  -90.0
    ! YMAX:   90.0
    ! ------------------------------------------------------------------
    DO N = 1,4

       ! Get next valid line
       CALL GetNextLine( IU_FILE, DUM, EOF, RC )
       IF ( RC /= HCO_SUCCESS .OR. EOF ) THEN
          ErrMsg= 'Error 2 reading ' // TRIM(GridFile)
          CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
          RETURN
       ENDIF

       ! Read integer after colon (this is the dimension size)
       LNG = LEN(TRIM(DUM))
       LOW = NextCharPos ( TRIM(DUM), COL, 1 )
       IF ( LOW < 0 .OR. LOW == LNG ) THEN
          ErrMsg = 'Cannot extract size information from ' // TRIM(DUM)
          CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
          RETURN
       ENDIF
       LOW = LOW + 1
       READ( DUM(LOW:LNG), * ) RG(N)

    ENDDO

    ! Pass to scalars
    XMIN = RG(1)
    XMAX = RG(2)
    YMIN = RG(3)
    YMAX = RG(4)

    ! Make sure values are in valid range
    IF ( XMIN >= XMAX ) THEN
       WRITE(ErrMsg,*) 'Lower lon must be smaller than upper lon: ', XMIN, XMAX
       CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
       RETURN
    ENDIF
    IF ( YMIN >= YMAX ) THEN
       WRITE(ErrMsg,*) 'Lower lat must be smaller than upper lat: ', YMIN, YMAX
       CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
       RETURN
    ENDIF

    ! Restrict latitude values to -90.0 and 90.0.
    IF ( YMIN < -90.0_hp ) THEN
       WRITE(ErrMsg,*) 'Lower latitude must be between -90 and 90 degN: ', YMIN
       CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
       RETURN
    ENDIF
    IF ( YMAX > 90.0_hp ) THEN
       WRITE(ErrMsg,*) 'Upper latitude must be between -90 and 90 degN: ', YMAX
       CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
       RETURN
    ENDIF

    ! ------------------------------------------------------------------
    ! Extract grid size (x,y,z)
    ! The grid sizes are expected to be provided in three separte lines:
    ! NX: 360
    ! NY: 180
    ! NZ: 1
    ! ------------------------------------------------------------------
    DO N = 1,3

       ! Get next valid line
       CALL GetNextLine( IU_FILE, DUM, EOF, RC )
       IF ( RC /= HCO_SUCCESS .OR. EOF ) THEN
          ErrMsg = 'Error 3 reading ' // TRIM(GridFile)
          CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
          RETURN
       ENDIF

       ! Read integer after colon (this is the dimension size)
       LNG = LEN(TRIM(DUM))
       LOW = NextCharPos ( TRIM(DUM), COL, 1 )
       IF ( LOW < 0 .OR. LOW == LNG ) THEN
          ErrMsg = 'Cannot extract size information from ' // TRIM(DUM)
          CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
          RETURN
       ENDIF
       LOW = LOW + 1
       READ( DUM(LOW:LNG), * ) SZ(N)

    ENDDO !N

    ! Grid dimensions
    NX = SZ(1)
    NY = SZ(2)
    NZ = SZ(3)

    ! ------------------------------------------------------------------
    ! Now that sizes are known, allocate all arrays
    ! ------------------------------------------------------------------
    ALLOCATE ( XMID     (NX,  NY,  1   ) )
    ALLOCATE ( YMID     (NX,  NY,  1   ) )
    ALLOCATE ( XEDGE    (NX+1,NY,  1   ) )
    ALLOCATE ( YEDGE    (NX,  NY+1,1   ) )
    ALLOCATE ( YSIN     (NX,  NY+1,1   ) )
    ALLOCATE ( AREA_M2  (NX,  NY,  1   ) )
    ALLOCATE ( AP       (          NZ+1) )
    ALLOCATE ( BP       (          NZ+1) )
    YSIN      = HCO_MISSVAL
    AREA_M2   = HCO_MISSVAL
    XMID      = HCO_MISSVAL
    YMID      = HCO_MISSVAL
    XEDGE     = HCO_MISSVAL
    YEDGE     = HCO_MISSVAL
    AP        = HCO_MISSVAL
    BP        = HCO_MISSVAL

    ! ------------------------------------------------------------------
    ! Check if grid box edges and/or midpoints are explicitly given.
    ! Those need be provided on one line, e.g.:
    ! YEDGE: -90.0 -89.0 -86.0 ... 86.0 89.0 90.0
    ! ------------------------------------------------------------------
    DO N = 1, 6 ! check for XEDGE, YEDGE, XMID, YMID

       ! Try to read line
       CALL GetNextLine( IU_FILE, DUM, EOF, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Error reading grid edges and/or midpoints in ' // TRIM(GridFile)
          CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
          RETURN
       ENDIF

       ! Exit loop here if end of file
       IF ( EOF ) EXIT

       ! Read XEDGES or YEDGES
       LNG = -1
       IF ( DUM(1:5) == 'XEDGE' .OR. DUM(1:5) == 'YEDGE' ) THEN
          LNG  = LEN(TRIM(DUM))
          STRT = 7 ! Start at string position 7 (e.g. 'XEDGE: XXX')
       ELSEIF ( DUM(1:4) == 'XMID' .OR. DUM(1:4) == 'YMID' ) THEN
          LNG = LEN(TRIM(DUM))
          STRT = 6 ! Start at string position 6 (e.g. 'XMID: XXX')
       ELSEIF ( DUM(1:2) == 'AP' .OR. DUM(1:2) == 'BP' ) THEN
          LNG = LEN(TRIM(DUM))
          STRT = 4 ! Start at string position 4 (e.g. 'AP: XXX')
       ENDIF

       IF ( LNG > 0 ) THEN

          LOW = -1
          UPP = -1
          I   = 0

          ! Walk through entire string
          DO J = STRT, LNG

             ! Need to evaluate if this is the last string character and/or
             ! whitespace character
             IF ( TRIM(DUM(J:J)) == HCO_SPC ) THEN

                ! If the lower substring bound is not set yet, assume that this
                ! is a lower substring bound, and continue search for upper bound
                IF ( LOW == -1 ) LOW = J

                ! Make sure the substring bounds are valid values
                IF ( (J-1) >= (LOW+1) ) THEN
                   UPP = J
                ELSE
                   LOW = J
                ENDIF

             ENDIF

             ! If this is the last character, set upper substring bound to J
             IF ( J == LNG ) UPP = J

             ! Read substring if both bounds are defined
             IF ( UPP > LOW ) THEN

                ! Read value
                READ( DUM(LOW:UPP), * ) DVAL

                ! Index to fill
                I = I + 1

                ! Pass to XEDGE
                IF ( TRIM(DUM(1:5)) == 'XEDGE' ) THEN
                   IF ( I > NX+1 ) THEN
                      WRITE(ErrMsg,*) 'More than ', NX+1, ' longitude edges found in ', TRIM(DUM)
                      CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
                      RETURN
                   ENDIF
                   XEDGE(I,:,1) = DVAL

                ! Pass to YEDGE
                ELSEIF ( TRIM(DUM(1:5)) == 'YEDGE' ) THEN
                   IF ( I > NY+1 ) THEN
                      WRITE(ErrMsg,*) 'More than ', NY+1, ' latitude edges found in ', TRIM(DUM)
                      CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
                      RETURN
                   ENDIF
                   YEDGE(:,I,1) = DVAL

                ! Pass to XMID
                ELSEIF ( TRIM(DUM(1:4)) == 'XMID' ) THEN
                   IF ( I > NX ) THEN
                      WRITE(ErrMsg,*) 'More than ', NX, ' latitude mid-points found in ', TRIM(DUM)
                      CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
                      RETURN
                   ENDIF
                   XMID(I,:,1) = DVAL

                ! Pass to YMID
                ELSEIF ( TRIM(DUM(1:4)) == 'YMID' ) THEN
                   IF ( I > NY ) THEN
                      WRITE(ErrMsg,*) 'More than ', NY, ' latitude mid-points found in ', TRIM(DUM)
                      CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
                      RETURN
                   ENDIF
                   YMID(:,I,1) = DVAL

                ! Pass to Ap
                ELSEIF ( TRIM(DUM(1:2)) == 'AP' ) THEN
                   IF ( I > (NZ+1) ) THEN
                      WRITE(ErrMsg,*) 'More than ', NZ+1, ' Ap values found in ', TRIM(DUM)
                      CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
                      RETURN
                   ENDIF
                   AP(I) = DVAL

                ! Pass to Bp
                ELSEIF ( TRIM(DUM(1:2)) == 'BP' ) THEN
                   IF ( I > (NZ+1) ) THEN
                      WRITE(ErrMsg,*) 'More than ', NZ+1, ' Bp values found in ', TRIM(DUM)
                      CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
                      RETURN
                   ENDIF
                   BP(I) = DVAL
                ENDIF

                ! Update bounds
                LOW = UPP
             ENDIF
          ENDDO

          ! Error check: all values must have been filled
          IF ( TRIM(DUM(1:5)) == 'XEDGE' .AND. I /= NX+1 ) THEN
             WRITE(ErrMsg,*) 'Error reading XEDGES: exactly ', NX+1, ' values must be given: ', TRIM(DUM)
             CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
             RETURN
          ENDIF
          IF ( TRIM(DUM(1:5)) == 'YEDGE' .AND. I /= NY+1 ) THEN
             WRITE(ErrMsg,*) 'Error reading YEDGES: exactly ', NY+1, ' values must be given: ', TRIM(DUM)
             CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
             RETURN
          ENDIF
          IF ( TRIM(DUM(1:4)) == 'XMID' .AND. I /= NX ) THEN
             WRITE(ErrMsg,*) 'Error reading XMID: exactly ', NX, ' values must be given: ', TRIM(DUM)
             CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
             RETURN
          ENDIF
          IF ( TRIM(DUM(1:4)) == 'YMID' .AND. I /= NY ) THEN
             WRITE(ErrMsg,*) 'Error reading YMID: exactly ', NY, ' values must be given: ', TRIM(DUM)
             CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
             RETURN
          ENDIF
          IF ( TRIM(DUM(1:2)) == 'AP' .AND. I /= NZ+1 ) THEN
             WRITE(ErrMsg,*) 'Error reading AP: exactly ', NZ+1, ' values must be given: ', TRIM(DUM)
             CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
             RETURN
          ENDIF
          IF ( TRIM(DUM(1:2)) == 'BP' .AND. I /= NZ+1 ) THEN
             WRITE(ErrMsg,*) 'Error reading BP: exactly ', NZ+1, ' values must be given: ', TRIM(DUM)
             CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
             RETURN
          ENDIF

       ENDIF
    ENDDO

    ! Error check: if AP is given, Bp must be given as well
    IF ( ALL(AP==HCO_MISSVAL) .AND. .NOT. ALL(BP==HCO_MISSVAL) ) THEN
       WRITE(ErrMsg,*) 'At least a few AP values are missing, please provide exactly ', &
                    NZ+1, 'AP and BP values.'
       CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
       RETURN
    ELSEIF ( .NOT. ALL(AP==HCO_MISSVAL) .AND. ALL(BP==HCO_MISSVAL) ) THEN
       WRITE(ErrMsg,*) 'At least a few BP values are missing, please provide exactly ', &
                    NZ+1, 'AP and BP values.'
       CALL HCO_Error( HcoState%Config%Err, ErrMsg, RC, THISLOC=ThisLoc )
       RETURN
    ENDIF

    ! ------------------------------------------------------------------
    ! Close file
    ! ------------------------------------------------------------------
    CLOSE( IU_FILE )

    ! ------------------------------------------------------------------
    ! Fill grid box values
    ! ------------------------------------------------------------------
    DLAT = ( YMAX - YMIN ) / NY

    ! Now fill values
    DO J = 1, NY
    DO I = 1, NX

       ! Set longitude and latitude edge values if not read from disk
       IF ( XEDGE(I,J,1) == HCO_MISSVAL ) THEN

          ! eventually get from mid-points
          IF ( XMID(I,J,1) /= HCO_MISSVAL ) THEN
             IF ( I > 1 ) THEN
                DLON         = XMID(I,J,1) - XMID(I-1,J,1)
             ELSE
                DLON         = XMID(I+1,J,1) - XMID(I,J,1)
             ENDIF
             XEDGE(I,J,1) = XMID(I,J,1) - DLON/2.0

          ! otherwise assume constant grid spacing
          ELSE
             DLON = ( XMAX - XMIN ) / NX
             XEDGE(I,J,1) = XMIN + ( (I-1) * DLON )
          ENDIF
       ELSE
          DLON = XEDGE(I+1,J,1) - XEDGE(I,J,1)
       ENDIF

       IF ( YEDGE(I,J,1) == HCO_MISSVAL ) THEN

          ! eventually get from mid-points
          IF ( YMID(I,J,1) /= HCO_MISSVAL ) THEN
             IF ( J > 1 ) THEN
                DLAT         = YMID(I,J,1) - YMID(I,J-1,1)
             ELSE
                DLAT         = YMID(I,J+1,1) - YMID(I,J,1)
             ENDIF
             YEDGE(I,J,1) = YMID(I,J,1) - DLAT/2.0

          ! otherwise assume constant grid spacing
          ELSE
             DLAT = ( YMAX - YMIN ) / NY
             YEDGE(I,J,1) = YMIN + ( (J-1) * DLAT )
          ENDIF
       ELSE
          DLAT = YEDGE(I,J+1,1) - YEDGE(I,J,1)
       ENDIF

       ! Set mid values
       IF ( XMID(I,J,1) == HCO_MISSVAL ) THEN
          XMID(I,J,1) = XEDGE(I,J,1) + ( DLON / 2.0_hp )
       ENDIF
       IF ( YMID(I,J,1) == HCO_MISSVAL ) THEN
          YMID(I,J,1) = YEDGE(I,J,1) + ( DLAT / 2.0_hp )
       ENDIF

       ! Get sine of latitude edges
       YDGR        = PI_180 * YEDGE(I,J,1)  ! radians
       YSN         = SIN( YDGR )            ! sine
       YSIN(I,J,1) = YSN

       ! Eventually set uppermost edge
       IF ( I == NX ) THEN
          IF ( XEDGE(I+1,J,1) == HCO_MISSVAL ) THEN
             XEDGE(I+1,J,1) = XMIN + I * DLON
          ENDIF
       ENDIF
       IF ( J == NY ) THEN
          IF ( YEDGE(I,J+1,1) == HCO_MISSVAL ) THEN
             YEDGE(I,J+1,1) = YMIN + J * DLAT
          ENDIF
          YDGR           = PI_180 * YEDGE(I,J+1,1)  ! radians
          YSN            = SIN( YDGR )              ! sine
          YSIN(I,J+1,1)  = YSN
       ENDIF

    ENDDO
    ENDDO

    ! Calculate grid box areas. Follow calculation from grid_mod.F90
    ! of GEOS-Chem.
    DO J = 1, NY

       ! delta latitude
       SIN_DELTA = YSIN(1,J+1,1) - YSIN(1,J,1)

       ! Grid box area.
       AM2 = DLON * PI_180 * HcoState%Phys%Re**2 * SIN_DELTA

       ! Pass to array
       AREA_M2(:,J,1) = AM2

    ENDDO

    ! Set grid dimensions
    HcoState%NX = NX
    HcoState%NY = NY
    HcoState%NZ = NZ

    ! Vertical grid definition
    IF ( ANY(AP/=HCO_MISSVAL) ) THEN
       CALL HCO_VertGrid_Define( HcoState%Config, &
                                 HcoState%Grid%zGrid, NZ, &
                                 Ap=Ap, Bp=Bp, RC=RC )
    ELSE
       CALL HCO_VertGrid_Define( HcoState%Config, &
                                 HcoState%Grid%zGrid, NZ, RC=RC )
    ENDIF
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Set pointers to grid variables
    HcoState%Grid%XMID%Val       => XMID   (:,:,1)
    HcoState%Grid%YMID%Val       => YMID   (:,:,1)
    HcoState%Grid%XEDGE%Val      => XEDGE  (:,:,1)
    HcoState%Grid%YEDGE%Val      => YEDGE  (:,:,1)
    HcoState%Grid%YSIN%Val       => YSIN   (:,:,1)
    HcoState%Grid%AREA_M2%Val    => AREA_M2(:,:,1)

    ! The pressure edges and grid box heights are obtained from
    ! an external file in ExtState_SetFields
    HcoState%Grid%PEDGE%Val      => NULL()
    HcoState%Grid%BXHEIGHT_M%Val => NULL()
    HcoState%Grid%ZSFC%Val       => NULL()
    HcoState%Grid%PSFC%Val       => NULL()

    ! Write grid information to log-file
    WRITE(Msg,*) 'HEMCO grid definitions:'
    CALL HCO_MSG(HcoState%Config%Err,MSG)

    WRITE(MSG,*) ' --> Number of longitude cells: ', NX
    CALL HCO_MSG(HcoState%Config%Err,MSG)
    WRITE(MSG,*) ' --> Number of latitude cells : ', NY
    CALL HCO_MSG(HcoState%Config%Err,MSG)
    WRITE(MSG,*) ' --> Number of levels         : ', NZ
    CALL HCO_MSG(HcoState%Config%Err,MSG)
    WRITE(MSG,*) ' --> Lon range [deg E]        : ', XMIN, XMAX
    CALL HCO_MSG(HcoState%Config%Err,MSG)
    WRITE(MSG,*) ' --> Lat range [deg N]        : ', YMIN, YMAX
    CALL HCO_MSG(HcoState%Config%Err,MSG)

    ! Cleanup
    IF ( ALLOCATED(AP) ) DEALLOCATE(AP)
    IF ( ALLOCATED(BP) ) DEALLOCATE(BP)

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Set_Grid
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_nnMatch
!
! !DESCRIPTION: Subroutine Get\_nnMatch returns the number of species
! found in both the HEMCO configuration and the species input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_nnMatch( HcoConfig, nnMatch, RC )
!
! !USES:
!
    USE HCO_Config_Mod, ONLY : Config_GetnSpecies
    USE HCO_Config_Mod, ONLY : Config_GetSpecNames
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(  OUT)  :: nnMatch   ! Number of HEMCO species that are
                                         ! also species in the atm. model
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ConfigObj), POINTER        :: HcoConfig ! Config object
    INTEGER,         INTENT(INOUT)  :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller   - Initial Version
!  18 Jan 2019 - R. Yantosca - Improve error trapping
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER            :: AS,     IDX
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! Get_nnMatch begins here
    !=================================================================

    ! Initialize
    RC      = HCO_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
     'Get_nnMatch (in module HEMCO/Interfaces/hcoi_standalone_mod.F90)'

    ! Extract number of HEMCO species and corresponding species names
    ! as read from the HEMCO config. file.
    nHcoSpec = Config_GetnSpecies ( HcoConfig )
    CALL Config_GetSpecNames( HcoConfig, &
                              HcoSpecNames, nHcoSpec, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Config_GetSpecNames"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Extract species to be used from input file
    CALL Model_GetSpecies( HcoConfig,                           &
                           nModelSpec,     ModelSpecNames,      &
                           ModelSpecIDs,   ModelSpecMW,         &
                           ModelSpecEmMW,  ModelSpecMolecRatio, &
                           ModelSpecK0,    ModelSpecCR,         &
                           ModelSpecPKA,   RC                    )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Model_GetSpecies"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! See how many species are also used in GEOS-Chem
    ALLOCATE(matchIDx(nHcoSpec),STAT=AS)
    IF ( AS/=0 ) THEN
       ErrMsg = 'Allocation error matchIDx'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    matchIDx(:) = -1
    CALL HCO_CharMatch( HcoSpecNames,   nHcoSpec,      &
                        ModelSpecNames, nModelSpec,    &
                        matchIDx,       nnMatch         )
    IF ( nnMatch == 0 ) THEN
       ErrMsg = 'HCO_CharMatch returned found matching species!'
       CALL HCO_Error(HcoConfig%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Get_nnMatch
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_Species
!
! !DESCRIPTION: Subroutine Register\_Species registers all species in the
!  HEMCO state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_Species( HcoState, RC )
!
! !USES:
!
    USE HCO_LogFile_Mod, ONLY : HCO_SPEC2LOG
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(HCO_STATE), POINTER :: HcoState
    INTEGER, INTENT(INOUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER :: CNT, I, IDX, CID

    !=================================================================
    ! REGISTER_SPECIES begins here
    !=================================================================

    ! Loop over all possible HEMCO species
    cnt = 0
    DO I = 1, nHcoSpec

       ! Skip if this HEMCO species is not used in GEOS-Chem
       IF ( MatchIDx(I) < 0 ) CYCLE

       ! increase counter: this is the index in HcoState%Spc!
       cnt                        = cnt + 1

       ! Set species name and GEOS-Chem tracer ID
       IDX                        = ModelSpecIDs(MatchIDx(I))
       HcoState%Spc(cnt)%SpcName  = HcoSpecNames(I)
       HcoState%Spc(cnt)%ModID    = IDX

       ! Molecular weights of species & emitted species.
       HcoState%Spc(cnt)%MW_g     = ModelSpecMW(IDX)
       HcoState%Spc(cnt)%EmMW_g   = ModelSpecEmMW(IDX)

       ! Emitted molecules per molecule of species.
       HcoState%Spc(cnt)%MolecRatio = ModelSpecMolecRatio(IDX)

       ! Set Henry coefficients
       HcoState%Spc(cnt)%HenryK0  = ModelSpecK0(IDX)
       HcoState%Spc(cnt)%HenryCR  = ModelSpecCR(IDX)
       HcoState%Spc(cnt)%HenryPKA = ModelSpecPKA(IDX)

       ! Logfile I/O
       CALL HCO_SPEC2LOG( HcoState, Cnt )

    ENDDO !I

    CALL HCO_MSG(HcoState%Config%Err,SEP1='-')

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Register_Species
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Define_Diagnostics
!
! !DESCRIPTION: Subroutine Define\_Diagnostics defines all diagnostics to be
!  used in this simulation.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Define_Diagnostics( HcoState, RC, SetDefault )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,   ONLY : GetExtNr
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_STATE),  POINTER                   :: HcoState
    LOGICAL,          INTENT(IN   ), OPTIONAL   :: SetDefault  ! Define default diagnostics?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)             :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller - Initial Version
!  05 Feb 2015 - C. Keller - Added SetDefault flag
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: SetDf
    INTEGER            :: I, N, ExtNr, HcoID

    ! Strings
    CHARACTER(LEN=31)  :: DiagnName
    CHARACTER(LEN=255) :: ErrMsg,   ThisLoc

    !=================================================================
    ! DEFINE_DIAGNOSTICS begins here
    !=================================================================

    ! Initialize
    RC      =  HCO_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
     'DEFINE_DIAGNOSTICS (in module HEMCO/Interfaces/hcoi_standalone_mod.F90'

    ! Get number of diagnostics currently defined in the default
    ! collection
    CALL DiagnCollection_Get( HcoState%Diagn, &
       HcoState%Diagn%HcoDiagnIDDefault, nnDiagn=N, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "DiagnCollection_Get"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    print*, '### Define_Diagnostics: NNDIAGN: ', N

    ! If there are no diagnostics defined yet, define some default
    ! diagnostics below. These are simply the overall emissions
    ! (across all extensions, categories, hierarchies) for each
    ! HEMCO species.
    IF ( PRESENT(SetDefault) ) THEN
       SetDf = SetDefault
    ELSE
       SetDf = ( N == 0 )
    ENDIF
    IF ( SetDf ) THEN

       ! Loop over all HEMCO species
       DO I = 1, HcoState%nSpc

          ! Get HEMCO ID
          HcoID = HcoState%Spc(I)%HcoID
          IF ( HcoID <= 0 ) CYCLE

          ! Create diagnostics
          DiagnName = 'HEMCO__EMIS_' // TRIM(HcoState%Spc(I)%SpcName)
          CALL Diagn_Create ( HcoState,                               &
                              cName     = DiagnName,                  &
                              ExtNr     = -1,                         &
                              Cat       = -1,                         &
                              Hier      = -1,                         &
                              HcoID     = HcoID,                      &
                              SpaceDim  = 3,                          &
                              LevIDx    = -1,                         &
                              OutUnit   = 'kg/m2/s',                  &
                              AutoFill  = 1,                          &
                              COL = HcoState%Diagn%HcoDiagnIDDefault, &
                              OkIfExist = .TRUE.,                     &
                              RC        = RC                           )

          ! Trap potential errors
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Error defining diagnostic: ' // TRIM( DiagnName )
             CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO !I
    ENDIF

    !--------------------------------------------------------------------------
    ! Define some additional diagnostics
    !--------------------------------------------------------------------------
    ExtNr = GetExtNr ( HcoState%Config%ExtList, 'LightNOx' )
    IF ( ExtNr > 0 ) THEN

       ! Loop over lighthing flash quantities
       DO I = 1, 3

          ! Pick the proper diagnostic name
          SELECT CASE( I )
             CASE( 1 )
                DiagnName = 'LIGHTNING_TOTAL_FLASHRATE'
             CASE( 2 )
                DiagnName = 'LIGHTNING_INTRACLOUD_FLASHRATE'
             CASE( 3 )
                DiagnName = 'LIGHTNING_CLOUDGROUND_FLASHRATE'
          END SELECT

          ! Define diagnostics ID
          N = 56000 + I

          ! Create diagnostic container
          CALL Diagn_Create( HcoState,                               &
                             cName     = TRIM( DiagnName ),          &
                             cID       = N,                          &
                             ExtNr     = ExtNr,                      &
                             Cat       = -1,                         &
                             Hier      = -1,                         &
                             HcoID     = -1,                         &
                             SpaceDim  = 2,                          &
                             LevIDx    = -1,                         &
                             OutUnit   = 'flashes/min/km2',          &
                             OutOper   = 'Mean',                     &
                             COL = HcoState%Diagn%HcoDiagnIDDefault, &
                             AutoFill  = 0,                          &
                             RC        = RC                           )

          ! Trap potential errors
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Error defining diagnostic: ' // TRIM( DiagnName )
             CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

       ! ----------------------------------------------------------
       ! Diagnostics for convective cloud top height.
       ! ----------------------------------------------------------

       ! Define diagnostics name and ID
       DiagnName = 'LIGHTNING_CLOUD_TOP'
       N         = 56004

       ! Create diagnostic container
       CALL Diagn_Create( HcoState,                               &
                          cName     = TRIM( DiagnName ),          &
                          cID       = N,                          &
                          ExtNr     = ExtNr,                      &
                          Cat       = -1,                         &
                          Hier      = -1,                         &
                          HcoID     = -1,                         &
                          SpaceDim  = 2,                          &
                          LevIDx    = -1,                         &
                          OutUnit   = '1',                        &
                          OutOper   = 'Mean',                     &
                          COL = HcoState%Diagn%HcoDiagnIDDefault, &
                          AutoFill  = 0,                          &
                          RC        = RC                           )

       ! Trap potential errors
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Error defining diagnostic: ' // TRIM( DiagnName )
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF ! Lightning NOx

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Define_Diagnostics
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Read_Time
!
! !DESCRIPTION: Subroutine READ\_TIME reads the time information for the
!  HEMCO standalone from an input file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Read_Time( HcoState, RC )
!
! !USES:
!
    USE inquireMod,      ONLY : findfreeLUN
    USE HCO_Extlist_Mod, ONLY : HCO_GetOpt, GetExtOpt, CoreNr
!
! !INPUT PARAMETERS:
!
    TYPE(HCO_State), POINTER       :: HcoState
!
! !INPUT/OUTPUT PARAMETERS
!
    INTEGER,         INTENT(INOUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: AS, IOS, IU_FILE
    INTEGER             :: I,  N,   LNG, LOW
    LOGICAL             :: EOF, FOUND

    ! Strings
    CHARACTER(LEN=  1)  :: COL
    CHARACTER(LEN=255)  :: ErrMsg, ThisLoc, DUM
    CHARACTER(LEN=255)  :: MyTimeFile

    !=================================================================
    ! READ_TIME begins here
    !=================================================================

    ! Initialize
    RC      = HCO_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
     'READ_TIME (in module HEMCO/Standalone/hcoi_standalone_mod.F90)'

    ! Try to get TimeFile from configuration file (in settings)
    CALL GetExtOpt ( HcoState%Config, CoreNr, 'TimeFile', &
                     OptValChar=MyTimeFile,   Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Hco_Run"!'
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    IF ( FOUND ) TimeFile = MyTimeFile

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Write colon character to local variable
    COL = HCO_GetOpt( HcoState%Config%ExtList, 'Colon' )

    ! Open time file
    OPEN( IU_FILE, FILE=TRIM(TimeFile), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       ErrMsg = 'Error 1 reading ' // TRIM(TimeFile)
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Read start and end of simulation
    DO N = 1,2

       CALL GetNextLine( IU_FILE, DUM, EOF, RC )
       IF ( RC /= HCO_SUCCESS .OR. EOF ) THEN
          ErrMsg = 'Error reading time in ' // TRIM(TimeFile)
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Remove 'BEGIN: ' or 'END: ' at the beginning
       LNG = LEN(TRIM(DUM))
       LOW = NextCharPos ( TRIM(DUM), COL, 1 )
       IF ( LOW < 0 .OR. LOW == LNG ) THEN
          ErrMsg = 'Cannot extract index after colon: ' // TRIM(DUM)
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       LOW = LOW + 1
       DUM = ADJUSTL(DUM(LOW:LNG))
       LNG = LEN(TRIM(DUM))

       ! Times have to be stored as:
       ! YYYY-MM-DD HH:MM:SS
       ! --> read year from position 1:4, month from 6:7, etc.
       IF ( LNG /= 19 ) THEN
          ErrMsg = 'Provided time stamp is not `YYYY-MM-DD HH:MM:SS`! ' // &
                   TRIM(DUM)
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       READ ( DUM( 1: 4), * ) YRS(N)
       READ ( DUM( 6: 7), * ) MTS(N)
       READ ( DUM( 9:10), * ) DYS(N)
       READ ( DUM(12:13), * ) HRS(N)
       READ ( DUM(15:16), * ) MNS(N)
       READ ( DUM(18:19), * ) SCS(N)

    ENDDO !I

    ! Get emission time step
    CALL GetNextLine( IU_FILE, DUM, EOF, RC )
    IF ( (RC /= HCO_SUCCESS) .OR. EOF ) THEN
       ErrMsg = 'Cannot read emission time step from ' // TRIM(TimeFile)
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Get index after colon
    LNG = LEN(TRIM(DUM))
    LOW = NextCharPos ( TRIM(DUM), COL, 1 )
    IF ( LOW < 0 .OR. LOW == LNG ) THEN
       ErrMsg = 'Cannot extract index after colon: ' // TRIM(DUM)
       CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    LOW = LOW + 1
    READ( DUM(LOW:LNG), * ) HcoState%TS_EMIS

    ! Set same chemical and dynamic time step
    HcoState%TS_CHEM = HcoState%TS_EMIS
    HcoState%TS_DYN  = HcoState%TS_EMIS

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Read_Time
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtState_SetFields
!
! !DESCRIPTION: Subroutine ExtState\_SetFields fills the ExtState data fields
! with data read through the HEMCO configuration file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtState_SetFields ( HcoState, ExtState, RC )
!
! !USES:
!
    USE HCO_ARR_MOD,        ONLY : HCO_ArrAssert
    USE HCO_GEOTOOLS_MOD,   ONLY : HCO_GetSUNCOS
    USE HCO_GEOTOOLS_MOD,   ONLY : HCO_CalcVertGrid
    USE HCOX_STATE_MOD,     ONLY : ExtDat_Set
    USE HCO_CLOCK_MOD,      ONLY : HcoClock_First
!
! !INPUT/OUTPUT PARAMETERS
!
    TYPE(HCO_STATE), POINTER       :: HcoState
    TYPE(EXT_STATE), POINTER       :: ExtState
    INTEGER,         INTENT(INOUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  28 Jul 2014 - C. Keller   - Initial Version
!  06 Oct 2014 - M. Sulprizio- Remove PCENTER. Now calculate from pressure edges
!  09 Jul 2015 - E. Lundgren - Add MODIS Chlorophyll-a (CHLR)
!  26 Oct 2016 - R. Yantosca - Don't nullify local ptrs in declaration stmts
!  15 Jan 2019 - R. Yantosca - Update met field names to be consistent with
!                              those used for the FlexGrid update
!  18 Jan 2019 - R. Yantosca - Improve error trapping
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: FIRST

    ! Strings
    CHARACTER(LEN=255) :: Name, ErrMsg, ThisLoc

    ! Pointers
    REAL(hp), POINTER  :: PSFC    (:,:  )
    REAL(hp), POINTER  :: ZSFC    (:,:  )
    REAL(hp), POINTER  :: TK      (:,:,:)
    REAL(hp), POINTER  :: BXHEIGHT(:,:,:)
    REAL(hp), POINTER  :: PEDGE   (:,:,:)

    !========================================================================
    ! ExtState_SetFields begins here
    !========================================================================

    ! Initialize
    RC       = HCO_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
     'ExtState_SetFields (in HEMCO/Interfaces/hcoi_standalone_mod.F90'

    ! Nullify pointers
    PSFC     => NULL()
    ZSFC     => NULL()
    TK       => NULL()
    BXHEIGHT => NULL()
    PEDGE    => NULL()

    ! Enter
    CALL HCO_ENTER( HcoState%Config%Err, ThisLoc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! First call?
    FIRST = HcoClock_First ( HcoState%Clock, .FALSE. )

    !------------------------------------------------------------------------
    ! %%%%% 2D fields %%%%%
    ! (1) Now use the same met field names as are specified in the
    !     the HEMCO_Config.rc file for the "FlexGrid" update
    ! (2) Not all extension fields are used for a given simulation type
    !------------------------------------------------------------------------

    !%%%%% 10-m winds %%%%%
    IF ( ExtState%U10M%DoUse ) THEN
       Name = 'U10M'
       CALL ExtDat_Set( HcoState,     ExtState%U10M,                         &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                    '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%V10M%DoUse ) THEN
       Name = 'V10M'
       CALL ExtDat_Set( HcoState,     ExtState%V10M,                         &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Albedo %%%%%
    IF ( ExtState%ALBD%DoUse ) THEN
       Name = 'ALBEDO'
       CALL ExtDat_Set( HcoState,     ExtState%ALBD,                         &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Land-water-ice flags  %%%%%
    IF ( ExtState%WLI%DoUse ) THEN
       Name = 'LWI'
       CALL ExtDat_Set( HcoState,     ExtState%WLI,                          &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Air and skin temperature %%%%%
    IF ( ExtState%T2M%DoUse ) THEN
       Name = 'T2M'
       CALL ExtDat_Set( HcoState,     ExtState%T2M,                          &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%TSKIN%DoUse ) THEN
       Name = 'TS'
       CALL ExtDat_Set( HcoState,     ExtState%TSKIN,                        &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Soil moisture %%%%%
    IF ( ExtState%GWETROOT%DoUse ) THEN
       Name = 'GWETROOT'
       CALL ExtDat_Set( HcoState,     ExtState%GWETROOT,                     &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%GWETTOP%DoUse ) THEN
       Name = 'GWETTOP'
       CALL ExtDat_Set( HcoState,     ExtState%GWETTOP,                      &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Snow fields %%%%%
    IF ( ExtState%SNOWHGT%DoUse ) THEN
       Name = 'SNOMAS'
       CALL ExtDat_Set( HcoState,     ExtState%SNOWHGT,                      &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%SNODP%DoUse ) THEN
       Name = 'SNODP'
       CALL ExtDat_Set( HcoState,     ExtState%SNODP,                        &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
      ENDIF
    ENDIF

    !%%%%% Friction velocity %%%%%
    IF ( ExtState%USTAR%DoUse ) THEN
       Name = 'USTAR'
       CALL ExtDat_Set( HcoState,     ExtState%USTAR,                        &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Roughness height %%%%%
    IF ( ExtState%Z0%DoUse ) THEN
       Name = 'Z0M'
       CALL ExtDat_Set( HcoState,     ExtState%Z0,                           &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg , RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Tropopause pressure %%%%%
    IF ( ExtState%TROPP%DoUse ) THEN
       Name = 'TROPPT'
       CALL ExtDat_Set( HcoState,     ExtState%TROPP,                        &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% PAR direct and diffuse %%%%%
    IF ( ExtState%PARDR%DoUse ) THEN
       Name = 'PARDR'
       CALL ExtDat_Set( HcoState,     ExtState%PARDR,                        &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%PARDF%DoUse ) THEN
       Name = 'PARDF'
       CALL ExtDat_Set( HcoState,     ExtState%PARDF,                        &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%RADSWG%DoUse ) THEN
       Name = 'SWGDN'
       CALL ExtDat_Set( HcoState,     ExtState%RADSWG,                       &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Cloud fraction @ surface %%%%%
    IF ( ExtState%CLDFRC%DoUse ) THEN
       Name = 'CLDTOT'
       CALL ExtDat_Set( HcoState,     ExtState%CLDFRC,                       &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Leaf area index %%%%%
    IF ( ExtState%LAI%DoUse ) THEN
       Name = 'LAI'
       CALL ExtDat_Set( HcoState,     ExtState%LAI,                          &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Flash density %%%%%
    IF ( ExtState%FLASH_DENS%DoUse ) THEN
       Name = 'FLASH_DENS'
       CALL ExtDat_Set( HcoState,     ExtState%FLASH_DENS,                   &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Convective depth %%%%%
    IF ( ExtState%CONV_DEPTH%DoUse ) THEN
       Name = 'CONV_DEPTH'
       CALL ExtDat_Set( HcoState,     ExtState%CONV_DEPTH,                   &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Fractional coverage fields %%%%%
    IF ( ExtState%FRCLND%DoUse ) THEN
       Name = 'FRCLND'
       CALL ExtDat_Set( HcoState,     ExtState%FRCLND,                       &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%FRLAND%DoUse ) THEN
       Name = 'FRLAND'
       CALL ExtDat_Set( HcoState,     ExtState%FRLAND,                       &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%FROCEAN%DoUse ) THEN
       Name = 'FROCEAN'
       CALL ExtDat_Set( HcoState,     ExtState%FROCEAN,                      &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%FRLAKE%DoUse ) THEN
       Name = 'FRLAKE'
       CALL ExtDat_Set( HcoState,     ExtState%FRLAKE,                       &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%FRLANDIC%DoUse ) THEN
       Name = 'FRLANDIC'
       CALL ExtDat_Set( HcoState,     ExtState%FRLANDIC,                     &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Solar zenith angle %%%%%
    IF ( ExtState%SZAFACT%DoUse ) THEN
       Name = 'SZAFACT'
       CALL ExtDat_Set( HcoState,     ExtState%SZAFACT,                      &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Photolysis values %%%%%
    IF ( ExtState%JNO2%DoUse ) THEN
       Name = 'JNO2'
       CALL ExtDat_Set( HcoState,     ExtState%JNO2,                         &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
      IF ( RC == HCO_SUCCESS ) THEN
         ErrMsg = 'Could not find quantity "' // TRIM( Name )             // &
                  '" for the HEMCO standalone simulation!'
         CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
         CALL HCO_Leave( HcoState%Config%Err, RC )
         RETURN
      ENDIF
   ENDIF

   IF ( ExtState%JOH%DoUse ) THEN
      Name = 'JOH'
      CALL ExtDat_Set( HcoState,     ExtState%JOH,                           &
                       TRIM( Name ), RC,       FIRST=FIRST                  )
      IF ( RC == HCO_SUCCESS ) THEN
         ErrMsg = 'Could not find quantity "' // TRIM( Name )             // &
                  '" for the HEMCO standalone simulation!'
         CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
         CALL HCO_Leave( HcoState%Config%Err, RC )
         RETURN
      ENDIF
   ENDIF

    !-----------------------------------------------------------------
    ! %%%%% 3D fields %%%%%
    ! (1) Now use the same met field names as are specified in the
    !     the HEMCO_Config.rc file for the "FlexGrid" update
    ! (2) Not all extension fields are used for a given simulation type
    !-----------------------------------------------------------------

    !%%%%% Cloud convection mass flux %%%%%
    IF ( ExtState%CNV_MFC%DoUse ) THEN
       Name = 'CMFMC'
       CALL ExtDat_Set( HcoState,     ExtState%CNV_MFC,                      &
                        TRIM( Name ), RC,       FIRST=FIRST,                 &
                        OnLevEdge=.TRUE.                                    )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Specific humidity %%%%%
    IF ( ExtState%SPHU%DoUse ) THEN
       Name = 'SPHU1'
       CALL ExtDat_Set( HcoState,     ExtState%SPHU,                         &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Temperature %%%%%
    IF ( ExtState%TK%DoUse ) THEN
       Name = 'TMPU1'
       CALL ExtDat_Set( HcoState,     ExtState%TK,                           &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Air mass, volume, density etc fields %%%%%
    IF ( ExtState%AIR%DoUse ) THEN
       Name = 'AIR'
       CALL ExtDat_Set( HcoState,     ExtState%AIR,                          &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC == HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%AIRVOL%DoUse ) THEN
       Name = 'AIRVOL'
       CALL ExtDat_Set( HcoState,     ExtState%AIRVOL,                       &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC == HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%AIRDEN%DoUse ) THEN
       Name = 'AIRDEN'
       CALL ExtDat_Set( HcoState,     ExtState%AIRDEN,                       &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC == HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Concentration fields %%%%%
    IF ( ExtState%O3%DoUse ) THEN
       Name = 'O3'
       CALL ExtDat_Set( HcoState,     ExtState%O3,                           &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC == HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%NO%DoUse ) THEN
       Name = 'NO'
       CALL ExtDat_Set( HcoState,     ExtState%NO,                           &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC == HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%NO2%DoUse ) THEN
       Name = 'NO2'
       CALL ExtDat_Set( HcoState,     ExtState%NO2,                          &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC == HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%HNO3%DoUse ) THEN
       Name = 'HNO3'
       CALL ExtDat_Set( HcoState,     ExtState%HNO3,                         &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC == HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Deposition fields (for soil NOx) %%%%%
    IF ( ExtState%DRY_TOTN%DoUse ) THEN
       Name = 'DRY_TOTN'
       CALL ExtDat_Set( HcoState,     ExtState%DRY_TOTN,                     &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC == HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    IF ( ExtState%WET_TOTN%DoUse ) THEN
       Name = 'WET_TOTN'
       CALL ExtDat_Set( HcoState,     ExtState%WET_TOTN,                     &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC == HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !%%%%% Fraction of PBL field (for sea exchange only) %%%%%
    IF ( ExtState%FRAC_OF_PBL%DoUse ) THEN
       Name = 'FRAC_OF_PBL'
       CALL ExtDat_Set( HcoState,     ExtState%FRAC_OF_PBL,                  &
                        TRIM( Name ), RC,       FIRST=FIRST                 )
       IF ( RC == HCO_SUCCESS ) THEN
          ErrMsg = 'Could not find quantity "' // TRIM( Name )            // &
                   '" for the HEMCO standalone simulation!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! ==> DRYCOEFF must be read from the configuration file in module
    !     hcox_soilnox_mod.F90.
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! Check for vertical grid update. This will try to read the
    ! vertical grid quantities from disk or calculate them from other
    ! quantities read from disk.
    !-----------------------------------------------------------------

    ! Eventually get temperature from disk
    IF ( ExtState%TK%DoUse ) TK => ExtState%TK%Arr%Val

    ! Attempt to calculate vertical grid quantities
    CALL HCO_CalcVertGrid( HcoState, PSFC, ZSFC, TK, BXHEIGHT, PEDGE, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Reset pointers
    PSFC     => NULL()
    ZSFC     => NULL()
    TK       => NULL()
    BXHEIGHT => NULL()
    PEDGE    => NULL()

    !-----------------------------------------------------------------
    ! If needed, calculate SUNCOS values
    !-----------------------------------------------------------------
    IF ( ExtState%SUNCOS%DoUse ) THEN
       IF ( FIRST ) THEN
          CALL HCO_ArrAssert( ExtState%SUNCOS%Arr, HcoState%NX, HcoState%NY, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'SUNCOS array is not the expected dimensions!'
             CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
             CALL HCO_Leave( HcoState%Config%Err, RC )
             RETURN
          ENDIF
       ENDIF

       CALL HCO_GetSUNCOS( HcoState, ExtState%SUNCOS%Arr%Val, 0, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Error encountered in routine "HCO_GetSuncos"!'
          CALL HCO_Error( HcoConfig%Err, ErrMsg, RC, ThisLoc )
          CALL HCO_Leave( HcoState%Config%Err, RC )
          RETURN
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! All done
    !-----------------------------------------------------------------

    ! Not first call any more
    FIRST = .FALSE.

    ! Leave w/ success
    CALL HCO_Leave( HcoState%Config%Err, RC )

  END SUBROUTINE ExtState_SetFields
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ExtState_UpdateFields
!
! !DESCRIPTION: Subroutine ExtState\_UpdateFields makes sure that all local
! variables that ExtState is pointing to are up to date. For the moment, this
! is just a placeholder routine as none of the ExtState fields is filled by
! local module fields. Content can be added to it if there are variables that
! need to be updated manually, e.g. not through netCDF input data.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ExtState_UpdateFields ( HcoState, ExtState, RC )
!
! !INPUT/OUTPUT PARAMETERS
!
    TYPE(HCO_STATE),  POINTER       :: HcoState
    TYPE(EXT_STATE),  POINTER       :: ExtState
    INTEGER,          INTENT(INOUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  28 Jul 2014 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!

    !=================================================================
    ! ExtState_UpdateFields begins here
    !=================================================================

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE ExtState_UpdateFields
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IsEndOfSimulation
!
! !DESCRIPTION: Function IsEndOfSimulation returns true if the passed date
! is beyond the end of the simulation date.
!\\
!\\
! !INTERFACE:
!
  FUNCTION IsEndOfSimulation( Yr, Mt, Dy, Hr, Mn, Sc ) RESULT ( IsEnd )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    INTEGER,          INTENT(IN   ) :: YR
    INTEGER,          INTENT(IN   ) :: MT
    INTEGER,          INTENT(IN   ) :: DY
    INTEGER,          INTENT(IN   ) :: HR
    INTEGER,          INTENT(IN   ) :: MN
    INTEGER,          INTENT(IN   ) :: SC
!
! !OUTPUT PARAMETERS
!
    LOGICAL                         :: IsEnd
!
! !REVISION HISTORY:
!  08 Sep 2014 - C. Keller - Initial Version
!  13 Jul 2015 - C. Keller - Bug fix: now save YYYYMMDD and hhmmss in different
!                            variables to avoid integer truncation errors.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER       :: THISYYYYMMDD
    INTEGER       :: THIShhmmss
    INTEGER, SAVE :: ENDYYYYMMDD = -1
    INTEGER, SAVE :: ENDhhmmss   = -1

    !=================================================================
    ! IsEndOfSimulation begins here
    !=================================================================

    ! Init
    IsEnd = .FALSE.

    ! Calculate simulation end datetime if not yet done so
    IF ( ENDYYYYMMDD < 0 ) THEN
       ENDYYYYMMDD = YRS(2)*10000 + MTS(2)*100 + DYS(2)
       ENDhhmmss   = HRS(2)*10000 + MNS(2)*100 + SCS(2)
    ENDIF

    ! Calculate current datetime
    THISYYYYMMDD = YR*10000 + MT*100 + DY
    THIShhmmss   = HR*10000 + MN*100 + SC

    ! Check if current datetime is beyond simulation end date
    IF ( THISYYYYMMDD > ENDYYYYMMDD ) THEN
       IsEnd = .TRUE.
    ELSEIF ( (THISYYYYMMDD == ENDYYYYMMDD) .AND. (THIShhmmss >= ENDhhmmss) ) THEN
       IsEnd = .TRUE.
    ENDIF

  END FUNCTION IsEndOfSimulation
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCOI_Sa_InitCleanup
!
! !DESCRIPTION: deallocates all local species arrays used during initialization.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCOI_SA_InitCleanup ( RC )
!
! !INPUT/OUTPUT PARAMETERS
!
    INTEGER, INTENT(INOUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  04 Feb 2016 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    !=================================================================
    ! HCOI_SA_InitCleanup begins here
    !=================================================================

    ! Deallocate local variables (not used anymore)
    IF ( ASSOCIATED(ModelSpecNames     ) ) DEALLOCATE(ModelSpecNames     )
    IF ( ASSOCIATED(ModelSpecIDs       ) ) DEALLOCATE(ModelSpecIDs       )
    IF ( ASSOCIATED(ModelSpecMW        ) ) DEALLOCATE(ModelSpecMW        )
    IF ( ASSOCIATED(ModelSpecEmMW      ) ) DEALLOCATE(ModelSpecEmMW      )
    IF ( ASSOCIATED(ModelSpecMolecRatio) ) DEALLOCATE(ModelSpecMolecRatio)
    IF ( ASSOCIATED(ModelSpecK0        ) ) DEALLOCATE(ModelSpecK0        )
    IF ( ASSOCIATED(ModelSpecCR        ) ) DEALLOCATE(ModelSpecCR        )
    IF ( ASSOCIATED(ModelSpecPKA       ) ) DEALLOCATE(ModelSpecPKA       )
    IF ( ASSOCIATED(matchIDx           ) ) DEALLOCATE(matchIDx           )
    IF ( ASSOCIATED(HcoSpecNames       ) ) DEALLOCATE(HcoSpecNames       )

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE HCOI_SA_InitCleanup
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Dry_Run
!
! !DESCRIPTION: Looks at the input arguments to determine if the user
!  has selected to do a GEOS-Chem dry-run.  If so, then the proper
!  fields of Input\_Opt will be populated accordingly, and the dry-run
!  log file will be opened.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Dry_Run( IsDryRun, RC )
!
! !USES:
!
    USE InquireMod, ONLY : FindFreeLUN
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: IsDryRun      ! Is it a dry-run?
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC            ! Success or failure?
!
! !REMARKS:
!  If in a "dry-run" mode, HEMCO will simply check whether files
!  are present (and possibly in the correct format) and go through
!  time-steps to check met fields and other IO issues.
!  No actual "compute" is performed.
!
! !REVISION HISTORY:
!  13 Nov 2019 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: nArg,   ArgLen

    ! Strings
    CHARACTER(LEN=255) :: ArgVal, ErrMsg, ThisLoc

    !=======================================================================
    ! Init_Dry_Run begins here!
    !=======================================================================

    ! Initialize
    RC      = HCO_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
        ' -> at Init_Dry_Run (in HEMCO/Interfaces/hcoi_standalone_mod.F90)'

    ! Enter
    CALL HCO_Enter( HcoState%Config%Err, ThisLoc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !=======================================================================
    ! Initialize dry-run fields of the HEMCO state object
    !=======================================================================
    IF ( IsDryRun ) THEN

       !--------------------------------------------------------------------
       ! If HEMCO is running in dry-run mode:
       !
       ! (1) Define dry-run parameters in HEMCO state
       ! (2) Print a warning to both to stdout and the HEMCO log file
       !--------------------------------------------------------------------

       ! Set parameters
       HcoState%Options%IsDryRun = IsDryRun

       ! Print dry-run header to stdout
       CALL Print_Dry_Run_Warning( 6 )

       ! Print dry-run header to the HEMCO log file
       CALL Print_Dry_Run_Warning( HcoState%Config%Err%LUN )

    ELSE

       !--------------------------------------------------------------------
       ! If this is a regular HEMCO standalone simuation,
       ! then set HEMCO dry-run parameters to default (off) values
       !--------------------------------------------------------------------
       HcoState%Options%IsDryRun = .FALSE.

    ENDIF

    ! Leave
    CALL HCO_Leave( HcoState%Config%Err, RC )

  END SUBROUTINE Init_Dry_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Dry_Run
!
! !DESCRIPTION: Looks at the input arguments to determine if the user
!  has selected to do a GEOS-Chem dry-run.  If so, then the proper
!  fields of Input\_Opt will be populated accordingly.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Dry_Run( RC )
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  13 Nov 2019 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255 ) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Cleanup_Dry_Run begins here!
    !=======================================================================

    ! Initialize
    RC     = HCO_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
       ' -> at Cleanup_Dry_Run (in HEMCO/Interfaces/hcoi_standalone_mod.F90)'

    ! Enter
    CALL HCO_Enter( HcoState%Config%Err, ThisLoc, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Only do the following for the dry-run simulation
    IF ( HcoState%Options%IsDryRun ) THEN

       ! Print dry-run header to stdout
       CALL Print_Dry_Run_Warning( 6 )

       ! Print dry-run header to the HEMCO log file
       CALL Print_Dry_Run_Warning( HcoState%Config%Err%LUN )

    ENDIF

    ! Leave
    CALL HCO_Leave( HcoState%Config%Err, RC )

  END SUBROUTINE Cleanup_Dry_Run
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_Dry_Run_Warning
!
! !DESCRIPTION: Prints the warning for the GEOS-Chem dry run to either
!  stdout (aka the GC log file) and the dry-run log file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_Dry_Run_Warning( U )
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: U

    !=================================================================
    ! Print warning info to the desired file
    !=================================================================
    WRITE( U, 100 )
    WRITE( U, 100 ) REPEAT( '!', 79 )
    WRITE( U, 100 ) '!!! HEMCO-STANDALONE IS IN DRY-RUN MODE!'
    WRITE( U, 100 ) '!!!'
    WRITE( U, 100 ) '!!! You will NOT get output for this run!'
    WRITE( U, 100 ) '!!! Use this command to validate a '                 // &
                     'HEMCO-STANDALONE run configuration:'
    WRITE( U, 100 ) '!!!    ./hemco_standalone.x -c CONFIG_FILE '         // &
                    '--dryrun > log'
    WRITE( U, 100 ) '!!!'
    WRITE( U, 100 ) '!!! REMOVE THE --dryrun ARGUMENT FROM THE COMMAND '  // &
                    'LINE'
    WRITE( U, 100 ) '!!! BEFORE RUNNING A HEMCO-STANDALONE PRODUCTION '   // &
                    'SIMULATION!'
    WRITE( U, 100 ) REPEAT( '!', 79 )
    WRITE( U, 120 ) '!!! Start Date       : ', YRS(1), MTS(1), DYS(1),       &
                                               HRS(1), MNS(1), SCS(1)
    WRITE( U, 120 ) '!!! End Date         : ', YRS(2), MTS(2), DYS(2),       &
                                               HRS(2), MNS(2), SCS(2)
    WRITE( U, 110 ) '!!! Meteorology      : ', TRIM(HcoState%Config%MetField)
    WRITE( U, 110 ) '!!! Grid Resolution  : ', TRIM(HcoState%Config%GridRes )
    WRITE( U, 100 ) REPEAT( '!', 79 )
    WRITE( U, 100 )

    ! Format statements
100 FORMAT( a                             )
110 FORMAT( a, a                          )
120 FORMAT( a, i4.4, 2(i2.2), 1x, 3(i2.2) )

  END SUBROUTINE Print_Dry_Run_Warning
!EOC
END MODULE HCOI_StandAlone_Mod
