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
  USE HCO_Error_MOD
  USE HCO_Diagn_Mod
  USE HCO_CharTools_Mod
  USE HCOX_State_Mod,      ONLY : Ext_State 
  USE HCO_State_Mod,       ONLY : HCO_State

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: HCOI_StandAlone_Run
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: HCOI_SA_Init
  PRIVATE :: HCOI_SA_Run
  PRIVATE :: HCOI_SA_Final
  PRIVATE :: Model_GetSpecies 
  PRIVATE :: Set_Grid
  PRIVATE :: Get_nnMatch 
  PRIVATE :: Register_Species
  PRIVATE :: Define_Diagnostics 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Default values for HEMCO input files: contain definitions of 
  ! species, grid, and time settings.
  CHARACTER(LEN=255)             :: GridFile  = 'HEMCO_sa_Grid'
  CHARACTER(LEN=255)             :: SpecFile  = 'HEMCO_sa_Spec'
  CHARACTER(LEN=255)             :: TimeFile  = 'HEMCO_sa_Time'

  ! HEMCO state 
  TYPE(HCO_State),       POINTER :: HcoState  => NULL()

  ! HEMCO extensions state
  TYPE(Ext_State),       POINTER :: ExtState  => NULL()

  ! Pointers used during initialization (for species matching)
  INTEGER                        :: nHcoSpec
  CHARACTER(LEN= 31),    POINTER :: HcoSpecNames(:) => NULL() 
  INTEGER                        :: nModelSpec
  CHARACTER(LEN= 31),    POINTER :: ModelSpecNames(:) => NULL()
  INTEGER,               POINTER :: ModelSpecIDs  (:) => NULL()
  REAL(hp),              POINTER :: ModelSpecMW   (:) => NULL()
  REAL(hp),              POINTER :: ModelSpecEmMW (:) => NULL()
  REAL(hp),              POINTER :: ModelSpecMolecRatio(:) => NULL()
  REAL(hp),              POINTER :: ModelSpecK0   (:) => NULL()
  REAL(hp),              POINTER :: ModelSpecCR   (:) => NULL()
  REAL(hp),              POINTER :: ModelSpecPKA  (:) => NULL()
  INTEGER,               POINTER :: matchidx(:) => NULL()

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
  SUBROUTINE HCOI_StandAlone_Run( ConfigFile )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)  :: ConfigFile
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL   :: am_I_Root
    INTEGER   :: RC

    !=================================================================
    ! HCOI_STANDALONE_RUN begins here!
    !=================================================================

    ! Treat as root CPU
    am_I_Root = .TRUE.

    ! Initialize the HEMCO standalone
    CALL HCOI_Sa_Init( am_I_Root, TRIM(ConfigFile), RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(*,*) 'Error in HCOI_SA_INIT'
       RETURN
    ENDIF

    ! Run the HEMCO standalone
    CALL HCOI_Sa_Run( am_I_Root, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       WRITE(*,*) 'Error in HCOI_SA_RUN'
       RETURN
    ENDIF

    ! Finalize the HEMCO standalone 
    CALL HCOI_Sa_Final( am_I_Root )
    
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
  SUBROUTINE HCOI_SA_Init( am_I_Root, ConfigFile, RC ) 
!
! !USES:
!
    USE HCO_Config_Mod,    ONLY : Config_ReadFile
    USE HCO_State_Mod,     ONLY : HcoState_Init
    USE HCO_Driver_Mod,    ONLY : HCO_Init
    USE HCOX_Driver_Mod,   ONLY : HCOX_Init
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN   ) :: am_I_Root  ! root CPU?
    CHARACTER(LEN=*),   INTENT(IN   ) :: ConfigFile ! Configuration file
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,            INTENT(INOUT) :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: nnMatch
    CHARACTER(LEN=255) :: LOC

    !=================================================================
    ! HCOI_SA_INIT begins here!
    !=================================================================

    ! Error handling 
    LOC = 'HCOI_SA_Init (hco_standalone_mod.F90)'

    !=================================================================
    ! Read HEMCO configuration file and save into buffer. This also
    ! sets the HEMCO error properties (verbose mode? log file name, 
    ! etc.) based upon the specifications in the configuration file.
    !=================================================================
    CALL Config_ReadFile( am_I_Root, TRIM(ConfigFile), 0, RC )
    IF( RC /= HCO_SUCCESS) RETURN 

    !=================================================================
    ! Open logfile 
    !=================================================================
    IF ( am_I_Root ) THEN
       CALL HCO_LogFile_Open( RC=RC ) 
       IF (RC /= HCO_SUCCESS) RETURN 
    ENDIF

    !=================================================================
    ! Initialize HEMCO state object and populate it 
    !=================================================================

    !-----------------------------------------------------------------
    ! Extract species to use in HEMCO 
    CALL Get_nnMatch( am_I_Root, nnMatch, RC )
    IF(RC /= HCO_SUCCESS) RETURN 

    !-----------------------------------------------------------------
    ! Initialize HCO state. Use only species that are used
    ! in GEOS-Chem and are also found in the HEMCO config. file.
    CALL HcoState_Init ( am_I_Root, HcoState, nnMatch, RC )
    IF(RC /= HCO_SUCCESS) RETURN 

    !-----------------------------------------------------------------
    ! Set grid
    CALL Set_Grid ( am_I_Root, RC )
    IF(RC /= HCO_SUCCESS) RETURN 

    !-----------------------------------------------------------------
    ! Register species
    CALL Register_Species( am_I_Root, RC )
    IF(RC /= HCO_SUCCESS) RETURN 

    !-----------------------------------------------------------------
    ! Read time information, incl. timesteps and simulation time(s) 
    CALL Read_Time( am_I_Root, RC )

    !=================================================================
    ! Set misc. parameter
    !=================================================================

    ! Set ESMF flag 
    HcoState%isESMF = .FALSE.  

    ! HEMCO configuration file
    HcoState%ConfigFile = ConfigFile

    ! Let HEMCO schedule the diagnostics output
    HcoState%Options%HcoWritesDiagn = .TRUE.

    !=================================================================
    ! Initialize HEMCO internal lists and variables. All data
    ! information is written into internal lists (ReadList) and 
    ! the HEMCO configuration file is removed from buffer in this
    ! step. Also initializes the HEMCO clock
    !=================================================================
    CALL HCO_Init( am_I_Root, HcoState, RC )
    IF(RC /= HCO_SUCCESS) RETURN 

    !=================================================================
    ! Initialize extensions.
    ! This initializes all (enabled) extensions and selects all met.
    ! fields needed by them. 
    !=================================================================
    CALL HCOX_Init( am_I_Root, HcoState, ExtState, RC )
    IF(RC /= HCO_SUCCESS) RETURN 

    !=================================================================
    ! Define diagnostics
    !=================================================================
    CALL Define_Diagnostics( am_I_Root, RC )
    IF( RC /= HCO_SUCCESS ) RETURN 

    !-----------------------------------------------------------------
    ! Leave 
    !-----------------------------------------------------------------

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
  SUBROUTINE HCOI_SA_RUN( am_I_Root, RC )
!
! !USES:
!
    USE HCO_FluxArr_Mod,       ONLY : HCO_FluxarrReset 
    USE HCO_Clock_Mod,         ONLY : HcoClock_Set
    USE HCO_Clock_Mod,         ONLY : HcoClock_Get
    USE HCO_Clock_Mod,         ONLY : HcoClock_Increase
    USE HCO_Driver_Mod,        ONLY : HCO_RUN
    USE HCOX_Driver_Mod,       ONLY : HCOX_RUN
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN   ) :: am_I_Root  ! root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT) :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: CNT
    INTEGER            :: YR, MT, DY, HR, MN, SC 
    CHARACTER(LEN=255) :: MSG, LOC

    !=================================================================
    ! HCOI_SA_RUN begins here!
    !=================================================================

    ! For error handling
    LOC = 'HCOI_SA_Run (hcoi_standalone_mod.F90)'

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

       !=================================================================
       ! Set HcoClock. On first call, use specified start date. Increase
       ! clock by one emission time step otherwise.
       !=================================================================
       IF ( CNT == 1 ) THEN
          CALL HcoClock_Set ( am_I_Root, HcoState, YRS(1), MTS(1), &
                              DYS(1),    HRS(1),   MNS(1), SCS(1), &
                              IsEmisTime=.TRUE.,   RC=RC)
          IF ( RC /= HCO_SUCCESS) RETURN 
       ELSE   
          CALL HcoClock_Increase ( am_I_Root,        HcoState,    &
                                   HcoState%TS_EMIS, .TRUE., RC=RC )
          IF ( RC /= HCO_SUCCESS) RETURN 
       ENDIF

       ! Get current time
       CALL HcoClock_Get ( cYYYY=YR, cMM=MT, cDD=DY, cH=HR, cM=MN, cS=SC, RC=RC )
       IF ( RC /= HCO_SUCCESS) RETURN 

       ! Leave loop if this is the end of the simulation
       IF ( IsEndOfSimulation(YR,MT,DY,HR,MN,SC) ) EXIT

       ! Write to logfile and standard output
       WRITE(MSG,100) YR, MT, DY, HR, MN, SC 
100    FORMAT( 'Calculate emissions at ', i4,'-',i2.2,'-',i2.2,' ', &
                 i2.2,':',i2.2,':',i2.2 )
       CALL HCO_MSG(MSG)
       WRITE(*,*) TRIM(MSG)
 
       ! ================================================================
       ! Reset all emission and deposition values
       ! ================================================================
       CALL HCO_FluxArrReset( HcoState, RC )
       IF ( RC/= HCO_SUCCESS) RETURN 
    
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
       ! ================================================================
       CALL HCO_Run( am_I_Root, HcoState, -1, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
   
       ! ================================================================
       ! Run HCO extensions
       ! ================================================================
   
       ! Set / update ExtState fields 
       CALL ExtState_SetFields ( am_I_Root, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       CALL ExtState_UpdateFields( am_I_Root, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
 
       ! Execute all enabled emission extensions. Emissions will be 
       ! added to corresponding flux arrays in HcoState.
       CALL HCOX_Run ( am_I_Root, HcoState, ExtState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       !=================================================================
       ! Update all autofill diagnostics 
       !=================================================================
       CALL HcoDiagn_AutoUpdate ( am_I_Root, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

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
  SUBROUTINE HCOI_SA_Final( am_I_Root )
!
! !USES:
!
    USE HCO_Driver_Mod,    ONLY : HCO_Final
    USE HCOX_Driver_Mod,   ONLY : HCOX_Final
    USE HCO_State_Mod,     ONLY : HcoState_Final
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN) :: am_I_Root
!
! !REVISION HISTORY: 
!  12 Sep 2013 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, RC
    CHARACTER(LEN=255) :: LOC
    CHARACTER(LEN= 31) :: RST='HEMCO_restart'

    !=================================================================
    ! HCOI_SA_FINAL begins here!
    !=================================================================

    ! Init
    LOC = 'HCOI_SA_FINAL (hco_standalone_mod.F90)'
 
    ! Cleanup HCO core
    CALL HCO_FINAL( am_I_Root, HcoState, .FALSE., RC )
    IF (RC /= HCO_SUCCESS) RETURN 

    ! Cleanup extensions and ExtState object
    ! This will also nullify all pointer to the met fields. 
    CALL HCOX_FINAL( am_I_Root, HcoState, ExtState, RC ) 
    IF (RC /= HCO_SUCCESS) RETURN 

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
  SUBROUTINE Model_GetSpecies( am_I_Root,                          &
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
    LOGICAL,            INTENT(IN ) :: am_I_Root
    INTEGER,            INTENT(OUT) :: nModelSpec
    CHARACTER(LEN= 31), POINTER     :: ModelSpecNames(:)
    INTEGER,            POINTER     :: ModelSpecIDs  (:)
    REAL(hp),           POINTER     :: ModelSpecMW   (:)
    REAL(hp),           POINTER     :: ModelSpecEmMW (:)
    REAL(hp),           POINTER     :: ModelSpecMolecRatio(:)
    REAL(hp),           POINTER     :: ModelSpecK0   (:)
    REAL(hp),           POINTER     :: ModelSpecCR   (:)
    REAL(hp),           POINTER     :: ModelSpecPKA  (:)
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
    CALL GetExtOpt ( CoreNr, 'SpecFile', OptValChar=MySpecFile, &
                     Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( FOUND ) SpecFile = MySpecFile

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Open spec file 
    OPEN( IU_FILE, FILE=TRIM(SpecFile), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       MSG = 'Error 1 reading ' // TRIM(SpecFile)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Get number of species
    nModelSpec = 0 
    DO 
       CALL GetNextLine( am_I_Root, IU_FILE, DUM, EOF, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
       IF ( EOF ) EXIT
       nModelSpec = nModelSpec + 1
    ENDDO

    ! Make sure we have one species
    IF ( nModelSpec == 0 ) THEN
       MSG = 'Species file ' // TRIM(SpecFile)      // &
             ' does not seem to have any content. ' // &
             'You must define at least one species.'
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
    ENDIF

    ! Go back to line one
    REWIND( IU_FILE )

    ! Get next valid line
!    CALL GetNextLine( am_I_Root, IU_FILE, DUM, EOF, RC )
!    IF ( RC /= HCO_SUCCESS .OR. EOF ) THEN
!       MSG = 'Error 2 reading ' // TRIM(SpecFile)
!       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
!       RETURN
!    ENDIF
!
!    LNG = LEN(TRIM(DUM))
!    LOW = NextCharPos ( TRIM(DUM), HCO_COL(), 1 )
!    IF ( LOW < 0 .OR. LOW == LNG ) THEN
!       MSG = 'Cannot extract index after colon: ' // TRIM(DUM)
!       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
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

       CALL GetNextLine( am_I_Root, IU_FILE, DUM, EOF, RC )
       IF ( RC /= HCO_SUCCESS .OR. EOF ) THEN
          WRITE(MSG,100) N, TRIM(SpecFile)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
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
                CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
                RETURN
             ENDIF
             UPP = NextCharPos( TRIM(DUM), HCO_SPC(), LOW )
             IF ( UPP < 0 ) UPP = LNG
          ENDDO
   
          UPP = UPP - 1 ! Don't read space
   
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
  SUBROUTINE SET_Grid( am_I_Root, RC ) 
!
! !USES:
!
!    USE Grid_Mod,        ONLY : DoGridComputation
    USE inquireMod,      ONLY : findFreeLUN
    USE HCO_ExtList_Mod, ONLY : GetExtOpt, CoreNr
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN   ) :: am_I_Root
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT) :: RC
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller - Initial Version
!  11 May 2015 - C. Keller - Now provide lon/lat edges instead of assuming
!                            global grid. 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    INTEGER               :: NX, NY, NZ
    INTEGER               :: I, J, N, LNG, LOW, UPP
    INTEGER               :: SZ(3)
    INTEGER               :: IU_FILE, IOS
    REAL(hp)              :: RG(4)
    REAL(hp)              :: XMIN, XMAX
    REAL(hp)              :: YMIN, YMAX
    REAL(hp)              :: DVAL
    REAL(hp)              :: DLON, DLAT
    REAL(hp)              :: PI_180, YDGR, YSN, SIN_DELTA, AM2 
    LOGICAL               :: FOUND,   EOF
    CHARACTER(LEN=255)    :: MSG, LOC 
    CHARACTER(LEN=255)    :: MyGridFile 
    CHARACTER(LEN=2047)   :: DUM

    !=================================================================
    ! SET_GRID begins here
    !=================================================================

    ! For error handling
    LOC = 'SET_GRID (hco_standalone_mod.F90)'

    ! Set PI_180
    PI_180 = HcoState%Phys%PI / 180.0_hp

    ! Try to get GridFile from configuration file (in settings)
    CALL GetExtOpt ( CoreNr, 'GridFile', OptValChar=MyGridFile, &
                     Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( FOUND ) GridFile = MyGridFile

    ! ------------------------------------------------------------------
    ! Open grid file
    ! ------------------------------------------------------------------

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Open grid file 
    OPEN( IU_FILE, FILE=TRIM(GridFile), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       MSG = 'Error 1 reading ' // TRIM(GridFile)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
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
       CALL GetNextLine( am_I_Root, IU_FILE, DUM, EOF, RC )
       IF ( RC /= HCO_SUCCESS .OR. EOF ) THEN
          MSG = 'Error 2 reading ' // TRIM(GridFile)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Read integer after colon (this is the dimension size)
       LNG = LEN(TRIM(DUM))
       LOW = NextCharPos ( TRIM(DUM), HCO_COL(), 1 )
       IF ( LOW < 0 .OR. LOW == LNG ) THEN
          MSG = 'Cannot extract size information from ' // TRIM(DUM)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
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
       WRITE(MSG,*) 'Lower lon must be smaller than upper lon: ', XMIN, XMAX
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    IF ( YMIN >= YMAX ) THEN
       WRITE(MSG,*) 'Lower lat must be smaller than upper lat: ', YMIN, YMAX
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    IF ( XMIN < -180.0_hp ) THEN
       WRITE(MSG,*) 'Lower longitude must be between -180 and 180 degE: ', XMIN
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF   
    IF ( XMAX > 180.0_hp ) THEN
       WRITE(MSG,*) 'Upper longitude must be between -180 and 180 degE: ', XMAX
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF   
    IF ( YMIN < -90.0_hp ) THEN
       WRITE(MSG,*) 'Lower latitude must be between -90 and 90 degN: ', YMIN
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF   
    IF ( YMAX > 90.0_hp ) THEN
       WRITE(MSG,*) 'Upper latitude must be between -90 and 90 degN: ', YMAX
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
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
       CALL GetNextLine( am_I_Root, IU_FILE, DUM, EOF, RC )
       IF ( RC /= HCO_SUCCESS .OR. EOF ) THEN
          MSG = 'Error 3 reading ' // TRIM(GridFile)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
 
       ! Read integer after colon (this is the dimension size)
       LNG = LEN(TRIM(DUM))
       LOW = NextCharPos ( TRIM(DUM), HCO_COL(), 1 )
       IF ( LOW < 0 .OR. LOW == LNG ) THEN
          MSG = 'Cannot extract size information from ' // TRIM(DUM)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       LOW = LOW + 1
       READ( DUM(LOW:LNG), * ) SZ(N)

    ENDDO !N

    ! Grid dimensions
    NX = SZ(1)
    NY = SZ(2)
    NZ = SZ(3)

    ! Close file
    CLOSE( IU_FILE )      

    ! ------------------------------------------------------------------
    ! Now that sizes are known, allocate all arrays
    ! ------------------------------------------------------------------
    ALLOCATE ( XMID     (NX,  NY,  1) )
    ALLOCATE ( YMID     (NX,  NY,  1) )
    ALLOCATE ( XEDGE    (NX+1,NY,  1) )
    ALLOCATE ( YEDGE    (NX,  NY+1,1) )
    ALLOCATE ( YSIN     (NX,  NY+1,1) )
    ALLOCATE ( AREA_M2  (NX,  NY,  1) )
    XMID      = 0.0_hp
    YMID      = 0.0_hp
    XEDGE     = 0.0_hp
    YEDGE     = 0.0_hp
    YSIN      = 0.0_hp
    AREA_M2   = 0.0_hp

    ! ------------------------------------------------------------------
    ! Fill grid box values
    ! ------------------------------------------------------------------
    DLON = ( XMAX - XMIN ) / NX
    DLAT = ( YMAX - YMIN ) / NY

    ! Now fill values
    DO J = 1, NY
    DO I = 1, NX

       ! Set longitude and latitude edge values
       XEDGE(I,J,1) = XMIN + ( (I-1) * DLON )
       YEDGE(I,J,1) = YMIN + ( (J-1) * DLAT )

       ! Set mid values
       XMID(I,J,1) = XEDGE(I,J,1) + ( DLON / 2.0_hp )
       YMID(I,J,1) = YEDGE(I,J,1) + ( DLAT / 2.0_hp )

       ! Get sine of latitude edges
       YDGR        = PI_180 * YEDGE(I,J,1)  ! radians       
       YSN         = SIN( YDGR )            ! sine
       YSIN(I,J,1) = YSN

       ! Eventually set uppermost edge
       IF ( I == NX ) THEN
          XEDGE(I+1,J,1) = XMIN + I * DLON
       ENDIF
       IF ( J == NY ) THEN
          YEDGE(I,J+1,1) = YMIN + J * DLAT
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

    ! Write grid information to log-file
    WRITE(MSG,*) 'HEMCO grid definitions:'
    CALL HCO_MSG(MSG)

    WRITE(MSG,*) ' --> Number of longitude cells: ', NX
    CALL HCO_MSG(MSG)
    WRITE(MSG,*) ' --> Number of latitude cells : ', NY
    CALL HCO_MSG(MSG)
    WRITE(MSG,*) ' --> Number of levels         : ', NZ
    CALL HCO_MSG(MSG)
    WRITE(MSG,*) ' --> Lon min, max [deg E]     : ', XMIN, XMAX
    CALL HCO_MSG(MSG)
    WRITE(MSG,*) ' --> Lon delta    [deg E]     : ', DLON
    CALL HCO_MSG(MSG)
    WRITE(MSG,*) ' --> Lat min, max [deg N]     : ', YMIN, YMAX
    CALL HCO_MSG(MSG)
    WRITE(MSG,*) ' --> Lat delta    [deg N]     : ', DLAT
    CALL HCO_MSG(MSG,SEP2="-")

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
  SUBROUTINE Get_nnMatch( am_I_Root, nnMatch, RC ) 
!
! !USES:
!
    USE HCO_Config_Mod, ONLY : Config_GetnSpecies
    USE HCO_Config_Mod, ONLY : Config_GetSpecNames
!
! !OUTPUT PARAMETERS:
!
    LOGICAL, INTENT(IN   )  :: am_I_Root ! Root CPU? 
    INTEGER, INTENT(  OUT)  :: nnMatch   ! Number of HEMCO species that are
                                         ! also species in the atm. model
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT)  :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER            :: AS, IDX
    CHARACTER(LEN=255) :: LOC

    !=================================================================
    ! Get_nnMatch begins here
    !=================================================================

    ! For error handling
    LOC = 'Get_nnMatch (hco_standalone_mod.F90)'
    
    ! Extract number of HEMCO species and corresponding species names 
    ! as read from the HEMCO config. file.
    nHcoSpec = Config_GetnSpecies ( ) 
    CALL Config_GetSpecNames( HcoSpecNames, nHcoSpec, RC )
    IF ( RC /= HCO_SUCCESS) RETURN 

    ! Extract species to be used from input file
    CALL Model_GetSpecies( am_I_Root,                           &
                           nModelSpec,     ModelSpecNames,      &
                           ModelSpecIDs,   ModelSpecMW,         &
                           ModelSpecEmMW,  ModelSpecMolecRatio, &
                           ModelSpecK0,    ModelSpecCR,         &
                           ModelSpecPKA,   RC                    )
    IF ( RC /= HCO_SUCCESS) RETURN

    ! See how many species are also used in GEOS-Chem
    ALLOCATE(matchIDx(nHcoSpec),STAT=AS)
    IF ( AS/=0 ) THEN 
       CALL HCO_ERROR('Allocation error matchIDx', RC, THISLOC=LOC )
       RETURN
    ENDIF
    matchIDx(:) = -1
    CALL HCO_CharMatch( HcoSpecNames,   nHcoSpec,      &
                        ModelSpecNames, nModelSpec,    &
                        matchIDx,       nnMatch         )
    IF ( nnMatch == 0 ) THEN
       CALL HCO_ERROR('No matching species!', RC, THISLOC=LOC )
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
  SUBROUTINE Register_Species( am_I_Root, RC )
!
! !USES:
!
    USE HCO_LogFile_Mod, ONLY : HCO_SPEC2LOG
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN   ) :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT) :: RC          ! Success or failure
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
       cnt = cnt + 1

       ! Set species name and GEOS-Chem tracer ID 
       IDX = ModelSpecIDs(MatchIDx(I))
       HcoState%Spc(cnt)%SpcName  = HcoSpecNames(I) 
       HcoState%Spc(cnt)%ModID    = IDX

       ! Molecular weights of species & emitted species.
       HcoState%Spc(cnt)%MW_g   = ModelSpecMW(IDX)
       HcoState%Spc(cnt)%EmMW_g = ModelSpecEmMW(IDX)

       ! Emitted molecules per molecule of species.
       HcoState%Spc(cnt)%MolecRatio = ModelSpecMolecRatio(IDX)

       ! Set Henry coefficients
       HcoState%Spc(cnt)%HenryK0  = ModelSpecK0(IDX)
       HcoState%Spc(cnt)%HenryCR  = ModelSpecCR(IDX)
       HcoState%Spc(cnt)%HenryPKA = ModelSpecPKA(IDX)

       ! Logfile I/O
       CALL HCO_SPEC2LOG( am_I_Root, HcoState, Cnt )

    ENDDO !I
    CALL HCO_MSG(SEP1='-')

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
  SUBROUTINE Define_Diagnostics( am_I_Root, RC )
!
! !USES:
!
    USE HCO_EXTLIST_MOD,   ONLY : GetExtNr
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN   ) :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER, INTENT(INOUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER           :: I, N, ExtNr, HcoID
    CHARACTER(LEN=31) :: DiagnName

    !=================================================================
    ! DEFINE_DIAGNOSTICS begins here
    !=================================================================

    ! Get number of diagnostics currently defined in the default
    ! collection
    CALL DiagnCollection_Get ( HcoDiagnIDDefault, nnDiagn=N, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If there are no diagnostics defined yet, define some default
    ! diagnostics below. These are simply the overall emissions 
    ! (across all extensions, categories, hierarchies) for each
    ! HEMCO species. 
    IF ( N == 0 ) THEN

       ! Loop over all HEMCO species
       DO I = 1, HcoState%nSpc 
   
          ! Get HEMCO ID
          HcoID = HcoState%Spc(I)%HcoID
          IF ( HcoID <= 0 ) CYCLE
   
          ! Create diagnostics
          DiagnName = 'EMIS_' // TRIM(HcoState%Spc(I)%SpcName) 
          CALL Diagn_Create ( am_I_Root,                         &
                              HcoState  = HcoState,              &
                              cName     = DiagnName,             &
                              ExtNr     = -1,                    &
                              Cat       = -1,                    &
                              Hier      = -1,                    &
                              HcoID     = HcoID,                 &
                              SpaceDim  = 2,                     &
                              LevIDx    = -1,                    &
                              OutUnit   = 'kg/m2/s',             &
                              AutoFill  = 1,                     &
                              COL       = HcoDiagnIDDefault,     &
                              RC        = RC                      )
          IF ( RC /= HCO_SUCCESS ) RETURN

       ENDDO !I
    ENDIF

    !--------------------------------------------------------------------------
    ! Define some additional diagnostics
    !--------------------------------------------------------------------------
    ExtNr = GetExtNr ( 'LightNOx' )
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
          CALL Diagn_Create( am_I_Root,                     & 
                             HcoState  = HcoState,          &
                             cName     = TRIM( DiagnName ), &
                             cID       = N,                 &
                             ExtNr     = ExtNr,             &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = -1,                &
                             SpaceDim  = 2,                 &
                             LevIDx    = -1,                &
                             OutUnit   = 'flashes/min/km2', &
                             OutOper   = 'Mean',            &
                             COL       = HcoDiagnIDDefault, &
                             AutoFill  = 0,                 &
                             RC        = RC                  ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDDO
 
       ! ---------------------------------------------------------- 
       ! Diagnostics for convective cloud top height.
       ! ---------------------------------------------------------- 

       ! Define diagnostics name and ID
       DiagnName = 'LIGHTNING_CLOUD_TOP'
       N         = 56004

       ! Create diagnostic container
       CALL Diagn_Create( am_I_Root,                     &
                          HcoState  = HcoState,          &
                          cName     = TRIM( DiagnName ), &
                          cID       = N,                 &
                          ExtNr     = ExtNr,             &
                          Cat       = -1,                &
                          Hier      = -1,                &
                          HcoID     = -1,                &
                          SpaceDim  = 2,                 &
                          LevIDx    = -1,                &
                          OutUnit   = '1',               &
                          OutOper   = 'Mean',            &
                          COL       = HcoDiagnIDDefault, &
                          AutoFill  = 0,                 &
                          RC        = RC                  )
       IF ( RC /= HCO_SUCCESS ) RETURN

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
  SUBROUTINE Read_Time( am_I_Root, RC ) 
!
! !USES:
!
    USE inquireMod,      ONLY : findfreeLUN
    USE HCO_Extlist_Mod, ONLY : GetExtOpt, CoreNr
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN   ) :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS
!
    INTEGER, INTENT(INOUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  13 Sep 2013 - C. Keller - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER             :: AS, IOS, IU_FILE
    INTEGER             :: I,  N,   LNG, LOW
    LOGICAL             :: EOF, FOUND
    CHARACTER(LEN=255)  :: MSG, LOC, DUM
    CHARACTER(LEN=255)  :: MyTimeFile 

    !=================================================================
    ! READ_TIME begins here
    !=================================================================

    ! For error handling
    LOC = 'READ_TIME (hcoi_standalone_mod.F90)'

    ! Try to get TimeFile from configuration file (in settings)
    CALL GetExtOpt ( CoreNr, 'TimeFile', OptValChar=MyTimeFile, &
                     Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( FOUND ) TimeFile = MyTimeFile

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Open time file 
    OPEN( IU_FILE, FILE=TRIM(TimeFile), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       MSG = 'Error 1 reading ' // TRIM(TimeFile)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
 
    ! Read start and end of simulation
    DO N = 1,2

       CALL GetNextLine( am_I_Root, IU_FILE, DUM, EOF, RC )
       IF ( RC /= HCO_SUCCESS .OR. EOF ) THEN
          MSG = 'Error reading time in ' // TRIM(TimeFile)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
     
       ! Remove 'BEGIN: ' or 'END: ' at the beginning 
       LNG = LEN(TRIM(DUM))
       LOW = NextCharPos ( TRIM(DUM), HCO_COL(), 1 )
       IF ( LOW < 0 .OR. LOW == LNG ) THEN
          MSG = 'Cannot extract index after colon: ' // TRIM(DUM)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       LOW = LOW + 1
       DUM = ADJUSTL(DUM(LOW:LNG))
       LNG = LEN(TRIM(DUM))
      
       ! Times have to be stored as:
       ! YYYY-MM-DD HH:MM:SS
       ! --> read year from position 1:4, month from 6:7, etc.
       IF ( LNG /= 19 ) THEN
          MSG = 'Provided time stamp is not `YYYY-MM-DD HH:MM:SS`! ' // TRIM(DUM)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
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
    CALL GetNextLine( am_I_Root, IU_FILE, DUM, EOF, RC )
    IF ( (RC /= HCO_SUCCESS) .OR. EOF ) THEN
       MSG = 'Cannot read emission time step from ' // TRIM(TimeFile)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Get index after colon 
    LNG = LEN(TRIM(DUM))
    LOW = NextCharPos ( TRIM(DUM), HCO_COL(), 1 )
    IF ( LOW < 0 .OR. LOW == LNG ) THEN
       MSG = 'Cannot extract index after colon: ' // TRIM(DUM)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
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
  SUBROUTINE ExtState_SetFields ( am_I_Root, RC )
!
! !USES:
!
    USE HCO_ARR_MOD,        ONLY : HCO_ArrAssert
    USE HCO_EMISLIST_MOD,   ONLY : HCO_GetPtr
    USE HCO_GEOTOOLS_MOD,   ONLY : HCO_GetSUNCOS
    USE HCOX_STATE_MOD,     ONLY : ExtDat_Set
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN   ) :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS
!
    INTEGER, INTENT(INOUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  28 Jul 2014 - C. Keller   - Initial Version
!  06 Oct 2014 - M. Sulprizio- Remove PCENTER. Now calculate from pressure edges
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER                       :: I,  J,  L 
    INTEGER                       :: NX, NY, NZ
    REAL(sp), POINTER             :: Ptr3D(:,:,:) => NULL()
    REAL(sp), POINTER             :: Ptr2D(:,:  ) => NULL()
    REAL(hp)                      :: P1, P2
    CHARACTER(LEN=255)            :: MSG
    LOGICAL                       :: ERR, FOUND
    LOGICAL, SAVE                 :: FIRST = .TRUE.
    CHARACTER(LEN=255), PARAMETER :: LOC = 'ExtState_SetFields (hcoi_standalone_mod.F90)'

    !=================================================================
    ! ExtState_SetFields begins here
    !=================================================================

    ! Enter
    CALL HCO_ENTER( LOC, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! Try to get pressure edges from external file.
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( am_I_Root, 'PEDGE', Ptr3D, RC, FOUND=FOUND )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If found, copy to pressure edge array
    IF ( FOUND ) THEN

       ! On first call, check array size dimensions. These are 
       ! expected to correspond to the HEMCO grid (NX, NY, NZ+1).
       IF ( FIRST ) THEN
          NX = SIZE(Ptr3D,1)
          NY = SIZE(Ptr3D,2)
          NZ = SIZE(Ptr3D,3)

          ! Make sure size dimensions are correct
          IF ( ( NX /=  HcoState%NX    ) .OR. &
               ( NY /=  HcoState%NY    ) .OR. &
               ( NZ /= (HcoState%NZ+1) )       ) THEN
             WRITE(MSG,*) 'Dimensions of field PEDGE do not correspond ', &
                          'to simulation grid: Expected dimensions: ',       &
                          HcoState%NX, HcoState%NY, HcoState%NZ+1,           &
                          '; dimensions of PEDGE: ', NX, NY, NZ
             CALL HCO_ERROR( MSG, RC )
             RETURN     
          ENDIF

          CALL HCO_ArrAssert( HcoState%Grid%PEDGE, NX, NY, NZ, RC ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF 
       HcoState%Grid%PEDGE%Val = Ptr3D
       Ptr3D => NULL()

    ! If not found, prompt a warning on the first call.
    ELSE
       IF ( FIRST ) THEN
          MSG = 'PEDGE is not a field in HEMCO configuration file '      // &
                '- cannot fill pressure edges. This may cause problems ' // &
                'in vertical interpolation and/or some of the extensions!'
          CALL HCO_WARNING ( MSG, RC, WarnLev=1 )
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! 2D fields 
    !-----------------------------------------------------------------

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%U10M, 'U10M', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%V10M, 'V10M', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%ALBD, 'ALBD', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%WLI, 'WLI', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%T2M, 'T2M', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%TSKIN, 'TSKIN', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%GWETROOT, 'GWETROOT', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%GWETTOP, 'GWETTOP', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%SNOWHGT, 'SNOWHGT', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%SNODP, 'SNODP', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%USTAR, 'USTAR', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%Z0, 'Z0', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%TROPP, 'TROPP', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%SZAFACT, 'SZAFACT', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%PARDR, 'PARDR', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%PARDF, 'PARDF', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%RADSWG, 'RADSWG', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%FRCLND, 'FRCLND', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%CLDFRC, 'CLDFRC', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%LAI, 'LAI', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%JNO2, 'JNO2', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%JO1D, 'JO1D', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%FRLAND, 'FRLAND', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%FROCEAN, 'FROCEAN', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%FRLAKE, 'FRLAKE', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%FRLANDIC, 'FRLANDIC', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! 3D fields 
    !-----------------------------------------------------------------

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%CNV_MFC, 'CNV_MFC', &
                      RC, FIRST,  OnLevEdge=.TRUE. )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%SPHU, 'SPHU', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%TK, 'TK', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%AIR, 'AIR', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%AIRVOL, 'AIRVOL', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%O3, 'O3', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%NO, 'NO', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%NO2, 'NO2', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%HNO3, 'HNO3', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%DRY_TOTN, 'DRY_TOTN', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%WET_TOTN, 'WET_TOTN', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL ExtDat_Set ( am_I_Root, HcoState, ExtState%FRAC_OF_PBL, 'FRAC_OF_PBL', RC, FIRST )
    IF ( RC /= HCO_SUCCESS ) RETURN

    !-----------------------------------------------------------------
    ! ==> DRYCOEFF must be read from the configuration file in module
    !     hcox_soilnox_mod.F90. 
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! If needed, calculate SUNCOS values
    !-----------------------------------------------------------------
    IF ( ExtState%SUNCOS%DoUse ) THEN
       IF ( FIRST ) THEN
          CALL HCO_ArrAssert( ExtState%SUNCOS%Arr, HcoState%NX, HcoState%NY, RC )
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF

       CALL HCO_GetSUNCOS( am_I_Root, HcoState, ExtState%SUNCOS%Arr%Val, 0, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN
    ENDIF

    !-----------------------------------------------------------------
    ! Grid box heights in meters. First try to get them from an 
    ! external file. If not defined in there, try to calculate them
    ! from grid box pressure edges and temperature via the hydrostatic
    ! equation.
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( am_I_Root, 'BXHEIGHT_M', Ptr3D, RC, FOUND=FOUND )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If found, copy to pressure edge array
    IF ( FOUND ) THEN

       ! On first call, check array size dimensions. These are 
       ! expected to correspond to the HEMCO grid (NX, NY, NZ).
       IF ( FIRST ) THEN
          NX = SIZE(Ptr3D,1)
          NY = SIZE(Ptr3D,2)
          NZ = SIZE(Ptr3D,3)

          ! Make sure size dimensions are correct
          IF ( ( NX /=  HcoState%NX  ) .OR. &
               ( NY /=  HcoState%NY  ) .OR. &
               ( NZ /= (HcoState%NZ) )       ) THEN
             WRITE(MSG,*) 'Dimensions of field BXHEIGHT_M do not correspond ', &
                          'to simulation grid: Expected dimensions: ',         &
                          HcoState%NX, HcoState%NY, HcoState%NZ,               &
                          '; dimensions of BXHEIGHT_M: ', NX, NY, NZ
             CALL HCO_ERROR( MSG, RC )
             RETURN     
          ENDIF

          CALL HCO_ArrAssert( HcoState%Grid%BXHEIGHT_M, NX, NY, NZ, RC ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF 
       HcoState%Grid%BXHEIGHT_M%Val = Ptr3D
       Ptr3D => NULL()

       ! Verbose
       IF ( FIRST .AND. Hco_IsVerb(1) ) THEN
          MSG = 'Grid box heights read from field BXHEIGHT_M'
          CALL HCO_MSG(MSG,SEP1=' ',SEP2=' ')
       ENDIF

    ! If not found, check if we can calculate it from pressure edges and temperature
    ELSEIF ( ExtState%TK%DoUse .AND. ASSOCIATED(HcoState%Grid%PEDGE%Val) ) THEN

       ! Make sure array is defined
       CALL HCO_ArrAssert( HcoState%Grid%BXHEIGHT_M, HcoState%NX,   &
                           HcoState%NY,              HcoState%NZ, RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Assume no error until otherwise
       ERR = .FALSE.

       ! Calculate box heights
!$OMP PARALLEL DO                                                      &
!$OMP DEFAULT( SHARED )                                                &
!$OMP PRIVATE( I, J, L, P1, P2 )                                       & 
!$OMP SCHEDULE( DYNAMIC )
       DO L = 1, HcoState%NZ
       DO J = 1, HcoState%NY
       DO I = 1, HcoState%NX
       
          ! Pressure at bottom and top edge [hPa]
          P1 = HcoState%Grid%PEDGE%Val(I,J,L)
          P2 = HcoState%Grid%PEDGE%Val(I,J,L+1)

          ! Box height
          IF ( P2 == 0.0_hp ) THEN
             ERR = .TRUE.
          ELSE
             HcoState%Grid%BXHEIGHT_M%Val(I,J,L) = HcoState%Phys%Rdg0 &
                                                 * ExtState%TK%Arr%Val(I,J,L) &
                                                 * LOG( P1 / P2 )
          ENDIF
       ENDDO !I
       ENDDO !J
       ENDDO !L
!$OMP END PARALLEL DO

       ! Error check
       IF ( ERR ) THEN
          MSG = 'Cannot calculate grid box heights - at least one pressure edge ' // &
                'value is zero! You can either provide an updated pressure edge ' // &
                'field (PEDGE) or add a field with the grid box heigths to your ' // &
                'configuration file (BXHEIGHT_M)'
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Verbose
       IF ( FIRST .AND. Hco_IsVerb(1) ) THEN
          MSG = 'Grid box heights approximated from hydrostatic equation'
          CALL HCO_MSG(MSG,SEP1=' ',SEP2=' ')
       ENDIF

    ! Otherwise, prompt a warning on the first call 
    ELSE
       IF ( FIRST ) THEN
          MSG = 'BXHEIGHT_M is not a field in HEMCO configuration file '   // &
                '- cannot fill grid box heights. This may cause problems ' // &
                'in some of the extensions (if so, they will crash)!'
          CALL HCO_WARNING ( MSG, RC, WarnLev=1 )
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! Surface geopotential height 
    !-----------------------------------------------------------------
    CALL HCO_GetPtr( am_I_Root, 'ZSFC', Ptr2D, RC, FOUND=FOUND )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! If found, copy to pressure edge array
    IF ( FOUND ) THEN

       ! On first call, check array size dimensions. These are 
       ! expected to correspond to the HEMCO grid (NX, NY, NZ).
       IF ( FIRST ) THEN
          NX = SIZE(Ptr2D,1)
          NY = SIZE(Ptr2D,2)

          ! Make sure size dimensions are correct
          IF ( ( NX /=  HcoState%NX  ) .OR. &
               ( NY /=  HcoState%NY  )       ) THEN
             WRITE(MSG,*) 'Dimensions of field ZSFC do not correspond ', &
                          'to simulation grid: Expected dimensions: ',         &
                          HcoState%NX, HcoState%NY,                            &
                          '; dimensions of ZSFC: ', NX, NY
             CALL HCO_ERROR( MSG, RC )
             RETURN     
          ENDIF

          CALL HCO_ArrAssert( HcoState%Grid%ZSFC, NX, NY, RC ) 
          IF ( RC /= HCO_SUCCESS ) RETURN
       ENDIF 
       HcoState%Grid%ZSFC%Val = Ptr2D
       Ptr2D => NULL()

       ! Verbose
       IF ( FIRST .AND. Hco_IsVerb(1) ) THEN
          MSG = 'Surface geopotential heights read from field ZSFC'
          CALL HCO_MSG(MSG,SEP1=' ',SEP2=' ')
       ENDIF

    ! If not found, check if we can calculate it from pressure edges and temperature
    ELSEIF ( ExtState%TK%DoUse .AND. ASSOCIATED(HcoState%Grid%PEDGE%Val) ) THEN

       ! Make sure array is defined
       CALL HCO_ArrAssert( HcoState%Grid%ZSFC, HcoState%NX, HcoState%NY, RC ) 
       IF ( RC /= HCO_SUCCESS ) RETURN

       ! Assume no error until otherwise
       ERR = .FALSE.

       ! Approximate geopotential heights
!$OMP PARALLEL DO                                                      &
!$OMP DEFAULT( SHARED )                                                &
!$OMP PRIVATE( I, J, P1, P2 )                                       & 
!$OMP SCHEDULE( DYNAMIC )
       DO J = 1, HcoState%NY
       DO I = 1, HcoState%NX
       
          ! Pressure at bottom and top edge [Pa]
          P1 = 101325.0
          P2 = HcoState%Grid%PEDGE%Val(I,J,1)

          ! Box height
          IF ( P2 == 0.0_hp ) THEN
             ERR = .TRUE.
          ELSE
             HcoState%Grid%ZSFC%Val(I,J) = HcoState%Phys%Rdg0 &
                                         * ExtState%TK%Arr%Val(I,J,1) &
                                         * LOG( P1 / P2 )
          ENDIF
       ENDDO !I
       ENDDO !J
!$OMP END PARALLEL DO

       ! Error check
       IF ( ERR ) THEN
          MSG = 'Cannot calculate surface geopotential heights - at least one ' // &
                'surface pressure value is zero! You can either provide an '    // &
                'updated pressure edge field (PEDGE) or add a field with the '  // &
                'surface geopotential height to your configuration file (ZSFC)'
          CALL HCO_ERROR ( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       ! Verbose
       IF ( FIRST .AND. Hco_IsVerb(1) ) THEN
          MSG = 'Geopotential heights approximated from hydrostatic equation'
          CALL HCO_MSG(MSG,SEP1=' ',SEP2=' ')
       ENDIF

    ! Otherwise, prompt a warning on the first call 
    ELSE
       IF ( FIRST ) THEN
          MSG = 'ZSFC is not a field in HEMCO configuration file '   // &
                '- cannot define surface geopotential height. This may cause ' // &
                'problems in some of the extensions (if so, they will crash)!'
          CALL HCO_WARNING ( MSG, RC, WarnLev=1 )
       ENDIF
    ENDIF

    !-----------------------------------------------------------------
    ! All done 
    !-----------------------------------------------------------------

    ! Not first call any more
    FIRST = .FALSE.

    ! Leave w/ success
    CALL HCO_LEAVE( RC )

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
  SUBROUTINE ExtState_UpdateFields ( am_I_Root, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS
!
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    INTEGER       :: THISDATETIME
    INTEGER, SAVE :: ENDDATETIME = -1

    !=================================================================
    ! IsEndOfSimulation begins here
    !=================================================================

    ! Init
    IsEnd = .FALSE.

    ! Calculate simulation end datetime if not yet done so
    IF ( ENDDATETIME < 0 ) THEN
       ENDDATETIME = YRS(2)*10000000000 + MTS(2)*100000000 + DYS(2)*1000000 &
                   + HRS(2)*10000       + MNS(2)*100       + SCS(2)
    ENDIF

    ! Calculate current datetime
    THISDATETIME = YR*10000000000 + MT*100000000 + DY*1000000 &
                 + HR*10000       + MN*100       + SC

    ! Check if current datetime is beyond simulation end date
    IF ( THISDATETIME >= ENDDATETIME ) IsEnd = .TRUE.

  END FUNCTION IsEndOfSimulation
!EOC
END MODULE HCOI_StandAlone_Mod
