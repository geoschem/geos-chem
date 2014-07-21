!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_standalone_mod.F90 
!
! !DESCRIPTION: Module HCO\_StandAlone\_Mod contains all wrapper routines 
! to run HEMCO in standalone mode, i.e. without any external model connected
! to it. All HEMCO input variables (grid, species, times) are read from disk.
! Note that for now, it is not possible to use HEMCO extensions that rely on
! meteorological variables when running HEMCO in standalone mode
!\\
!\\
! Subroutine HCO\_StandAlone\_Run will execute the standalone version of 
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
! \item HEMCO\_sa\_Time: contains the time definitions. The first row must
!  contain the emission time step (e.g. TS\_EMIS: 3600.0). The second row must
!  contain the number of desired time steps (e.g. NSTEPS: 1), and the following
!  rows (3-END) contain the dates and times of all time steps in format 
!  YYYY-MM-DD HH:MM:SS (e.g. 2013-07-01 00:00:00).
! 
! \end{itemize}
!
! The file names of the species, grid, and time input files can be provided
! in the settings section of the HEMCO configuration file 
! (e.g. 'SpecFile: myHEMCO\_sa\_Spec'). Otherwise, the default file names
! (HEMCO\_sa\_Spec, HEMCO\_sa\_Grid, HEMCO\_sa\_Time) will be used.
!
! !INTERFACE:
!
MODULE HCO_StandAlone_Mod
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
  PUBLIC  :: HCO_StandAlone_Run
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
  PRIVATE :: Read_Time
!
! !REVISION HISTORY:
!  20 Aug 2013 - C. Keller   - Initial version. 
!  14 Jul 2014 - R. Yantosca - Now use F90 free-format indentation
!  14 Jul 2014 - R. Yantosca - Cosmetic changes in ProTeX headers
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

  ! Grid
  REAL(df), ALLOCATABLE, TARGET  :: XMID   (:,:,:)
  REAL(df), ALLOCATABLE, TARGET  :: YMID   (:,:,:)
  REAL(df), ALLOCATABLE, TARGET  :: XEDGE  (:,:,:)
  REAL(df), ALLOCATABLE, TARGET  :: YEDGE  (:,:,:)
  REAL(df), ALLOCATABLE, TARGET  :: YSIN   (:,:,:)
  REAL(df), ALLOCATABLE, TARGET  :: AREA_M2(:,:,:)

  ! Simulation times
  INTEGER                        :: NSTEPS
  INTEGER,  ALLOCATABLE          :: YRS(:)
  INTEGER,  ALLOCATABLE          :: MTS(:)
  INTEGER,  ALLOCATABLE          :: DYS(:)
  INTEGER,  ALLOCATABLE          :: HRS(:)
  INTEGER,  ALLOCATABLE          :: MNS(:)
  INTEGER,  ALLOCATABLE          :: SCS(:)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HCO_StandAlone_Run
!
! !DESCRIPTION: Subroutine HCO\_StandAlone\_Run runs the standalone version
! of HEMCO. All input variables are taken from input files.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HCO_StandAlone_Run( ConfigFile )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=255), INTENT(IN)  :: ConfigFile
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
    ! HCO_STANDALONE_RUN begins here!
    !=================================================================

    ! Treat as root CPU
    am_I_Root = .TRUE.

    ! Initialize the HEMCO standalone
    CALL HCOI_Sa_Init( am_I_Root, ConfigFile, RC )
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
    
  END SUBROUTINE HCO_StandAlone_Run
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
    USE HCOI_GC_Diagn_Mod, ONLY : HCOI_Diagn_Init
!
! !INPUT PARAMETERS:
!
    LOGICAL,            INTENT(IN   ) :: am_I_Root  ! root CPU?
    CHARACTER(LEN=255), INTENT(IN   ) :: ConfigFile ! Configuration file
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
    CALL Config_ReadFile( am_I_Root, ConfigFile, RC )
    IF( RC /= HCO_SUCCESS) RETURN 

    !=================================================================
    ! Open logfile 
    !=================================================================
    IF ( am_I_Root ) THEN
       CALL HCO_LogFile_Open( RC=RC ) 
       IF (RC /= HCO_SUCCESS) RETURN 
    ELSE
       ! If this is not the root CPU, always disable verbose mode.
       CALL HCO_Verbose_Set( .FALSE. )
    ENDIF

    !=================================================================
    ! Initialize HEMCO state object and populate it 
    !=================================================================

    !-----------------------------------------------------------------
    ! Extract species to use in HEMCO 
    CALL Get_nnMatch( nnMatch, RC )
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

    !-----------------------------------------------------------------
    ! Set pointers to met fields.
    !-----------------------------------------------------------------
    ! tbd

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
    ! HEMCO routines 
    USE HCO_FluxArr_Mod,       ONLY : HCO_FluxarrReset 
    USE HCO_Clock_Mod,         ONLY : HcoClock_Set
    USE HCO_Driver_Mod,        ONLY : HCO_RUN
    USE HCOX_Driver_Mod,       ONLY : HCOX_RUN
    USE HCOI_GC_Diagn_Mod,     ONLY : HCOI_DIAGN_WRITEOUT
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
    INTEGER            :: N
    CHARACTER(LEN=255) :: MSG, LOC

    !=================================================================
    ! HCOI_SA_RUN begins here!
    !=================================================================

    ! For error handling
    LOC = 'HCOI_SA_Run (hcoi_standalone_mod.F90)'

    ! Do for every time step
    DO N = 1, NSTEPS

       ! Write to output
       WRITE(MSG,100) YRS(N), MTS(N), DYS(N), HRS(N), MNS(N), SCS(N)
100    FORMAT( 'Calculate emissions at ', i4,'-',i2.2,'-',i2.2,' ', &
                 i2.2,':',i2.2,':',i2.2 )
       CALL HCO_Msg(MSG)

       !=================================================================
       ! Set HcoClock 
       !=================================================================
       CALL HcoClock_Set ( am_I_Root, HcoState, YRS(N), MTS(N), &
                           DYS(N),    HRS(N),   MNS(N), SCS(N), RC=RC)
       IF ( RC/= HCO_SUCCESS) RETURN 
   
       !=================================================================
       ! Output diagnostics 
       !=================================================================
       CALL HCOI_Diagn_WriteOut ( am_I_Root, HcoState, .FALSE., RC )
       IF ( RC/= HCO_SUCCESS) RETURN 
    
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
       CALL HCO_Run( am_I_Root, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN 
   
       ! ================================================================
       ! Run HCO extensions
       ! ================================================================
   
       ! Execute all enabled emission extensions. Emissions will be 
       ! added to corresponding flux arrays in HcoState.
       CALL HCOX_Run ( am_I_Root, HcoState, ExtState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

       !=================================================================
       ! Update all autofill diagnostics 
       !=================================================================
       CALL HCO_Diagn_Update ( am_I_Root, HcoState, RC )
       IF ( RC /= HCO_SUCCESS ) RETURN

    ENDDO !N

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
    USE HCOI_GC_Diagn_Mod, ONLY : HCOI_Diagn_WriteOut
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

    !=================================================================
    ! HCOI_SA_FINAL begins here!
    !=================================================================

    ! Init
    LOC = 'HCOI_SA_FINAL (hco_standalone_mod.F90)'

    ! Write out all diagnostics
    CALL HCOI_DIAGN_WRITEOUT ( am_I_Root, HcoState, .TRUE., RC, &
                               UsePrevTime=.FALSE. )
    IF (RC /= HCO_SUCCESS) RETURN 
 
    ! Cleanup diagnostics
    CALL Diagn_Cleanup()
 
    ! Cleanup extensions and ExtState object
    ! This will also nullify all pointer to the met fields. 
    CALL HCOX_FINAL( ExtState ) 

    ! Cleanup HCO core
    CALL HCO_FINAL()

    ! Remove pointer references to grid arrays first.
    HcoState%Grid%XMID       => NULL() 
    HcoState%Grid%YMID       => NULL() 
    HcoState%Grid%XEDGE      => NULL()
    HcoState%Grid%YEDGE      => NULL()
    HcoState%Grid%YSIN       => NULL()
    HcoState%Grid%AREA_M2    => NULL()

    ! Deallocate module arrays/pointers
    IF ( ALLOCATED( XMID    ) ) DEALLOCATE ( XMID    )
    IF ( ALLOCATED( YMID    ) ) DEALLOCATE ( YMID    )
    IF ( ALLOCATED( XEDGE   ) ) DEALLOCATE ( XEDGE   )
    IF ( ALLOCATED( YEDGE   ) ) DEALLOCATE ( YEDGE   )
    IF ( ALLOCATED( YSIN    ) ) DEALLOCATE ( YSIN    )
    IF ( ALLOCATED( AREA_M2 ) ) DEALLOCATE ( AREA_M2 )
    IF ( ALLOCATED( YRS     ) ) DEALLOCATE ( YRS     )
    IF ( ALLOCATED( MTS     ) ) DEALLOCATE ( MTS     )
    IF ( ALLOCATED( DYS     ) ) DEALLOCATE ( DYS     )
    IF ( ALLOCATED( HRS     ) ) DEALLOCATE ( HRS     )
    IF ( ALLOCATED( MNS     ) ) DEALLOCATE ( MNS     )
    IF ( ALLOCATED( SCS     ) ) DEALLOCATE ( SCS     )

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
  SUBROUTINE Model_GetSpecies( nModelSpec,     ModelSpecNames,     &
                               ModelSpecIDs,   ModelSpecMW,        &
                               ModelSpecEmMW,  ModelSpecMolecRatio,&
                               ModelSpecK0,    ModelSpecCR,        &
                               ModelSpecPKA,   RC                   )
!
! !USES:
!
    USE inquireMod,      ONLY : findfreeLUN
    USE HCO_EXTLIST_Mod, ONLY : GetExtOpt
!
! !OUTPUT PARAMETERS:
!
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
    LOGICAL             :: FOUND
    CHARACTER(LEN=255)  :: MSG, LOC 
    CHARACTER(LEN=255)  :: MySpecFile 
    CHARACTER(LEN=2047) :: DUM

    !=================================================================
    ! Model_GetSpecies begins here
    !=================================================================

    ! For error handling
    LOC = 'Model_GetSpecies (hcoi_standalone_mod.F90)'

    ! Try to get SpecFile from configuration file (in settings)
    CALL GetExtOpt ( 0, 'SpecFile', OptValChar=MySpecFile, &
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
    READ( IU_FILE, '(A)', IOSTAT=IOS ) DUM
    IF ( IOS /= 0 ) THEN
       MSG = 'Error 2 reading ' // TRIM(SpecFile)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    LNG = LEN(TRIM(DUM))
    LOW = NextCharPos ( TRIM(DUM), HCO_COL(), 1 )
    IF ( LOW < 0 .OR. LOW == LNG ) THEN
       MSG = 'Cannot extract index after colon: ' // TRIM(DUM)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    LOW = LOW + 1
    READ ( DUM(LOW:LNG), * ) nModelSpec
      
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

       ! Read line
       READ( IU_FILE, '(A)', IOSTAT=IOS ) DUM
       IF ( IOS /= 0 ) THEN
          WRITE(MSG,100) N, TRIM(SpecFile)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       LNG = LEN(TRIM(DUM))

       ! Start reading line from beginning 
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
          IF ( I == 1 ) THEN
             READ( DUM(LOW:UPP), * ) ModelSpecIDs(N)
          ELSEIF ( I == 2 ) THEN
             READ( DUM(LOW:UPP), * ) ModelSpecNames(N)
          ELSEIF ( I == 3 ) THEN
             READ( DUM(LOW:UPP), * ) ModelSpecMW(N)
          ELSEIF ( I == 4 ) THEN
             READ( DUM(LOW:UPP), * ) ModelSpecEmMW(N)
          ELSEIF ( I == 5 ) THEN
             READ( DUM(LOW:UPP), * ) ModelSpecMolecRatio(N)
          ELSEIF ( I == 6 ) THEN
             READ( DUM(LOW:UPP), * ) ModelSpecK0(N)
          ELSEIF ( I == 7 ) THEN
             READ( DUM(LOW:UPP), * ) ModelSpecCR(N)
          ELSEIF ( I == 8 ) THEN
             READ( DUM(LOW:UPP), * ) ModelSpecPKA(N)
          ENDIF

          ! Continue from upper position (+1 to skip space). The
          ! while loop at the beginning of the do-loop will advance
          ! low by another one position, so the next character
          ! search will start at position UPP + 2, which is exactly
          ! what we want (UPP is the position BEFORE the space!).
          LOW = UPP + 1

       ENDDO !I
    ENDDO !N

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
! !DESCRIPTION: SUBROUTINE SET\_GRID sets the grid when running HEMCO
!  in standalone mode. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_Grid( am_I_Root, RC ) 
!
! !USES:
!
    USE Grid_Mod,        ONLY : DoGridComputation
    USE inquireMod,      ONLY : findFreeLUN
    USE HCO_ExtList_Mod, ONLY : GetExtOpt
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    INTEGER               :: NX, NY, NZ, JSP, JNP
    INTEGER               :: I, N, M, LNG, LOW, UPP
    INTEGER               :: SZ(3)
    INTEGER               :: IU_FILE, IOS
    LOGICAL               :: FOUND
    CHARACTER(LEN=255)    :: MSG, LOC 
    CHARACTER(LEN=255)    :: MyGridFile 
    CHARACTER(LEN=2047)   :: DUM

    ! Arrays
    REAL(df)              :: DVAL
    REAL(df), ALLOCATABLE :: DLON(:,:,:)
    REAL(df), ALLOCATABLE :: DLAT(:,:,:)
    REAL(df), ALLOCATABLE :: YMID_R(:,:,:)
    REAL(df), ALLOCATABLE :: YEDGE_R(:,:,:)
    REAL(df), ALLOCATABLE :: YMID_R_W(:,:,:)
    REAL(df), ALLOCATABLE :: YEDGE_R_W(:,:,:)

    !=================================================================
    ! SET_GRID begins here
    !=================================================================

    ! For error handling
    LOC = 'SET_GRID (hco_standalone_mod.F90)'

    ! Try to get GridFile from configuration file (in settings)
    CALL GetExtOpt ( 0, 'GridFile', OptValChar=MyGridFile, &
                     Found=FOUND, RC=RC )
    IF ( RC /= HCO_SUCCESS ) RETURN
    IF ( FOUND ) GridFile = MyGridFile

    ! Read grid information from file

    ! Find a free file LUN
    IU_FILE = findFreeLUN()

    ! Open grid file 
    OPEN( IU_FILE, FILE=TRIM(GridFile), STATUS='OLD', IOSTAT=IOS )
    IF ( IOS /= 0 ) THEN
       MSG = 'Error 1 reading ' // TRIM(GridFile)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Get grid size (x,y,z)
    DO N = 1, 3
       READ( IU_FILE, '(A)', IOSTAT=IOS ) DUM
       IF ( IOS /= 0 ) THEN
          MSG = 'Error 2 reading ' // TRIM(GridFile)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       LNG = LEN(TRIM(DUM))
 
       ! Read integer after colon (this is the dimension size)
       LOW = NextCharPos ( TRIM(DUM), HCO_COL(), 1 )
       IF ( LOW < 0 .OR. LOW == LNG ) THEN
          MSG = 'Cannot extract size information from ' // TRIM(DUM)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       LOW = LOW + 1
       READ( DUM(LOW:LNG), * ) SZ(N) 
    ENDDO

    ! Grid dimensions
    NX = SZ(1)
    NY = SZ(2)
    NZ = SZ(3)

    ! Now that sizes are known, allocate all arrays
    ALLOCATE ( DLON(NX,NY,NZ) )
    DLON = 0.0_df
    ALLOCATE ( DLAT(NX,NY,NZ) )
    DLAT = 0.0_df
    ALLOCATE ( XMID(NX,NY,NZ) )
    XMID = 0.0_df
    ALLOCATE ( YMID(NX,NY,NZ) )
    YMID = 0.0_df
    ALLOCATE ( XEDGE(NX+1,NY,NZ) )
    XEDGE = 0.0_df
    ALLOCATE ( YEDGE(NX,NY+1,NZ) )
    YEDGE = 0.0_df
    ALLOCATE ( YSIN(NX,NY+1,NZ) )
    YSIN = 0.0_df
    ALLOCATE ( YMID_R(NX,NY,NZ) )
    YMID_R = 0.0_df
    ALLOCATE ( YEDGE_R(NX,NY+1,NZ) )
    YEDGE_R = 0.0_df
    ALLOCATE ( AREA_M2(NX,NY,NZ) )
    AREA_M2 = 0.0_df
    ALLOCATE ( YMID_R_W(NX,NY,NZ) )
    YMID_R_W = 0.0_df
    ALLOCATE ( YEDGE_R_W(NX,NY+1,NZ) )
    YEDGE_R_W = 0.0_df

    ! Read delta lon (dx) or lat (dy) from file. This can be
    ! specified for every grid point or as single value (applied
    ! to all grid boxes).
    DO N=1,2
       READ( IU_FILE, '(A)', IOSTAT=IOS ) DUM
       IF ( IOS /= 0 ) THEN
          MSG = 'Error 3 reading ' // TRIM(GridFile)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       LNG = LEN(TRIM(DUM))
          
       ! Get position after colon
       LOW = NextCharPos ( TRIM(DUM), HCO_COL(), 1 )
       IF ( LOW < 0 .OR. LOW == LNG ) THEN 
          MSG = 'Cannot extract grid space from ' // TRIM(DUM)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       LOW = LOW + 1
 
       ! Read data for each grid cell
       M = 1 ! Index in dlon/dlat
       DO
          ! Get index up to next space.
          UPP = NextCharPos ( TRIM(DUM), HCO_SPC(), LOW )

          ! If no space found in word after index LOW, read until end of word.
          IF ( UPP < 0 ) UPP = LNG

          ! Read value and pass to DLON / DLAT.
          ! Ignore if only space.
          IF ( DUM(LOW:UPP) /= HCO_SPC() ) THEN
             READ( DUM(LOW:UPP), * ) DVAL
             IF ( N == 1 ) THEN
                DLON(M,:,:) = DVAL
             ELSE
                DLAT(:,M,:) = DVAL
             ENDIF
             M = M+1 ! Advance index in DLON/DLAT
          ENDIF

          ! Continue at next position
          LOW = UPP + 1
          IF ( LOW >= LNG ) EXIT
       ENDDO

       ! Check if we have to fill up dlon
       IF ( M == 2 ) THEN
          IF ( N == 1 ) THEN
             DLON(:,:,:) = DVAL
          ELSE
             DLAT(:,:,:) = DVAL
          ENDIF
       ELSEIF ( M == 1 ) THEN
          MSG = 'no grid space(s) found: ' // TRIM(DUM)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ELSEIF ( N == 1 .AND. M /= (NX+1) ) THEN
          MSG = 'lon grid space error: ' // TRIM(DUM)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ELSEIF ( N == 2 .AND. M /= (NY+1) ) THEN
          MSG = 'lat grid space error: ' // TRIM(DUM)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
    ENDDO

    ! Close file
    CLOSE( IU_FILE )      

    ! Set grid dimensions
    HcoState%NX = NX
    HcoState%NY = NY
    HcoState%NZ = NZ 

    ! For now, assume global grid with south pole at index 1 and north 
    ! pole at last index position
    JSP = 1
    JNP = NY

    Call DoGridComputation( am_I_Root = am_I_Root, &
                            I1        = 1,         &
                            I2        = NX,        & 
                            J1        = 1,         &
                            J2        = NY,        & 
                            JSP       = JSP,       &
                            JNP       = JNP,       &
                            L1        = 1,         &
                            L2        = NZ,        &
                            DLON      = DLON,      &
                            DLAT      = DLAT,      &
                            I_LO      = 1,         &
                            J_LO      = 1,         &
                            IOFF      = 0,         &
                            JOFF      = 0,         &
                            XMD       = XMID,      &
                            XDG       = XEDGE,     &
                            YMD       = YMID,      &
                            YDG       = YEDGE,     &
                            YSN       = YSIN,      &
                            YMDR      = YMID_R,    &
                            YDGR      = YEDGE_R,   &
                            YMDRW     = YMID_R_W,  &
                            YDGRW     = YEDGE_R_W, &
                            AM2       = AREA_M2,   &
                            RC        = RC          )

    ! Set pointers to grid variables
    HcoState%Grid%XMID       => XMID   (:,:,1)
    HcoState%Grid%YMID       => YMID   (:,:,1)
    HcoState%Grid%XEDGE      => XEDGE  (:,:,1)
    HcoState%Grid%YEDGE      => YEDGE  (:,:,1)
    HcoState%Grid%YSIN       => YSIN   (:,:,1)
    HcoState%Grid%AREA_M2    => AREA_M2(:,:,1)

    ! Cleanup
    DEALLOCATE( YMID_R, YEDGE_R, YMID_R_W, YEDGE_R_W )

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
! !DESCRIPTION: Subroutine Get\_nnMatch returns the number of HEMCO species
! that are also used in the atmospheric model. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_nnMatch( nnMatch, RC ) 
!
! !USES:
!
    USE HCO_Config_Mod, ONLY : Config_GetnSpecies
    USE HCO_Config_Mod, ONLY : Config_GetSpecNames
!
! !OUTPUT PARAMETERS:
!
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
    CALL Model_GetSpecies( nModelSpec,     ModelSpecNames,      &
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

       ! Register for output (through diagnostics)
       CALL Diagn_Create ( am_I_Root,                         &
                           HcoState,                          &
                           cName     = TRIM(HcoSpecNames(I)), &
                           ExtNr     = -1,                    &
                           Cat       = -1,                    &
                           Hier      = -1,                    &
                           HcoID     = CNT,                   &
                           SpaceDim  = 2,                     &
                           LevIDx    = -1,                    &
                           OutUnit   = 'kg/m2/s',             &
                           WriteFreq = 'Hourly',              &
                           AutoFill  = 1,                     &
                           cID       = cID,                   &
                           RC        = RC                      )
       IF ( RC /= HCO_SUCCESS ) RETURN

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
    USE HCO_Extlist_Mod, ONLY : GetExtOpt
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
    LOGICAL             :: FOUND
    CHARACTER(LEN=255)  :: MSG, LOC, DUM
    CHARACTER(LEN=255)  :: MyTimeFile 

    !=================================================================
    ! READ_TIME begins here
    !=================================================================

    ! For error handling
    LOC = 'READ_TIME (hcoi_standalone_mod.F90)'

    ! Try to get TimeFile from configuration file (in settings)
    CALL GetExtOpt ( 0, 'TimeFile', OptValChar=MyTimeFile, &
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

    ! Get emission time step
    READ( IU_FILE, '(A)', IOSTAT=IOS ) DUM
    IF ( IOS /= 0 ) THEN
       MSG = 'Error 2 reading ' // TRIM(TimeFile)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    LNG = LEN(TRIM(DUM))
 
    ! Get index after colon 
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
     
    ! Get number of simulation steps
    READ( IU_FILE, '(A)', IOSTAT=IOS ) DUM
    IF ( IOS /= 0 ) THEN
       MSG = 'Error 2 reading ' // TRIM(TimeFile)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    LNG = LEN(TRIM(DUM))
    LOW = NextCharPos ( TRIM(DUM), HCO_COL(), 1 )
    IF ( LOW < 0 .OR. LOW == LNG ) THEN
       MSG = 'Cannot extract index after colon: ' // TRIM(DUM)
       CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
       RETURN
    ENDIF
    LOW = LOW + 1
    READ ( DUM(LOW:LNG), * ) NSTEPS

    ! Allocate time arrays
    ALLOCATE ( YRS(NSTEPS), MTS(NSTEPS), DYS(NSTEPS), STAT=AS )
    IF ( AS/= 0 ) THEN
       CALL HCO_ERROR( 'Time alloc error 1', RC, THISLOC=LOC )
       RETURN
    ENDIF
    ALLOCATE ( HRS(NSTEPS), MNS(NSTEPS), SCS(NSTEPS), STAT=AS )
    IF ( AS/= 0 ) THEN
       CALL HCO_ERROR( 'Time alloc error 2', RC, THISLOC=LOC )
       RETURN
    ENDIF

    ! Set emission times
    DO N = 1, NSTEPS

       READ( IU_FILE, '(A)', IOSTAT=IOS ) DUM
       IF ( IOS /= 0 ) THEN
          MSG = 'Error reading time in ' // TRIM(TimeFile)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF
       LNG = LEN(TRIM(DUM))

       ! Times have to be stored as:
       ! YYYY-MM-DD HH:MM:SS
       ! --> read year from position 1:4, month from 6:7, etc.
       IF ( LNG < 19 ) THEN
          MSG = 'Provided time stamp is too short! ' // TRIM(DUM)
          CALL HCO_ERROR( MSG, RC, THISLOC=LOC )
          RETURN
       ENDIF

       READ ( DUM( 1: 4), * ) YRS(N) 
       READ ( DUM( 6: 7), * ) MTS(N) 
       READ ( DUM( 9:10), * ) DYS(N) 
       READ ( DUM(12:13), * ) HRS(N) 
       READ ( DUM(15:16), * ) MNS(N) 
       READ ( DUM(18:19), * ) SCS(N) 
    ENDDO

    ! Return w/ success
    RC = HCO_SUCCESS

  END SUBROUTINE Read_Time
!EOC
END MODULE HCO_StandAlone_Mod
