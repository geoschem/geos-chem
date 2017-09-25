!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: state_diag_mod.F90
!
! !DESCRIPTION: Module STATE\_DIAG\_MOD contains the derived type
!  used to define the Diagnostics State object for GEOS-Chem.
!\\
!\\
!  This module also contains the routines that allocate and deallocate memory
!  to the Diagnostics State object.  The Diagnostics State object is not
!  defined in this module.  It must be be declared as variable in the top-level
!  driver routine, and then passed to lower-level routines as an argument.
!\\
!\\
! !INTERFACE:
!
MODULE State_Diag_Mod
!
! USES:
!
  USE ErrCode_Mod
  USE Precision_Mod
  USE Registry_Mod, ONLY : MetaRegItem

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_State_Diag
  PUBLIC :: Get_Metadata_State_Diag
  PUBLIC :: Lookup_State_Diag
  PUBLIC :: Print_State_Diag
  PUBLIC :: Cleanup_State_Diag
!
! !PRIVATE MEMBER FUNCTIONS
!
  PRIVATE :: Register_DiagField
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Diagnostics State
  !=========================================================================
  TYPE, PUBLIC :: DgnState

     !----------------------------------------------------------------------
     ! Diagnostic arrays
     !----------------------------------------------------------------------
     REAL(f4),  POINTER :: F_Adv_EW(:,:,:,:)
     REAL(f4),  POINTER :: F_Adv_NS(:,:,:,:)
     REAL(f4),  POINTER :: F_Adv_UP(:,:,:,:)

     ! ewl 
     REAL(f4),  POINTER :: SpeciesConc(:,:,:,:) ! Spc Conc for diag output
     REAL(f4),  POINTER :: DryDepFlux (:,:,:,:) ! Dry deposition flux
     REAL(f4),  POINTER :: DryDepVel  (:,:,:,:) ! Dry deposition velocity

     !----------------------------------------------------------------------
     ! Registry of variables contained within State_Diag
     !----------------------------------------------------------------------
     CHARACTER(LEN=4)           :: State     = 'DIAG'   ! Name of this state
     TYPE(MetaRegItem), POINTER :: Registry  => NULL()  ! Registry object

  END TYPE DgnState

  !=========================================================================
  ! Derived type for Diagnostics Item (unique item in HISTORY.rc)
  !=========================================================================
  TYPE, PUBLIC :: DgnItem
     CHARACTER(LEN=255)      :: name 
     CHARACTER(LEN=255)      :: category 
     LOGICAL                 :: isWildcard
     LOGICAL                 :: isSpecies
     CHARACTER(LEN=255)      :: wildcard
     CHARACTER(LEN=255)      :: species
     TYPE(DgnItem), POINTER  :: next
  END TYPE DgnItem

  !=========================================================================
  ! Derived type for Diagnostics List (unique items in HISTORY.rc)
  !=========================================================================
  TYPE, PUBLIC :: DgnList
     TYPE(DgnItem), POINTER  :: head
     INTEGER                 :: numItems
  END TYPE DgnList
!
! !REMARKS:
!  TBD
!
! !REVISION HISTORY: 
!  05 Jul 2017 - R. Yantosca - Initial version
!  22 Sep 2017 - E. Lundgren - Fill in content to allocate State_Diag
!
!EOC
!------------------------------------------------------------------------------
!BOC
!
! !MODULE INTERFACES:
!
  INTERFACE Register_DiagField
     MODULE PROCEDURE Register_DiagField_R4_2D
     MODULE PROCEDURE Register_DiagField_R4_3D
     MODULE PROCEDURE Register_DiagField_R4_4D
  END INTERFACE Register_DiagField
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_State_Diag
!
! !DESCRIPTION: Subroutine INIT\_STATE\_DIAG allocates all fields of
!  the diagnostics state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_State_Diag( am_I_Root, IM, JM, LM, Input_Opt, &
                              State_Chm, State_Diag, RC )
!
! !USES:
!
    USE Diagnostics_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Chm_Mod, ONLY : ChmState
    USe Registry_Mod,  ONLY : Registry_AddField
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    INTEGER,        INTENT(IN)    :: IM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: JM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: LM          ! # longitudes on this PET
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostic State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REMARKS:
!  For consistency, maybe this should be moved to a different module.
!
! !REVISION HISTORY: 
!  05 Jul 2017 - R. Yantosca -  Initial version
!  22 Sep 2017 - E. Lundgren -  Fill in content
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255)     :: ErrMsg, ThisLoc, Desc, Units, historyConfigFile
    CHARACTER(LEN=255)     :: diagID
    INTEGER                :: nSpecies, nAdvect, nDryDep, nKppSpc, nWetDep
    LOGICAL                :: EOF, Found
    TYPE(DgnList), POINTER :: DiagList

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      =  GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Init_State_Diag (in Headers/state_diag_mod.F90)'
    Desc    = ''
    Units   = ''
    Found   = .False.

    ! Number of species per category
    nSpecies = State_Chm%nSpecies
    nDrydep  = State_Chm%nAdvect
    nDryDep  = State_Chm%nDryDep
    nKppSpc  = State_Chm%nKppSpc
    nWetDep  = State_Chm%nWetDep

    ! Start with these
    State_Diag%SpeciesConc => NULL()
    State_Diag%DryDepFlux  => NULL()
    State_Diag%DryDepVel   => NULL()

    ! placeholder
    State_Diag%F_Adv_EW => NULL()
    State_Diag%F_Adv_NS => NULL()
    State_Diag%F_Adv_UP => NULL()

    !---------------------------------------------------------
    ! Initialize the list of unique entries in HISTORY.rc
    !---------------------------------------------------------
    historyConfigFile = 'HISTORY.rc' ! eventually use Input_Opt variable
    CALL Init_DiagList( am_I_Root, historyConfigFile, DiagList, RC )

    !--------------------------------------------
    ! Species Concentration [v/v dry] 
    !--------------------------------------------
    diagID = 'SpeciesConc'
    CALL Check_DiagList( am_I_Root, DiagList, diagID, Found, RC )
    IF ( Found ) THEN
       PRINT *, 'Allocating ' // TRIM(diagID) // ' diagnostic'
       ALLOCATE( State_Diag%SpeciesConc( IM, JM, LM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( 'State_Diag%' // TRIM(diagID), 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%SpeciesConc = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%SpeciesConc, &
                                State_Chm, State_Diag, RC )
    ENDIF

    !--------------------------------------------
    ! Dry deposition flux
    !--------------------------------------------
    diagID = 'DryDepFlux'
    CALL Check_DiagList( am_I_Root, DiagList, diagID, Found, RC )
    IF ( Found ) THEN
       PRINT *, 'Allocating ' // TRIM(diagID) // ' diagnostic'
       ALLOCATE( State_Diag%DryDepFlux( IM, JM, LM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'State_Diag%' // TRIM(diagID), 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepFlux = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepFlux, &
                                State_Chm, State_Diag, RC )
    ENDIF

    !--------------------------------------------
    ! Dry deposition velocity
    !-------------------------------------------- 
    diagID = 'DryDepVel'
    CALL Check_DiagList( am_I_Root, DiagList, diagID, Found, RC )
    IF ( Found ) THEN
       PRINT *, 'Allocating ' // TRIM(diagID) // 'diagnostic'
       ALLOCATE( State_Diag%DryDepVel( IM, JM, LM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'State_Diag%' // TRIM(diagID), 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepVel = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepVel, &
                                State_Chm, State_Diag, RC )
    ENDIF

    !=======================================================================
    ! Print information about the registered fields (short format)
    !=======================================================================
    CALL Print_State_Diag( am_I_Root, State_Diag, RC, ShortFormat=.TRUE.)
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Print_State_Met"'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Init_State_Diag
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_State_Diag
!
! !DESCRIPTION: Subroutine CLEANUP\_STATE\_DIAG deallocates all fields 
!  of the meteorology state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_State_Diag( am_I_Root, State_Diag, RC )
!
! !USES:
!
    USE Registry_Mod, ONLY : Registry_Destroy
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Return code
!
! !REVISION HISTORY: 
!  05 Jul 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !========================================================================
    ! Initialize
    !========================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> Cleanup_State_Diag (in Headers/state_diag_mod.F90)'

    !=======================================================================
    ! Deallocate module variables
    !=======================================================================

    ! For now, just nullify
    State_Diag%F_Adv_EW => NULL()
    State_Diag%F_Adv_NS => NULL()
    State_Diag%F_Adv_UP => NULL()

    IF ( ASSOCIATED(State_Diag%SpeciesConc)) DEALLOCATE( State_Diag%SpeciesConc)
    IF ( ASSOCIATED(State_Diag%DryDepFlux )) DEALLOCATE( State_Diag%DryDepFlux )
    IF ( ASSOCIATED(State_Diag%DryDepVel  )) DEALLOCATE( State_Diag%DryDepVel  )

    !=======================================================================
    ! Destroy the registry of fields for this module
    !=======================================================================
    CALL Registry_Destroy( am_I_Root, State_Diag%Registry, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not destroy registry object State_Diag%Registry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Cleanup_State_Diag
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Print_State_Diag
!
! !DESCRIPTION: Print information about all the registered variables
!  contained within the State\_Diag object.  This is basically a wrapper for
!  routine REGISTRY\_PRINT in registry\_mod.F90.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_State_Diag( am_I_Root, State_Diag, RC, ShortFormat )
!
! !USES:
!
    USE Registry_Mod, ONLY : Registry_Print
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)  :: am_I_Root    ! Root CPU?  
    TYPE(DgnState), INTENT(IN)  :: State_Diag   ! Meteorology State object
    LOGICAL,        OPTIONAL    :: ShortFormat  ! Print truncated info
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC           ! Success/failure?
!
! !REVISION HISTORY:
!  29 Jun 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Print_State_Diag (in Headers/state_diag_mod.F90)'

    !=======================================================================
    ! Print info about registered variables
    !=======================================================================

    ! Header line
    IF ( am_I_Root ) THEN
       WRITE( 6, 10 )
10     FORMAT( /, 'Registered variables contained within the State_Diag object:' )
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF

    ! Print registry info in truncated format
    CALL Registry_Print( am_I_Root   = am_I_Root,                            &
                         Registry    = State_Diag%Registry,                  &
                         ShortFormat = ShortFormat,                          &
                         RC          = RC                                   )

    ! Trap error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "Registry_Print"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Print_State_Diag
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Lookup_State_Diag
!
! !DESCRIPTION: Return metadata and/or a pointer to the data for any
!  variable contained within the State\_Diag object by searching for its name.
!  This is basically a wrapper for routine REGISTRY\_LOOKUP in 
!  registry\_mod.F90.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Lookup_State_Diag( am_I_Root,  State_Diag,   Variable,          & 
                                RC,         Description,  Dimensions,        &
                                KindVal,    MemoryInKb,   Rank,              & 
                                Units,      OnLevelEdges, Ptr2d_4,           &
                                Ptr3d_4,    Ptr2d_I,      Ptr3d_I           )
!
! !USES:
!
    USE Registry_Mod, ONLY : Registry_Lookup, To_Uppercase
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)  :: am_I_Root       ! Is this the root CPU? 
    TYPE(DgnState),   INTENT(IN)  :: State_Diag      ! Diagnostics State
    CHARACTER(LEN=*), INTENT(IN)  :: Variable        ! Variable name
!
! !OUTPUT PARAMETERS:
!
    ! Required outputs
    INTEGER,          INTENT(OUT) :: RC              ! Success or failure?

    ! Optional outputs
    CHARACTER(LEN=255),  OPTIONAL :: Description     ! Description of data
    INTEGER,             OPTIONAL :: Dimensions(3)   ! Dimensions of data
    INTEGER,             OPTIONAL :: KindVal         ! Numerical KIND value
    REAL(fp),            OPTIONAL :: MemoryInKb      ! Memory usage
    INTEGER,             OPTIONAL :: Rank            ! Size of data
    CHARACTER(LEN=255),  OPTIONAL :: Units           ! Units of data
    LOGICAL,             OPTIONAL :: OnLevelEdges    ! =T if data is defined
                                                     !    on level edges
                                                     ! =F if on centers

    ! Pointers to data
    REAL(f4),   POINTER, OPTIONAL :: Ptr2d_4(:,:  )  ! 2D 4-byte data
    REAL(f4),   POINTER, OPTIONAL :: Ptr3d_4(:,:,:)  ! 3D 4-byte data
    INTEGER,    POINTER, OPTIONAL :: Ptr2d_I(:,:  )  ! 2D integer data
    INTEGER,    POINTER, OPTIONAL :: Ptr3d_I(:,:,:)  ! 3D integer data
!
! !REMARKS:
!  We keep the StateName variable private to this module. Users only have
!  to supply the name of each module variable.
!
! !REVISION HISTORY:
!  05 Jul 2017 - R. Yantosca - Initial version
!  25 Sep 2017 - E. Lundgren - Use all caps when looking up name in registry
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, name_AllCaps

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Lookup_State_Diag (in Headers/state_diag_mod.F90)'

    ! Convert to all uppercase when looking in registry
    name_AllCaps = To_Uppercase( Variable )

    !=======================================================================
    ! Look up a variable; Return metadata and/or a pointer to the data
    !=======================================================================
    CALL Registry_Lookup( am_I_Root    = am_I_Root,                          &
                          Registry     = State_Diag%Registry,                &
                          State        = State_Diag%State,                   &
                          Variable     = name_allCaps,                       &
                          Description  = Description,                        &
                          Dimensions   = Dimensions,                         &
                          KindVal      = KindVal,                            &
                          MemoryInKb   = MemoryInKb,                         &
                          Rank         = Rank,                               &
                          Units        = Units,                              &
                          OnLevelEdges = OnLevelEdges,                       &
                          Ptr2d_4      = Ptr2d_4,                            &
                          Ptr3d_4      = Ptr3d_4,                            &
                          Ptr2d_I      = Ptr2d_I,                            &
                          Ptr3d_I      = Ptr3d_I,                            &
                          RC           = RC                                 )

    ! Trap error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not find variable "' // TRIM( Variable ) // &
               '" in the State_Diag registry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Lookup_State_Diag
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Metadata_State_Diag
!
! !DESCRIPTION: Subroutine GET\_METADATA\_STATE\_DIAG retrieves basic 
!  information about each State_Diag field.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Metadata_State_Diag( am_I_Root,  Name,  Desc,  Units, &
                                      SpeciesSet, Rank,  Type,  VLoc,  &
                                      Found,      RC    )
!
! !USES:
!
    USE Registry_Params_Mod
    USE Registry_Mod, ONLY : To_Uppercase
!
! !INPUT PARAMETERS:
! 
    LOGICAL,             INTENT(IN)  :: am_I_Root  ! Is this the root CPU?
    CHARACTER(LEN=*),    INTENT(IN)  :: Name       ! Sate_Met field name
!
! !OUTPUT PARAMETERS:
!
    CHARACTER(LEN=255),  OPTIONAL    :: Desc       ! Long name string
    CHARACTER(LEN=255),  OPTIONAL    :: Units      ! Units string
    CHARACTER(LEN=255),  OPTIONAL    :: SpeciesSet ! Max species wildcard
    INTEGER,             OPTIONAL    :: Rank       ! # of dimensions
    INTEGER,             OPTIONAL    :: Type       ! Desc of data type
    INTEGER,             OPTIONAL    :: VLoc       ! Vertical placement
    LOGICAL,             OPTIONAL    :: Found      ! Item found?
    INTEGER,             OPTIONAL    :: RC         ! Return code
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  20 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, Name_AllCaps
    LOGICAL            :: isDesc, isUnits, isRank, isType, isVLoc, isSpecies
    
    !=======================================================================
    ! Initialize
    !=======================================================================

    ! Assume success
    RC    =  GC_SUCCESS
    ThisLoc = ' -> at Get_Metadata_State_Diag (in Headers/state_diag_mod.F90)'
    Found = .True.

    ! Optional arguments present?
    isDesc    = PRESENT( Desc  )
    isUnits   = PRESENT( Units )
    isRank    = PRESENT( Rank  )
    isType    = PRESENT( Type  )
    isVLoc    = PRESENT( VLoc  )
    isSpecies = PRESENT( SpeciesSet )

    ! Set defaults for optional arguments. Assume type and vertical 
    ! location are real (flexible precision) and center unless specified 
    ! otherwise
    IF ( isUnits ) Units = ''
    IF ( isDesc  ) Desc  = ''              
    IF ( isRank  ) Rank  = -1              ! Initialize # dims as bad value 
    IF ( isType  ) Type  = KINDVAL_F4      ! Assume real*4 for diagnostics
    IF ( isVLoc  ) VLoc  = VLocationCenter ! Assume vertically centered
    IF ( isSpecies ) SpeciesSet = ''         ! Assume not species-dependent

    ! Convert name to uppercase
    Name_AllCaps = To_Uppercase( TRIM( Name ) )

    !=======================================================================
    ! Values for Retrieval (string comparison slow but happens only once)
    !=======================================================================
    SELECT CASE ( TRIM( Name_AllCaps ) )

       CASE ( 'SPECIESCONC' )
          IF ( isDesc  )   Desc       = 'Dry mixing ratio of species'
          IF ( isUnits )   Units      = 'mol mol-1 dry'
          IF ( isRank  )   Rank       = 3
          IF ( isSpecies ) SpeciesSet = 'ALL'

       CASE ( 'DRYDEPFLUX' )
          IF ( isDesc  )   Desc       = 'Dry deposition flux of species'
          IF ( isUnits )   Units      = 'molec cm-2 s-1'
          IF ( isRank  )   Rank       = 3
          IF ( isSpecies ) SpeciesSet = 'DRY'

       CASE ( 'DRYDEPVEL' )
          IF ( isDesc  )   Desc       = 'Dry deposition velocity of species'
          IF ( isUnits )   Units      = 'cm s-1'
          IF ( isRank  )   Rank       = 3
          IF ( isSpecies ) SpeciesSet = 'DRY'

       CASE DEFAULT
          ! Need to add better error handling
          Found = .False.
          PRINT *, 'WARNING: Metadata not found for State_Diag field: ' // TRIM( Name )
          RETURN

    END SELECT

    ! Set VLoc to undefined if variable is 2d
    IF ( isVLoc .AND. Rank == 2 ) THEN
       VLoc = VLocationNone
    ENDIF

   END SUBROUTINE Get_Metadata_State_Diag
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_DiagField_R4_2D
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R4_2D( am_I_Root,  Name,       Ptr2Data,  &
                                       State_Chm,  State_Diag, RC )
!
! !USES:
!
    USE Species_Mod,   ONLY : Species
    USE State_Chm_Mod, ONLY : ChmState
    USE Registry_Mod,  ONLY : Registry_AddField
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: Name            ! diagnostics name
    TYPE(ChmState),    INTENT(IN)    :: State_Chm       ! Obj for chem state
    REAL*4,            POINTER       :: Ptr2Data(:,:)   ! pointer to data
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),    INTENT(INOUT) :: State_Diag      ! Obj for diag state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC              ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: ErrMsg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, speciesSet
    CHARACTER(LEN=255)     :: thisSpcDiagName, thisSpcDiagDesc
    INTEGER                :: N, D
    INTEGER                :: rank, type,  vloc
    LOGICAL                :: found
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize
    RC             = GC_SUCCESS
    ThisLoc        = ' -> at Register_DiagField_R4_2D ' // &
                     '(in Headers/state_diag_mod.F90)'

    CALL Get_Metadata_State_Diag( am_I_Root,   Name,           desc=desc,   &
                                  units=units, rank=rank,      type=type,   &
                                  vloc=vloc,   speciesSet=speciesSet,       &
                                  Found=Found, RC=RC     )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in Get_DiagField_Metadata'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Check that metadata consistent with data pointer
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // TRIM(Name)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! If tied to all species then register each one
    IF ( speciesSet == 'ALL' ) THEN       

       DO N = 1, State_Chm%nSpecies
          SpcInfo  => State_Chm%SpcData(N)%Info
          thisSpcDiagName = TRIM( Name ) // '_' // TRIM( SpcInfo%Name )
          thisSpcDiagDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcDiagName,      &
                                  Description  = thisSpcDiagDesc,      &
                                  Units        = units,                &
                                  Data1d_4     = Ptr2Data(:,N),    &
                                  RC           = RC                   )
          SpcInfo => NULL()
          CALL GC_CheckVar( 'State_Diag%' // TRIM(Name), 1, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in Registry_AddField'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If tied to advected species then register advected species only
    ELSEIF ( speciesSet == 'ADV' ) THEN       

       DO N = 1, State_Chm%nAdvect
          SpcInfo  => State_Chm%SpcData(N)%Info
          thisSpcDiagName = TRIM( Name ) // '_' // TRIM( SpcInfo%Name )
          thisSpcDiagDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcDiagName,      &
                                  Description  = thisSpcDiagDesc,      &
                                  Units        = units,                &
                                  Data1d_4     = Ptr2Data(:,N),        &
                                  RC           = RC                   )
          SpcInfo => NULL()
          CALL GC_CheckVar( 'State_Diag%' // TRIM(Name), 1, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in Registry_AddField'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If tied to drydep species then register drydep species only
    ELSEIF ( speciesSet == 'DRY' ) THEN       
       DO D = 1, State_Chm%nDryDep

          ! Get species number and pointer to the species database
          N =  State_Chm%Map_DryDep(D)
          SpcInfo  => State_Chm%SpcData(N)%Info

          ! Set the registry name for this diagnostic and species
          thisSpcDiagName = TRIM( Name ) // '_' // TRIM( SpcInfo%Name )
          thisSpcDiagDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcDiagName,      &
                                  Description  = thisSpcDiagDesc,      &
                                  Units        = units,                &
                                  Data1d_4     = Ptr2Data(:,N),        &
                                  RC           = RC                   )
          SpcInfo => NULL()
          CALL GC_CheckVar( 'State_Diag%' // TRIM(Name), 1, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in Registry_AddField'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If tied to wetdep species then register drydep species only
    ELSEIF ( speciesSet == 'WET' ) THEN       
       DO D = 1, State_Chm%nWetDep

          ! Get species number and pointer to the species database
          N =  State_Chm%Map_WetDep(D)
          SpcInfo  => State_Chm%SpcData(N)%Info

          ! Set the registry name for this diagnostic and species
          thisSpcDiagName = TRIM( Name ) // '_' // TRIM( SpcInfo%Name )
          thisSpcDiagDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcDiagName,      &
                                  Description  = thisSpcDiagDesc,      &
                                  Units        = units,                &
                                  Data1d_4     = Ptr2Data(:,N),        &
                                  RC           = RC                   )
          SpcInfo => NULL()
          CALL GC_CheckVar( 'State_Diag%' // TRIM(Name), 1, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in Registry_AddField'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If not tied to species then simply add the single field
    ELSEIF ( speciesSet == '' ) THEN
       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Diag%Registry,  &
                               State        = State_Diag%State,     &
                               Variable     = Name,                 &
                               Description  = desc,                 &
                               Units        = units,                &
                               Data2d_4     = Ptr2Data,             &
                               RC           = RC                   )
       CALL GC_CheckVar( 'State_Diag%' // TRIM(Name), 1, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in Registry_AddField'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ELSE
       ! TODO: Add better error handling
       PRINT *, "SpeciesSet for diagnostic " // TRIM(Name) // &
                " is not definied in Register_DiagField_R4_2D: ", SpeciesSet
       RETURN

    ENDIF

  END SUBROUTINE Register_DiagField_R4_2D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_DiagField_R4_3D
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R4_3D( am_I_Root,  Name,       Ptr2Data,  &
                                       State_Chm,  State_Diag, RC )
!
! !USES:
!
    USE Species_Mod,   ONLY : Species
    USE State_Chm_Mod, ONLY : ChmState
    USE Registry_Mod,  ONLY : Registry_AddField
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: Name            ! diagnostics name
    TYPE(ChmState),    INTENT(IN)    :: State_Chm       ! Obj for chem state
    REAL*4,            POINTER       :: Ptr2Data(:,:,:) ! pointer to data
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),    INTENT(INOUT) :: State_Diag      ! Obj for diag state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC              ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: ErrMsg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, speciesSet
    CHARACTER(LEN=255)     :: thisSpcDiagName, thisSpcDiagDesc
    INTEGER                :: N, D
    INTEGER                :: rank, type,  vloc
    LOGICAL                :: found
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize
    RC             = GC_SUCCESS
    ThisLoc        = ' -> at Register_DiagField_R4_3D ' // &
                     '(in Headers/state_diag_mod.F90)'

    CALL Get_Metadata_State_Diag( am_I_Root,   Name,           desc=desc,   &
                                  units=units, rank=rank,      type=type,   &
                                  vloc=vloc,   speciesSet=speciesSet,       &
                                  Found=Found, RC=RC     )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in Get_DiagField_Metadata'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Check that metadata consistent with data pointer
    IF ( rank /= 2 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // TRIM(Name)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! If tied to all species then register each one
    IF ( speciesSet == 'ALL' ) THEN       

       DO N = 1, State_Chm%nSpecies
          SpcInfo  => State_Chm%SpcData(N)%Info
          thisSpcDiagName = TRIM( Name ) // '_' // TRIM( SpcInfo%Name )
          thisSpcDiagDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcDiagName,      &
                                  Description  = thisSpcDiagDesc,      &
                                  Units        = units,                &
                                  Data2d_4     = Ptr2Data(:,:,N),    &
                                  RC           = RC                   )
          SpcInfo => NULL()
          CALL GC_CheckVar( 'State_Diag%' // TRIM(Name), 1, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in Registry_AddField'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If tied to advected species then register advected species only
    ELSEIF ( speciesSet == 'ADV' ) THEN       

       DO N = 1, State_Chm%nAdvect
          SpcInfo  => State_Chm%SpcData(N)%Info
          thisSpcDiagName = TRIM( Name ) // '_' // TRIM( SpcInfo%Name )
          thisSpcDiagDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcDiagName,      &
                                  Description  = thisSpcDiagDesc,      &
                                  Units        = units,                &
                                  Data2d_4     = Ptr2Data(:,:,N),    &
                                  RC           = RC                   )
          SpcInfo => NULL()
          CALL GC_CheckVar( 'State_Diag%' // TRIM(Name), 1, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in Registry_AddField'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If tied to drydep species then register drydep species only
    ELSEIF ( speciesSet == 'DRY' ) THEN       
       DO D = 1, State_Chm%nDryDep

          ! Get species number and pointer to the species database
          N =  State_Chm%Map_DryDep(D)
          SpcInfo  => State_Chm%SpcData(N)%Info

          ! Set the registry name for this diagnostic and species
          thisSpcDiagName = TRIM( Name ) // '_' // TRIM( SpcInfo%Name )
          thisSpcDiagDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcDiagName,      &
                                  Description  = thisSpcDiagDesc,      &
                                  Units        = units,                &
                                  Data2d_4     = Ptr2Data(:,:,N),    &
                                  RC           = RC                   )
          SpcInfo => NULL()
          CALL GC_CheckVar( 'State_Diag%' // TRIM(Name), 1, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in Registry_AddField'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If tied to wetdep species then register drydep species only
    ELSEIF ( speciesSet == 'WET' ) THEN       
       DO D = 1, State_Chm%nWetDep

          ! Get species number and pointer to the species database
          N =  State_Chm%Map_WetDep(D)
          SpcInfo  => State_Chm%SpcData(N)%Info

          ! Set the registry name for this diagnostic and species
          thisSpcDiagName = TRIM( Name ) // '_' // TRIM( SpcInfo%Name )
          thisSpcDiagDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcDiagName,      &
                                  Description  = thisSpcDiagDesc,      &
                                  Units        = units,                &
                                  Data2d_4     = Ptr2Data(:,:,N),    &
                                  RC           = RC                   )
          SpcInfo => NULL()
          CALL GC_CheckVar( 'State_Diag%' // TRIM(Name), 1, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in Registry_AddField'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If not tied to species then simply add the single field
    ELSEIF ( speciesSet == '' ) THEN
       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Diag%Registry,  &
                               State        = State_Diag%State,     &
                               Variable     = Name,                 &
                               Description  = desc,                 &
                               Units        = units,                &
                               Data3d_4     = Ptr2Data,             &
                               RC           = RC                   )
       CALL GC_CheckVar( 'State_Diag%' // TRIM(Name), 1, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in Registry_AddField'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ELSE
       ! TODO: Add better error handling
       PRINT *, "SpeciesSet for diagnostic " // TRIM(Name) // &
                " is not definied in Register_DiagField_R4_3D: ", SpeciesSet
       RETURN

    ENDIF

  END SUBROUTINE Register_DiagField_R4_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Register_DiagField_R4_4D
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Register_DiagField_R4_4D( am_I_Root,  Name,       Ptr2Data,  &
                                       State_Chm,  State_Diag, RC )
!
! !USES:
!
    USE Species_Mod,   ONLY : Species
    USE State_Chm_Mod, ONLY : ChmState
    USE Registry_Mod,  ONLY : Registry_AddField
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: Name            ! diagnostics name
    TYPE(ChmState),    INTENT(IN)    :: State_Chm       ! Obj for chem state
    REAL*4,            POINTER       :: Ptr2Data(:,:,:,:) ! pointer to data
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState),    INTENT(INOUT) :: State_Diag      ! Obj for diag state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC              ! Success/failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  20 Sep 2017 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!   
    CHARACTER(LEN=255)     :: ErrMsg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, speciesSet
    CHARACTER(LEN=255)     :: thisSpcDiagName, thisSpcDiagDesc
    INTEGER                :: N, D
    INTEGER                :: rank, type,  vloc
    LOGICAL                :: found
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize
    RC             = GC_SUCCESS
    ThisLoc        = ' -> at Register_DiagField_R4_4D ' // &
                     '(in Headers/state_diag_mod.F90)'

    CALL Get_Metadata_State_Diag( am_I_Root,   Name,           desc=desc,   &
                                  units=units, rank=rank,      type=type,   &
                                  vloc=vloc,   speciesSet=speciesSet,       &
                                  Found=Found, RC=RC     )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in Get_DiagField_Metadata'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Check that metadata consistent with data pointer
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // TRIM(Name)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! If tied to all species then register each one
    IF ( speciesSet == 'ALL' ) THEN       

       DO N = 1, State_Chm%nSpecies
          SpcInfo  => State_Chm%SpcData(N)%Info
          thisSpcDiagName = TRIM( Name ) // '_' // TRIM( SpcInfo%Name )
          thisSpcDiagDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcDiagName,      &
                                  Description  = thisSpcDiagDesc,      &
                                  Units        = units,                &
                                  Data3d_4     = Ptr2Data(:,:,:,N),    &
                                  RC           = RC                   )
          SpcInfo => NULL()
          CALL GC_CheckVar( 'State_Diag%' // TRIM(Name), 1, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in Registry_AddField'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If tied to advected species then register advected species only
    ELSEIF ( speciesSet == 'ADV' ) THEN       

       DO N = 1, State_Chm%nAdvect
          SpcInfo  => State_Chm%SpcData(N)%Info
          thisSpcDiagName = TRIM( Name ) // '_' // TRIM( SpcInfo%Name )
          thisSpcDiagDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcDiagName,      &
                                  Description  = thisSpcDiagDesc,      &
                                  Units        = units,                &
                                  Data3d_4     = Ptr2Data(:,:,:,N),    &
                                  RC           = RC                   )
          SpcInfo => NULL()
          CALL GC_CheckVar( 'State_Diag%' // TRIM(Name), 1, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in Registry_AddField'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If tied to drydep species then register drydep species only
    ELSEIF ( speciesSet == 'DRY' ) THEN       
       DO D = 1, State_Chm%nDryDep

          ! Get species number and pointer to the species database
          N =  State_Chm%Map_DryDep(D)
          SpcInfo  => State_Chm%SpcData(N)%Info

          ! Set the registry name for this diagnostic and species
          thisSpcDiagName = TRIM( Name ) // '_' // TRIM( SpcInfo%Name )
          thisSpcDiagDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcDiagName,      &
                                  Description  = thisSpcDiagDesc,      &
                                  Units        = units,                &
                                  Data3d_4     = Ptr2Data(:,:,:,N),    &
                                  RC           = RC                   )
          SpcInfo => NULL()
          CALL GC_CheckVar( 'State_Diag%' // TRIM(Name), 1, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in Registry_AddField'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If tied to wetdep species then register drydep species only
    ELSEIF ( speciesSet == 'WET' ) THEN       
       DO D = 1, State_Chm%nWetDep

          ! Get species number and pointer to the species database
          N =  State_Chm%Map_WetDep(D)
          SpcInfo  => State_Chm%SpcData(N)%Info

          ! Set the registry name for this diagnostic and species
          thisSpcDiagName = TRIM( Name ) // '_' // TRIM( SpcInfo%Name )
          thisSpcDiagDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcDiagName,      &
                                  Description  = thisSpcDiagDesc,      &
                                  Units        = units,                &
                                  Data3d_4     = Ptr2Data(:,:,:,N),    &
                                  RC           = RC                   )
          SpcInfo => NULL()
          CALL GC_CheckVar( 'State_Diag%' // TRIM(Name), 1, RC )
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in Registry_AddField'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ELSE
       ! TODO: Add better error handling
       PRINT *, "SpeciesSet for diagnostic " // TRIM(Name) // &
                " is not definied in Register_DiagField_R4_4D: ", SpeciesSet
       RETURN

    ENDIF

  END SUBROUTINE Register_DiagField_R4_4D
!EOC
END MODULE State_Diag_Mod
