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
  USE Diagnostics_Mod
  USE ErrCode_Mod
  USE Precision_Mod
  USE Registry_Mod
  USE Species_Mod,   ONLY : Species
  USE State_Chm_Mod, ONLY : ChmState

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_State_Diag
  PUBLIC :: Cleanup_State_Diag
  PUBLIC :: Get_Metadata_State_Diag
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
     REAL(f4),  POINTER :: SpeciesConc(:,:,:,:) ! Spc Conc for diag output
     REAL(f4),  POINTER :: DryDepFlux (:,:,:,:) ! Dry deposition flux
     REAL(f4),  POINTER :: DryDepVel  (:,:,:,:) ! Dry deposition velocity

     !----------------------------------------------------------------------
     ! Registry of variables contained within State_Diag
     !----------------------------------------------------------------------
     CHARACTER(LEN=4)           :: State     = 'DIAG'   ! Name of this state
     TYPE(MetaRegItem), POINTER :: Registry  => NULL()  ! Registry object

  END TYPE DgnState
!
! !REMARKS:
!  TBD
!
! !REVISION HISTORY: 
!  05 Jul 2017 - R. Yantosca - Initial version
!  22 Sep 2017 - E. Lundgren - Fill in content to allocate State_Diag; add 
!                              subroutines to get metadata and register fields
!  26 Sep 2017 - E. Lundgren - Remove Lookup_State_Diag and Print_State_Diag
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
                              State_Chm, Diag_List,  State_Diag, RC )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    INTEGER,        INTENT(IN)    :: IM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: JM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: LM          ! # longitudes on this PET
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry state object
    TYPE(DgnList),  INTENT(IN)    :: Diag_List   ! Diagnostics list object

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
    CHARACTER(LEN=255)     :: ErrMsg, ThisLoc, Desc, Units
    CHARACTER(LEN=255)     :: diagID
    INTEGER                :: nSpecies, nAdvect, nDryDep, nKppSpc, nWetDep
    LOGICAL                :: EOF, Found

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

    !--------------------------------------------
    ! Species Concentration [v/v dry] 
    !--------------------------------------------
    diagID = 'SpeciesConc'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
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
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
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
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
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
    IF ( am_I_Root ) THEN
       WRITE( 6, 10 )
10     FORMAT( /, 'Registered variables contained within the State_Diag object:' )
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF
    CALL Registry_Print( am_I_Root   = am_I_Root,            &
                         Registry    = State_Diag%Registry,  &
                         ShortFormat = .TRUE.,               &
                         RC          = RC                      )

    ! Trap error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered!'
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
! !IROUTINE: Get_Metadata_State_Diag
!
! !DESCRIPTION: Subroutine GET\_METADATA\_STATE\_DIAG retrieves basic 
!  information about each State_Diag field.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Metadata_State_Diag( am_I_Root,  Name,  Desc,  Units, &
                                      PerSpecies, Rank,  Type,  VLoc,  &
                                      Found,      RC    )
!
! !USES:
!
    USE Registry_Params_Mod
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
    CHARACTER(LEN=255),  OPTIONAL    :: PerSpecies ! Max species wildcard
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
    isSpecies = PRESENT( PerSpecies )

    ! Set defaults for optional arguments. Assume type and vertical 
    ! location are real (flexible precision) and center unless specified 
    ! otherwise
    IF ( isUnits ) Units = ''
    IF ( isDesc  ) Desc  = ''              
    IF ( isRank  ) Rank  = -1              ! Initialize # dims as bad value 
    IF ( isType  ) Type  = KINDVAL_F4      ! Assume real*4 for diagnostics
    IF ( isVLoc  ) VLoc  = VLocationCenter ! Assume vertically centered
    IF ( isSpecies ) PerSpecies = ''         ! Assume not species-dependent

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
          IF ( isSpecies ) PerSpecies = 'ALL'

       CASE ( 'DRYDEPFLUX' )
          IF ( isDesc  )   Desc       = 'Dry deposition flux of species'
          IF ( isUnits )   Units      = 'molec cm-2 s-1'
          IF ( isRank  )   Rank       = 3
          IF ( isSpecies ) PerSpecies = 'DRY'

       CASE ( 'DRYDEPVEL' )
          IF ( isDesc  )   Desc       = 'Dry deposition velocity of species'
          IF ( isUnits )   Units      = 'cm s-1'
          IF ( isRank  )   Rank       = 3
          IF ( isSpecies ) PerSpecies = 'DRY'

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
    CHARACTER(LEN=255)     :: desc, units, perSpecies
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
                                  vloc=vloc,   perSpecies=perSpecies,       &
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
    IF ( perSpecies == 'ALL' ) THEN       

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
    ELSEIF ( perSpecies == 'ADV' ) THEN       

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
    ELSEIF ( perSpecies == 'DRY' ) THEN       
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
    ELSEIF ( perSpecies == 'WET' ) THEN       
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
    ELSEIF ( perSpecies == '' ) THEN
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
       PRINT *, "PerSpecies for diagnostic " // TRIM(Name) // &
                " is not definied in Register_DiagField_R4_2D: ", PerSpecies
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
    CHARACTER(LEN=255)     :: desc, units, perSpecies
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
                                  vloc=vloc,   perSpecies=perSpecies,       &
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
    IF ( perSpecies == 'ALL' ) THEN       

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
    ELSEIF ( perSpecies == 'ADV' ) THEN       

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
    ELSEIF ( perSpecies == 'DRY' ) THEN       
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
    ELSEIF ( perSpecies == 'WET' ) THEN       
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
    ELSEIF ( perSpecies == '' ) THEN
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
       PRINT *, "PerSpecies for diagnostic " // TRIM(Name) // &
                " is not definied in Register_DiagField_R4_3D: ", PerSpecies
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
    CHARACTER(LEN=255)     :: desc, units, perSpecies
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
                                  vloc=vloc,   perSpecies=perSpecies,       &
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
    IF ( perSpecies == 'ALL' ) THEN       

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
    ELSEIF ( perSpecies == 'ADV' ) THEN       

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
    ELSEIF ( perSpecies == 'DRY' ) THEN       
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
    ELSEIF ( perSpecies == 'WET' ) THEN       
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
       PRINT *, "PerSpecies for diagnostic " // TRIM(Name) // &
                " is not definied in Register_DiagField_R4_4D: ", PerSpecies
       RETURN

    ENDIF

  END SUBROUTINE Register_DiagField_R4_4D
!EOC
END MODULE State_Diag_Mod
