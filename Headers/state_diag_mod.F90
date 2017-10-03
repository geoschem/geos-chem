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
     REAL(f4),  POINTER :: JValues    (:,:,:,:) ! Photolysis rates

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
!                              subroutines to get metadata and interface to
!                              register fields
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
    CHARACTER(LEN=255)     :: ErrMsg, ThisLoc
    CHARACTER(LEN=255)     :: diagID
    INTEGER                :: nSpecies, nAdvect, nDryDep, nKppSpc, nWetDep
    LOGICAL                :: EOF, Found

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      =  GC_SUCCESS
    ThisLoc = ' -> at Init_State_Diag (in Headers/state_diag_mod.F90)'
    Found   = .FALSE.

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
    State_Diag%JValues     => NULL()

#if defined( NC_DIAG )
    !--------------------------------------------
    ! Species Concentration
    !--------------------------------------------
    diagID = 'SpeciesConc'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       PRINT *, 'Allocating ' // TRIM(diagID) // ' diagnostic'
       ALLOCATE( State_Diag%SpeciesConc( IM, JM, LM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( 'State_Diag%SpeciesConc', 0, RC )
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
       CALL GC_CheckVar( 'State_Diag%DryDepFlux', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepFlux = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepFlux, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! Dry deposition velocity
    !-------------------------------------------- 
    diagID = 'DryDepVel'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       PRINT *, 'Allocating ' // TRIM(diagID) // ' diagnostic'
       ALLOCATE( State_Diag%DryDepVel( IM, JM, LM, nDryDep ), STAT=RC )
       CALL GC_CheckVar( 'State_Diag%DryDepVel', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%DryDepVel = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%DryDepVel, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    !--------------------------------------------
    ! J-Values
    !-------------------------------------------- 
    ! TODO: Mapping needs work
    diagID = 'JValues'
    CALL Check_DiagList( am_I_Root, Diag_List, diagID, Found, RC )
    IF ( Found ) THEN
       PRINT *, 'Allocating ' // TRIM(diagID) // ' diagnostic'
       ALLOCATE( State_Diag%JValues( IM, JM, LM, nSpecies ), STAT=RC )
       CALL GC_CheckVar( 'State_Diag%JValues', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       State_Diag%JValues = 0.0_f4
       CALL Register_DiagField( am_I_Root, diagID, State_Diag%JValues, &
                                State_Chm, State_Diag, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
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
#endif

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
    IF ( ASSOCIATED(State_Diag%JValues    )) DEALLOCATE( State_Diag%JValues    )

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
  SUBROUTINE Get_Metadata_State_Diag( am_I_Root, Name,   Found,      RC,   &
                                      Desc,      Units,  PerSpecies, Rank, &
                                      Type,      VLoc )
!
! !USES:
!
    USE Registry_Params_Mod
!
! !INPUT PARAMETERS:
! 
    LOGICAL,             INTENT(IN)  :: am_I_Root  ! Is this the root CPU?
    CHARACTER(LEN=*),    INTENT(IN)  :: Name       ! State_Diag field name
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,             INTENT(OUT)           :: Found      ! Item found?
    INTEGER,             INTENT(OUT)           :: RC         ! Return code
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: Desc       ! Long name string
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: Units      ! Units string
    CHARACTER(LEN=255),  INTENT(OUT), OPTIONAL :: PerSpecies ! Max spc wildcard
    INTEGER,             INTENT(OUT), OPTIONAL :: Rank       ! # of dimensions
    INTEGER,             INTENT(OUT), OPTIONAL :: Type       ! Desc of data type
    INTEGER,             INTENT(OUT), OPTIONAL :: VLoc       ! Vert placement
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
    Found = .TRUE.

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
    IF ( isSpecies ) PerSpecies = ''       ! Assume not per species

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

       CASE ( 'JVALUES' )
          IF ( isDesc  )   Desc       = 'Photolysis rate' !TODO: append to this?
          IF ( isUnits )   Units      = 's-1'
          IF ( isRank  )   Rank       = 3
          IF ( isSpecies ) PerSpecies = 'ALL' ! TODO: fix species mapping

       CASE DEFAULT
          ! Need to add better error handling
          Found = .False.
          ErrMsg = 'WARNING: Metadata not found for State_Diag field: ' &
                   // TRIM( Name )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN

    END SELECT

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
  SUBROUTINE Register_DiagField_R4_2D( am_I_Root, metadataID, Ptr2Data, &
                                       State_Chm, State_Diag, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID      ! Name
    REAL*4,            POINTER       :: Ptr2Data(:,:)   ! pointer to data
    TYPE(ChmState),    INTENT(IN)    :: State_Chm       ! Obj for chem state
    TYPE(DgnState),    INTENT(IN)    :: State_Diag      ! Obj for diag state
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
    CHARACTER(LEN=255)     :: ErrMsg, ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    INTEGER                :: N, D, nSpecies
    INTEGER                :: rank, type, vloc
    LOGICAL                :: found
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_2D (in Headers/state_diag_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Diag field'

    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,   &
                                  desc=desc,   units=units, rank=rank,    &
                                  type=type,   vloc=vloc,                 &
                                  perSpecies=perSpecies                 )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Check that metadata consistent with data pointer
    IF ( ( ( perSpecies == '' ) .AND. ( rank /= 2 ) )  &
         .OR. ( ( perSpecies /= '' ) .AND. ( rank /= 1 ) ) ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
          
    IF ( perSpecies /= '' ) THEN

       ! TODO: add more cases as needed
       SELECT CASE ( perSpecies )
          CASE ( 'ALL' )
             nSpecies = State_Chm%nSpecies
          CASE ( 'ADV' )
             nSpecies = State_Chm%nAdvect
          CASE ( 'DRY' )
             nSpecies = State_Chm%nDryDep
          CASE ( 'WET' )
             nSpecies = State_Chm%nWetDep
          CASE DEFAULT
             ErrMsg = 'Handling of perSpecies ' // TRIM(perSpecies) // &
                      ' is not implemented for this combo of data type and size'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
       END SELECT

       DO N = 1, nSpecies          
          ! TODO: add more cases as needed
          SELECT CASE ( perSpecies )
             CASE ( 'ALL', 'ADV' )
                D = N
             CASE ( 'DRY' )
                D =  State_Chm%Map_DryDep(N)
             CASE ( 'WET' )
                D =  State_Chm%Map_WetDep(N)
          END SELECT
          SpcInfo  => State_Chm%SpcData(D)%Info
          thisSpcName = TRIM( metadataID ) // '__' // TRIM( SpcInfo%Name )
          thisSpcDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcName,          &
                                  Description  = thisSpcDesc,          &
                                  Units        = units,                &
                                  Data1d_4     = Ptr2Data(:,N),        &
                                  RC           = RC                   )
          SpcInfo => NULL()
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = ErrMsg_reg // ' where perSpecies is ' // TRIM(perSpecies)
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If not tied to species then simply add the single field
    ELSE
       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Diag%Registry,  &
                               State        = State_Diag%State,     &
                               Variable     = MetadataID,           &
                               Description  = desc,                 &
                               Units        = units,                &
                               Data2d_4     = Ptr2Data,             &
                               RC           = RC                   )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = ErrMsg_reg // ' where diagnostics is not tied to species'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
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
  SUBROUTINE Register_DiagField_R4_3D( am_I_Root, metadataID, Ptr2Data,  &
                                       State_Chm, State_Diag, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root       ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID      ! Name
    REAL*4,            POINTER       :: Ptr2Data(:,:,:) ! pointer to data
    TYPE(ChmState),    INTENT(IN)    :: State_Chm       ! Obj for chem state
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
    CHARACTER(LEN=255)     :: ErrMsg, ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    INTEGER                :: N, D, nSpecies
    INTEGER                :: rank, type,  vloc
    LOGICAL                :: found
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_3D (in Headers/state_diag_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Diag field'

    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,   &
                                  desc=desc,   units=units, rank=rank,    &
                                  type=type,   vloc=vloc,                 &
                                  perSpecies=perSpecies                 )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Check that metadata consistent with data pointer
    IF ( ( ( perSpecies == '' ) .AND. ( rank /= 3 ) )  &
         .OR. ( ( perSpecies /= '' ) .AND. ( rank /= 2 ) ) ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( perSpecies /= '' ) THEN

       ! TODO: add more cases as needed
       SELECT CASE ( perSpecies )
          CASE ( 'ALL' )
             nSpecies = State_Chm%nSpecies
          CASE ( 'ADV' )
             nSpecies = State_Chm%nAdvect
          CASE ( 'DRY' )
             nSpecies = State_Chm%nDryDep
          CASE ( 'WET' )
             nSpecies = State_Chm%nWetDep
          CASE DEFAULT
             ErrMsg = 'Handling of perSpecies ' // TRIM(perSpecies) // &
                      ' is not implemented for this combo of data type and size'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
       END SELECT

       DO N = 1, nSpecies          
          ! TODO: add more cases as needed
          SELECT CASE ( perSpecies )
             CASE ( 'ALL', 'ADV' )
                D = N
             CASE ( 'DRY' )
                D =  State_Chm%Map_DryDep(N)
             CASE ( 'WET' )
                D =  State_Chm%Map_WetDep(N)
          END SELECT
          SpcInfo  => State_Chm%SpcData(D)%Info
          thisSpcName = TRIM( metadataID ) // '__' // TRIM( SpcInfo%Name )
          thisSpcDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
          CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                                  Registry     = State_Diag%Registry,  &
                                  State        = State_Diag%State,     &
                                  Variable     = thisSpcName,          &
                                  Description  = thisSpcDesc,          &
                                  Units        = units,                &
                                  Data2d_4     = Ptr2Data(:,:,N),      &
                                  RC           = RC                   )
          SpcInfo => NULL()
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = ErrMsg_reg // ' where perSpecies is ' // TRIM(perSpecies)
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDDO

    ! If not tied to species then simply add the single field
    ELSE
       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Diag%Registry,  &
                               State        = State_Diag%State,     &
                               Variable     = metadataID,           &
                               Description  = desc,                 &
                               Units        = units,                &
                               Data3d_4     = Ptr2Data,             &
                               RC           = RC                   )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = ErrMsg_reg // ' where diagnostics is not tied to species'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
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
  SUBROUTINE Register_DiagField_R4_4D( am_I_Root, metadataID, Ptr2Data,  &
                                       State_Chm, State_Diag, RC )
!
! !USES:
!
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN)    :: am_I_Root         ! Root CPU?
    CHARACTER(LEN=*),  INTENT(IN)    :: metadataID        ! Name
    REAL*4,            POINTER       :: Ptr2Data(:,:,:,:) ! pointer to data
    TYPE(ChmState),    INTENT(IN)    :: State_Chm         ! Obj for chem state
    TYPE(DgnState),    INTENT(IN)    :: State_Diag        ! Obj for diag state
!
! !INPUT/OUTPUT PARAMETERS:
!
!
! !OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(OUT)   :: RC                ! Success/failure
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
    CHARACTER(LEN=255)     :: ErrMsg, ErrMsg_reg, ThisLoc
    CHARACTER(LEN=255)     :: desc, units, perSpecies
    CHARACTER(LEN=255)     :: thisSpcName, thisSpcDesc
    INTEGER                :: N, D, nSpecies
    INTEGER                :: rank, type,  vloc
    LOGICAL                :: found
    TYPE(Species), POINTER :: SpcInfo

    ! Initialize
    RC = GC_SUCCESS
    ThisLoc = ' -> at Register_DiagField_R4_4D (in Headers/state_diag_mod.F90)'
    ErrMsg_reg = 'Error encountered while registering State_Diag field'

    CALL Get_Metadata_State_Diag( am_I_Root,   metadataID,  Found,  RC,   &
                                  desc=desc,   units=units, rank=rank,    &
                                  type=type,   vloc=vloc,                 &
                                  perSpecies=perSpecies                 )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( ErrMsg_reg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Check that metadata consistent with data pointer
    IF ( rank /= 3 ) THEN
       ErrMsg = 'Data dims and metadata rank do not match for ' // &
                TRIM(metadataID)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
    
    ! Assume always tied to a species
    ! TODO: add more cases as needed
    SELECT CASE ( perSpecies )
       CASE ( 'ALL' )
          nSpecies = State_Chm%nSpecies
       CASE ( 'ADV' )
          nSpecies = State_Chm%nAdvect
       CASE ( 'DRY' )
          nSpecies = State_Chm%nDryDep
       CASE ( 'WET' )
          nSpecies = State_Chm%nWetDep
       CASE DEFAULT
          ErrMsg = 'Handling of perSpecies ' // TRIM(perSpecies) // &
                   ' is not implemented for this combo of data type and size'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
    END SELECT

    DO N = 1, nSpecies          
       ! TODO: add more cases as needed
       SELECT CASE ( perSpecies )
          CASE ( 'ALL', 'ADV' )
             D = N
          CASE ( 'DRY' )
             D =  State_Chm%Map_DryDep(N)
          CASE ( 'WET' )
             D =  State_Chm%Map_WetDep(N)
       END SELECT
       SpcInfo  => State_Chm%SpcData(D)%Info
       thisSpcName = TRIM( metadataID ) // '__' // TRIM( SpcInfo%Name )
       thisSpcDesc = TRIM( Desc ) // ' ' // TRIM( SpcInfo%Name )
       CALL Registry_AddField( am_I_Root    = am_I_Root,            &
                               Registry     = State_Diag%Registry,  &
                               State        = State_Diag%State,     &
                               Variable     = thisSpcName,          &
                               Description  = thisSpcDesc,          &
                               Units        = units,                &
                               Data3d_4     = Ptr2Data(:,:,:,N),    &
                               RC           = RC                   )
       SpcInfo => NULL()
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = ErrMsg_reg // ' where perSpecies is ' // TRIM(perSpecies)
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDDO

  END SUBROUTINE Register_DiagField_R4_4D
!EOC
END MODULE State_Diag_Mod
