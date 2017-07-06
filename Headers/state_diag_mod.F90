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
  USE Precision_Mod
  USE Registry_Mod, ONLY : MetaRegItem

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Init_State_Diag
  PUBLIC :: Cleanup_State_Diag
  PUBLIC :: Lookup_State_Diag
  PUBLIC :: Print_State_Diag
!
! !PUBLIC DATA MEMBERS:
!
  !=========================================================================
  ! Derived type for Meteorology State
  !=========================================================================
  TYPE, PUBLIC :: DgnState

     !----------------------------------------------------------------------
     ! Diagnostic arrays
     !----------------------------------------------------------------------
     REAL(f4),          POINTER :: F_Adv_EW(:,:,:,:)
     REAL(f4),          POINTER :: F_Adv_NS(:,:,:,:)
     REAL(f4),          POINTER :: F_Adv_UP(:,:,:,:)

     !----------------------------------------------------------------------
     ! Registry of variables contained within State_Met
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
!EOC
!------------------------------------------------------------------------------
!BOC
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
  SUBROUTINE Init_State_Diag( am_I_Root, IM, JM, LM, Input_Opt, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod, ONLY : OptInput
    USe Registry_Mod,  ONLY : Registry_AddField
!
! !INPUT PARAMETERS:
! 
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
    INTEGER,        INTENT(IN)    :: IM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: JM          ! # longitudes on this PET
    INTEGER,        INTENT(IN)    :: LM          ! # longitudes on this PET
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, Desc, Units

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      =  GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Init_State_Diag (in Headers/state_diag_mod.F90)'
    Desc    = ''
    Units   = ''

    !=======================================================================
    ! Allocate module variables
    !=======================================================================

    ! For now, just nullify
    State_Diag%F_Adv_EW => NULL()
    State_Diag%F_Adv_NS => NULL()
    State_Diag%F_Adv_UP => NULL()
    
  END SUBROUTINE Init_State_Diag
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_State_Met
!
! !DESCRIPTION: Subroutine CLEANUP\_STATE\_MET deallocates all fields 
!  of the meteorology state object.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_State_Diag( am_I_Root, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
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
!  contained within the State_Diag object.  This is basically a wrapper for
!  routine REGISTRY_PRINT in registry_mod.F90.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Print_State_Diag( am_I_Root, State_Diag, RC, ShortFormat )
!
! !USES:
!
    USE ErrCode_Mod
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
    PRINT*
    PRINT*, 'Registered variables contained within the State_Diag object:'
    PRINT*, REPEAT( '=', 79 )

    ! Print registry info in truncated format
    CALL Registry_Print( am_I_Root   = am_I_Root,            &
                         Registry    = State_Diag%Registry,  &
                         ShortFormat = ShortFormat,          &
                         RC          = RC                   )

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
!  variable contained within the State_Diag object by searching for its name.
!  This is basically a wrapper for routine REGISTRY_LOOKUP in registry_mod.F90.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Lookup_State_Diag( am_I_Root,  State_Diag,  Variable,  & 
                                RC,         Description, KindVal,   & 
                                MemoryInKb, Rank,        Units,     &    
                                Ptr2d_4,    Ptr3d_4,     Ptr2d_I,   &
                                Ptr3d_I                            )
!
! !USES:
!
    USE ErrCode_Mod
    USE Registry_Mod, ONLY : Registry_Lookup
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
    INTEGER,             OPTIONAL :: KindVal         ! Numerical KIND value
    REAL(fp),            OPTIONAL :: MemoryInKb      ! Memory usage
    INTEGER,             OPTIONAL :: Rank            ! Size of data
    CHARACTER(LEN=255),  OPTIONAL :: Units           ! Units of data

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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! Initialize
    !=======================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Lookup_State_Diag (in Headers/state_met_mod.F90)'

    !=======================================================================
    ! Look up a variable; Return metadata and/or a pointer to the data
    !=======================================================================
    CALL Registry_Lookup( am_I_Root   = am_I_Root,            &
                          Registry    = State_Diag%Registry,  &
                          State       = State_Diag%State,     &
                          Variable    = Variable,             &
                          Description = Description,          &
                          KindVal     = KindVal,              &
                          MemoryInKb  = MemoryInKb,           &
                          Rank        = Rank,                 &
                          Units       = Units,                &
                          Ptr2d_4     = Ptr2d_4,              &
                          Ptr3d_4     = Ptr3d_4,              &
                          Ptr2d_I     = Ptr2d_I,              &
                          Ptr3d_I     = Ptr3d_I,              &
                          RC          = RC                   )

    ! Trap error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Could not find variable ' // TRIM( Variable ) // &
               ' in the State_Diag registry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Lookup_State_Diag
!EOC
END MODULE State_Diag_Mod
