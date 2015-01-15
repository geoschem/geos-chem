#if defined( DEVEL )
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diagnostics_mod.F90
!
! !DESCRIPTION: Module diagnostics\_mod.F90 is the first crack of the new
! diagnostics package for GEOS-Chem.
!
! !INTERFACE:
!
MODULE Diagnostics_Mod
!
! !USES:
!
  USE Precision_Mod
  USE HCO_Error_Mod
  USE HCO_Diagn_Mod
 
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Diagnostics_Init
  PUBLIC :: Diagnostics_Write
  PUBLIC :: Diagnostics_Final
!
! !DEFINED PARAMETERS:
!
  ! Prefix of restart file. This file will hold all diagnostics that are 
  ! written out at the end of a simulation (either because their output 
  ! frequency is set to 'End' or the run finishes and these diagnostics
  ! haven't reached the end of their output interval yet).
  CHARACTER(LEN=31), PARAMETER :: RST = 'GEOSCHEM_Restart'
  CHARACTER(LEN=31), PARAMETER :: DGN = 'GEOSCHEM_Diagnostics'
!
! !REVISION HISTORY:
!  09 Jan 2015 - C. Keller   - Initial version. 
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagnostics_Init
!
! !DESCRIPTION: Subroutine Diagnostics\_Init initializes the GEOS-Chem 
! diagnostics collection. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagnostics_Init( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
    USE GRID_MOD,           ONLY : AREA_M2
    USE TIME_MOD,           ONLY : GET_TS_CHEM
    USE DRYDEP_MOD,         ONLY : DIAGINIT_DRYDEP
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(MetState),   INTENT(IN   ) :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   ) :: State_Chm  ! Chemistry state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT) :: Input_Opt  ! Input opts
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  09 Jan 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(sp)           :: TS
    REAL(fp), POINTER  :: AM2(:,:) => NULL()
    CHARACTER(LEN=255) :: LOC = 'Diagnostics_Init (diagnostics_mod.F90)'

    !=======================================================================
    ! Diagnostics_Init begins here 
    !=======================================================================

    ! Define collection variables
    AM2    => AREA_M2(:,:,1)
    TS     =  GET_TS_CHEM() * 60.0_sp

    !-----------------------------------------------------------------------
    ! Create diagnostics collection for GEOS-Chem.  This will keep the
    ! GEOS-Chem diagostics separate from the HEMCO diagnostics.
    !-----------------------------------------------------------------------
    CALL DiagnCollection_Create( am_I_Root,                             &
                                 NX        = IIPAR,                     &
                                 NY        = JJPAR,                     &
                                 NZ        = LLPAR,                     &
                                 TS        = TS,                        &
                                 AM2       = AM2,                       &
                                 PREFIX    = DGN,                       &
                                 COL       = Input_Opt%DIAG_COLLECTION, &
                                 OVERWRITE = .TRUE.,                    & 
                                 RC        = RC         )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'Cannot overwrite collection', LOC ) 
    ENDIF

    ! Cleanup
    AM2 => NULL()

    !-----------------------------------------------------------------------
    ! Add diagnostics to collection 
    ! (Add calls to additional subroutines for other diagnostics here!)
    !-----------------------------------------------------------------------

    ! Drydep diagnostic (ND44)
    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )

    ! Leave w/ success
    RC = GIGC_SUCCESS

  END SUBROUTINE Diagnostics_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagnostics_Write
!
! !DESCRIPTION: Subroutine Diagnostics\_Write writes the GEOS-Chem diagnostics
! to disk. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagnostics_Write( am_I_Root, Input_Opt, LAST, RC ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCO_STATE_MOD,      ONLY : HCO_STATE
    USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoState, SetHcoTime
    USE HCOIO_Diagn_Mod,    ONLY : HCOIO_Diagn_WriteOut
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN ) :: Input_Opt  ! Input Options object
    LOGICAL,        INTENT(IN ) :: LAST       ! Is this the last call? 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  09 Jan 2015 - C. Keller   - Initial version
!  15 Jan 2015 - R. Yantosca - Now accept Input_Opt via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(HCO_STATE), POINTER :: HcoState => NULL()
    CHARACTER(LEN=255)       :: LOC = 'Diagnostics_Write (diagnostics_mod.F90)'

    !=======================================================================
    ! Diagnostics_Write begins here 
    !=======================================================================

    ! Make sure HEMCO time is in sync
    CALL SetHcoTime( am_I_Root, .FALSE., RC )
    IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( 'Cannot update HEMCO time', LOC ) 

    ! Get pointer to HEMCO state object.
    CALL GetHcoState( HcoState )
    IF ( .NOT. ASSOCIATED(HcoState) ) THEN
       CALL ERROR_STOP( 'Cannot get HEMCO state object', LOC )
    ENDIF

    !-----------------------------------------------------------------------
    ! Write diagnostics to the diagnostics file
    !-----------------------------------------------------------------------
    CALL HCOIO_DIAGN_WRITEOUT( am_I_Root,                                &
                               HcoState,                                 &
                               WriteAll    = .FALSE.,                    &
                               UsePrevTime = .FALSE.,                    & 
                               COL         = Input_Opt%DIAG_COLLECTION,  & 
                               RC          = RC                         )

    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'Diagnostics write error', LOC ) 
    ENDIF

    !-----------------------------------------------------------------------
    ! Last call: write out all diagnostics. 
    ! Use current time stamp and write into restart file
    !-----------------------------------------------------------------------
    IF ( LAST ) THEN
       CALL HCOIO_DIAGN_WRITEOUT( am_I_Root,                                & 
                                  HcoState,                                 &
                                  WriteAll    = .TRUE.,                     &
                                  UsePrevTime = .FALSE.,                    &
                                  PREFIX      = RST,                        &  
                                  COL         = Input_Opt%DIAG_COLLECTION,  &
                                  RC          = RC                         )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'Diagnostics write error at end of run', LOC ) 
       ENDIF

    ENDIF

    ! Free pointer
    HcoState => NULL()

    ! Leave w/ success
    RC = GIGC_SUCCESS

  END SUBROUTINE Diagnostics_Write
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diagnostics_Final
!
! !DESCRIPTION: Subroutine Diagnostics\_Final finalizes the GEOS-Chem 
! diagnostics collection. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagnostics_Final( am_I_Root, Input_Opt, RC ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN ) :: Input_Opt  ! Input Options objec
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  09 Jan 2015 - C. Keller   - Initial version 
!  15 Jan 2015 - R. Yantosca - Now accept Input_Opt, am_I_Root, RC arguments
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Finalize diagnostics
    CALL DiagnCollection_Cleanup( COL = Input_Opt%DIAG_COLLECTION )

    ! Return with success
    RC = GIGC_SUCCESS

  END SUBROUTINE Diagnostics_Final
!EOC
END MODULE Diagnostics_Mod
#endif
