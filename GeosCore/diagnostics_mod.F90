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
! !PRIVATE MEMBER FUNCTIONS:
!
!
! !MODULE PARAMETER
!
  ! Prefix of restart file. This file will hold all diagnostics that are 
  ! written out at the end of a simulation (either because their output 
  ! frequency is set to 'End' or the run finishes and these diagnostics
  ! haven't reached the end of their output interval yet).
  CHARACTER(LEN= 31), PARAMETER           :: RST = 'GEOSCHEM_Restart'
  CHARACTER(LEN= 31), PARAMETER           :: DGN = 'GEOSCHEM_Diagnostics'
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

    ! Define diagnostics in these modules:
    USE DRYDEP_MOD,         ONLY : DIAGINIT_DRYDEP
    USE UCX_MOD,            ONLY : DIAGINIT_UCX
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(INOUT)  :: Input_Opt  ! Input opts
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
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

    ! ----------------------------------------------------------------------
    ! Create collection
    ! ----------------------------------------------------------------------
    CALL DiagnCollection_Create( am_I_Root,            &
                                 NX        = IIPAR,    &
                                 NY        = JJPAR,    &
                                 NZ        = LLPAR,    &
                                 TS        = TS,       & 
                                 AM2       = AM2,      & 
                                 PREFIX    = DGN,      &
                                 COL       = GCDiagNr, &
                                 OVERWRITE = .FALSE.,  & 
                                 RC        = RC         )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'Cannot overwrite collection', LOC ) 
    ENDIF

    ! Cleanup
    AM2 => NULL()

    ! ----------------------------------------------------------------------
    ! Add diagnostics to collection 
    ! ----------------------------------------------------------------------
    IF ( Input_Opt%LDRYD ) CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
    IF ( Input_Opt%LUCX  ) CALL DIAGINIT_UCX   ( am_I_Root, Input_Opt, RC )

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
! to disk. If the variable RESTART is set to true, all GEOS-Chem diagnostics
! are passed to the restart file.
!\\
!\\
! The Diagnostics\_Write routine is called from main.F, at the end of the time
! loop and during cleanup (to write the restart file).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagnostics_Write ( am_I_Root, RESTART, RC ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE HCO_STATE_MOD,      ONLY : HCO_STATE
    USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoState, SetHcoTime
    USE HCOIO_Diagn_Mod,    ONLY : HCOIO_Diagn_WriteOut
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
    LOGICAL,          INTENT(IN   )  :: RESTART    ! Write restart file? 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  09 Jan 2015 - C. Keller   - Initial version 
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

    ! RESTART: write out all diagnostics. Use current time stamp and save into
    ! restart file.
    IF ( RESTART ) THEN
       CALL HCOIO_DIAGN_WRITEOUT ( am_I_Root, HcoState, WriteAll=.TRUE., RC=RC, &
                                   UsePrevTime=.FALSE., PREFIX=RST, COL=GCDiagNr )
       IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( 'Restart write error', LOC ) 

    ! Not restart: write out regular diagnostics. Use current time stamp.
    ELSE
       CALL HCOIO_DIAGN_WRITEOUT ( am_I_Root, HcoState, WriteAll=.FALSE., RC=RC, &
                                   UsePrevTime=.FALSE., COL=GCDiagNr )
       IF ( RC /= HCO_SUCCESS ) CALL ERROR_STOP( 'Diagnostics write error', LOC ) 
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
  SUBROUTINE Diagnostics_Final ( ) 
!
! !REVISION HISTORY: 
!  09 Jan 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC

    CALL DiagnCollection_Cleanup( COL = GCDiagNr )

  END SUBROUTINE Diagnostics_Final
!EOC
END MODULE Diagnostics_Mod
