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
  PUBLIC  :: Diagnostics_Init
  PUBLIC  :: Diagnostics_Write
  PUBLIC  :: Diagnostics_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: DiagInit_Drydep
  PRIVATE :: DiagInit_Tracer_Conc
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
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GRID_MOD,           ONLY : AREA_M2
    USE TIME_MOD,           ONLY : GET_TS_CHEM
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
    !
    ! NOTE: Right now there is only 1 GEOS-Chem collection, but eventually
    ! we may want to add more (i.e. hourly, instantaneous, monthly, etc.)
    !-----------------------------------------------------------------------

!    ! Rn/Pb/Be source diagnostic (ND01)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! Rn/Pb/Be decay diagnostic (ND02)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! Hg emissions/prod/loss diagnostic (ND03)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! Sulfate prod/loss diagnostic (ND05)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! Carbon aerosol sources diagnostic (ND07)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! Boundary layer fraction diagnostic (ND12)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! Cloud convection mass flux diagnostic (ND14)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! Boundary-layer mixing mass flux diagnostic (ND15)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! Areal fraction of precip diagnostic (ND16)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! Rainout fraction diagnostic (ND17)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! Washout fraction diagnostic (ND18)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! CH4 loss diagnostic (ND19)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! Optical depths diagnostic (ND21)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! Photolysis rates (J-values) diagnostic (ND22)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! E/W transport flux diagnostic (ND24)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! N/S transport flux diagnostic (ND25)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! U/D transport flux diagnostic (ND26)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! Strat influx (NOx, Ox, HNO3 diagnostic (ND27)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
!    ! Land map diagnostic (ND30)
!    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
!    ENDIF
!
    ! Drydep diagnostic (ND44)
    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
    ENDIF

    ! Tracer concentration diagnostics (ND45)
    CALL DIAGINIT_TRACER_CONC( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_TRACER_CONC', LOC ) 
    ENDIF

    ! Leave with success
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
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_drydep
!
! !DESCRIPTION: Subroutine DIAGINIT\_DRYDEP initializes the dry deposition 
! diagnostics arrays. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Drydep( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_DIAGN_MOD,      ONLY : Diagn_Create
    USE HCO_ERROR_MOD
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  13 Jan 2015 - C. Keller   - Initial version 
!  15 Jan 2015 - R. Yantosca - Init drydep velocity & flux diagnostics
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: cId,      Collection, D, N
    CHARACTER(LEN=15)  :: OutOper,  WriteFreq
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGINIT_DRYDEP (diagnostics_mod.F)' 
    
    !=======================================================================
    ! DIAGINIT_DRYDEP begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND44 diagnostic is turned off
    IF ( Input_Opt%ND44 <= 0 ) RETURN
      
    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND44_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND44_OUTPUT_FREQ
      
    ! Loop over # of depositing species
    DO D = 1, Input_Opt%NUMDEP
         
       ! Corresponding GEOS-Chem tracer number
       N = Input_Opt%NTRAIND(D)

       ! If this tracer number N is scheduled for output in input.geos, 
       ! then define the diagnostic containers for drydep velocity & flux
       IF ( ANY( Input_Opt%TINDEX(44,:) == N ) ) THEN

          !----------------------------------------------------------------
          ! Create containers for drydep velocity [m/s]
          !----------------------------------------------------------------

          ! Diagnostic name
          DiagnName = 'DEPVEL_' // TRIM( Input_Opt%DEPNAME(D) )

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cName     = TRIM( DiagnName ), &
                             AutoFill  = 0,                 &
                             ExtNr     = -1,                &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = -1,                &
                             SpaceDim  =  2,                &
                             LevIDx    = -1,                &
                             OutUnit   = 's-1',             &
                             OutOper   = TRIM( OutOper   ), &
                             WriteFreq = TRIM( WriteFreq ), &
                             cID       = cId,               &
                             RC        = RC )

          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
          
          !----------------------------------------------------------------
          ! Create containers for drydep flux [molec/cm2/s]
          !----------------------------------------------------------------

          ! Diagnostic name
          DiagnName = 'DEPFLUX_' // TRIM( Input_Opt%DEPNAME(D) )

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        &
                             cName     = TRIM( DiagnName ), &
                             AutoFill  = 0,                 &
                             ExtNr     = -1,                &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = -1,                &
                             SpaceDim  =  2,                &
                             LevIDx    = -1,                &
                             OutUnit   = 'cm-2 s-1',        &
                             OutOper   = TRIM( OutOper   ), &
                             WriteFreq = TRIM( WriteFreq ), &
                             cID       = cId,               &
                             RC        = RC )

          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP ( MSG, LOC ) 
          ENDIF
       ENDIF
    ENDDO

  END SUBROUTINE DiagInit_Drydep
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_tracer_conc
!
! !DESCRIPTION: Subroutine DIAGINIT\_TRACER\_CONC initializes the tracer
!  concentration diagnostic (aka ND45).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Tracer_Conc( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE Error_Mod,          ONLY : Error_Stop
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE HCO_Error_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jan 2015 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: cId,      Collection, N
    CHARACTER(LEN=15)  :: OutOper,  WriteFreq
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGINIT_TRACER_CONC (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_TRACER_CONC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND45 diagnostic is turned off
    IF ( Input_Opt%ND45 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND45_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND45_OUTPUT_FREQ
      
    ! Loop over # of depositing species
    DO N = 1, Input_Opt%N_TRACERS
         
       ! If this tracer number N is scheduled for output in input.geos, 
       ! then define the diagnostic containers for drydep velocity & flux
       IF ( ANY( Input_Opt%TINDEX(45,:) == N ) ) THEN

          !----------------------------------------------------------------
          ! Create containers for drydep velocity [m/s]
          !----------------------------------------------------------------

          ! Diagnostic name
          DiagnName = 'TRACER_CONC_' // TRIM( Input_Opt%TRACER_NAME(N) )

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cName     = TRIM( DiagnName ), &
                             AutoFill  = 0,                 &
                             ExtNr     = -1,                &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = -1,                &
                             SpaceDim  =  2,                &
                             LevIDx    = -1,                &
                             OutUnit   = 'v/v',             &
                             OutOper   = TRIM( OutOper   ), &
                             WriteFreq = TRIM( WriteFreq ), &
                             cId       = cId,               &
                             RC        = RC )

          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
       ENDIF
    ENDDO

  END SUBROUTINE DiagInit_Tracer_Conc
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_pb210_emiss
!
! !DESCRIPTION: Subroutine DIAGINIT\_PB210\_EMISS initializes the 210Pb emissions diagnostic (aka part of ND01).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Pb210_Emiss( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE Error_Mod,          ONLY : Error_Stop
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE HCO_Error_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  21 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: cId,      Collection, N
    CHARACTER(LEN=15)  :: OutOper,  WriteFreq
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGINIT_PB210_EMISS (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_TRACER_CONC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND01 diagnostic is turned off
    IF ( Input_Opt%ND01 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND01_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND01_OUTPUT_FREQ
        
    ! Corresponding GEOS-Chem tracer number
    N = ! WHAT TO PUT HERE??? Study TINDEX / tracer storage in Input_Mod
  
    ! If this tracer number is scheduled for output in input.geos, 
    ! then define the diagnostic container for 210Pb emissions.
    IF ( ANY ( Input_Opt%TINDEX(1,:) ) == N ) THEN

       !----------------------------------------------------------------
       ! Create containers for 210Pb emissions [kg/s]
       !----------------------------------------------------------------

       ! Diagnostic name
       DiagnName = 'TRACER_CONC_' // TRIM( Input_Opt%TRACER_NAME(N) )

       ! Create container
       CALL Diagn_Create( am_I_Root,                     &
                          Col       = Collection,        & 
                          cName     = TRIM( DiagnName ), &
                          AutoFill  = 0,                 &
                          ExtNr     = -1,                &
                          Cat       = -1,                &
                          Hier      = -1,                &
                          HcoID     = -1,                &
                          SpaceDim  =  2,                &
                          LevIDx    = -1,                &
                          OutUnit   = 'kg/s',             &
                          OutOper   = TRIM( OutOper   ), &
                          WriteFreq = TRIM( WriteFreq ), &
                          cId       = cId,               &
                          RC        = RC )

       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF
    ENDIF

  END SUBROUTINE DiagInit_Tracer_Conc
!EOC
END MODULE Diagnostics_Mod
#endif
