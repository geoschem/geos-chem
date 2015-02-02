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
  PRIVATE :: DiagInit_Pb_Emiss


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

    ! Pb emissions diagnostic (ND01)
    CALL DIAGINIT_PB_EMISS( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_PB_EMISS', LOC ) 
    ENDIF

    ! Rn/Pb/Be decay diagnostic (ND02)
    CALL DIAGINIT_RN_DECAY( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_RN_DECAY', LOC ) 
    ENDIF

!    ! Hg emissions/prod/loss diagnostic (ND03)
!    CALL DIAGINIT_HG_SOURCE( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_HG_SOURCE', LOC ) 
!    ENDIF
!
!    ! Sulfate prod/loss diagnostic (ND05)
!    CALL DIAGINIT_SULFATE_PL( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_SULFATE_PL', LOC ) 
!    ENDIF
!
!    ! Carbon aerosol sources diagnostic (ND07)
!    CALL DIAGINIT_C_AERSRC( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_C_AERSRC', LOC ) 
!    ENDIF
!
    ! Boundary layer fraction diagnostic (ND12)
    CALL DIAGINIT_BL_FRAC( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_BL_FRAC', LOC ) 
    ENDIF

    ! Cloud convection mass flux diagnostic (ND14)
    CALL DIAGINIT_CLDCONV_FLX( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_CLDCONV_FLX', LOC ) 
    ENDIF

    ! Boundary-layer mixing mass flux diagnostic (ND15)
    CALL DIAGINIT_BLMIX_FLX( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
    ENDIF

    ! Areal fraction of precip diagnostic (ND16)
    CALL DIAGINIT_PRECIP_FRAC( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_PRECIP_FRAC', LOC ) 
    ENDIF

    ! Rainout fraction diagnostic (ND17)
    CALL DIAGINIT_RAIN_FRAC( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_RAIN_FRAC', LOC ) 
    ENDIF

    ! Washout fraction diagnostic (ND18)
    CALL DIAGINIT_WASH_FRAC( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_WASH_FRAC', LOC ) 
    ENDIF

    ! CH4 loss diagnostic (ND19)
    CALL DIAGINIT_CH4_LOSS( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_CH4_LOSS', LOC ) 
    ENDIF

!    ! Optical depths diagnostic (ND21)
!    CALL DIAGINIT_CLD_OD( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_CLD_OD', LOC ) 
!    ENDIF
!
!    ! Photolysis rates (J-values) diagnostic (ND22)
!    CALL DIAGINIT_PHOTOLYSIS( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_PHOTOLYSIS', LOC ) 
!    ENDIF
!
    ! E/W transport flux diagnostic (ND24)
    CALL DIAGINIT_EW_FLX( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_EW_FLUX', LOC ) 
    ENDIF

    ! N/S transport flux diagnostic (ND25)
    CALL DIAGINIT_NS_FLX( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_NS_FLX', LOC ) 
    ENDIF

    ! U/D transport flux diagnostic (ND26)
    CALL DIAGINIT_VERT_FLX( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_VERT_FLX', LOC ) 
    ENDIF

    ! Strat influx (NOx, Ox, HNO3 diagnostic (ND27)
    CALL DIAGINIT_STRAT_FLX( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_STRAT_FLX', LOC ) 
    ENDIF

    ! Land map diagnostic (ND30)
    CALL DIAGINIT_LANDMAP( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_LANDMAP', LOC ) 
    ENDIF

!    ! Surface Pressure diagnostic (ND31)
!    CALL DIAGINIT_PRESSURE( am_I_Root, Input_Opt, RC )
!    IF ( RC /= GIGC_SUCCESS ) THEN
!       CALL ERROR_STOP( 'Error in DIAGINIT_PRESSURE', LOC )
!    ENDIF

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
    ! Tropopause Height diagnostic (ND55)
    CALL DIAGINIT_TROP_HEIGHT( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_TROP_HEIGHT', LOC )
    ENDIF

    ! Lightning Flashes diagnostic (ND56)
    CALL DIAGINIT_LIGHTNING( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_LIGHTNING', LOC )
    ENDIF

    ! Potential Temperature (THETA) diagnostic (ND57)
    CALL DIAGINIT_THETA( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_THETA', LOC )
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
    ! Write out diagnostics to file using current time stamp.
    ! If last, save into restart file. Else, write out regular diagnostics.
    !-----------------------------------------------------------------------

    ! DEBUGGING - LL 2/2/15
    PRINT *, 'Attempting to write diagnostics to netcdf'
    ! End of debugging

    IF ( LAST ) THEN
       CALL HCOIO_DIAGN_WRITEOUT( am_I_Root,                                & 
                                  HcoState,                                 &
                                  WriteAll    = .TRUE.,                     &
                                  UsePrevTime = .FALSE.,                    &
                                  PREFIX      = RST,                        &   
                                  COL         = Input_Opt%DIAG_COLLECTION,  &
                                  RC          = RC                         )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'Restart write error', LOC )
       ENDIF
    ELSE
       CALL HCOIO_DIAGN_WRITEOUT( am_I_Root,                                &
                                  HcoState,                                 &
                                  WriteAll    = .FALSE.,                    &
                                  UsePrevTime = .FALSE.,                    & 
                                  COL         = Input_Opt%DIAG_COLLECTION,  & 
                                  RC          = RC                         )

       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'Diagnostics write error', LOC ) 
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
                             SpaceDim  =  3,                &
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
! !IROUTINE: diaginit_pb_emiss
!
! !DESCRIPTION: Subroutine DIAGINIT\_PB\_EMISS initializes the Pb emissions 
!  diagnostic (aka ND01). Other ND01 tracer emissions (Rn and Be7) are 
!  handled within HEMCO.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Pb_Emiss( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE Error_Mod,          ONLY : Error_Stop
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE HCO_Error_Mod
    USE TRACERID_MOD,       ONLY : IDTPB   
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
    INTEGER            :: cId,      Collection
    CHARACTER(LEN=15)  :: OutOper,  WriteFreq
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGINIT_PB_EMISS (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_PB_EMISS begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND01 diagnostic is turned off
    IF ( Input_Opt%ND01 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND01_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND01_OUTPUT_FREQ
    
    ! If the tracer number for lead is scheduled for output in input.geos, 
    ! then define the diagnostic container for 210Pb emissions.
    IF ( ANY ( Input_Opt%TINDEX(1,:) == IDTPB ) ) THEN

       !----------------------------------------------------------------
       ! Create containers for Pb emissions [kg/s]
       !----------------------------------------------------------------

       ! Diagnostic name
       DiagnName = 'EMISS_PB'

       ! Create container
       CALL Diagn_Create( am_I_Root,                     &
                          Col       = Collection,        & 
                          cName     = TRIM( DiagnName ), &
                          AutoFill  = 0,                 &
                          ExtNr     = -1,                &
                          Cat       = -1,                &
                          Hier      = -1,                &
                          HcoID     = -1,                &
                          SpaceDim  =  3,                &
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

  END SUBROUTINE DiagInit_Pb_Emiss
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_rn_decay
!
! !DESCRIPTION: Subroutine DIAGINIT\_RN\_DECAY initializes the Rn/Pb/Be7
!  decay diagnostic (aka ND02).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Rn_Decay( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE Error_Mod,          ONLY : Error_Stop
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE HCO_Error_Mod
    USE TRACERID_MOD,       ONLY : IDTPb, IDTRn, IDTBe7   
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
!  23 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId, Collection, M
    INTEGER, PARAMETER   :: NumTracers = 3        ! Does a var exist for this?
    INTEGER              :: TracersN ( NumTracers ) 
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_RN_DECAY (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_RN_DECAY begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND02 diagnostic is turned off
    IF ( Input_Opt%ND02 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND02_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND02_OUTPUT_FREQ
    
    ! Assign array of tracer numbers corresponding to decaying species
    TracersN = (/ IDTPb, IDTRN, IDTBe7 /)

    ! Loop over # of radon decay diagnostics
    DO M = 1, NumTracers

       ! If the tracer number is scheduled for output in input.geos, 
       ! then define the diagnostic container for that tracer.
       IF ( ANY ( Input_Opt%TINDEX(2,:) == TracersN( M ) ) ) THEN

          !----------------------------------------------------------------
          ! Create containers for Rn/Pb/Be7 decay [kg/s]
          !----------------------------------------------------------------

          ! Diagnostic name
          DiagnName = 'DECAY_' //                           &
                      TRIM( Input_Opt%TRACER_NAME( TracersN( M ) ) )

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cName     = TRIM( DiagnName ), &
                             AutoFill  = 0,                 &
                             ExtNr     = -1,                &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = -1,                &
                             SpaceDim  =  3,                &
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
    ENDDO

  END SUBROUTINE DiagInit_Rn_Decay
!EOC

!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: diaginit_hg_source (LL in progress)
!!
!! !DESCRIPTION: Subroutine DIAGINIT\_HG\_SOURCE initializes the mercury
!!  emissions, production, and loss diagnostic (aka ND03).
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE DiagInit_Hg_Source( am_I_Root, Input_Opt, RC )
!!
!! !USES:
!!
!    USE Error_Mod,          ONLY : Error_Stop
!    USE GIGC_ErrCode_Mod
!    USE GIGC_Input_Opt_Mod, ONLY : OptInput
!    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
!    USE HCO_Error_Mod
!    USE TRACERID_MOD,        ONLY : ??? 
!!
!! !INPUT PARAMETERS:
!!
!    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
!    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
!! 
!! !REVISION HISTORY: 
!!  23 Jan 2015 - E. Lundgren - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    INTEGER              :: cId,      Collection, N
!!    INTEGER, PARAMETER   :: NUMTRC = 3
!!    INTEGER              :: TRCN ( NUMTRC )
!    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
!    CHARACTER(LEN=60)    :: DiagnName
!    CHARACTER(LEN=255)   :: MSG
!    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_HG_SOURCE (diagnostics_mod.F90)' 
!
!    !=======================================================================
!    ! DIAGINIT_HG_SOURCE begins here!
!    !=======================================================================
!      
!    ! Assume successful return
!    RC = GIGC_SUCCESS
!
!    ! Skip if ND03 diagnostic is turned off
!    IF ( Input_Opt%ND03 <= 0 ) RETURN
!
!    ! Get diagnostic parameters from the Input_Opt object
!    Collection = Input_Opt%DIAG_COLLECTION
!    OutOper    = Input_Opt%ND03_OUTPUT_TYPE
!    WriteFreq  = Input_Opt%ND03_OUTPUT_FREQ
!    
!    ! Assign array of tracer numbers corresponding to decaying species
!    TRCN(1) = IDTPB
!    TRCN(2) = IDTRN
!    TRCN(3) = IDTBE7
!
!    ! Loop over # of mercury diagnostics
!    DO N = 1, NumDiagND03
!
!       ! If the tracer number is scheduled for output in input.geos, 
!       ! then define the diagnostic container for that tracer.
!       IF ( ANY ( Input_Opt%TINDEX(3,:) == TRCN(N) ) ) THEN
!
!          !----------------------------------------------------------------
!          ! Create containers for mercury sources [units?]
!          !----------------------------------------------------------------
!
!          ! Diagnostic name
!          DiagnName = 'HG_SOURCE_' // TRIM( Input_Opt%TRACER_NAME(N) )
!
!          ! Create container
!          CALL Diagn_Create( am_I_Root,                     &
!                             Col       = Collection,        & 
!                             cName     = TRIM( DiagnName ), &
!                             AutoFill  = 0,                 &
!                             ExtNr     = -1,                &
!                             Cat       = -1,                &
!                             Hier      = -1,                &
!                             HcoID     = -1,                &
!                             SpaceDim  =  2,                &
!                             LevIDx    = -1,                &
!                             OutUnit   = 'kg' ,             &
!                             OutOper   = TRIM( OutOper   ), &
!                             WriteFreq = TRIM( WriteFreq ), &
!                             cId       = cId,               &
!                             RC        = RC )
!
!          IF ( RC /= HCO_SUCCESS ) THEN
!             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
!             CALL ERROR_STOP( MSG, LOC ) 
!          ENDIF
!       ENDIF
!    ENDDO
!
!  END SUBROUTINE DiagInit_Hg_Source
!!EOC

!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: diaginit_sulfate_pl (LL in progress)
!!
!! !DESCRIPTION: Subroutine DIAGINIT\_SULFATE\_PL initializes the sulfate
!!  chemistry production and loss diagnostic (aka ND05).
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE DiagInit_Sulfate_PL( am_I_Root, Input_Opt, RC )
!!
!! !USES:
!!
!    USE Error_Mod,          ONLY : Error_Stop
!    USE GIGC_ErrCode_Mod
!    USE GIGC_Input_Opt_Mod, ONLY : OptInput
!    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
!    USE HCO_Error_Mod
!!
!! !INPUT PARAMETERS:
!!
!    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
!    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
!! 
!! !REVISION HISTORY: 
!!  26 Jan 2015 - E. Lundgren - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    INTEGER              :: cId,      Collection, N
!    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
!    CHARACTER(LEN=60)    :: DiagnName
!    CHARACTER(LEN=255)   :: MSG
!    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_SULFATE_PL (diagnostics_mod.F90)' 
!
!    !=======================================================================
!    ! DIAGINIT_SULFATE_PL begins here!
!    !=======================================================================
!      
!    ! Assume successful return
!    RC = GIGC_SUCCESS
!
!    ! Skip if ND05 diagnostic is turned off
!    IF ( Input_Opt%ND05 <= 0 ) RETURN
!
!    ! Get diagnostic parameters from the Input_Opt object
!    Collection = Input_Opt%DIAG_COLLECTION
!    OutOper    = Input_Opt%ND05_OUTPUT_TYPE
!    WriteFreq  = Input_Opt%ND05_OUTPUT_FREQ
!    
!    ! Loop over # of diagnostics (only ones with units kg S)
!    DO N = 1, NumDiagND05
!
!       ! If the tracer number is scheduled for output in input.geos, 
!       ! then define the diagnostic container for that tracer.
!       IF ( ANY ( Input_Opt%TINDEX(5,:) == TRCN(N) ) ) THEN
!
!          !----------------------------------------------------------------
!          ! Create containers for sulfate PL [kg S] (8 of these so need a loop)
!          !----------------------------------------------------------------
!
!          ! Diagnostic name
!          DiagnName = 'PL_SUL_' // TRIM( Input_Opt%TRACER_NAME(N) )
!
!          ! Create container
!          CALL Diagn_Create( am_I_Root,                     &
!                             Col       = Collection,        & 
!                             cName     = TRIM( DiagnName ), &
!                             AutoFill  = 0,                 &
!                             ExtNr     = -1,                &
!                             Cat       = -1,                &
!                             Hier      = -1,                &
!                             HcoID     = -1,                &
!                             SpaceDim  =  3,                &
!                             LevIDx    = -1,                &
!                             OutUnit   = 'kg S' ,           &
!                             OutOper   = TRIM( OutOper   ), &
!                             WriteFreq = TRIM( WriteFreq ), &
!                             cId       = cId,               &
!                             RC        = RC )
!
!          IF ( RC /= HCO_SUCCESS ) THEN
!             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
!             CALL ERROR_STOP( MSG, LOC ) 
!          ENDIF
!
!          !----------------------------------------------------------------
!          ! Create containers for L(OH) by DMS [kg OH] - just 1
!          !----------------------------------------------------------------
!
!          ! Diagnostic name
!          DiagnName = 'PL_SUL_' // TRIM( Input_Opt%TRACER_NAME(N) )
!
!          ! Create container
!          CALL Diagn_Create( am_I_Root,                     &
!                             Col       = Collection,        & 
!                             cName     = TRIM( DiagnName ), &
!                             AutoFill  = 0,                 &
!                             ExtNr     = -1,                &
!                             Cat       = -1,                &
!                             Hier      = -1,                &
!                             HcoID     = -1,                &
!                             SpaceDim  =  3,                &
!                             LevIDx    = -1,                &
!                             OutUnit   = 'kg OH' ,          &
!                             OutOper   = TRIM( OutOper   ), &
!                             WriteFreq = TRIM( WriteFreq ), &
!                             cId       = cId,               &
!                             RC        = RC )
!
!          IF ( RC /= HCO_SUCCESS ) THEN
!             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
!             CALL ERROR_STOP( MSG, LOC ) 
!
!          !----------------------------------------------------------------
!          ! Create containers for L(NO3) by DMS [kg NO3] - just 1
!          !----------------------------------------------------------------
!
!          ! Diagnostic name
!          DiagnName = 'PL_SUL_' // TRIM( Input_Opt%TRACER_NAME(N) )
!
!          ! Create container
!          CALL Diagn_Create( am_I_Root,                     &
!                             Col       = Collection,        & 
!                             cName     = TRIM( DiagnName ), &
!                             AutoFill  = 0,                 &
!                             ExtNr     = -1,                &
!                             Cat       = -1,                &
!                             Hier      = -1,                &
!                             HcoID     = -1,                &
!                             SpaceDim  =  3,                &
!                             LevIDx    = -1,                &
!                             OutUnit   = 'kg NO3' ,         &
!                             OutOper   = TRIM( OutOper   ), &
!                             WriteFreq = TRIM( WriteFreq ), &
!                             cId       = cId,               &
!                             RC        = RC )
!
!          IF ( RC /= HCO_SUCCESS ) THEN
!             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
!             CALL ERROR_STOP( MSG, LOC ) 
!
!       ENDIF
!    ENDDO
!!
!  END SUBROUTINE DiagInit_Sulfate_PL
!!EOC

!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: diaginit_C_AerSrc (LL in progress)
!!
!! !DESCRIPTION: Subroutine DIAGINIT\_C\_AERSRC initializes the sources of
!!  black carbon and organic carbon diagnostic (aka ND07).
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE DiagInit_C_AerSrc( am_I_Root, Input_Opt, RC )
!!
!! !USES:
!!
!    USE Error_Mod,          ONLY : Error_Stop
!    USE GIGC_ErrCode_Mod
!    USE GIGC_Input_Opt_Mod, ONLY : OptInput
!    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
!    USE HCO_Error_Mod
!!
!! !INPUT PARAMETERS:
!!
!    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
!    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
!! 
!! !REVISION HISTORY: 
!!  26 Jan 2015 - E. Lundgren - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    INTEGER              :: cId,      Collection, N
!    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
!    CHARACTER(LEN=60)    :: DiagnName
!    CHARACTER(LEN=255)   :: MSG
!    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_C_AERSRC (diagnostics_mod.F90)' 
!
!    !=======================================================================
!    ! DIAGINIT_C_AERSRC begins here!
!    !=======================================================================
!      
!    ! Assume successful return
!    RC = GIGC_SUCCESS
!
!    ! Skip if ND07 diagnostic is turned off
!    IF ( Input_Opt%ND07 <= 0 ) RETURN
!
!    ! Get diagnostic parameters from the Input_Opt object
!    Collection = Input_Opt%DIAG_COLLECTION
!    OutOper    = Input_Opt%ND07_OUTPUT_TYPE
!    WriteFreq  = Input_Opt%ND07_OUTPUT_FREQ
!    
!    ! Loop over # of diagnostics (NOTE: diags 10-23 only if SOA simulation)
!    DO N = 1, NumDiagND07
!
!       ! If the tracer number is scheduled for output in input.geos, 
!       ! then define the diagnostic container for that tracer.
!       IF ( ANY ( Input_Opt%TINDEX(7,:) == TRCN(N) ) ) THEN
!
!          !----------------------------------------------------------------
!          ! Create containers for carbon aerosol sources [kg]
!          !----------------------------------------------------------------
!
!          ! Diagnostic name
!          DiagnName = 'C_AERSRC_' // TRIM( Input_Opt%TRACER_NAME(N) )
!
!          ! Create container
!          CALL Diagn_Create( am_I_Root,                     &
!                             Col       = Collection,        & 
!                             cName     = TRIM( DiagnName ), &
!                             AutoFill  = 0,                 &
!                             ExtNr     = -1,                &
!                             Cat       = -1,                &
!                             Hier      = -1,                &
!                             HcoID     = -1,                &
!                             SpaceDim  =  2,                &
!                             LevIDx    = -1,                &
!                             OutUnit   = 'kg' ,             &
!                             OutOper   = TRIM( OutOper   ), &
!                             WriteFreq = TRIM( WriteFreq ), &
!                             cId       = cId,               &
!                             RC        = RC )
!
!          IF ( RC /= HCO_SUCCESS ) THEN
!             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
!             CALL ERROR_STOP( MSG, LOC ) 
!          ENDIF
!       ENDIF
!    ENDDO
!
!  END SUBROUTINE DiagInit_C_AerSrc
!!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_bl_frac (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGINIT\_BL\_FRAC initializes the distribution 
!  of surface emissions in the boundary layer diagnostics (aka ND12).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_BL_Frac( am_I_Root, Input_Opt, RC )
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
!  26 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_BL_FRAC (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_BL_FRAC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND12 diagnostic is turned off
    IF ( Input_Opt%ND12 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND12_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND12_OUTPUT_FREQ

    ! Check if certain tracer(s) listed for ND12 in input.geos???

    !----------------------------------------------------------------
    ! Create containers for fraction of BL occupied by level L [.]
    !----------------------------------------------------------------

    ! Diagnostic name
    DiagnName = 'BL_FRAC'

    ! Create container
    CALL Diagn_Create( am_I_Root,                     &
                       Col       = Collection,        & 
                       cName     = TRIM( DiagnName ), &
                       AutoFill  = 0,                 &
                       ExtNr     = -1,                &
                       Cat       = -1,                &
                       Hier      = -1,                &
                       HcoID     = -1,                &
                       SpaceDim  =  3,                &
                       LevIDx    = -1,                &
                       OutUnit   = '.' ,              &
                       OutOper   = TRIM( OutOper   ), &
                       WriteFreq = TRIM( WriteFreq ), &
                       cId       = cId,               &
                       RC        = RC )

    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP( MSG, LOC ) 
    ENDIF

  END SUBROUTINE DiagInit_BL_Frac
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_cldconv_flx (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGINIT\_CLDCONV\_FLX initializes the upward
!  mass flux due to wet convection diagnostic (aka ND14).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_CldConv_Flx( am_I_Root, Input_Opt, RC )
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
!  26 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_CLDCONV_FLX (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_CLDCONV_FLX begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND14 diagnostic is turned off
    IF ( Input_Opt%ND14 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND14_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND14_OUTPUT_FREQ

    ! Check if certain tracer(s) listed for ND14 in input.geos???
    ! NEED TO ADD LOOP OVER TRACERS
    
    ! Loop over tracers
    DO N = 1, Input_Opt%N_TRACERS
         
       ! If this tracer number N is scheduled for output in input.geos, 
       ! then define the diagnostic containers for cloud convection mass change
       IF ( ANY( Input_Opt%TINDEX(14,:) == N ) ) THEN

          !----------------------------------------------------------------
          ! Create container for mass change due to cloud convection [kg/s]
          !----------------------------------------------------------------
      
          ! Diagnostic name
          DiagnName = 'CLDCONV_FLX_' // TRIM( Input_Opt%TRACER_NAME(N) )
      
          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cName     = TRIM( DiagnName ),  &
                             AutoFill  = 0,                 &
                             ExtNr     = -1,                &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = -1,                &
                             SpaceDim  =  3,                &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/s' ,           &
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

  END SUBROUTINE DiagInit_CldConv_Flx
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_blmix_flx (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGINIT\_BLMIX\_FLX initializes the upward
!  mass flux from boundary-layer mixing diagnostic (aka ND15).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_BLMix_Flx( am_I_Root, Input_Opt, RC )
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
!  26 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_BLMIX_FLX (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_BLMIX_FLX begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND15 diagnostic is turned off
    IF ( Input_Opt%ND15 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND15_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND15_OUTPUT_FREQ
    
    ! Check if certain tracer(s) listed for ND15 in input.geos???

    ! Loop over tracers
    DO N = 1, Input_Opt%N_TRACERS
         
       ! If this tracer number N is scheduled for output in input.geos, 
       ! then define the diagnostic containers for BL mixing upward flux
       IF ( ANY( Input_Opt%TINDEX(15,:) == N ) ) THEN

          !----------------------------------------------------------------
          ! Create containers for boundary-layer mixing upward mass flux [kg/s]
          !----------------------------------------------------------------
      
          ! Diagnostic name
          DiagnName = 'BLMIX_FLX_' // TRIM( Input_Opt%TRACER_NAME(N) )
      
          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cName     = TRIM( DiagnName ),  &
                             AutoFill  = 0,                 &
                             ExtNr     = -1,                &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = -1,                &
                             SpaceDim  =  3,                &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/s' ,             &
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

  END SUBROUTINE DiagInit_BLMix_Flx
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_precip_frac (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGINIT\_PRECIP\_FRAC initializes the areal
!  fraction of precipitation diagnostic (aka ND16).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Precip_Frac( am_I_Root, Input_Opt, RC )
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
!  26 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, M
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_PRECIP_FRAC (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_PRECIP_FRAC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND16 diagnostic is turned off
    IF ( Input_Opt%ND16 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND16_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND16_OUTPUT_FREQ
    
    DO M = 1, 2          
       SELECT CASE ( M )
          CASE ( 1 )
             !----------------------------------------------------------------
             ! Name container for fraction of grid box with rainout and
             ! washout (large-scale precipitation) [.]
             !----------------------------------------------------------------
             DiagnName = 'PRECIP_FRAC_LS' 
          CASE ( 2 )
             !----------------------------------------------------------------
             ! Name container for fraction of grid box with rainout and
             ! washout (convective precipitation) [.]
             !----------------------------------------------------------------
             DiagnName = 'PRECIP_FRAC_CONV'
          CASE DEFAULT
             IF ( RC /= HCO_SUCCESS ) THEN
                MSG = 'ND16 diagnostic name not defined.'
                CALL ERROR_STOP( MSG, LOC ) 
             ENDIF             
       END SELECT


       ! Create container
       CALL Diagn_Create( am_I_Root,                     &
                          Col       = Collection,        & 
                          cName     = TRIM( DiagnName ), &
                          AutoFill  = 0,                 &
                          ExtNr     = -1,                &
                          Cat       = -1,                &
                          Hier      = -1,                &
                          HcoID     = -1,                &
                          SpaceDim  =  3,                &
                          LevIDx    = -1,                &
                          OutUnit   = '.' ,              &
                          OutOper   = TRIM( OutOper   ), &
                          WriteFreq = TRIM( WriteFreq ), &
                          cId       = cId,               &
                          RC        = RC )
   
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot create ND16 diagnostic: ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF
   ENDDO

  END SUBROUTINE DiagInit_Precip_Frac
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_rain_frac (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGINIT\_RAIN\_FRAC initializes the rainout
!  fraction in precipitation diagnostic (aka ND17).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Rain_Frac( am_I_Root, Input_Opt, RC )
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
!  26 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N, M
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_RAIN_FRAC (diagnostics_mod.F90)' !
    !=======================================================================
    ! DIAGINIT_RAIN_FRAC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND03 diagnostic is turned off
    IF ( Input_Opt%ND17 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND17_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND17_OUTPUT_FREQ
    
    ! Check if certain tracer(s) listed for ND17 in input.geos???

    ! Loop over tracers
    DO N = 1, Input_Opt%N_TRACERS
         
       ! If this tracer number N is scheduled for output in input.geos, 
       ! then define the diagnostic containers for rainout fraction
       IF ( ANY( Input_Opt%TINDEX(17,:) == N ) ) THEN

          DO M = 1, 2          
             SELECT CASE ( M )
                CASE ( 1 )
                !-----------------------------------------------------
                ! Name container for the fraction of soluble tracer
                ! lost to rainout (large-scale precipitation) [.]
                !-----------------------------------------------------
                DiagnName = 'RAIN_FRAC_LS_'                    &
                     // TRIM( Input_Opt%TRACER_NAME(N) )
                CASE ( 2 )
               !------------------------------------------------------
               ! Name container for the fraction of soluble tracer
               ! lost to rainout (convective precipitation) [.]
               !------------------------------------------------------
               ! Diagnostic name
               DiagnName = 'RAIN_FRAC_CONV_'                   &
                    // TRIM( Input_Opt%TRACER_NAME(N) )
                CASE DEFAULT
                   IF ( RC /= HCO_SUCCESS ) THEN
                      MSG = 'ND17 diagnostic name not defined.'
                      CALL ERROR_STOP( MSG, LOC ) 
                   ENDIF             
             END SELECT

      
             ! Create container
             CALL Diagn_Create( am_I_Root,                     &
                                Col       = Collection,        & 
                                cName     = TRIM( DiagnName ), &
                                AutoFill  = 0,                 &
                                ExtNr     = -1,                &
                                Cat       = -1,                &
                                Hier      = -1,                &
                                HcoID     = -1,                &
                                SpaceDim  =  3,                &
                                LevIDx    = -1,                &
                                OutUnit   = '.' ,              &
                                OutOper   = TRIM( OutOper   ), &
                                WriteFreq = TRIM( WriteFreq ), &
                                cId       = cId,               &
                                RC        = RC )
         
             IF ( RC /= HCO_SUCCESS ) THEN
                MSG = 'Cannot create ND17 diagnostic: ' // TRIM(DiagnName)
                CALL ERROR_STOP( MSG, LOC ) 
             ENDIF 
          ENDDO     
       ENDIF
    ENDDO     

  END SUBROUTINE DiagInit_Rain_Frac
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_wash_frac (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGINIT\_WASH\_FRAC initializes the washout
!  fraction diagnostic (aka ND18).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Wash_Frac( am_I_Root, Input_Opt, RC )
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
!  26 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N, M
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_WASH_FRAC (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_WASH_FRAC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND18 diagnostic is turned off
    IF ( Input_Opt%ND18 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND18_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND18_OUTPUT_FREQ
 
    ! Check if certain tracer(s) listed for ND18 in input.geos???

    ! Loop over tracers
    DO N = 1, Input_Opt%N_TRACERS
         
       ! If this tracer number N is scheduled for output in input.geos, 
       ! then define the diagnostic containers for washout fraction
       IF ( ANY( Input_Opt%TINDEX(18,:) == N ) ) THEN

          DO M = 1, 2          
             SELECT CASE ( M )
                CASE ( 1 )
                   !-----------------------------------------------------
                   ! Name container for the fraction of soluble tracer
                   ! lost to washout (large-scale precipitation) [.]
                   !-----------------------------------------------------
                   DiagnName = 'WASH_FRAC_LS'                  &
                        // TRIM( Input_Opt%TRACER_NAME(N) )
                CASE ( 2 )
                   !------------------------------------------------------
                   ! Name container for the fraction of soluble tracer
                   ! lost to washout (convective precipitation) [.]
                   !------------------------------------------------------
                   DiagnName = 'WASH_FRAC_CONV_'               &
                        // TRIM( Input_Opt%TRACER_NAME(N) )
                CASE DEFAULT
                   IF ( RC /= HCO_SUCCESS ) THEN
                      MSG = 'ND18 diagnostic name not defined.'
                      CALL ERROR_STOP( MSG, LOC ) 
                   ENDIF             
             END SELECT

             ! Create container
             CALL Diagn_Create( am_I_Root,                     &
                                Col       = Collection,        & 
                                cName     = TRIM( DiagnName ), &
                                AutoFill  = 0,                 &
                                ExtNr     = -1,                &
                                Cat       = -1,                &
                                Hier      = -1,                &
                                HcoID     = -1,                &
                                SpaceDim  =  3,                &
                                LevIDx    = -1,                &
                                OutUnit   = '.' ,              &
                                OutOper   = TRIM( OutOper   ), &
                                WriteFreq = TRIM( WriteFreq ), &
                                cId       = cId,               &
                                RC        = RC )
         
             IF ( RC /= HCO_SUCCESS ) THEN
                MSG = 'Cannot create ND18 diagnostic: ' // TRIM(DiagnName)
                CALL ERROR_STOP( MSG, LOC ) 
             ENDIF
          ENDDO
       ENDIF
    ENDDO

  END SUBROUTINE DiagInit_Wash_Frac
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_CH4_Loss (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGINIT\_CH4\_LOSS initializes the methane
!  loss diagnostic (aka ND19).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_CH4_Loss( am_I_Root, Input_Opt, RC )
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
!  26 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_CH4_LOSS (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_CH4_LOSS begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND19 diagnostic is turned off
    IF ( Input_Opt%ND19 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND19_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND19_OUTPUT_FREQ
    
    ! Check if certain tracer(s) listed for ND19 in input.geos???

    !----------------------------------------------------------------
    ! Create containers for CH4 removal by OH [kg CH4]
    !----------------------------------------------------------------

    ! Diagnostic name
    DiagnName = 'CH4_LOSS'

    ! Create container
    CALL Diagn_Create( am_I_Root,                     &
                       Col       = Collection,        & 
                       cName     = TRIM( DiagnName ), &
                       AutoFill  = 0,                 &
                       ExtNr     = -1,                &
                       Cat       = -1,                &
                       Hier      = -1,                &
                       HcoID     = -1,                &
                       SpaceDim  =  3,                &
                       LevIDx    = -1,                &
                       OutUnit   = 'kg CH4' ,         &
                       OutOper   = TRIM( OutOper   ), &
                       WriteFreq = TRIM( WriteFreq ), &
                       cId       = cId,               &
                       RC        = RC )

    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP( MSG, LOC ) 
    ENDIF

  END SUBROUTINE DiagInit_CH4_Loss
!EOC

!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: diaginit_Cld_OD (LL in progress)
!!
!! !DESCRIPTION: Subroutine DIAGINIT\_CLD\_OD initializes the cloud
!!  optical depths and cloud fractions diagnostic (aka ND21).
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE DiagInit_Cld_OD( am_I_Root, Input_Opt, RC )
!!
!! !USES:
!!
!    USE Error_Mod,          ONLY : Error_Stop
!    USE GIGC_ErrCode_Mod
!    USE GIGC_Input_Opt_Mod, ONLY : OptInput
!    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
!    USE HCO_Error_Mod
!!
!! !INPUT PARAMETERS:
!!
!    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
!    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
!! 
!! !REVISION HISTORY: 
!!  26 Jan 2015 - E. Lundgren - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    INTEGER              :: cId,      Collection, N
!    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
!    CHARACTER(LEN=60)    :: DiagnName
!    CHARACTER(LEN=255)   :: MSG
!    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_CLD_OD (diagnostics_mod.F90)' 
!
!    !=======================================================================
!    ! DIAGINIT_CLD_OD begins here!
!    !=======================================================================
!      
!    ! Assume successful return
!    RC = GIGC_SUCCESS
!
!    ! Skip if ND21 diagnostic is turned off
!    IF ( Input_Opt%ND21 <= 0 ) RETURN
!
!    ! Get diagnostic parameters from the Input_Opt object
!    Collection = Input_Opt%DIAG_COLLECTION
!    OutOper    = Input_Opt%ND21_OUTPUT_TYPE
!    WriteFreq  = Input_Opt%ND21_OUTPUT_FREQ
!    
!    ! Loop over # of cloud OD and fractions diagnostics
!    DO N = 1, NumDiagND21
!
!       ! If the tracer number is scheduled for output in input.geos, 
!       ! then define the diagnostic container for that tracer.
!       IF ( ANY ( Input_Opt%TINDEX(21,:) == TRCN(N) ) ) THEN
!
!          !----------------------------------------------------------------
!          ! Create containers for ...  [???] - lots of diags. Revisit!
!          !----------------------------------------------------------------
!
!          ! Diagnostic name
!          DiagnName = 'CLD_OD_' // TRIM( Input_Opt%TRACER_NAME(N) )
!
!          ! Create container
!          CALL Diagn_Create( am_I_Root,                     &
!                             Col       = Collection,        & 
!                             cName     = TRIM( DiagnName ), &
!                             AutoFill  = 0,                 &
!                             ExtNr     = -1,                &
!                             Cat       = -1,                &
!                             Hier      = -1,                &
!                             HcoID     = -1,                &
!                             SpaceDim  =  3,                &
!                             LevIDx    = -1,                &
!                             OutUnit   = 'kg' ,             &
!                             OutOper   = TRIM( OutOper   ), &
!                             WriteFreq = TRIM( WriteFreq ), &
!                             cId       = cId,               &
!                             RC        = RC )
!
!          IF ( RC /= HCO_SUCCESS ) THEN
!             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
!             CALL ERROR_STOP( MSG, LOC ) 
!          ENDIF
!       ENDIF
!    ENDDO
!
!  END SUBROUTINE DiagInit_Cld_OD
!!EOC

!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: diaginit_photolysis (LL in progress)
!!
!! !DESCRIPTION: Subroutine DIAGINIT\_PHOTOLYSIS initializes the photolysis
!!  rates (J-values) diagnostic (aka ND22).
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE DiagInit_Photolysis( am_I_Root, Input_Opt, RC )
!!
!! !USES:
!!
!    USE Error_Mod,          ONLY : Error_Stop
!    USE GIGC_ErrCode_Mod
!    USE GIGC_Input_Opt_Mod, ONLY : OptInput
!    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
!    USE HCO_Error_Mod
!!
!! !INPUT PARAMETERS:
!!
!    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
!    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!!
!! !INPUT/OUTPUT PARAMETERS:
!!
!    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
!! 
!! !REVISION HISTORY: 
!!  26 Jan 2015 - E. Lundgren - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    INTEGER              :: cId,      Collection, N
!    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
!    CHARACTER(LEN=60)    :: DiagnName
!    CHARACTER(LEN=255)   :: MSG
!    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_PHOTOLYSIS (diagnostics_mod.F90)' 
!
!    !=======================================================================
!    ! DIAGINIT_PHOTOLYSIS begins here!
!    !=======================================================================
!      
!    ! Assume successful return
!    RC = GIGC_SUCCESS
!
!    ! Skip if ND22 diagnostic is turned off
!    IF ( Input_Opt%ND22 <= 0 ) RETURN
!
!    ! Get diagnostic parameters from the Input_Opt object
!    Collection = Input_Opt%DIAG_COLLECTION
!    OutOper    = Input_Opt%ND22_OUTPUT_TYPE
!    WriteFreq  = Input_Opt%ND22_OUTPUT_FREQ
!    
!    ! Loop over # of photolysis rate diagnostics
!    DO N = 1, NumDiagND22 (NO2, HNO3, H2O2, CH2O, O3, OH source from O3)
!
!       ! If the tracer number is scheduled for output in input.geos, 
!       ! then define the diagnostic container for that tracer.
!       IF ( ANY ( Input_Opt%TINDEX(22,:) == TRCN(N) ) ) THEN
!
!          !----------------------------------------------------------------
!          ! Create containers for J-Values [1/s] - Revisit this!!!
!          !----------------------------------------------------------------
!
!          ! Diagnostic name
!          DiagnName = 'JVALUES_' // TRIM( Input_Opt%TRACER_NAME(N) )
!
!          ! Create container
!          CALL Diagn_Create( am_I_Root,                  &
!                          Col       = Collection,        & 
!                          cName     = TRIM( DiagnName ), &
!                          AutoFill  = 0,                 &
!                          ExtNr     = -1,                &
!                          Cat       = -1,                &
!                          Hier      = -1,                &
!                          HcoID     = -1,                &
!                          SpaceDim  =  3,                &
!                          LevIDx    = -1,                &
!                          OutUnit   = '1/s' ,            &
!                          OutOper   = TRIM( OutOper   ), &
!                          WriteFreq = TRIM( WriteFreq ), &
!                          cId       = cId,               &
!                          RC        = RC )
!
!          IF ( RC /= HCO_SUCCESS ) THEN
!             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
!             CALL ERROR_STOP( MSG, LOC ) 
!          ENDIF
!       ENDIF
!    ENDDO
!
!    ! If simulation is SOA, also output J(GLYX) and J(MGLY)...revisit!!!
!
!  END SUBROUTINE DiagInit_Photolysis
!!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_ew_flx (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGINIT\_EW\_FLX initializes the zonal
!  (east/west) horizontal mass transport flux diagnostic (aka ND24).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_EW_Flx( am_I_Root, Input_Opt, RC )
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
!  26 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_EW_FLX (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_EW_FLX begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND24 diagnostic is turned off
    IF ( Input_Opt%ND24 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND24_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND24_OUTPUT_FREQ
 
    ! Loop over tracers
    DO N = 1, Input_Opt%N_TRACERS
         
       ! If this tracer number N is scheduled for output in input.geos, 
       ! then define the diagnostic containers for E/W flux
       IF ( ANY( Input_Opt%TINDEX(24,:) == N ) ) THEN

          !----------------------------------------------------------------
          ! Create container for east/west mass flux by transport [kg/s]
          !----------------------------------------------------------------
      
          ! Diagnostic name
          DiagnName = 'EW_FLX_' // TRIM( Input_Opt%TRACER_NAME(N) )
      
          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cName     = TRIM( DiagnName ), &
                             AutoFill  = 0,                 &
                             ExtNr     = -1,                &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = -1,                &
                             SpaceDim  =  3,                &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/s' ,           &
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

  END SUBROUTINE DiagInit_EW_Flx
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_ns_flx (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGINIT\_NS\_FLX initializes the meridional 
!  (north/south) horizontal mass transport flux diagnostic (aka ND25).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_NS_Flx( am_I_Root, Input_Opt, RC )
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
!  26 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_NS_FLX (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_NS_FLX begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND25 diagnostic is turned off
    IF ( Input_Opt%ND25 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND25_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND25_OUTPUT_FREQ
 
    ! Loop over tracers
    DO N = 1, Input_Opt%N_TRACERS
         
       ! If this tracer number N is scheduled for output in input.geos, 
       ! then define the diagnostic containers for N/S flux
       IF ( ANY( Input_Opt%TINDEX(25,:) == N ) ) THEN

          !----------------------------------------------------------------
          ! Create container for the north/south mass flux by transport [kg/s]
          !----------------------------------------------------------------
      
          ! Diagnostic name
          DiagnName = 'NS_FLX_' // TRIM( Input_Opt%TRACER_NAME(N) )
      
          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cName     = TRIM( DiagnName ), &
                             AutoFill  = 0,                 &
                             ExtNr     = -1,                &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = -1,                &
                             SpaceDim  =  3,                &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/s' ,           &
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

  END SUBROUTINE DiagInit_NS_Flx
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_vert_flx (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGINIT\_VERT\_FLX initializes the vertical
!  (up/down) mass transport flux diagnostic (aka ND26).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Vert_Flx( am_I_Root, Input_Opt, RC )
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
!  26 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_VERT_FLX (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_VERT_FLX begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND26 diagnostic is turned off
    IF ( Input_Opt%ND26 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND26_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND26_OUTPUT_FREQ
 
    ! Loop over tracers
    DO N = 1, Input_Opt%N_TRACERS
         
       ! If this tracer number N is scheduled for output in input.geos, 
       ! then define the diagnostic containers for N/S flux
       IF ( ANY( Input_Opt%TINDEX(26,:) == N ) ) THEN

          !----------------------------------------------------------------
          ! Create container for the up/down mass flux by transport [kg/s]
          !----------------------------------------------------------------
      
          ! Diagnostic name
          DiagnName = 'VERT_FLX_' // TRIM( Input_Opt%TRACER_NAME(N) )
      
          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cName     = TRIM( DiagnName ), &
                             AutoFill  = 0,                 &
                             ExtNr     = -1,                &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = -1,                &
                             SpaceDim  =  3,                &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/s' ,           &
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

  END SUBROUTINE DiagInit_Vert_Flx
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_strat_flx (LL in progress...obsolete diag???)
!
! !DESCRIPTION: Subroutine DIAGINIT\_STRAT\_FLX initializes the stratospheric
!  influx diagnostics (aka ND27).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Strat_Flx( am_I_Root, Input_Opt, RC )
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
!  26 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_STRAT_FLX (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_STRAT_FLX begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND27 diagnostic is turned off
    IF ( Input_Opt%ND27 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND27_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND27_OUTPUT_FREQ
    
    IF ( Input_Opt%SIM_TYPE == 3 ) THEN     ! NOx-Ox-HC simulation (full chem)

       ! Check if certain tracer(s) listed for ND27 in input.geos???

       !----------------------------------------------------------------
       ! Create container for NOx, Ox, and NHO3 influx from stratosphere [kg/s]
       !----------------------------------------------------------------

       ! Diagnostic name
       DiagnName = 'STRAT_FLX_FULL'

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
                          OutUnit   = 'kg/s' ,           &
                          OutOper   = TRIM( OutOper   ), &
                          WriteFreq = TRIM( WriteFreq ), &
                          cId       = cId,               &
                          RC        = RC )

       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF   
     
    ELSEIF ( Input_Opt%SIM_TYPE == 13 ) THEN      ! H2/HD simulation
  
       !----------------------------------------------------------------
       ! Create container for H2, HD influx from stratosphere [kg/s]
       !----------------------------------------------------------------

       ! Check if certain tracer(s) listed for ND27 in input.geos???

       ! Diagnostic name
       DiagnName = 'STRAT_FLX_H2HD'

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
                          OutUnit   = 'kg/s' ,           &
                          OutOper   = TRIM( OutOper   ), &
                          WriteFreq = TRIM( WriteFreq ), &
                          cId       = cId,               &
                          RC        = RC )

       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF   

    ELSE
  
       !----------------------------------------------------------------
       ! Create container for Ox influx from stratosphere [kg/s]
       !----------------------------------------------------------------

       ! Check if certain tracer(s) listed for ND27 in input.geos???

       ! Diagnostic name
       DiagnName = 'STRAT_FLX_Ox'

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
                          OutUnit   = 'kg/s' ,           &
                          OutOper   = TRIM( OutOper   ), &
                          WriteFreq = TRIM( WriteFreq ), &
                          cId       = cId,               &
                          RC        = RC )

       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF   

     ENDIF

  END SUBROUTINE DiagInit_Strat_Flx
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_landmap (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGINIT\_LANDMAP initializes the land map  
!  diagnostic (aka ND30).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_LandMap( am_I_Root, Input_Opt, RC )
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
!  26 Jan 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_LANDMAP (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_LANDMAP begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND30 diagnostic is turned off
    IF ( Input_Opt%ND30 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND30_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND30_OUTPUT_FREQ
   
    ! Check if certain tracer(s) listed for ND30 in input.geos???

    !----------------------------------------------------------------
    ! Create container for the GMAO land-water indices [.]
    !----------------------------------------------------------------

    ! Diagnostic name
    DiagnName = 'LANDMAP'

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
                       OutUnit   = '.' ,              &
                       OutOper   = TRIM( OutOper   ), &
                       WriteFreq = TRIM( WriteFreq ), &
                       cId       = cId,               &
                       RC        = RC )

    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP( MSG, LOC ) 
    ENDIF   

  END SUBROUTINE DiagInit_LandMap

!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_pressure
!
! !DESCRIPTION: Subroutine DIAGINIT\_PRESSURE initializes the diagnostics
!  for Pressure (aka ND31).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Pressure( am_I_Root, Input_Opt, RC )
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
!  28 Jan 2015 - M. Yannetti - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_PRESSURE (diagnostics_mod.F90)'

    !=======================================================================
    ! DIAGINIT_PRESSURE begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND31 diagnostic is turned off

      IF ( Input_Opt%ND31 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND31_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND31_OUTPUT_FREQ

    ! TODO Check if certain tracer(s) listed for ND31 in input.geos 

    !----------------------------------------------------------------
    ! Create containers for Pressure
    !----------------------------------------------------------------

    ! Diagnostic name
    DiagnName = 'PRESSURE'

    ! Create container
    CALL Diagn_Create( am_I_Root,                     &
                       Col       = Collection,        &
                       cName     = TRIM( DiagnName ), &
                       AutoFill  = 0,                 &
                       ExtNr     = -1,                &
                       Cat       = -1,                &
                       Hier      = -1,                &
                       HcoID     = -1,                &
                       SpaceDim  =  3,                &
                       LevIDx    = -1,                &
                       OutUnit   = '.' ,              &
                       OutOper   = TRIM( OutOper   ), &
                       WriteFreq = TRIM( WriteFreq ), &
                       cId       = cId,               &
                       RC        = RC )

    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP( MSG, LOC )
    ENDIF

  END SUBROUTINE DiagInit_Pressure

!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_trop_height
!
! !DESCRIPTION: Subroutine DIAGINIT\_TROP\_HEIGHT initializes the diagnostics
!  for Tropospheric Height (aka ND55).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Trop_Height( am_I_Root, Input_Opt, RC )
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
!  28 Jan 2015 - M. Yannetti - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_TROP_HEIGHT (diagnostics_mod.F90)'

    !=======================================================================
    ! DIAGINIT_TROP_HEIGHT begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND55 diagnostic is turned off

      IF ( Input_Opt%ND55 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND55_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND55_OUTPUT_FREQ

    ! TODO Check if certain tracer(s) listed for ND55 in input.geos 

    !----------------------------------------------------------------
    ! Create containers for Tropospheric Height
    !----------------------------------------------------------------

    ! Diagnostic name
    DiagnName = 'TROP_HEIGHT'

    ! Create container
    CALL Diagn_Create( am_I_Root,                     &
                       Col       = Collection,        &
                       cName     = TRIM( DiagnName ), &
                       AutoFill  = 0,                 &
                       ExtNr     = -1,                &
                       Cat       = -1,                &
                       Hier      = -1,                &
                       HcoID     = -1,                &
                       SpaceDim  =  3,                &
                       LevIDx    = -1,                &
                       OutUnit   = '.' ,              &
                       OutOper   = TRIM( OutOper   ), &
                       WriteFreq = TRIM( WriteFreq ), &
                       cId       = cId,               &
                       RC        = RC )

    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP( MSG, LOC )
    ENDIF

  END SUBROUTINE DiagInit_Trop_Height

!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_lightning
!
! !DESCRIPTION: Subroutine DIAGINIT\_LIGHTNING initializes the diagnostics
!  for Lightning Flashrate (aka ND56).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Lightning( am_I_Root, Input_Opt, RC )
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
!  28 Jan 2015 - M. Yannetti - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N, M
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_LIGHTNING (diagnostics_mod.F90)'

    !=======================================================================
    ! DIAGINIT_LIGHTNING begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND56 diagnostic is turned off

    IF ( Input_Opt%ND56 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND56_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND56_OUTPUT_FREQ

    ! TODO Check if certain tracer(s) listed for ND56 in input.geos 

    !----------------------------------------------------------------
    ! Create containers for Lightning
    !----------------------------------------------------------------

    ! Loop over ND56 diagnostic tracers
    DO M = 1, 3

    ! Define quantities
      N = Input_Opt%TINDEX(56,M)

      SELECT CASE ( N )
         CASE ( 1 )
           DiagnName = 'LIGHTNING_TOTAL_FLASHRATE'
         CASE ( 2 )
           DiagnName = 'LIGHTNING_INTRACLOUD_FLASHRATE'
         CASE ( 3 )
           DiagnName = 'LIGHTNING_CLOUDGROUND_FLASHRATE'
         CASE DEFAULT
           MSG = 'Lightning index N must not exceed 3!'
           CALL ERROR_STOP ( MSG, LOC )
       END SELECT

       ! Create container
       CALL Diagn_Create( am_I_Root,                     &
                       Col       = Collection,        &
                       cName     = TRIM( DiagnName ), &
                       AutoFill  = 0,                 &
                       ExtNr     = -1,                &
                       Cat       = -1,                &
                       Hier      = -1,                &
                       HcoID     = -1,                &
                       SpaceDim  =  3,                &
                       LevIDx    = -1,                &
                       OutUnit   = '.' ,              &
                       OutOper   = TRIM( OutOper   ), &
                       WriteFreq = TRIM( WriteFreq ), &
                       cId       = cId,               &
                       RC        = RC )

       IF ( RC /= HCO_SUCCESS ) THEN
           MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
           CALL ERROR_STOP( MSG, LOC )
       ENDIF
    ENDDO
  END SUBROUTINE DiagInit_Lightning

!EOC


!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_theta
!
! !DESCRIPTION: Subroutine DIAGINIT\_THETA initializes the diagnostics
!  for Potential Temperature (aka ND57).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Theta( am_I_Root, Input_Opt, RC )
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
!  28 Jan 2015 - M. Yannetti - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N
    CHARACTER(LEN=15)    :: OutOper,  WriteFreq
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_THETA (diagnostics_mod.F90)'

    !=======================================================================
    ! DIAGINIT_THETA begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND57 diagnostic is turned off

      IF ( Input_Opt%ND57 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND57_OUTPUT_TYPE
    WriteFreq  = Input_Opt%ND57_OUTPUT_FREQ

    ! TODO Check if certain tracer(s) listed for ND57 in input.geos 

    !----------------------------------------------------------------
    ! Create containers for Potential Temperature
    !----------------------------------------------------------------

    ! Diagnostic name
    DiagnName = 'THETA'

    ! Create container
    CALL Diagn_Create( am_I_Root,                     &
                       Col       = Collection,        &
                       cName     = TRIM( DiagnName ), &
                       AutoFill  = 0,                 &
                       ExtNr     = -1,                &
                       Cat       = -1,                &
                       Hier      = -1,                &
                       HcoID     = -1,                &
                       SpaceDim  =  3,                &
                       LevIDx    = -1,                &
                       OutUnit   = '.' ,              &
                       OutOper   = TRIM( OutOper   ), &
                       WriteFreq = TRIM( WriteFreq ), &
                       cId       = cId,               &
                       RC        = RC )

    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP( MSG, LOC )
    ENDIF

  END SUBROUTINE DiagInit_Theta

!EOC


END MODULE Diagnostics_Mod
#endif
