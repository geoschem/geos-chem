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
  USE GIGC_ErrCode_Mod
  USE Error_Mod,          ONLY : Error_Stop
  USE GIGC_Input_Opt_Mod, ONLY : OptInput
  USE GIGC_State_Met_Mod, ONLY : MetState
  USE GIGC_State_Chm_Mod, ONLY : ChmState

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Diagnostics_Write

  ! All of the following is currently only active for NETCDF=y
#if defined( NETCDF )
  PUBLIC  :: Diagnostics_Init
  PUBLIC  :: Diagnostics_Final
  PUBLIC  :: DiagnUpdate_NTracers_3D ! added by ewl in early 2015 NewDiag
  PUBLIC  :: CalcDobsonColumn        ! added by ck

!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: DiagInit_Drydep        ! ND44 init
  PRIVATE :: DiagInit_Tracer_Conc   ! ND45 init
  PRIVATE :: DiagInit_Pb_Emiss      ! ND01 init
  PRIVATE :: DiagInit_Rn_Decay      ! ND02 init
  PRIVATE :: DiagInit_BL_Frac       ! ND12 init
  PRIVATE :: DiagInit_CldConv_Flx   ! ND14 init
  PRIVATE :: DiagInit_BlMix_Flx     ! ND15 init
  PRIVATE :: DiagInit_Precip_Frac   ! ND16 init
  PRIVATE :: DiagInit_Rain_Frac     ! ND17 init
  PRIVATE :: DiagInit_Wash_Frac     ! ND18 init
  PRIVATE :: DiagInit_CH4_Loss      ! ND19 init
  PRIVATE :: DiagInit_EW_Flx        ! ND24 init
  PRIVATE :: DiagInit_NS_Flx        ! ND25 init
  PRIVATE :: DiagInit_Vert_Flx      ! ND26 init
  PRIVATE :: DiagInit_Strat_Flx     ! ND27 init
  PRIVATE :: DiagInit_Landmap       ! ND30 init
  PRIVATE :: DiagInit_GridBox       ! General init
  PRIVATE :: DiagInit_Conv_Loss     ! added by ck
  PRIVATE :: DiagInit_Wetdep_Loss   ! added by ck
  PRIVATE :: DiagInit_Tracer_Emis   ! added by ck
  PRIVATE :: DiagInit_Dobson        ! added by ck
!
! !DEFINED PARAMETERS:
!
  ! Prefix of restart file. This file will hold all diagnostics that are 
  ! written out at the end of a simulation (either because their output 
  ! frequency is set to 'End' or the run finishes and these diagnostics
  ! haven't reached the end of their output interval yet).
  CHARACTER(LEN=31), PARAMETER :: RST = 'GEOSCHEM_Restart'
  CHARACTER(LEN=31), PARAMETER :: DGN = 'GEOSCHEM_Diagnostics'

  ! Toggle to enable species diagnostics. This will write out species 
  ! concentrations (in addition to the tracers). Not recommended unless
  ! you have a good reason (ckeller, 8/11/2015).
  LOGICAL, PARAMETER, PUBLIC   :: DiagnSpec = .FALSE.
#endif
!
! !REVISION HISTORY:
!  09 Jan 2015 - C. Keller   - Initial version. 
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
#if defined( NETCDF )
!EOC
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
    USE GRID_MOD,           ONLY : AREA_M2
    USE TIME_MOD,           ONLY : GET_TS_CHEM
#if defined( DEVEL )
    USE TENDENCIES_MOD,     ONLY : TEND_INIT
#endif
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
!  25 Mar 2015 - C. Keller   - Moved UCX initialization to UCX_mod.F
!  06 Nov 2015 - C. Keller   - Added argument OutTimeStamp
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: CollectionID
    INTEGER            :: DeltaYMD, DeltaHMS 
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

    ! Define output write frequency. In ESMF environment, make sure 
    ! diagnostics is always passed to MAPL history!
    CALL DiagnCollection_GetDefaultDelta ( am_I_Root, deltaYMD, deltaHMS, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    CALL DiagnCollection_Create( am_I_Root,                      &
                                 NX           = IIPAR,           &
                                 NY           = JJPAR,           &
                                 NZ           = LLPAR,           &
                                 TS           = TS,              &
                                 AM2          = AM2,             &
                                 PREFIX       = DGN,             &
                                 deltaYMD     = deltaYMD,        &
                                 deltaHMS     = deltaHMS,        &
                                 COL          = CollectionID,    &
                                 OutTimeStamp = HcoDiagnEnd,     &
                                 RC           = RC                )
    IF ( RC /= HCO_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in creating diagnostics collection '//TRIM(DGN), LOC ) 
    ENDIF

    ! Cleanup
    AM2 => NULL()

    ! Save collection ID in Input_Opt%DIAG_COLLECTION for easy future 
    ! reference
    Input_Opt%DIAG_COLLECTION = CollectionID

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

    ! Surface Pressure diagnostic (ND31)
    CALL DIAGINIT_PEDGE( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_PEDGE', LOC )
    ENDIF

    ! Trop. Column Sum of Tracer (ND33)
    CALL DIAGINIT_COLUMNT( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_COLUMNT', LOC )
    ENDIF

    ! UCX diagnostics are now initialized in ucx_mod.F. The UCX diagnostics
    ! currently only include the PSC state, which is written into the HEMCO
    ! restart file for now. Initialize this outside this module to make sure
    ! that this field is also diagnosed if we are not using the DEVEL compiler
    ! switch (ckeller, 3/25/2015). 
!    ! UCX diagnostics
!    IF ( Input_Opt%LUCX ) THEN
!       CALL DIAGINIT_UCX( am_I_Root, Input_Opt, RC )
!       IF ( RC /= GIGC_SUCCESS ) THEN
!          CALL ERROR_STOP( 'Error in DIAGINIT_UCX', LOC ) 
!       ENDIF
!    ENDIF

    ! Convective scavenging loss (ND38)
    CALL DIAGINIT_CONV_LOSS( am_I_Root, Input_Opt, State_Met, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_CONV_LOSS', LOC ) 
    ENDIF

    ! Wetdeposition scavenging loss (ND39)
    CALL DIAGINIT_WETDEP_LOSS( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_WETDEP_LOSS', LOC ) 
    ENDIF

    ! Drydep diagnostic (ND44)
    CALL DIAGINIT_DRYDEP( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_DRYDEP', LOC ) 
    ENDIF

    ! Tracer concentration diagnostics (ND45)
    CALL DIAGINIT_TRACER_CONC( am_I_Root, Input_Opt, State_Chm, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_TRACER_CONC', LOC ) 
    ENDIF
    ! Tropopause diagnostics (ND55)
    CALL DIAGINIT_TR_PAUSE( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_TR_PAUSE', LOC )
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

    ! Tracer Mixing Ratio (IJ-MAX) Diagnostic (ND71)
    CALL DIAGINIT_TRACER_MIXING( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_TRACER_MIXING', LOC )
    ENDIF

    ! Grid box quantities (ND68)
    CALL DIAGINIT_GRIDBOX( am_I_Root, Input_Opt, State_Met, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_GRIDBOX', LOC )
    ENDIF

    ! Tracer emission diagnostics (NEW)
    IF ( .FALSE. ) THEN
       CALL DIAGINIT_TRACER_EMIS( am_I_Root, Input_Opt, State_Met, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DIAGINIT_TRACER_EMIS', LOC ) 
       ENDIF
    ENDIF

    CALL DIAGINIT_DOBSON( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_DOBSON', LOC ) 
    ENDIF

#if defined( DEVEL )
    ! Initialize tendencies
    CALL TEND_INIT( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in TENDENCIES_INIT', LOC ) 
    ENDIF
#endif

    ! Leave with success
    RC = GIGC_SUCCESS

  END SUBROUTINE Diagnostics_Init
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
#if defined( DEVEL )
    USE TENDENCIES_MOD,     ONLY : TEND_CLEANUP
#endif
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

!    ! Finalize diagnostics
!    CALL DiagnCollection_Cleanup( COL = Input_Opt%DIAG_COLLECTION )

#if defined( DEVEL )
    CALL TEND_CLEANUP()
#endif

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
    INTEGER            :: cID,      Collection, D, N
    CHARACTER(LEN=15)  :: OutOper
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
          DiagnName = 'DEPVEL_' // TRIM( Input_Opt%TRACER_NAME(N) )
          cID       = 44000 + N

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cID       = cID,               &
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
                             OkIfExist = .TRUE.,            &
                             RC        = RC )

          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
          
          !----------------------------------------------------------------
          ! Create containers for drydep flux [kg/m2/s]
          !----------------------------------------------------------------

          ! Diagnostic name
          DiagnName = 'DEPFLUX_' // TRIM( Input_Opt%TRACER_NAME(N) )
          cID       = 44500 + N

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        &
                             cID       = cID,               &
                             cName     = TRIM( DiagnName ), &
                             AutoFill  = 0,                 &
                             ExtNr     = -1,                &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = -1,                &
                             SpaceDim  =  2,                &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg m-2 s-1',      &
                             OutOper   = TRIM( OutOper   ), &
                             OkIfExist = .TRUE.,            &
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
  SUBROUTINE DiagInit_Tracer_Conc( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE COMODE_LOOP_MOD,    ONLY : NAMEGAS, NTSPEC, NCSURBAN
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState),   INTENT(IN   ) :: State_Chm  ! Chemistry state 
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
    CHARACTER(LEN=15)  :: OutOper
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
                             RC        = RC )

          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
       ENDIF
    ENDDO

    ! To also write out all species concentrations (not tracers)
    IF ( DiagnSpec ) THEN

       ! Loop over species
       DO N = 1, NTSPEC(NCSURBAN)
   
          !----------------------------------------------------------------
          ! Create containers for species concentrations in molec/cm3
          !----------------------------------------------------------------
   
          ! Diagnostic name
          IF ( TRIM(NAMEGAS(N)) == '' ) CYCLE
          DiagnName = 'SPECIES_CONC_' // TRIM(NAMEGAS(N)) 
   
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
                             OutUnit   = 'molec/cm3',       &
                             OutOper   = TRIM( OutOper   ), &
                             RC        = RC )
   
          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC )
          ENDIF
   
       ENDDO
   
    ENDIF ! DiagnSpec
   
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
    CHARACTER(LEN=15)  :: OutOper
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
    
    ! If the tracer number for lead is scheduled for output in input.geos, 
    ! then define the diagnostic container for 210Pb emissions.
    IF ( ANY ( Input_Opt%TINDEX(1,:) == IDTPB ) ) THEN

       !----------------------------------------------------------------
       ! Create containers for Pb emissions [kg/s]
       !----------------------------------------------------------------

       ! Diagnostic name
       DiagnName = 'EMISS_' // TRIM( Input_Opt%TRACER_NAME( IDTPB ) )

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
                          OutUnit   = 'kg/s',            &
                          OutOper   = TRIM( OutOper   ), &
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
    CHARACTER(LEN=15)    :: OutOper
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
                             OutUnit   = 'kg/s',            &
                             OutOper   = TRIM( OutOper   ), &
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
!    CHARACTER(LEN=15)    :: OutOper
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
!    CHARACTER(LEN=15)    :: OutOper
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
!    CHARACTER(LEN=15)    :: OutOper
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
    CHARACTER(LEN=15)    :: OutOper
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
    CHARACTER(LEN=15)    :: OutOper
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
    CHARACTER(LEN=15)    :: OutOper
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
    CHARACTER(LEN=15)    :: OutOper
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
    CHARACTER(LEN=15)    :: OutOper
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
    CHARACTER(LEN=15)    :: OutOper
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
    CHARACTER(LEN=15)    :: OutOper
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
!    CHARACTER(LEN=15)    :: OutOper
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
!    CHARACTER(LEN=15)    :: OutOper
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
! !IROUTINE: diaginit_ew_flx
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
    CHARACTER(LEN=15)    :: OutOper
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
    CHARACTER(LEN=15)    :: OutOper
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
    CHARACTER(LEN=15)    :: OutOper
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
    CHARACTER(LEN=15)    :: OutOper
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
    CHARACTER(LEN=15)    :: OutOper 
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
! !IROUTINE: diaginit_pedge
!
! !DESCRIPTION: Subroutine DIAGINIT\_PEDGE initializes the diagnostics
!  for Pressure (aka ND31).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Pedge( am_I_Root, Input_Opt, RC )
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
    CHARACTER(LEN=15)    :: OutOper 
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_PEDGE (diagnostics_mod.F90)'

    !=======================================================================
    ! DIAGINIT_PEDGE begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND31 diagnostic is turned off

      IF ( Input_Opt%ND31 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND31_OUTPUT_TYPE

    ! TODO Check if certain tracer(s) listed for ND31 in input.geos 

    !----------------------------------------------------------------
    ! Create containers for Pressure
    !----------------------------------------------------------------

    ! Diagnostic name
    DiagnName = 'PEDGE'

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
                       OutUnit   = 'mb',              &
                       OutOper   = TRIM( OutOper   ), &
                       cId       = cId,               &
                       RC        = RC )

    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP( MSG, LOC )
    ENDIF

  END SUBROUTINE DiagInit_Pedge

!EOC


!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_columnt
!
! !DESCRIPTION: Subroutine DIAGINIT\_COLUMNT initializes the diagnostics
!  for Trop. Column Sum (aka ND33).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_ColumnT( am_I_Root, Input_Opt, RC )
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
!  29 Jan 2015 - M. Yannetti - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N, M
    CHARACTER(LEN=15)    :: OutOper 
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_COLUMNT (diagnostics_mod.F90)'

    !=======================================================================
    ! DIAGINIT_COLUMNT begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND33 diagnostic is turned off

    IF ( Input_Opt%ND33 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND33_OUTPUT_TYPE

    ! TODO Check if certain tracer(s) listed for ND33 in input.geos 

    !----------------------------------------------------------------
    ! Create containers for Column Tracer
    !----------------------------------------------------------------

    ! Diagnostic name
    DiagnName = 'COLUMN-T'

    DO M = 1, Input_Opt%TMAX(33)
      N  = Input_Opt%TINDEX(33,M)
      IF ( N > Input_Opt%N_TRACERS ) CYCLE


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
                       OutUnit   = 'kg' ,             &
                       OutOper   = TRIM( OutOper   ), &
                       cId       = cId,               &
                       RC        = RC )

      IF ( RC /= HCO_SUCCESS ) THEN
         MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
         CALL ERROR_STOP( MSG, LOC )
      ENDIF
    ENDDO
  END SUBROUTINE DiagInit_ColumnT

!EOC


!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_tr_pause
!
! !DESCRIPTION: Subroutine DIAGINIT\_TR\_PAUSE initializes the diagnostics
!  for Tropopause height/level/pressure. (aka ND55).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Tr_Pause( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE Error_Mod,          ONLY : Error_Stop
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE HCO_Error_Mod
    USE CMN_DIAG_MOD,       ONLY : PD55
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
    CHARACTER(LEN=15)    :: OutOper 
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=40)    :: THEUNIT
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_TR_PAUSE (diagnostics_mod.F90)'

    !=======================================================================
    ! DIAGINIT_TR_PAUSE begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND55 diagnostic is turned off

      IF ( Input_Opt%ND55 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND55_OUTPUT_TYPE

    ! TODO Check if certain tracer(s) listed for ND55 in input.geos 

    !----------------------------------------------------------------
    ! Create containers for Tropopause Diags
    !----------------------------------------------------------------

    ! Diagnostic name
    DiagnName = 'TR-PAUSE'

    ! Loop over ND55 diagnostic tracers
    DO M = 1, Input_Opt%TMAX(55)

    ! Define quantities
      N = Input_Opt%TINDEX(55,M)
      IF ( N > PD55 ) CYCLE

      ! Pick the appropriate unit string
      SELECT CASE ( N )
        CASE ( 1 )
          THEUNIT = 'unitless'
        CASE ( 2 )
          THEUNIT = 'km'
        CASE ( 3 )
          THEUNIT = 'mb'
      END SELECT

    ! pulled from diag3... not sure where these vars are... (mdy)
!     ARRAY(:,:,1) = AD55(:,:,N) / SCALEDYN

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
                       OutUnit   = TRIM( THEUNIT   ), &
                       OutOper   = TRIM( OutOper   ), &
                       cId       = cId,               &
                       RC        = RC )

      IF ( RC /= HCO_SUCCESS ) THEN
         MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
         CALL ERROR_STOP( MSG, LOC )
      ENDIF

  ENDDO

  END SUBROUTINE DiagInit_Tr_Pause

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
    CHARACTER(LEN=15)    :: OutOper 
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
                       OutUnit   = 'flashes/min/km2', &
                       OutOper   = TRIM( OutOper   ), &
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
    CHARACTER(LEN=15)    :: OutOper 
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

    ! TODO Check if certain tracer(s) listed for ND57 in input.geos 

    !----------------------------------------------------------------
    ! Create containers for Potential Temperature
    !----------------------------------------------------------------

    ! Diagnostic name
    DiagnName = 'THETA'

    ! another import from diag3... (mdy)
!    ARRAY(:,:,1:LD57) = AD57(:,:,1:LD57) / SCALEDIAG

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
                       OutUnit   = 'K' ,              &
                       OutOper   = TRIM( OutOper   ), &
                       cId       = cId,               &
                       RC        = RC )

    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP( MSG, LOC )
    ENDIF

  END SUBROUTINE DiagInit_Theta

!EOC


!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_tracer_mixing
!
! !DESCRIPTION: Subroutine DIAGINIT\_TRACER\_MIXING initializes the diagnostics
!  for Tracer Mixing Ratio at Surface (aka ND71).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Tracer_Mixing( am_I_Root, Input_Opt, RC )
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
!  30 Jan 2015 - M. Yannetti - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N, M
    CHARACTER(LEN=15)    :: OutOper 
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGINIT_TRACER_MIXING (diagnostics_mod.F90)'

    !=======================================================================
    ! DIAGINIT_TRACER_MIXING begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND71 diagnostic is turned off

      IF ( Input_Opt%ND71 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND71_OUTPUT_TYPE

    ! TODO Check if certain tracer(s) listed for ND71 in input.geos 

    !----------------------------------------------------------------
    ! Create containers for Tracer Mixing Ratio
    !----------------------------------------------------------------

    ! Diagnostic name
    DiagnName = 'IJ-MAX'
    DO M = 1, Input_Opt%TMAX(71)
       N = Input_Opt%TINDEX(71,M)
       IF ( N > Input_Opt%N_TRACERS ) CYCLE

!       ! line from diag3
!       ARRAY(:,:,1) = AD71(:,:,N)/AD71_COUNT

       ! Create container
       CALL Diagn_Create( am_I_Root,                  &
                       Col       = Collection,        &
                       cName     = TRIM( DiagnName ), &
                       AutoFill  = 0,                 &
                       ExtNr     = -1,                &
                       Cat       = -1,                &
                       Hier      = -1,                &
                       HcoID     = -1,                &
                       SpaceDim  =  3,                &
                       LevIDx    = -1,                &
                       OutUnit   = 'v/v' ,            &
                       OutOper   = TRIM( OutOper   ), &
                       cId       = cId,               &
                       RC        = RC )

       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC )
       ENDIF
    ENDDO

  END SUBROUTINE DiagInit_Tracer_Mixing

!EOC
  
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagnupdate_ntracers_3D
!
! !DESCRIPTION: Subroutine DIAGNUPDATE\_NTRACERS\_3D updates a generic set of
!  3D diagnostics that are distinguished from each other by tracer.
!  
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnUpdate_NTracers_3D( am_I_Root, DiagnPrefix, DiagnNum, &
                                     DiagnArray, Input_Opt, RC )
!
! !USES:
!
    USE Error_Mod,          ONLY : Error_Stop
    USE GIGC_ErrCode_Mod    ! Is this needed?
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Update
    USE HCO_Error_Mod
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR ! Works for all cases???
!
! !INPUT PARAMETERS:
!
    LOGICAL,           INTENT(IN) :: am_I_Root          ! Is this the root CPU?
    CHARACTER(LEN=60), INTENT(IN) :: DiagnPrefix         ! Diag name prefix
    INTEGER,           INTENT(IN) :: DiagnNum            ! Diagn # (eg. 24)
    REAL(fp),          INTENT(IN) :: DiagnArray(:,:,:,:) ! data
    TYPE(OptInput),    INTENT(IN) :: Input_Opt          ! Input options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,           INTENT(INOUT) :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  6 Feb 2015 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: cId,      Collection, N
    CHARACTER(LEN=15)    :: OutOper 
    CHARACTER(LEN=255)   :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGNUPDATE_NTRACERS_3D (diagnostics_mod.F90)' 
    INTEGER              :: HCRC
    REAL(fp), TARGET     :: DiagnArray_tracer(IIPAR,JJPAR,LLPAR)
    REAL(fp), POINTER    :: Ptr3D(:,:,:)

    !=======================================================================
    ! DIAGNUPDATE_NTRACERS_3D begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Assume this is called if Input_Opt%NDXX > 0

    ! Loop over tracers
    DO N = 1, Input_Opt%N_TRACERS

       ! If this tracer number N is scheduled for output in input.geos, 
       ! then update its diagnostic container for transport flux
       IF ( ANY( Input_Opt%TINDEX(DiagnNum,:) == N ) ) THEN
         
          !----------------------------------------------------------------
          ! Update diagnostic container
          !----------------------------------------------------------------
      
          ! Diagnostic name
          DiagnName = TRIM( DiagnPrefix )                              & 
                      // TRIM( Input_Opt%TRACER_NAME( N ) )

          ! Assign temporary 3D array
          DiagnArray_tracer = DiagnArray(:,:,:,N)

          ! Point to the array
          Ptr3D => DiagnArray_tracer

          ! Create container
          CALL Diagn_Update( am_I_Root,                                &
                             cName     = TRIM( DiagnName ),            &
                             Array3D   = Ptr3D,                        &
                             RC        = RC )

          ! Free the point before error handling
          Ptr3D => NULL()

          ! Stop with error if the diagnostics update was unsuccessful.
          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot update diagnostic: ' // TRIM( DiagnName )
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF  
       ENDIF
    ENDDO

  END SUBROUTINE DiagnUpdate_NTracers_3D
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_gridbox
!
! !DESCRIPTION: Subroutine DIAGINIT\_GRIDBOX initializes the grid box quantity 
!  diagnostic (aka ND68).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DIAGINIT_GRIDBOX( am_I_Root, Input_Opt, State_Met, RC )
!
! !USES:
!
    USE CMN_GCTM_MOD,       ONLY : XNUMOLAIR
    USE Error_Mod,          ONLY : Error_Stop
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE HCO_Error_Mod
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN   ) :: State_Met  ! Met state
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
    INTEGER            :: cId, Collection, N
    INTEGER            :: SpaceDim 
    REAL(hp)           :: ScaleFact
    CHARACTER(LEN=15)  :: OutOper, OutUnit
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGINIT_GRIDBOX (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_GRIDBOX begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND68 diagnostic is turned off
    IF ( Input_Opt%ND68 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND68_OUTPUT_TYPE
      
    ! There are six diagnostics 
    DO N = 1,6

       !----------------------------------------------------------------
       ! Create containers
       !----------------------------------------------------------------

       SELECT CASE ( N )
          CASE ( 1 )
             DiagnName = 'BOXHEIGHT'
             OutUnit   = 'm'
             ScaleFact = 1.0_hp
             SpaceDim  = 3
          CASE ( 2 )
             DiagnName = 'AIRMASS'
             OutUnit   = 'kg'
             ScaleFact = 1.0_hp
             SpaceDim  = 3
          CASE ( 3 )
             DiagnName = 'AIRDENSITY'
             OutUnit   = 'kg m-3'
             ScaleFact = 1.0_hp
             SpaceDim  = 3
          CASE ( 4 )
             IF ( .NOT. ASSOCIATED(State_Met%AVGW) ) CYCLE
             DiagnName = 'AVGW'
             OutUnit   = 'v/v'
             ScaleFact = 1.0_hp
             SpaceDim  = 3
          CASE ( 5 )
             DiagnName = 'TROPP_PRESSURE'
             OutUnit   = 'hPa'
             ScaleFact = 1.0_hp
             SpaceDim  = 2
          CASE ( 6 )
             DiagnName = 'TROPP_LEVEL'
             OutUnit   = 'count'
             ScaleFact = 1.0_hp
             SpaceDim  = 2
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
                          SpaceDim  =  SpaceDim,         &
                          LevIDx    = -1,                &
                          OutUnit   = TRIM( OutUnit   ), &
                          OutOper   = TRIM( OutOper   ), &
                          ScaleFact = ScaleFact,         &
                          RC        = RC )

       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF
    ENDDO

  END SUBROUTINE DiagInit_GridBox
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DiagInit_Tracer_Emis
!
! !DESCRIPTION: Subroutine DiagInit\_Tracer\_Emis initializes diagnostics for 
!  total species emissions diagnostics. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Tracer_Emis( am_I_Root, Input_Opt, State_Met, RC ) 
!
! !USES:
!
    USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoID
    USE TRACERID_MOD
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt   ! Input Options object
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  05 Mar 2015 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: ID
    INTEGER            :: cId,      Collection, N
    CHARACTER(LEN=15)  :: OutOper
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGINIT_TRACER_EMIS (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_TRACER_CONC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    ! Use same output frequency and operations as for tracer concentrations.
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND45_OUTPUT_TYPE
 
    ! Loop over # of species 
    DO N = 1, Input_Opt%N_TRACERS
       
       ! HEMCO ID
       ID = GetHcoID( TrcID = N )
 
       ! Restrict diagnostics to these species
       IF ( N /= IDTNO    .AND. N /= IDTCO     .AND. &
            N /= IDTALK4  .AND. N /= IDTISOP   .AND. &
            N /= IDTHNO3  .AND. N /= IDTACET   .AND. &
            N /= IDTMEK   .AND. N /= IDTALD2   .AND. &
            N /= IDTPRPE  .AND. N /= IDTC3H8   .AND. &
            N /= IDTC2H6  .AND. N /= IDTDMS    .AND. &
            N /= IDTSO2   .AND. N /= IDTSO4    .AND. &
            N /= IDTNH3   .AND. N /= IDTBCPI   .AND. &
            N /= IDTOCPI  .AND. N /= IDTBCPO   .AND. &
            N /= IDTOCPO  .AND. N /= IDTDST1   .AND. &
            N /= IDTDST2  .AND. N /= IDTDST3   .AND. &
            N /= IDTDST4  .AND. N /= IDTSALA   .AND. &
            N /= IDTSALC  .AND. N /= IDTBr2    .AND. &
            N /= IDTBrO   .AND. N /= IDTCH2Br2 .AND. &
            N /= IDTCH3Br .AND. N /= IDTO3             ) THEN
          ID = -1
       ENDIF
 
       ! If this is an emission tracer, add diagnostics for emissions. 
       IF ( ID > 0 ) THEN 

          !----------------------------------------------------------------
          ! Create container for emission flux (kg/s) 
          !----------------------------------------------------------------

          ! Diagnostic name
          DiagnName = 'TRACER_EMIS_' // TRIM( Input_Opt%TRACER_NAME(N) )

          ! Diagnostics ID
          cID = 10000 + ID

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cID       = cID,               &
                             cName     = TRIM( DiagnName ), &
                             AutoFill  = 0,                 &
                             ExtNr     = -1,                &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = ID,                &
                             SpaceDim  =  3,                &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/s',            &
                             OutOper   = TRIM( OutOper   ), &
                             RC        = RC                  )

          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
       ENDIF
    ENDDO

  END SUBROUTINE DiagInit_Tracer_Emis
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_conv_loss
!
! !DESCRIPTION: Subroutine DIAGINIT\_CONV\_LOSS initializes the convective 
!  scavenging loss diagnostic (aka ND38). For now, this is a 2D column 
!  diagnostics, i.e. we only write out the loss due to convection of the
!  entire column.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DIAGINIT_CONV_LOSS( am_I_Root, Input_Opt, State_Met, RC )
!
! !USES:
!
    USE WETSCAV_MOD,    ONLY : GET_WETDEP_NSOL, GET_WETDEP_IDWETD
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN   ) :: State_Met  ! Met state
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Mar 2015 - C. Keller   - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: N, NN, M
    INTEGER            :: cId, Collection
    CHARACTER(LEN=15)  :: OutOper, OutUnit
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGINIT_CONV_LOSS (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_CONV_LOSS begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND39 diagnostic is turned off
    IF ( Input_Opt%ND38 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    ! Use same output frequency and operations as for tracer concentrations.
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND38_OUTPUT_TYPE

    ! Get number of soluble species
    M = GET_WETDEP_NSOL()

    ! Loop over # of species 
    DO N = 1, M

       ! Get GEOS-Chem tracer number
       NN = GET_WETDEP_IDWETD( N )

       ! Check if this is a species asked in input.geos
       IF ( ANY( Input_Opt%TINDEX(38,:) == NN ) ) THEN

          !----------------------------------------------------------------
          ! Create container for convective loss (kg/s) 
          !----------------------------------------------------------------

          ! Diagnostic name
          DiagnName = 'CONV_LOSS_' // TRIM( Input_Opt%TRACER_NAME(NN) )

          ! Define diagnostics ID
          cID = 38000 + NN

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cID       = cID,               &
                             cName     = TRIM( DiagnName ), &
                             AutoFill  = 0,                 &
                             ExtNr     = -1,                &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = -1,                &
                             SpaceDim  =  2,                & ! 2D for now!!
                             LevIDx    = -1,                & ! sum over all vert. levels
                             OutUnit   = 'kg/m2/s',         &
                             OutOper   = TRIM( OutOper   ), &
                             RC        = RC                  )

          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
       ENDIF
    ENDDO !N

  END SUBROUTINE DIAGINIT_CONV_LOSS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_wetdep_loss
!
! !DESCRIPTION: Subroutine DIAGINIT\_WETDEP\_LOSS initializes the wet 
!  deposition scavenging loss diagnostic (aka ND39).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DIAGINIT_WETDEP_LOSS( am_I_Root, Input_Opt,     &
                                   State_Met, State_Chm, RC )
!
! !USES:
!
    USE Species_Mod, ONLY : Species   
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Mar 2015 - C. Keller   - Initial version
!   3 Sep 2015 - R. Yantosca - Now get wetdep species from State_Chm
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: N, NN, M
    INTEGER            :: cId, Collection
    CHARACTER(LEN=15)  :: OutOper, OutUnit
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGINIT_WETDEP_LOSS (diagnostics_mod.F90)' 

    ! Pointers
    TYPE(Species), POINTER :: ThisSpc

    !=======================================================================
    ! DIAGINIT_WETDEP_LOSS begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND39 diagnostic is turned off
    IF ( Input_Opt%ND39 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    ! Use same output frequency and operations as for tracer concentrations.
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND39_OUTPUT_TYPE

    ! Get number of soluble species
    M = State_Chm%nWetDep

    ! Loop over all species
    DO N = 1, State_Chm%nSpecies

       ! Get info about the Nth species from the species database
       ThisSpc => State_Chm%SpcData(N)%Info

       ! Skip if this is not a wet-depositing species
       IF ( .not. ThisSpc%Is_WetDep ) CYCLE

       ! Wetdep species index
       NN = ThisSpc%WetDepId

       ! Check if this is a species asked in input.geos
       IF ( ANY( Input_Opt%TINDEX(39,:) == NN ) ) THEN

          !----------------------------------------------------------------
          ! Create container for wetdep loss (kg/s) 
          !----------------------------------------------------------------

          ! Diagnostic name
          DiagnName = 'WETDEP_LOSS_' // TRIM( Input_Opt%TRACER_NAME(NN) )

          ! Define diagnostics ID
          cID = 39000 + NN

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cID       = cID,               &
                             cName     = TRIM( DiagnName ), &
                             AutoFill  = 0,                 &
                             ExtNr     = -1,                &
                             Cat       = -1,                &
                             Hier      = -1,                &
                             HcoID     = -1,                &
                             SpaceDim  =  2,                &
                             LevIDx    = -1,                &
                             OutUnit   = 'kg/m2/s',         &
                             OutOper   = TRIM( OutOper   ), &
                             RC        = RC                  )

          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
       ENDIF

       ! Free pointer
       ThisSpc => NULL()

    ENDDO !N

  END SUBROUTINE DIAGINIT_WETDEP_LOSS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diaginit_dobson
!
! !DESCRIPTION: Subroutine DIAGINIT\_DOBSON initializes the O3 dobson column
! diagnostics. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagInit_Dobson( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE TRACERID_MOD, ONLY : IDTO3
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
!  07 Jul 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: Collection, N
    CHARACTER(LEN=15)  :: OutOper
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGINIT_DOBSON (diagnostics_mod.F)' 
    
    !=======================================================================
    ! DIAGINIT_DOBSON begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Nothing to do if O3 is not a tracer
    IF ( IDTO3 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = 'Mean'
     
    ! troposphere and total column
    DO N = 1, 2
 
       ! Diagnostic name
       IF ( N == 1 ) THEN
          DiagnName = 'O3_COLUMN' 
       ELSEIF ( N == 2 ) THEN
          DiagnName = 'O3_TROPCOLUMN' 
       ENDIF

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
                          OutUnit   = 'dobson',          &
                          OutOper   = TRIM( OutOper   ), &
                          OkIfExist = .TRUE.,            &
                          RC        = RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL ERROR_STOP( 'Cannot create diagnostics '//TRIM(DiagnName), LOC ) 
      ENDIF  
 
   ENDDO
   
  END SUBROUTINE DiagInit_Dobson
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CalcDobsonColumn 
!
! !DESCRIPTION: Subroutine CalcDobsonColumn calculates total ozone column in
! dobsons and adds them to the GEOS-Chem diagnostics. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CalcDobsonColumn( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
    USE CMN_GCTM_MOD,       ONLY : AIRMW, AVO,   g0
    USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
    USE TRACERID_MOD,       ONLY : IDTO3
    USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_TROP
    USE PRESSURE_MOD,       ONLY : GET_PEDGE
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root  ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
    TYPE(MetState), INTENT(IN   ) :: State_Met  ! Met state
    TYPE(ChmState), INTENT(IN   ) :: State_Chm  ! Chemistry state 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  07 Jul 2015 - C. Keller   - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DiagnCont), POINTER :: DgnPtr => NULL()

    INTEGER                  :: I, J, L

    REAL(fp)                 :: constant
    REAL(fp)                 :: DP, O3vv, DU
    REAL(fp)                 :: TROPO3(IIPAR,JJPAR)
    REAL(fp)                 :: TOTO3 (IIPAR,JJPAR)

    LOGICAL, SAVE            :: FIRST      = .TRUE.
    LOGICAL, SAVE            :: TropDiagn  = .FALSE.
    LOGICAL, SAVE            :: TotDiagn   = .FALSE.

    CHARACTER(LEN=255)       :: MSG
    CHARACTER(LEN=255)       :: LOC = 'CalcDobsonColumn (diagnostics_mod.F)' 
    
    !=======================================================================
    ! CalcDobsonColumn begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Nothing to do if O3 is not a tracer
    IF ( IDTO3 <= 0 ) RETURN

    ! On first call, check if any of the two diagnostics is defined
    IF ( FIRST ) THEN
       CALL DiagnCont_Find( -1, -1, -1, -1, -1, 'O3_COLUMN', &
                            -1, TotDiagn, DgnPtr )
       DgnPtr => NULL()

       CALL DiagnCont_Find( -1, -1, -1, -1, -1, 'O3_TROPCOLUMN', &
                            -1, TropDiagn, DgnPtr ) 

       DgnPtr => NULL()
       FIRST = .FALSE.
    ENDIF

    ! Nothing to do if none of the diagnostics exist
    IF ( .NOT. TotDiagn .AND. .NOT. TropDiagn ) RETURN 

    ! Initialize values
    TROPO3 = 0.0_fp
    TOTO3  = 0.0_fp

    ! Constant
    constant = 0.01_fp * AVO / ( g0 * ( AIRMW/1000.0_fp) )

    ! Do for all levels
    !$OMP PARALLEL DO &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, DP, O3vv, DU )
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Pressure difference in hPa
       DP = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

       ! Ozone in v/v
       O3vv = State_Chm%Tracers(I,J,L,IDTO3) * Input_Opt%TCVV(IDTO3) / State_Met%AD(I,J,L)

       ! Calculate O3 in DU for this grid box 
       DU = O3vv * DP * constant / 2.69e16_fp

       ! Add to totals
       IF ( TotDiagn ) THEN
          TOTO3(I,J) = TOTO3(I,J) + DU
       ENDIF
       IF ( TropDiagn ) THEN
          IF ( ITS_IN_THE_TROP(I,J,L,State_Met) ) THEN
             TROPO3(I,J) = TROPO3(I,J) + DU
          ENDIF
       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Update diagnostics
    IF ( TotDiagn ) THEN
       CALL Diagn_Update( am_I_Root, cName = 'O3_COLUMN',     Array2D = TOTO3,  RC=RC )
    ENDIF
    IF ( TropDiagn ) THEN
       CALL Diagn_Update( am_I_Root, cName = 'O3_TROPCOLUMN', Array2D = TROPO3, RC=RC )
    ENDIF

    ! Return w/ success
    RC = GIGC_SUCCESS
   
  END SUBROUTINE CalcDobsonColumn
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: TotalsToLogfile 
!
! !DESCRIPTION: Subroutine TotalsToLogfile is a helper routine to print the 
! monthly emission totals to the GEOS-Chem logfile.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE TotalsToLogfile( am_I_Root, Input_Opt, RC ) 
!
! !USES:
!
    USE HCO_ERROR_MOD
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE TIME_MOD,           ONLY : GET_YEAR, GET_MONTH 
    USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoID
    USE HCO_DIAGN_MOD,      ONLY : Diagn_TotalGet
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input opts
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  15 Mar 2015 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER       :: I, cID, HCRC
    INTEGER       :: YEAR, MONTH 
    INTEGER, SAVE :: SAVEMONTH = -999
    REAL(sp)      :: TOTAL
    LOGICAL       :: FOUND

    !=================================================================
    ! TotalsToLogfile begins here!
    !=================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Don't do anything if not root
    IF ( .NOT. Am_I_Root ) RETURN

    ! Get this month
    MONTH = GET_MONTH()

    ! Print totals if it's a new month, but not on first call
    IF ( (SAVEMONTH /= MONTH) .AND. (SAVEMONTH > 0) ) THEN 

       ! Get year/month of previous month
       YEAR = GET_YEAR()
       IF ( MONTH == 1 ) THEN
          YEAR  = YEAR -1
          MONTH = 12
       ELSE
          MONTH = MONTH - 1
       ENDIF

       ! Print header
       WRITE(6,*    ) ''
       WRITE(6,'(a)') REPEAT( '-', 79 )
       WRITE(6,100  ) MONTH, YEAR

       ! Loop over all tracers
       DO I = 1, Input_Opt%N_TRACERS

          ! Only if it's a HEMCO species...
          cID = GetHcoID( TrcID=I )
          IF ( cID <= 0 ) CYCLE
          
          ! Define diagnostics ID
          cID = 10000 + cID

          ! Get the total [kg] 
          CALL Diagn_TotalGet( am_I_Root,                         & 
                             cID     = cID,                       &
                             Found   = Found,                     &
                             Total   = Total,                     &
                             COL     = Input_Opt%DIAG_COLLECTION, &
                             Reset   = .TRUE.,                    &
                             RC      = HCRC                        )
          IF ( .NOT. FOUND ) CYCLE

          ! Only if there have been any emissions...
          IF ( Total == 0.0_fp ) CYCLE
   
          ! Convert total from kg to Mg
          Total = Total / 1000.0_fp

          ! Write out
          WRITE(6,101) TRIM(Input_Opt%TRACER_NAME(I)), Total
       ENDDO
    ENDIF
    
    ! Update month counter
    SAVEMONTH = MONTH

100 FORMAT( 'Emissions in month ', i2, ' of year ', i4, ':')
101 FORMAT( a9, ': ', e11.3, ' Mg/month')

  END SUBROUTINE TotalsToLogfile 
!EOC
#endif
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
  SUBROUTINE Diagnostics_Write ( am_I_Root, Input_Opt, RESTART, RC ) 
!
! !USES:
!
    USE HCO_STATE_MOD,      ONLY : HCO_STATE
    USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoState
    USE HCOI_GC_MAIN_MOD,   ONLY : HCOI_GC_WriteDiagn 
    USE HCOIO_Diagn_Mod,    ONLY : HCOIO_Diagn_WriteOut
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
    LOGICAL,          INTENT(IN   )  :: RESTART    ! Write restart file? 
    TYPE(OptInput),   INTENT(IN )    :: Input_Opt  ! Input Options object
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

    ! Write HEMCO diagnostics
    CALL HCOI_GC_WriteDiagn( am_I_Root, Input_Opt, RESTART, RC )
    IF ( RC /= GIGC_SUCCESS ) RETURN

    ! Write netCDF GEOS-Chem diagnostics
#if defined( NETCDF )
    ! Get pointer to HEMCO state object.
    CALL GetHcoState( HcoState )
    IF ( .NOT. ASSOCIATED(HcoState) ) THEN
       CALL ERROR_STOP( 'Cannot get HEMCO state object', LOC )
    ENDIF

    !-----------------------------------------------------------------------
    ! Eventually write out emission totals to GEOS-Chem logfile
    !-----------------------------------------------------------------------
    CALL TotalsToLogfile( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN 
       CALL ERROR_STOP ('Error in TotalsToLogfile', LOC ) 
    ENDIF

    !-----------------------------------------------------------------------
    ! RESTART: write out all diagnostics. Use current time stamp and save into
    ! restart file.
    ! GEOS-Chem restart file is currently not defined ...
    !-----------------------------------------------------------------------
    IF ( RESTART ) THEN
!       CALL HCOIO_DIAGN_WRITEOUT( am_I_Root,                                &
!                                  HcoState,                                 &
!                                  ForceWrite  = .TRUE.,                     &
!                                  UsePrevTime = .FALSE.,                    & 
!                                  PREFIX      = RST,                        &
!                                  COL         = Input_Opt%DIAG_COLLECTION,  &
!                                  OnlyIfFirst = .TRUE.,                     &
!                                  RC          = RC                         )
!   
!       IF ( RC /= HCO_SUCCESS ) THEN
!          CALL ERROR_STOP( 'Diagnostics restart write error', LOC ) 
!       ENDIF

    !-----------------------------------------------------------------------
    ! Not restart: write out regular diagnostics. Use current time stamp.
    !-----------------------------------------------------------------------
    ELSE
       CALL HCOIO_DIAGN_WRITEOUT( am_I_Root,                                & 
                                  HcoState,                                 &
                                  ForceWrite  = .FALSE.,                    &
                                  UsePrevTime = .FALSE.,                    &
                                  COL         = Input_Opt%DIAG_COLLECTION,  &
                                  RC          = RC                         )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'Diagnostics write error', LOC ) 
       ENDIF

    ENDIF

    ! Free pointer
    HcoState => NULL()
#endif

    ! Leave w/ success
    RC = GIGC_SUCCESS

  END SUBROUTINE Diagnostics_Write
!EOC
END MODULE Diagnostics_Mod
