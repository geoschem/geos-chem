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

  ! All of the following is currently only active in development mode:
#if defined( DEVEL )
  PUBLIC  :: Diagnostics_Init
  PUBLIC  :: Diagnostics_Final
  PUBLIC  :: CalcDobsonColumn
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: DiagInit_Drydep
  PRIVATE :: DiagInit_Conv_Loss
  PRIVATE :: DiagInit_Wetdep_Loss
  PRIVATE :: DiagInit_Tracer_Conc
  PRIVATE :: DiagInit_Tracer_Emis
  PRIVATE :: DiagInit_GridBox
  PRIVATE :: DiagInit_Dobson
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
#if defined( DEVEL )
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
    USE TENDENCIES_MOD,     ONLY : TEND_INIT
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

    ! Save number of wet-depositing species for later use
    nWetDep = State_Chm%nWetDep

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

    ! Initialize tendencies
    CALL TEND_INIT( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in TENDENCIES_INIT', LOC ) 
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
    USE TENDENCIES_MOD,     ONLY : TEND_CLEANUP
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
    CALL TEND_CLEANUP()

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
    USE TRACER_MOD,         ONLY : XNUMOLAIR
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
!  24 Feb 2015 - C. Keller   - Initial version
!  19 Oct 2015 - C. Keller   - Rename AIRDEN to AIRDENSITY to
!                              avoid name conflict with GEOS-5 model.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: cId, Collection, N
    INTEGER            :: SpaceDim 
    REAL(hp)           :: ScaleFact
    CHARACTER(LEN=31)  :: OutOper, OutUnit
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGINIT_GRIDBOX (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGINIT_TRACER_CONC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Skip if ND68 diagnostic is turned off
    IF ( Input_Opt%ND68 <= 0 ) RETURN

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND68_OUTPUT_TYPE
      
    ! There are four diagnostics 
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

  END SUBROUTINE DIAGINIT_GRIDBOX
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

    ! Write new GEOS-Chem diagnostics
#if defined( DEVEL )
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
