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
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: DiagInit_Drydep
  PRIVATE :: DiagInit_Conv_Loss
  PRIVATE :: DiagInit_Wetdep_Loss
  PRIVATE :: DiagInit_Tracer_Conc
  PRIVATE :: DiagInit_Tracer_Emis
  PRIVATE :: DiagInit_GridBox
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
#endif
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
#if defined( DEVEL )
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
!    USE UCX_MOD,            ONLY : DIAGINIT_UCX
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: CollectionID
    REAL(sp)           :: TS
    REAL(fp), POINTER  :: AM2(:,:) => NULL()
    CHARACTER(LEN=15)  :: WriteFreq
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
#if defined(ESMF_)
    WriteFreq = 'Always' 
#else
    WriteFreq = 'Daily' 
#endif

    CALL DiagnCollection_Create( am_I_Root,                   &
                                 NX        = IIPAR,           &
                                 NY        = JJPAR,           &
                                 NZ        = LLPAR,           &
                                 TS        = TS,              &
                                 AM2       = AM2,             &
                                 PREFIX    = DGN,             &
                                 WriteFreq = TRIM(WriteFreq), &
                                 COL       = CollectionID,    &
                                 RC        = RC         )
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
    CALL DIAGINIT_WETDEP_LOSS( am_I_Root, Input_Opt, State_Met, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_WETDEP_LOSS', LOC ) 
    ENDIF

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

    ! Grid box quantities (ND68)
    CALL DIAGINIT_GRIDBOX( am_I_Root, Input_Opt, State_Met, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_GRIDBOX', LOC )
    ENDIF

    ! Tracer emission diagnostics (NEW)
    CALL DIAGINIT_TRACER_EMIS( am_I_Root, Input_Opt, State_Met, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DIAGINIT_TRACER_EMIS', LOC ) 
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
  SUBROUTINE DiagInit_Tracer_Conc( am_I_Root, Input_Opt, RC )
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: cId, Collection, N
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
    DO N = 1, 4 

       !----------------------------------------------------------------
       ! Create containers
       !----------------------------------------------------------------

       SELECT CASE ( N )
          CASE ( 1 )
             DiagnName = 'BOXHEIGHT'
             OutUnit   = 'm'
             ScaleFact = 1.0_hp
          CASE ( 2 )
             DiagnName = 'AIRMASS'
             OutUnit   = 'kg'
             ScaleFact = 1.0_hp
          CASE ( 3 )
             DiagnName = 'AIRDENS'
             OutUnit   = 'molecules air m-3'
             ScaleFact = XNUMOLAIR 
          CASE ( 4 )
             IF ( .NOT. ASSOCIATED(State_Met%AVGW) ) CYCLE
             DiagnName = 'AVGW'
             OutUnit   = 'v/v'
             ScaleFact = 1.0_hp
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
  SUBROUTINE DIAGINIT_WETDEP_LOSS( am_I_Root, Input_Opt, State_Met, RC )
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
    CHARACTER(LEN=255) :: LOC = 'DIAGINIT_WETDEP_LOSS (diagnostics_mod.F90)' 

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
    M = GET_WETDEP_NSOL()

    ! Loop over # of species 
    DO N = 1, M 

       ! Get GEOS-Chem tracer number
       NN = GET_WETDEP_IDWETD( N )

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
    ENDDO !N

  END SUBROUTINE DIAGINIT_WETDEP_LOSS
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
