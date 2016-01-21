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
  USE CMN_SIZE_Mod,       ONLY : IIPAR, JJPAR, LLPAR, NNPAR
  USE Error_Mod,          ONLY : Error_Stop
  USE HCO_Error_Mod
  USE GIGC_ErrCode_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Diagnostics_Write

  ! All of the following is currently only active for NETCDF=y
#if defined( NETCDF )
  PUBLIC  :: Diagnostics_Init
  PUBLIC  :: Diagnostics_Final
  PUBLIC  :: DiagnUpdate_NTracers_3D    ! utility subroutine
  PUBLIC  :: DiagnUpdate_Met            ! ND31, 55, 56, 57, 66, 67, 68, and 69
  PUBLIC  :: DiagnUpdate_Transport_Flux ! ND24, ND25, ND26
#if defined( DEVEL )
  PUBLIC  :: CalcDobsonColumn           ! added by ck - no existing bpch diag
#endif
!
! !PRIVATE MEMBER FUNCTIONS:
!
! Benchmark diagnostic groups
  !PRIVATE :: DiagnInit_Sulfate_ProdLoss ! ND05 (implementation not complete)
  !PRIVATE :: DiagnInit_Carbon_Sources   ! ND07 (implementation not complete)
  !PRIVATE :: DiagnInit_Cloud_Properties ! ND21 (implementation not complete)
  !PRIVATE :: DiagnInit_Photolysis_Rates ! ND22 (implementation not complete)
  PRIVATE :: DiagnInit_Transport_Flux    ! ND24, 25, and 26
  PRIVATE :: DiagnInit_Conv_Loss        ! ND38
  PRIVATE :: DiagnInit_Wetdep_Loss      ! ND39
  PRIVATE :: DiagnInit_Drydep           ! ND44
  PRIVATE :: DiagnInit_Tracer_Conc      ! ND45
  PRIVATE :: DiagnInit_Met              ! ND31, 55, 56, 57, 66, 67, 68, and 69

! Specialty simulation diagnostic groups
  PRIVATE :: DiagnInit_Pb_Emiss      ! ND01 init
  PRIVATE :: DiagnInit_Rn_Decay      ! ND02 init
  PRIVATE :: DiagnInit_CH4_Loss      ! ND19 init
  !   to do: ND03 (Hg emissions, P/L)
  !          ND43 (POPs emissions prod/loss)
  !          ND60 (wetland fraction)
  !          ND61 (TOMAS aerosol rates)
  !          ND61 (TOMAS 3D rates)
  !          ND72 (RRTMD radiative output)
  
! Optional diagnostic groups
  PRIVATE :: DiagnInit_BL_Frac          ! ND12 init
  PRIVATE :: DiagnInit_CldConv_Flx      ! ND14 init
  PRIVATE :: DiagnInit_BlMix_Flx        ! ND15 init
  PRIVATE :: DiagnInit_Precip_Frac      ! ND16 init
  PRIVATE :: DiagnInit_Rain_Frac        ! ND17 init
  PRIVATE :: DiagnInit_Wash_Frac        ! ND18 init
  PRIVATE :: DiagnInit_Landmap          ! ND30 init
  PRIVATE :: DiagnInit_Tracer_Emis      ! added by ck
  !PRIVATE :: DiagnInit_HrlyMax_SurfConc ! ND71 init (impementation commented
                                         ! out pending decision on if to keep)
  !   to do: ND37 (updraft scav fraction)
  !          ND62 (instantaneous column maps)

#if defined( DEVEL )
  PRIVATE :: DiagnInit_Dobson        ! added by ck
#endif
!
! !DEFINED PARAMETERS:
!
  ! Prefix of restart file. This file will hold all diagnostics that are 
  ! written out at the end of a simulation (either because their output 
  ! frequency is set to 'End' or the run finishes and these diagnostics
  ! haven't reached the end of their output interval yet).
  CHARACTER(LEN=31), PARAMETER :: RST = 'GEOSCHEM_Restart'
  CHARACTER(LEN=31), PARAMETER :: DGN = 'GEOSCHEM_Diagnostics_Hrly'

  ! Toggle to enable species diagnostics. This will write out species 
  ! concentrations (in addition to the tracers). Not recommended unless
  ! you have a good reason (ckeller, 8/11/2015).
  LOGICAL, PARAMETER, PUBLIC   :: DiagnSpec = .FALSE.
#endif
#if defined( NETCDF )
  ! Initialize GEOS-Chem diagnostics container id (ewl, 1/20/16)
  ! If we add containers, we will need to modify this to have a different
  ! container id variable per collection.
  INTEGER, SAVE, PRIVATE       :: cId = 0                
#endif
!
! !REVISION HISTORY:
!  09 Jan 2015 - C. Keller   - Initial version. 
!  14 Jan 2016 - E. Lundgren - Add several GEOS-Chem diagnostics
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
! diagnostics collections and populates them with diagnostics containers. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diagnostics_Init( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE GRID_MOD,           ONLY : AREA_M2
    USE HCO_DIAGN_MOD
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

    ! Temporarily manually set to hourly. Eventually this will be set
    ! per collection in input.geos or a new diagnostics input file
    !  (ewl, 1/13/16)
    CALL DiagnCollection_Create( am_I_Root,                      &
                                 NX           = IIPAR,           &
                                 NY           = JJPAR,           &
                                 NZ           = LLPAR,           &
                                 TS           = TS,              &
                                 AM2          = AM2,             &
                                 PREFIX       = DGN,             &
!                                 deltaYMD     = deltaYMD,        &
!                                 deltaHMS     = deltaHMS,        &
                                 deltaYMD     = 00000000,        &
                                 deltaHMS     = 010000,          &
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
    ! Create diagnostics containers and add to collections by initializing
    ! each diagnostic group
    ! (Add calls to additional subroutines for other diagnostics here!)
    !
    ! NOTE: Right now there is only one GEOS-Chem collection, but eventually
    ! we may want to add more (i.e. hourly, instantaneous, monthly, etc.)
    !-----------------------------------------------------------------------

!   ! Sulfate prod/loss diagnostic (ND05)
!   IF ( Input_Opt%ND05 > 0 ) THEN
!      CALL DiagnInit_Sulfate_ProdLoss( am_I_Root, Input_Opt, RC )
!      IF ( RC /= GIGC_SUCCESS ) THEN
!         CALL ERROR_STOP( 'Error in DiagnInit_Sulfate_ProdLoss', LOC ) 
!      ENDIF
!   ENDIF
!
!   ! Carbon sources diagnostic (ND07)
!   IF ( Input_Opt%ND07 > 0 ) THEN
!      CALL DiagnInit_Carbon_Sources( am_I_Root, Input_Opt, RC )
!      IF ( RC /= GIGC_SUCCESS ) THEN
!         CALL ERROR_STOP( 'Error in DiagnInit_Carbon_Sources', LOC ) 
!      ENDIF
!   ENDIF
!
!   ! Cloud properties (ND21)
!   IF ( Input_Opt%ND21 > 0 ) THEN
!      CALL DiagnInit_Cloud_Properties( am_I_Root, Input_Opt, RC )
!      IF ( RC /= GIGC_SUCCESS ) THEN
!         CALL ERROR_STOP( 'Error in DiagnInit_Cloud_Properties', LOC ) 
!      ENDIF
!   ENDIF
!
!   ! Photolysis rates diagnostic (ND22)
!   IF ( Input_Opt%ND22 > 0 ) THEN
!      CALL DiagnInit_Photolysis_Rates( am_I_Root, Input_Opt, RC )
!      IF ( RC /= GIGC_SUCCESS ) THEN
!         CALL ERROR_STOP( 'Error in DiagnInit_Photolysis_Rates', LOC ) 
!      ENDIF
!   ENDIF

    ! Transport fluxes diagnostic (ND24, ND25, and ND26)
    ! For now, use ND24 as indicator of whether to turn on diag group
    IF ( Input_Opt%ND24 > 0 ) THEN
       CALL DiagnInit_Transport_Flux( am_I_Root, Input_Opt, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Transport_Flux', LOC ) 
       ENDIF
    ENDIF

    ! Convective scavenging loss (ND38)
    IF ( Input_Opt%ND38 > 0 ) THEN
       CALL DiagnInit_Conv_Loss( am_I_Root, Input_Opt, State_Met, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Conv_Loss', LOC ) 
       ENDIF
    ENDIF

    ! Wetdeposition scavenging loss (ND39)
    IF ( Input_Opt%ND39 > 0 ) THEN
       CALL DiagnInit_WetDep_Loss( am_I_Root, Input_Opt, State_Met,  &
                                  State_Chm, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_WetDep_Loss', LOC ) 
       ENDIF
    ENDIF

    ! Drydep diagnostic (ND44)
    IF ( Input_Opt%ND44 > 0 ) THEN
       CALL DiagnInit_DryDep( am_I_Root, Input_Opt, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_DryDep', LOC ) 
       ENDIF
    ENDIF

    ! Tracer concentration diagnostics (ND45)
    IF ( Input_Opt%ND45 > 0 ) THEN
       CALL DiagnInit_Tracer_Conc( am_I_Root, Input_Opt, State_Chm, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Tracer_Conc', LOC ) 
       ENDIF
    ENDIF

    ! Meteorology state diagnostic (ND31, 55, 56, 57, 66, 67, 68, and 69)
    ! For now, use ND68 as indicator for entire diagnostic group
    IF ( Input_Opt%ND68 > 0 ) THEN
       CALL DiagnInit_Met( am_I_Root, Input_Opt, State_Met, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Met', LOC )
       ENDIF
    ENDIF

    ! Pb emissions diagnostic (ND01)
    IF ( Input_Opt%ND01 > 0 ) THEN
       CALL DiagnInit_Pb_Emiss( am_I_Root, Input_Opt, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Pb_Emiss', LOC ) 
       ENDIF
    ENDIF

    ! Rn/Pb/Be decay diagnostic (ND02)
    IF ( Input_Opt%ND02 > 0 ) THEN
       CALL DiagnInit_Rn_Decay( am_I_Root, Input_Opt, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Rn_Decay', LOC ) 
       ENDIF
    ENDIF

    ! Boundary layer fraction diagnostic (ND12)
    IF ( Input_Opt%ND12 > 0 ) THEN
       CALL DiagnInit_BL_Frac( am_I_Root, Input_Opt, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_BL_Frac', LOC ) 
       ENDIF
    ENDIF

    ! Cloud convection mass flux diagnostic (ND14)
    IF ( Input_Opt%ND14 > 0 ) THEN
       CALL DiagnInit_CldConv_Flx( am_I_Root, Input_Opt, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_CldConv_Flx', LOC ) 
       ENDIF
    ENDIF

    ! Boundary-layer mixing mass flux diagnostic (ND15)
    IF ( Input_Opt%ND15 > 0 ) THEN
       CALL DiagnInit_BLMix_Flx( am_I_Root, Input_Opt, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_BLMix_Flx', LOC ) 
       ENDIF
    ENDIF

    ! Areal fraction of precip diagnostic (ND16)
    IF ( Input_Opt%ND16 > 0 ) THEN
       CALL DiagnInit_Precip_Frac( am_I_Root, Input_Opt, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Precip_Frac', LOC ) 
       ENDIF
    ENDIF

    ! Rainout fraction diagnostic (ND17)
    IF ( Input_Opt%ND17 > 0 ) THEN
       CALL DiagnInit_Rain_Frac( am_I_Root, Input_Opt, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Rain_Frac', LOC ) 
       ENDIF
    ENDIF

    ! Washout fraction diagnostic (ND18)
    IF ( Input_Opt%ND18 > 0 ) THEN
       CALL DiagnInit_Wash_Frac( am_I_Root, Input_Opt, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Wash_Frac', LOC ) 
       ENDIF
    ENDIF

    ! CH4 loss diagnostic (ND19)
    IF ( Input_Opt%ND19 > 0 ) THEN
       CALL DiagnInit_CH4_Loss(am_I_Root, Input_Opt, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_CH4_Loss', LOC ) 
       ENDIF
    ENDIF


    ! Land map diagnostic (ND30)
    IF ( Input_Opt%ND30 > 0 ) THEN
       CALL DiagnInit_LandMap( am_I_Root, Input_Opt, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_LandMap', LOC ) 
       ENDIF
    ENDIF

    ! UCX diagnostics are now initialized in ucx_mod.F. The UCX diagnostics
    ! currently only include the PSC state, which is written into the HEMCO
    ! restart file for now. Initialize this outside this module to make sure
    ! that this field is also diagnosed if we are not using the DEVEL compiler
    ! switch (ckeller, 3/25/2015). 
!    ! UCX diagnostics
!    IF ( Input_Opt%LUCX ) THEN
!       CALL DiagnInit_UCX( am_I_Root, Input_Opt, RC )
!       IF ( RC /= GIGC_SUCCESS ) THEN
!          CALL ERROR_STOP( 'Error in DiagnInit_UCX', LOC ) 
!       ENDIF
!    ENDIF

    !! Houly-maximum tracer mixing ratio (IJ-MAX) at surface (ND71)
    !IF ( Input_Opt%ND71 > 0 ) THEN
    !   CALL DiagnInit_HrlyMax_SurfConc( am_I_Root, Input_Opt, RC )
    !   IF ( RC /= GIGC_SUCCESS ) THEN
    !      CALL ERROR_STOP( 'Error in DiagnInit_HrlyMax_SurfConc', LOC )
    !   ENDIF
    !ENDIF

    ! Tracer emission diagnostics (NEW) (added by Christoph)
    ! NOTE: Currently this diagnostic must be initialized last since since
    ! container ids start at 10000 for routine TotalsToLogFile (ewl, 1/20/16)
    IF ( .FALSE. ) THEN
       CALL DiagnInit_Tracer_Emis( am_I_Root, Input_Opt, State_Met, RC )
       IF ( RC /= GIGC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Tracer_Emis', LOC ) 
       ENDIF
    ENDIF

#if defined( DEVEL )
    CALL DiagnInit_Dobson( am_I_Root, Input_Opt, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DiagnInit_Dobson', LOC ) 
    ENDIF

    ! Initialize tendencies
    CALL Tend_Init( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in Tendencies_Init', LOC ) 
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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
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
! !IROUTINE: diagninit_drydep
!
! !DESCRIPTION: Subroutine DIAGNINIT\_DRYDEP initializes the dry deposition 
! diagnostics arrays. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_DryDep( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE TRACERID_Mod,       ONLY : IDTISOPN, IDTMMN
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
    INTEGER            :: Collection, D, N, M
    LOGICAL            :: First_ISOPN, First_MMN
    CHARACTER(LEN=15)  :: OutOper
    CHARACTER(LEN=60)  :: DiagnName, DryDepName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DiagnInit_DruDep (diagnostics_mod.F)' 
    
    !=======================================================================
    ! DIAGNINIT_DRYDEP begins here!
    !=======================================================================

    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND44_OUTPUT_TYPE

    ! Special cases if ISOPN or MMN since there are two 
    ! corresponding dry deposition species each, ISOPND and ISOPNB for
    ! ISOPN, and MACRN and MAVKN for MMN. Therefore, track when ISOPN and
    ! MMN are first encountered. (ewl, 1/20/16)
    First_ISOPN = .TRUE.
    First_MMN = .TRUE.

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

          IF ( N == IDTISOPN .AND. First_ISOPN ) THEN
             DryDepName = 'ISOPND'
             First_ISOPN = .FALSE.
          ELSE IF ( N == IDTISOPN .AND. .NOT. First_ISOPN ) THEN
             DryDepName = 'ISOPNB'
          ELSE IF ( N == IDTMMN .AND. First_MMN ) THEN
             DryDepName = 'MACRN'
             First_MMN = .FALSE.
          ELSE IF ( N == IDTMMN .AND. .NOT. First_MMN ) THEN
             DryDepName = 'MAVKN'
          ELSE
             DryDepName =  Input_Opt%TRACER_NAME(N)
          ENDIF
             
          ! Diagnostic container name and id
          DiagnName = 'DRYDEP_VEL_' // TRIM( DryDepName )
          cId = cId + 1
   
          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        &
                             cId       = cId,               & 
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
   
          ! Diagnostic container name and id
          DiagnName = 'DRYDEP_FLX_' // TRIM( DryDepName )
          cId       = cId + 1
   
          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        &
                             cId       = cId,               &
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

  END SUBROUTINE DiagnInit_DryDep
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_tracer_conc
!
! !DESCRIPTION: Subroutine DIAGNINIT\_TRACER\_CONC initializes the tracer
!  concentration diagnostic (aka ND45).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_Tracer_Conc( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE COMODE_LOOP_MOD,    ONLY : NAMEGAS, NTSPEC, NCSURBAN
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create

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
!  13 Jan 2016 - E. Lundgren - Define diagnostic ID
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
    CHARACTER(LEN=255) :: LOC = 'DIAGNINIT_TRACER_CONC (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_TRACER_CONC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    !----------------------------------------------------------------
    ! Create containers for drydep velocity [m/s]
    !----------------------------------------------------------------
    
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
    
          ! Diagnostic container name and id
          DiagnName = 'TRACER_CONC_' // TRIM( Input_Opt%TRACER_NAME(N) )
          cId = cId + 1
    
          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cId       = cId,               &
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
   
          ! Diagnostic container name and id
          IF ( TRIM(NAMEGAS(N)) == '' ) CYCLE
          DiagnName = 'SPECIES_CONC_' // TRIM(NAMEGAS(N)) 
          cId = cId + 1
  
          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        &
                             cId       = cId,               &
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
   
   END SUBROUTINE DiagnInit_Tracer_Conc
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_pb_emiss
!
! !DESCRIPTION: Subroutine DIAGNINIT\_PB\_EMISS initializes the Pb emissions 
!  diagnostic (aka ND01). Other ND01 tracer emissions (Rn and Be7) are 
!  handled within HEMCO.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_Pb_Emiss( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create

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
    INTEGER            :: Collection
    CHARACTER(LEN=15)  :: OutOper
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGNINIT_PB_EMISS (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_PB_EMISS begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND01_OUTPUT_TYPE
    
    ! If the tracer number for lead is scheduled for output in input.geos, 
    ! then define the diagnostic container for 210Pb emissions.
    IF ( ANY ( Input_Opt%TINDEX(1,:) == IDTPB ) ) THEN

       !----------------------------------------------------------------
       ! Create containers for Pb emissions [kg/s]
       !----------------------------------------------------------------

       ! Diagnostic container name and id
       DiagnName = 'EMISS_' // TRIM( Input_Opt%TRACER_NAME( IDTPB ) )
       cId = cId + 1

       ! Create container
       CALL Diagn_Create( am_I_Root,                     &
                          Col       = Collection,        & 
                          cId       = cId,               &
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
                          RC        = RC )

       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF
    ENDIF

  END SUBROUTINE DiagnInit_Pb_Emiss
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_rn_decay
!
! !DESCRIPTION: Subroutine DIAGNINIT\_RN\_DECAY initializes the Rn/Pb/Be7
!  decay diagnostic (aka ND02).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_Rn_Decay( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
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
    INTEGER              :: Collection, M
    INTEGER, PARAMETER   :: NumTracers = 3        ! Does a var exist for this?
    INTEGER              :: TracersN ( NumTracers ) 
    CHARACTER(LEN=15)    :: OutOper
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGNINIT_RN_DECAY (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_RN_DECAY begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

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

          ! Diagnostic container name and id
          DiagnName = 'DECAY_' //                           &
                      TRIM( Input_Opt%TRACER_NAME( TracersN( M ) ) )
          cId = cId + 1
          
          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cId       = cId,               &
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
                             RC        = RC )

          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
       ENDIF
    ENDDO

  END SUBROUTINE DiagnInit_Rn_Decay
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_bl_frac (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGNINIT\_BL\_FRAC initializes the distribution 
!  of surface emissions in the boundary layer diagnostics (aka ND12).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_BL_Frac( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create

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
    INTEGER              :: Collection, N
    CHARACTER(LEN=15)    :: OutOper
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGNINIT_BL_FRAC (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_BL_FRAC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND12_OUTPUT_TYPE

    ! Check if certain tracer(s) listed for ND12 in input.geos???

    !----------------------------------------------------------------
    ! Create containers for fraction of BL occupied by level L [.]
    !----------------------------------------------------------------

    ! Diagnostic container name and id
    DiagnName = 'BL_FRAC'
    cId = cId + 1

    ! Create container
    CALL Diagn_Create( am_I_Root,                     &
                       Col       = Collection,        & 
                       cId       = cId,               &
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
                       RC        = RC )

    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP( MSG, LOC ) 
    ENDIF

  END SUBROUTINE DiagnInit_BL_Frac
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_cldconv_flx (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGNINIT\_CLDCONV\_FLX initializes the upward
!  mass flux due to wet convection diagnostic (aka ND14).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_CldConv_Flx( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create

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
    INTEGER              :: Collection, N
    CHARACTER(LEN=15)    :: OutOper
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGNINIT_CLDCONV_FLX (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_CLDCONV_FLX begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

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
      
          ! Diagnostic container name and id
          DiagnName = 'CLDCONV_FLX_' // TRIM( Input_Opt%TRACER_NAME(N) )
          cId = cId + 1

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cId       = cId,               &
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
                             RC        = RC )
      
          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF      
       ENDIF
    ENDDO   

  END SUBROUTINE DiagnInit_CldConv_Flx
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_blmix_flx (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGNINIT\_BLMIX\_FLX initializes the upward
!  mass flux from boundary-layer mixing diagnostic (aka ND15).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_BLMix_Flx( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create

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
    INTEGER              :: Collection, N
    CHARACTER(LEN=15)    :: OutOper
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGNINIT_BLMIX_FLX (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_BLMIX_FLX begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

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
      
          ! Diagnostic container name and id
          DiagnName = 'BLMIX_FLX_' // TRIM( Input_Opt%TRACER_NAME(N) )
          cId = cId + 1

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cId       = cId,               &
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
                             RC        = RC )
      
          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
       ENDIF
    ENDDO   

  END SUBROUTINE DiagnInit_BLMix_Flx
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_precip_frac (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGNINIT\_PRECIP\_FRAC initializes the areal
!  fraction of precipitation diagnostic (aka ND16).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_Precip_Frac( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create

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
    INTEGER              :: Collection, M
    CHARACTER(LEN=15)    :: OutOper
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGNINIT_PRECIP_FRAC (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_PRECIP_FRAC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

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

       ! Diagnostic container id
       cId = cId + 1

       ! Create container
       CALL Diagn_Create( am_I_Root,                     &
                          Col       = Collection,        & 
                          cId       = cId,               &
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
                          RC        = RC )
   
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot create ND16 diagnostic: ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF
   ENDDO

  END SUBROUTINE DiagnInit_Precip_Frac
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_rain_frac (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGNINIT\_RAIN\_FRAC initializes the rainout
!  fraction in precipitation diagnostic (aka ND17).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_Rain_Frac( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create

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
    INTEGER              :: Collection, N, M
    CHARACTER(LEN=15)    :: OutOper
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGNINIT_RAIN_FRAC (diagnostics_mod.F90)' !
    !=======================================================================
    ! DIAGNINIT_RAIN_FRAC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

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

             ! Diagnostic container id
             cId = cId + 1
      
             ! Create container
             CALL Diagn_Create( am_I_Root,                     &
                                Col       = Collection,        & 
                                cId       = cId,               &
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
                                RC        = RC )
         
             IF ( RC /= HCO_SUCCESS ) THEN
                MSG = 'Cannot create ND17 diagnostic: ' // TRIM(DiagnName)
                CALL ERROR_STOP( MSG, LOC ) 
             ENDIF 
          ENDDO     
       ENDIF
    ENDDO     

  END SUBROUTINE DiagnInit_Rain_Frac
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_wash_frac (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGNINIT\_WASH\_FRAC initializes the washout
!  fraction diagnostic (aka ND18).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_Wash_Frac( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create

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
    INTEGER              :: Collection, N, M
    CHARACTER(LEN=15)    :: OutOper
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGNINIT_WASH_FRAC (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_WASH_FRAC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

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

             ! Diagnostic container id
             cId = cId + 1

             ! Create container
             CALL Diagn_Create( am_I_Root,                     &
                                Col       = Collection,        & 
                                cId       = cId,               &
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
                                RC        = RC )
         
             IF ( RC /= HCO_SUCCESS ) THEN
                MSG = 'Cannot create ND18 diagnostic: ' // TRIM(DiagnName)
                CALL ERROR_STOP( MSG, LOC ) 
             ENDIF
          ENDDO
       ENDIF
    ENDDO

  END SUBROUTINE DiagnInit_Wash_Frac
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_CH4_Loss (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGNINIT\_CH4\_LOSS initializes the methane
!  loss diagnostic (aka ND19).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_CH4_Loss( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create

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
    INTEGER              :: Collection, N
    CHARACTER(LEN=15)    :: OutOper
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGNINIT_CH4_LOSS (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_CH4_LOSS begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND19_OUTPUT_TYPE
    
    ! Check if certain tracer(s) listed for ND19 in input.geos???

    !----------------------------------------------------------------
    ! Create containers for CH4 removal by OH [kg CH4]
    !----------------------------------------------------------------

    ! Diagnostic container name and id
    DiagnName = 'CH4_LOSS'
    cId = cId + 1

    ! Create container
    CALL Diagn_Create( am_I_Root,                     &
                       Col       = Collection,        & 
                       cId       = cId,               &
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
                       RC        = RC )

    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP( MSG, LOC ) 
    ENDIF

  END SUBROUTINE DiagnInit_CH4_Loss
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_transport_flux
!
! !DESCRIPTION: Subroutine DIAGNINIT\_TRANSPORT\_FLUX initializes the zonal
!  (east/west), meridional (north/south), and vertical mass transport 
!  flux diagnostics (aka ND24, ND25, and ND26).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_Transport_Flux( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create

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
!  19 Jan 2016 - E. Lundgren - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER              :: Collection, N, M
    CHARACTER(LEN=15)    :: OutOper
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=30)    :: NamePrefix
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DiagnInit_Transport_Flux (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_TRANSPORT_FLUX begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    ! For now, use ND24 for entire transport diagnostic group
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND24_OUTPUT_TYPE

    ! Loop over 3 types of transport fluxes
    DO M = 1, 3

       ! Loop over tracers
       DO N = 1, Input_Opt%N_TRACERS
         
          ! If this tracer number N is scheduled for output
          ! then define the diagnostic containers for mass flux
          ! NOTE: for now, use tracers indicated for ND24 only until 
          ! consolidation into one group in new netcdf diag input file
          IF ( ANY( Input_Opt%TINDEX(24,:) == N ) ) THEN

             !----------------------------------------------------------------
             ! Create container for mass flux by transport [kg/s]
             !----------------------------------------------------------------

             SELECT CASE ( M )
                CASE ( 1 )
                   NamePrefix = 'TRANSPORT_FLX_EW_'   ! ND24
                CASE ( 2 )
                   NamePrefix = 'TRANSPORT_FLX_NS_'   ! ND25
                CASE ( 3 )
                   NamePrefix = 'TRANSPORT_FLX_VERT_' ! ND26
             END SELECT

             ! Diagnostic container name and ID
             DiagnName = TRIM( NamePrefix ) // TRIM( Input_Opt%TRACER_NAME(N) )
             cId = cId + 1
   
             ! Create container
             CALL Diagn_Create( am_I_Root,                     &
                                Col       = Collection,        & 
                                cId       = cId,               &
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
                                RC        = RC )
         
             IF ( RC /= HCO_SUCCESS ) THEN
                MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
                CALL ERROR_STOP( MSG, LOC ) 
             ENDIF  
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE DiagnInit_Transport_Flux
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_landmap (LL in progress)
!
! !DESCRIPTION: Subroutine DIAGNINIT\_LANDMAP initializes the land map  
!  diagnostic (aka ND30).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_LandMap( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
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
    INTEGER              :: Collection, N
    CHARACTER(LEN=15)    :: OutOper 
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGNINIT_LANDMAP (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_LANDMAP begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND30_OUTPUT_TYPE
   
    !----------------------------------------------------------------
    ! Create container for the GMAO land-water indices [.]
    !----------------------------------------------------------------

    ! Diagnostic container name and id
    DiagnName = 'LANDMAP'
    cId = cId + 1

    ! Create container
    CALL Diagn_Create( am_I_Root,                     &
                       Col       = Collection,        & 
                       cId       = cId,               &
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
                       RC        = RC )

    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP( MSG, LOC ) 
    ENDIF   

  END SUBROUTINE DiagnInit_LandMap

!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_met
!
! !DESCRIPTION: Subroutine DIAGNINIT\_MET initializes the meteorology state
!  diagnostics (aka ND31, ND55, ND56, ND57, ND66, ND67, ND68, and ND69).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_Met( am_I_Root, Input_Opt, State_Met, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE PHYSCONSTANTS,      ONLY : XNUMOLAIR
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
!  15 Jan 2016 - E. Lundgren - Revise for all MET-related diagnostics
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: Collection, N, Num2D, Num3D
    INTEGER            :: SpaceDim 
    REAL(hp)           :: ScaleFact
    CHARACTER(LEN=15)  :: OutOper, OutUnit
    CHARACTER(LEN=30)  :: NameSuffix, DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DiagnInit_Met (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_MET begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND68_OUTPUT_TYPE

    ! Set number of 3D and 2D MET diagnostics
    Num3d = 12
    Num2d = 28

    !----------------------------------------------------------------
    ! Create containers
    !----------------------------------------------------------------
  
    ! Loop over 3D diagnostics within Met category
    DO N = 1, Num3D

       SELECT CASE ( N )
          CASE ( 1 )
             NameSuffix = 'THETA'           ! ND57, trcr 1
             OutUnit    = 'K'
          CASE ( 2 )
             NameSuffix = 'UWND'            ! ND66, trcr 1 
             OutUnit    = 'm/s'
          CASE ( 3 )
             NameSuffix = 'VWND'            ! ND66, trcr 2
             OutUnit    = 'm/s'
          CASE ( 4 )
             NameSuffix = 'TMPU'            ! ND66, trcr 3
             OutUnit    = 'K'
          CASE ( 5 )
             NameSuffix = 'SPHU'            ! ND66, trcr 4
             OutUnit    = 'g H2O/kg air'
          CASE ( 6 )
             NameSuffix = 'CLDMAS'          ! ND66, trcr 5
             OutUnit    = 'kg/m2/s'
          CASE ( 7 )
             NameSuffix = 'DTRAIN'          ! ND66, trcr 6
             OutUnit    = 'kg/m2/s'
          CASE ( 8 )
             NameSuffix = 'BXHEIGHT'        ! ND68, trcr 1
             OutUnit    = 'm'
          CASE ( 9 )
             NameSuffix = 'DRYAIRMASS'      ! ND68, trcr 2
             OutUnit    = 'kg'
          CASE ( 10 )
             IF ( .NOT. ASSOCIATED(State_Met%AVGW) ) CYCLE
             NameSuffix = 'AVGW'            ! ND68, trcr 3
             OutUnit    = 'v/v'
          CASE ( 11 )
             NameSuffix = 'NAIR'            ! ND68, trcr 4
             OutUnit    = 'molec dry air/m3'
          CASE ( 12 )
             NameSuffix = 'PEDGE'           ! ND31, trcr 1
             OutUnit    = 'hPa'
          CASE DEFAULT
             IF ( N < Num3D ) THEN
                MSG = 'Num3D is less than number of named 3D MET diagnostics'
             ELSE
                MSG = 'Undefined 3D diagnostic in MET diagnostic group'
             ENDIF
             CALL ERROR_STOP( MSG, LOC ) 
       END SELECT

       ! Diagnostic container info
       DiagnName = 'MET_' // NameSuffix
       cId       = cId + 1
       SpaceDim  = 3

       ! Create container
       CALL Diagn_Create( am_I_Root,                     &
                          Col       = Collection,        & 
                          cId       = cId,               &
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
          MSG = 'Cannot create MET diagnostic ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF
    ENDDO

    ! Loop over 2D diagnostics within Met category
    DO N = 1, Num2D

       SELECT CASE ( N )
          CASE ( 1 )
             NameSuffix = 'TR-PAUSE_Level'  ! ND55, trcr 1
             OutUnit    = 'Level index'
          CASE ( 2 )
             NameSuffix = 'TR-PAUSE_km'     ! ND55, trcr 2
             OutUnit    = 'km'
          CASE ( 3 )
             NameSuffix = 'TR-PAUSE_mb'     ! ND55, trcr 3
             OutUnit    = 'mb'
          CASE ( 4 )
             NameSuffix = 'LFLASH_TOTAL'    ! ND56, trcr 1
             OutUnit    = 'flashes/min/km2'
          CASE ( 5 )
             NameSuffix = 'LFLASH_INTRA'    ! ND56, trcr 2
             OutUnit    = 'flashes/min/km2'
          CASE ( 6 )
             NameSuffix = 'LFLASH_CLDGND'  ! ND56, trcr 3
             OutUnit    = 'flashes/min/km2'
          CASE ( 7 )
             NameSuffix = 'HFLUX'           ! ND67, trcr 1
             OutUnit    = 'W/m2'
          CASE ( 8 )
             NameSuffix = 'RADSWG'          ! ND67, trcr 2
             OutUnit    = 'W/m2'
          CASE ( 9 )
             NameSuffix = 'PREACC'          ! ND67, trcr 3
             OutUnit    = 'mm/day'
          CASE ( 10 )
             NameSuffix = 'PRECON'          ! ND67, trcr 4
             OutUnit    = 'mm/day'
          CASE ( 11 )
             NameSuffix = 'TS'              ! ND67, trcr 5
             OutUnit    = 'K'
          CASE ( 12 )
             NameSuffix = 'RADSWT'          ! ND67, trcr 6
             OutUnit    = 'W/m2'
          CASE ( 13 )
             NameSuffix = 'USTAR'           ! ND67, trcr 7
             OutUnit    = 'm/s'
          CASE ( 14 )
             NameSuffix = 'Z0'              ! ND67, trcr 8
             OutUnit    = 'm'
          CASE ( 15 )
             NameSuffix = 'PBL'             ! ND67, trcr 9
             OutUnit    = 'hPa'
          CASE ( 16 )
             NameSuffix = 'CLDFRC'          ! ND67, trcr 10
             OutUnit    = '0-1'
          CASE ( 17 )
             NameSuffix = 'U10M'            ! ND67, trcr 11
             OutUnit    = 'm/s'
          CASE ( 18 )
             NameSuffix = 'V10M'            ! ND67, trcr 12
             OutUnit    = 'm/s'
          CASE ( 19 )
             NameSuffix = 'PS-PBL'          ! ND67, trcr 13
             OutUnit    = 'hPa'
          CASE ( 20 )
             NameSuffix = 'ALBD'            ! ND67, trcr 14
             OutUnit    = 'unitless'
          CASE ( 21 )
             NameSuffix = 'PHIS'            ! ND67, trcr 15
             OutUnit    = 'm'
          CASE ( 22 )
             NameSuffix = 'CLTOP'           ! ND67, trcr 16
             OutUnit    = 'levels'
          CASE ( 23 )
             NameSuffix = 'TROPP'           ! ND67, trcr 17
             OutUnit    = 'hPa'
          CASE ( 24 )
             NameSuffix = 'SLP'             ! ND67, trcr 18
             OutUnit    = 'hPa'
          CASE ( 25 )
             NameSuffix = 'TSKIN'           ! ND67, trcr 19
             OutUnit    = 'K'
          CASE ( 26 )
             NameSuffix = 'PARDF'           ! ND67, trcr 20
             OutUnit    = 'W/2'
          CASE ( 27 )
             NameSuffix = 'PARDR'           ! ND67, trcr 21
             OutUnit    = 'W/m2'
          CASE ( 28 )
             NameSuffix = 'GWETTOP'         ! ND67, trcr 22
             OutUnit    = 'unitless'
          CASE DEFAULT
             IF ( N < Num2D ) THEN
                MSG = 'Num2D is less than number of named 2D MET diagnostics'
             ELSE
                MSG = 'Undefined 2D diagnostic in MET diagnostic group'
             ENDIF
             CALL ERROR_STOP( MSG, LOC ) 
       END SELECT

       ! Diagnostic container info
       DiagnName = 'MET_' // NameSuffix
       cId       = cId + 1
       SpaceDim  = 2

       ! Create container
       CALL Diagn_Create( am_I_Root,                     &
                          Col       = Collection,        & 
                          cId       = cId,               &
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
          MSG = 'Cannot create MET diagnostic ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF
    ENDDO

  END SUBROUTINE DiagnInit_Met
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DiagnInit_Tracer_Emis
!
! !DESCRIPTION: Subroutine DiagnInit\_Tracer\_Emis initializes diagnostics for 
!  total species emissions diagnostics. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_Tracer_Emis( am_I_Root, Input_Opt, State_Met, RC ) 
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
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
    INTEGER            :: Collection, N
    CHARACTER(LEN=15)  :: OutOper
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGNINIT_TRACER_EMIS (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_TRACER_EMIS begins here!
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
       ID = GetHcoID( TrcId = N )
 
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

          ! Diagnostic container name and id
          DiagnName = 'TRACER_EMIS_' // TRIM( Input_Opt%TRACER_NAME(N) )
          cId = 10000 + ID

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cId       = cId,               &
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

  END SUBROUTINE DiagnInit_Tracer_Emis
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_conv_loss
!
! !DESCRIPTION: Subroutine DIAGNINIT\_CONV\_LOSS initializes the convective 
!  scavenging loss diagnostic (aka ND38). For now, this is a 2D column 
!  diagnostics, i.e. we only write out the loss due to convection of the
!  entire column.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DIAGNINIT_CONV_LOSS( am_I_Root, Input_Opt, State_Met, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
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
    INTEGER            :: Collection
    CHARACTER(LEN=15)  :: OutOper, OutUnit
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGNINIT_CONV_LOSS (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_CONV_LOSS begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

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
          cId = cId + 1

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cId       = cId,               &
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

  END SUBROUTINE DIAGNINIT_CONV_LOSS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_wetdep_loss
!
! !DESCRIPTION: Subroutine DIAGNINIT\_WETDEP\_LOSS initializes the wet 
!  deposition scavenging loss diagnostic (aka ND39).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DIAGNINIT_WETDEP_LOSS( am_I_Root, Input_Opt,     &
                                   State_Met, State_Chm, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
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
    INTEGER            :: Collection
    CHARACTER(LEN=15)  :: OutOper, OutUnit
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGNINIT_WETDEP_LOSS (diagnostics_mod.F90)' 

    ! Pointers
    TYPE(Species), POINTER :: ThisSpc

    !=======================================================================
    ! DIAGNINIT_WETDEP_LOSS begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

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
          cId = cId + 1

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             Col       = Collection,        & 
                             cId       = cId,               &
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

  END SUBROUTINE DIAGNINIT_WETDEP_LOSS
!EOC
#if defined( DEVEL )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_dobson
!
! !DESCRIPTION: Subroutine DIAGNINIT\_DOBSON initializes the O3 dobson column
! diagnostics. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_Dobson( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
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
    CHARACTER(LEN=255) :: LOC = 'DIAGNINIT_DOBSON (diagnostics_mod.F)' 
    
    !=======================================================================
    ! DIAGNINIT_DOBSON begins here!
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
 
       ! Diagnostic container name and id
       IF ( N == 1 ) THEN
          DiagnName = 'O3_COLUMN' 
       ELSEIF ( N == 2 ) THEN
          DiagnName = 'O3_TROPCOLUMN' 
       ENDIF
       cId = cId + 1

       ! Create container
       CALL Diagn_Create( am_I_Root,                     &
                          Col       = Collection,        & 
                          cId       = cId,               & 
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
   
  END SUBROUTINE DiagnInit_Dobson
!EOC
#endif
#if defined( DEVEL )
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
    USE PHYSCONSTANTS,      ONLY : AIRMW, AVO,   g0
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE HCO_Diagn_Mod,      ONLY : Diagn_Update
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
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagnupdate_ntracers_3D
!
! !DESCRIPTION: Subroutine DIAGNUPDATE\_NTRACERS\_3D updates a generic set of
!  3D diagnostics that are distinguished from each other by chemical tracer.
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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_Update
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
    INTEGER              :: Collection, N
    CHARACTER(LEN=15)    :: OutOper 
    CHARACTER(LEN=255)   :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DiagnUpdate_NTracers_3D (diagnostics_mod.F90)' 
    REAL(fp), TARGET     :: DiagnArray_tracer(IIPAR,JJPAR,LLPAR)
    REAL(fp), POINTER    :: Ptr3D(:,:,:)

    !=======================================================================
    ! DIAGNUPDATE_NTRACERS_3D begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

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

          ! Free the pointer before error handling
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
! !IROUTINE: diagnupdate_met
!
! !DESCRIPTION: Subroutine DIAGNUPDATE\_MET updates the meteorology state
!  diagnostics (aka ND31, ND55, ND56, ND57, ND66, ND67, ND68, and ND69).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnUpdate_Met( am_I_Root, Input_Opt, State_Met, RC )
!
! !USES:
!
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE HCO_Diagn_Mod,      ONLY : Diagn_Update
    USE PHYSCONSTANTS,      ONLY : XNUMOLAIR
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
!  15 Jan 2016 - E. Lundgren - Revise for all MET-related diagnostics
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: Collection, N, Num2D, Num3D
    CHARACTER(LEN=30)  :: NameSuffix, DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DiagnUpdate_Met (diagnostics_mod.F90)'
    REAL(fp), POINTER  :: Ptr3D(:,:,:) => NULL()
    REAL(fp), POINTER  :: Ptr2D(:,:)   => NULL()
    ! Need to define Arr2D and Arr3D

    !=======================================================================
    ! DIAGNUPDATE_MET begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    ! This is not currently used, but keep for possible later use/editing
    Collection = Input_Opt%DIAG_COLLECTION

    ! Set number of 3D and 2D MET diagnostics
    Num3d = 12
    Num2d = 28

    !----------------------------------------------------------------
    ! Update containers
    !----------------------------------------------------------------
  
    ! Loop over 3D diagnostics within Met category
    DO N = 1, Num3D

       SELECT CASE ( N )
          CASE ( 1 )
             NameSuffix = 'THETA'           ! ND57, trcr 1
             Ptr3D => State_Met%T  ! placeholder
          CASE ( 2 )
             NameSuffix = 'UWND'            ! ND66, trcr 1 
             Ptr3D => State_Met%T  ! placeholder
          CASE ( 3 )
             NameSuffix = 'VWND'            ! ND66, trcr 2
             Ptr3D => State_Met%T  ! placeholder
          CASE ( 4 )
             NameSuffix = 'TMPU'            ! ND66, trcr 3
             Ptr3D => State_Met%T  ! placeholder
          CASE ( 5 )
             NameSuffix = 'SPHU'            ! ND66, trcr 4
             Ptr3D => State_Met%T  ! placeholder
          CASE ( 6 )
             NameSuffix = 'CLDMAS'          ! ND66, trcr 5
             Ptr3D => State_Met%T  ! placeholder
          CASE ( 7 )
             NameSuffix = 'DTRAIN'          ! ND66, trcr 6
             Ptr3D => State_Met%T  ! placeholder
          CASE ( 8 )
             NameSuffix = 'BXHEIGHT'        ! ND68, trcr 1
             Ptr3D => State_Met%BXHEIGHT 
          CASE ( 9 )
             NameSuffix = 'DRYAIRMASS'      ! ND68, trcr 2
             Ptr3D => State_Met%AD
          CASE ( 10 )
             IF ( .NOT. ASSOCIATED(State_Met%AVGW) ) CYCLE
             NameSuffix = 'AVGW'            ! ND68, trcr 3
             Ptr3D => State_Met%AVGW
          CASE ( 11 )
             NameSuffix = 'NAIR'            ! ND68, trcr 4
!             Arr3D = State_Met%AIRDEN * XNUMOLAIR
             Ptr3D => State_Met%T  ! placeholder
          CASE ( 12 )
             NameSuffix = 'PEDGE'           ! ND31, trcr 1
             Ptr3D => State_Met%PEDGE
          CASE DEFAULT
             IF ( N < Num3D ) THEN
                MSG = 'Num3D is less than number of named 3D MET diagnostics'
             ELSE
                MSG = 'Undefined 3D diagnostic in MET diagnostic group'
             ENDIF
             CALL ERROR_STOP( MSG, LOC ) 
       END SELECT

       ! Diagnostic name
       DiagnName = 'MET_' // NameSuffix

       ! Update diagnostics
       IF ( ASSOCIATED(Ptr3D) ) THEN 
          CALL Diagn_Update( am_I_Root,                    &
                             cName   = TRIM( DiagnName ),  &
                             Array3D = Ptr3D,              &
                             RC      = RC )

          ! Free the pointer
          Ptr3D => NULL()
 
          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot update MET diagnostic ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
       ENDIF
    ENDDO

    ! Loop over 2D diagnostics within Met category
    DO N = 1, Num2D

       SELECT CASE ( N )
          CASE ( 1 )
             NameSuffix = 'TR-PAUSE_Level'  ! ND55, trcr 1
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 2 )
             NameSuffix = 'TR-PAUSE_km'     ! ND55, trcr 2
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 3 )
             NameSuffix = 'TR-PAUSE_mb'     ! ND55, trcr 3
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 4 )
             NameSuffix = 'LFLASH_TOTAL'    ! ND56, trcr 1
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 5 )
             NameSuffix = 'LFLASH_INTRA'    ! ND56, trcr 2
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 6 )
             NameSuffix = 'LFLASH_CLDGND'   ! ND56, trcr 3
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 7 )
             NameSuffix = 'HFLUX'           ! ND67, trcr 1
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 8 )
             NameSuffix = 'RADSWG'          ! ND67, trcr 2
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 9 )
             NameSuffix = 'PREACC'          ! ND67, trcr 3
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 10 )
             NameSuffix = 'PRECON'          ! ND67, trcr 4
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 11 )
             NameSuffix = 'TS'              ! ND67, trcr 5
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 12 )
             NameSuffix = 'RADSWT'          ! ND67, trcr 6
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 13 )
             NameSuffix = 'USTAR'           ! ND67, trcr 7
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 14 )
             NameSuffix = 'Z0'              ! ND67, trcr 8
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 15 )
             NameSuffix = 'PBL'             ! ND67, trcr 9
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 16 )
             NameSuffix = 'CLDFRC'          ! ND67, trcr 10
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 17 )
             NameSuffix = 'U10M'            ! ND67, trcr 11
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 18 )
             NameSuffix = 'V10M'            ! ND67, trcr 12
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 19 )
             NameSuffix = 'PS-PBL'          ! ND67, trcr 13
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 20 )
             NameSuffix = 'ALBD'            ! ND67, trcr 14
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 21 )
             NameSuffix = 'PHIS'            ! ND67, trcr 15
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 22 )
             NameSuffix = 'CLTOP'           ! ND67, trcr 16
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 23 )
             NameSuffix = 'TROPP'           ! ND67, trcr 17
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 24 )
             NameSuffix = 'SLP'             ! ND67, trcr 18
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 25 )
             NameSuffix = 'TSKIN'           ! ND67, trcr 19
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 26 )
             NameSuffix = 'PARDF'           ! ND67, trcr 20
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 27 )
             NameSuffix = 'PARDR'           ! ND67, trcr 21
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE ( 28 )
             NameSuffix = 'GWETTOP'         ! ND67, trcr 22
             Ptr2D => State_Met%TROPP  ! placeholder
          CASE DEFAULT
             IF ( N < Num2D ) THEN
                MSG = 'Num2D is less than number of named 2D MET diagnostics'
             ELSE
                MSG = 'Undefined 2D diagnostic in MET diagnostic group'
             ENDIF
             CALL ERROR_STOP( MSG, LOC ) 
       END SELECT

       ! Diagnostic info
       DiagnName = 'MET_' // NameSuffix

       ! Update diagnostics
       IF ( ASSOCIATED(Ptr2D) ) THEN 
          CALL Diagn_Update( am_I_Root,                   &
                             cName   = TRIM( DiagnName ), &
                             Array2D = Ptr2D,             &
                             RC      = RC )

          ! Free the pointer
          Ptr2D => NULL()
 
          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot update 2D MET diagnostic ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
       ENDIF
    ENDDO

  END SUBROUTINE DiagnUpdate_Met
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: diagnupdate_transport_flux
!
! !DESCRIPTION: Subroutine DIAGNUPDATE\_TRANSPORT\_FLUX updates the transport 
!  mass flux diagnostics (zonal ND24, meridional ND25, and vertical ND26) 
!  that are written to netCDF. This routine is called within various
!  TPCORE routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnUpdate_Transport_Flux( am_I_Root, FLUX_EW, FLUX_NS,  &
                                         FLUX_VERT, Input_Opt, RC )

!
! !USES:
!
   USE GIGC_Input_Opt_Mod, ONLY : OptInput    
   USE HCO_Diagn_Mod,      ONLY : Diagn_Update
!
! !INPUT PARAMETERS:
!
   LOGICAL,            INTENT(IN)  :: am_I_Root        ! Are we on root CPU?
   REAL(fp),           INTENT(IN)  :: FLUX_EW(:,:,:,:) ! east-west mass flux
   REAL(fp),           INTENT(IN)  :: FLUX_NS(:,:,:,:) ! north-south flux
   REAL(fp),           INTENT(IN)  :: FLUX_VERT(:,:,:,:) ! up-down flux
   TYPE(OptInput),     INTENT(IN)  :: Input_Opt   ! Input options object
!
! !INPUT/OUTPUT PARAMETERS:
!
   INTEGER,           INTENT(INOUT) :: RC          ! Success or failure
!
! !REVISION HISTORY:
!  12 Feb 2015 - E. Lundgren   - Initial version
!  19 Jan 2016 - E. Lundgren   - Updated to include all transport flux diags
! 
! !REMARKS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   INTEGER             :: M, NumLevels
   CHARACTER(LEN=60)   :: DiagnPrefix
   CHARACTER(LEN=255)  :: MSG
   CHARACTER(LEN=255)  :: LOC = 'DiagnUpdate_Transport_Flux' &
                                 // ' (Transport_mod.F)'

   ! Array to hold input mass flux array with reverse-order levels
   ! and only levels that have data (ie. less than LLPAR if input.geos
   ! levels for the diagnostic are less than LLPAR).
   REAL(fp)   :: FlxArray_ReverseLevels(IIPAR,JJPAR,LLPAR,NNPAR)

   !=================================================================
   ! DIAGNUPDATE_TRANSPORT_FULX begins here!
   !=================================================================
   
    ! Assume successful return
    RC = GIGC_SUCCESS

    ! Loop over diagnostics
    DO M = 1, 3

       ! Set diagnostic name suffix and reverse-order the levels such that 
       ! index 1 corresponds to surface. 
       SELECT CASE ( M )
          CASE ( 1 )
             DiagnPrefix = 'TRANSPORT_FLX_EW_'
             NumLevels = Input_Opt%LD24
             FlxArray_ReverseLevels(:,:,1:NumLevels,:) = &
                   FLUX_EW(:,:,LLPAR:LLPAR-NumLevels+1:-1,:)
          CASE ( 2 )
             DiagnPrefix = 'TRANSPORT_FLX_NS_'
             NumLevels = Input_Opt%LD24
             FlxArray_ReverseLevels(:,:,1:NumLevels,:) = &
                   FLUX_NS(:,:,LLPAR:LLPAR-NumLevels+1:-1,:)
          CASE ( 3 )
             DiagnPrefix = 'TRANSPORT_FLX_VERT_'
             NumLevels = Input_Opt%LD24
             FlxArray_ReverseLevels(:,:,1:NumLevels,:) = &
                   FLUX_VERT(:,:,LLPAR:LLPAR-NumLevels+1:-1,:)
       END SELECT

       ! Update the diagnostic container using routine in diagnostics module
       ! For now, pass diagnostic number 24 which can be temporarily used to
       ! to indicate which tracers to include for entire transport flux diag 
       ! group
       CALL DiagnUpdate_NTracers_3D( am_I_Root, DiagnPrefix, 24,  &
                                     FlxArray_ReverseLevels, Input_Opt, &
                                     RC )

    ENDDO

   END SUBROUTINE DiagnUpdate_Transport_Flux
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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE HCO_Diagn_Mod,      ONLY : Diagn_TotalGet
    USE TIME_MOD,           ONLY : GET_YEAR, GET_MONTH 
    USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoID
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
    INTEGER       :: I, cId, HCRC
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
          cId = GetHcoID( TrcId=I )
          IF ( cId <= 0 ) CYCLE
          
          ! Define diagnostics ID
          cId = 10000 + cId

          ! Get the total [kg] 
          CALL Diagn_TotalGet( am_I_Root,                         & 
                             cId     = cId,                       &
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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
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

! Not yet working diagnostics code

!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: diagninit_cloud_properties (LL in progress)
!!
!! !DESCRIPTION: Subroutine DIAGNINIT\_CLOUD\_PROPERTIES initializes the cloud
!!  optical depths and cloud fractions diagnostic (aka ND21).
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE DiagnInit_Cloud_Properties( am_I_Root, Input_Opt, RC )
!!
!! !USES:
!!
!    USE Error_Mod,          ONLY : Error_Stop
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
!    INTEGER              :: Collection, N
!    CHARACTER(LEN=15)    :: OutOper
!    CHARACTER(LEN=60)    :: DiagnName
!    CHARACTER(LEN=255)   :: MSG
!    CHARACTER(LEN=255)   :: LOC = 'DiagnInit_Cloud_Properties (diagnostics_mod.F90)' 
!
!    !=======================================================================
!    ! DIAGNINIT_CLOUD_PROPERTIES begins here!
!    !=======================================================================
!      
!    ! Assume successful return
!    RC = GIGC_SUCCESS
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
!                             cId       = cId,               &
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
!                             RC        = RC )
!
!          IF ( RC /= HCO_SUCCESS ) THEN
!             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
!             CALL ERROR_STOP( MSG, LOC ) 
!          ENDIF
!       ENDIF
!    ENDDO
!
!  END SUBROUTINE DiagnInit_Cloud_Properties
!!EOC

!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: diagninit_photolysis_rates (LL in progress)
!!
!! !DESCRIPTION: Subroutine DIAGNINIT\_PHOTOLYSIS\_RATES initializes the photolysis
!!  rates (J-values) diagnostic (aka ND22).
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE DiagnInit_Photolysis_Rates( am_I_Root, Input_Opt, RC )
!!
!! !USES:
!!
!    USE Error_Mod,          ONLY : Error_Stop
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
!    INTEGER              :: Collection, N
!    CHARACTER(LEN=15)    :: OutOper
!    CHARACTER(LEN=60)    :: DiagnName
!    CHARACTER(LEN=255)   :: MSG
!    CHARACTER(LEN=255)   :: LOC = 'DiagnInit_Photolysis_Rates (diagnostics_mod.F90)' 
!
!    !=======================================================================
!    ! DIAGNINIT_PHOTOLYSIS begins here!
!    !=======================================================================
!      
!    ! Assume successful return
!    RC = GIGC_SUCCESS
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
!                          cId       = cId,               &
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
!  END SUBROUTINE DiagnInit_Photolysis_Rates
!!EOC

!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: diagninit_hg_source (LL in progress)
!!
!! !DESCRIPTION: Subroutine DIAGNINIT\_HG\_SOURCE initializes the mercury
!!  emissions, production, and loss diagnostic (aka ND03).
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE DiagnInit_Hg_Source( am_I_Root, Input_Opt, RC )
!!
!! !USES:
!!
!    USE Error_Mod,          ONLY : Error_Stop
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
!    INTEGER              :: Collection, N
!!    INTEGER, PARAMETER  :: NUMTRC = 3
!!    INTEGER             :: TRCN ( NUMTRC )
!    CHARACTER(LEN=15)    :: OutOper
!    CHARACTER(LEN=60)    :: DiagnName
!    CHARACTER(LEN=255)   :: MSG
!    CHARACTER(LEN=255)   :: LOC = 'DIAGNINIT_HG_SOURCE (diagnostics_mod.F90)' 
!
!    !=======================================================================
!    ! DIAGNINIT_HG_SOURCE begins here!
!    !=======================================================================
!      
!    ! Assume successful return
!    RC = GIGC_SUCCESS
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
!                             cId       = cId,               &
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
!                             RC        = RC )
!
!          IF ( RC /= HCO_SUCCESS ) THEN
!             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
!             CALL ERROR_STOP( MSG, LOC ) 
!          ENDIF
!       ENDIF
!    ENDDO
!
!  END SUBROUTINE DiagnInit_Hg_Source
!!EOC

!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: diagninit_sulfate_prodloss (LL in progress)
!!
!! !DESCRIPTION: Subroutine DIAGNINIT\_SULFATE\_PRODLOSS initializes the sulfate
!!  chemistry production and loss diagnostic (aka ND05).
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE DiagnInit_Sulfate_ProdLoss( am_I_Root, Input_Opt, RC )
!!
!! !USES:
!!
!    USE GIGC_Input_Opt_Mod, ONLY : OptInput
!    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
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
!    INTEGER              :: Collection, N
!    CHARACTER(LEN=15)    :: OutOper
!    CHARACTER(LEN=60)    :: DiagnName
!    CHARACTER(LEN=255)   :: MSG
!    CHARACTER(LEN=255)   :: LOC = 'DiagnInit_Sulfate_ProdLoss (diagnostics_mod.F90)' 
!
!    !=======================================================================
!    ! DIAGNINIT_SULFATE_PRODLOSS begins here!
!    !=======================================================================
!      
!    ! Assume successful return
!    RC = GIGC_SUCCESS
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
!                             cId       = cId,               &
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
!                             cId       = cId,               &
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
!                             cId       = cId,               &
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
!                             RC        = RC )
!
!          IF ( RC /= HCO_SUCCESS ) THEN
!             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
!             CALL ERROR_STOP( MSG, LOC ) 
!
!       ENDIF
!    ENDDO
!!
!  END SUBROUTINE DiagnInit_Sulfate_ProdLoss
!!EOC

!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: diagninit_Carbon_Sources (LL in progress)
!!
!! !DESCRIPTION: Subroutine DIAGNINIT\_CARBON\_SOURCES initializes the sources of
!!  black carbon and organic carbon diagnostic (aka ND07).
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE DiagnInit_Carbon_Sources( am_I_Root, Input_Opt, RC )
!!
!! !USES:
!!
!    USE Error_Mod,          ONLY : Error_Stop
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
!    INTEGER              :: Collection, N
!    CHARACTER(LEN=15)    :: OutOper
!    CHARACTER(LEN=60)    :: DiagnName
!    CHARACTER(LEN=255)   :: MSG
!    CHARACTER(LEN=255)   :: LOC = 'DiagnInit_Carbon_Sources (diagnostics_mod.F90)' 
!
!    !=======================================================================
!    ! DIAGNINIT_CARBON_SOURCES begins here!
!    !=======================================================================
!      
!    ! Assume successful return
!    RC = GIGC_SUCCESS
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
!                             cId       = cId,               &
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
!                             RC        = RC )
!
!          IF ( RC /= HCO_SUCCESS ) THEN
!             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
!             CALL ERROR_STOP( MSG, LOC ) 
!          ENDIF
!       ENDIF
!    ENDDO
!
!  END SUBROUTINE DiagnInit_Carbon_Sources
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: diagninit_hrlymax_surfconc
!!
!! !DESCRIPTION: Subroutine DIAGNINIT\_HRLYMAX\_SURFCONC initializes the 
!!  diagnostics for hourly maximum tracer mixing ratio at the surface (aka ND71).
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE DiagnInit_HrlyMax_SurfConc( am_I_Root, Input_Opt, RC )
!!
!! !USES:
!!
!    USE GIGC_Input_Opt_Mod, ONLY : OptInput
!    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
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
!!  30 Jan 2015 - M. Yannetti - Initial version
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    INTEGER              :: Collection, N, M
!    CHARACTER(LEN=15)    :: OutOper 
!    CHARACTER(LEN=60)    :: DiagnName
!    CHARACTER(LEN=255)   :: MSG
!    CHARACTER(LEN=255)   :: LOC = 'DIAGNINIT_HRLYMAX_SURFCONC (diagnostics_mod.F90)'
!
!    !=======================================================================
!    ! DIAGNINIT_HRLYMAX_SURFCONC begins here!
!    !=======================================================================
!
!    ! Assume successful return
!    RC = GIGC_SUCCESS
!
!    ! Get diagnostic parameters from the Input_Opt object
!    Collection = Input_Opt%DIAG_COLLECTION
!    OutOper    = Input_Opt%ND71_OUTPUT_TYPE
!
!    ! TODO Check if certain tracer(s) listed for ND71 in input.geos 
!
!    !----------------------------------------------------------------
!    ! Create containers for Tracer Mixing Ratio
!    !----------------------------------------------------------------
!
!    ! Diagnostic name
!    DiagnName = 'IJ-MAX'
!    DO M = 1, Input_Opt%TMAX(71)
!       N = Input_Opt%TINDEX(71,M)
!       IF ( N > Input_Opt%N_TRACERS ) CYCLE
!
!!       ! line from diag3
!!       ARRAY(:,:,1) = AD71(:,:,N)/AD71_COUNT
!
!       ! Create container
!       CALL Diagn_Create( am_I_Root,                  &
!                       Col       = Collection,        &
!                       cId       = cId,               &
!                       cName     = TRIM( DiagnName ), &
!                       AutoFill  = 0,                 &
!                       ExtNr     = -1,                &
!                       Cat       = -1,                &
!                       Hier      = -1,                &
!                       HcoID     = -1,                &
!                       SpaceDim  =  3,                &
!                       LevIDx    = -1,                &
!                       OutUnit   = 'v/v' ,            &
!                       OutOper   = TRIM( OutOper   ), &
!                       RC        = RC )
!
!       IF ( RC /= HCO_SUCCESS ) THEN
!          MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
!          CALL ERROR_STOP( MSG, LOC )
!       ENDIF
!    ENDDO
!
!  END SUBROUTINE DiagnInit_HrlyMax_SurfConc
!!EOC


END MODULE Diagnostics_Mod
