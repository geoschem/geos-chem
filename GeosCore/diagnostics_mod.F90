!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  USE CMN_SIZE_Mod
  USE ErrCode_Mod
  USE Error_Mod,          ONLY : Error_Stop
  USE HCO_Error_Mod
  USE HCO_INTERFACE_MOD,  ONLY : HcoState
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Diagnostics_Write

  ! All of the following is currently only active for NETCDF=y
#if defined( NC_DIAG )
  PUBLIC  :: Diagnostics_Init
  PUBLIC  :: Diagnostics_Final
  PUBLIC  :: DiagnUpdate_Met            ! ND31, 55, 57, 66, 67, and 68
  PUBLIC  :: DiagnUpdate_Transport_Flux ! ND24, ND25, ND26
#if defined( DIAG_DEVEL )
  PUBLIC  :: CalcDobsonColumn           ! added by ck - no existing bpch diag
#endif
!
! !PRIVATE MEMBER FUNCTIONS:
!
! Benchmark diagnostic groups
  PRIVATE :: DiagnInit_WetScav          ! ND38 and 39 (add 37)
  PRIVATE :: DiagnInit_DryDep           ! ND44
  PRIVATE :: DiagnInit_Species_Conc     ! ND45
  PRIVATE :: DiagnInit_Met              ! ND31, 55, 57, 66, 67, and 68

  !PRIVATE :: DiagnInit_Sulfate_ProdLoss ! ND05 (implementation not complete)
  !PRIVATE :: DiagnInit_Carbon_Sources   ! ND07 (implementation not complete)
  !PRIVATE :: DiagnInit_Cloud_Properties ! ND21 (implementation not complete)
  !PRIVATE :: DiagnInit_Photolysis_Rates ! ND22 (implementation not complete)

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
  PRIVATE :: DiagnInit_Transport_Flux   ! ND24, 25, and 26
  PRIVATE :: DiagnInit_BL_Frac          ! ND12 init
  PRIVATE :: DiagnInit_CldConv_Flx      ! ND14 init
  PRIVATE :: DiagnInit_BlMix_Flx        ! ND15 init
  PRIVATE :: DiagnInit_Precip_Frac      ! ND16 init
  PRIVATE :: DiagnInit_Rain_Frac        ! ND17 init
  PRIVATE :: DiagnInit_Wash_Frac        ! ND18 init
  PRIVATE :: DiagnInit_Landmap          ! ND30 init
#if defined( DIAG_DEVEL )
  PRIVATE :: DiagnInit_Species_Emis     ! added by ck
#endif
  !PRIVATE :: DiagnInit_HrlyMax_SurfConc ! ND71 init (impementation commented
                                         ! out pending decision on if to keep)
  !   to do: ND37 (updraft scav fraction)
  !          ND62 (instantaneous column maps)

#if defined( DIAG_DEVEL )
  PRIVATE :: DiagnInit_Dobson          ! added by ck
#endif

#if defined( DIAG_DEVEL )
  ! Diagnostics for testing FlexChem (may be removed later)
  PRIVATE :: DiagnInit_KPP_Rates
  PRIVATE :: DiagnInit_KPP_Spec
#endif
!
! !DEFINED PARAMETERS:
!
  ! Prefix of diagnostics file. This file will hold all diagnostics that 
  ! are written out at the end of a simulation (either because their
  ! output frequency is set to 'End' or the run finishes and these
  ! diagnostics haven't reached the end of their output interval yet).
  CHARACTER(LEN=31), PARAMETER :: DGN = 'GEOSCHEM_Diagnostics_Hrly'

  ! Toggle to enable species diagnostics. This will write out chemical species 
  ! concentrations (in addition to the advected species). Not recommended unless
  ! you have a good reason (ckeller, 8/11/2015).
  LOGICAL, PARAMETER, PUBLIC   :: DiagnSpec = .FALSE.

  ! Also write out species tendencies
  LOGICAL, PARAMETER, PUBLIC   :: DiagnSpecTend = .FALSE.

  ! Initialize GEOS-Chem diagnostics container id (ewl, 1/20/16)
  ! If we add containers, we will need to modify this to have a different
  ! container id variable per collection.
  INTEGER, SAVE, PRIVATE       :: cId = 0            

#endif ! netcdf 1
!
! !REVISION HISTORY:
!  09 Jan 2015 - C. Keller   - Initial version. 
!  14 Jan 2016 - E. Lundgren - Add several GEOS-Chem diagnostics
!  29 Jan 2016 - E. Lundgren - Update diagnostics for recent HEMCO updates
!  20 Jul 2016 - R. Yantosca - Replace #ifdef DEVEL with #ifdef DIAG_DEVEL
!  29 Nov 2016 - R. Yantosca - grid_mod.F90 is now gc_grid_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Species ID flags
  INTEGER :: id_Rn, id_Pb, id_Be7

CONTAINS
#if defined( NC_DIAG )
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
    USE GC_GRID_MOD,        ONLY : AREA_M2
    USE HCO_DIAGN_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE TIME_MOD,           ONLY : GET_TS_CHEM
#if defined( USE_TEND )
    USE TENDENCIES_MOD,     ONLY : TEND_INIT
#endif
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   ) :: am_I_Root  ! Are we on the root CPU?
    TYPE(MetState),   INTENT(IN   ) :: State_Met  ! Met state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm  ! Chemistry state 
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
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  22 Jun 2016 - R. Yantosca - Now use Ind_() to define id_Rn, id_Pb, id_Be7
!  22 Jun 2016 - R. Yantosca - Remove reference to Species type
!  01 Jul 2016 - R. Yantosca - Pass State_Chm to several other routines
!  19 Jul 2016 - R. Yantosca - Now bracket tendency calls with #ifdef USE_TEND
!  19 Jul 2016 - R. Yantosca - Add missing nAdvect variable
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N
    INTEGER            :: CollectionID
    INTEGER            :: DeltaYMD,    DeltaHMS 
    REAL(sp)           :: TS

    ! Pointers
    REAL(fp), POINTER  :: AM2(:,:)

    ! Strings
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'Diagnostics_Init (diagnostics_mod.F90)'

    !=======================================================================
    ! Diagnostics_Init begins here 
    !=======================================================================

    ! Assume successful return
    RC      =  GC_SUCCESS

    ! Define collection variables
    AM2     => AREA_M2(:,:,1)
    TS      =  GET_TS_CHEM() * 60.0_sp

    !-----------------------------------------------------------------------
    ! Store local diagnostic species ID flags
    !-----------------------------------------------------------------------

    ! Get species ID flags
    id_Rn   =  Ind_('Rn' )
    id_Pb   =  Ind_('Pb' )
    id_Be7  =  Ind_('Be7')

    ! Double check that HcoState and also the diagnostics bundle are
    ! available
    IF ( .NOT. ASSOCIATED(HcoState) ) THEN
       MSG = 'Cannot initialize diagnostics - HEMCO state '     // &
             'is not yet initialized. The HEMCO initialization '// &
             'routines must be called first.'
       CALL ERROR_STOP ( MSG, LOC )
       RETURN
    ENDIF
    IF ( .NOT. ASSOCIATED(HcoState%Diagn) ) THEN
       MSG = 'Cannot initialize diagnostics - HEMCO diagnostics '       // &
             'bundle is not yet initialized. The HEMCO initialization ' // &
             'routines must be called first.'
       CALL ERROR_STOP ( MSG, LOC )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Create diagnostics collection for GEOS-Chem.  This will keep the
    ! GEOS-Chem diagostics separate from the HEMCO diagnostics.
    !-----------------------------------------------------------------------

    ! Define output write frequency. In ESMF environment, make sure 
    ! diagnostics is always passed to MAPL history!
    CALL DiagnCollection_GetDefaultDelta ( am_I_Root, HcoState, &
       deltaYMD, deltaHMS, RC )
    IF ( RC /= HCO_SUCCESS ) RETURN

    ! Temporarily manually set to hourly. Eventually this will be set
    ! per collection in input.geos or a new diagnostics input file
    !  (ewl, 1/13/16)
    CALL DiagnCollection_Create( am_I_Root, HcoState%Diagn,      &
                                 NX           = IIPAR,           &
                                 NY           = JJPAR,           &
                                 NZ           = LLPAR,           &
                                 TS           = TS,              &
                                 AM2          = AM2,             &
                                 PREFIX       = DGN,             &
                                 deltaYMD     = deltaYMD,        &
                                 deltaHMS     = deltaHMS,        &
                                 !deltaYMD     = 00000000,        &
                                 !deltaHMS     = 010000,          &
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
!      IF ( RC /= GC_SUCCESS ) THEN
!         CALL ERROR_STOP( 'Error in DiagnInit_Sulfate_ProdLoss', LOC ) 
!      ENDIF
!   ENDIF
!
!   ! Carbon sources diagnostic (ND07)
!   IF ( Input_Opt%ND07 > 0 ) THEN
!      CALL DiagnInit_Carbon_Sources( am_I_Root, Input_Opt, RC )
!      IF ( RC /= GC_SUCCESS ) THEN
!         CALL ERROR_STOP( 'Error in DiagnInit_Carbon_Sources', LOC ) 
!      ENDIF
!   ENDIF
!
!   ! Cloud properties (ND21)
!   IF ( Input_Opt%ND21 > 0 ) THEN
!      CALL DiagnInit_Cloud_Properties( am_I_Root, Input_Opt, RC )
!      IF ( RC /= GC_SUCCESS ) THEN
!         CALL ERROR_STOP( 'Error in DiagnInit_Cloud_Properties', LOC ) 
!      ENDIF
!   ENDIF
!
!   ! Photolysis rates diagnostic (ND22)
!   IF ( Input_Opt%ND22 > 0 ) THEN
!      CALL DiagnInit_Photolysis_Rates( am_I_Root, Input_Opt, RC )
!      IF ( RC /= GC_SUCCESS ) THEN
!         CALL ERROR_STOP( 'Error in DiagnInit_Photolysis_Rates', LOC ) 
!      ENDIF
!   ENDIF

    ! Transport fluxes diagnostic (ND24, ND25, and ND26)
    ! For now, use ND24 as indicator of whether to turn on diag group
    IF ( Input_Opt%ND24 > 0 ) THEN
       CALL DiagnInit_Transport_Flux( am_I_Root, Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Transport_Flux', LOC ) 
       ENDIF
    ENDIF

    ! Wet scavenging diagnostic group (ND38 and ND39) (add others?)
    ! Assume both ND38 and ND39 are on if ND38 is on
    IF ( Input_Opt%ND38 > 0 ) THEN
       CALL DiagnInit_WetScav( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_WetScav', LOC ) 
       ENDIF
    ENDIF

    ! Drydep diagnostic (ND44)
    IF ( Input_Opt%ND44 > 0 ) THEN
       CALL DiagnInit_DryDep( am_I_Root, Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_DryDep', LOC ) 
       ENDIF
    ENDIF

    ! Species concentration diagnostics (ND45)
    IF ( Input_Opt%ND45 > 0 ) THEN
       CALL DiagnInit_Species_Conc( am_I_Root, Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Species_Conc', LOC ) 
       ENDIF
    ENDIF

    ! Meteorology state diagnostic (ND31, 55, 57, 66, 67, and 68)
    ! For now, use ND68 as indicator for entire diagnostic group
    IF ( Input_Opt%ND68 > 0 ) THEN
       CALL DiagnInit_Met( am_I_Root, Input_Opt, State_Met, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Met', LOC )
       ENDIF
    ENDIF

    ! Pb emissions diagnostic (ND01)
    IF ( Input_Opt%ND01 > 0 ) THEN
       CALL DiagnInit_Pb_Emiss( am_I_Root, Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Pb_Emiss', LOC ) 
       ENDIF
    ENDIF

    ! Rn/Pb/Be decay diagnostic (ND02)
    IF ( Input_Opt%ND02 > 0 ) THEN
       CALL DiagnInit_Rn_Decay( am_I_Root, Input_Opt, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Rn_Decay', LOC ) 
       ENDIF
    ENDIF

    ! Boundary layer fraction diagnostic (ND12)
    IF ( Input_Opt%ND12 > 0 ) THEN
       CALL DiagnInit_BL_Frac( am_I_Root, Input_Opt, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_BL_Frac', LOC ) 
       ENDIF
    ENDIF

    ! Cloud convection mass flux diagnostic (ND14)
    IF ( Input_Opt%ND14 > 0 ) THEN
       CALL DiagnInit_CldConv_Flx( am_I_Root, Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_CldConv_Flx', LOC ) 
       ENDIF
    ENDIF

    ! Boundary-layer mixing mass flux diagnostic (ND15)
    IF ( Input_Opt%ND15 > 0 ) THEN
       CALL DiagnInit_BLMix_Flx( am_I_Root, Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_BLMix_Flx', LOC ) 
       ENDIF
    ENDIF

    ! Areal fraction of precip diagnostic (ND16)
    IF ( Input_Opt%ND16 > 0 ) THEN
       CALL DiagnInit_Precip_Frac( am_I_Root, Input_Opt, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Precip_Frac', LOC ) 
       ENDIF
    ENDIF

    ! Rainout fraction diagnostic (ND17)
    IF ( Input_Opt%ND17 > 0 ) THEN
       CALL DiagnInit_Rain_Frac( am_I_Root, Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Rain_Frac', LOC ) 
       ENDIF
    ENDIF

    ! Washout fraction diagnostic (ND18)
    IF ( Input_Opt%ND18 > 0 ) THEN
       CALL DiagnInit_Wash_Frac( am_I_Root, Input_Opt, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Wash_Frac', LOC ) 
       ENDIF
    ENDIF

    ! CH4 loss diagnostic (ND19)
    IF ( Input_Opt%ND19 > 0 ) THEN
       CALL DiagnInit_CH4_Loss(am_I_Root, Input_Opt, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_CH4_Loss', LOC ) 
       ENDIF
    ENDIF


    ! Land map diagnostic (ND30)
    IF ( Input_Opt%ND30 > 0 ) THEN
       CALL DiagnInit_LandMap( am_I_Root, Input_Opt, RC )
       IF ( RC /= GC_SUCCESS ) THEN
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
!       CALL DiagnInit_UCX( am_I_Root, Input_Opt, State_Chm, RC )
!       IF ( RC /= GC_SUCCESS ) THEN
!          CALL ERROR_STOP( 'Error in DiagnInit_UCX', LOC ) 
!       ENDIF
!    ENDIF

    !! Houly-maximum species mixing ratio (IJ-MAX) at surface (ND71)
    !IF ( Input_Opt%ND71 > 0 ) THEN
    !   CALL DiagnInit_HrlyMax_SurfConc( am_I_Root, Input_Opt, RC )
    !   IF ( RC /= GC_SUCCESS ) THEN
    !      CALL ERROR_STOP( 'Error in DiagnInit_HrlyMax_SurfConc', LOC )
    !   ENDIF
    !ENDIF

#if defined( DIAG_DEVEL )
    CALL DiagnInit_Dobson( am_I_Root, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DiagnInit_Dobson', LOC ) 
    ENDIF
#endif

#if defined( USE_TEND )
    ! Initialize tendencies
    CALL Tend_Init( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in Tendencies_Init', LOC ) 
    ENDIF
#endif

#if defined( DIAG_DEVEL )
    ! Species emission diagnostics (NEW) (added by Christoph)
    ! NOTE: Currently this diagnostic must be initialized last since since
    ! container ids start at 10000 for routine TotalsToLogFile (ewl, 1/20/16)
    IF ( .FALSE. ) THEN
       CALL DiagnInit_Species_Emis( am_I_Root, Input_Opt, State_Met, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP( 'Error in DiagnInit_Species_Emis', LOC ) 
       ENDIF
    ENDIF
#endif

#if defined( DIAG_DEVEL )
    ! KPP diagnostics
    CALL DiagnInit_KPP_Rates( am_I_Root, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DiagnInit_KPP_Rates', LOC ) 
    ENDIF

    CALL DiagnInit_KPP_Spec( am_I_Root, Input_Opt, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL ERROR_STOP( 'Error in DiagnInit_KPP_Spec', LOC ) 
    ENDIF
#endif

    ! Leave with success
    RC = GC_SUCCESS

  END SUBROUTINE Diagnostics_Init
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  SUBROUTINE Diagnostics_Final( am_I_Root, RC ) 
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
#if defined( USE_TEND )
    USE TENDENCIES_MOD,     ONLY : TEND_CLEANUP
#endif
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN ) :: am_I_Root  ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  09 Jan 2015 - C. Keller   - Initial version 
!  15 Jan 2015 - R. Yantosca - Now accept Input_Opt, am_I_Root, RC arguments
!  19 Jul 2016 - R. Yantosca - Now bracket tendency calls with #ifdef USE_TEND
!EOP
!------------------------------------------------------------------------------
!BOC

!    ! Finalize diagnostics
!    CALL DiagnCollection_Cleanup( COL = Input_Opt%DIAG_COLLECTION )

#if defined( USE_TEND )
    CALL TEND_CLEANUP()
#endif

    ! Return with success
    RC = GC_SUCCESS

  END SUBROUTINE Diagnostics_Final
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
  SUBROUTINE DiagnInit_Transport_Flux( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  19 Jan 2016 - E. Lundgren - Initial version
!  01 Jul 2016 - R. Yantosca - Use State_Chm%nAdvect to replace N_TRACERS
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: Collection, N, M, NA

    ! Strings
    CHARACTER(LEN=15)      :: OutOper
    CHARACTER(LEN=60)      :: DiagnName
    CHARACTER(LEN=30)      :: NamePrefix
    CHARACTER(LEN=255)     :: MSG
    CHARACTER(LEN=255)     :: LOC='DiagnInit_Transport_Flux (diagnostics_mod.F90)' 

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !=======================================================================
    ! DIAGNINIT_TRANSPORT_FLUX begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    ! For now, use ND24 for entire transport diagnostic group
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%TRANSPORT_OUTPUT_TYPE

    ! Loop over 3 types of transport fluxes
    DO M = 1, 3

       ! Loop over advected species
       DO NA = 1, State_Chm%nAdvect
         
          ! Species ID
          N = State_Chm%Map_Advect(NA)

          ! Entry in the species database
          SpcInfo => State_Chm%SpcData(N)%Info

          ! If this species number N is scheduled for output
          ! then define the diagnostic containers for mass flux
          ! NOTE: for now, use species indicated for ND24 only until 
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
             DiagnName = TRIM( NamePrefix ) // TRIM( SpcInfo%Name )
             cId = cId + 1
   
             ! Create container
             CALL Diagn_Create( am_I_Root,  HcoState,          &
                                Col       = Collection,        & 
                                cId       = cId,               &
                                cName     = TRIM( DiagnName ), &
                                SpaceDim  =  3,                &
                                OutUnit   = 'kg/s' ,           &
                                OutOper   = TRIM( OutOper   ), &
                                RC        = RC )
         
             IF ( RC /= HCO_SUCCESS ) THEN
                MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
                CALL ERROR_STOP( MSG, LOC ) 
             ENDIF  
          ENDIF

          ! Free pointer
          SpcInfo => NULL()
       ENDDO
    ENDDO

  END SUBROUTINE DiagnInit_Transport_Flux
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_wetscav
!
! !DESCRIPTION: Subroutine DIAGNINIT\_WETSCAV initializes the wet
!  scavenging loss diagnostics, including loss in moist convection (ND38)
!  and loss in aerosol wet deposition (ND39).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_WetScav( am_I_Root, Input_Opt, State_Met,  &
                                State_Chm, RC )
!
! !USES:
!
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species   
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Met state
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  21 Jan 2016 - E. Lundgren - Initial version
!  01 Jul 2016 - R. Yantosca - Use State_Chm%Map_WETDEP to get the wetdep ID's
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: N, NW
    INTEGER                :: Collection

    ! Strings
    CHARACTER(LEN=15)      :: OutOper
    CHARACTER(LEN=60)      :: DiagnName
    CHARACTER(LEN=255)     :: MSG
    CHARACTER(LEN=255)     :: LOC = 'DiagnInit_WetScav (diagnostics_mod.F90)' 

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !=======================================================================
    ! DIAGNINIT_WETSCAV begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    ! Use same output frequency and operations as for species concentrations.
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%WETSCAV_OUTPUT_TYPE ! Assume same for ND38 and ND39

    ! Loop over # of soluble species 
    DO NW = 1, State_Chm%nWetDep

       ! Get the species ID from the wetdep ID
       N = State_Chm%Map_WetDep(NW)

       ! Get info about the Nth species from the species database
       SpcInfo => State_Chm%SpcData(N)%Info
          
       ! Check if this species is turned on for this diagnostic
       ! NOTE: Use ND38 species for entire group
       IF ( ANY( Input_Opt%TINDEX(38,:) == N ) ) THEN

          !----------------------------------------------------------------
          ! Create container for convective loss (kg/s) (ND38) 
          !----------------------------------------------------------------

          ! Diagnostic container name and id
          DiagnName = 'WETSCAV_CONVLOSS_' // TRIM( SpcInfo%Name )
          cId = cId + 1

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             HcoState  = HcoState,          & 
                             Col       = Collection,        & 
                             cId       = cId,               &
                             cName     = TRIM( DiagnName ), &
                             SpaceDim  =  2,                & ! 2D for now!!
                             OutUnit   = 'kg/s',            &
                             OutOper   = TRIM( OutOper   ), &
                             RC        = RC                  )

          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF

          !----------------------------------------------------------------
          ! Create container for wetdep loss (kg/s) (ND39)
          !----------------------------------------------------------------

          ! Skip if this is not a wet-depositing species
          IF ( .not. SpcInfo%Is_WetDep ) CYCLE

          ! Diagnostic container name and id
          DiagnName = 'WETSCAV_DEPLOSS_' // TRIM( SpcInfo%Name )
          cId = cId + 1

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             HcoState  = HcoState,          & 
                             Col       = Collection,        &
                             cID       = cID,               &
                             cName     = TRIM( DiagnName ), &
                             SpaceDim  =  2,                &
                             OutUnit   = 'kg/s',            &
                             OutOper   = TRIM( OutOper   ), &
                             RC        = RC                  )

          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
       ENDIF

       ! Free pointer
       SpcInfo => NULL()
    ENDDO

  END SUBROUTINE DiagnInit_WetScav
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_drydep
!
! !DESCRIPTION: Subroutine DIAGNINIT\_DRYDEP initializes the dry deposition 
! diagnostics, including drydep velocity and drydep flux. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_DryDep( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMTERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  13 Jan 2015 - C. Keller   - Initial version 
!  15 Jan 2015 - R. Yantosca - Init drydep velocity & flux diagnostics
!  01 Jul 2016 - R. Yantosca - Use State_Chm%Map_Drydep to get the drydep ID's
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: Collection, D, N, M

    ! Strings
    CHARACTER(LEN=15)      :: OutOper
    CHARACTER(LEN=60)      :: Prefix, Units, DiagnName
    CHARACTER(LEN=255)     :: MSG
    CHARACTER(LEN=255)     :: LOC = 'DiagnInit_DryDep (diagnostics_mod.F)' 
    
    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !=======================================================================
    ! DIAGNINIT_DRYDEP begins here!
    !=======================================================================

    ! Assume successful return
    RC = GC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%DRYDEP_OUTPUT_TYPE

    ! Loop over 3 types of drydep diagnostics
    DO M = 1, 3

       SELECT CASE ( M )
          CASE ( 1 )
             Prefix = 'DRYDEP_VEL_'
             Units  = 'cm s-1'
          CASE ( 2 )
             Prefix = 'DRYDEP_FLX_CHEM_'
             Units  = 'molec cm-2 s-1'
          CASE ( 3 )
             Prefix = 'DRYDEP_FLX_MIX_'
             Units  = 'molec cm-2 s-2'
       END SELECT

       ! Loop over # of depositing species
       DO D = 1, State_Chm%nDryDep
            
          ! Corresponding GEOS-Chem species number
          N = State_Chm%Map_DryDep(D)
       
          ! Entry in the species database
          SpcInfo => State_Chm%SpcData(N)%Info

          ! If this species number N is scheduled for output in input.geos, 
          ! then define the diagnostic containers for drydep velocity & flux
          IF ( ANY( Input_Opt%TINDEX(44,:) == N ) ) THEN
       
             !----------------------------------------------------------------
             ! Create diagnostic container
             !----------------------------------------------------------------
       
             ! Diagnostic container name and id
             DiagnName = TRIM( Prefix ) // TRIM( SpcInfo%Name )
             cId = cId + 1
       
             ! Create container
             CALL Diagn_Create( am_I_Root,                     &
                                HcoState  = HcoState,          &
                                Col       = Collection,        &
                                cId       = cId,               & 
                                cName     = TRIM( DiagnName ), &
                                SpaceDim  =  2,                &
                                OutUnit   = TRIM( Units ),     &
                                OutOper   = TRIM( OutOper ),   &
                                OkIfExist = .TRUE.,            &
                                RC        = RC )
       
             IF ( RC /= HCO_SUCCESS ) THEN
                MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
                CALL ERROR_STOP( MSG, LOC ) 
             ENDIF
          ENDIF

          ! Free pointer
          SpcInfo => NULL()
       ENDDO
    ENDDO

  END SUBROUTINE DiagnInit_DryDep
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_species_conc
!
! !DESCRIPTION: Subroutine DIAGNINIT\_SPECIES\_CONC initializes the species
!  concentration diagnostic (aka ND45).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_Species_Conc( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN   ) :: State_Chm  ! Chemistry state 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  20 Jan 2015 - R. Yantosca - Initial version
!  13 Jan 2016 - E. Lundgren - Define diagnostic ID
!  06 Jun 2016 - M. Sulprizio- Replace NTSPEC with State_Chm%nSpecies and
!                              NAMEGAS with ThisSpc%Name from species database
!  01 Jul 2016 - R. Yantosca - Use module variable nAdvect to replace N_TRACERS
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: Collection, N, NA
    CHARACTER(LEN=15)  :: OutOper
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGNINIT_SPECIES_CONC (diagnostics_mod.F90)' 

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    !=======================================================================
    ! DIAGNINIT_SPECIES_CONC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GC_SUCCESS

    ! Initialize pointer
    SpcInfo => NULL()

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%SPECIES_CONC_OUTPUT_TYPE
      
    ! Loop over all species
    DO N = 1, State_Chm%nSpecies
         
       ! Get info about this species from the species database
       SpcInfo => State_Chm%SpcData(N)%Info

       ! Only save out chemical species if DiagnSpec option is on
       IF ( .not. DiagnSpec .and. .not. SpcInfo%Is_Advected ) THEN
          CYCLE
       ENDIF
       
       ! If this species number N is scheduled for output in input.geos, 
       ! then define the diagnostic containers for drydep velocity & flux
       IF ( ANY( Input_Opt%TINDEX(45,:) == N ) ) THEN
    
          !----------------------------------------------------------------
          ! Create containers for species concentrations
          !----------------------------------------------------------------
    
          ! Diagnostic container name and id
          DiagnName = 'SPECIES_CONC_' // TRIM( SpcInfo%Name )
          cId = cId + 1
    
          ! Create container
          CALL Diagn_Create( am_I_Root,  HcoState,          &
                             Col       = Collection,        & 
                             cId       = cId,               &
                             cName     = TRIM( DiagnName ), &
                             SpaceDim  =  3,                &
                             OutUnit   = 'v/v',             &
                             OutOper   = TRIM( OutOper   ), &
                             RC        = RC )
    
          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
       ENDIF

       ! Free pointer
       SpcInfo => NULL()
    ENDDO
   
   END SUBROUTINE DiagnInit_Species_Conc
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_met
!
! !DESCRIPTION: Subroutine DIAGNINIT\_MET initializes the meteorology state
!  diagnostics (aka ND31, ND55, ND57, ND66, ND67, and ND68).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_Met( am_I_Root, Input_Opt, State_Met, RC )
!
! !USES:
!
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants,      ONLY : XNUMOLAIR
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN   ) :: State_Met  ! Met state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
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
    ! Scalars
    INTEGER            :: Collection, N, Num2D, Num3D
    INTEGER            :: SpaceDim 

    ! Strings
    CHARACTER(LEN=15)  :: OutOper, Units
    CHARACTER(LEN=30)  :: NameSuffix, DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DiagnInit_Met (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_MET begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%MET_OUTPUT_TYPE

    ! Set number of 3D and 2D MET diagnostics
    Num3d = 12
    Num2d = 26

    !----------------------------------------------------------------
    ! Create containers
    !----------------------------------------------------------------
  
    ! Loop over 3D diagnostics within Met category
    DO N = 1, Num3D

       SELECT CASE ( N )
          CASE ( 1 )
             NameSuffix = 'THETA'        ! ND57, trcr 1: potential temperature
             Units      = 'K'
          CASE ( 2 )
             NameSuffix = 'UWND'         ! ND66, trcr 1: Zonal wind
             Units      = 'm/s'
          CASE ( 3 )
             NameSuffix = 'VWND'         ! ND66, trcr 2: Meridional wind
             Units      = 'm/s'
          CASE ( 4 )
             NameSuffix = 'TMPU'         ! ND66, trcr 3: Temperature
             Units      = 'K'
          CASE ( 5 )
             NameSuffix = 'SPHU'         ! ND66, trcr 4: Specific humidity 
             Units      = 'g H2O/kg air'
          CASE ( 6 )
             NameSuffix = 'CLDMAS'       ! ND66, trcr 5: Convective mass flux
             Units      = 'kg/m2/s'      
          CASE ( 7 )
             NameSuffix = 'DTRAIN'       ! ND66, trcr 6: Detrainment flux
             Units      = 'kg/m2/s'
          CASE ( 8 )
             NameSuffix = 'BXHEIGHT'     ! ND68, trcr 1
             Units      = 'm'
          CASE ( 9 )
             NameSuffix = 'DRYAIRMASS'   ! ND68, trcr 2
             Units      = 'kg'
          CASE ( 10 )
             NameSuffix = 'AVGW'         ! ND68, trcr 3
             Units      = 'v/v'
          CASE ( 11 )
             NameSuffix = 'NAIR'         ! ND68, trcr 4
             Units      = 'molec/m3'
          CASE ( 12 )
             NameSuffix = 'PEDGE'        ! ND31, trcr 1
             Units      = 'hPa'
          CASE DEFAULT
             IF ( N < Num3D ) THEN
                MSG = 'Num3D is less than number of named 3D MET diagnostics'
             ELSE
                MSG = 'Undefined 3D diagnostic case in MET diagnostic group'
             ENDIF
             CALL ERROR_STOP( MSG, LOC ) 
       END SELECT

       ! Diagnostic container info
       DiagnName = 'MET_3D_' // NameSuffix
       cId       = cId + 1
       SpaceDim  = 3

       ! Create container
       CALL Diagn_Create( am_I_Root,  HcoState,          &
                          Col       = Collection,        & 
                          cId       = cId,               &
                          cName     = TRIM( DiagnName ), &
                          SpaceDim  =  SpaceDim,         &
                          OutUnit   = TRIM( Units ),     &
                          OutOper   = TRIM( OutOper ),   &
                          RC        = RC )

       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot create 3D MET diagnostic ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF
    ENDDO

    ! Loop over 2D diagnostics within Met category
    DO N = 1, Num2D

       SELECT CASE ( N )
          CASE ( 1 )
             NameSuffix = 'TRPAUSE_LVL'  ! ND55, trcr 1
             Units      = 'level'
          CASE ( 2 )
             NameSuffix = 'TRPAUSE_HGHT' ! ND55, trcr 2
             Units      = 'km'
          CASE ( 3 )
             NameSuffix = 'TRPAUSE_PRESS' ! ND55, trcr 3: Tropopause Pressure
             Units      = 'mb'
          CASE ( 4 )
             NameSuffix = 'HFLUX'        ! ND67, trcr 1: Sens heat flux
             Units      = 'W/m2'
          CASE ( 5 )
             NameSuffix = 'RADSWG'       ! ND67, trcr 2: SW rad @ sfc
             Units      = 'W/m2'
          CASE ( 6 )
             NameSuffix = 'PREACC'       ! ND67, trcr 3: Tot prec [kg/m2/s]???
             Units      = 'mm/day'
          CASE ( 7 )
             NameSuffix = 'PRECON'       ! ND67, trcr 4: Sfc conv prec[kg/m2/s]?
             Units      = 'mm/day'
          CASE ( 8 )
             NameSuffix = 'TS'           ! ND67, trcr 5: T @ 2m height
             Units      = 'K'
          CASE ( 9 )
             NameSuffix = 'RADSWT'       ! ND67, trcr 6
             Units      = 'W/m2'
          CASE ( 10 )
             NameSuffix = 'USTAR'        ! ND67, trcr 7: Friction vel
             Units      = 'm/s'
          CASE ( 11 )
             NameSuffix = 'Z0'           ! ND67, trcr 8: Roughness height
             Units      = 'm'
          CASE ( 12 )
             NameSuffix = 'PBL'          ! ND67, trcr 9: PBL height [m]?
             Units      = 'hPa'
          CASE ( 13 )
             NameSuffix = 'CLDFRC'       ! ND67, trcr 10: Column cld fraction
             Units      = '0-1'
          CASE ( 14 )
             NameSuffix = 'U10M'         ! ND67, trcr 11: U-wind @ 10m
             Units      = 'm/s'
          CASE ( 15 )
             NameSuffix = 'V10M'         ! ND67, trcr 12: V-wind @ 10m
             Units      = 'm/s'
          CASE ( 16 )
             NameSuffix = 'PS-PBL'       ! ND67, trcr 13
             Units      = 'hPa'
          CASE ( 17 )
             NameSuffix = 'ALBD'         ! ND67, trcr 14: Surface albedo
             Units      = 'unitless'
          CASE ( 18 )
             NameSuffix = 'PHIS'         ! ND67, trcr 15: Surface geopotential
             Units      = 'm'
          CASE ( 19 )
             NameSuffix = 'CLDTOP'       ! ND67, trcr 16: Cloud top level
             Units      = 'level'
          CASE ( 20 )
             NameSuffix = 'TROPP'        ! ND67, trcr 17: T'pause pressure
             Units      = 'hPa'
          CASE ( 21 )
             NameSuffix = 'SLP'          ! ND67, trcr 18: Sea level pressure
             Units      = 'hPa'
          CASE ( 22 )
             NameSuffix = 'TSKIN'        ! ND67, trcr 19: Surface skin temp
             Units      = 'K'
          CASE ( 23 )
             NameSuffix = 'PARDF'        ! ND67, trcr 20: Diffuse PAR
             Units      = 'W/m2'
          CASE ( 24 )
             NameSuffix = 'PARDR'        ! ND67, trcr 21: Direct PAR
             Units      = 'W/m2'
          CASE ( 25 )
             NameSuffix = 'GWETTOP'      ! ND67, trcr 22: Topsoil wetness [frac]
             Units      = 'unitless'
          CASE ( 26 )
             NameSuffix = 'EFLUX'        ! ND67, trcr 23: Latent heat flux
             Units      = 'W/m2'
          CASE DEFAULT
             IF ( N < Num2D ) THEN
                MSG = 'Num2D is less than number of named 2D MET diagnostics'
             ELSE
                MSG = 'Undefined 2D diagnostic case in MET diagnostic group'
             ENDIF
             CALL ERROR_STOP( MSG, LOC ) 
        END SELECT

       ! Diagnostic container info
       DiagnName = 'MET_2D_' // NameSuffix
       cId       = cId + 1
       SpaceDim  = 2

       ! Create container
       CALL Diagn_Create( am_I_Root,  HcoState,          &
                          Col       = Collection,        & 
                          cId       = cId,               &
                          cName     = TRIM( DiagnName ), &
                          SpaceDim  =  SpaceDim,         &
                          OutUnit   = TRIM( Units ),     &
                          OutOper   = TRIM( OutOper ),   &
                          RC        = RC )

       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot create 2D MET diagnostic ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF
    ENDDO

  END SUBROUTINE DiagnInit_Met
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_pb_emiss
!
! !DESCRIPTION: Subroutine DIAGNINIT\_PB\_EMISS initializes the Pb emissions 
!  diagnostic (aka ND01). Other ND01 species emissions (Rn and Be7) are 
!  handled within HEMCO.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_Pb_Emiss( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  21 Jan 2015 - E. Lundgren - Initial version
!  16 Jun 2016 - K. Travis   - Now define species ID's with the Ind_ function
!  01 Jul 2016 - R. Yantosca - Remove reference to Input_Opt%TRACER_NAME
!  19 Jul 2016 - R. Yantosca - Now pass State_Chm as an argument
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: Collection
    INTEGER            :: id_Pb

    ! Strings
    CHARACTER(LEN=15)  :: OutOper
    CHARACTER(LEN=60)  :: DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGNINIT_PB_EMISS (diagnostics_mod.F90)' 

    !=======================================================================
    ! DIAGNINIT_PB_EMISS begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GC_SUCCESS

    ! Pb210 species ID
    id_Pb      = Ind_('Pb')

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND01_OUTPUT_TYPE
    
    ! If the species number for lead is scheduled for output in input.geos, 
    ! then define the diagnostic container for 210Pb emissions.
    IF ( ANY ( Input_Opt%TINDEX(1,:) == id_Pb ) ) THEN

       !----------------------------------------------------------------
       ! Create containers for Pb emissions [kg/s]
       !----------------------------------------------------------------
       
       ! Diagnostic container name and id
       DiagnName = 'EMISS_' // TRIM( State_Chm%SpcData(id_Pb)%Info%Name )
       cId = cId + 1

       ! Create container
       CALL Diagn_Create( am_I_Root,  HcoState,          &
                          Col       = Collection,        & 
                          cId       = cId,               &
                          cName     = TRIM( DiagnName ), &
                          SpaceDim  =  3,                &
                          OutUnit   = 'kg/s',            &
                          OutOper   = TRIM( OutOper ),   &
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
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  23 Jan 2015 - E. Lundgren - Initial version
!  01 Jul 2016 - R. Yantosca - Remove references to Input_Opt%TRACER_NAME
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    INTEGER, PARAMETER   :: NumSpecies = 3        ! Does a var exist for this
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER              :: Collection, M
    
    ! Strings
    CHARACTER(LEN=15)    :: OutOper
    CHARACTER(LEN=60)    :: DiagnName
    CHARACTER(LEN=255)   :: MSG
    CHARACTER(LEN=255)   :: LOC = 'DIAGNINIT_RN_DECAY (diagnostics_mod.F90)' 

    ! Arrays
    INTEGER              :: SpeciesN   (NumSpecies) 
    CHARACTER(LEN=3)     :: SpeciesName(NumSpecies)


    !=======================================================================
    ! DIAGNINIT_RN_DECAY begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND02_OUTPUT_TYPE
    
    ! Assign array of species numbers corresponding to decaying species
    SpeciesN    = (/ id_Pb, id_Rn, id_Be7 /)
    SpeciesName = (/ 'Pb',  'Rn',  'Be7'  /)

    ! Loop over # of radon decay diagnostics
    DO M = 1, NumSpecies

       ! If the species number is scheduled for output in input.geos, 
       ! then define the diagnostic container for that species.
       IF ( ANY ( Input_Opt%TINDEX(2,:) == SpeciesN( M ) ) ) THEN

          !----------------------------------------------------------------
          ! Create containers for Rn/Pb/Be7 decay [kg/s]
          !----------------------------------------------------------------

          ! Diagnostic container name and id
          DiagnName = 'DECAY_' // TRIM( SpeciesName(M) )
          cId = cId + 1
          
          ! Create container
          CALL Diagn_Create( am_I_Root,  HcoState,          &
                             Col       = Collection,        & 
                             cId       = cId,               &
                             cName     = TRIM( DiagnName ), &
                             SpaceDim  =  3,                &
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
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
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
    RC = GC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND12_OUTPUT_TYPE

    ! Check if certain species listed for ND12 in input.geos???

    !----------------------------------------------------------------
    ! Create containers for fraction of BL occupied by level L [.]
    !----------------------------------------------------------------

    ! Diagnostic container name and id
    DiagnName = 'BL_FRAC'
    cId = cId + 1

    ! Create container
    CALL Diagn_Create( am_I_Root,  HcoState,          &
                       Col       = Collection,        & 
                       cId       = cId,               &
                       cName     = TRIM( DiagnName ), &
                       SpaceDim  =  3,                &
                       OutUnit   = '.' ,              &
                       OutOper   = TRIM( OutOper ),   &
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
  SUBROUTINE DiagnInit_CldConv_Flx( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  26 Jan 2015 - E. Lundgren - Initial version
!  01 Jul 2016 - R. Yantosca - Use State_Chm%nAdvect to replace N_TRACERS
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: Collection, N, NA
    CHARACTER(LEN=15)      :: OutOper
    CHARACTER(LEN=60)      :: DiagnName
    CHARACTER(LEN=255)     :: MSG
    CHARACTER(LEN=255)     :: LOC='DIAGNINIT_CLDCONV_FLX (diagnostics_mod.F90)'

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !=======================================================================
    ! DIAGNINIT_CLDCONV_FLX begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND14_OUTPUT_TYPE

    ! Check if certain species listed for ND14 in input.geos???
    
    ! Loop over advected species
    DO NA = 1, State_Chm%nAdvect
         
       ! Species ID 
       N = State_Chm%Map_Advect(NA)

       ! Entry in the species database
       SpcInfo => State_Chm%SpcData(N)%Info

       ! If this species number N is scheduled for output in input.geos, 
       ! then define the diagnostic containers for cloud convection mass change
       IF ( ANY( Input_Opt%TINDEX(14,:) == N ) ) THEN

          !----------------------------------------------------------------
          ! Create container for mass change due to cloud convection [kg/s]
          !----------------------------------------------------------------
      
          ! Diagnostic container name and id
          DiagnName = 'CLDCONV_FLX_' // TRIM( SpcInfo%Name )
          cId = cId + 1

          ! Create container
          CALL Diagn_Create( am_I_Root,  HcoState,          &
                             Col       = Collection,        & 
                             cId       = cId,               &
                             cName     = TRIM( DiagnName ), &
                             SpaceDim  =  3,                &
                             OutUnit   = 'kg/s' ,           &
                             OutOper   = TRIM( OutOper ),   &
                             RC        = RC )
      
          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF      
       ENDIF

       ! Free pointer
       SpcInfo => NULL()
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
  SUBROUTINE DiagnInit_BLMix_Flx( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  26 Jan 2015 - E. Lundgren - Initial version
!  01 Jul 2016 - R. Yantosca - Use State_Chm%nAdvect to replace N_TRACERS
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: Collection, N, NA

    ! Strings
    CHARACTER(LEN=15)      :: OutOper
    CHARACTER(LEN=60)      :: DiagnName
    CHARACTER(LEN=255)     :: MSG
    CHARACTER(LEN=255)     :: LOC = 'DIAGNINIT_BLMIX_FLX (diagnostics_mod.F90)' 

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !=======================================================================
    ! DIAGNINIT_BLMIX_FLX begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND15_OUTPUT_TYPE
    
    ! Check if certain species listed for ND15 in input.geos???

    ! Loop over only the advected species
    DO NA = 1, State_Chm%nAdvect

       ! Species ID
       N = State_Chm%Map_Advect(NA)

       ! Entry in the species database
       SpcInfo => State_Chm%SpcData(N)%Info
         
       ! If this species number N is scheduled for output in input.geos, 
       ! then define the diagnostic containers for BL mixing upward flux
       IF ( ANY( Input_Opt%TINDEX(15,:) == N ) ) THEN

          !----------------------------------------------------------------
          ! Create containers for boundary-layer mixing upward mass flux [kg/s]
          !----------------------------------------------------------------
      
          ! Diagnostic container name and id
          DiagnName = 'BLMIX_FLX_' // TRIM( SpcInfo%Name )
          cId = cId + 1

          ! Create container
          CALL Diagn_Create( am_I_Root,  HcoState,          &
                             Col       = Collection,        & 
                             cId       = cId,               &
                             cName     = TRIM( DiagnName ), &
                             SpaceDim  =  3,                &
                             OutUnit   = 'kg/s' ,           &
                             OutOper   = TRIM( OutOper ),   &
                             RC        = RC )
      
          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
       ENDIF

       ! Free pointer
       SpcInfo => NULL()
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
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
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
    RC = GC_SUCCESS

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
       CALL Diagn_Create( am_I_Root,  HcoState,          &
                          Col       = Collection,        & 
                          cId       = cId,               &
                          cName     = TRIM( DiagnName ), &
                          SpaceDim  =  3,                &
                          OutUnit   = '.' ,              &
                          OutOper   = TRIM( OutOper ),   &
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
  SUBROUTINE DiagnInit_Rain_Frac( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  26 Jan 2015 - E. Lundgren - Initial version
!  01 Jul 2016 - R. Yantosca - Use module variable nAdvect to replace N_TRACERS
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!  19 Jul 2016 - R. Yantosca - Bug fix: now declare N, nAdvect variables
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: Collection, M, NA, N

    ! Strings
    CHARACTER(LEN=15)      :: OutOper
    CHARACTER(LEN=60)      :: DiagnName
    CHARACTER(LEN=255)     :: MSG
    CHARACTER(LEN=255)     :: LOC='DIAGNINIT_RAIN_FRAC (diagnostics_mod.F90)'

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !=======================================================================
    ! DIAGNINIT_RAIN_FRAC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND17_OUTPUT_TYPE
    
    ! Check if certain species listed for ND17 in input.geos???

    ! Loop over only the advected species
    DO NA = 1, State_Chm%nAdvect
         
       ! Species ID
       N = State_Chm%Map_Advect(NA)

       ! Entry in the species database
       SpcInfo => State_Chm%SpcData(N)%Info

       ! If this species number N is scheduled for output in input.geos, 
       ! then define the diagnostic containers for rainout fraction
       IF ( ANY( Input_Opt%TINDEX(17,:) == N ) ) THEN

          DO M = 1, 2          
             SELECT CASE ( M )
                CASE ( 1 )
                !-----------------------------------------------------
                ! Name container for the fraction of soluble species
                ! lost to rainout (large-scale precipitation) [.]
                !-----------------------------------------------------
                DiagnName = 'RAIN_FRAC_LS_' // TRIM( SpcInfo%Name )
                CASE ( 2 )
               !------------------------------------------------------
               ! Name container for the fraction of soluble species
               ! lost to rainout (convective precipitation) [.]
               !------------------------------------------------------
               ! Diagnostic name
               DiagnName = 'RAIN_FRAC_CONV_' // TRIM( SpcInfo%Name )
                CASE DEFAULT
                   IF ( RC /= HCO_SUCCESS ) THEN
                      MSG = 'ND17 diagnostic name not defined.'
                      CALL ERROR_STOP( MSG, LOC ) 
                   ENDIF             
             END SELECT

             ! Diagnostic container id
             cId = cId + 1
      
             ! Create container
             CALL Diagn_Create( am_I_Root,  HcoState,          &
                                Col       = Collection,        & 
                                cId       = cId,               &
                                cName     = TRIM( DiagnName ), &
                                SpaceDim  =  3,                &
                                OutUnit   = '.' ,              &
                                OutOper   = TRIM( OutOper ),   &
                                RC        = RC )
         
             IF ( RC /= HCO_SUCCESS ) THEN
                MSG = 'Cannot create ND17 diagnostic: ' // TRIM(DiagnName)
                CALL ERROR_STOP( MSG, LOC ) 
             ENDIF 
          ENDDO     
       ENDIF

       ! Free pointer
       SpcInfo => NULL()
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
  SUBROUTINE DiagnInit_Wash_Frac( am_I_Root, Input_Opt, State_Chm, RC )
!
! !USES:
!
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  26 Jan 2015 - E. Lundgren - Initial version
!  01 Jul 2016 - R. Yantosca - Use module variable nAdvect to replace N_TRACERS
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: Collection, M, NA, N

    ! Strings
    CHARACTER(LEN=15)      :: OutOper
    CHARACTER(LEN=60)      :: DiagnName
    CHARACTER(LEN=255)     :: MSG
    CHARACTER(LEN=255)     :: LOC='DiagnInit_Wash_Frac (diagnostics_mod.F90)' 

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !=======================================================================
    ! DIAGNINIT_WASH_FRAC begins here!
    !=======================================================================
      
    ! Assume successful return
    RC         = GC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND18_OUTPUT_TYPE
 
    ! Check if certain species listed for ND18 in input.geos???

    ! Loop over advected species
    DO NA = 1, State_Chm%nAdvect
         
       ! Species ID
       N = State_Chm%Map_Advect(NA)

       ! Entry in the species database
       SpcInfo => State_Chm%SpcData(N)%Info

       ! If this species number N is scheduled for output in input.geos, 
       ! then define the diagnostic containers for washout fraction
       IF ( ANY( Input_Opt%TINDEX(18,:) == N ) ) THEN

          DO M = 1, 2          
             SELECT CASE ( M )
                CASE ( 1 )
                   !-----------------------------------------------------
                   ! Name container for the fraction of soluble species
                   ! lost to washout (large-scale precipitation) [.]
                   !-----------------------------------------------------
                   DiagnName = 'WASH_FRAC_LS' // TRIM( SpcInfo%Name )
                CASE ( 2 )
                   !------------------------------------------------------
                   ! Name container for the fraction of soluble species
                   ! lost to washout (convective precipitation) [.]
                   !------------------------------------------------------
                   DiagnName = 'WASH_FRAC_CONV_' // TRIM( SpcInfo%Name )
                CASE DEFAULT
                   IF ( RC /= HCO_SUCCESS ) THEN
                      MSG = 'ND18 diagnostic name not defined.'
                      CALL ERROR_STOP( MSG, LOC ) 
                   ENDIF             
             END SELECT

             ! Diagnostic container id
             cId = cId + 1

             ! Create container
             CALL Diagn_Create( am_I_Root,  HcoState,          &
                                Col       = Collection,        & 
                                cId       = cId,               &
                                cName     = TRIM( DiagnName ), &
                                SpaceDim  =  3,                &
                                OutUnit   = '.' ,              &
                                OutOper   = TRIM( OutOper ),   &
                                RC        = RC )
         
             IF ( RC /= HCO_SUCCESS ) THEN
                MSG = 'Cannot create ND18 diagnostic: ' // TRIM(DiagnName)
                CALL ERROR_STOP( MSG, LOC ) 
             ENDIF
          ENDDO
       ENDIF

       ! Free pointer
       SpcInfo => NULL()
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
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
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
    RC = GC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%ND19_OUTPUT_TYPE
    
    ! Check if certain species listed for ND19 in input.geos???

    !----------------------------------------------------------------
    ! Create containers for CH4 removal by OH [kg CH4]
    !----------------------------------------------------------------

    ! Diagnostic container name and id
    DiagnName = 'CH4_LOSS'
    cId = cId + 1

    ! Create container
    CALL Diagn_Create( am_I_Root,  HcoState,          &
                       Col       = Collection,        & 
                       cId       = cId,               &
                       cName     = TRIM( DiagnName ), &
                       SpaceDim  =  3,                &
                       OutUnit   = 'kg CH4' ,         &
                       OutOper   = TRIM( OutOper ),   &
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
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
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
    RC = GC_SUCCESS

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
    CALL Diagn_Create( am_I_Root,  HcoState,          &
                       Col       = Collection,        & 
                       cId       = cId,               &
                       cName     = TRIM( DiagnName ), &
                       SpaceDim  =  2,                &
                       OutUnit   = '.' ,              &
                       OutOper   = TRIM( OutOper ),   &
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
! !IROUTINE: diagnupdate_met
!
! !DESCRIPTION: Subroutine DIAGNUPDATE\_MET updates the meteorology state
!  diagnostics (aka ND31, ND55, ND57, ND66, ND67, and ND68). Only the subset
!  of Met diagnostics which save State_Met data or are calculated from
!  State_Met data are updated in this subroutine. All other Met diagnostics
!  are updated within Met-field read modules.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnUpdate_Met( am_I_Root, Input_Opt, State_Met, RC )
!
! !USES:
!
    USE HCO_Diagn_Mod,      ONLY : Diagn_Update
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants,      ONLY : XNUMOLAIR, SCALE_HEIGHT
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN   ) :: State_Met  ! Met state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  21 Jan 2015 - E. Lundgren - Initial version
!  15 Jan 2016 - E. Lundgren - Revise for all MET-related diagnostics
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, N, Collection, Num2D, Num3D, HCRC
    CHARACTER(LEN=30)  :: NameSuffix, DiagnName
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DiagnUpdate_Met (diagnostics_mod.F90)'
    REAL(fp), TARGET   :: Temp2D( IIPAR, JJPAR )
    REAL(fp), TARGET   :: Temp3D( IIPAR, JJPAR, LLPAR )
    REAL(fp), POINTER  :: Ptr2D(:,:)
    REAL(fp), POINTER  :: Ptr3D(:,:,:)

    !=======================================================================
    ! DIAGNUPDATE_MET begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GC_SUCCESS

    ! Initialize pointers
    Ptr2D => NULL()
    Ptr3D => NULL()

    ! Get diagnostic parameters from the Input_Opt object
    ! This is not currently used, but keep for possible later use/editing
    Collection = Input_Opt%DIAG_COLLECTION

    ! For now, only update if GEOS-FP or MERRA2 since diagnostics
    ! for those are stored in State_Met and are therefore accessible here.
#if defined( GEOS_FP ) || defined ( MERRA2 )

    ! Set number of 3D and 2D MET diagnostics that are updated 
    Num2d = 23
    Num3d = 14

    !----------------------------------------------------------------
    ! Update 2D containers
    !----------------------------------------------------------------
  
    ! Loop over 2D diagnostics within Met category
    DO N = 1, Num2D

       ! Zero the temporary array in case it is used
       Temp2D = 0.0e+0_fp

       SELECT CASE ( N )
          CASE ( 1 )
             NameSuffix = 'HFLUX'   ! ND67, trcr 1: Sens heat flux [W/m2]
             Ptr2D => State_Met%HFLUX  
          CASE ( 2 )
             NameSuffix = 'RADSWG'  ! ND67, trcr 2: SW rad @ sfc [W/m2]
             Ptr2D => State_Met%SWGDN  
          CASE ( 3 )
             NameSuffix = 'PREACC'  ! ND67, trcr 3: Tot prec [kg/m2/s]
             Ptr2D => State_Met%PRECTOT 
          CASE ( 4 )
             NameSuffix = 'PRECON'  ! ND67, trcr 4: Sfc conv prec[kg/m2/s]
             Ptr2D => State_Met%PRECCON 
          CASE ( 5 )
             NameSuffix = 'TS'      ! ND67, trcr 5: T @ 2m height [K]
             Ptr2D => State_Met%TS   
          CASE ( 6 )
             NameSuffix = 'RADSWT'  ! ND67, trcr 6
             Ptr2D => Temp2D        ! BPCH is set to 0 for GEOS-FP/MERRA2
          CASE ( 7 )
             NameSuffix = 'USTAR'   ! ND67, trcr 7: Friction vel [m/s]
             Ptr2D => State_Met%USTAR   
          CASE ( 8 )
             NameSuffix = 'Z0'      ! ND67, trcr 8: Roughness height [m]
             Ptr2D => State_Met%Z0    
          CASE ( 9 )
             NameSuffix = 'PBL'     ! ND67, trcr 9: PBL height [m]
             Ptr2D => State_Met%PBLH  
          CASE ( 10 )
             NameSuffix = 'CLDFRC'  ! ND67, trcr 10: Column cld fraction
             Ptr2D => State_Met%CLDFRC  
          CASE ( 11 )
             NameSuffix = 'U10M'    ! ND67, trcr 11: U-wind @ 10m [m/s]
             Ptr2D => State_Met%U10M   
          CASE ( 12 )
             NameSuffix = 'V10M'    ! ND67, trcr 12: V-wind @ 10m [m/s]
             Ptr2D => State_Met%V10M  
          CASE ( 13 )
             NameSuffix = 'PS-PBL'  ! ND67, trcr 13: Boundary layer top pressure
             Temp2D = State_Met%PEDGE(:,:,1) *      &             ! [hPa]
                      EXP( -State_Met%PBLH(:,:) / SCALE_HEIGHT ) 
             Ptr2D => Temp2D
          CASE ( 14 )
             NameSuffix = 'ALBD'    ! ND67, trcr 14: Sfc albedo [unitless]
             Ptr2D => State_Met%ALBD   
          CASE ( 15 )
             NameSuffix = 'PHIS'    ! ND67, trcr 15: Sfc geopotential [m]
             Ptr2D => State_Met%PHIS
          CASE ( 16 )
             NameSuffix = 'CLDTOP'  ! ND67, trcr 16: Max cloud top level
             Temp2D = FLOAT( State_Met%CLDTOPS )
             Ptr2D => Temp2D
          CASE ( 17 )
             NameSuffix = 'TROPP'   ! ND67, trcr 17: T'pause pressure
             Ptr2D => State_Met%TROPP !              [hPa]
          CASE ( 18 )
             NameSuffix = 'SLP'     ! ND67, trcr 18: Sea level prs [hPa]
             Ptr2D => State_Met%SLP 
          CASE ( 19 )
             NameSuffix = 'TSKIN'   ! ND67, trcr 19: Sfc skin temp [K]
             Ptr2D => State_Met%TSKIN   
          CASE ( 20 )
             NameSuffix = 'PARDF'   ! ND67, trcr 20: Diffuse PAR [W/m2]
             Ptr2D => State_Met%PARDF 
          CASE ( 21 )
             NameSuffix = 'PARDR'   ! ND67, trcr 21: Direct PAR [W/m2]
             Ptr2D => State_Met%PARDR   
          CASE ( 22 )
             NameSuffix = 'GWETTOP' ! ND67, trcr 22: Topsoil wetness
             Ptr2D => State_Met%GWETTOP !            [frac]
          CASE ( 23 )
             NameSuffix = 'EFLUX'   ! ND67, trcr 23: Latent heat flx
             Ptr2D => State_Met%EFLUX !              [W/m2] 
          CASE DEFAULT
             IF ( N > Num2D ) THEN
                MSG = 'Undefined 2D diagnostic in MET diagnostic group'
             ENDIF
             CALL ERROR_STOP( MSG, LOC ) 
       END SELECT

       ! Diagnostic info
       DiagnName = 'MET_2D_' // NameSuffix

       ! Update diagnostics
       IF ( ASSOCIATED(Ptr2D) ) THEN 
          CALL Diagn_Update( am_I_Root, HcoState,         &
                             cName   = TRIM( DiagnName ), &
                             Array2D = Ptr2D,             &
                             RC      = HCRC )

          ! Free the pointer
          Ptr2D => NULL()
 
          IF ( HCRC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot update 2D MET diagnostic ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
       ENDIF
    ENDDO

    !----------------------------------------------------------------
    ! Update 3D containers
    !----------------------------------------------------------------
    ! Loop over 3D diagnostics within Met category
    DO N = 1, Num3D

       ! Zero the temporary array in case it is used
       Temp3D = 0.0e+0_fp

       SELECT CASE ( N ) 
          CASE ( 1 )
             NameSuffix = 'PEDGE'        ! ND31, trcr 1
             Ptr3D => State_Met%PEDGE
          CASE ( 2 )
             NameSuffix = 'THETA'        ! Potential temp, ND57, trcr 1
             !$OMP PARALLEL DO &
             !$OMP DEFAULT( SHARED ) &
             !$OMP PRIVATE( I, J, L )
             DO L = 1, LLPAR
             DO J = 1, JJPAR
             DO I = 1, IIPAR
                Temp3D(I,J,L) = State_Met%T(I,J,L)               &
                                * ( State_MET%PEDGE(I,J,1)       & 
                                / State_Met%PMID(I,J,L) )**0.286
             ENDDO
             ENDDO
             ENDDO
             !$OMP END PARALLEL DO
             Ptr3D => Temp3D
          CASE ( 3 )
             NameSuffix = 'UWND'         ! ND66, trcr 1: Zonal wind [m/s]
             Ptr3D => State_Met%U 
          CASE ( 4 )
             NameSuffix = 'VWND'         ! ND66, trcr 2: Meridional wind [m/s]
             Ptr3D => State_Met%V 
          CASE ( 5 )
             NameSuffix = 'TMPU'         ! ND66, trcr 3: Temperature [K]
             Ptr3D => State_Met%TMPU1
          CASE ( 6 )
             NameSuffix = 'TMPU'         ! ND66, trcr 3: Temperature [K]
             Ptr3D => State_Met%TMPU2
          CASE ( 7 )
             NameSuffix = 'SPHU'         ! ND66, trcr 4: Specific humidity
             Ptr3D => State_Met%SPHU1                   ! [g H2O/kg air]
          CASE ( 8 )
             NameSuffix = 'SPHU'         ! ND66, trcr 4: Specific humidity
             Ptr3D => State_Met%SPHU2                   ! [g H2O/kg air]
          CASE ( 9 )
             NameSuffix = 'CLDMAS'       ! ND66, trcr 5: Convective mass flux
             Ptr3D => State_Met%CMFMC                  ! [kg/m2/s]
          CASE ( 10 )
             NameSuffix = 'DTRAIN'       ! ND66, trcr 6: Detrainment flux 
             Ptr3D => State_Met%DTRAIN                 ! [kg/m2/s]            
          CASE ( 11 )
             NameSuffix = 'BXHEIGHT'     ! ND68, trcr 1: grid box height [m]
             Ptr3D => State_Met%BXHEIGHT 
          CASE ( 12 )
             NameSuffix = 'DRYAIRMASS'   ! ND68, trcr 2: grid box dry air mass
             Ptr3D => State_Met%AD       !               [kg dry air]
          CASE ( 13 )
             NameSuffix = 'AVGW'         ! ND68, trcr 3: H2O vapor mixing ratio
             Ptr3D => State_Met%AVGW     !               [v/v dry air]
          CASE ( 14 )
             NameSuffix = 'NAIR'         ! ND68, trcr 4: Air number density
             Temp3D = State_Met%AIRDEN * XNUMOLAIR  !    [molec dry air/m3]
             Ptr3D => Temp3D
          CASE DEFAULT
             IF ( N < Num3D ) THEN
                MSG = 'Num3D is less than number of named 3D MET diagnostics'
             ELSE
                MSG = 'Undefined 3D diagnostic in MET diagnostic group'
             ENDIF
             CALL ERROR_STOP( MSG, LOC ) 
       END SELECT

       ! Diagnostic name
       DiagnName = 'MET_3D_' // NameSuffix

       ! Update diagnostics
       CALL Diagn_Update( am_I_Root, HcoState,          &
                          cName   = TRIM( DiagnName ),  &
                          Array3D = Ptr3D,              &
                          RC      = RC )

       ! Free the pointer
       Ptr3D => NULL()
 
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot update MET diagnostic ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF
    ENDDO

#endif

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
  SUBROUTINE DiagnUpdate_Transport_Flux( am_I_Root, FLUX_EW,   FLUX_NS,       &
                                         FLUX_VERT, Input_Opt, State_Chm, RC )

!
! !USES:
!
   USE HCO_Diagn_Mod,      ONLY : Diagn_Update
   USE Input_Opt_Mod,      ONLY : OptInput    
   USE Species_Mod,        ONLY : Species
   USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
   LOGICAL,        INTENT(IN)    :: am_I_Root          ! Are we on root CPU?
   REAL(fp),       INTENT(IN)    :: FLUX_EW(:,:,:,:)   ! zonal mass flux
   REAL(fp),       INTENT(IN)    :: FLUX_NS(:,:,:,:)   ! meridional flux
   REAL(fp),       INTENT(IN)    :: FLUX_VERT(:,:,:,:) ! vertical mass flux
   TYPE(OptInput), INTENT(IN)    :: Input_Opt          ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
   TYPE(ChmState), INTENT(INOUT) :: State_Chm          ! Chemistry State object
!
!
! !OUTPUT PARAMETERS:
!
   INTEGER,        INTENT(OUT)   :: RC                 ! Success or failure
!
! !REMARKS:
!
! !REVISION HISTORY:
!  12 Feb 2015 - E. Lundgren - Initial version
!  19 Jan 2016 - E. Lundgren - Updated to include all transport flux diags
!  30 Jun 2016 - R. Yantosca - Now get advected species IDs from State_Chm
!  20 Jul 2016 - R. Yantosca - Remove references to NNPAR
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
   ! Scalars
   INTEGER                :: N, M, NumLevels, NA

   ! Strings
   CHARACTER(LEN=60)      :: DiagnPrefix, DiagnName
   CHARACTER(LEN=255)     :: MSG
   CHARACTER(LEN=255)     :: LOC = 'DiagnUpdate_Transport_Flux'  &
                               // ' (Transport_mod.F)'

   ! Array to hold input mass flux array with reverse-order levels
   ! and only levels that have data (ie. less than LLPAR if input.geos
   ! levels for the diagnostic are less than LLPAR).
   REAL(fp), TARGET       :: DiagnArray( IIPAR, JJPAR, Input_Opt%LD24,    &
                                                       State_Chm%nAdvect )

   ! Pointers
   REAL(fp),      POINTER :: Ptr3D(:,:,:)
   TYPE(Species), POINTER :: SpcInfo

   !=================================================================
   ! DIAGNUPDATE_TRANSPORT_FULX begins here!
   !=================================================================
   
    ! Assume successful return
    RC        =  GC_SUCCESS

    ! Get number of levels for transport diagnostics
    NumLevels =  Input_Opt%LD24

    ! Initialize pointers
    Ptr3d     => NULL()
    SpcInfo   => NULL()

    ! Loop over diagnostics
    DO M = 1, 3

       ! Set diagnostic name suffix and reverse-order the levels such that 
       ! index 1 corresponds to surface. 
       SELECT CASE ( M )
          CASE ( 1 )
             DiagnPrefix = 'TRANSPORT_FLX_EW_'
             DiagnArray  = FLUX_EW(:,:,LLPAR:LLPAR-NumLevels+1:-1,:)
          CASE ( 2 )
             DiagnPrefix = 'TRANSPORT_FLX_NS_'
             DiagnArray  = FLUX_NS(:,:,LLPAR:LLPAR-NumLevels+1:-1,:)
          CASE ( 3 )
             DiagnPrefix = 'TRANSPORT_FLX_VERT_'
             DiagnArray  = FLUX_VERT(:,:,LLPAR:LLPAR-NumLevels+1:-1,:)
       END SELECT

       ! Loop over only the advected species
       DO NA = 1, State_Chm%nAdvect
      
          ! Get the species ID from the advected species ID
          N = State_Chm%Map_Advect(NA)

          ! If this species number N is scheduled for output in input.geos, 
          ! then update its diagnostic container for transport flux
          IF ( ANY( Input_Opt%TINDEX(24,:) == N ) ) THEN
            
             !----------------------------------------------------------------
             ! Update diagnostic container
             !----------------------------------------------------------------
         
             ! Point to the corresponding entry in the species database
             SpcInfo   => State_Chm%SpcData(N)%Info

             ! Diagnostic container name
             DiagnName =  TRIM( DiagnPrefix ) // TRIM( SpcInfo%Name )
       
             ! Point to the array with levels corrected
             Ptr3D     => DiagnArray(:,:,:,NA)
       
             ! Create container
             CALL Diagn_Update( am_I_Root,                        &
                                HcoState  = HcoState,             &
                                cName     = TRIM( DiagnName ),    &
                                Array3D   = Ptr3D,                &
                                RC        = RC )
       
             ! Free the pointers
             Ptr3D   => NULL()
             SpcInfo => NULL()
       
             ! Stop with error if the diagnostics update was unsuccessful.
             IF ( RC /= HCO_SUCCESS ) THEN
                MSG = 'Cannot update diagnostic: ' // TRIM( DiagnName )
                CALL ERROR_STOP( MSG, LOC ) 
             ENDIF 
          ENDIF
       ENDDO
    ENDDO

   END SUBROUTINE DiagnUpdate_Transport_Flux
!EOC
#if defined( DIAG_DEVEL )
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DiagnInit_Species_Emis
!
! !DESCRIPTION: Subroutine DiagnInit\_Species\_Emis initializes diagnostics for 
!  total species emissions diagnostics. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_Species_Emis( am_I_Root, Input_Opt,     &
                                    State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE HCO_INTERFACE_MOD,  ONLY : GetHcoID
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt   ! Input Options object
    TYPE(MetState),   INTENT(IN   )  :: State_Met  ! Met state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)    :: RC         ! Failure or success
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  05 Mar 2015 - C. Keller   - Initial version
!  16 Jun 2016 - K. Travis   - Now define species ID's with the Ind_ function 
!  01 Jul 2016 - R. Yantosca - Use State_Chm%nAdvect to replace N_TRACERS
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: ID
    INTEGER                :: Collection, N, NA

    ! Strings
    CHARACTER(LEN=15)      :: OutOper
    CHARACTER(LEN=60)      :: DiagnName
    CHARACTER(LEN=255)     :: MSG
    CHARACTER(LEN=255)     :: LOC='DiagnInit_Species_Emis (diagnostics_mod.F90)'

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !=======================================================================
    ! DIAGNINIT_SPECIES_EMIS begins here!
    !=======================================================================
      
    ! Assume successful return
    RC = GC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    ! Use same output frequency and operations as for species concentrations.
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = Input_Opt%SPECIES_EMIS_OUTPUT_TYPE
 
    ! Loop over # of species 
    DO NA = 1, State_Chm%nAdvect

       ! Species ID
       N = State_Chm%Map_Advect(NA)

       ! Entry in the species database
       SpcInfo => State_Chm%SpcData(N)%Info
       
       ! HEMCO ID
       ID = GetHcoID( SpcId = N )
 
       ! Restrict diagnostics to these species
       IF ( N /= Ind_('NO'   ) .AND. N /= Ind_('CO'    ) .AND. &
            N /= Ind_('ALK4' ) .AND. N /= Ind_('ISOP'  ) .AND. &
            N /= Ind_('HNO3' ) .AND. N /= Ind_('ACET'  ) .AND. &
            N /= Ind_('MEK'  ) .AND. N /= Ind_('ALD2'  ) .AND. &
            N /= Ind_('PRPE' ) .AND. N /= Ind_('C3H8'  ) .AND. &
            N /= Ind_('C2H6' ) .AND. N /= Ind_('DMS'   ) .AND. &
            N /= Ind_('SO2'  ) .AND. N /= Ind_('SO4'   ) .AND. &
            N /= Ind_('NH3'  ) .AND. N /= Ind_('BCPI'  ) .AND. &
            N /= Ind_('OCPI' ) .AND. N /= Ind_('BCPO'  ) .AND. &
            N /= Ind_('OCPO' ) .AND. N /= Ind_('DST1'  ) .AND. &
            N /= Ind_('DST2' ) .AND. N /= Ind_('DST3'  ) .AND. &
            N /= Ind_('DST4' ) .AND. N /= Ind_('SALA'  ) .AND. &
            N /= Ind_('SALC' ) .AND. N /= Ind_('Br2'   ) .AND. &
            N /= Ind_('BrO'  ) .AND. N /= Ind_('CH2Br2') .AND. &
            N /= Ind_('CH3Br') .AND. N /= Ind_('O3'    )        ) THEN
          ID = -1
       ENDIF
 
       ! If this is an emission species, add diagnostics for emissions. 
       IF ( ID > 0 ) THEN 

          !----------------------------------------------------------------
          ! Create container for emission flux (kg/s) 
          !----------------------------------------------------------------

          ! Diagnostic container name and id
          DiagnName = 'SPECIES_EMIS_' // TRIM( SpcInfo%Name )
          cId = 10000 + ID

          ! Create container
          CALL Diagn_Create( am_I_Root,                     &
                             HcoState  = HcoState,          & 
                             Col       = Collection,        & 
                             cId       = cId,               &
                             cName     = TRIM( DiagnName ), &
                             HcoID     = ID,                &
                             SpaceDim  =  3,                &
                             OutUnit   = 'kg/s',            &
                             OutOper   = TRIM( OutOper ),   &
                             RC        = RC                  )

          IF ( RC /= HCO_SUCCESS ) THEN
             MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
             CALL ERROR_STOP( MSG, LOC ) 
          ENDIF
       ENDIF

       ! Free pointer
       SpcInfo => NULL()
       
    ENDDO

  END SUBROUTINE DiagnInit_Species_Emis
!EOC
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
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : Ind_
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  07 Jul 2015 - C. Keller   - Initial version 
!  16 Jun 2016 - K. Travis   - Now define species ID's with the Ind_ function 
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
    RC = GC_SUCCESS

    ! Nothing to do if O3 is not a species
    IF ( Ind_('O3') <= 0 ) RETURN

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
                          HcoState  = HcoState,          & 
                          Col       = Collection,        & 
                          cId       = cId,               & 
                          cName     = TRIM( DiagnName ), &
                          SpaceDim  =  2,                &
                          OutUnit   = 'dobson',          &
                          OutOper   = TRIM( OutOper ),   &
                          OkIfExist = .TRUE.,            &
                          RC        = RC )
      IF ( RC /= HCO_SUCCESS ) THEN
         CALL ERROR_STOP( 'Cannot create diagnostics '//TRIM(DiagnName), LOC ) 
      ENDIF  
 
   ENDDO
   
  END SUBROUTINE DiagnInit_Dobson
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
    USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_TROP
    USE HCO_Diagn_Mod,      ONLY : Diagn_Update
    USE HCO_Diagn_Mod,      ONLY : DiagnCont_Find
    USE HCO_TYPES_MOD,      ONLY : DiagnCont
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants,      ONLY : AIRMW, AVO,   g0
    USE PRESSURE_MOD,       ONLY : GET_PEDGE
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root  ! Is this the root CPU?!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
    TYPE(MetState), INTENT(IN   ) :: State_Met  ! Met state
    TYPE(ChmState), INTENT(IN   ) :: State_Chm  ! Chemistry state 
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
! 
! !REVISION HISTORY: 
!  07 Jul 2015 - C. Keller   - Initial version 
!  27 Feb 2016 - C. Keller   - Deactivate OMP statement
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  16 Jun 2016 - K. Travis   - Now define species ID's with the Ind_ function 
!  22 Jun 2016 - M. Yannetti - Replace TCVV with species db MW and phys constant
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    TYPE(DiagnCont), POINTER :: DgnPtr

    INTEGER                  :: I, J, L, id_O3

    REAL(fp)                 :: constant
    REAL(fp)                 :: DP, O3vv, DU, MW_O3
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
    RC = GC_SUCCESS

    ! Initialize
    DgnPtr => NULL()
    id_O3  = Ind_('O3')
    MW_O3  = State_Chm%SpcData(id_O3)%Info%emMW_g

    ! Nothing to do if O3 is not a species
    IF ( id_O3 <= 0 ) RETURN

    ! On first call, check if any of the two diagnostics is defined
    IF ( FIRST ) THEN
       CALL DiagnCont_Find( HcoState%Diagn, -1, -1, -1, -1, -1, &
                           'O3_COLUMN', -1, TotDiagn, DgnPtr )
       DgnPtr => NULL()

       CALL DiagnCont_Find( HcoState%Diagn, -1, -1, -1, -1, -1, &
                           'O3_TROPCOLUMN', -1, TropDiagn, DgnPtr ) 
                            

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
    ! The following OMP loop does produce non-reproducable results, 
    ! presumably due to the call of GET_PEDGE or ITS_IN_THE_TROP? 
    ! Deactivate for now (ckeller, 2/27/16) 
#if defined(0)
    !$OMP PARALLEL DO &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, DP, O3vv, DU )
#endif
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR

       ! Pressure difference in hPa
       DP = GET_PEDGE(I,J,L) - GET_PEDGE(I,J,L+1)

       ! Ozone in v/v
       O3vv = State_Chm%Species(I,J,L,id_O3) &
            * (AIRMW / MW_O3) &
            / State_Met%AD(I,J,L)

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
#if defined(0)
    !$OMP END PARALLEL DO
#endif

    ! Update diagnostics
    IF ( TotDiagn ) THEN
       CALL Diagn_Update( am_I_Root, HcoState, cName = 'O3_COLUMN',     &
                          Array2D = TOTO3,  RC=RC )
    ENDIF
    IF ( TropDiagn ) THEN
       CALL Diagn_Update( am_I_Root, HcoState, cName = 'O3_TROPCOLUMN', &
                          Array2D = TROPO3, RC=RC )
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
! !IROUTINE: diagninit_KPP_rates
!
! !DESCRIPTION: Subroutine DIAGNINIT\_KPP_RATES initializes the KPP-based
! diagnostics arrays. It is associated (at this point) with the
! ND45 diagnostic.
!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_KPP_Rates( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE gckpp_Parameters
    USE gckpp_Monitor
    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
    USE Input_Opt_Mod,      ONLY : OptInput
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
!  09 Nov 2016 - M. Long     - Initial version, designed to output
!                              ALL reaction rates (including Phot. & HET)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: cID,      Collection, D, N
    CHARACTER(LEN=15)  :: OutOper
    CHARACTER(LEN=155) :: DiagnName, DiagnNameBase
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGNINIT_KPP_RATES (diagnostics_mod.F90)' 
    
    !=======================================================================
    ! DIAGNINIT_KPP_RATES begins here!
    !=======================================================================

    ! Assume successful return
    RC = GC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = 'Instantaneous'

    ! Loop over # of rates
    DO D = 1, NREACT
       !----------------------------------------------------------------
       ! Create containers for reaction rates (mcl/cc/s)
       ! CURRENTLY initializes ALL rates.
       !----------------------------------------------------------------

       ! Diagnostic name
       DiagnNameBase = TRIM( ADJUSTL(EQN_NAMES(D)) )
       write(DiagnName,'(a)')'RR_' // TRIM(DiagnNameBase)
       cId = cId + 1

       ! Create container
       CALL Diagn_Create( am_I_Root, HcoState,           &
                          Col       = Collection,        & 
                          cId       = cId,               &
                          cName     = TRIM( DiagnName ), &
                          SpaceDim  =  2,                &
                          OutUnit   = 'cm-3 s-1',        &
                          OutOper   = TRIM( OutOper   ), &
                          RC        = RC )
       
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF
       
    ENDDO
 
    ! Diagnostic name
    DiagnName = 'T'
    cId = cId + 1

    ! Create container
    CALL Diagn_Create( am_I_Root, HcoState,           &
                       Col       = Collection,        & 
                       cId       = cId,               &
                       cName     = TRIM( DiagnName ), &
                       SpaceDim  =  3,                &
                       OutUnit   = 'K',               &
                       OutOper   = TRIM( OutOper   ), &
                       RC        = RC )
    
    IF ( RC /= HCO_SUCCESS ) THEN
       MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
       CALL ERROR_STOP( MSG, LOC ) 
    ENDIF

    !--------------------------------------
    ! Sticking coefficients from HetRates
    !--------------------------------------
    ! Create container
    cId = cId + 1
    CALL Diagn_Create( am_I_Root, HcoState,&
         Col       = Collection,           & 
         cId       = cId,                  &
         cName     = TRIM( 'STK_N2O5' ),   &
         SpaceDim  =  3,                   &
         OutUnit   = '1',                  &
         OutOper   = TRIM( OutOper   ),    &
         RC        = RC )
    
    cId = cId + 1
    CALL Diagn_Create( am_I_Root, HcoState,&
         Col       = Collection,           & 
         cId       = cId,                  &
         cName     = TRIM( 'STK_HBr' ),    &
         SpaceDim  =  3,                   &
         OutUnit   = '1',                  &
         OutOper   = TRIM( OutOper   ),    &
         RC        = RC )
    
    cId = cId + 1
    CALL Diagn_Create( am_I_Root, HcoState,&
         Col       = Collection,           & 
         cId       = cId,                  &
         cName     = TRIM( 'STK_HBrICE' ), &
         SpaceDim  =  3,                   &
         OutUnit   = '1',                  &
         OutOper   = TRIM( OutOper   ),    &
         RC        = RC )
    

  END SUBROUTINE DiagnInit_KPP_Rates
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: diagninit_KPP_Spec
!
! !DESCRIPTION: Subroutine DIAGNINIT\_KPP initializes the KPP-based
! diagnostics arrays. It is associated (at this point) with the
! ND45 diagnostic.
!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DiagnInit_KPP_Spec( am_I_Root, Input_Opt, RC )
!
! !USES:
!
    USE gckpp_Parameters
    USE gckpp_Monitor
    USE Input_Opt_Mod,      ONLY : OptInput
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
!  09 Nov 2016 - M. Long     - Initial version, designed to output
!                              ALL reaction rates (including Phot. & HET)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: Collection, D, N
    CHARACTER(LEN=15)  :: OutOper
    CHARACTER(LEN=155) :: DiagnName, DiagnNameBase
    CHARACTER(LEN=255) :: MSG
    CHARACTER(LEN=255) :: LOC = 'DIAGNINIT_KPP_SPEC (diagnostics_mod.F90)' 
    
    !=======================================================================
    ! DIAGNINIT_KPP_SPEC begins here!
    !=======================================================================

    ! Assume successful return
    RC = GC_SUCCESS

    ! Get diagnostic parameters from the Input_Opt object
    Collection = Input_Opt%DIAG_COLLECTION
    OutOper    = 'Instantaneous'

    ! Loop over # of rates
    DO D = 1, NSPEC
       !----------------------------------------------------------------
       ! Create containers for reaction rates (mcl/cc/s)
       ! CURRENTLY initializes ALL rates.
       !----------------------------------------------------------------

       ! Diagnostic name
       DiagnNameBase = TRIM( ADJUSTL(SPC_NAMES(D)) )
       write(DiagnName,'(a)')'SPC_' // TRIM(DiagnNameBase)
       cId = cId + 1

       ! Create container
       CALL Diagn_Create( am_I_Root, HcoState,           &
                          Col       = Collection,        & 
                          cId       = cId,               &
                          cName     = TRIM( DiagnName ), &
                          SpaceDim  =  2,                &
                          OutUnit   = 'cm-3',            &
                          OutOper   = TRIM( OutOper   ), &
                          RC        = RC )
       
       IF ( RC /= HCO_SUCCESS ) THEN
          MSG = 'Cannot create diagnostics: ' // TRIM(DiagnName)
          CALL ERROR_STOP( MSG, LOC ) 
       ENDIF
       
    ENDDO
 
  END SUBROUTINE DiagnInit_KPP_Spec
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
  SUBROUTINE TotalsToLogfile( am_I_Root, Input_Opt, State_Chm, RC ) 
!
! !USES:
!
    USE HCO_INTERFACE_MOD,  ONLY : GetHcoID
    USE HCO_Diagn_Mod,      ONLY : Diagn_TotalGet
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE Species_Mod,        ONLY : Species
    USE TIME_MOD,           ONLY : GET_YEAR, GET_MONTH 
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input Options object
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)    :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  15 Mar 2015 - C. Keller   - Initial version 
!  18 Aug 2016 - M. Sulprizio- Use State_Chm%nAdvect to replace N_TRACERS
!EOP
!------------------------------------------------------------------------------
!BOC
    INTEGER       :: N, NA, cId, HCRC
    INTEGER       :: YEAR, MONTH 
    INTEGER, SAVE :: SAVEMONTH = -999
    REAL(sp)      :: TOTAL
    LOGICAL       :: FOUND

    ! Pointers
    TYPE(Species), POINTER :: SpcInfo

    !=================================================================
    ! TotalsToLogfile begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

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

       ! Loop over advected species
       DO NA = 1, State_Chm%nAdvect

          N = State_Chm%Map_Advect(NA)

          ! Entry in the species database
          SpcInfo => State_Chm%SpcData(N)%Info

          ! Only if it's a HEMCO species...
          cId = GetHcoID( SpcId=N )
          IF ( cId <= 0 ) CYCLE
          
          ! Define diagnostics ID
          cId = 10000 + cId

          ! Get the total [kg] 
          CALL Diagn_TotalGet( am_I_Root, HcoState%Diagn,         & 
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
          WRITE(6,101) TRIM( SpcInfo%Name ), Total
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
!                  GEOS-Chem Global Chemical Transport Model                  !
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
  SUBROUTINE Diagnostics_Write ( am_I_Root, Input_Opt, State_Chm, RESTART, RC ) 
!
! !USES:
!
    USE HCO_INTERFACE_MOD,  ONLY : HcoState, GetHcoID
    USE HCOI_GC_MAIN_MOD,   ONLY : HCOI_GC_WriteDiagn 
    USE HCOIO_Diagn_Mod,    ONLY : HCOIO_Diagn_WriteOut
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! Are we on the root CPU?
    LOGICAL,          INTENT(IN   )  :: RESTART    ! Write restart file? 
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input Options object
    TYPE(ChmState),   INTENT(IN   )  :: State_Chm  ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT  )  :: RC         ! Failure or success
!
! !REVISION HISTORY: 
!  09 Jan 2015 - C. Keller   - Initial version
!  15 Jan 2015 - R. Yantosca - Now accept Input_Opt via the arg list
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)       :: LOC = 'Diagnostics_Write (diagnostics_mod.F90)'

    !=======================================================================
    ! Diagnostics_Write begins here 
    !=======================================================================

    ! Write HEMCO diagnostics
    CALL HCOI_GC_WriteDiagn( am_I_Root, Input_Opt, RESTART, RC )
    IF ( RC /= GC_SUCCESS ) RETURN

    ! Write netCDF GEOS-Chem diagnostics
#if defined( NC_DIAG )

    !-----------------------------------------------------------------------
    ! Eventually write out emission totals to GEOS-Chem logfile
    !-----------------------------------------------------------------------
    CALL TotalsToLogfile( am_I_Root, Input_Opt, State_Chm, RC )
    IF ( RC /= GC_SUCCESS ) THEN 
       CALL ERROR_STOP ('Error in TotalsToLogfile', LOC ) 
    ENDIF

    !-----------------------------------------------------------------------
    ! RESTART: write out all diagnostics. Use current time stamp and save into
    ! restart file.
    ! GEOS-Chem restart file is currently not defined ...
    !-----------------------------------------------------------------------
    IF ( RESTART ) THEN

       ! Do nothing. GEOS-Chem restart is written by call to routine
       ! WRITE_GC_RESTART in main.F. That routine is implemented in
       ! GeosCore/restart_mod.F (ewl, 2/5/16)

       ! Force regular diagnostics to be written out (ckeller, 11/1/16)
       CALL HCOIO_DIAGN_WRITEOUT( am_I_Root,                                & 
                                  HcoState,                                 &
                                  ForceWrite  = .TRUE.,                     &
                                  UsePrevTime = .FALSE.,                    &
                                  COL         = Input_Opt%DIAG_COLLECTION,  &
                                  RC          = RC                         )
       IF ( RC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP( 'Diagnostics write error A', LOC ) 
       ENDIF

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
          CALL ERROR_STOP( 'Diagnostics write error B', LOC ) 
       ENDIF

    ENDIF
#endif

    ! Leave w/ success
    RC = GC_SUCCESS

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
!    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
!    USE HCO_Error_Mod
!    USE Input_Opt_Mod,      ONLY : OptInput
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
!    RC = GC_SUCCESS
!
!    ! Get diagnostic parameters from the Input_Opt object
!    Collection = Input_Opt%DIAG_COLLECTION
!    OutOper    = Input_Opt%ND21_OUTPUT_TYPE
!    
!    ! Loop over # of cloud OD and fractions diagnostics
!    DO N = 1, NumDiagND21
!
!       ! If the species number is scheduled for output in input.geos, 
!       ! then define the diagnostic container for that species.
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
!          CALL Diagn_Create( am_I_Root,  HcoState,          &
!                             Col       = Collection,        & 
!                             cId       = cId,               &
!                             cName     = TRIM( DiagnName ), &
!                             SpaceDim  =  3,                &
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
!    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
!    USE HCO_Error_Mod
!    USE Input_Opt_Mod,      ONLY : OptInput
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
!    RC = GC_SUCCESS
!
!    ! Get diagnostic parameters from the Input_Opt object
!    Collection = Input_Opt%DIAG_COLLECTION
!    OutOper    = Input_Opt%ND22_OUTPUT_TYPE
!    
!    ! Loop over # of photolysis rate diagnostics
!    DO N = 1, NumDiagND22 (NO2, HNO3, H2O2, CH2O, O3, OH source from O3)
!
!       ! If the species number is scheduled for output in input.geos, 
!       ! then define the diagnostic container for that species.
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
!          CALL Diagn_Create( am_I_Root, HcoState,        &
!                          Col       = Collection,        & 
!                          cId       = cId,               &
!                          cName     = TRIM( DiagnName ), &
!                          SpaceDim  =  3,                &
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
!    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
!    USE HCO_Error_Mod
!    USE Input_Opt_Mod,      ONLY : OptInput
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
!    RC = GC_SUCCESS
!
!    ! Get diagnostic parameters from the Input_Opt object
!    Collection = Input_Opt%DIAG_COLLECTION
!    OutOper    = Input_Opt%ND03_OUTPUT_TYPE
!    
!    ! Assign array of species numbers corresponding to decaying species
!    TRCN(1) = IDTPB
!    TRCN(2) = IDTRN
!    TRCN(3) = IDTBE7
!
!    ! Loop over # of mercury diagnostics
!    DO N = 1, NumDiagND03
!
!       ! If the species number is scheduled for output in input.geos, 
!       ! then define the diagnostic container for that species.
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
!          CALL Diagn_Create( am_I_Root,  HcoState,          &
!                             Col       = Collection,        & 
!                             cId       = cId,               &
!                             cName     = TRIM( DiagnName ), &
!                             SpaceDim  =  2,                &
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
!    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
!    USE Input_Opt_Mod,      ONLY : OptInput
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
!    RC = GC_SUCCESS
!
!    ! Get diagnostic parameters from the Input_Opt object
!    Collection = Input_Opt%DIAG_COLLECTION
!    OutOper    = Input_Opt%ND05_OUTPUT_TYPE
!    
!    ! Loop over # of diagnostics (only ones with units kg S)
!    DO N = 1, NumDiagND05
!
!       ! If the species number is scheduled for output in input.geos, 
!       ! then define the diagnostic container for that species.
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
!          CALL Diagn_Create( am_I_Root,  HcoState,          &
!                             Col       = Collection,        & 
!                             cId       = cId,               &
!                             cName     = TRIM( DiagnName ), &
!                             SpaceDim  =  3,                &
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
!          CALL Diagn_Create( am_I_Root,  HcoState,          &
!                             Col       = Collection,        & 
!                             cId       = cId,               &
!                             cName     = TRIM( DiagnName ), &
!                             SpaceDim  =  3,                &
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
!          CALL Diagn_Create( am_I_Root,  HcoState,          &
!                             Col       = Collection,        & 
!                             cId       = cId,               &
!                             cName     = TRIM( DiagnName ), &
!                             SpaceDim  =  3,                &
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
!    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
!    USE HCO_Error_Mod
!    USE Input_Opt_Mod,      ONLY : OptInput
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
!    RC = GC_SUCCESS
!
!    ! Get diagnostic parameters from the Input_Opt object
!    Collection = Input_Opt%DIAG_COLLECTION
!    OutOper    = Input_Opt%ND07_OUTPUT_TYPE
!    
!    ! Loop over # of diagnostics (NOTE: diags 10-23 only if SOA simulation)
!    DO N = 1, NumDiagND07
!
!       ! If the species number is scheduled for output in input.geos, 
!       ! then define the diagnostic container for that species.
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
!          CALL Diagn_Create( am_I_Root,  HcoState,          &
!                             Col       = Collection,        & 
!                             cId       = cId,               &
!                             cName     = TRIM( DiagnName ), &
!                             SpaceDim  =  2,                &
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
!!  diagnostics for hourly maximum species mixing ratio at the surface (aka
!!  ND71).
!!\\
!!\\
!! !INTERFACE:
!!
!  SUBROUTINE DiagnInit_HrlyMax_SurfConc( am_I_Root, Input_Opt, State_Chm, RC )
!!
!! !USES:
!!
!    USE HCO_Diagn_Mod,      ONLY : Diagn_Create
!    USE Input_Opt_Mod,      ONLY : OptInput
!    USE State_Chm_Mod,      ONLY : ChmState
!!
!! !INPUT PARAMETERS:
!!
!    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?!
!    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
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
!    RC = GC_SUCCESS
!
!    ! Get diagnostic parameters from the Input_Opt object
!    Collection = Input_Opt%DIAG_COLLECTION
!    OutOper    = Input_Opt%ND71_OUTPUT_TYPE
!
!    ! TODO Check if certain species listed for ND71 in input.geos 
!
!    !----------------------------------------------------------------
!    ! Create containers for Species Mixing Ratio
!    !----------------------------------------------------------------
!
!    ! Diagnostic name
!    DiagnName = 'IJ-MAX'
!    DO M = 1, Input_Opt%TMAX(71)
!       N = Input_Opt%TINDEX(71,M)
!       IF ( N > State_Chm%nAdvect ) CYCLE
!
!!       ! line from diag3
!!       ARRAY(:,:,1) = AD71(:,:,N)/AD71_COUNT
!
!       ! Create container
!       CALL Diagn_Create( am_I_Root, HcoState,        &
!                       Col       = Collection,        &
!                       cId       = cId,               &
!                       cName     = TRIM( DiagnName ), &
!                       SpaceDim  =  3,                &
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
