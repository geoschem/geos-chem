!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: mixing_mod.F90
!
! !DESCRIPTION: Module mixing\_mod.F90 is a wrapper module for the PBL mixing
! in GEOS-Chem.
!\\
!\\
! !INTERFACE:
!
MODULE Mixing_Mod
!
! !USES:
!
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: DO_MIXING
  PUBLIC :: DO_TEND
  PUBLIC :: Do_Emis_Above_Sfc
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Add_Emis_To_Box
  PRIVATE :: Check_For_Negative_Conc
  PRIVATE :: Get_Emis_Top
  PRIVATE :: Init_Tend_Species_Ids
  PRIVATE :: Is_ChemGrid_Only
  PRIVATE :: Tend_Budget
  PRIVATE :: Tend_Convert_Units
!
! !PRIVATE TYPES:
!
  LOGICAL :: ids_are_set = .FALSE.
  INTEGER :: id_ACET,  id_ALD2,  id_ALK4,  id_Br,    id_Br2
  INTEGER :: id_BrNO3, id_BrO,   id_C2H6,  id_C3H8,  id_CH2O
  INTEGER :: id_CO2,   id_HBr,   id_HNO3,  id_HOBr,  id_MACR
  INTEGER :: id_O3,    id_PRPE,  id_RCHO
!
! !REVISION HISTORY:
!  04 Mar 2015 - C. Keller   - Initial version.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_mixing
!
! !DESCRIPTION: Subroutine DO\_MIXING performs the PBL mixing.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Mixing( Input_Opt,  State_Chm, State_Diag, &
                        State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE Pbl_Mix_Mod,    ONLY : Do_Full_Pbl_Mixing
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE Vdiff_Mod,      ONLY : Do_Vdiff
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt   ! Input Options
    TYPE(GrdState),   INTENT(IN   )  :: State_Grid  ! Grid State
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState),   INTENT(INOUT)  :: State_Met   ! Meteorology State
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm   ! Chemistry State
    TYPE(DgnState),   INTENT(INOUT)  :: State_Diag  ! Diagnostics State
    INTEGER,          INTENT(INOUT)  :: RC          ! Failure or success
!
! !REMARKS
!  While all dry deposition rates are calculated either in DO_PBL_MIX2 or
!  DO_TEND, settling of aerosols is still computed in the dust & seasalt
!  modules.  The order of operations for full PBL mixing (TURBDAY) or no
!  PBL mixing stay the same, so TURBDAY and no-mixing stay zero-diff.
!
! !REVISION HISTORY:
!  04 Mar 2015 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! DO_MIXING begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at DO_MIXING (in module GeosCore/mixing_mod.F90)'

    !========================================================================
    ! Non-local PBL mixing (VDIFF)
    !========================================================================
    IF ( Input_Opt%LTURB .and. Input_Opt%LNLPBL ) THEN

       ! Add emissions in levels 2 and higher directly to the species array.
       ! NOTE: Do_Emis_Above_Sfc is only called for non-local PBL mixing.
       CALL Do_Emis_Above_Sfc( Input_Opt,  State_Chm, State_Diag,            &
                               State_Grid, State_Met, RC                    )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Do_Emis_Above_Sfc"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Call VDIFF to perform the non-local PBL mixing.  The non-local term
       ! (surface emissions and dry deposition flux) enters VDIFF via the
       ! State_Chm%SurfaceFlux array (defined in Compute_Sflx_for_Vdiff,
       ! located in module GeosCore/hco_interface_gc_mod.F90).
       CALL Do_Vdiff( Input_Opt,  State_Chm, State_Diag,                     &
                      State_Grid, State_Met, RC                             )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Do_Vdiff"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Exit successfully
       RETURN
    ENDIF

    !========================================================================
    ! Full PBL mixing (TURBDAY), or no PBL mixing
    !========================================================================

    ! DO_TEND applies all emissions (level 1 and up) and dry deposition
    ! to the species array.  It is not called when VDIFF is used; in
    ! that case emissions in levels 2 and higher are applied by
    ! Do_Emis_Above_Sfc, and level-1 emissions and drydep are passed
    ! to VDIFF via State_Chm%SurfaceFlux.  (GCHP also calls DO_TEND
    ! directly when turbulence is handled outside of GEOS-Chem.)
    CALL Do_Tend( Input_Opt,  State_Chm, State_Diag,                         &
                  State_Grid, State_Met, RC                                 )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "DO_TEND"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! If PBL mixing is on, fully mix the PBL (TURBDAY)
    IF ( Input_Opt%LTURB ) THEN
       CALL Do_Full_Pbl_Mixing( Input_Opt,  State_Chm, State_Diag,           &
                                State_Grid, State_Met, RC                   )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Do_Full_Pbl_Mixing"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE DO_MIXING
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_tend
!
! !DESCRIPTION: Subroutine DO\_TEND adds the species tendencies (dry deposition
!  and emissions) to the species array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_TEND( Input_Opt, State_Chm, State_Diag, State_Grid,          &
                      State_Met, RC,        DT                              )
!
! !USES:
!
    USE Diagnostics_Mod,      ONLY : Compute_Budget_Diagnostics
    USE ErrCode_Mod
    USE ERROR_MOD,            ONLY : SAFE_DIV
    USE GET_NDEP_MOD,         ONLY : SOIL_DRYDEP
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_GetDiagn
    USE HCO_Utilities_GC_Mod, ONLY : GetHcoValEmis, GetHcoValDep, InquireHco
    USE HCO_Utilities_GC_Mod, ONLY : LoadHcoValEmis, LoadHcoValDep
    USE Input_Opt_Mod,        ONLY : OptInput
    USE PhysConstants,        ONLY : AVO
    USE Species_Mod,          ONLY : Species
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Chm_Mod,        ONLY : Ind_
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE TIME_MOD,             ONLY : GET_TS_DYN, GET_TS_CONV, GET_TS_CHEM
    USE Timers_Mod,           ONLY : Timer_End, Timer_Start
    USE UnitConv_Mod
#ifdef MODEL_CLASSIC
    use hco_utilities_gc_mod, only: TMP_MDL ! danger
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN   )           :: Input_Opt    ! Input opts
    TYPE(MetState),   INTENT(IN   )           :: State_Met    ! Met state
    TYPE(GrdState),   INTENT(IN   )           :: State_Grid   ! Grid state
    REAL(fp),         INTENT(IN   ), OPTIONAL :: DT           ! Time step [s]
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT)           :: State_Chm    ! Chemistry state
    TYPE(DgnState),   INTENT(INOUT)           :: State_Diag   ! Diags State
    INTEGER,          INTENT(INOUT)           :: RC           ! Success/Failure
!
! !REVISION HISTORY:
!  04 Mar 2015 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                 :: found,      chem_grid_only
    LOGICAL                 :: dryd_spec,  emis_spec,      error
    INTEGER                 :: I,          J,              L
    INTEGER                 :: err_status, max_level,      N
    INTEGER                 :: D,          NN,             NA
    INTEGER                 :: S,          dryDepId,       previous_units
    INTEGER                 :: pbl_top,    dryd_top,       emis_top
    REAL(fp)                :: TS,         air_sea_freq,   freq,
    REAL(fp)                :: RKT,        frac,           flux
    REAL(fp)                :: AREA_M2,    MWkg,           denom

    ! PARANOX loss fluxes (kg/m2/s). These are obtained from the
    ! HEMCO PARANOX extension via the diagnostics module.
    REAL(fp)                :: PNOxLoss
    REAL(f4),      POINTER  :: PNOxLoss_O3  (:,:)
    REAL(f4),      POINTER  :: PNOxLoss_HNO3(:,:)
    REAL(f4),      POINTER  :: Ptr2D        (:,:)

    ! Pointers and objects
    TYPE(Species), POINTER  :: SpcInfo
    REAL(fp),      POINTER  :: DepFreq(:,:,:)  ! IM, JM, nDryDep

    ! Strings
    CHARACTER(LEN=255)      :: thisLoc
    CHARACTER(LEN=512)      :: errMsg, errorMsg

#ifdef ADJOINT
    LOGICAL                 :: IS_ADJ
#endif

    !========================================================================
    ! DO_TEND begins here!
    !========================================================================

    ! Assume success
    RC            =  GC_SUCCESS
    errMsg        =  ''
    thisLoc       =  ' -> at DO_TEND (in module GeosCore/mixing_mod.F90)'

    ! Special case that there is no dry deposition and emissions
    IF ( .NOT. Input_Opt%LDRYD .AND. .NOT. Input_Opt%DoEmissions ) RETURN

    ! Initialize
    error         =  .FALSE.
    DepFreq       => State_Chm%DryDepFreq
    SpcInfo       => NULL()
    PNOxLoss_O3   => NULL()
    PNOxLoss_HNO3 => NULL()
    Ptr2D         => NULL()

    ! Look up species indices (first call only)
    CALL Init_Tend_Species_Ids()

    !------------------------------------------------------------------------
    ! Emissions/dry deposition budget diagnostics - Part 1 of 2
    !------------------------------------------------------------------------
    CALL Tend_Budget(                                                        &
         Input_Opt = Input_Opt,                                              &
         State_Chm = State_Chm,
        State_Diag,                     &
                      State_Grid, State_Met, .TRUE.,     RC                 )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Tend_Budget"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Convert species to [kg/m2] (removes area dependency, ewl 9/30/15)
    !------------------------------------------------------------------------
    CALL Tend_Convert_Units( Input_Opt,  State_Chm, State_Grid,              &
                             State_Met,  KG_SPECIES_PER_M2, RC,              &
                             previous_units                                 )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Tend_Convert_Units"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Get time step [s]
    IF ( PRESENT( dt ) ) THEN
       time_step = dt
    ELSE
       time_step = Get_Ts_Dyn()
    ENDIF
#ifdef ADJOINT
    IF ( Input_Opt%Is_Adjoint ) THEN
       time_step = time_step * -1
    ENDIF
#endif

    !=======================================================================
    ! Get pointers to the PARANOX loss fluxes. These are
    ! stored in diagnostics 'PARANOX_O3_DEPOSITION_FLUX' and
    ! 'PARANOX_HNO3_DEPOSITION_FLUX'. The call below links pointers
    ! PNOxLoss_O3 and PNOxLoss_HNO3 to the data values stored in the
    ! respective diagnostics. The pointers will remain unassociated if
    ! the diagnostics do not exist.
    ! This is only needed if non-local PBL scheme is not being used.
    ! Otherwise, PARANOX fluxes are applied in vdiff_mod.F90.
    !  (ckeller, 4/10/2015)
    !
    ! If using HEMCO Intermediate grid feature, then the call needs to be
    ! refreshed at every time step for regridding. (hplin, 6/21/20)
    !=======================================================================
    IF ( .NOT. Input_Opt%LNLPBL ) THEN

       ! O3 flux from PARANOX
       CALL HCO_GC_GetDiagn(                                                 &
            Input_Opt      = Input_Opt,                                      &
            State_Grid     = State_Grid,                                     &
            DiagnName      = 'PARANOX_O3_DEPOSITION_FLUX',                   &
            StopIfNotFound = .FALSE.,                                        &
            Ptr2D          = Ptr2D,                                          &
            RC             = RC                                             )
       IF ( ASSOCIATED( Ptr2D ) ) THEN
          ALLOCATE( PNOxLoss_O3( State_Grid%NX, State_Grid%NY ), STAT=RC )
          PNOxLoss_O3 = Ptr2D
       ENDIF
       Ptr2D => NULL()

       ! HNO3 flux from PARANOX
       CALL HCO_GC_GetDiagn(                                                 &
            Input_Opt      = Input_Opt,                                      &
            State_Grid     = State_Grid,                                     &
            DiagnName      = 'PARANOX_HNO3_DEPOSITION_FLUX',                 &
            StopIfNotFound = .FALSE.,                                        &
            Ptr2D          = Ptr2D,                                          &
            RC             = RC                                             )
       IF( ASSOCIATED( Ptr2D )) THEN
          ALLOCATE( PNOxLoss_HNO3( State_Grid%NX, State_Grid%NY ), STAT=RC )
          PNOxLoss_HNO3(:,:) = Ptr2D(:,:)
       ENDIF
       Ptr2D => NULL()
    ENDIF

    !========================================================================
    ! Do for every advected species and grid box
    !
    ! NOTE: For GEOS-Chem Classic HEMCO "Intermediate Grid" feature, where
    ! HEMCO is running on a distinct grid from the model, the on-demand
    ! regridding is most optimized when contiguous accesses to GetHcoValEmis
    ! and GetHcoValDep are performed for the given species. Therefore,
    ! it is most optimal to call in the following fashion (IJKN)
    !
    !  Emis(1,1,1,1) -> Emis(1,1,2,1) -> ... -> Dep(1,1,1,1) -> Dep(1,1,2,1)
    !
    ! By switching emis/dep or the species # LAST, as either of these
    ! changing WILL trigger a new regrid and thrashing of the old buffer.
    !
    ! Therefore, the loop below has been adjusted to run serially for each
    ! species, and parallelizing the inner I, J loop instead.  Also, moved
    ! some non-I,J specific variables outside of the loop for optimization
    ! (Haipeng Lin 6/27/20).
    !========================================================================
    DO NA = 1, State_Chm%nAdvect

       ! Get the species ID from the advected species ID
       N = State_Chm%Map_Advect(NA)

       ! Get info about this species from the species database
       SpcInfo => State_Chm%SpcData(N)%Info

       ! Molecular weight in kg
       MWkg = SpcInfo%MW_g * 1.e-3_fp

       !--------------------------------------------------------------------
       ! Check if we need to do dry deposition for this species
       !--------------------------------------------------------------------

       ! Initialize
       is_hco_drydep = .FALSE.
       id_drydep     = -1

       ! Only if dry deposition is turned on and we do want to consider
       ! processes below the PBL...
       IF ( Input_Opt%LDRYD ) THEN

          ! Get dry deposition ID (used by drydep_mod.F90) for this species.
          id_drydep = SpcInfo%DryDepId

          ! Check if this is a HEMCO drydep species
          is_hco_drydep = ( id_drydep > 0 )
          IF ( .NOT. is_hco_drydep ) THEN
             CALL InquireHco ( N, Dep=is_hco_drydep )
          ENDIF

          ! Special case for O3 or HNO3: include PARANOX loss
          IF ( N == id_O3   .AND. ASSOCIATED( PNOxLoss_O3   ) ) THEN
             is_hco_drydep = .TRUE.
          ENDIF
          IF ( N == id_HNO3 .AND. ASSOCIATED( PNOxLoss_HNO3 ) ) THEN
             is_hco_drydep = .TRUE.
          ENDIF
       ENDIF

       ! Set emissions top level:
       ! This is the top of atmosphere unless concentration build-up
       ! in stratosphere wants to be avoided.
       chem_grid_only  = .FALSE.

       ! Set emissions to zero above chemistry grid for the following VOCs
       IF ( N == id_MACR .OR. N == id_RCHO .OR. &
            N == id_ACET .OR. N == id_ALD2 .OR. &
            N == id_ALK4 .OR. N == id_C2H6 .OR. &
            N == id_C3H8 .OR. N == id_CH2O .OR. &
            N == id_PRPE                         ) THEN
          ChemGridOnly = .TRUE.
       ENDIF

       ! Bry concentrations become prescribed in lin. strat. chemistry.
       ! Therefore avoid any emissions of these compounds above the
       ! chemistry grid (lin. strat. chem. applies above chemistry grid
       ! only).
       IF ( Input_Opt%LINEAR_CHEM ) THEN
          IF ( N == id_BrO  .OR. N == id_Br2   .OR. &
               N == id_Br   .OR. N == id_HOBr  .OR. &
               N == id_HBr  .OR. N == id_BrNO3       ) THEN
             ChemGridOnly = .TRUE.
          ENDIF
       ENDIF

       !--------------------------------------------------------------------
       ! Check if we need to do emissions for this species
       !--------------------------------------------------------------------
       IF ( Input_Opt%DoEmissions ) THEN
          CALL InquireHco ( N, Emis=EmisSpec )
       ELSE
          EmisSpec = .FALSE.
       ENDIF

       ! If there is emissions for this species, it must be loaded into
       ! memory first.  This is achieved by attempting to retrieve a grid
       ! box while NOT in a parallel loop. Failure to load this will result
       ! in severe performance issues!! (Haipeng Lin, 9/27/20)
       IF ( EmisSpec ) THEN
          CALL LoadHcoValEmis( Input_Opt, State_Grid, N )
       ENDIF

       IF ( DryDepSpec ) THEN
          CALL LoadHcoValDep( Input_Opt, State_Grid, N )
       ENDIF

       !--------------------------------------------------------------------
       ! Can go to next species if this species does not have
       ! dry deposition and/or emissions
       !--------------------------------------------------------------------
       IF ( .NOT. DryDepSpec .AND. .NOT. EmisSpec ) CYCLE

       ! Loop over all grid boxes
       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( I,        J,          L,            max_level         )&
       !$OMP PRIVATE( pbl_top,  found       air_sea_freq, freq              )&
       !$OMP PRIVATE( RKT,      frac,       flux,         area_m2           )&
       !$OMP PRIVATE( dryd_top, emis_top,   PNOxLoss,     denom,            )&
       !$OMP PRIVATE( S,        err_status                                  )&
       !$OMP COLLAPSE( 2                                                    )&
       !$OMP SCHEDULE( GUIDED                                               )&
       !$OMP REDUCTION( .OR. : error                                        )
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Continue to end of loop if an error has occured.  Then we can
          ! exit outside the parallel loop.  This is thread-safe.
          IF ( error ) CYCLE

          !-----------------------------------------------------------------
          ! Define various quantities before computing tendencies
          !-----------------------------------------------------------------

          ! Get the PBL top level at this grid box
          pbl_top = MAX( 1, FLOOR( State_Met%PBL_TOP_L(I,J) ) )

          ! Top level for dry deposition
          IF ( Input_Opt%PBL_DRYDEP ) THEN
             dryd_top = pbl_top
          ELSE
             dryd_top = 1
          ENDIF

          ! Restrict emissions of some species to the chemistry grid
          chem_grid_only = Is_ChemGrid_Only( Input_Opt, N )

          ! Top level for emissions
          emis_top = Get_Emis_Top( State_Grid, State_Met, chem_grid_only, I, J )

          ! Loop over surface levels (to the top level for
          ! dry deposition or emissions, whichever is higher)
          DO L = 1, MAX( dryd_top, emis_top )

             !--------------------------------------------------------------
             ! Apply dry deposition frequencies to all levels below the
             ! PBL top.
             !--------------------------------------------------------------
             IF ( dryd_spec .AND. ( L <= dryd_top ) ) THEN

                ! Init
                freq = 0.0_fp

                ! Get the drydep frequency [s-1] from drydep_mod.F90.
                IF ( id_drydep > 0 ) THEN
                   freq = DepFreq(I,J,id_drydep)
                ENDIF

                ! Then add the air-to-sea deposition frequency [s-1]
                ! from the HEMCO SeaFlux extension if it exists.
                CALL GetHcoValDep(                                           &
                     I          = I,                                         &
                     J          = J,                                         &
                     L          = 1,                                         &
                     N          = N,                                         &
                     Input_Opt  = Input_Opt,                                 &
                     State_Grid = State_Grid,                                &
                     found      = found,                                     &
                     dep        = air_sea_freq                              )
                IF ( found ) freq = freq + air_sea_freq

                ! Get PARANOX deposition loss [kg/m2/s], at surface only
                PNOxLoss = 0.0_fp
                IF ( L == 1 ) THEN
                   IF ( N == id_O3   .AND. ASSOCIATED( PNOxLoss_O3  ) ) THEN
                      PNOxLoss = PNOxLoss_O3(I,J)
                   ENDIF
                   IF ( N == id_HNO3 .AND. ASSOCIATED( PNOxLoss_HNO3) ) THEN
                      PNOxLoss = PNOxLoss_HNO3(I,J)
                   ENDIF
                ENDIF

                ! Apply dry deposition
                IF ( freq > 0.0_fp .OR. PNOxLoss > 0.0_fp ) THEN

                   ! Compute exponential loss term
                   RKT  = freq * ts
                   frac = EXP( -RKT )

                   ! Loss in kg/m2
                   flux = ( 1.0_fp - frac ) * State_Chm%Species(N)%Conc(I,J,L)

                   ! Apply dry deposition
                   State_Chm%Species(N)%Conc(I,J,L) = frac *                 &
                      State_Chm%Species(N)%Conc(I,J,L)

#ifdef ADJOINT
                   if (Input_Opt%Is_Adjoint) then
                      State_Chm%SpeciesAdj(I,J,L,N) = frac *                 &
                         State_Chm%SpeciesAdj(I,J,L,N)
                   endif
#endif

                   ! Apply PARANOX loss [kg/m2/s]
                   IF ( PNOxLoss > 0 ) THEN
                      State_Chm%Species(N)%Conc(I,J,L) =                     &
                      State_Chm%Species(N)%Conc(I,J,L) - ( PNOxLoss * TS )
                      flux = flux + ( PNOxLoss * TS )
                   ENDIF

                   ! Loss in [molec/cm2/s]
                   denom = ( MWkg * ts * 1.0e+4_fp ) / AVO
                   flux  = SAFE_DIV( flux, denom, 0.0_fp )  ! molec/cm2/s

                   ! Eventually add to SOIL_DRYDEP
                   IF ( Input_Opt%LSOILNOX ) THEN
                      CALL SOIL_DRYDEP( I, J, N, flux, State_Chm )
                   ENDIF

                   !--------------------------------------------------------
                   ! HISTORY: Archive drydep flux loss from mixing
                   ! Units = molec/cm2/s
                   !
                   ! NOTE: we don't need to multiply by the ratio of
                   ! TS_CONV / TS_CHEM, as the updating frequency for
                   ! HISTORY is determined by the "frequency" setting in
                   ! the "HISTORY.rc" input file.
                   !--------------------------------------------------------
                   IF ( ( State_Diag%Archive_DryDepMix .or.                  &
                          State_Diag%Archive_DryDep        )   .and.         &
                          DryDepID > 0                       ) THEN
                      S = State_Diag%Map_DryDepMix%id2slot(DryDepID)
                      IF ( S > 0 ) THEN
                         State_Diag%DryDepMix(I,J,S) = flux
                      ENDIF
                   ENDIF

                ENDIF ! apply drydep
             ENDIF ! L <= PBLTOP

             !--------------------------------------------------------------
             ! Apply emissions from HEMCO [kg/m2/s] x TS
             !--------------------------------------------------------------
             IF ( emisSpec .AND. ( L <= emis_top ) ) THEN
                CALL Add_Emis_To_Box( Input_Opt, State_Grid, State_Chm,      &
                                      N,         I,          J,              &
                                      L,         ts                         )
             ENDIF

             ! Exit this loop if negatives are found.  We will trap the 
             ! error outside of the parallel loop in a thread-safe manner.
             CALL Check_For_Negative_Conc(                                   &
                  I           = I,                                           &
                  J           = J,                                           &
                  L           = L,                                           &
                  N           = N,                                           &
                  State_Chm   = State_Chm,                                   &
                  is_negative = is_negative                                 )
             IF ( is_negative ) THEN
                error      = .TRUE.
                err_status = 1
                CYCLE
             ENDIF

          ENDDO !L
       ENDDO !J
       ENDDO !I
!$OMP END PARALLEL DO

       !====================================================================
       ! Exit if any routine in the loop encountered an error
       ! This is a thread-safe implementation.
       !====================================================================
       IF ( error ) THEN

          ! Construct error message
          SELECT CASE( err_status )
             CASE( 1 )
                errMsg = 'Error encountered in ""!'
          END SELECT
          CALL GC_Error( errMsg, RC, thisLoc )

          ! Free pointers and return
          IF ( ASSOCIATED( PNOxLoss_O3   ) ) DEALLOCATE( PNOxLoss_O3   )
          IF ( ASSOCIATED( PNOxLoss_HNO3 ) ) DEALLOCATE( PNOxLoss_HNO3 )
          PNOxLoss_O3   => NULL()
          PNOxLoss_HNO3 => NULL()
          SpcInfo       => NULL()
          RETURN

       ENDIF


    ENDDO !N

#if defined( ADJOINT )  && defined ( DEBUG )
    IF (Input_Opt%is_adjoint .and. Input_Opt%IS_FD_SPOT_THIS_PET) THEN
       WRITE(*,*) ' SpcAdj(IFD,JFD) before unit converstion: ',  &
            State_Chm%SpeciesAdj(Input_Opt%IFD, Input_Opt%JFD, &
            Input_Opt%LFD, Input_Opt%NFD)
       WRITE(*,*) ' Spc(IFD,JFD) before unit converstion: ',  &
            State_Chm%Species(Input_Opt%NFD)%Conc(Input_Opt%IFD, Input_Opt%JFD, Input_Opt%LFD)
    ENDIF

#endif

    ! Convert species back to the original units
    CALL Tend_Convert_Units( Input_Opt, State_Chm,      State_Grid,          &
                             State_Met, previous_units, RC                  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Tend_Convert_Units" #2!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       IF ( ASSOCIATED( PNOxLoss_O3   ) ) DEALLOCATE( PNOxLoss_O3   )
       IF ( ASSOCIATED( PNOxLoss_HNO3 ) ) DEALLOCATE( PNOxLoss_HNO3 )
       PNOxLoss_O3   => NULL()
       PNOxLoss_HNO3 => NULL()
       SpcInfo       => NULL()
       RETURN
    ENDIF

    ! Start mixing timer again
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_Start( "Boundary layer mixing", RC )
    ENDIF

#if defined( ADJOINT )  && defined ( DEBUG )
    IF (Input_Opt%is_adjoint .and. Input_Opt%IS_FD_SPOT_THIS_PET) THEN
       WRITE(*,*) ' SpcAdj(IFD,JFD) after unit converstion: ',  &
            State_Chm%SpeciesAdj(Input_Opt%IFD, Input_Opt%JFD, &
            Input_Opt%LFD, Input_Opt%NFD)
       WRITE(*,*) ' Spc(IFD,JFD) after unit converstion: ',  &
            State_Chm%Species(Input_Opt%NFD)%Conc(Input_Opt%IFD, Input_Opt%JFD, Input_Opt%LFD)
    ENDIF

#endif

    !========================================================================
    ! Emissions/dry deposition budget diagnostics - Part 2 of 2
    !========================================================================
    IF ( State_Diag%Archive_BudgetEmisDryDep ) THEN

       ! Timestep for diagnostics [s]
       DT_Tend = DBLE( TS )

       CALL Tend_Budget( Input_Opt,  State_Chm, State_Diag,                  &
                         State_Grid, State_Met, .FALSE.,   RC,               &
                         timeStep=DT_Tend                                   )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Tend_Budget" #2!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          IF ( ASSOCIATED( PNOxLoss_O3   ) ) DEALLOCATE( PNOxLoss_O3   )
          IF ( ASSOCIATED( PNOxLoss_HNO3 ) ) DEALLOCATE( PNOxLoss_HNO3 )
          PNOxLoss_O3   => NULL()
          PNOxLoss_HNO3 => NULL()
          SpcInfo       => NULL()
          RETURN
       ENDIF
    ENDIF

    !========================================================================
    ! Cleanup and quit
    !========================================================================

    ! Deallocate pointers
    IF ( ASSOCIATED( PNOxLoss_O3   ) ) DEALLOCATE( PNOxLoss_O3   )
    IF ( ASSOCIATED( PNOxLoss_HNO3 ) ) DEALLOCATE( PNOxLoss_HNO3 )

    ! Free pointers
    DepFreq       => NULL()
    PNOxLoss_O3   => NULL()
    PNOxLoss_HNO3 => NULL()
    SpcInfo       => NULL()

  END SUBROUTINE Do_Tend
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_emis_above_sfc
!
! !DESCRIPTION: Adds HEMCO emissions in levels 2 and higher to the species
!  array.  This is called only for non-local PBL mixing (VDIFF), before
!  VDIFF.  Level-1 emissions and dry deposition are passed to VDIFF via
!  State\_Chm\%SurfaceFlux (see Compute\_Sflx\_For\_Vdiff).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Emis_Above_Sfc( Input_Opt,  State_Chm, State_Diag,           &
                                State_Grid, State_Met, RC                   )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Utilities_GC_Mod, ONLY : InquireHco, LoadHcoValEmis
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE Time_Mod,             ONLY : Get_Ts_Dyn
    USE UnitConv_Mod,         ONLY : KG_SPECIES_PER_M2
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid state
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Met state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Emissions above the surface are counted in the BudgetEmisDryDep*
!  diagnostics, just as DO_TEND does for full PBL mixing.
!
! !REVISION HISTORY:
!  16 Oct 2026 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: chem_grid_only, has_emissions
    INTEGER            :: emis_top,       I
    INTEGER            :: J,              L
    INTEGER            :: L,              N
    INTEGER            :: NA,             previous_units
    REAL(fp)           :: time_step
    REAL(f8)           :: dt_tend

    ! Strings
    CHARACTER(LEN=255) :: thisLoc
    CHARACTER(LEN=512) :: errMsg

    !========================================================================
    ! Do_Emis_Above_Sfc begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    error   = .FALSE.
    errMsg  = ''
    thisLoc = &
       ' -> at Do_Emis_Above_Sfc (in module GeosCore/mixing_mod.F90)'

    ! Exit if emissions are turned off
    IF ( .not. Input_Opt%DoEmissions ) RETURN

    ! Look up species indices (first call only)
    CALL Init_Tend_Species_Ids()

    ! Emissions time step [s] (same as the VDIFF time step)
    time_step = Get_Ts_Dyn()
#ifdef ADJOINT
    IF ( Input_Opt%Is_Adjoint ) time_step = -time_step
#endif

    !========================================================================
    ! Emissions/dry deposition budget diagnostics - Part 1 of 2
    !========================================================================
    CALL Tend_Budget(                                                        &
         Input_Opt   = Input_Opt,                                            &
         State_Chm   = State_Chm,                                            &
         State_Diag  = State_Diag,                                           &
         State_Grid  = State_Grid,                                           &
         State_Met   = State_Met,                                            &
         before      = .TRUE.                                                &
         time_step   = DBLE( time_step ),                                    &
         RC          = RC                                                   )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Tend_Budget" (1 of 2)!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Convert species to [kg/m2]
    !========================================================================
    CALL Tend_Convert_Units(                                                 &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         new_units      = KG_SPECIES_PER_M2,                                 &
         previous_units = previous_units,                                    &
         RC             = RC                                                )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Tend_Convert_Units" (to kg/m2)!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Add emissions in levels 2 and higher to the species array.
    !
    ! Keep one species in the HEMCO regrid buffer at a time (loop over
    ! species outermost), as required for the HEMCO intermediate grid.
    ! errMsg and RC are SHARED: every thread that finds a negative
    ! concentration writes the same values.
    !========================================================================
    DO NA = 1, State_Chm%nAdvect

       ! Skip species without HEMCO emissions
       N = State_Chm%Map_Advect(NA)
       CALL InquireHco( N, emis=has_emissions )
       IF ( .not. has_emissions ) CYCLE

       ! Load emissions into memory outside of the parallel loop
       ! (Haipeng Lin, 9/27/20)
       CALL LoadHcoValEmis( Input_Opt, State_Grid, N )

       ! Determine if emissions of the current species is restricted
       ! to the chemistry grid (or can be applied to the entire grid).
       ! This may be needed to prevent unphysical species concentrations
       ! near the top of the atmosphere.
       chem_grid_only = Is_ChemGrid_Only( Input_Opt, N )

       !---------------------------------------------------------------------
       ! Loop over surface grid boxes
       !---------------------------------------------------------------------
       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( I, J, L, emis_top, is_negative                        )&
       !$OMP COLLAPSE( 2                                                    )&
       !$OMP SCHEDULE( STATIC                                               )&
       !$OMP REDUCTION( .OR. : error                                        )
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Continue to the end of this loop if an error is encountered
          ! on any thread.  Then we'll trap the error outside of the
          ! parallel loop, as this is thread-safe.
          IF ( error ) CYCLE

          ! Find the top level for emissions
          CALL Get_Emis_Top(                                                 &
               I              = I,                                           &
               J              = J,                                           &
               chem_grid_only = chem_grid_only,                              &
               State_Grid     = State_Grid,                                  &
               State_Met      = State_Met,                                   &
               emis_top       = emis_top                                    )

          ! Add emissions in the column to the species array
          DO L = 2, emis_top
             CALL Add_Emis_To_Box(                                           &
                  I           = I,                                           &
                  J           = J,                                           &
                  L           = L,                                           &
                  N           = N,                                           &
                  time_step   = time_step,                                   &
                  Input_Opt   = Input_Opt,                                   &
                  State_Grid  = State_Grid,                                  &
                  State_Chm   = State_Chm                                   )

             ! Check for negative concentrations
             CALL Check_For_Negative_Conc(                                   &
                  I           = I,                                           &
                  J           = J,                                           &
                  L           = L,                                           &
                  N           = N,                                           &
                  State_Chm   = State_Chm,                                   &
                  is_negative = is_negative                                 )
             IF ( is_negative ) THEN
                error = .TRUE.
             ENDIF

          ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       !---------------------------------------------------------------------
       ! Throw an error if any threads returned with error status in the 
       ! loop above.  This is the thread-safe implementation.
       !---------------------------------------------------------------------
       IF ( is_negative ) THEN
          errMsg = 'Negative species concentations encountered! '         // &
                   'This may be fixed by increasing the background '      // &
                   'concentration or by shortening the transport '        // &
                   'time step.'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDDO

    !========================================================================
    ! Convert species back to the original units
    !========================================================================
    CALL Tend_Convert_Units(                                                 &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         previous_units = previous_units,                                    &
         RC             = RC                                                )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Tend_Convert_Units" (from kg/m2)!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Emissions/dry deposition budget diagnostics - Part 2 of 2
    !========================================================================
    CALL Tend_Budget(                                                        &
         Input_Opt   = Input_Opt,                                            &
         State_Chm   = State_Chm,                                            &
         State_Diag  = State_Diag,                                           &
         State_Grid  = State_Grid,                                           &
         State_Met   = State_Met,                                            &
         before      = .FALSE.                                               &
         time_step   = DBLE( time_step ),                                    &
         RC          = RC                                                   )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Tend_Budget" (2 of 2)!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Do_Emis_Above_Sfc
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: add_emis_to_box
!
! !DESCRIPTION: Adds the HEMCO emission flux of species N at grid box
!  (I,J,L), multiplied by the time step, to the species array.  Species
!  must be in [kg/m2].  Safe to call from within a parallel loop.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Add_Emis_To_Box( Input_Opt, State_Grid, State_Chm,              &
                              N,         I,          J,                      &
                              L,         time_step                          )
!
! !USES:
!
    USE HCO_Utilities_GC_Mod, ONLY : GetHcoValEmis
#ifdef MODEL_CLASSIC
    USE HCO_Utilities_GC_Mod, ONLY : TMP_MDL
#endif
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Grid_Mod,       ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid state
    INTEGER,        INTENT(IN)    :: N           ! Species index
    INTEGER,        INTENT(IN)    :: I, J, L     ! Grid box indices
    REAL(fp),       INTENT(IN)    :: time_step   ! Time step [s]
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state
!
! !REVISION HISTORY:
!  16 Oct 2026 - R. Yantosca - Initial version, from DO_TEND
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL  :: found
    REAL(fp) :: emis, flux

    !========================================================================
    ! Add_Emis_To_Box begins here!
    !========================================================================

    ! Initialize
    emis = 0.0_fp
    flux = 0.0_fp

    ! Get HEMCO emissions [kg/m2/s].  For the intermediate grid, read the
    ! regrid buffer directly; this assumes adjacent calls are for the same
    ! species (hplin, 10/10/20).
#ifdef MODEL_CLASSIC
    IF ( Input_Opt%LIMGRID ) THEN
       found = .TRUE.
       emis  = TMP_MDL(I,J,L)
    ELSE
#endif
       CALL GetHcoValEmis(                                                   &
            I          = I,                                                  &
            J          = J,                                                  &
            L          = L,                                                  &
            N          = N,                                                  &
            Input_Opt  = Input_Opt,                                          &
            State_Grid = State_Grid,                                         &
            SkipCheck  =.TRUE.,                                              &
            found      = found,                                              &
            emis       = emis                                               )
#ifdef MODEL_CLASSIC
    ENDIF
#endif

    ! Negative fluxes are allowed (ckeller, 4/12/17)
    IF ( .not. found ) RETURN

    ! Flux: [kg/m2] = [kg/m2/s] x [s]
    flux = emis * time_step

#ifdef ADJOINT
    IF ( I == Input_Opt%IFD .and. J == Input_Opt%JFD .and.                   &
         L == Input_Opt%LFD .and. N == Input_Opt%NFD ) THEN
       WRITE( *, * ) ' GetHcoVal(IFD,JFD) = ', emis, ' FLUX = ', flux
       IF ( Input_Opt%Is_Adjoint ) THEN
          WRITE( *, * ) ' SpeciesAdj(FD) = ', State_Chm%SpeciesAdj(I,J,L,N)
       ENDIF
    ENDIF
#endif

    ! Add to species concentration array
    State_Chm%Species(N)%Conc(I,J,L) = State_Chm%Species(N)%Conc(I,J,L) + flux

  END SUBROUTINE Add_Emis_To_Box
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_negative_conc
!
! !DESCRIPTION: Flags (or, for TOMAS, resets) a negative concentration of
!  species N at grid box (I,J,L).  CO2 is exempt.  Safe to call from
!  within a parallel loop.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Check_For_Negative_Conc( State_Chm, I, J, L, N, is_negative )
!
! !USES:
!
    USE ErrCode_Mod
    USE State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L     ! Grid box indices
    INTEGER,        INTENT(IN)    :: N           ! Species index
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry state
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,        INTENT(OUT)   :: isNegative  ! Was a negative value found?
!   
! !REMARKS:
!  Thread-safe: writes only Conc(I,J,L) and the private flag isNegative.
!  The caller sets RC and the error message after the parallel loop.
!
! !REVISION HISTORY:
!  16 Oct 2026 - R. Yantosca - Initial version, from DO_TEND
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !========================================================================
    ! Check_For_Negative_Conc begins here!
    !========================================================================
    
    ! Initialize
    isNegative = .FALSE.
    
    ! Return if the species concentration is not negative
    IF ( .not. ( State_Chm%Species(N)%Conc(I,J,L) < 0.0_fp ) ) RETURN

#ifdef TOMAS

    ! For TOMAS simulations only, reset negatives to a 
    ! small positive value while we look for the root cause.
    !  -- Betty Croft, Bob Yantosca (21 Jan 2022)
    !$OMP CRITICAL
    PRINT*, 'Found negative ', N, State_Chm%Species(N)%Conc(I,J,L)
    !$OMP END CRITICAL
    State_Chm%Species(N)%Conc(I,J,L) = 1e-26_fp

#else

    ! CO2 may be negative (e.g. from the balanced biosphere),
    ! so don't throw an error in this case.
    IF ( N == id_CO2 ) RETURN

    ! Otherwise print a warning and set is_negative=T so that 
    ! we can trap the error outside of any parallel loops
    !$OMP CRITICAL
    PRINT*, 'WARNING: Negative concentration for species ',                  &
         TRIM( State_Chm%SpcData(N)%Info%Name ), ' at (I,J,L) = ', I, J, L
    !$OMP END CRITICAL
    isNegative = .TRUE.

#endif

  END SUBROUTINE Check_For_Negative_Conc
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_chemgrid_only
!
! !DESCRIPTION: Returns .TRUE. if emissions of species N must be confined
!  to the chemistry grid (selected VOCs always; Bry species when linear
!  stratospheric chemistry is on).
!\\
!\\
! !INTERFACE:
!
  FUNCTION Is_ChemGrid_Only( Input_Opt, N ) RESULT( chem_grid_only )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt     ! Input options
    INTEGER,        INTENT(IN) :: N             ! Species index
!
! !RETURN VALUE:
!
    LOGICAL                    :: chem_grid_only  ! Restrict to chem grid?
!
! !REVISION HISTORY:
!  16 Oct 2026 - R. Yantosca - Initial version, from DO_TEND
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !========================================================================
    ! Is_ChemGrid_Only begins here!
    !========================================================================

    ! Avoid build-up of these VOCs in the stratosphere
    chem_grid_only = ( N == id_MACR .or. N == id_RCHO .or.                   &
                       N == id_ACET .or. N == id_ALD2 .or.                   &
                       N == id_ALK4 .or. N == id_C2H6 .or.                   &
                       N == id_C3H8 .or. N == id_CH2O .or.                   &
                       N == id_PRPE                                         )

    ! Bry is prescribed above the chemistry grid in linear strat chemistry
    IF ( Input_Opt%LINEAR_CHEM ) THEN
       chem_grid_only = chem_grid_only                   .or.                &
                        ( N == id_BrO  .or. N == id_Br2  .or.                &
                          N == id_Br   .or. N == id_HOBr .or.                &
                          N == id_HBr  .or. N == id_BrNO3                   )
    ENDIF

  END FUNCTION Is_ChemGrid_Only
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_emis_top
!
! !DESCRIPTION: Returns the top level for emissions in column (I,J): the
!  chemistry grid top if chemGridOnly is .TRUE., otherwise the model top.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Emis_Top( State_Grid, State_Met, chem_grid_only,            &
                           I,          J,         emis_top                  )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid       ! Grid state
    TYPE(MetState), INTENT(IN)  :: State_Met        ! Met state
    LOGICAL,        INTENT(IN)  :: chem_grid_only   ! Restrict to chem grid?
    INTEGER,        INTENT(IN)  :: I, J             ! Lon and lat indices
!
! !OUTPUT PARAMETERS
!
    INTEGER,        INTENT(OUT) :: emis_top         ! Top emission level
!
! !REVISION HISTORY:
!  16 Oct 2026 - R. Yantosca - Initial version, from DO_TEND
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !========================================================================
    ! Get_Emis_Top begins here!
    !========================================================================

    ! Emit to the top of the chemistry grid
    IF ( chem_grid_only ) THEN
       emis_top = MIN( State_Grid%NZ, State_Met%ChemGridLev(I,J) )
       RETURN
    ENDIF

    ! Otherwise emit to the top of the atmosphere
    emis_top = State_Grid%NZ

  END SUBROUTINE Get_Emis_Top
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tend_species_ids
!
! !DESCRIPTION: Looks up the species indices used by DO_TEND and
!  Do_Emis_Above_Sfc.  Only does work on the first call.  Must be called
!  outside of parallel loops.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Tend_Species_Ids()
!
! !USES:
!
    USE State_Chm_Mod, ONLY : Ind_
!
! !REVISION HISTORY:
!  16 Oct 2026 - R. Yantosca - Initial version, from DO_TEND
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    !========================================================================
    ! Init_Tend_Species_Ids begins here!
    !========================================================================

    ! Exit if this has been done previously
    IF ( ids_are_set ) RETURN

    ! Initialize
    id_ACET     = Ind_( 'ACET'  )
    id_ALD2     = Ind_( 'ALD2'  )
    id_ALK4     = Ind_( 'ALK4'  )
    id_Br       = Ind_( 'Br'    )
    id_Br2      = Ind_( 'Br2'   )
    id_BrNO3    = Ind_( 'BrNO3' )
    id_BrO      = Ind_( 'BrO'   )
    id_C2H6     = Ind_( 'C2H6'  )
    id_C3H8     = Ind_( 'C3H8'  )
    id_CH2O     = Ind_( 'CH2O'  )
    id_CO2      = Ind_( 'CO2'   )
    id_HBr      = Ind_( 'HBr'   )
    id_HNO3     = Ind_( 'HNO3'  )
    id_HOBr     = Ind_( 'HOBr'  )
    id_MACR     = Ind_( 'MACR'  )
    id_O3       = Ind_( 'O3'    )
    id_PRPE     = Ind_( 'PRPE'  )
    id_RCHO     = Ind_( 'RCHO'  )
    ids_are_set = .TRUE.

  END SUBROUTINE Init_Tend_Species_Ids
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tend_budget
!
! !DESCRIPTION: Computes the BudgetEmisDryDep* diagnostics.  Call with
!  beforeOp=.TRUE. before applying tendencies (stores initial column
!  masses), and with beforeOp=.FALSE. and timeStep afterwards (stores the
!  change in column mass [kg/s]).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Budget( Input_Opt,  State_Chm, State_Diag,                 &
                          State_Grid, State_Met, beforeOp,   RC,             &
                          timeStep                                          )
!
! !USES:
!
    USE Diagnostics_Mod, ONLY : Compute_Budget_Diagnostics
    USE ErrCode_Mod
    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Diag_Mod,  ONLY : DgnState
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Met_Mod,   ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)           :: Input_Opt   ! Input options
    TYPE(GrdState), INTENT(IN)           :: State_Grid  ! Grid state
    TYPE(MetState), INTENT(IN)           :: State_Met   ! Met state
    LOGICAL,        INTENT(IN)           :: beforeOp    ! Before tendencies?
    REAL(f8),       INTENT(IN), OPTIONAL :: timeStep    ! Time step [s]
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT)        :: State_Chm   ! Chemistry state
    TYPE(DgnState), INTENT(INOUT)        :: State_Diag  ! Diagnostics state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)          :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  16 Oct 2026 - R. Yantosca - Initial version, from DO_TEND
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: errMsg, thisLoc

    RC      = GC_SUCCESS
    thisLoc = ' -> at Tend_Budget (in module GeosCore/mixing_mod.F90)'

    IF ( .not. State_Diag%Archive_BudgetEmisDryDep ) RETURN

    IF ( beforeOp ) THEN

       ! Get initial column masses (full, trop, PBL, levels)
       CALL Compute_Budget_Diagnostics(                                      &
            Input_Opt   = Input_Opt,                                         &
            State_Chm   = State_Chm,                                         &
            State_Diag  = State_Diag,                                        &
            State_Grid  = State_Grid,                                        &
            State_Met   = State_Met,                                         &
            isFull      = State_Diag%Archive_BudgetEmisDryDepFull,           &
            diagFull    = NULL(),                                            &
            mapDataFull = State_Diag%Map_BudgetEmisDryDepFull,               &
            isTrop      = State_Diag%Archive_BudgetEmisDryDepTrop,           &
            diagTrop    = NULL(),                                            &
            mapDataTrop = State_Diag%Map_BudgetEmisDryDepTrop,               &
            isPBL       = State_Diag%Archive_BudgetEmisDryDepPBL,            &
            diagPBL     = NULL(),                                            &
            mapDataPBL  = State_Diag%Map_BudgetEmisDryDepPBL,                &
            isLevs      = State_Diag%Archive_BudgetEmisDryDepLevs,           &
            diagLevs    = NULL(),                                            &
            mapDataLevs = State_Diag%Map_BudgetEmisDryDepLevs,               &
            colMass     = State_Diag%BudgetColumnMass,                       &
            before_op   = .TRUE.,                                            &
            RC          = RC                                                )
       errMsg = 'Emissions/dry deposition budget diagnostics error 1'

    ELSE

       ! Compute change in column masses and store in diagnostics [kg/s]
       CALL Compute_Budget_Diagnostics(                                      &
            Input_Opt   = Input_Opt,                                         &
            State_Chm   = State_Chm,                                         &
            State_Diag  = State_Diag,                                        &
            State_Grid  = State_Grid,                                        &
            State_Met   = State_Met,                                         &
            isFull      = State_Diag%Archive_BudgetEmisDryDepFull,           &
            diagFull    = State_Diag%BudgetEmisDryDepFull,                   &
            mapDataFull = State_Diag%Map_BudgetEmisDryDepFull,               &
            isTrop      = State_Diag%Archive_BudgetEmisDryDepTrop,           &
            diagTrop    = State_Diag%BudgetEmisDryDepTrop,                   &
            mapDataTrop = State_Diag%Map_BudgetEmisDryDepTrop,               &
            isPBL       = State_Diag%Archive_BudgetEmisDryDepPBL,            &
            diagPBL     = State_Diag%BudgetEmisDryDepPBL,                    &
            mapDataPBL  = State_Diag%Map_BudgetEmisDryDepPBL,                &
            isLevs      = State_Diag%Archive_BudgetEmisDryDepLevs,           &
            diagLevs    = State_Diag%BudgetEmisDryDepLevs,                   &
            mapDataLevs = State_Diag%Map_BudgetEmisDryDepLevs,               &
            colMass     = State_Diag%BudgetColumnMass,                       &
            timeStep    = timeStep,                                          &
            RC          = RC                                                )
       errMsg = 'Emissions/dry deposition budget diagnostics error 2'

    ENDIF

    IF ( RC /= GC_SUCCESS ) CALL GC_Error( errMsg, RC, thisLoc )

  END SUBROUTINE Tend_Budget
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: tend_convert_units
!
! !DESCRIPTION: Converts the advected species to new\_units, halting the
!  "Boundary layer mixing" timer so that the conversion is timed separately.
!  Pass previous\_units on the way in to save the original units; omit it
!  on the way back.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Convert_Units( Input_Opt, State_Chm, State_Grid,           &
                                 State_Met, new_units, RC,                   &
                                 previous_units                             )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE Timers_Mod,     ONLY : Timer_End, Timer_Start
    USE UnitConv_Mod,   ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)            :: Input_Opt      ! Input options
    TYPE(GrdState), INTENT(IN)            :: State_Grid     ! Grid state
    TYPE(MetState), INTENT(IN)            :: State_Met      ! Met state
    INTEGER,        INTENT(IN)            :: new_units      ! Target units
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT)         :: State_Chm      ! Chemistry state
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)           :: RC             ! Success?
    INTEGER,        INTENT(OUT), OPTIONAL :: previous_units ! Original units
!
! !REVISION HISTORY:
!  16 Oct 2026 - R. Yantosca - Initial version, from DO_TEND
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: errMsg, thisLoc

    RC      = GC_SUCCESS
    thisLoc = ' -> at Tend_Convert_Units (in module GeosCore/mixing_mod.F90)'

#if defined( ADJOINT ) && defined( DEBUG )
    IF ( Input_Opt%Is_Adjoint .and. Input_Opt%IS_FD_SPOT_THIS_PET ) THEN
       WRITE( *, * ) ' SpcAdj(IFD,JFD) before unit conversion: ',           &
            State_Chm%SpeciesAdj( Input_Opt%IFD, Input_Opt%JFD,              &
                                  Input_Opt%LFD, Input_Opt%NFD )
       WRITE( *, * ) ' Spc(IFD,JFD) before unit conversion: ',              &
            State_Chm%Species(Input_Opt%NFD)%Conc( Input_Opt%IFD,            &
                                                   Input_Opt%JFD,            &
                                                   Input_Opt%LFD )
    ENDIF
#endif

    ! Halt mixing timer (so that unit conv can be timed separately)
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End( "Boundary layer mixing", RC )
    ENDIF

    ! Convert units (previous_units may be absent; it is OPTIONAL there too)
    CALL Convert_Spc_Units(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         mapping        = State_Chm%Map_Advect,                              &
         new_units      = new_units,                                         &
         previous_units = previous_units,                                    &
         RC             = RC                                                )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Unit conversion error!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Start mixing timer again
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_Start( "Boundary layer mixing", RC )
    ENDIF

#if defined( ADJOINT ) && defined( DEBUG )
    IF ( Input_Opt%Is_Adjoint .and. Input_Opt%IS_FD_SPOT_THIS_PET ) THEN
       WRITE( *, * ) ' SpcAdj(IFD,JFD) after unit conversion: ',            &
            State_Chm%SpeciesAdj( Input_Opt%IFD, Input_Opt%JFD,              &
                                  Input_Opt%LFD, Input_Opt%NFD )
       WRITE( *, * ) ' Spc(IFD,JFD) after unit conversion: ',               &
            State_Chm%Species(Input_Opt%NFD)%Conc( Input_Opt%IFD,            &
                                                   Input_Opt%JFD,            &
                                                   Input_Opt%LFD )
    ENDIF
#endif

  END SUBROUTINE Tend_Convert_Units
!EOC
END MODULE MIXING_MOD
