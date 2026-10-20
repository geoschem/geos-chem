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
  PUBLIC :: Do_Mixing
  PUBLIC :: Do_Tend
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
  SUBROUTINE Do_Mixing( Input_Opt,  State_Chm, State_Diag,                   &
                        State_Grid, State_Met, RC                           )
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
!  While all dry deposition rates are calculated either in DO_VDIFF or
!  DO_TEND, settling of aerosols is still computed in the dust & seasalt
!  modules.  The order of operations for full PBL mixing (TURBDAY) or no
!  PBL mixing stays the same, so TURBDAY and no-mixing stay zero-diff.
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
       CALL Do_Emis_Above_Sfc(                                               &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Diag = State_Diag,                                         &
            State_Grid = State_Grid,                                         &
            State_Met  = State_Met,                                          &
            RC         = RC                                                 )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Do_Emis_Above_Sfc"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Call VDIFF to perform the non-local PBL mixing.  The non-local term
       ! (surface emissions and dry deposition flux) enters VDIFF via the
       ! State_Chm%SurfaceFlux array (defined in Compute_Sflx_for_Vdiff,
       ! located in module GeosCore/hco_interface_gc_mod.F90).
       CALL Do_Vdiff(                                                        &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Diag = State_Diag,                                         &
            State_Grid = State_Grid,                                         &
            State_Met  = State_Met,                                          &
            RC         = RC                                                 )
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
    CALL Do_Tend(                                                            &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Diag = State_Diag,                                            &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         RC         = RC                                                    )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "DO_TEND"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! If PBL mixing is on, fully mix the PBL (TURBDAY)
    IF ( Input_Opt%LTURB ) THEN
       CALL Do_Full_Pbl_Mixing(                                              &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Diag = State_Diag,                                         &
            State_Grid = State_Grid,                                         &
            State_Met  = State_Met,                                          &
            RC         = RC                                                 )
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
  SUBROUTINE Do_Tend( Input_Opt, State_Chm, State_Diag, State_Grid,          &
                      State_Met, RC,        DT                              )
!
! !USES:
!
    USE ErrCode_Mod
    USE Error_Mod,            ONLY : Safe_Div
    USE Get_Ndep_Mod,         ONLY : Soil_DryDep
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_GetDiagn
    USE HCO_Utilities_GC_Mod, ONLY : GetHcoValDep
    USE HCO_Utilities_GC_Mod, ONLY : InquireHco
    USE HCO_Utilities_GC_Mod, ONLY : LoadHcoValDep
    USE HCO_Utilities_GC_Mod, ONLY : LoadHcoValEmis
    USE Input_Opt_Mod,        ONLY : OptInput
    USE PhysConstants,        ONLY : AVO
    USE Species_Mod,          ONLY : Species
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE Time_Mod,             ONLY : Get_Ts_Dyn
    USE UnitConv_Mod
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
    LOGICAL                :: chem_grid_only, found
    LOGICAL                :: error,          has_dep_flux
    LOGICAL                :: has_emis_flux,  is_negative
    INTEGER                :: dryd_top,       emis_top
    INTEGER                :: I,              id_drydep
    INTEGER                :: J,              L
    INTEGER                :: N,              NA
    INTEGER                :: pbl_top,        previous_units
    INTEGER                :: S
    REAL(fp)               :: air_sea_freq,   denom
    REAL(fp)               :: flux,           frac
    REAL(fp)               :: freq,           lnox_kgm2
    REAL(fp)               :: lnox_kgm2s,     mw_kg
    REAL(fp)               :: RkT,            time_step

    ! Pointers and objects
    REAL(f4),      POINTER :: lnox_kgm2s_o3  (:,:)
    REAL(f4),      POINTER :: lnox_kgm2s_hno3(:,:)
    REAL(f4),      POINTER :: ptr2D          (:,:)
    TYPE(Species), POINTER :: SpcInfo

    ! Strings
    CHARACTER(LEN=255)     :: thisLoc
    CHARACTER(LEN=512)     :: errMsg

    !========================================================================
    ! Do_Tend begins here!
    !========================================================================

    ! Assume success
    RC            =  GC_SUCCESS
    errMsg        =  ''
    thisLoc       =  ' -> at DO_TEND (in module GeosCore/mixing_mod.F90)'

    ! Special case that there is no dry deposition and emissions
    IF ( .NOT. Input_Opt%LDRYD .AND. .NOT. Input_Opt%DoEmissions ) RETURN

    ! Initialize
    error           = .FALSE.
    SpcInfo         => NULL()
    lnox_kgm2s_o3   => NULL()
    lnox_kgm2s_hno3 => NULL()
    ptr2d           => NULL()

    ! Look up species indices (first call only)
    CALL Init_Tend_Species_Ids()

    !------------------------------------------------------------------------
    ! Emissions/dry deposition budget diagnostics - Part 1 of 2
    !------------------------------------------------------------------------
    CALL Tend_Budget(                                                        &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Diag = State_Diag,                                            &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         before     = .TRUE.,                                                &
         RC         = RC                                                    )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Tend_Budget"!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Convert species to [kg/m2] (removes area dependency, ewl 9/30/15)
    !------------------------------------------------------------------------
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

    ! Get time step [s]
    IF ( PRESENT( dt ) ) THEN
       time_step = dt
    ELSE
       time_step = Get_Ts_Dyn()
    ENDIF
#ifdef ADJOINT
    IF ( Input_Opt%Is_Adjoint ) THEN
       time_step = time_step * (-1)
    ENDIF
#endif

    !========================================================================
    ! Get pointers to the PARANOX loss fluxes. These are
    ! stored in diagnostics 'PARANOX_O3_DEPOSITION_FLUX' and
    ! 'PARANOX_HNO3_DEPOSITION_FLUX'.  The call below links pointers
    ! lnox_kgm2s_o3 and lnox_kgm2s_hno3 to the data values stored in the
    ! respective diagnostics. The pointers will remain unassociated if
    ! the diagnostics do not exist.  This is only needed if non-local PBL
    ! scheme is not being used. Otherwise, PARANOX fluxes are applied
    ! in vdiff_mod.F90. (CHristoph Keller, 4/10/2015)
    !
    ! If using HEMCO Intermediate grid feature, then the call needs to be
    ! refreshed at every time step for regridding. (Haipeng Lin, 6/21/20)
    !========================================================================
    IF ( .NOT. Input_Opt%LNLPBL ) THEN

       ! O3 flux from PARANOX
       CALL HCO_GC_GetDiagn(                                                 &
            Input_Opt      = Input_Opt,                                      &
            State_Grid     = State_Grid,                                     &
            DiagnName      = 'PARANOX_O3_DEPOSITION_FLUX',                   &
            StopIfNotFound = .FALSE.,                                        &
            ptr2d          = ptr2d,                                          &
            RC             = RC                                             )
       IF ( ASSOCIATED( ptr2d ) ) THEN
          ALLOCATE( lnox_kgm2s_o3( State_Grid%NX, State_Grid%NY ), STAT=RC )
          lnox_kgm2s_o3 = ptr2d
       ENDIF
       ptr2d => NULL()

       ! HNO3 flux from PARANOX
       CALL HCO_GC_GetDiagn(                                                 &
            Input_Opt      = Input_Opt,                                      &
            State_Grid     = State_Grid,                                     &
            DiagnName      = 'PARANOX_HNO3_DEPOSITION_FLUX',                 &
            StopIfNotFound = .FALSE.,                                        &
            ptr2d          = ptr2d,                                          &
            RC             = RC                                             )
       IF( ASSOCIATED( ptr2d ) ) THEN
          ALLOCATE( lnox_kgm2s_hno3( State_Grid%NX, State_Grid%NY ), STAT=RC )
          lnox_kgm2s_hno3 = ptr2d
       ENDIF
       ptr2d => NULL()
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
    DO NA = 1, State_Chm%nAdvect                   ! Advected species ID

       ! Initialize various species quantities
       N              =  State_Chm%Map_Advect(NA)   ! Species ID
       SpcInfo        => State_Chm%SpcData(N)%Info  ! Species database entry
       mw_kg          =  SpcInfo%MW_g * 1.0e-3_fp   ! Species mol. wt. [kg/mol]
       has_dep_flux   = .FALSE.                     ! Does it have deposition?
       has_emis_flux  = .FALSE.                     ! Does it have emissions?
       id_drydep      = -1                          ! Drydep species index

       ! Should we restrict emissions of species N to the chemistry grid
       ! (to avoid building up concentration near the top-of-atmosphere)?
       chem_grid_only = Is_ChemGrid_Only( Input_Opt, N )

       !---------------------------------------------------------------------
       ! Check if this species has a surface deposition flux
       !---------------------------------------------------------------------
       IF ( Input_Opt%LDRYD ) THEN

          ! Check if this is one of the drydep species from drydep_mod.F90
          id_drydep    = SpcInfo%DryDepId
          has_dep_flux = ( id_drydep > 0 )

          ! Also check if the species has air-to-sea flux from HEMCO
          IF ( .not. has_dep_flux ) THEN
             CALL InquireHco ( N, dep=has_dep_flux )
          ENDIF

          ! Special case for O3 or HNO3: also include PARANOX loss
          IF ( N == id_O3   .and. ASSOCIATED( lnox_kgm2s_o3   ) ) THEN
             has_dep_flux = .TRUE.
          ENDIF
          IF ( N == id_HNO3 .and. ASSOCIATED( lnox_kgm2s_hno3 ) ) THEN
             has_dep_flux = .TRUE.
          ENDIF
       ENDIF

       !---------------------------------------------------------------------
       ! Check if this species has an emissions flux
       !---------------------------------------------------------------------
       IF ( Input_Opt%DoEmissions ) THEN
          CALL InquireHco ( N, emis=has_emis_flux )
       ENDIF

       !---------------------------------------------------------------------
       ! Skip unless this species has drydep and/or emissions fluxes
       !---------------------------------------------------------------------
       IF ( .not. has_dep_flux .and. .not. has_emis_flux ) CYCLE

       !---------------------------------------------------------------------
       ! Emissions or deposition (from HEMCO) must be loaded into memory
       ! first.  This is achieved by attempting to retrieve a grid box
       ! while NOT in a parallel loop. Failure to load this will result
       ! in severe performance issues!! (Haipeng Lin, 9/27/20)
       !---------------------------------------------------------------------
       IF ( has_emis_flux ) THEN
          CALL LoadHcoValEmis( Input_Opt, State_Grid, N )
       ENDIF
       IF ( has_dep_flux ) THEN
          CALL LoadHcoValDep( Input_Opt, State_Grid, N )
       ENDIF

       !=====================================================================
       ! Loop over surface grid boxes
       !=====================================================================
       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( I,            J,          pbl_top, dryd_top           )&
       !$OMP PRIVATE( emis_top,     L,          freq,    found              )&
       !$OMP PRIVATE( air_sea_freq, lnox_kgm2s, RkT,     frac               )&
       !$OMP PRIVATE( flux,         lnox_kgm2,  denom,   S                  )&
       !$OMP PRIVATE( is_negative                                           )&
       !$OMP COLLAPSE( 2                                                    )&
       !$OMP SCHEDULE( GUIDED                                               )&
       !$OMP REDUCTION( .OR. : error                                        )
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Continue to end of loop if an error has occured.  Then we can
          ! exit outside the parallel loop.  This is thread-safe.
          IF ( error ) CYCLE

          !------------------------------------------------------------------
          ! Define various quantities before computing tendencies
          !------------------------------------------------------------------

          ! Get the PBL top level at this grid box
          pbl_top = MAX( 1, FLOOR( State_Met%PBL_TOP_L(I,J) ) )

          ! Top level for dry deposition
          IF ( Input_Opt%PBL_DRYDEP ) THEN
             dryd_top = pbl_top
          ELSE
             dryd_top = 1
          ENDIF

          ! Find the top level for emissions
          CALL Get_Emis_Top(                                                 &
               I              = I,                                           &
               J              = J,                                           &
               chem_grid_only = chem_grid_only,                              &
               State_Grid     = State_Grid,                                  &
               State_Met      = State_Met,                                   &
               emis_top       = emis_top                                    )

          ! Loop up to the top level for drydep or emissions
          ! (whichever is higher)
          DO L = 1, MAX( dryd_top, emis_top )

             !---------------------------------------------------------------
             ! If the species has a dry deposition flux,
             ! then apply it to all levels below the PBL top
             !---------------------------------------------------------------
             IF ( has_dep_flux .AND. ( L <= dryd_top ) ) THEN

                ! Get the drydep frequency [s-1] from drydep_mod.F90.
                freq = 0.0_fp
                IF ( id_drydep > 0 ) THEN
                   freq = State_Chm%DryDepFreq(I,J,id_drydep)
                ENDIF

                ! Then add the air-to-sea deposition frequency [s-1]
                ! from the HEMCO SeaFlux extension if it exists.
                CALL GetHcoValDep(                                           &
                     I          = I,                                         &
                     J          = J,                                         &
                     L          = 1,                                         &
                     TrcId      = N,                                         &
                     Input_Opt  = Input_Opt,                                 &
                     State_Grid = State_Grid,                                &
                     found      = found,                                     &
                     dep        = air_sea_freq                              )
                IF ( found ) freq = freq + air_sea_freq

                ! For O3 and HNO3, also get the PARANOX deposition
                ! loss [kg/m2/s] at the surface (if it exists).
                lnox_kgm2s = 0.0_fp
                IF ( L == 1 ) THEN
                   IF ( N == id_O3   .and. ASSOCIATED( lnox_kgm2s_o3   ) ) THEN
                      lnox_kgm2s = lnox_kgm2s_o3(I,J)
                   ENDIF
                   IF ( N == id_HNO3 .and. ASSOCIATED( lnox_kgm2s_hno3 ) ) THEN
                      lnox_kgm2s = lnox_kgm2s_hno3(I,J)
                   ENDIF
                ENDIF

                !------------------------------------------------------------
                ! Proceed to apply dry depostion flux if it is nonzero
                !------------------------------------------------------------
                IF ( freq > 0.0_fp .OR. lnox_kgm2s > 0.0_fp ) THEN

                   ! Compute exponential loss term
                   RkT  = freq * time_step
                   frac = EXP( -RkT )

                   ! Loss in kg/m2
                   flux = ( 1.0_fp - frac ) * State_Chm%Species(N)%Conc(I,J,L)

                   ! Apply dry deposition
                   State_Chm%Species(N)%Conc(I,J,L) = frac *                 &
                      State_Chm%Species(N)%Conc(I,J,L)

#ifdef ADJOINT
                   IF ( Input_Opt%Is_Adjoint ) THEN
                      State_Chm%SpeciesAdj(I,J,L,N) = frac *                 &
                         State_Chm%SpeciesAdj(I,J,L,N)
                   ENDIF
#endif

                   ! Apply PARANOX loss [kg/m2/s]
                   IF ( lnox_kgm2s > 0.0_fp ) THEN
                      lnox_kgm2 = lnox_kgm2s * time_step

                      State_Chm%Species(N)%Conc(I,J,L) =                     &
                      State_Chm%Species(N)%Conc(I,J,L) - lnox_kgm2

                      flux = flux + lnox_kgm2
                   ENDIF

                   ! Loss in [molec/cm2/s]
                   denom = ( mw_kg * time_step * 1.0e+4_fp ) / AVO
                   flux  = Safe_Div( flux, denom, 0.0_fp )  ! [molec/cm2/s]

                   ! Eventually add to SOIL_DRYDEP
                   IF ( Input_Opt%LSOILNOX ) THEN
                      CALL SOIL_DRYDEP( I, J, N, flux, State_Chm )
                   ENDIF

                   !---------------------------------------------------------
                   ! HISTORY: Archive drydep flux loss from mixing
                   ! Units = molec/cm2/s
                   !
                   ! NOTE: we don't need to multiply by the ratio of
                   ! TS_CONV / TS_CHEM, as the updating frequency for
                   ! HISTORY is determined by the "frequency" setting in
                   ! the "HISTORY.rc" input file.
                   !---------------------------------------------------------
                   IF ( ( State_Diag%Archive_DryDepMix .or.                  &
                          State_Diag%Archive_DryDep        )   .and.         &
                          id_drydep > 0                       ) THEN
                      S = State_Diag%Map_DryDepMix%id2slot(id_drydep)
                      IF ( S > 0 ) THEN
                         State_Diag%DryDepMix(I,J,S) = flux
                      ENDIF
                   ENDIF

                ENDIF ! apply drydep
             ENDIF    ! L <= PBLTOP

             !---------------------------------------------------------------
             ! Apply emissions from HEMCO [kg/m2/s], converted to [kg/m2]
             !---------------------------------------------------------------
             IF ( has_emis_flux .AND. ( L <= emis_top ) ) THEN
                CALL Add_Emis_To_Box(                                        &
                     I           = I,                                        &
                     J           = J,                                        &
                     L           = L,                                        &
                     N           = N,                                        &
                     time_step   = time_step,                                &
                     Input_Opt   = Input_Opt,                                &
                     State_Grid  = State_Grid,                               &
                     State_Chm   = State_Chm                                )
             ENDIF

             !---------------------------------------------------------------
             ! Exit this loop if negatives are found.  We will trap the
             ! error outside of the parallel loop in a thread-safe manner.
             !---------------------------------------------------------------
             CALL Check_For_Negative_Conc(                                   &
                  I           = I,                                           &
                  J           = J,                                           &
                  L           = L,                                           &
                  N           = N,                                           &
                  State_Chm   = State_Chm,                                   &
                  is_negative = is_negative                                 )
             IF ( is_negative ) THEN
                error      = .TRUE.
                EXIT
             ENDIF

          ENDDO !L
       ENDDO    !J
       ENDDO    !I
       !$OMP END PARALLEL DO

       !---------------------------------------------------------------------
       ! Exit if a negative value was found on any thread in the loop.
       ! This is a thread-safe implementation.
       !---------------------------------------------------------------------
       IF ( error ) THEN
          errMsg = 'Negative species concentrations encountered! '         // &
                   'This may be fixed by increasing the background '      // &
                   'concentration or by shortening the transport '        // &
                   'time step.'
          CALL GC_Error( errMsg, RC, thisLoc )

          ! Free pointers and return
          IF ( ASSOCIATED( lnox_kgm2s_o3   ) ) DEALLOCATE( lnox_kgm2s_o3   )
          IF ( ASSOCIATED( lnox_kgm2s_hno3 ) ) DEALLOCATE( lnox_kgm2s_hno3 )
          lnox_kgm2s_o3   => NULL()
          lnox_kgm2s_hno3 => NULL()
          SpcInfo         => NULL()
          RETURN
       ENDIF


    ENDDO !N

    !========================================================================
    ! Convert species from [kg/m2] back to original units
    !========================================================================
    CALL Tend_Convert_Units(                                                 &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         new_units      = previous_units,                                    &
         RC             = RC                                                )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Tend_Convert_Units" (from kg/m2)!'
       CALL GC_Error( errMsg, RC, thisLoc )
       IF ( ASSOCIATED( lnox_kgm2s_o3   ) ) DEALLOCATE( lnox_kgm2s_o3   )
       IF ( ASSOCIATED( lnox_kgm2s_hno3 ) ) DEALLOCATE( lnox_kgm2s_hno3 )
       lnox_kgm2s_o3   => NULL()
       lnox_kgm2s_hno3 => NULL()
       SpcInfo         => NULL()
       RETURN
    ENDIF

    !========================================================================
    ! Emissions/dry deposition budget diagnostics (part 2 of 2)
    !========================================================================
    CALL Tend_Budget(                                                        &
         Input_Opt    = Input_Opt,                                           &
         State_Chm    = State_Chm,                                           &
         State_Diag   = State_Diag,                                          &
         State_Grid   = State_Grid,                                          &
         State_Met    = State_Met,                                           &
         before       = .FALSE.,                                             &
         time_step    = DBLE( time_step ),                                   &
         RC           = RC                                                  )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Tend_Budget" #2!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       IF ( ASSOCIATED( lnox_kgm2s_o3   ) ) DEALLOCATE( lnox_kgm2s_o3   )
       IF ( ASSOCIATED( lnox_kgm2s_hno3 ) ) DEALLOCATE( lnox_kgm2s_hno3 )
       lnox_kgm2s_o3   => NULL()
       lnox_kgm2s_hno3 => NULL()
       SpcInfo         => NULL()
       RETURN
    ENDIF

    !========================================================================
    ! Cleanup and quit
    !========================================================================

    ! Deallocate pointers
    IF ( ASSOCIATED( lnox_kgm2s_o3   ) ) DEALLOCATE( lnox_kgm2s_o3   )
    IF ( ASSOCIATED( lnox_kgm2s_hno3 ) ) DEALLOCATE( lnox_kgm2s_hno3 )

    ! Free pointers
    lnox_kgm2s_o3   => NULL()
    lnox_kgm2s_hno3 => NULL()
    ptr2d           => NULL()
    SpcInfo         => NULL()

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
    USE HCO_Utilities_GC_Mod, ONLY : InquireHco
    USE HCO_Utilities_GC_Mod, ONLY : LoadHcoValEmis
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
    LOGICAL            :: chem_grid_only, error
    LOGICAL            :: has_emis_flux,  is_negative
    INTEGER            :: emis_top,       I
    INTEGER            :: J,              L
    INTEGER            :: N,              NA
    INTEGER            :: previous_units
    REAL(fp)           :: time_step

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
    ! Emissions/dry deposition budget diagnostics (part 1 of 2)
    !========================================================================
    CALL Tend_Budget(                                                        &
         Input_Opt   = Input_Opt,                                            &
         State_Chm   = State_Chm,                                            &
         State_Diag  = State_Diag,                                           &
         State_Grid  = State_Grid,                                           &
         State_Met   = State_Met,                                            &
         before      = .TRUE.,                                               &
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
    !========================================================================
    DO NA = 1, State_Chm%nAdvect

       ! Skip species without HEMCO emissions
       N = State_Chm%Map_Advect(NA)
       CALL InquireHco( N, emis=has_emis_flux )
       IF ( .not. has_emis_flux ) CYCLE

       ! Load emissions into memory outside of the parallel loop
       ! to avoid a performance penalty (Haipeng Lin, 9/27/20)
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
                EXIT
             ENDIF

          ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       !---------------------------------------------------------------------
       ! Throw an error if any threads returned with error status in
       ! the loop above.  This is the thread-safe implementation.
       !---------------------------------------------------------------------
       IF ( error ) THEN
          errMsg = 'Negative species concentrations encountered! '         // &
                   'This may be fixed by increasing the background '      // &
                   'concentration or by shortening the transport '        // &
                   'time step.'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDDO

    !========================================================================
    ! Convert species from [kg/m2] back to original units
    !========================================================================
    CALL Tend_Convert_Units(                                                 &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         new_units  = previous_units,                                        &
         RC         = RC                                                    )
    IF ( RC /= GC_SUCCESS ) THEN
       errMsg = 'Error encountered in "Tend_Convert_Units" (from kg/m2)!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Emissions/dry deposition budget diagnostics (part 2 of 2)
    !========================================================================
    CALL Tend_Budget(                                                        &
         Input_Opt   = Input_Opt,                                            &
         State_Chm   = State_Chm,                                            &
         State_Diag  = State_Diag,                                           &
         State_Grid  = State_Grid,                                           &
         State_Met   = State_Met,                                            &
         before      = .FALSE.,                                              &
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
    USE HCO_Utilities_GC_Mod, ONLY : Tmp_MDL
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
       emis  = Tmp_MDL(I,J,L)
    ELSE
#endif
       CALL GetHcoValEmis(                                                   &
            I          = I,                                                  &
            J          = J,                                                  &
            L          = L,                                                  &
            TrcId      = N,                                                  &
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
    ! Debug printout (adjoint only)
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
! !IROUTINE: check_for_negative_conc
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
    LOGICAL,        INTENT(OUT)   :: is_negative ! Was a negative value found?
!
! !REMARKS:
!  Thread-safe: writes only Conc(I,J,L) and the private flag is_negative.
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
    is_negative = .FALSE.

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
    is_negative = .TRUE.

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
!  masses), and with before=.FALSE. and time_step afterwards (stores the
!  change in column mass [kg/s]).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tend_Budget( Input_Opt,  State_Chm, State_Diag,                 &
                          State_Grid, State_Met, before,                     &
                          RC,         time_step                             )
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
    LOGICAL,        INTENT(IN)           :: before      ! Before tendencies?
    REAL(f8),       INTENT(IN), OPTIONAL :: time_step   ! Time step [s]
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

    !========================================================================
    ! Tend_Budget begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Tend_Budget (in module GeosCore/mixing_mod.F90)'

    ! Exit if budget diagnostics are disabled
    IF ( .not. State_Diag%Archive_BudgetEmisDryDep ) RETURN

    ! Compute diagnostics
    IF ( before ) THEN

       !---------------------------------------------------------------------
       ! Get initial column masses (full, trop, PBL, levels)
       !---------------------------------------------------------------------
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
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error encountered in "Compute_Budget_Diagnostics" '   // &
                   '(part 1 of 2)!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ELSE

       !---------------------------------------------------------------------
       ! Compute change in column masses and store in diagnostics [kg/s]
       !---------------------------------------------------------------------
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
            timeStep    = time_step,                                         &
            RC          = RC                                                )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error encountered in "Compute_Budget_Diagnostics" '   // &
                   '(part 2 of 2)!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

    ENDIF

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
  SUBROUTINE Tend_Convert_Units( Input_Opt,     State_Chm, State_Grid,       &
                                 State_Met,     new_units, RC,               &
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
    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !========================================================================
    ! Tend_Convert_Units begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Tend_Convert_Units (in module GeosCore/mixing_mod.F90)'

#if defined( ADJOINT ) && defined( DEBUG )
    ! Debug output (adjoint only)
    IF ( Input_Opt%Is_Adjoint .and. Input_Opt%IS_FD_SPOT_THIS_PET ) THEN
       WRITE( *, * ) ' SpcAdj(IFD,JFD) before unit conversion: ',            &
            State_Chm%SpeciesAdj( Input_Opt%IFD, Input_Opt%JFD,              &
                                  Input_Opt%LFD, Input_Opt%NFD )
       WRITE( *, * ) ' Spc(IFD,JFD) before unit conversion: ',               &
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
    ! Debug output (adjoint only)
    IF ( Input_Opt%Is_Adjoint .and. Input_Opt%IS_FD_SPOT_THIS_PET ) THEN
       WRITE( *, * ) ' SpcAdj(IFD,JFD) after unit conversion: ',             &
            State_Chm%SpeciesAdj( Input_Opt%IFD, Input_Opt%JFD,              &
                                  Input_Opt%LFD, Input_Opt%NFD )
       WRITE( *, * ) ' Spc(IFD,JFD) after unit conversion: ',                &
            State_Chm%Species(Input_Opt%NFD)%Conc( Input_Opt%IFD,            &
                                                   Input_Opt%JFD,            &
                                                   Input_Opt%LFD )
    ENDIF
#endif

  END SUBROUTINE Tend_Convert_Units
!EOC
END MODULE MIXING_MOD
