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
!  (A) While all dry deposition rates are calculated either in
!      DO_PBL_MIX2 or DO_TEND, settling of aerosols is still
!      computed in the dust/seasalt modules.
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
    ! Scalars
    LOGICAL            :: OnlyAbovePBL
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! DO_MIXING begins here!
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at DO_MIXING (in module GeosCore/mixing_mod.F90)'

    !-----------------------------------------------------------------------
    ! Do non-local PBL mixing. This will apply the species tendencies
    ! (emission fluxes and dry deposition rates) below the PBL.
    ! This is done for all species with defined emissions / dry
    ! deposition rates, including dust.
    !
    ! Set OnlyAbovePBL flag (used below by DO_TEND) to indicate that
    ! fluxes within the PBL have already been applied.
    ! ----------------------------------------------------------------------
    IF ( Input_Opt%LTURB .AND. Input_Opt%LNLPBL ) THEN

       !--------------------------------------------------------------------
       ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
       !
       ! Initialize the diagnostic array for the History Component.  This will
       ! prevent leftover values from being carried over to this timestep.
       ! (For example, if on the last iteration, the PBL height was higher than
       ! it is now, then we will have stored drydep fluxes up to that height,
       ! so we need to zero these out.)
       !--------------------------------------------------------------------

       ! Non-local mixing
       CALL Do_Vdiff( Input_Opt,  State_Chm, State_Diag,                     &
                      State_Grid, State_Met, RC                             )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountred in "DO_PBL_MIX_2"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Fluxes in PBL have been applied (non-local PBL mixing)
       OnlyAbovePBL = .TRUE.

    ELSE

       ! Fluxes in PBL have not been applied (full PBL mixing)
       OnlyAbovePBL = .FALSE.

    ENDIF

    !-----------------------------------------------------------------------
    ! Apply tendencies. This will apply dry deposition rates and
    ! emission fluxes below the PBL if it has not yet been done
    ! via the non-local PBL mixing. It also adds the emissions above
    ! the PBL to the species array. Emissions of some species may be
    ! capped at the tropopause to avoid build-up in stratosphere.
    !-----------------------------------------------------------------------

    ! Apply tendencies
    CALL DO_TEND( Input_Opt, State_Chm,  State_Diag, &
                  State_Grid, State_Met, OnlyAbovePBL, RC )

    ! Trap potential error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountred in "DO_TEND"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Do full pbl mixing. This fully mixes the updated species
    ! concentrations within the PBL.
    !
    ! Now also archive concentrations and calculate turbulence
    ! tendencies (ckeller, 7/15/2015)
    !-----------------------------------------------------------------------
    IF ( Input_Opt%LTURB .AND. .NOT. Input_Opt%LNLPBL ) THEN

       ! Full PBL mixing
       CALL Do_Full_Pbl_Mixing( Input_Opt,  State_Chm, State_Diag,            &
                                State_Grid, State_Met, RC                    )

       ! Trap potential error
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountred in "DO_PBL_MIX"!'
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
  SUBROUTINE DO_TEND( Input_Opt, State_Chm,    State_Diag, State_Grid, &
                      State_Met, OnlyAbovePBL, RC,         DT )
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
    LOGICAL,          INTENT(IN   )           :: OnlyAbovePBL ! Only above PBL?
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
    INTEGER                 :: I, J, L, L1, L2, N, D, NN, NA, nAdvect, S
    INTEGER                 :: DRYDEPID, previous_units
    INTEGER                 :: PBL_TOP, DRYD_TOP, EMIS_TOP
    REAL(fp)                :: TS, TMP, FRQ, RKT, FRAC, FLUX, AREA_M2
    REAL(fp)                :: MWkg, DENOM
    LOGICAL                 :: FND
    LOGICAL                 :: PBL_DRYDEP, LINEAR_CHEM, ChemGridOnly
    LOGICAL                 :: LEMIS,      LDRYD
    LOGICAL                 :: DryDepSpec, EmisSpec
    REAL(f8)                :: DT_Tend

    ! PARANOX loss fluxes (kg/m2/s). These are obtained from the
    ! HEMCO PARANOX extension via the diagnostics module.
    REAL(fp)                :: PNOXLOSS
    REAL(f4), POINTER       :: Ptr2D        (:,:) => NULL()
    REAL(f4), POINTER       :: PNOXLOSS_O3  (:,:)
    REAL(f4), POINTER       :: PNOXLOSS_HNO3(:,:)

    ! SAVEd scalars (defined on first call only)
    LOGICAL,           SAVE :: FIRST = .TRUE.
    INTEGER,           SAVE :: id_MACR,  id_RCHO,  id_ACET, id_ALD2
    INTEGER,           SAVE :: id_ALK4,  id_C2H6,  id_C3H8, id_CH2O
    INTEGER,           SAVE :: id_PRPE,  id_O3,    id_HNO3, id_BrO
    INTEGER,           SAVE :: id_Br2,   id_Br,    id_HOBr, id_HBr
    INTEGER,           SAVE :: id_BrNO3, id_CO2

    ! Pointers and objects
    TYPE(Species), POINTER  :: SpcInfo
    REAL(fp),      POINTER  :: DepFreq(:,:,:  )  ! IM, JM, nDryDep

    ! Strings
    CHARACTER(LEN=255)      :: ErrMsg, ErrorMsg, ThisLoc

#ifdef ADJOINT
    LOGICAL                 :: IS_ADJ
#endif

    !=================================================================
    ! DO_TEND begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ErrMsg  = ''
    ThisLoc = ' -> at DO_TEND (in module GeosCore/mixing_mod.F90)'

    ! Special case that there is no dry deposition and emissions
    IF ( .NOT. Input_Opt%LDRYD .AND. .NOT. Input_Opt%DoEmissions ) RETURN

    ! Initialize
    LINEAR_CHEM       = Input_Opt%LINEAR_CHEM
    LEMIS             = Input_Opt%DoEmissions
    LDRYD             = Input_Opt%LDRYD
    PBL_DRYDEP        = Input_Opt%PBL_DRYDEP
    nAdvect           = State_Chm%nAdvect

    ! Initialize pointer
    SpcInfo           => NULL()
    DepFreq           => State_Chm%DryDepFreq

    PNOxLoss_O3       => NULL()
    PNOxLoss_HNO3     => NULL()

    !------------------------------------------------------------------------
    ! Emissions/dry deposition budget diagnostics - Part 1 of 2
    !------------------------------------------------------------------------
    IF ( State_Diag%Archive_BudgetEmisDryDep ) THEN

       ! Get initial column masses (full, trop, PBL)
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
          ErrMsg = 'Emissions/dry deposition budget diagnostics error 1'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

#if defined( ADJOINT )  && defined ( DEBUG )
    IF (Input_Opt%is_adjoint .and. Input_Opt%IS_FD_SPOT_THIS_PET) THEN
       WRITE(*,*) ' SpcAdj(IFD,JFD) before unit converstion: ',  &
            State_Chm%SpeciesAdj(Input_Opt%IFD, Input_Opt%JFD, &
            Input_Opt%LFD, Input_Opt%NFD)
       WRITE(*,*) ' Spc(IFD,JFD) before unit converstion: ',  &
            State_Chm%Species(Input_Opt%NFD)%Conc(Input_Opt%IFD, Input_Opt%JFD, Input_Opt%LFD)
    ENDIF
#endif

    ! Halt mixing timer (so that unit conv can be timed separately)
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End( "Boundary layer mixing", RC )
    ENDIF

    ! DO_TEND previously operated in units of kg. The species arrays are in
    ! v/v for mixing, hence needed to convert before and after.
    ! Now use units kg/m2 as State_Chm%SPECIES units in DO_TEND to
    ! remove area-dependency (ewl, 9/30/15)
    CALL Convert_Spc_Units(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         mapping        = State_Chm%Map_Advect,                              &
         new_units      = KG_SPECIES_PER_M2,                                 &
         previous_units = previous_units,                                    &
         RC             = RC                                                )
    
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
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

    ! Trap potential error
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Get time step [s]
    IF ( PRESENT(DT) ) THEN
       TS = DT
    ELSE
       TS = GET_TS_DYN()
    ENDIF
#ifdef ADJOINT
    if (Input_Opt%Is_Adjoint) then
       TS = TS * -1
    endif
#endif

    ! First-time setup
    IF ( FIRST ) THEN

       ! Define species indices on the first call
       id_MACR = Ind_('MACR' )
       id_RCHO = Ind_('RCHO' )
       id_ACET = Ind_('ACET' )
       id_ALD2 = Ind_('ALD2' )
       id_ALK4 = Ind_('ALK4' )
       id_C2H6 = Ind_('C2H6' )
       id_C3H8 = Ind_('C3H8' )
       id_CH2O = Ind_('CH2O' )
       id_CO2  = Ind_('CO2'  )
       id_PRPE = Ind_('PRPE' )
       id_O3   = Ind_('O3'   )
       id_HNO3 = Ind_('HNO3' )
       id_BrO  = Ind_('BrO'  )
       id_Br2  = Ind_('Br2'  )
       id_Br   = Ind_('Br'   )
       id_HOBr = Ind_('HOBr' )
       id_HBr  = Ind_('HBr'  )
       id_BrNO3= Ind_('BrNO3')

       FIRST = .FALSE.
    ENDIF

    ! On first call, get pointers to the PARANOX loss fluxes. These are
    ! stored in diagnostics 'PARANOX_O3_DEPOSITION_FLUX' and
    ! 'PARANOX_HNO3_DEPOSITION_FLUX'. The call below links pointers
    ! PNOXLOSS_O3 and PNOXLOSS_HNO3 to the data values stored in the
    ! respective diagnostics. The pointers will remain unassociated if
    ! the diagnostics do not exist.
    ! This is only needed if non-local PBL scheme is not being used.
    ! Otherwise, PARANOX fluxes are applied in vdiff_mod.F90.
    !  (ckeller, 4/10/2015)
    !
    ! If using HEMCO Intermediate grid feature, then the call needs to be
    ! refreshed at every time step for regridding. (hplin, 6/21/20)
    IF ( .NOT. Input_Opt%LNLPBL ) THEN
      CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, 'PARANOX_O3_DEPOSITION_FLUX', &
                            .FALSE.,   RC, Ptr2D = Ptr2D          )
      IF( ASSOCIATED( Ptr2D )) THEN
        ALLOCATE ( PNOxLoss_O3( State_Grid%NX, State_Grid%NY ), STAT=RC )
        PNOxLoss_O3(:,:) = Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, 'PARANOX_HNO3_DEPOSITION_FLUX',&
                            .FALSE.,   RC, Ptr2D = Ptr2D        )
      IF( ASSOCIATED( Ptr2D )) THEN
        ALLOCATE ( PNOxLoss_HNO3( State_Grid%NX, State_Grid%NY ), STAT=RC )
        PNOxLoss_HNO3(:,:) = Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()
    ENDIF

    !=======================================================================
    ! Do for every advected species and grid box
    !=======================================================================
    ! Note: For GEOS-Chem Classic HEMCO "Intermediate Grid" feature,
    ! where HEMCO is running on a distinct grid from the model, the
    ! on-demand regridding is most optimized when contiguous accesses
    ! to GetHcoValEmis and GetHcoValDep are performed for the given species.
    ! Therefore, it is most optimal to call in the following fashion (IJKN)
    !    Emis(1,1,1,1) -> Emis(1,1,2,1) -> ... -> Dep(1,1,1,1) -> Dep(1,1,2,1)
    ! By switching emis/dep or the species # LAST, as either of these changing
    ! WILL trigger a new regrid and thrashing of the old buffer.
    !
    ! Therefore, the loop below has been adjusted to run serially for each
    ! species, and parallelizing the inner I, J loop instead (hplin, 6/27/20)
    ! Also, moved some non-I,J specific variables outside of the loop for optimization

    DO NA = 1, nAdvect

       ! Initialize PRIVATE error-handling variables
       ErrorMsg  = ''

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
       DryDepSpec = .FALSE.
       DryDepID   = -1

       ! Only if dry deposition is turned on and we do want to consider
       ! processes below the PBL...
       IF ( LDRYD .AND. .NOT. OnlyAbovePBL ) THEN

          ! Get dry deposition ID (used by drydep_mod.F90) for this species.
          ! This is now stored in the species database object. (bmy, 7/6/16)
          DryDepID = SpcInfo%DryDepId

          ! Check if this is a HEMCO drydep species
          DryDepSpec = ( DryDepId > 0 )
          IF ( .NOT. DryDepSpec ) THEN
             CALL InquireHco ( N, Dep=DryDepSpec )
          ENDIF

          ! Special case for O3 or HNO3: include PARANOX loss
          IF ( N == id_O3   .AND. ASSOCIATED(PNOXLOSS_O3  ) )    &
               DryDepSpec = .TRUE.
          IF ( N == id_HNO3 .AND. ASSOCIATED(PNOXLOSS_HNO3) )    &
               DryDepSpec = .TRUE.
       ENDIF

       ! Set emissions top level:
       ! This is the top of atmosphere unless concentration build-up
       ! in stratosphere wants to be avoided.
       ChemGridOnly = .FALSE.

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
       IF ( LINEAR_CHEM ) THEN
          IF ( N == id_BrO  .OR. N == id_Br2   .OR. &
               N == id_Br   .OR. N == id_HOBr  .OR. &
               N == id_HBr  .OR. N == id_BrNO3       ) THEN
             ChemGridOnly = .TRUE.
          ENDIF
       ENDIF

       !--------------------------------------------------------------------
       ! Check if we need to do emissions for this species
       !--------------------------------------------------------------------
       IF ( LEMIS ) THEN
          CALL InquireHco ( N, Emis=EmisSpec )
       ELSE
          EmisSpec = .FALSE.
       ENDIF

       ! If there is emissions for this species, it must be loaded into memory first.
       ! This is achieved by attempting to retrieve a grid box while NOT in a parallel
       ! loop. Failure to load this will result in severe performance issues!! (hplin, 9/27/20)
       IF ( EmisSpec ) THEN
          CALL LoadHcoValEmis ( Input_Opt, State_Grid, N )
       ENDIF

       IF ( DryDepSpec ) THEN
          CALL LoadHcoValDep ( Input_Opt, State_Grid, N )
       ENDIF

       !--------------------------------------------------------------------
       ! Can go to next species if this species does not have
       ! dry deposition and/or emissions
       !--------------------------------------------------------------------
       IF ( .NOT. DryDepSpec .AND. .NOT. EmisSpec ) CYCLE

!$OMP PARALLEL DO                                                           &
!$OMP DEFAULT( SHARED                                                     ) &
!$OMP PRIVATE( I,        J,            L,          L1,       L2           ) &
!$OMP PRIVATE( PBL_TOP,  FND,          TMP                                ) &
!$OMP PRIVATE( FRQ,      RKT,          FRAC,       FLUX,     Area_m2      ) &
!$OMP PRIVATE( DRYD_TOP, EMIS_TOP,     PNOXLOSS,   DENOM                  ) &
!$OMP PRIVATE( S,        ErrorMsg                                         )

       ! Loop over all grid boxes
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          !-----------------------------------------------------------------
          ! Define various quantities before computing tendencies
          !-----------------------------------------------------------------

          ! Get PBL_TOP at this grid box
          PBL_TOP = MAX( 1, FLOOR( State_Met%PBL_TOP_L(I,J) ) )

          ! Determine lower level L1 to be used:
          ! If specified so, apply emissions only above the PBL_TOP.
          ! This will also disable dry deposition.
          IF ( OnlyAbovePBL ) THEN
             L1 = PBL_TOP + 1
          ELSE
             L1 = 1
          ENDIF

          ! Set dry deposition top level based on PBL_DRYDEP flag of
          ! Input_Opt.
          IF ( PBL_DRYDEP ) THEN
             DRYD_TOP = PBL_TOP
          ELSE
             DRYD_TOP = 1
          ENDIF

          ! Restrict to chemistry grid
          IF ( ChemGridOnly ) THEN
             EMIS_TOP = State_Met%ChemGridLev(I,J)
             EMIS_TOP = MIN(State_Grid%NZ,EMIS_TOP)
          ELSE
             EMIS_TOP = State_Grid%NZ
          ENDIF

          ! L2 is the upper level index to loop over
          L2 = MAX(DRYD_TOP, EMIS_TOP)

          ! This should not happen:
          IF ( L2 < L1 ) CYCLE

          ! Loop over selected vertical levels
          DO L = L1, L2

             !--------------------------------------------------------------
             ! Apply dry deposition frequencies to all levels below the
             ! PBL top.
             !--------------------------------------------------------------
             IF ( DryDepSpec .AND. ( L <= DRYD_TOP ) ) THEN

                ! Init
                FRQ = 0.0_fp

                ! Dry deposition frequency from drydep_mod.F90. This is
                ! stored in State_Chm%DryDepFreq. Units are [s-1].
                IF ( DRYDEPID > 0 ) THEN
                   FRQ = DepFreq(I,J,DRYDEPID)
                ENDIF

                ! Dry deposition frequency from HEMCO. HEMCO calculates
                ! dry deposition frequencies for air-sea exchange and
                ! from ship NOx plume parameterization (PARANOx). The
                ! units are [s-1].
                CALL GetHcoValDep ( Input_Opt, State_Grid, N, I, J, 1, FND, TMP )

                ! Add to dry dep frequency from drydep_mod.F90
                IF ( FND ) FRQ = FRQ + TMP

                ! Get PARANOX deposition loss. Apply to surface level only.
                ! PNOXLOSS is in kg/m2/s. (ckeller, 4/10/15)
                PNOXLOSS = 0.0_fp
                IF ( L == 1 ) THEN
                   IF ( N == id_O3 .AND. ASSOCIATED(PNOXLOSS_O3) ) THEN
                      PNOXLOSS = PNOXLOSS_O3(I,J)
                   ENDIF
                   IF ( N == id_HNO3 .AND. ASSOCIATED(PNOXLOSS_HNO3) ) THEN
                      PNOXLOSS = PNOXLOSS_HNO3(I,J)
                   ENDIF
                ENDIF

                ! Apply dry deposition
                IF ( FRQ > 0.0_fp .OR. PNOXLOSS > 0.0_fp ) THEN

                   ! Compute exponential loss term
                   RKT  = FRQ * TS
                   FRAC = EXP(-RKT)

                   ! Loss in kg/m2
                   FLUX = ( 1.0_fp - FRAC ) * State_Chm%Species(N)%Conc(I,J,L)

                   ! Apply dry deposition
                   State_Chm%Species(N)%Conc(I,J,L) = FRAC *    &
                                            State_Chm%Species(N)%Conc(I,J,L)

#ifdef ADJOINT
                   if (Input_Opt%Is_Adjoint) then
                      State_Chm%SpeciesAdj(I,J,L,N) = FRAC *  &
                           State_Chm%SpeciesAdj(I,J,L,N)
                   endif
#endif
                   ! Eventually add PARANOX loss. PNOXLOSS is in kg/m2/s.
                   ! Make sure PARANOx loss is applied to tracers. (ckeller,
                   ! 3/29/16).
                   IF ( PNOXLOSS > 0 ) THEN
                      State_Chm%Species(N)%Conc(I,J,L) = &
                         State_Chm%Species(N)%Conc(I,J,L) - ( PNOXLOSS * TS )
                      FLUX = FLUX + ( PNOXLOSS * TS )
                   ENDIF

                   ! Loss in [molec/cm2/s]
                   ! Added a safe_div due to small parallelization error
                   ! (mdy, 5/15)
                   !
                   ! NOTE: The original computation was:
                   !   FLUX = FLUX / MWkg * AVO / TS / ( AREA_M2 * 1.0e4_fp ) ]
                   ! so the denominator as we had it was wrong.
                   ! Now corrected (elundgren, bmy, 6/12/15)
                   DENOM = ( MWkg * TS * 1.0e+4_fp ) / AVO
                   FLUX  = SAFE_DIV( FLUX, DENOM, 0.0e+0_fp )  ! molec/cm2/s

                   ! Eventually add to SOIL_DRYDEP
                   IF ( Input_Opt%LSOILNOX ) THEN
                      CALL SOIL_DRYDEP( I, J, N, FLUX, State_Chm )
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
                         State_Diag%DryDepMix(I,J,S) = Flux
                      ENDIF
                   ENDIF

                ENDIF ! apply drydep
             ENDIF ! L <= PBLTOP

             !--------------------------------------------------------------
             ! Apply emissions.
             ! These are always taken from HEMCO
             !--------------------------------------------------------------
             IF ( EmisSpec .AND. ( L <= EMIS_TOP ) ) THEN

                ! Get HEMCO emissions. Units are [kg/m2/s].
                ! Fix hplin: for intermediate grid, pass SkipCheck in a tight loop. Note that this assumes that adjacent
                ! calls to GetHcoValEmis are from the same species ID, or there will be big trouble. (hplin, 10/10/20)

#ifdef MODEL_CLASSIC
                IF ( Input_Opt%LIMGRID ) THEN
                  FND = .true.
                  TMP = TMP_MDL(I,J,L) ! this is a kludge for the tight loop optimization
                ELSE
#endif
                  CALL GetHcoValEmis ( Input_Opt, State_Grid, N, I, J, L, FND, TMP, SkipCheck=.true. )
#ifdef MODEL_CLASSIC
                ENDIF
#endif

                ! Add emissions (if any)
                ! Bug fix: allow negative fluxes. (ckeller, 4/12/17)
                !IF ( FND .AND. (TMP > 0.0_fp) ) THEN
                IF ( FND ) THEN

                   ! Flux: [kg/m2] = [kg m-2 s-1 ] x [s]
                   FLUX = TMP * TS
#ifdef ADJOINT
                   IF ( I .eq. Input_Opt%IFD .and. J .eq. Input_Opt%JFD .and. &
                        L .eq. Input_Opt%LFD .and. N .eq. Input_Opt%NFD) THEN
                      WRITE(*,*) ' GetHcoVal(IFD,JFD) = ', TMP,  ' FLUX = ', FLUX
                      IF ( Input_Opt%is_adjoint ) THEN
                         WRITE(*,*) ' SpeciesAdj(FD) = ', State_Chm%SpeciesAdj(I,J,L,N)
                      ENDIF
                   ENDIF
#endif

                   ! Add to species array
                   State_Chm%Species(N)%Conc(I,J,L) = &
                         State_Chm%Species(N)%Conc(I,J,L) + FLUX
                ENDIF
             ENDIF

             ! Check for negative concentrations
             IF ( State_Chm%Species(N)%Conc(I,J,L) < 0.0_fp ) THEN
#ifdef TOMAS
                ! For TOMAS simulations only, look for negative and reset
                ! to small positive.  This prevents the run from dying,
                ! while we look for the root cause of the issue.
                !  -- Betty Croft, Bob Yantosca (21 Jan 2022)
                print *, 'Found negative ', N, State_Chm%Species(N)%Conc(I,J,L)
                State_Chm%Species(N)%Conc(I,J,L) = 1e-26_fp
#else

                IF ( N /= id_CO2 ) THEN
                   Print*, 'WARNING: Negative concentration for species ',    &
                            TRIM( SpcInfo%Name), ' at (I,J,L) = ', I, J, L
                   ErrorMsg = 'Negative species concentations encountered.'// &
                            ' This may be fixed by increasing the'        //  &
                            ' background concentration or by shortening'  //  &
                            ' the transport time step.'
                   RC = GC_FAILURE
                ENDIF
#endif
             ENDIF

          ENDDO !L
       ENDDO !J
       ENDDO !I
!$OMP END PARALLEL DO

       ! Exit with error condition
       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error( ErrorMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Nullify pointer
       SpcInfo  => NULL()

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

    ! Halt mixing timer (so that unit conv can be timed separately)
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End( "Boundary layer mixing", RC )
    ENDIF

    ! Convert State_Chm%Species back to original units
    CALL Convert_Spc_Units(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         mapping    = State_Chm%Map_Advect,                                  &
         new_units  = previous_units,                                        &
         RC         = RC                                                    )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
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

    !------------------------------------------------------------------------
    ! Emissions/dry deposition budget diagnostics - Part 2 of 2
    !------------------------------------------------------------------------
    IF ( State_Diag%Archive_BudgetEmisDryDep ) THEN

       ! Timestep for diagnostics [s]
       DT_Tend = DBLE( TS )

       ! Compute change in column masses (after emis/dryd - before emis/dryd)
       ! and store in diagnostic arrays.  Units are [kg/s].
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
            timeStep    = DT_Tend,                                           &
            RC          = RC                                                )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Emissions/dry deposition budget diagnostics error 2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Free pointers
    DepFreq => NULL()

  IF ( ASSOCIATED( PNOxLoss_O3 ) )   DEALLOCATE( PNOxLoss_O3 )
  IF ( ASSOCIATED( PNOxLoss_HNO3 ) ) DEALLOCATE( PNOxLoss_HNO3 )

  PNOxLoss_O3 => NULL()
  PNOxLoss_HNO3 => NULL()

  END SUBROUTINE DO_TEND

!EOC
END MODULE MIXING_MOD
