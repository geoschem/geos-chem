!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
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
MODULE MIXING_MOD
!
! !USES:
!
  USE PRECISION_MOD

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: INIT_MIXING
  PUBLIC :: DO_MIXING
  PUBLIC :: DO_TEND
!
! !REVISION HISTORY:
!  04 Mar 2015 - C. Keller   - Initial version.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_mixing
!
! !DESCRIPTION: Subroutine INIT\_MIXING initialized the pbl mixing wrapper
! module.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_MIXING ( Input_Opt,  State_Chm, State_Diag, &
                            State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,   ONLY : OptInput
    USE PBL_MIX_MOD,     ONLY : COMPUTE_PBL_HEIGHT
    USE PBL_MIX_MOD,     ONLY : DO_PBL_MIX
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Diag_Mod,  ONLY : DgnState
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Met_Mod,   ONLY : MetState
    USE VDIFF_MOD,       ONLY : DO_PBL_MIX_2
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt   ! Input Options
    TYPE(GrdState),   INTENT(IN   )  :: State_Grid  ! Grid State
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm   ! Chemistry State
    TYPE(DgnState),   INTENT(INOUT)  :: State_Diag  ! Diagnostics State
    TYPE(MetState),   INTENT(INOUT)  :: State_Met   ! Meteorology State
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
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=======================================================================
    ! INIT_MIXING begins here!
    !=======================================================================

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at INIT_MIXING (in module GeosCore/mixing_mod.F90)'

    !-----------------------------------------------------------------------
    ! Initialize PBL mixing scheme
    !-----------------------------------------------------------------------
    IF ( Input_Opt%LNLPBL ) THEN

       ! Initialize non-local PBL mixing scheme
       CALL DO_PBL_MIX_2( .FALSE.,    Input_Opt,  State_Chm, &
                          State_Diag, State_Grid, State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "DO_PBL_MIX" at initialization!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ELSE

       ! Initialize full PBL mixing scheme
       CALL DO_PBL_MIX( .FALSE.,    Input_Opt, State_Chm,  &
                        State_Diag, State_Grid, State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "DO_PBL_MIX" at initialization!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

#if !defined( ESMF_ ) && !defined( MODEL_WRF )
    !-----------------------------------------------------------------------
    ! Compute the various PBL quantities with the initial met fields.
    ! This is needed so that HEMCO won't be passed a zero PBL height
    ! (bmy, 10/26/16)
    !
    ! In ESMF mode this routine should not be called during the init
    ! stage: the required met quantities are not yet defined.
    ! (ckeller, 11/23/16)
    !
    ! In GC-WRF, which uses the same entry-point as the GEOS-5 GCM, the
    ! required met quantities are also not defined until GIGC_Chunk_Run,
    ! so also skip this here (hplin, 8/9/18)
    !-----------------------------------------------------------------------
    CALL COMPUTE_PBL_HEIGHT( State_Grid, State_Met, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "COMPUTE_PBL_HEIGHT" at initialization!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#endif

  END SUBROUTINE INIT_MIXING
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
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
  SUBROUTINE DO_MIXING( Input_Opt,  State_Chm, State_Diag, &
                        State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PBL_MIX_MOD,        ONLY : DO_PBL_MIX
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_MOd,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE VDIFF_MOD,          ONLY : DO_PBL_MIX_2
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
       IF ( State_Diag%Archive_DryDepMix .or.  &
            State_Diag%Archive_DryDep        ) THEN
          State_Diag%DryDepMix = 0.0_f4
       ENDIF

       ! Non-local mixing
       CALL DO_PBL_MIX_2( Input_Opt%LTURB, Input_Opt,  State_Chm, &
                          State_Diag,      State_Grid, State_Met, RC )

       ! Trap potential error
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

       ! Initialize the diagnostic array for the History Component. This
       ! will prevent leftover values from being carried over to this
       ! timestep. (For example, if on the last iteration, the PBL height
       ! was higher than it is now, then we will have stored drydep fluxes
       ! up to that height, so we need to zero these out.)
       IF ( State_Diag%Archive_DryDepMix .or.  &
            State_Diag%Archive_DryDep        ) THEN
          State_Diag%DryDepMix = 0.0_f4
       ENDIF

       ! Full PBL mixing
       CALL DO_PBL_MIX( Input_Opt%LTURB, Input_Opt,  State_Chm, &
                        State_Diag,      State_Grid, State_Met, RC )

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
!                  Harvard-NASA Emissions Component (HEMCO)                   !
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
    USE Diagnostics_Mod,    ONLY : Compute_Column_Mass
    USE Diagnostics_Mod,    ONLY : Compute_Budget_Diagnostics
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : SAFE_DIV
    USE GET_NDEP_MOD,       ONLY : SOIL_DRYDEP
    USE HCO_INTERFACE_MOD,  ONLY : GetHcoVal, GetHcoDiagn
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants,      ONLY : AVO
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_DYN, GET_TS_CONV, GET_TS_CHEM
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
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
    INTEGER                 :: I, J, L, L1, L2, N, D, NN, NA, nAdvect
    INTEGER                 :: DRYDEPID
    INTEGER                 :: PBL_TOP, DRYD_TOP, EMIS_TOP
    REAL(fp)                :: TS, TMP, FRQ, RKT, FRAC, FLUX, AREA_M2
    REAL(fp)                :: MWkg, DENOM
    LOGICAL                 :: FND
    LOGICAL                 :: PBL_DRYDEP, LSCHEM, ChemGridOnly
    LOGICAL                 :: LEMIS,      LDRYD
    LOGICAL                 :: DryDepSpec, EmisSpec
    CHARACTER(LEN=63)       :: OrigUnit
    CHARACTER(LEN=255)      :: MSG

    ! PARANOX loss fluxes (kg/m2/s). These are obtained from the
    ! HEMCO PARANOX extension via the diagnostics module.
    REAL(fp)                :: PNOXLOSS
    REAL(f4), POINTER, SAVE :: PNOXLOSS_O3  (:,:) => NULL()
    REAL(f4), POINTER, SAVE :: PNOXLOSS_HNO3(:,:) => NULL()

    ! SAVEd scalars (defined on first call only)
    LOGICAL,           SAVE :: FIRST = .TRUE.
    INTEGER,           SAVE :: id_MACR,  id_RCHO,  id_ACET, id_ALD2
    INTEGER,           SAVE :: id_ALK4,  id_C2H6,  id_C3H8, id_CH2O
    INTEGER,           SAVE :: id_PRPE,  id_O3,    id_HNO3, id_BrO
    INTEGER,           SAVE :: id_Br2,   id_Br,    id_HOBr, id_HBr
    INTEGER,           SAVE :: id_BrNO3

    ! Pointers and objects
    TYPE(Species), POINTER  :: SpcInfo
    REAL(fp),      POINTER  :: DEPSAV       (:,:,:  )  ! IM, JM, nDryDep

    ! Temporary save for total ch4 (Xueying Yu, 12/08/2017)
    LOGICAL                 :: ITS_A_CH4_SIM
    REAL(fp)                :: total_ch4_pre_soil_absorp(State_Grid%NX, &
                                                         State_Grid%NY, &
                                                         State_Grid%NZ)

    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! DO_TEND begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ErrMsg  = ''
    ThisLoc = ' -> at DO_TEND (in module GeosCore/mixing_mod.F90)'

    ! Special case that there is no dry deposition and emissions
    IF ( .NOT. Input_Opt%LDRYD .AND. .NOT. Input_Opt%LEMIS ) RETURN

    ! Initialize
    LSCHEM            = Input_Opt%LSCHEM
    LEMIS             = Input_Opt%LEMIS
    LDRYD             = Input_Opt%LDRYD
    PBL_DRYDEP        = Input_Opt%PBL_DRYDEP
    ITS_A_CH4_SIM     = Input_Opt%ITS_A_CH4_SIM
    nAdvect           = State_Chm%nAdvect

    ! Initialize pointer
    SpcInfo           => NULL()
    DEPSAV            => State_Chm%DryDepSav

    !----------------------------------------------------------
    ! Emissions/dry deposition budget diagnostics - Part 1 of 2
    !----------------------------------------------------------
    IF ( State_Diag%Archive_BudgetEmisDryDep ) THEN
       ! Get initial column masses
       CALL Compute_Column_Mass( Input_Opt,                               &
                                 State_Chm, State_Grid, State_Met,        &
                                 State_Chm%Map_Advect,                    &
                                 State_Diag%Archive_BudgetEmisDryDepFull, &
                                 State_Diag%Archive_BudgetEmisDryDepTrop, &
                                 State_Diag%Archive_BudgetEmisDryDepPBL,  &
                                 State_Diag%BudgetMass1,                  &
                                 RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Emissions/dry deposition budget diagnostics error 1'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! DO_TEND previously operated in units of kg. The species arrays are in
    ! v/v for mixing, hence needed to convert before and after.
    ! Now use units kg/m2 as State_Chm%SPECIES units in DO_TEND to
    ! remove area-dependency (ewl, 9/30/15)
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'kg/m2', RC, OrigUnit=OrigUnit )

    ! Trap potential error
    IF ( RC /= GC_SUCCESS ) THEN
       MSG = 'Unit conversion error!'
       CALL GC_Error( MSG, RC, 'DO_TEND in mixing_mod.F90' )
       RETURN
    ENDIF

    ! Get time step [s]
    IF ( PRESENT(DT) ) THEN
       TS = DT
    ELSE
       TS = GET_TS_DYN()
    ENDIF

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
       id_PRPE = Ind_('PRPE' )
       id_O3   = Ind_('O3'   )
       id_HNO3 = Ind_('HNO3' )
       id_BrO  = Ind_('BrO'  )
       id_Br2  = Ind_('Br2'  )
       id_Br   = Ind_('Br'   )
       id_HOBr = Ind_('HOBr' )
       id_HBr  = Ind_('HBr'  )
       id_BrNO3= Ind_('BrNO3')

       ! On first call, get pointers to the PARANOX loss fluxes. These are
       ! stored in diagnostics 'PARANOX_O3_DEPOSITION_FLUX' and
       ! 'PARANOX_HNO3_DEPOSITION_FLUX'. The call below links pointers
       ! PNOXLOSS_O3 and PNOXLOSS_HNO3 to the data values stored in the
       ! respective diagnostics. The pointers will remain unassociated if
       ! the diagnostics do not exist.
       ! This is only needed if non-local PBL scheme is not being used.
       ! Otherwise, PARANOX fluxes are applied in vdiff_mod.F90.
       !  (ckeller, 4/10/2015)
       IF ( .NOT. Input_Opt%LNLPBL ) THEN
          CALL GetHcoDiagn( 'PARANOX_O3_DEPOSITION_FLUX'  , &
                            .FALSE.,   RC, Ptr2D = PNOXLOSS_O3          )
          CALL GetHcoDiagn( 'PARANOX_HNO3_DEPOSITION_FLUX', &
                            .FALSE.,   RC, Ptr2D = PNOXLOSS_HNO3        )
       ENDIF
       FIRST = .FALSE.
    ENDIF

    !-----------------------------------------------------------------------
    ! For tagged CH4 simulations
    ! Save the total CH4 concentration before apply soil absorption
    !-----------------------------------------------------------------------
    IF ( ITS_A_CH4_SIM .and. Input_Opt%LSPLIT ) THEN
       total_ch4_pre_soil_absorp = State_Chm%Species(:,:,:,1)
    ENDIF

    !=======================================================================
    ! Do for every advected species and grid box
    !=======================================================================
!$OMP PARALLEL DO                                                           &
!$OMP DEFAULT( SHARED                                                     ) &
!$OMP PRIVATE( I,        J,            L,          L1,       L2           ) &
!$OMP PRIVATE( N,        PBL_TOP,      FND,        TMP,      DryDepId     ) &
!$OMP PRIVATE( FRQ,      RKT,          FRAC,       FLUX,     Area_m2      ) &
!$OMP PRIVATE( MWkg,     ChemGridOnly, DryDepSpec, EmisSpec, DRYD_TOP     ) &
!$OMP PRIVATE( EMIS_TOP, PNOXLOSS,     DENOM,      SpcInfo,  NA           )
    DO NA = 1, nAdvect

       ! Get the species ID from the advected species ID
       N = State_Chm%Map_Advect(NA)

       ! Get info about this species from the species database
       SpcInfo => State_Chm%SpcData(N)%Info

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
             CALL GetHcoVal ( N, 1, 1, 1, DryDepSpec, dep = TMP )
          ENDIF

          ! Special case for O3 or HNO3: include PARANOX loss
          IF ( N == id_O3   .AND. ASSOCIATED(PNOXLOSS_O3  ) )    &
               DryDepSpec = .TRUE.
          IF ( N == id_HNO3 .AND. ASSOCIATED(PNOXLOSS_HNO3) )    &
               DryDepSpec = .TRUE.
       ENDIF

       !--------------------------------------------------------------------
       ! Check if we need to do emissions for this species
       !--------------------------------------------------------------------
       IF ( LEMIS ) THEN
          CALL GetHcoVal ( N, 1, 1, 1, EmisSpec, emis = TMP )
       ELSE
          EmisSpec = .FALSE.
       ENDIF

       !--------------------------------------------------------------------
       ! Can go to next species if this species does not have
       ! dry deposition and/or emissions
       !--------------------------------------------------------------------
       IF ( .NOT. DryDepSpec .AND. .NOT. EmisSpec ) CYCLE

       ! Loop over all grid boxes
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          !-----------------------------------------------------------------
          ! Define various quantities before computing tendencies
          !-----------------------------------------------------------------

          ! Get PBL_TOP at this grid box
          PBL_TOP = MAX( 1, FLOOR( State_Met%PBL_TOP_L(I,J) ) )

          ! Molecular weight in kg
          MWkg = SpcInfo%emMW_g * 1.e-3_fp

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
          IF ( LSCHEM ) THEN
             IF ( N == id_BrO  .OR. N == id_Br2   .OR. &
                  N == id_Br   .OR. N == id_HOBr  .OR. &
                  N == id_HBr  .OR. N == id_BrNO3       ) THEN
                ChemGridOnly = .TRUE.
             ENDIF
          ENDIF

          ! For non-UCX runs, never emit above the chemistry grid.
          ! (ckeller, 6/18/15)
          ! Exclude all specialty simulations (ewl, 3/17/16)
          IF ( Input_Opt%ITS_A_FULLCHEM_SIM .AND.  &
               .NOT. Input_Opt%LUCX ) THEN
             ChemGridOnly = .TRUE.
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
                ! stored in DEPSAV. Units are [s-1].
                IF ( DRYDEPID > 0 ) THEN
                   FRQ = DEPSAV(I,J,DRYDEPID)
                ENDIF

                ! Dry deposition frequency from HEMCO. HEMCO calculates
                ! dry deposition frequencies for air-sea exchange and
                ! from ship NOx plume parameterization (PARANOx). The
                ! units are [s-1].
                CALL GetHcoVal ( N, I, J, 1, FND, dep=TMP )

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
                   FLUX = ( 1.0_fp - FRAC ) * State_Chm%Species(I,J,L,N)

                   ! Apply dry deposition
                   State_Chm%Species(I,J,L,N) = FRAC *    &
                                                State_Chm%Species(I,J,L,N)

                   ! Eventually add PARANOX loss. PNOXLOSS is in kg/m2/s.
                   ! Make sure PARANOx loss is applied to tracers. (ckeller,
                   ! 3/29/16).
                   IF ( PNOXLOSS > 0 ) THEN
                      State_Chm%Species(I,J,L,N) = &
                         State_Chm%Species(I,J,L,N) - ( PNOXLOSS * TS )
                      FLUX = FLUX + ( PNOXLOSS * TS )
                   ENDIF

!                   IF (AREA_M2 .eq. 0.0_fp) THEN
!                     PRINT*, "FLUX: ", FLUX
!                     PRINT*, "MWkg: ", MWkg
!                     PRINT*, "AVO: ", AVO
!                     PRINT*, "TS: ", TS
!                     PRINT*, "AREA_M2: ", AREA_M2
!                     CALL FLUSH(6)
!                   ENDIF

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
                      CALL SOIL_DRYDEP( I, J, L, N, FLUX, State_Chm )
                   ENDIF

                   !--------------------------------------------------------
                   ! HISTORY: Archive drydep flux loss from mixing
                   ! Units = molec/cm2/s
                   !
                   ! NOTE: we don't need to multiply by the ratio of
                   ! TS_CONV / TS_CHEM, as the updating frequency for
                   ! HISTORY is determined by the "frequency" setting in
                   ! the "HISTORY.rc"input file.  The old bpch diagnostics
                   ! archived the drydep due to chemistry every chemistry
                   ! timestep = 2X the dynamic timestep.  So in order to
                   ! avoid double-counting the drydep flux from mixing,
                   ! you had to multiply by TS_CONV / TS_CHEM.
                   !
                   ! ALSO NOTE: When comparing History output to bpch
                   ! output, you must use an updating frequency equal to
                   ! the dynamic timestep so that the drydep fluxes due to
                   ! mixing will be equivalent w/ the bpch output.  It is
                   ! also recommended to turn off chemistry so as to be
                   ! able to compare the drydep fluxes due to mixing in
                   ! bpch vs. History as an "apples-to-apples" comparison.
                   !
                   !    -- Bob Yantosca (yantosca@seas.harvard.edu)
                   !--------------------------------------------------------
                   IF ( ( State_Diag%Archive_DryDepMix .or.        &
                          State_Diag%Archive_DryDep        ) .and. &
                          DryDepID > 0 ) THEN
                      State_Diag%DryDepMix(I,J,DryDepId) = Flux
                   ENDIF

                ENDIF ! apply drydep
             ENDIF ! L <= PBLTOP

             !--------------------------------------------------------------
             ! Apply emissions.
             ! These are always taken from HEMCO
             !--------------------------------------------------------------
             IF ( EmisSpec .AND. ( L <= EMIS_TOP ) ) THEN

                ! Get HEMCO emissions. Units are [kg/m2/s].
                CALL GetHcoVal ( N, I, J, L, FND, emis=TMP )

                ! Add emissions (if any)
                ! Bug fix: allow negative fluxes. (ckeller, 4/12/17)
                !IF ( FND .AND. (TMP > 0.0_fp) ) THEN
                IF ( FND ) THEN

                   ! Flux: [kg/m2] = [kg m-2 s-1 ] x [s]
                   FLUX = TMP * TS

                   ! Add to species array
                   State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N) &
                                              + FLUX
                ENDIF
             ENDIF

             !--------------------------------------------------------------
             ! Special handling for tagged CH4 simulations
             !
             ! Tagged CH4 species are split off into a separate loop to
             ! ensure we remove soil absorption from NA=1 (total CH4) first
             !--------------------------------------------------------------
             IF ( ITS_A_CH4_SIM .and. Input_Opt%LSPLIT ) THEN

                ! If we are in the chemistry grid
                IF ( L <= EMIS_TOP ) THEN

                   ! Total CH4 species
                   IF ( NA == 1 ) THEN

                      ! Get soil absorption from HEMCO. Units are [kg/m2/s].
                      ! CH4_SAB is species #15
                      CALL GetHcoVal ( 15, I, J, L, FND, emis=TMP )

                      ! Remove soil absorption from total CH4 emissions
                      IF ( FND ) THEN

                         ! Flux: [kg/m2] = [kg m-2 s-1 ] x [s]
                         FLUX = TMP * TS

                         ! Apply soil absorption as loss
                         State_Chm%Species(I,J,L,N) =                       &
                         State_Chm%Species(I,J,L,N) - FLUX
                      ENDIF

                   ENDIF

                ENDIF

             ENDIF

             ! Prevent negative concentrations. (ckeller, 3/29/16)
             State_Chm%Species(I,J,L,N) = MAX(State_Chm%Species(I,J,L,N),0.0_fp)

          ENDDO !L
       ENDDO !J
       ENDDO !I

       ! Nullify pointer
       SpcInfo  => NULL()

    ENDDO !N
!$OMP END PARALLEL DO

    !--------------------------------------------------------------
    ! Special handling for tagged CH4 simulations
    !--------------------------------------------------------------
    IF ( ITS_A_CH4_SIM .and. Input_Opt%LSPLIT ) THEN

!$OMP PARALLEL DO               &
!$OMP DEFAULT( SHARED         ) &
!$OMP PRIVATE( I, J, L, N, NA )
       DO NA = 1, nAdvect

          ! Get the species ID from the advected species ID
          N = State_Chm%Map_Advect(NA)

          ! Loop over all grid boxes
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Tagged CH4 tracers
             IF ( NA >= 2 .and. NA <= nAdvect-1 ) THEN

                ! Apply soil absorption to each tagged CH4 species
                ! (Xueying Yu, 12/08/2017)
                State_Chm%Species(I,J,L,N) =                    &
                     SAFE_DIV(State_Chm%Species(I,J,L,N),       &
                              total_ch4_pre_soil_absorp(I,J,L), &
                              0.e+0_fp) *                       &
                     State_Chm%Species(I,J,L,1)

             ENDIF

             ! Prevent negative concentrations. (ckeller, 3/29/16)
             State_Chm%Species(I,J,L,N) = MAX(State_Chm%Species(I,J,L,N),0.0_fp)

          ENDDO
          ENDDO
          ENDDO
       ENDDO
!$OMP END PARALLEL DO

    ENDIF

    ! Convert State_Chm%Species back to original units
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       MSG = 'Unit conversion error!'
       CALL GC_Error( MSG, RC, 'DO_TEND in mixing_mod.F90' )
       RETURN
    ENDIF

    !----------------------------------------------------------
    ! Emissions/dry deposition budget diagnostics - Part 2 of 2
    !----------------------------------------------------------
    IF ( State_Diag%Archive_BudgetEmisDryDep ) THEN
       ! Get final column masses and compute diagnostics
       CALL Compute_Column_Mass( Input_Opt,                                   &
                                 State_Chm, State_Grid, State_Met,            &
                                 State_Chm%Map_Advect,                        &
                                 State_Diag%Archive_BudgetEmisDryDepFull,     &
                                 State_Diag%Archive_BudgetEmisDryDepTrop,     &
                                 State_Diag%Archive_BudgetEmisDryDepPBL,      &
                                 State_Diag%BudgetMass2,                      &
                                 RC )
       CALL Compute_Budget_Diagnostics( State_Grid,                           &
                                     State_Chm%Map_Advect,                    &
                                     TS,                                      &
                                     State_Diag%Archive_BudgetEmisDryDepFull, &
                                     State_Diag%Archive_BudgetEmisDryDepTrop, &
                                     State_Diag%Archive_BudgetEmisDryDepPBL,  &
                                     State_Diag%BudgetEmisDryDepFull,         &
                                     State_Diag%BudgetEmisDryDepTrop,         &
                                     State_Diag%BudgetEmisDryDepPBL,          &
                                     State_Diag%BudgetMass1,                  &
                                     State_Diag%BudgetMass2,                  &
                                     RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Emissions/dry deposition budget diagnostics error 2'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  ! Nullify pointers
  NULLIFY( DEPSAV )

  END SUBROUTINE DO_TEND
!EOC
END MODULE MIXING_MOD
