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
! !PRIVATE MEMBER FUNCTIONS:
!
! !REVISION HISTORY:
!  04 Mar 2015 - C. Keller   - Initial version. 
!EOP
!------------------------------------------------------------------------------
!BOC
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
  SUBROUTINE INIT_MIXING ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PBL_MIX_MOD,        ONLY : COMPUTE_PBL_HEIGHT
    USE PBL_MIX_MOD,        ONLY : DO_PBL_MIX
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE VDIFF_MOD,          ONLY : DO_PBL_MIX_2
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input opts
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState),   INTENT(INOUT)  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm  ! Chemistry state 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS
!  (A) While all dry deposition rates are calculated either in
!      DO_PBL_MIX2 or DO_TEND, settling of aerosols is still 
!      computed in the dust/seasalt modules.
!
! !REVISION HISTORY: 
!  04 Mar 2015 - C. Keller   - Initial version
!  26 Oct 2016 - R. Yantosca - Now also call COMPUTE_PBL_HEIGHT so that we
!                              populate PBL quantities w/ the initial met
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! DO_MIXING begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS
   
    ! ------------------------------------------------------------------
    ! Initialize PBL mixing scheme
    ! ------------------------------------------------------------------
    IF ( Input_Opt%LNLPBL ) THEN
       CALL DO_PBL_MIX_2( am_I_Root, .FALSE. ,  Input_Opt, &
                          State_Met, State_Chm, RC          )
    ELSE 
       CALL DO_PBL_MIX( am_I_Root, .FALSE.,   Input_Opt, &
                        State_Met, State_Chm, RC          )
    ENDIF

    ! Compute the various PBL quantities with the initial met fields.
    ! This is needed so that HEMCO won't be passed a zero PBL height
    ! (bmy, 10/26/16)
    CALL COMPUTE_PBL_HEIGHT( State_Met )

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
  SUBROUTINE DO_MIXING ( am_I_Root, Input_Opt, State_Met, State_Chm, RC ) 
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : GC_Error
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PBL_MIX_MOD,        ONLY : DO_PBL_MIX
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE UnitConv_Mod
    USE VDIFF_MOD,          ONLY : DO_PBL_MIX_2
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )  :: am_I_Root  ! root CPU?
    TYPE(OptInput),   INTENT(IN   )  :: Input_Opt  ! Input opts
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState),   INTENT(INOUT)  :: State_Met  ! Met state
    TYPE(ChmState),   INTENT(INOUT)  :: State_Chm  ! Chemistry state 
    INTEGER,          INTENT(INOUT)  :: RC         ! Failure or success
!
! !REMARKS
!  (A) While all dry deposition rates are calculated either in
!      DO_PBL_MIX2 or DO_TEND, settling of aerosols is still 
!      computed in the dust/seasalt modules.
!
! !REVISION HISTORY: 
!  04 Mar 2015 - C. Keller   - Initial version 
!  12 Aug 2015 - E. Lundgren - Input tracer units are now [kg/kg] and 
!                              are converted to [v/v] for mixing
!  30 Sep 2014 - E. Lundgren - Move unit conversion for DO_TEND to DO_TEND
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  08 Aug 2016 - R. Yantosca - Remove temporary tracer-removal code
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL :: OnlyAbovePBL

    !=================================================================
    ! DO_MIXING begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Convert [kg/kg dry air] to [v/v dry air] for mixing (ewl, 8/12/15)
     CALL ConvertSpc_KgKgDry_to_VVDry( am_I_Root, State_Chm, RC )  

    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'DO_MIXING in mixing_mod.F')
       RETURN
    ENDIF  

    ! ------------------------------------------------------------------
    ! Do non-local PBL mixing. This will apply the species tendencies
    ! (emission fluxes and dry deposition rates) below the PBL.
    ! This is done for all species with defined emissions / dry
    ! deposition rates, including dust.
    ! Set OnlyAbovePBL flag (used below by DO_TEND) to indicate that
    ! fluxes within the PBL have already been applied. 
    ! ------------------------------------------------------------------
    IF ( Input_Opt%LTURB .AND. Input_Opt%LNLPBL ) THEN
       CALL DO_PBL_MIX_2( am_I_Root, Input_Opt%LTURB, Input_Opt, &
                          State_Met, State_Chm,       RC          )
       OnlyAbovePBL = .TRUE.
    ELSE
       OnlyAbovePBL = .FALSE.
    ENDIF

    ! ------------------------------------------------------------------
    ! Apply tendencies. This will apply dry deposition rates and 
    ! emission fluxes below the PBL if it has not yet been done
    ! via the non-local PBL mixing. It also adds the emissions above 
    ! the PBL to the species array. Emissions of some species may be 
    ! capped at the tropopause to avoid build-up in stratosphere.
    ! ------------------------------------------------------------------

    ! Apply tendencies
    CALL DO_TEND ( am_I_Root, Input_Opt,    State_Met, &
                   State_Chm, OnlyAbovePBL, RC          )

    ! ------------------------------------------------------------------
    ! Do full pbl mixing. This fully mixes the updated species 
    ! concentrations within the PBL. 
    ! 
    ! Now also archive concentrations and calculate turbulence 
    ! tendencies (ckeller, 7/15/2015)
    ! ------------------------------------------------------------------
    IF ( Input_Opt%LTURB .AND. .NOT. Input_Opt%LNLPBL ) THEN
       CALL DO_PBL_MIX( am_I_Root, Input_Opt%LTURB, Input_Opt, & 
                        State_Met, State_Chm,       RC )
    ENDIF

    ! Convert species conc back to [kg/kg dry air] after mixing (ewl, 8/12/15)
    CALL ConvertSpc_VVDry_to_KgKgDry( am_I_Root, State_Chm, RC )

    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'DO_MIXING in mixing_mod.F')
       RETURN
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
  SUBROUTINE DO_TEND ( am_I_Root, Input_Opt,    State_Met, &
                       State_Chm, OnlyAbovePBL, RC,    DT   ) 
!
! !USES:
!
    USE CHEMGRID_MOD,       ONLY : GET_CHEMGRID_LEVEL
    USE CMN_DIAG_MOD,       ONLY : ND44
    USE CMN_SIZE_MOD,       ONLY : IIPAR,   JJPAR,   LLPAR
    USE DRYDEP_MOD,         ONLY : DEPSAV
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : ERROR_STOP, SAFE_DIV
    USE GET_NDEP_MOD,       ONLY : SOIL_DRYDEP
    USE HCO_INTERFACE_MOD,  ONLY : GetHcoVal, GetHcoDiagn
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
    USE PhysConstants,      ONLY : AVO
    USE Species_Mod,        ONLY : Species
    USE State_Met_Mod,      ONLY : MetState
    USE State_Chm_Mod,      ONLY : ChmState
    USE TIME_MOD,           ONLY : GET_TS_DYN, GET_TS_CONV, GET_TS_CHEM
    USE State_Chm_Mod,      ONLY : Ind_
    USE UnitConv_Mod
#if defined( BPCH_DIAG )
    USE DIAG_MOD,           ONLY : AD44
#endif
#if defined( NC_DIAG )
    USE HCO_ERROR_MOD
    USE HCO_INTERFACE_MOD,  ONLY : HcoState, GetHcoID
    USE HCO_DIAGN_MOD,      ONLY : Diagn_Update
#endif
#if defined( USE_TEND )
    USE TENDENCIES_MOD
#endif
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN   )           :: am_I_Root    ! root CPU?
    TYPE(OptInput),   INTENT(IN   )           :: Input_Opt    ! Input opts
    TYPE(MetState),   INTENT(IN   )           :: State_Met    ! Met state
    LOGICAL,          INTENT(IN   )           :: OnlyAbovePBL ! Only above PBL?
    REAL(fp),         INTENT(IN   ), OPTIONAL :: DT           ! Time step [s]
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT)           :: State_Chm    ! Chemistry state 
    INTEGER,          INTENT(INOUT)           :: RC           ! Failure or success
!
! !REVISION HISTORY: 
!  04 Mar 2015 - C. Keller   - Initial version 
!  09 Mar 2015 - R. Yantosca - Bug fix: Use the drydep ID number instead of the
!                              tracer number to index the AD44 drydep array
!  09 Mar 2015 - R. Yantosca - Bug fix: Remove an IF ( L1 == 1 ) block where
!                              we define DRYDEPID.  This isn't needed here.
!  10 Apr 2015 - C. Keller   - Now exchange PARANOX loss fluxes via HEMCO 
!                              diagnostics.
!  12 Jun 2015 - R. Yantosca - Bug fix in SAFE_DIV: the denominator was
!                              arranged wrongly.  Now corrected.
!  18 Jun 2015 - C. Keller   - Now restrict all emissions to chemistry grid
!                              if UCX=false.
!  30 Sep 2015 - E. Lundgren - Now convert locally to kg/m2 for area-independent
!                              compatibility between tracer units and flux
!  22 Mar 2016 - C. Keller   - Bug fix: make sure drydep velocities are written
!                              to diagnostics if emissions are zero.
!  16 Mar 2016 - E. Lundgren - Exclude specialty simulations in restriction of
!                              all emissions to chemistry grid if UCX=false
!  29 Feb 2016 - C. Keller   - Make sure PARANOx fluxes are applied to tracers.
!  29 Apr 2016 - R. Yantosca - Don't initialize pointers in declaration stmts
!  25 May 2016 - E. Lundgren - Replace input_opt%TRACER_MW_KG with species
!                              database field emMW_g (emitted species g/mol)
!  16 Jun 2016 - C. Miller   - Now define species ID's with the Ind_ function
!  17 Jun 2016 - R. Yantosca - Only define species ID's on the first call
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!  01 Jul 2016 - R. Yantosca - Now rename species DB object ThisSpc to SpcInfo
!  19 Jul 2016 - R. Yantosca - Now bracket tendency calls with #ifdef USE_TEND
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

#if defined( NC_DIAG )
    ! For netCDF diagnostics
    INTEGER                 :: cID, spc_id, HCRC
    CHARACTER(LEN=30)       :: DiagnName
    REAL(fp), TARGET        :: DryDepFlux( IIPAR, JJPAR, State_Chm%nDryDep ) 
    REAL(fp), TARGET        :: EMIS( IIPAR, JJPAR, LLPAR, State_Chm%nAdvect) 
    REAL(fp)                :: TOTFLUX(State_Chm%nAdvect)
    REAL(fp), POINTER       :: Ptr2D(:,:  )
    REAL(fp), POINTER       :: Ptr3D(:,:,:)
#endif

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

    !=================================================================
    ! DO_TEND begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Special case that there is no dry deposition and emissions
    IF ( .NOT. Input_Opt%LDRYD .AND. .NOT. Input_Opt%LEMIS ) RETURN

    ! Initialize
    LSCHEM     = Input_Opt%LSCHEM
    LEMIS      = Input_Opt%LEMIS 
    LDRYD      = Input_Opt%LDRYD 
    PBL_DRYDEP = Input_Opt%PBL_DRYDEP
    nAdvect    = State_Chm%nAdvect

    ! Initialize pointer
    SpcInfo    => NULL()

    ! DO_TEND previously operated in units of kg. The species arrays are in
    ! v/v for mixing, hence needed to convert before and after.
    ! Now use units kg/m2 as State_Chm%SPECIES units in DO_TEND to 
    ! remove area-dependency (ewl, 9/30/15)
    ! v/v --> kg/m2
    CALL ConvertSpc_VVDry_to_KgKgDry( am_I_Root,            State_Chm, RC )
    CALL ConvertSpc_KgKgDry_to_Kgm2 ( am_I_Root, State_Met, State_Chm, RC )

    ! Get time step [s]
    IF ( PRESENT(DT) ) THEN
       TS = DT
    ELSE
       TS = GET_TS_DYN() * 60.0_fp
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
       ! Otherwise, PARANOX fluxes are applied in vdiff_mod.F.
       !  (ckeller, 4/10/2015) 
       IF ( .NOT. Input_Opt%LNLPBL ) THEN
          CALL GetHcoDiagn( am_I_Root, 'PARANOX_O3_DEPOSITION_FLUX'  , &
                            .FALSE.,   RC, Ptr2D = PNOXLOSS_O3          ) 
          CALL GetHcoDiagn( am_I_Root, 'PARANOX_HNO3_DEPOSITION_FLUX', &
                            .FALSE.,   RC, Ptr2D = PNOXLOSS_HNO3        ) 
       ENDIF
       FIRST = .FALSE.
    ENDIF

#if defined( NC_DIAG )
    ! Initialize local diagnostic variables
    DryDepFlux  = 0.0_fp
    EMIS        = 0.0_fp
    TOTFLUX     = 0.0_fp
#endif

#if defined( USE_TEND )
    ! Archive concentrations for tendencies (ckeller, 7/15/2015) 
      CALL TEND_STAGE1( am_I_Root, Input_Opt, State_Met, &
                        State_Chm, 'FLUX', RC )
#endif

    ! Do for every advected species and grid box
!$OMP PARALLEL DO                                                    &
!$OMP DEFAULT ( SHARED )                                             &
!$OMP PRIVATE( I, J, L, L1, L2, N, D, PBL_TOP, FND, TMP, DRYDEPID  ) &
!$OMP PRIVATE( FRQ, RKT, FRAC, FLUX, AREA_M2,   MWkg, ChemGridOnly ) & 
!$OMP PRIVATE( DryDepSpec, EmisSpec, DRYD_TOP,  EMIS_TOP, PNOXLOSS ) &
!$OMP PRIVATE( DENOM, SpcInfo, NA                                  )
    DO NA = 1, nAdvect

       ! Get the species ID from the advected species ID
       N = State_Chm%Map_Advect(NA)

       ! Get info about this species from the species database
       SpcInfo => State_Chm%SpcData(N)%Info

       !----------------------------------------------------------------
       ! Check if we need to do dry deposition for this species 
       !----------------------------------------------------------------

       ! Initialize
       DryDepSpec = .FALSE.
       DryDepID   = -1

       ! Only if dry deposition is turned on and we do want to consider
       ! processes below the PBL...
       IF ( LDRYD .AND. .NOT. OnlyAbovePBL ) THEN

          ! Get dry deposition ID (used by drydep_mod.F) for this species.
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

       !----------------------------------------------------------------
       ! Check if we need to do emissions for this species 
       !----------------------------------------------------------------
       IF ( LEMIS ) THEN
          CALL GetHcoVal ( N, 1, 1, 1, EmisSpec, emis = TMP )
       ELSE
          EmisSpec = .FALSE.
       ENDIF

       !----------------------------------------------------------------
       ! Can go to next species if this species does not have 
       ! dry deposition and/or emissions
       !----------------------------------------------------------------
       IF ( .NOT. DryDepSpec .AND. .NOT. EmisSpec ) CYCLE

       ! Loop over all grid boxes
       DO J = 1, JJPAR    
       DO I = 1, IIPAR    

          !------------------------------------------------------------
          ! Define various quantities before computing tendencies
          !------------------------------------------------------------

          ! Get PBL_TOP at this grid box
          PBL_TOP = State_Met%PBL_TOP_L(I,J)

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

          ! Set emissions to zero above chemistry grid for the following 
          ! VOCs (adopted from aeic_mod.F).
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
             EMIS_TOP = GET_CHEMGRID_LEVEL( I, J, State_Met )
             EMIS_TOP = MIN(LLPAR,EMIS_TOP)
          ELSE
             EMIS_TOP = LLPAR
          ENDIF

          ! L2 is the upper level index to loop over
          L2 = MAX(DRYD_TOP, EMIS_TOP)

          ! This should not happen:
          IF ( L2 < L1 ) CYCLE

          ! Loop over selected vertical levels 
          DO L = L1, L2

             !----------------------------------------------------------
             ! Apply dry deposition frequencies to all levels below the
             ! PBL top.
             !----------------------------------------------------------
             IF ( DryDepSpec .AND. ( L <= DRYD_TOP ) ) THEN

                ! Init
                FRQ = 0.0_fp
 
                ! Dry deposition frequency from drydep_mod.F. This is 
                ! stored in DEPSAV. Units are [s-1].
                IF ( DRYDEPID > 0 ) THEN
                   FRQ = DEPSAV(I,J,DRYDEPID)
                ENDIF

                ! Dry deposition frequency from HEMCO. HEMCO calculates
                ! dry deposition frequencies for air-sea exchange and 
                ! from ship NOx plume parameterization (PARANOx). The
                ! units are [s-1].
                CALL GetHcoVal ( N, I, J, 1, FND, dep=TMP )

                ! Add to dry dep frequency from drydep_mod.F
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

#if defined( NC_DIAG )
                   ! Prior to 1/22/16, archive deposition flux in kg/m2/s:
                   !DepFluxes(I,J,N) = DepFluxes(I,J,N) + ( FLUX / TS )
#endif
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
                      CALL SOIL_DRYDEP( I, J, L, N, FLUX )
                   ENDIF

                   !========================================================
                   ! ND44: Dry deposition diagnostic [molec/cm2/s]
                   !========================================================
                   IF ( ND44 > 0 .and. DryDepID > 0 ) THEN
#if defined( BPCH_DIAG )
                      ! For bpch diagnostic, store data in global AD44 array
                      AD44(I,J,DryDepID,1) = AD44(I,J,DryDepID,1) + FLUX &
                                             * GET_TS_CONV() / GET_TS_CHEM() 
#endif
#if defined( NC_DIAG )
                      ! For netcdf diagnostic, store data in local array
                      ! Now use same units as bpch for comparison (ewl, 1/22/16)
                      DryDepFlux(I,J,DryDepID) = DryDepFlux(I,J,DryDepID) &
                                                + FLUX
#endif
                   ENDIF

                ENDIF ! apply drydep
             ENDIF ! L <= PBLTOP

             !----------------------------------------------------------
             ! Apply emissions.
             ! These are always taken from HEMCO
             !----------------------------------------------------------
             IF ( EmisSpec .AND. ( L <= EMIS_TOP ) ) THEN
    
                ! Get HEMCO emissions. Units are [kg/m2/s].
                CALL GetHcoVal ( N, I, J, L, FND, emis=TMP )
           
                ! Add emissions (if any)
                IF ( FND .AND. (TMP > 0.0_fp) ) THEN

                   ! Flux: [kg/m2] = [kg m-2 s-1 ] x [s]
                   FLUX = TMP * TS

                   ! Add to species array
                   State_Chm%Species(I,J,L,N) = State_Chm%Species(I,J,L,N) & 
                                              + FLUX 

#if defined( NC_DIAG )
                   ! Update new species emissions diagnostics
                   EMIS(I,J,L,N) = TMP                 ! kg/m2/s
                   TOTFLUX(N)    = TOTFLUX(N) + FLUX   ! kg/m2/s
#endif
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

#if defined( USE_TEND )
      ! Calculate tendencies and write to diagnostics (ckeller, 7/15/2015)
      CALL TEND_STAGE2( am_I_Root, Input_Opt, State_Met, &
                        State_Chm, 'FLUX', TS, RC )
#endif


#if defined( NC_DIAG )
    !========================================================      
    ! ND44: Dry deposition diagnostic [molec/cm2/s] (netcdf) 
    !========================================================
    IF ( Input_Opt%ND44 > 0 ) THEN

       ! Assume success
       HCRC = HCO_SUCCESS

       ! Loop over all drydep species
       DO D = 1, State_Chm%nDryDep

          ! Get the species ID from the drydep ID
          spc_id = State_Chm%Map_DryDep(D)
      
          ! If this species number is scheduled for output in input.geos, 
          ! then archive the latest depvel data into the diagnostic structure
          IF ( ANY( Input_Opt%TINDEX(44,:) == spc_id ) ) THEN
       
             ! Define diagnostic name
             DiagnName = 'DRYDEP_FLX_MIX_' // &
                         TRIM( State_Chm%SpcData(spc_id)%Info%Name )

             ! Point to data
             Ptr2D => DryDepFlux(:,:,D)

             ! Update diagnostic container
             CALL Diagn_Update( am_I_Root, HcoState,                 &
                                cName   = TRIM( DiagnName),          &
                                Array2D = Ptr2D,                     &
                                COL     = Input_Opt%DIAG_COLLECTION, &
                                RC      = HCRC                          )
       
             ! Free the pointer
             Ptr2D => NULL()
       
             ! Stop with error if the diagnostics were unsuccessful
             IF ( HCRC /= HCO_SUCCESS ) THEN
                CALL ERROR_STOP( 'Cannot update drydep flux' //       &
                                 ' diagnostic: ' // TRIM( DiagnName), &
                                 'DO_TEND (mixing_mod.F)')
             ENDIF
          ENDIF
       ENDDO
    ENDIF

    ! For now, always output species emissions diagnostic since there
    ! is no logical switch for it in Input_Mod, nor an entry in input.geos
    DO NA = 1, nAdvect

       ! Get the species ID from the advected species count
       N = State_Chm%Map_Advect(NA)
 
       ! Skip if there are no emissions
       IF ( TOTFLUX(N) == 0.0_fp ) CYCLE

       ! Only if HEMCO species is defined
       ! NOTE: For other netcdf diagnostics we look at Input_Opt%TINDEX
       ! to see which speciess are turned on for the diagnostic
       cID = GetHcoID( SpcID=N )
       IF ( cID > 0 ) THEN 

          ! Entry in the species database
          SpcInfo => State_Chm%SpcData(N)%Info

          ! Define diagnostics name
          DiagnName = 'SPECIES_EMIS_' // TRIM( SpcInfo%Name )

          ! Point to species slice 
          Ptr3D => EMIS(:,:,:,N)

          ! Update the emissions diagnostics
          CALL Diagn_Update( am_I_Root, HcoState,                 &
                             cName   = TRIM( DiagnName ),         &
                             Array3D = Ptr3D,                     &
                             Total   = TOTFLUX(N),                &
                             COL     = Input_Opt%DIAG_COLLECTION, &
                             RC      = HCRC                        )

          ! Free pointers
          Ptr3D   => NULL()
          SpcInfo => NULL()

          ! Error check
          IF ( HCRC /= HCO_SUCCESS ) THEN
             CALL ERROR_STOP ('Error updating diagnostics: '// DiagnName, &
                              'DO_TEND (mixing_mod.F90)' )
          ENDIF
       ENDIF
    ENDDO

#endif

    ! Convert State_Chm%SPECIES back: kg/m2 --> v/v (ewl, 9/30/15)
    CALL ConvertSpc_Kgm2_to_KgKgDry ( am_I_Root, State_Met, State_Chm, RC )
    CALL ConvertSpc_KgKgDry_to_VVDry( am_I_Root,            State_Chm, RC )

  END SUBROUTINE DO_TEND 
!EOC
END MODULE MIXING_MOD 
