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
! !IROUTINE: INIT_MIXING 
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
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE PBL_MIX_MOD,        ONLY : DO_PBL_MIX
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
!  04 Mar 2015 - C. Keller    - Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! DO_MIXING begins here!
    !=================================================================

    ! Assume success
    RC = GIGC_SUCCESS
   
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

  END SUBROUTINE INIT_MIXING 
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DO_MIXING 
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
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE ERROR_MOD,          ONLY : GIGC_ERROR
    USE PBL_MIX_MOD,        ONLY : DO_PBL_MIX
    USE VDIFF_MOD,          ONLY : DO_PBL_MIX_2
    USE UNITCONV_MOD
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
!  04 Mar 2015 - C. Keller    - Initial version 
!  12 Aug 2015 - E. Lundgren  - Input tracer units are now [kg/kg] and 
!                               are converted to [v/v] for mixing
!  30 Sep 2014 - E. Lundgren  - Move unit conversion for DO_TEND to DO_TEND
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    LOGICAL   :: OnlyAbovePBL

    !=================================================================
    ! DO_MIXING begins here!
    !=================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Convert [kg/kg dry air] to [v/v dry air] for mixing (ewl, 8/12/15)
    CALL Convert_KgKgDry_to_VVDry( am_I_Root, Input_Opt, &
                                   State_Chm, RC )  
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL GIGC_Error('Unit conversion error', RC, &
                       'DO_MIXING in mixing_mod.F')
       RETURN
    ENDIF  

    ! ------------------------------------------------------------------
    ! Do non-local PBL mixing. This will apply the tracer tendencies
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
    ! the PBL to the tracer array. Emissions of some species may be 
    ! capped at the tropopause to avoid build-up in stratosphere.
    ! ------------------------------------------------------------------

    ! Apply tendencies
    CALL DO_TEND ( am_I_Root, Input_Opt,    State_Met, &
                   State_Chm, OnlyAbovePBL, RC          )

    ! ------------------------------------------------------------------
    ! Do full pbl mixing. This fully mixes the updated tracer 
    ! concentrations within the PBL. 
    ! 
    ! Now also archive concentrations and calculate turbulence 
    ! tendencies (ckeller, 7/15/2015)
    ! ------------------------------------------------------------------
    IF ( Input_Opt%LTURB .AND. .NOT. Input_Opt%LNLPBL ) THEN
       CALL DO_PBL_MIX( am_I_Root, Input_Opt%LTURB, Input_Opt, & 
                        State_Met, State_Chm,       RC )
    ENDIF

    ! Convert tracer conc back to [kg/kg dry air] after mixing (ewl, 8/12/15)
    CALL Convert_VVDry_to_KgKgDry( am_I_Root, Input_Opt, &
                                   State_Chm, RC )
    IF ( RC /= GIGC_SUCCESS ) THEN
       CALL GIGC_Error('Unit conversion error', RC, &
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
! !IROUTINE: DO_TEND 
!
! !DESCRIPTION: Subroutine DO\_TEND adds the tracer tendencies (dry deposition
!  and emissions) to the tracer array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DO_TEND ( am_I_Root, Input_Opt,    State_Met, &
                       State_Chm, OnlyAbovePBL, RC,    DT   ) 
!
! !USES:
!
    USE GIGC_ErrCode_Mod
#if defined( DEVEL )
    USE TENDENCIES_MOD
#endif
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE ERROR_MOD,          ONLY : ERROR_STOP, SAFE_DIV
    USE CMN_SIZE_MOD,       ONLY : IIPAR,   JJPAR,   LLPAR
    USE TRACERID_MOD,       ONLY : IDTMACR, IDTRCHO, IDTACET, IDTALD2
    USE TRACERID_MOD,       ONLY : IDTALK4, IDTC2H6, IDTC3H8, IDTCH2O
    USE TRACERID_MOD,       ONLY : IDTPRPE, IDTO3,   IDTHNO3
    USE TRACERID_MOD,       ONLY : IDTBrO,  IDTBr2,  IDTBr,   IDTHOBr
    USE TRACERID_MOD,       ONLY : IDTHBr,  IDTBrNO3 
    USE UNITCONV_MOD
    USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
    USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoVal, GetHcoDiagn
    USE TIME_MOD,           ONLY : GET_TS_DYN
    USE CHEMGRID_MOD,       ONLY : GET_CHEMGRID_LEVEL
    USE DRYDEP_MOD,         ONLY : DEPSAV
    USE GET_NDEP_MOD,       ONLY : SOIL_DRYDEP
    USE PHYSCONSTANTS,      ONLY : AVO
    USE CMN_DIAG_MOD,       ONLY : ND44
#if defined( BPCH )
    USE DIAG_MOD,           ONLY : AD44
#endif
#if defined( NETCDF )
    USE HCO_ERROR_MOD
    USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoID
    USE HCO_DIAGN_MOD,      ONLY : Diagn_Update
    USE TRACERID_MOD,       ONLY : IDTISOPN, IDTMMN
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, L1, L2, N, D, NN
    INTEGER            :: DRYDEPID
    INTEGER            :: PBL_TOP, DRYD_TOP, EMIS_TOP
    REAL(fp)           :: TS, TMP, FRQ, RKT, FRAC, FLUX, AREA_M2
    REAL(fp)           :: MWkg, DENOM
    LOGICAL            :: FND
    LOGICAL            :: PBL_DRYDEP, LSCHEM, ChemGridOnly
    LOGICAL            :: LEMIS,      LDRYD
    LOGICAL            :: DryDepSpec, EmisSpec

    ! For diagnostics
#if defined( NETCDF )
    INTEGER            :: cID, trc_id, HCRC
    CHARACTER(LEN=30)  :: DiagnName
    REAL(fp), TARGET   :: DryDepFlux( IIPAR, JJPAR, Input_Opt%NUMDEP ) 
    REAL(fp), POINTER  :: Ptr2D(:,:)   => NULL()

    ! tracer emissions diagnostic (does not correspond to any bpch diags)
    REAL(fp), TARGET   :: EMIS( IIPAR, JJPAR, LLPAR, Input_Opt%N_TRACERS ) 
    REAL(fp), POINTER  :: Ptr3D(:,:,:) => NULL()
    REAL(fp)           :: TOTFLUX( Input_Opt%N_TRACERS )
#endif

    ! PARANOX loss fluxes (kg/m2/s). These are obtained from the 
    ! HEMCO PARANOX extension via the diagnostics module.
    REAL(fp)                :: PNOXLOSS
    REAL(f4), POINTER, SAVE :: PNOXLOSS_O3  (:,:) => NULL()
    REAL(f4), POINTER, SAVE :: PNOXLOSS_HNO3(:,:) => NULL()

    ! First call?
    LOGICAL,           SAVE :: FIRST = .TRUE.

    !=================================================================
    ! DO_TEND begins here!
    !=================================================================

    ! Assume success
    RC = GIGC_SUCCESS

    ! Special case that there is no dry deposition and emissions
    IF ( .NOT. Input_Opt%LDRYD .AND. .NOT. Input_Opt%LEMIS ) RETURN

    ! Shadow variables
    LSCHEM     = Input_Opt%LSCHEM
    LEMIS      = Input_Opt%LEMIS 
    LDRYD      = Input_Opt%LDRYD 
    PBL_DRYDEP = Input_Opt%PBL_DRYDEP

    ! DO_TEND previously operated in units of kg. The tracer arrays are in
    ! v/v for mixing, hence needed to convert before and after.
    ! CALL CONVERT_UNITS( 2, Input_Opt%N_TRACERS, Input_Opt%TCVV, &
    !                     State_Met%AD, State_Chm%Tracers ) 
    ! Now use units kg/m2 as State_Chm%TRACERS units in DO_TEND to 
    ! remove area-dependency (ewl, 9/30/15)
    ! v/v --> kg/m2
    CALL Convert_VVDry_to_KgKgDry( am_I_Root, Input_Opt,  &
                                   State_Chm, RC )

    CALL Convert_KgKgDry_to_Kgm2( am_I_Root, Input_Opt,   &
                                  State_Met, State_Chm, RC )

    ! Get time step [s]
    IF ( PRESENT(DT) ) THEN
       TS = DT
    ELSE
       TS = GET_TS_DYN() * 60.0_fp
    ENDIF

    ! On first call, get pointers to the PARANOX loss fluxes. These are
    ! stored in diagnostics 'PARANOX_O3_DEPOSITION_FLUX' and 
    ! 'PARANOX_HNO3_DEPOSITION_FLUX'. The call below links pointers 
    ! PNOXLOSS_O3 and PNOXLOSS_HNO3 to the data values stored in the
    ! respective diagnostics. The pointers will remain unassociated if
    ! the diagnostics do not exist.
    ! This is only needed if non-local PBL scheme is not being used. 
    ! Otherwise, PARANOX fluxes are applied in vdiff_mod.F.
    !  (ckeller, 4/10/2015) 
    IF ( FIRST ) THEN
       IF ( .NOT. Input_Opt%LNLPBL ) THEN
          CALL GetHcoDiagn( am_I_Root, 'PARANOX_O3_DEPOSITION_FLUX'  , &
                            .FALSE.,   RC, Ptr2D = PNOXLOSS_O3          ) 
          CALL GetHcoDiagn( am_I_Root, 'PARANOX_HNO3_DEPOSITION_FLUX', &
                            .FALSE.,   RC, Ptr2D = PNOXLOSS_HNO3        ) 
       ENDIF
       FIRST = .FALSE.
    ENDIF

#if defined( NETCDF )
    ! Initialize local diagnostic variables
    DryDepFlux  = 0.0_fp
    EMIS        = 0.0_fp
    TOTFLUX     = 0.0_fp
#endif

    ! Archive concentrations for tendencies (ckeller, 7/15/2015) 
#if defined( DEVEL )
      CALL TEND_STAGE1( am_I_Root, Input_Opt, State_Met, &
                        State_Chm, 'FLUX', .FALSE., RC )
#endif

    ! Do for every tracer and grid box
!$OMP PARALLEL DO                                                    &
!$OMP DEFAULT ( SHARED )                                             &
!$OMP PRIVATE( I, J, L, L1, L2, N, D, PBL_TOP, FND, TMP, DRYDEPID ) &
!$OMP PRIVATE( FRQ, RKT, FRAC, FLUX, AREA_M2,   MWkg, ChemGridOnly ) & 
!$OMP PRIVATE( DryDepSpec, EmisSpec, DRYD_TOP,  EMIS_TOP, PNOXLOSS ) &
!$OMP PRIVATE( DENOM                                               )
    DO N = 1, Input_Opt%N_TRACERS

       !----------------------------------------------------------------
       ! Check if we need to do dry deposition for this species 
       !----------------------------------------------------------------

       ! Initialize
       DryDepSpec = .FALSE.
       DryDepID   = -1

       ! Only if dry deposition is turned on and we do want to consider
       ! processes below the PBL...
       IF ( LDRYD .AND. .NOT. OnlyAbovePBL ) THEN

          ! Get dry deposition index DRYDEPID. This is the ID used by
          ! drydep_mod.F90 for this species. 
          DO D = 1, Input_Opt%NUMDEP
             IF ( Input_Opt%NTRAIND(D) == N ) THEN
                DRYDEPID = D
                EXIT
             ENDIF
          ENDDO ! D

          ! Check if this is a HEMCO drydep species 
          DryDepSpec = ( DRYDEPID > 0 )
          IF ( .NOT. DryDepSpec ) THEN
             CALL GetHcoVal ( N, 1, 1, 1, DryDepSpec, dep = TMP )
          ENDIF

          ! Special case for O3 or HNO3: include PARANOX loss
          IF ( N == IDTO3   .AND. ASSOCIATED(PNOXLOSS_O3  ) ) DryDepSpec = .TRUE. 
          IF ( N == IDTHNO3 .AND. ASSOCIATED(PNOXLOSS_HNO3) ) DryDepSpec = .TRUE. 
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
          MWkg = Input_Opt%TRACER_MW_KG(N)

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
          IF ( N == IDTMACR .OR. N == IDTRCHO .OR. &
               N == IDTACET .OR. N == IDTALD2 .OR. & 
               N == IDTALK4 .OR. N == IDTC2H6 .OR. & 
               N == IDTC3H8 .OR. N == IDTCH2O .OR. & 
               N == IDTPRPE                         ) THEN 
             ChemGridOnly = .TRUE. 
          ENDIF

          ! Bry concentrations become prescribed in lin. strat. chemistry.
          ! Therefore avoid any emissions of these compounds above the 
          ! chemistry grid (lin. strat. chem. applies above chemistry grid
          ! only).
          IF ( LSCHEM ) THEN
             IF ( N == IDTBrO  .OR. N == IDTBr2   .OR. &
                  N == IDTBr   .OR. N == IDTHOBr  .OR. & 
                  N == IDTHBr  .OR. N == IDTBrNO3       ) THEN
                ChemGridOnly = .TRUE.
             ENDIF
          ENDIF

          ! For non-UCX runs, never emit above the chemistry grid.
          ! (ckeller, 6/18/15)
          IF ( .NOT. Input_Opt%LUCX ) ChemGridOnly = .TRUE.

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
                   IF ( N == IDTO3 .AND. ASSOCIATED(PNOXLOSS_O3) ) THEN
                      PNOXLOSS = PNOXLOSS_O3(I,J)
                   ENDIF
                   IF ( N == IDTHNO3 .AND. ASSOCIATED(PNOXLOSS_HNO3) ) THEN
                      PNOXLOSS = PNOXLOSS_HNO3(I,J)
                   ENDIF
                ENDIF

                ! Apply dry deposition
                IF ( FRQ > 0.0_fp .OR. PNOXLOSS > 0.0_fp ) THEN

                   ! Compute exponential loss term
                   RKT  = FRQ * TS
                   FRAC = EXP(-RKT)

                   ! Apply dry deposition
                   State_Chm%Tracers(I,J,L,N) = FRAC *    &
                                                State_Chm%Tracers(I,J,L,N)

                   ! Loss in kg/m2
                   FLUX = ( 1.0_fp - FRAC ) * State_Chm%Tracers(I,J,L,N) 

                   ! Eventually add PARANOX loss. PNOXLOSS is in kg/m2/s. 
                   IF ( PNOXLOSS > 0 ) THEN
                      FLUX = FLUX + ( PNOXLOSS * TS )
                   ENDIF 

#if defined( NETCDF )
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
#if defined( BPCH )
                   IF ( ND44 > 0 ) THEN
                      ! For bpch diagnostic, store data in global AD44 array
                      AD44(I,J,DryDepID,1) = AD44(I,J,DryDepID,1) + FLUX
                   ENDIF
#endif
#if defined( NETCDF )
                   IF ( ND44 > 0 ) THEN
                      ! For netcdf diagnostic, store data in local array
                      ! Now use same units as bpch for comparison (ewl, 1/22/16)
                      DryDepFlux(I,J,DryDepID) = DryDepFlux(I,J,DryDepID) &
                                                + FLUX
                   ENDIF
#endif

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

                   ! Add to tracer array
                   State_Chm%Tracers(I,J,L,N) = State_Chm%Tracers(I,J,L,N) & 
                                              + FLUX 

#if defined( NETCDF )
                   ! Update new tracer emissions diagnostics
                   EMIS(I,J,L,N) = TMP                 ! kg/m2/s
                   TOTFLUX(N)    = TOTFLUX(N) + FLUX   ! kg/m2/s
#endif
                ENDIF
             ENDIF

          ENDDO !L
       ENDDO !J
       ENDDO !I
    ENDDO !N
!$OMP END PARALLEL DO

      ! Calculate tendencies and write to diagnostics (ckeller, 7/15/2015)
#if defined( DEVEL )
      CALL TEND_STAGE2( am_I_Root, Input_Opt, State_Met, &
                        State_Chm, 'FLUX', .FALSE., TS, RC )
#endif


#if defined( NETCDF )
    !========================================================      
    ! ND44: Dry deposition diagnostic [molec/cm2/s] (netcdf) 
    !========================================================
    IF ( Input_Opt%ND44 > 0 ) THEN

       ! Assume success
       HCRC = HCO_SUCCESS

       ! Loop over all depositing species
       DO D = 1, Input_Opt%NUMDEP
       
          ! Get the corresponding GEOS-Chem tracer number
          trc_id = Input_Opt%NTRAIND( D )
       
          ! If this tracer number is scheduled for output in input.geos, 
          ! then archive the latest depvel data into the diagnostic structure
          IF ( ANY( Input_Opt%TINDEX(44,:) == trc_id ) ) THEN
       
             ! Define diagnostic name
             DiagnName = 'DRYDEP_FLX_MIX_' // TRIM( Input_OPt%DEPNAME(D) )

             ! Point to data
             Ptr2D => DryDepFlux(:,:,D)

             ! Update diagnostic container
             CALL Diagn_Update( am_I_Root,                           &
                                cName     = TRIM( DiagnName),        &
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

    ! For now, always output tracer emissions diagnostic since there
    ! is no logical switch for it in Input_Mod, nor an entry in input.geos
    DO N = 1, Input_Opt%N_TRACERS
 
       ! Skip if there are no emissions
       IF ( TOTFLUX(N) == 0.0_fp ) CYCLE

       ! Only if HEMCO tracer is defined
       ! NOTE: For other netcdf diagnostics we look at Input_Opt%TINDEX
       ! to see which tracers are turned on for the diagnostic
       cID = GetHcoID( TrcID=N )
       IF ( cID > 0 ) THEN 

          ! Define diagnostics name
          DiagnName = 'TRACER_EMIS_' // TRIM( Input_Opt%TRACER_NAME(N) )

          ! Point to species slice 
          Ptr3D => EMIS(:,:,:,N)

          ! Update the emissions diagnostics
          CALL Diagn_Update( am_I_Root,                           &
                             cName   = TRIM( DiagnName ),         &
                             Array3D = Ptr3D,                     &
                             Total   = TOTFLUX(N),                &
                             COL     = Input_Opt%DIAG_COLLECTION, &
                             RC      = HCRC                        )

          ! Free pointer
          Ptr3D => NULL()

          ! Error check
          IF ( HCRC /= HCO_SUCCESS ) THEN
             CALL ERROR_STOP ('Error updating diagnostics: '// DiagnName, &
                              'DO_TEND (mixing_mod.F90)' )
          ENDIF
       ENDIF
    ENDDO

#endif

    ! Convert State_Chm%TRACERS back: kg/m2 --> v/v (ewl, 9/30/15)
    CALL Convert_Kgm2_to_KgKgDry( am_I_Root, Input_Opt,   &
                                  State_Met, State_Chm, RC )
    CALL Convert_KgKgDry_to_VVDry( am_I_Root, Input_Opt,          &
                                   State_Chm, RC )

  END SUBROUTINE DO_TEND 
!EOC
END MODULE MIXING_MOD 
