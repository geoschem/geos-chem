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
    USE DAO_MOD,            ONLY : CONVERT_UNITS
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
       CALL DO_PBL_MIX( .FALSE., Input_Opt, State_Met, State_Chm )
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
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE PBL_MIX_MOD,        ONLY : DO_PBL_MIX
    USE VDIFF_MOD,          ONLY : DO_PBL_MIX_2
    USE DAO_MOD,            ONLY : CONVERT_UNITS
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
    LOGICAL       :: OnlyAbovePBL

    !=================================================================
    ! DO_MIXING begins here!
    !=================================================================

    ! Assume success
    RC = GIGC_SUCCESS

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

    ! DO_TEND operates in units of kg. The tracer arrays enters as
    ! v/v, hence need to convert before and after.
    ! v/v --> kg
    CALL CONVERT_UNITS( 2, Input_Opt%N_TRACERS, Input_Opt%TCVV, &
                        State_Met%AD, State_Chm%Tracers ) 

    ! Apply tendencies
    CALL DO_TEND ( am_I_Root, Input_Opt,    State_Met, &
                   State_Chm, OnlyAbovePBL, RC          )

    ! Convert back: kg --> v/v
    CALL CONVERT_UNITS( 1, Input_Opt%N_TRACERS, Input_Opt%TCVV, &
                        State_Met%AD, State_Chm%Tracers ) 

    ! ------------------------------------------------------------------
    ! Do full pbl mixing. This fully mixes the updated tracer 
    ! concentrations within the PBL. 
    ! ------------------------------------------------------------------
    IF ( Input_Opt%LTURB .AND. .NOT. Input_Opt%LNLPBL ) THEN
       CALL DO_PBL_MIX( Input_Opt%LTURB, Input_Opt, State_Met, State_Chm )
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
    USE GIGC_Input_Opt_Mod, ONLY : OptInput
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE GIGC_State_Chm_Mod, ONLY : ChmState
    USE ERROR_MOD,          ONLY : ERROR_STOP
    USE CMN_SIZE_MOD,       ONLY : IIPAR,   JJPAR,   LLPAR
    USE TRACERID_MOD,       ONLY : IDTMACR, IDTRCHO, IDTACET, IDTALD2
    USE TRACERID_MOD,       ONLY : IDTALK4, IDTC2H6, IDTC3H8, IDTCH2O
    USE TRACERID_MOD,       ONLY : IDTPRPE
    USE TRACERID_MOD,       ONLY : IDTBrO,  IDTBr2,  IDTBr,   IDTHOBr
    USE TRACERID_MOD,       ONLY : IDTHBr,  IDTBrNO3 
    USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
    USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoVal
    USE TIME_MOD,           ONLY : GET_TS_DYN
    USE CHEMGRID_MOD,       ONLY : GET_CHEMGRID_LEVEL
    USE DRYDEP_MOD,         ONLY : NTRAIND, DEPSAV
    USE GET_NDEP_MOD,       ONLY : SOIL_DRYDEP
    USE CMN_GCTM_MOD,       ONLY : AVO
    USE CMN_DIAG_MOD,       ONLY : ND44
#if !defined( NO_BPCH )
    USE DIAG_MOD,           ONLY : AD44
#endif
#if defined( DEVEL )
    USE HCO_ERROR_MOD
    USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoID
    USE HCO_DIAGN_MOD,      ONLY : Diagn_Update
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER            :: I, J, L, L1, L2, N, NN
    INTEGER            :: DRYDEPID
    INTEGER            :: PBL_TOP, DRYD_TOP, EMIS_TOP
    REAL(fp)           :: TS, TMP, FRQ, RKT, FRAC, FLUX, AREA_M2
    REAL(fp)           :: MWkg 
    LOGICAL            :: FND
    LOGICAL            :: PBL_DRYDEP, LSCHEM, ChemGridOnly
    LOGICAL            :: LEMIS,      LDRYD
    LOGICAL            :: DryDepSpec, EmisSpec

    ! For diagnostics
#if defined( DEVEL )
    INTEGER            :: cID, HCRC
    REAL(fp), POINTER  :: Ptr3D(:,:,:) => NULL()
    REAL(fp), TARGET   :: EMIS(IIPAR,JJPAR,LLPAR,Input_Opt%N_TRACERS) 
    REAL(fp)           :: TOTFLUX(Input_Opt%N_TRACERS)
#endif

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

    ! Get time step [s]
    IF ( PRESENT(DT) ) THEN
       TS = DT
    ELSE
       TS = GET_TS_DYN() * 60.0_fp
    ENDIF

    ! Init diagnostics
#if defined( DEVEL )
    EMIS    = 0.0_fp
    TOTFLUX = 0.0_fp
#endif

    ! Do for every tracer and grid box
!$OMP PARALLEL DO                                                      &
!$OMP DEFAULT ( SHARED )                                               &
!$OMP PRIVATE( I, J, L, L1, L2, N, NN, PBL_TOP, FND, TMP, DRYDEPID   ) &
!$OMP PRIVATE( FRQ, RKT, FRAC, FLUX, AREA_M2,   MWkg, ChemGridOnly   ) & 
!$OMP PRIVATE( DryDepSpec, EmisSpec, DRYD_TOP,  EMIS_TOP             )
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
          DO NN = 1, Input_Opt%NUMDEP
             IF ( NTRAIND(NN) == N ) THEN
                DRYDEPID = NN
                EXIT
             ENDIF
          ENDDO !NN

          ! Check if this is a HEMCO drydep species 
          DryDepSpec = ( DRYDEPID > 0 )
          IF ( .NOT. DryDepSpec ) THEN
             CALL GetHcoVal ( N, 1, 1, 1, DryDepSpec, dep = TMP )
          ENDIF
       ENDIF

       !----------------------------------------------------------------
       ! Check if we need to do emisisons for this species 
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

          ! Grid box area
          AREA_M2 = State_Met%AREA_M2(I,J,1)

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

                ! Apply dry deposition
                IF ( FRQ > 0.0_fp ) THEN

                   ! Compute exponential loss term
                   RKT  = FRQ * TS
                   FRAC = EXP(-RKT)

                   ! Apply dry deposition
                   State_Chm%Tracers(I,J,L,N) = FRAC * State_Chm%Tracers(I,J,L,N)

                   ! Loss in [molec/cm2/s]
                   FLUX = ( 1.0_fp - FRAC ) * State_Chm%Tracers(I,J,L,N)     ! kg
                   FLUX = FLUX / MWkg * AVO / TS / ( AREA_M2 * 1.0e4_fp ) ! molec/cm2/s

                   ! Diagnostics. These are the same as DRYFLX.
                   ! Diagnostics are in molec/cm2/s.
#if !defined( NO_BPCH )
                   ! ND44 diagnostics
                   IF ( ND44 > 0 .and. DryDepID > 0 ) THEN
                      AD44(I,J,DryDepID,1) = AD44(I,J,DryDepID,1) + FLUX
                   ENDIF
#endif
                   ! Eventually add to SOIL_DRYDEP
                   IF ( Input_Opt%LSOILNOX ) THEN
                      CALL SOIL_DRYDEP( I, J, L, N, FLUX )
                   ENDIF
                ENDIF

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

                   ! Flux: [kg] = [kg m-2 s-1 ] x [s] x [m2]
                   FLUX = TMP * TS * AREA_M2

                   ! Add to tracer array
                   State_Chm%Tracers(I,J,L,N) = State_Chm%Tracers(I,J,L,N) & 
                                              + FLUX 

                   ! Update diagnostics
#if defined( DEVEL )
                   EMIS(I,J,L,N) = TMP
                   TOTFLUX(N)    = TOTFLUX(N) + FLUX 
#endif
                ENDIF
             ENDIF

          ENDDO !L
       ENDDO !J
       ENDDO !I
    ENDDO !N
!$OMP END PARALLEL DO

#if defined( DEVEL )
    !-------------------------------------------------------------------
    ! Update diagnostics that will get saved to netCDF files.
    ! These are defined in diagnostics_mod.F90
    ! 
    ! NOTE: For now, this is only activated by compiling with DEVEL=y,
    ! but in the future this will replace the bpch diagnostics!
    !-------------------------------------------------------------------
    DO N = 1, Input_Opt%N_TRACERS

       ! Only if HEMCO tracer is defined
       cID = GetHcoID( TrcID=N )
       IF ( cID <= 0 ) CYCLE

       ! Define diagnostics name
       cID = 10000 + cID

       ! Point to species slice 
       Ptr3D => EMIS(:,:,:,N)

       ! Update the emissions diagnostics
       CALL Diagn_Update( am_I_Root,                           &
                          cID     = cID,                       &
                          Array3D = Ptr3D,                     &
                          Total   = TOTFLUX(N),                &
                          COL     = Input_Opt%DIAG_COLLECTION, &
                          RC      = HCRC                        )

       ! Free pointer
       Ptr3D => NULL()

       ! Error check
       IF ( HCRC /= HCO_SUCCESS ) THEN
          CALL ERROR_STOP ('Error updating diagnostics: '// &
                           'TRACER_EMIS_'//TRIM(Input_Opt%TRACER_NAME(N)), &
                           'DO_TEND (mixing_mod.F90)' )
       ENDIF
    ENDDO
#endif

  END SUBROUTINE DO_TEND 
!EOC
END MODULE MIXING_MOD 
