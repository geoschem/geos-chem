!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fullchem_mod.F90
!
! !DESCRIPTION: Contines arrays and routines for the GEOS_Chem "fullchem"
!  mechanism, which is implemented in KPP-generated Fortran code.
!\\
!\\
! !INTERFACE:
!
MODULE FullChem_Mod
!
! !USES:
!
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Do_FullChem
  PUBLIC  :: Init_FullChem
  PUBLIC  :: Cleanup_FullChem
!
! !REVISION HISTORY:
!  14 Dec 2015 - M.S. Long   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Species ID flags (and logicals to denote if species are present)
  INTEGER               :: id_OH,  id_HO2,  id_O3P,  id_O1D, id_CH4
  INTEGER               :: id_PCO, id_LCH4, id_NH3
#ifdef MODEL_GEOS
  INTEGER               :: id_O3
  INTEGER               :: id_A3O2, id_ATO2, id_B3O2, id_BRO2
  INTEGER               :: id_ETO2, id_LIMO2, id_MO2, id_PIO2, id_PO2
  INTEGER               :: id_PRN1, id_R4N1, id_R4O2, id_TRO2, id_XRO2
  INTEGER               :: id_IHOO1, id_IHOO4, id_IHCOO, id_ICHOO, id_IHPOO1
  INTEGER               :: id_IHPOO2, id_IHPOO3, id_IEPOXAOO, id_IEPOXBOO
  INTEGER               :: id_C4HVP1, id_C4HVP2, id_HPALD1OO, id_HPALD2OO
  INTEGER               :: id_ISOPNOO1, id_ISOPNOO2, id_INO2B, id_INO2D
  INTEGER               :: id_IDHNBOO, id_IDHNDOO1, id_IDHNDOO2
  INTEGER               :: id_IHPNBOO, id_IHPNDOO, id_ICNOO, id_IDNOO
#endif
  INTEGER               :: id_SALAAL, id_SALCAL, id_SO4, id_SALC ! MSL
  LOGICAL               :: ok_OH, ok_HO2, ok_O1D, ok_O3P
  LOGICAL               :: Failed2x

  ! Diagnostic flags
  LOGICAL               :: Do_Diag_OH_HO2_O1D_O3P
#ifdef MODEL_GEOS
  LOGICAL               :: Archive_O3concAfterchem
  LOGICAL               :: Archive_RO2concAfterchem
#endif

  ! SAVEd scalars
  INTEGER,  SAVE        :: PrevDay   = -1
  INTEGER,  SAVE        :: PrevMonth = -1

  ! Arrays
  INTEGER,  ALLOCATABLE :: PL_Kpp_ID (:      )
  REAL(f4), ALLOCATABLE :: JvCountDay(:,:,:  )
  REAL(f4), ALLOCATABLE :: JvCountMon(:,:,:  )
  REAL(f4), ALLOCATABLE :: JvSumDay  (:,:,:,:)
  REAL(f4), ALLOCATABLE :: JvSumMon  (:,:,:,:)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: do_fullchem
!
! !DESCRIPTION: Driver subroutine for the KPP fullchem mechanism.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_FullChem( Input_Opt,  State_Chm, State_Diag,                 &
                          State_Grid, State_Met, RC                         )
!
! !USES:
!
    USE CMN_FJX_MOD
    USE ErrCode_Mod
    USE ERROR_MOD
    USE FAST_JX_MOD,              ONLY : PHOTRATE_ADJ, FAST_JX
    USE FAST_JX_MOD,              ONLY : RXN_O3_1, RXN_NO2
    USE fullchem_AutoReduceFuncs, ONLY : fullchem_AR_KeepHalogensActive
    USE fullchem_AutoReduceFuncs, ONLY : fullchem_AR_SetKeepActive
    USE fullchem_AutoReduceFuncs, ONLY : fullchem_AR_UpdateKppDiags
    USE fullchem_AutoReduceFuncs, ONLY : fullchem_AR_SetIntegratorOptions
    USE fullchem_HetStateFuncs,   ONLY : fullchem_SetStateHet
    USE fullchem_SulfurChemFuncs, ONLY : fullchem_ConvertAlkToEquiv
    USE fullchem_SulfurChemFuncs, ONLY : fullchem_ConvertEquivToAlk
    USE fullchem_SulfurChemFuncs, ONLY : fullchem_HetDropChem
    USE GcKpp_Monitor,            ONLY : SPC_NAMES, FAM_NAMES, EQN_NAMES
    USE GcKpp_Parameters
    USE GcKpp_Integrator,         ONLY : Integrate
    USE GcKpp_Function
    USE GcKpp_Global
    USE GcKpp_Rates,              ONLY : UPDATE_RCONST, RCONST
    USE GcKpp_Util,               ONLY : Get_OHreactivity
    USE Input_Opt_Mod,            ONLY : OptInput
    USE PhysConstants,            ONLY : AVO, AIRMW
    USE PRESSURE_MOD
    USE Species_Mod,              ONLY : Species
    USE State_Chm_Mod,            ONLY : ChmState
    USE State_Chm_Mod,            ONLY : Ind_
    USE State_Diag_Mod,           ONLY : DgnState
    USE State_Diag_Mod,           ONLY : DgnMap
    USE State_Grid_Mod,           ONLY : GrdState
    USE State_Met_Mod,            ONLY : MetState
    USE TIME_MOD,                 ONLY : GET_TS_CHEM
    USE TIME_MOD,                 ONLY : Get_Day
    USE TIME_MOD,                 ONLY : Get_Month
    USE TIME_MOD,                 ONLY : Get_Year
    USE Timers_Mod
    USE UnitConv_Mod,             ONLY : Convert_Spc_Units
    USE UCX_MOD,                  ONLY : CALC_STRAT_AER
    USE UCX_MOD,                  ONLY : SO4_PHOTFRAC
    USE UCX_MOD,                  ONLY : UCX_NOX
    USE UCX_MOD,                  ONLY : UCX_H2SO4PHOT
#ifdef TOMAS
#ifdef BPCH_DIAG
    USE TOMAS_MOD,                ONLY : H2SO4_RATE
#endif
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met  ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure
!
! !REVISION HISTORY:
!  14 Dec 2015 - M.S. Long   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                :: prtDebug,   IsLocNoon, Size_Res, Failed2x
    INTEGER                :: I,          J,         L,        N
    INTEGER                :: NA,         F,         SpcID,    KppID
    INTEGER                :: P,          MONTH,     YEAR,     Day
    INTEGER                :: WAVELENGTH, IERR,      S,        Thread
    REAL(fp)               :: SO4_FRAC,   T,         TIN
    REAL(fp)               :: TOUT,       SR,        LWC

    ! Strings
    CHARACTER(LEN=63)      :: OrigUnit
    CHARACTER(LEN=255)     :: ErrMsg,   ThisLoc

    ! SAVEd scalars
    LOGICAL,  SAVE         :: FIRSTCHEM = .TRUE.
    INTEGER,  SAVE         :: CH4_YEAR  = -1

    ! For

#ifdef MODEL_CLASSIC
#ifndef NO_OMP
    INTEGER, EXTERNAL      :: OMP_GET_THREAD_NUM
#endif
#endif

    ! Arrays
    INTEGER                :: ICNTRL (20)
    INTEGER                :: ISTATUS(20)
    REAL(dp)               :: RCNTRL (20)
    REAL(dp)               :: RSTATE (20)
    REAL(fp)               :: Before(State_Grid%NX, State_Grid%NY,           &
                                     State_Grid%NZ, State_Chm%nAdvect       )

    ! For tagged CO saving
    REAL(fp)               :: LCH4, PCO_TOT, PCO_CH4, PCO_NMVOC

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    ! For testing purposes
    LOGICAL                :: DO_PHOTCHEM

    ! OH reactivity and KPP reaction rate diagnostics
    REAL(fp)               :: OHreact
    REAL(dp)               :: Vloc(NVAR),     Aout(NREACT)
#ifdef MODEL_GEOS
    REAL(f4)               :: NOxTau, NOxConc, NOx_weight, NOx_tau_weighted
    REAL(f4)               :: TROP_NOx_Tau
    REAL(f4)               :: TROPv_NOx_tau(State_Grid%NX,State_Grid%NY)
    REAL(f4)               :: TROPv_NOx_mass(State_Grid%NX,State_Grid%NY)
    REAL(dp)               :: localC(NSPEC)
#endif
#ifdef MODEL_WRF
    REAL(dp)               :: localC(NSPEC)
#endif

    ! Grid box integration time diagnostic
    REAL(fp)               :: TimeStart, TimeEnd

    ! Objects
    TYPE(DgnMap), POINTER :: mapData => NULL()
!
! !DEFINED PARAMETERS
!
    ! Defines the slot in which the H-value from the KPP integrator is stored
    ! This should be the same as the value of Nhnew in gckpp_Integrator.F90
    ! (assuming Rosenbrock solver).  Define this locally in order to break
    ! a compile-time dependency.  -- Bob Yantosca (05 May 2022)
    INTEGER,    PARAMETER :: Nhnew = 3

    !========================================================================
    ! Do_FullChem begins here!
    ! NOTE: FlexChem timer is started in DO_CHEMISTRY (the calling routine)
    !========================================================================

    ! Initialize
    RC        =  GC_SUCCESS
    ErrMsg    =  ''
    ThisLoc   =  ' -> at Do_FullChem (in module GeosCore/FullChem_mod.F90)'
    SpcInfo   => NULL()
    prtDebug  =  ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    Day       =  Get_Day()    ! Current day
    Month     =  Get_Month()  ! Current month
    Year      =  Get_Year()   ! Current year
    Thread    =  1
    Failed2x  = .FALSE.

    ! Set a switch that allows you to toggle off photolysis for testing
    ! (default value : TRUE)
    DO_PHOTCHEM = .TRUE.

    ! Print information the first time that DO_FULLCHEM is called
    CALL PrintFirstTimeInfo( Input_Opt, State_Chm, FirstChem, Do_PhotChem )

    ! Zero diagnostic archival arrays to make sure that we don't have any
    ! leftover values from the last timestep near the top of the chemgrid
    IF (State_Diag%Archive_Loss           ) State_Diag%Loss           = 0.0_f4
    IF (State_Diag%Archive_Prod           ) State_Diag%Prod           = 0.0_f4
    IF (State_Diag%Archive_Jval           ) State_Diag%Jval           = 0.0_f4
    IF (State_Diag%Archive_JvalO3O1D      ) State_Diag%JvalO3O1D      = 0.0_f4
    IF (State_Diag%Archive_JvalO3O3P      ) State_Diag%JvalO3O3P      = 0.0_f4
    IF (State_Diag%Archive_JNoon          ) State_Diag%JNoon          = 0.0_f4
    IF (State_Diag%Archive_ProdCOfromCH4  ) State_Diag%ProdCOfromCH4  = 0.0_f4
    IF (State_Diag%Archive_ProdCOfromNMVOC) State_Diag%ProdCOfromNMVOC= 0.0_f4
    IF (State_Diag%Archive_OHreactivity   ) State_Diag%OHreactivity   = 0.0_f4
    IF (State_Diag%Archive_RxnRate        ) State_Diag%RxnRate        = 0.0_f4
    IF (State_Diag%Archive_SatDiagnRxnRate) State_Diag%SatDiagnRxnRate= 0.0_f4
    IF (State_Diag%Archive_KppDiags) THEN
       IF (State_Diag%Archive_KppIntCounts) State_Diag%KppIntCounts   = 0.0_f4
       IF (State_Diag%Archive_KppJacCounts) State_Diag%KppJacCounts   = 0.0_f4
       IF (State_Diag%Archive_KppTotSteps ) State_Diag%KppTotSteps    = 0.0_f4
       IF (State_Diag%Archive_KppAccSteps ) State_Diag%KppAccSteps    = 0.0_f4
       IF (State_Diag%Archive_KppRejSteps ) State_Diag%KppRejSteps    = 0.0_f4
       IF (State_Diag%Archive_KppLuDecomps) State_Diag%KppLuDecomps   = 0.0_f4
       IF (State_Diag%Archive_KppSubsts   ) State_Diag%KppSubsts      = 0.0_f4
       IF (State_Diag%Archive_KppSmDecomps) State_Diag%KppSmDecomps   = 0.0_f4
       IF (State_Diag%Archive_KppAutoReducerNVAR)                            &
                                      State_Diag%KppAutoReducerNVAR   = 0.0_f4
       IF (State_Diag%Archive_KppcNONZERO)  State_Diag%KppcNONZERO    = 0.0_f4
    ENDIF
    
    ! Also zero satellite diagnostic archival arrays
    IF ( State_Diag%Archive_SatDiagnLoss ) State_Diag%SatDiagnLoss    = 0.0_f4
    IF ( State_Diag%Archive_SatDiagnProd ) State_Diag%SatDiagnProd    = 0.0_f4
    IF ( State_Diag%Archive_SatDiagnJval ) THEN
       State_Diag%SatDiagnJval = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_SatDiagnJvalO3O1D ) THEN
       State_Diag%SatDiagnJvalO3O1D = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_SatDiagnJvalO3O3P ) THEN
       State_Diag%SatDiagnJvalO3O3P = 0.0_f4
    ENDIF

    ! Keep track of the boxes where it is local noon in the JNoonFrac
    ! diagnostic. When time-averaged, this will be the fraction of time
    ! that local noon occurred at a grid box. (bmy, 4/2/19)
    IF ( State_Diag%Archive_JNoonFrac ) THEN
       WHERE( State_Met%IsLocalNoon )
          State_Diag%JNoonFrac = 1.0_f4
       ELSEWHERE
          State_Diag%JNoonFrac = 0.0_f4
       ENDWHERE
    ENDIF

#if defined( MODEL_GEOS )
    IF ( State_Diag%Archive_NoxTau     ) State_Diag%NoxTau(:,:,:) = 0.0_f4
    IF ( State_Diag%Archive_TropNOxTau ) THEN
       State_Diag%TropNOxTau(:,:) = 0.0_f4
       TROPv_NOx_mass(:,:) = 0.0_f4
       TROPv_NOx_tau(:,:)  = 0.0_f4
    ENDIF
#endif

    !========================================================================
    ! Zero out certain species:
    !    - isoprene oxidation counter species (dkh, bmy, 6/1/06)
    !    - isoprene-NO3 oxidation counter species (hotp, 6/1/10)
    !    - if SOA or SOA_SVPOA, aromatic oxidation counter species
    !      (dkh, 10/06/06)
    !    - if SOA_SVPOA, LNRO2H and LNRO2N for NAP (hotp 6/25/09
    !========================================================================
    DO N = 1, State_Chm%nSpecies

       ! Get info about this species from the species database
       SpcInfo => State_Chm%SpcData(N)%Info

       ! isoprene oxidation counter species
       IF ( TRIM( SpcInfo%Name ) == 'LISOPOH' .or. &
            TRIM( SpcInfo%Name ) == 'LISOPNO3' ) THEN
          State_Chm%Species(N)%Conc(:,:,:) = 0.0_fp
       ENDIF

       ! aromatic oxidation counter species
       IF ( Input_Opt%LSOA .or. Input_Opt%LSVPOA ) THEN
          SELECT CASE ( TRIM( SpcInfo%Name ) )
             CASE ( 'LBRO2H', 'LBRO2N', 'LTRO2H', 'LTRO2N', &
                    'LXRO2H', 'LXRO2N', 'LNRO2H', 'LNRO2N' )
                State_Chm%Species(N)%Conc(:,:,:) = 0.0_fp
          END SELECT
       ENDIF

       ! Temporary fix for CO2
       ! CO2 is a dead species and needs to be set to zero to
       ! match the old SMVGEAR code (mps, 6/14/16)
       IF ( TRIM( SpcInfo%Name ) == 'CO2' ) THEN
          State_Chm%Species(N)%Conc(:,:,:) = 0.0_fp
       ENDIF

       ! Free pointer
       SpcInfo => NULL()

    ENDDO

    !========================================================================
    ! Convert species to [molec/cm3] (ewl, 8/16/16)
    !========================================================================
    CALL Convert_Spc_Units( Input_Opt,            State_Chm,   State_Grid,   &
                            State_Met,           'molec/cm3', RC,            &
                            OrigUnit=OrigUnit                               )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'fullchem_mod.F90')
       RETURN
    ENDIF

    !========================================================================
    ! Call photolysis routine to compute J-Values
    !========================================================================
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End  ( "=> FlexChem",           RC )
       CALL Timer_Start( "=> FAST-JX photolysis", RC )
    ENDIF

    ! Do Photolysis
    WAVELENGTH = 0
    CALL FAST_JX( WAVELENGTH, Input_Opt,  State_Chm,                         &
                  State_Diag, State_Grid, State_Met, RC                     )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "FAST_JX"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !### Debug
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### Do_FullChem: after FAST_JX' )
    ENDIF

    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End  ( "=> FAST-JX photolysis", RC )
       CALL Timer_Start( "=> FlexChem",           RC ) ! ended in Do_Chemistry
    ENDIF

#if defined( MODEL_GEOS ) || defined( MODEL_WRF )
    ! Init diagnostics
    IF ( ASSOCIATED(State_Diag%KppError) ) THEN
       State_Diag%KppError(:,:,:) = 0.0
    ENDIF
#endif

    !=======================================================================
    ! Archive concentrations before chemistry
    !=======================================================================
    IF ( State_Diag%Archive_ConcBeforeChem ) THEN
       ! Point to mapping obj specific to ConcBeforeChem diagnostic collection
       mapData => State_Diag%Map_ConcBeforeChem

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( N, S   )
       DO S = 1, mapData%nSlots
          N = mapData%slot2id(S)
          State_Diag%ConcBeforeChem(:,:,:,S) = State_Chm%Species(N)%Conc(:,:,:)
       ENDDO
       !$OMP END PARALLEL DO

       ! Free pointer
       mapData => NULL()
    ENDIF

    !========================================================================
    ! Set up integration convergence conditions and timesteps
    ! (cf. M. J. Evans)
    !
    ! NOTE: ATOL and RTOL are defined in gckpp_Global.F90 so they
    ! are probably only used as INTENT(IN).  Therefore, it is
    ! probably safe to define them here outside the OpenMP loop.
    ! (bmy, 3/28/16)
    !
    ! The ICNTRL vector specifies options for the solver.  We can
    ! define ICNTRL outside of the "SOLVE CHEMISTRY" parallel loop
    ! below because it is passed to KPP as INTENT(IN).
    ! (bmy, 3/28/16)
    !
    ! ICNTRL now needs to be updated within the TimeLoop for the AR solver
    ! (hplin, 4/13/22)
    !========================================================================

    !%%%%% TIMESTEPS %%%%%

    ! mje Set up conditions for the integration
    ! mje chemical timestep and convert it to seconds.
    DT        = GET_TS_CHEM() ! [s]
    T         = 0d0
    TIN       = T
    TOUT      = T + DT

    !%%%%% CONVERGENCE CRITERIA %%%%%

    ! Absolute tolerance
    ATOL      = 1e-2_dp

    ! Relative tolerance
    RTOL      = 0.5e-2_dp

    !=======================================================================
    ! %%%%% SOLVE CHEMISTRY -- This is the main KPP solver loop %%%%%
    !=======================================================================
100 format('No. of function calls:', i6, /,                                 &
           'No. of jacobian calls:', i6, /,                                 &
           'No. of steps:         ', i6, /,                                 &
           'No. of accepted steps:', i6, /,                                 &
           'No. of rejected steps ', i6, /,                                 &
           '       (except at very beginning)',          /,                 &
           'No. of LU decompositions:             ', i6, /,                 &
           'No. of forward/backward substitutions:', i6, /,                 &
           'No. of singular matrix decompositions:', i6, /,                 &
            /,                                                              &
           'Texit, the time corresponding to the      ',        /,          &
           '       computed Y upon return:            ', f11.4, /,          &
           'Hexit, last accepted step before exit:    ', f11.4, /,          &
           'Hnew, last predicted step (not yet taken):', f11.4 )

    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_Start( "  -> FlexChem loop", RC )
    ENDIF

    !------------------------------------------------------------------------
    ! Always consider halogens as "fast" species for auto-reduce
    !------------------------------------------------------------------------
    IF ( FIRSTCHEM .and. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
       IF ( Input_Opt%AutoReduce_Is_KeepActive ) THEN
          CALL fullchem_AR_KeepHalogensActive( Input_Opt%amIRoot )
       ENDIF
    ENDIF

    !========================================================================
    ! MAIN LOOP: Compute reaction rates and call chemical solver
    !
    ! Variables not listed here are held THREADPRIVATE in gckpp_Global.F90
    ! !$OMP COLLAPSE(3) vectorizes the loop and !$OMP DYNAMIC(24) sends
    ! 24 boxes at a time to each core... then when that core is finished,
    ! it gets another chunk of 24 boxes.  This should lead to better
    ! load balancing, and will spread the sunrise/sunset boxes across
    ! more cores.
    !========================================================================
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I,        J,        L,       N                           )&
    !$OMP PRIVATE( ICNTRL                                                   )&
    !$OMP PRIVATE( SO4_FRAC, IERR,     RCNTRL,  ISTATUS,   RSTATE           )&
    !$OMP PRIVATE( SpcID,    KppID,    F,       P,         Vloc             )&
    !$OMP PRIVATE( Aout,     Thread,   RC,      S,         LCH4             )&
    !$OMP PRIVATE( OHreact,  PCO_TOT,  PCO_CH4, PCO_NMVOC, SR               )&
    !$OMP PRIVATE( SIZE_RES, LWC                                            )&
#ifdef MODEL_GEOS
    !$OMP PRIVATE( NOxTau,     NOxConc, localC                              )&
    !$OMP PRIVATE( NOx_weight, NOx_tau_weighted                             )&
#endif
#ifdef MODEL_WRF
    !$OMP PRIVATE( localC                                                   )&
#endif
    !$OMP COLLAPSE( 3                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 24                                             )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !=====================================================================
       ! Initialize private loop variables for each (I,J,L)
       ! Other private variables will be assigned in Set_Kpp_GridBox_Values
       !=====================================================================
       IERR      = 0                        ! KPP success or failure flag
       ISTATUS   = 0.0_dp                   ! Rosenbrock output
       ICNTRL    = 0                        ! Rosenbrock input (integer)
       RCNTRL    = 0.0_fp                   ! Rosenbrock input (real)
       RSTATE    = 0.0_dp                   ! Rosenbrock output
       SO4_FRAC  = 0.0_fp                   ! Frac of SO4 avail for photolysis
       P         = 0                        ! GEOS-Chem photolyis species ID
       LCH4      = 0.0_fp                   ! P/L diag: Methane loss rate
       PCO_TOT   = 0.0_fp                   ! P/L diag: Total P(CO)
       PCO_CH4   = 0.0_fp                   ! P/L diag: P(CO) from CH4
       PCO_NMVOC = 0.0_fp                   ! P/L diag: P(CO) from NMVOC
       SR        = 0.0_fp                   ! Enhancement to O2 catalysis rate
       LWC       = 0.0_fp                   ! Liquid water content
       SIZE_RES  = .FALSE.                  ! Size resolved calculation?
       C         = 0.0_dp                   ! KPP species conc's
       RCONST    = 0.0_dp                   ! KPP rate constants
       PHOTOL    = 0.0_dp                   ! Photolysis array for KPP
       K_CLD     = 0.0_dp                   ! Sulfur in-cloud rxn het rates
       K_MT      = 0.0_dp                   ! Sulfur sea salt rxn het rates
       CFACTOR   = 1.0_dp                   ! KPP conversion factor
#ifdef MODEL_CLASSIC
#ifndef NO_OMP
       Thread    = OMP_GET_THREAD_NUM() + 1 ! OpenMP thread number
#endif
#endif
#if defined( MODEL_GEOS ) || defined( MODEL_WRF )
       localC    = 0.0_dp                   ! Local backup array for C
#endif

       ! Per discussions for Lin et al., force keepActive throughout the atmosphere
       ! if keepActive option is enabled. (hplin, 2/9/22)
       !keepActive = .true.
       CALL fullchem_AR_SetKeepActive( option=.TRUE. )

       ! Start measuring KPP-related routine timing for this grid box
       IF ( State_Diag%Archive_KppTime ) THEN
          call cpu_time(TimeStart)
       ENDIF

       !=====================================================================
       ! Get photolysis rates (daytime only)
       !
       ! NOTE: The ordering of the photolysis reactions here is
       ! the order in the Fast-J definition file FJX_j2j.dat.
       ! I've assumed that these are the same as in the text files
       ! but this may have been changed.  This needs to be checked
       ! through more thoroughly -- M. Long (3/28/16)
       !
       ! ALSO NOTE: We moved this section above the test to see if grid
       ! box (I,J,L) is in the chemistry grid.  This will ensure that
       ! J-value diagnostics are defined for all levels in the column.
       ! This modification was validated by a geosfp_4x5_standard
       ! difference test. (bmy, 1/18/18)
       !
       ! Update SUNCOSmid threshold from 0 to cos(98 degrees) since
       ! fast-jx allows for SZA down to 98 degrees. This is important in
       ! the stratosphere-mesosphere where sunlight still illuminates at
       ! high altitudes if the sun is below the horizon at the surface
       ! (update submitted by E. Fleming (NASA), 10/11/2018)
       !=====================================================================
       IF ( State_Met%SUNCOSmid(I,J) > -0.1391731e+0_fp ) THEN

          ! Start timer
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( TimerName = "  -> Photolysis rates",          &
                               InLoop    = .TRUE.,                           &
                               ThreadNum = Thread,                           &
                               RC        = RC                               )
          ENDIF

          ! Get the fraction of H2SO4 that is available for photolysis
          ! (this is only valid for UCX-enabled mechanisms)
          SO4_FRAC = SO4_PHOTFRAC( I, J, L, State_Chm )

          ! Adjust certain photolysis rates:
          ! (1) H2SO4 + hv -> SO2 + OH + OH   (UCX-based mechanisms)
          ! (2) O3    + hv -> O2  + O         (UCX-based mechanisms)
          ! (2) O3    + hv -> OH  + OH        (trop-only mechanisms)
          CALL PHOTRATE_ADJ( Input_Opt, State_Diag, State_Met, I,            &
                             J,         L,          SO4_FRAC,  IERR         )

          ! Loop over the FAST-JX photolysis species
          DO N = 1, JVN_

             IF ( DO_PHOTCHEM ) THEN

                ! Copy photolysis rate from FAST_JX into KPP PHOTOL array
                PHOTOL(N) = ZPJ(L,N,I,J)

             ENDIF

             !===============================================================
             ! HISTORY (aka netCDF diagnostics)
             !
             ! Instantaneous photolysis rates [s-1] (aka J-values)
             ! and noontime photolysis rates [s-1]
             !
             !    NOTE: Attach diagnostics here instead of in module
             !    fast_jx_mod.F90 so that we can get the adjusted photolysis
             !    rates (output from routne PHOTRATE_ADJ above).
             !
             ! The mapping between the GEOS-Chem photolysis species and
             ! the FAST-JX photolysis species is contained in the lookup
             ! table in input file FJX_j2j.dat.

             ! Some GEOS-Chem photolysis species may have multiple
             ! branches for photolysis reactions.  These will be
             ! represented by multiple entries in the FJX_j2j.dat
             ! lookup table.
             !
             !    NOTE: For convenience, we have stored the GEOS-Chem
             !    photolysis species index (range: 1..State_Chm%nPhotol)
             !    for each of the FAST-JX photolysis species (range;
             !    1..JVN_) in the GC_PHOTO_ID array (located in module
             !    CMN_FJX_MOD.F90).
             !===============================================================

             ! GC photolysis species index
             P = GC_Photo_Id(N)

             ! If this FAST_JX photolysis species maps to a valid
             ! GEOS-Chem photolysis species (for this simulation)...
             IF ( P > 0 .and. P <= State_Chm%nPhotol ) THEN

                ! Archive the instantaneous photolysis rate
                ! (summing over all reaction branches)
                IF ( State_Diag%Archive_Jval ) THEN
                   S = State_Diag%Map_Jval%id2slot(P)
                   IF ( S > 0 ) THEN
                      State_Diag%Jval(I,J,L,S) =                             &
                      State_Diag%Jval(I,J,L,S) + PHOTOL(N)
                   ENDIF
                ENDIF

                ! Satellite diagnostics
                ! Archive the instantaneous photolysis rate
                ! (summing over all reaction branches)
                IF ( State_Diag%Archive_SatDiagnJval ) THEN
                   S = State_Diag%Map_SatDiagnJval%id2slot(P)
                   IF ( S > 0 ) THEN
                      State_Diag%SatDiagnJval(I,J,L,S) =                     &
                      State_Diag%SatDiagnJval(I,J,L,S) + PHOTOL(N)
                   ENDIF
                ENDIF

                ! Archive the noontime photolysis rate
                ! (summing over all reaction branches)
                IF ( State_Met%IsLocalNoon(I,J) ) THEN
                   IF ( State_Diag%Archive_JNoon ) THEN
                      S = State_Diag%Map_JNoon%id2slot(P)
                      IF ( S > 0 ) THEN
                         State_Diag%JNoon(I,J,L,S) =                         &
                         State_Diag%JNoon(I,J,L,S) + PHOTOL(N)
                      ENDIF
                   ENDIF
                ENDIF

             ELSE IF ( P == State_Chm%nPhotol+1 ) THEN

                ! J(O3_O1D).  This used to be stored as the nPhotol+1st
                ! diagnostic in Jval, but needed to be broken off
                ! to facilitate cleaner diagnostic indexing (bmy, 6/3/20)
                IF ( State_Diag%Archive_JvalO3O1D ) THEN
                   State_Diag%JvalO3O1D(I,J,L) =                             &
                   State_Diag%JvalO3O1D(I,J,L) + PHOTOL(N)
                ENDIF

                ! J(O3_O1D) for satellite diagnostics
                IF ( State_Diag%Archive_SatDiagnJvalO3O1D ) THEN
                   State_Diag%SatDiagnJvalO3O1D(I,J,L) =                     &
                   State_Diag%SatDiagnJvalO3O1D(I,J,L) + PHOTOL(N)
                ENDIF

             ELSE IF ( P == State_Chm%nPhotol+2 ) THEN

                ! J(O3_O3P).  This used to be stored as the nPhotol+2nd
                ! diagnostic in Jval, but needed to be broken off
                ! to facilitate cleaner diagnostic indexing (bmy, 6/3/20)
                IF ( State_Diag%Archive_JvalO3O3P ) THEN
                   State_Diag%JvalO3O3P(I,J,L) =                             &
                   State_Diag%JvalO3O3P(I,J,L) + PHOTOL(N)
                ENDIF

                ! J(O3_O3P) for satellite diagnostics
                IF ( State_Diag%Archive_SatDiagnJvalO3O3P ) THEN
                   State_Diag%SatDiagnJvalO3O3P(I,J,L) =                     &
                   State_Diag%SatDiagnJvalO3O3P(I,J,L) + PHOTOL(N)
                ENDIF

             ENDIF
          ENDDO

          ! Stop timer
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( TimerName = "  -> Photolysis rates",            &
                             InLoop    = .TRUE.,                             &
                             ThreadNum = Thread,                             &
                             RC        = RC                                 )
          ENDIF
       ENDIF

       !=====================================================================
       ! Test if we need to do the chemistry for box (I,J,L),
       ! otherwise move onto the next box.
       !=====================================================================

       ! If we are not in the troposphere don't do the chemistry!
       IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

       ! Skipping buffer zone (lzh, 08/10/2014)
       IF ( State_Grid%NestedGrid ) THEN
          IF ( J <=                 State_Grid%SouthBuffer ) CYCLE
          IF ( J >  State_Grid%NY - State_Grid%NorthBuffer ) CYCLE
          IF ( I <=                 State_Grid%EastBuffer  ) CYCLE
          IF ( I >  State_Grid%NX - State_Grid%WestBuffer  ) CYCLE
       ENDIF

       !=====================================================================
       ! Initialize the KPP "C" vector of species concentrations [molec/cm3]
       !=====================================================================
       DO N = 1, NSPEC
          SpcID = State_Chm%Map_KppSpc(N)
          C(N)  = 0.0_dp
          IF ( SpcId > 0 ) C(N) = State_Chm%Species(SpcID)%Conc(I,J,L)
       ENDDO

       !=====================================================================
       ! Start KPP main timer
       !=====================================================================
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start( TimerName = "  -> KPP",                          &
                            InLoop    = .TRUE.,                              &
                            ThreadNum = Thread,                              &
                            RC        = RC                                  )
       ENDIF

       !=====================================================================
       ! CHEMISTRY MECHANISM INITIALIZATION (#1)
       !
       ! Populate KPP global variables and arrays in gckpp_global.F90
       !
       ! NOTE: This has to be done before Set_Sulfur_Chem_Rates, so that
       ! the NUMDEN and SR_TEMP KPP variables will be populated first.
       ! Otherwise this can lead to differences in output that are evident
       ! when running with different numbers of OpenMP cores.
       ! See https://github.com/geoschem/geos-chem/issues/1157
       !    -- Bob Yantosca (08 Mar 2022)
       !=====================================================================

       ! Start timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start( TimerName = "  -> Init KPP",                     &
                            InLoop    = .TRUE.,                              &
                            ThreadNum = Thread,                              &
                            RC        =  RC                                 )
       ENDIF

       ! Copy values into the various KPP global variables
       CALL Set_Kpp_GridBox_Values( I          = I,                          &
                                    J          = J,                          &
                                    L          = L,                          &
                                    Input_Opt  = Input_Opt,                  &
                                    State_Chm  = State_Chm,                  &
                                    State_Grid = State_Grid,                 &
                                    State_Met  = State_Met,                  &
                                    RC         = RC                         )

       ! Stop timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( TimerName =  "  -> Init KPP",                      &
                          InLoop    = .TRUE.,                                &
                          ThreadNum = Thread,                                &
                          RC        =  RC                                   )
       ENDIF

       !=====================================================================
       ! CHEMISTRY MECHANISM INITIALIZATION (#2)
       !
       ! Update reaction rates [1/s] for sulfur chemistry in cloud and on
       ! seasalt.  These will be passed to the KPP chemical solver.
       !
       ! NOTE: This has to be done before fullchem_SetStateHet so that
       ! State_Chm%HSO3_aq and State_Chm%SO3_aq will be populated first.
       ! These are copied into State_Het%HSO3_aq and State_Het%SO3_aq.
       ! See https://github.com/geoschem/geos-chem/issues/1157
       !    -- Bob Yantosca (08 Mar 2022)
       !=====================================================================

       ! Start timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start( TimerName = "     RCONST",                       &
                            InLoop    = .TRUE.,                              &
                            ThreadNum = Thread,                              &
                            RC        = RC                                  )
       ENDIF

       ! Compute sulfur chemistry reaction rates [1/s]
       ! If size_res = T, we'll call fullchem_HetDropChem below.
       CALL Set_Sulfur_Chem_Rates( I          = I,                           &
                                   J          = J,                           &
                                   L          = L,                           &
                                   Input_Opt  = Input_Opt,                   &
                                   State_Chm  = State_Chm,                   &
                                   State_Diag = State_Diag,                  &
                                   State_Grid = State_Grid,                  &
                                   State_Met  = State_Met,                   &
                                   size_res   = size_res,                    &
                                   RC         = RC                          )

       ! Stop timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( TimerName = "     RCONST",                         &
                          InLoop    = .TRUE.,                                &
                          ThreadNum = Thread,                                &
                          RC        = RC                                    )
       ENDIF

       !=====================================================================
       ! CHEMISTRY MECHANISM INITIALIZATION (#3)
       !
       ! Populate the various fields of the State_Het object.
       !
       ! NOTE: This has to be done after fullchem_SetStateHet so that
       ! State_Chm%HSO3_aq and State_Chm%SO3_aq will be populated first.
       ! These are copied into State_Het%HSO3_aq and State_Het%SO3_aq.
       ! See https://github.com/geoschem/geos-chem/issues/1157
       !    -- Bob Yantosca (08 Mar 2022)
       !=====================================================================

       ! Start timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start( TimerName = "  -> Init KPP",                     &
                            InLoop    = .TRUE.,                              &
                            ThreadNum = Thread,                              &
                            RC        =  RC                                 )
       ENDIF

       ! Populate fields of the State_Het object
       CALL fullchem_SetStateHet( I         = I,                             &
                                  J         = J,                             &
                                  L         = L,                             &
                                  Input_Opt = Input_Opt,                     &
                                  State_Chm = State_Chm,                     &
                                  State_Met = State_Met,                     &
                                  H         = State_Het,                     &
                                  RC        = RC                            )

       ! Stop timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( TimerName =  "  -> Init KPP",                      &
                          InLoop    = .TRUE.,                                &
                          ThreadNum = Thread,                                &
                          RC        =  RC                                   )
       ENDIF

       !=====================================================================
       ! CHEMISTRY MECHANISM INITIALIZATION (#5)
       !
       ! Call Het_Drop_Chem (formerly located in sulfate_mod.F90) to
       ! estimate the in-cloud sulfate production rate in heterogeneous
       ! cloud droplets based on the Yuen et al., 1996 parameterization.
       ! Code by Becky Alexander (2011) with updates by Mike Long and Bob
       ! Yantosca (2021).
       !
       ! We will only call Het_Drop_Chem if:
       ! (1) It is at least 0.01% cloudy
       ! (2) We are doing a size-resolved computation
       ! (3) The grid box is over water
       ! (4) The temperature is above -5C
       ! (5) Liquid water content is nonzero
       !=====================================================================
       IF ( State_Met%CLDF(I,J,L) > 1.0e-4_fp ) THEN

          ! Start timer
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( TimerName = "     RCONST",                    &
                               InLoop    = .TRUE.,                           &
                               ThreadNum = Thread,                           &
                               RC        =  RC                              )
          ENDIF

          ! Liquid water content (same formula from the old sulfate_mod.F90)
          LWC = ( State_Met%QL(I,J,L) * State_Met%AIRDEN(I,J,L)              &
              *   1.0e-3_fp           / State_Met%CLDF(I,J,L)               )


          ! Eexecute fullchem_HetDropChem if criteria are satisfied
          ! NOTE: skip if LWC is very small, which will blow up equations!
          IF ( ( size_res                                           )  .and. &
               ( State_Met%IsWater(I,J)                             )  .and. &
               ( TEMP                    > 268.15_fp                )  .and. &
               ( LWC                     > 1.0e-20_fp               ) ) THEN

             CALL fullchem_HetDropChem( I         = I,                       &
                                        J         = J,                       &
                                        L         = L,                       &
                                        SR        = SR,                      &
                                        Input_Opt = Input_Opt,               &
                                        State_Met = State_Met,               &
                                        State_Chm = State_Chm               )

             ! Add result as an enhancement to O2 metal catalysis rate
             ! as a 1st order reaction
             K_CLD(3) = K_CLD(3) + SR
          ENDIF

          ! Stop timer
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( TimerName = "     RCONST",                       &
                             InLoop    = .TRUE.,                              &
                             ThreadNum = Thread,                              &
                             RC        =  RC                                 )
          ENDIF
       ENDIF

       !=====================================================================
       ! Start KPP main timer and prepare arrays
       !=====================================================================

       ! Zero out dummy species index in KPP
       DO F = 1, NFAM
          KppID = PL_Kpp_Id(F)
          IF ( KppID > 0 ) C(KppID) = 0.0_dp
       ENDDO

       !=====================================================================
       ! Update reaction rates
       !=====================================================================

       ! Start timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start( TimerName = "     RCONST",                       &
                            InLoop    = .TRUE.,                              &
                            ThreadNum = Thread,                              &
                            RC        =  RC                                 )
       ENDIF

       ! Update the array of rate constants
       CALL Update_RCONST( )

       ! Stop timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( TimerName = "     RCONST",                         &
                          InLoop    = .TRUE.,                                &
                          ThreadNum = Thread,                                &
                          RC        =  RC                                   )
       ENDIF

       !=====================================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Archive KPP reaction rates [s-1]
       ! See gckpp_Monitor.F90 for a list of chemical reactions
       !
       ! NOTE: In KPP 2.5.0+, VAR and FIX are now private to the integrator
       ! and point to C.  Therefore, pass C(1:NVAR) instead of VAR and
       ! C(NVAR+1:NSPEC) instead of FIX to routine FUN.
       !=====================================================================
       IF ( State_Diag%Archive_RxnRate                                  .or. &
            State_Diag%Archive_SatDiagnRxnRate                        ) THEN
  
          ! Get equation rates (Aout)
          CALL Fun( V       = C(1:NVAR),                                     &
                    F       = C(NVAR+1:NSPEC),                               &
                    RCT     = RCONST,                                        &
                    Vdot    = Vloc,                                          &
                    Aout    = Aout                                          )

          ! Archive the RxnRate diagnostic collection
          IF ( State_Diag%Archive_RxnRate ) THEN
             DO S = 1, State_Diag%Map_RxnRate%nSlots
                N = State_Diag%Map_RxnRate%slot2Id(S)
                State_Diag%RxnRate(I,J,L,S) = Aout(N)
             ENDDO
          ENDIF

          ! Archive the SatDiagnRxnRate diagnostic collection
          IF ( State_Diag%Archive_SatDiagnRxnRate ) THEN
             DO S = 1, State_Diag%Map_SatDiagnRxnRate%nSlots
                N = State_Diag%Map_SatDiagnRxnRate%slot2Id(S)
                State_Diag%SatDiagnRxnRate(I,J,L,S) = Aout(N)
             ENDDO
          ENDIF
       ENDIF

       !=====================================================================
       ! Set options for the KPP integrator in vectors ICNTRL and RCNTRL
       ! This now needs to be done within the parallel loop
       !=====================================================================
       CALL fullchem_AR_SetIntegratorOptions( Input_Opt, State_Chm,          &
                                              State_Met, FirstChem,          &
                                              I,         J,         L,       &
                                              ICNTRL,    RCNTRL             )

       !=====================================================================
       ! Integrate the box forwards
       !=====================================================================

       ! Start timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start( TimerName = "     Integrate 1",                  &
                            InLoop    =  .TRUE.,                             &
                            ThreadNum = Thread,                              &
                            RC        = RC                                  )
       ENDIF

       ! Call the KPP integrator
       CALL Integrate( TIN,    TOUT,    ICNTRL,                              &
                       RCNTRL, ISTATUS, RSTATE, IERR                        )

       ! Stop timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( TimerName = "     Integrate 1",                    &
                          InLoop    = .TRUE.,                                &
                          ThreadNum = Thread,                                &
                          RC        = RC                                    )
       ENDIF

       ! Print grid box indices to screen if integrate failed
       IF ( IERR < 0 ) THEN
          WRITE(6,*) '### INTEGRATE RETURNED ERROR AT: ', I, J, L
       ENDIF

#if defined( MODEL_GEOS ) || defined( MODEL_WRF )
       ! Print grid box indices to screen if integrate failed
       IF ( IERR < 0 ) THEN
          WRITE(6,*) '### INTEGRATE RETURNED ERROR AT: ', I, J, L
          IF ( ASSOCIATED(State_Diag%KppError) ) THEN
             State_Diag%KppError(I,J,L) = State_Diag%KppError(I,J,L) + 1.0
          ENDIF
       ENDIF
#endif

       !=====================================================================
       ! HISTORY: Archive KPP solver diagnostics
       !=====================================================================
       IF ( State_Diag%Archive_KppDiags ) THEN

          ! # of integrator calls
          IF ( State_Diag%Archive_KppIntCounts ) THEN
             State_Diag%KppIntCounts(I,J,L) = ISTATUS(1)
          ENDIF

          ! # of times Jacobian was constructed
          IF ( State_Diag%Archive_KppJacCounts ) THEN
             State_Diag%KppJacCounts(I,J,L) = ISTATUS(2)
          ENDIF

          ! # of internal timesteps
          IF ( State_Diag%Archive_KppTotSteps ) THEN
             State_Diag%KppTotSteps(I,J,L) = ISTATUS(3)
          ENDIF

          ! # of accepted internal timesteps
          IF ( State_Diag%Archive_KppAccSteps ) THEN
             State_Diag%KppAccSteps(I,J,L) = ISTATUS(4)
          ENDIF

          ! # of rejected internal timesteps
          IF ( State_Diag%Archive_KppRejSteps ) THEN
             State_Diag%KppRejSteps(I,J,L) = ISTATUS(5)
          ENDIF

          ! # of LU-decompositions
          IF ( State_Diag%Archive_KppLuDecomps ) THEN
             State_Diag%KppLuDecomps(I,J,L) = ISTATUS(6)
          ENDIF

          ! # of forward and backwards substitutions
          IF ( State_Diag%Archive_KppSubsts ) THEN
             State_Diag%KppSubsts(I,J,L) = ISTATUS(7)
          ENDIF

          ! # of singular-matrix decompositions
          IF ( State_Diag%Archive_KppSmDecomps ) THEN
             State_Diag%KppSmDecomps(I,J,L) = ISTATUS(8)
          ENDIF

          ! Update autoreduce solver statistics
          ! (only if the autoreduction is turned on)
          IF ( Input_Opt%Use_AutoReduce ) THEN
             CALL fullchem_AR_UpdateKppDiags( I, J, L, RSTATE, State_Diag )
          ENDIF
       ENDIF

       !=====================================================================
       ! Try another time if it failed
       !=====================================================================
       IF ( IERR < 0 ) THEN

#if defined( MODEL_GEOS ) || defined( MODEL_WRF )
          ! Save a copy of the C vector (GEOS and WRF only)
          localC    = C
#endif

          ! Reset first time step and start concentrations
          ! Retry the integration with non-optimized settings
          RCNTRL(3) = 0.0_dp
          C         = 0.0_dp

          ! Disable auto-reduce solver for the second iteration for safety
          IF ( Input_Opt%Use_AutoReduce ) THEN
             RCNTRL(12) = -1.0_dp ! without using ICNTRL
          ENDIF

          ! Update rates again
          CALL Update_RCONST( )

          ! Start timer
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( TimerName = "     Integrate 2",               &
                               InLoop    =  .TRUE.,                          &
                               ThreadNum = Thread,                           &
                               RC        = RC                               )
          ENDIF

          ! Integrate again
          CALL Integrate( TIN,    TOUT,    ICNTRL,                           &
                          RCNTRL, ISTATUS, RSTATE, IERR                     )

          ! Stop timer
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( TimerName = "     Integrate 2",                 &
                             InLoop    =  .TRUE.,                            &
                             ThreadNum = Thread,                             &
                             RC        = RC                                 )
          ENDIF

          !==================================================================
          ! HISTORY: Archive KPP solver diagnostics
          ! This time, add to the existing value
          !==================================================================
          IF ( State_Diag%Archive_KppDiags ) THEN

             ! # of integrator calls
             IF ( State_Diag%Archive_KppIntCounts ) THEN
                State_Diag%KppIntCounts(I,J,L) =                             &
                State_Diag%KppIntCounts(I,J,L) + ISTATUS(1)
             ENDIF

             ! # of times Jacobian was constructed
             IF ( State_Diag%Archive_KppJacCounts ) THEN
                State_Diag%KppJacCounts(I,J,L) =                             &
                State_Diag%KppJacCounts(I,J,L) + ISTATUS(2)
             ENDIF

             ! # of internal timesteps
             IF ( State_Diag%Archive_KppTotSteps ) THEN
                State_Diag%KppTotSteps(I,J,L) =                              &
                State_Diag%KppTotSteps(I,J,L) + ISTATUS(3)
             ENDIF

             ! # of accepted internal timesteps
             IF ( State_Diag%Archive_KppAccSteps ) THEN
                State_Diag%KppAccSteps(I,J,L) =                              &
                State_Diag%KppAccSteps(I,J,L) + ISTATUS(4)
             ENDIF

             ! # of rejected internal timesteps
             IF ( State_Diag%Archive_KppRejSteps ) THEN
                State_Diag%KppRejSteps(I,J,L) =                              &
                State_Diag%KppRejSteps(I,J,L) + ISTATUS(5)
             ENDIF

             ! # of LU-decompositions
             IF ( State_Diag%Archive_KppLuDecomps ) THEN
                State_Diag%KppLuDecomps(I,J,L) =                             &
                State_Diag%KppLuDecomps(I,J,L) + ISTATUS(6)
             ENDIF

             ! # of forward and backwards substitutions
             IF ( State_Diag%Archive_KppSubsts ) THEN
                State_Diag%KppSubsts(I,J,L) =                                &
                State_Diag%KppSubsts(I,J,L) + ISTATUS(7)
             ENDIF

             ! # of singular-matrix decompositions
!             IF ( State_Diag%Archive_KppSmDecomps ) THEN
!                State_Diag%KppSmDecomps(I,J,L) =                             &
!                State_Diag%KppSmDecomps(I,J,L) + ISTATUS(8)
!             ENDIF
          ENDIF

          !==================================================================
          ! Exit upon the second failure
          !==================================================================
          IF ( IERR < 0 ) THEN
             WRITE(6,     '(a   )' ) '## INTEGRATE FAILED TWICE !!! '
             WRITE(ERRMSG,'(a,i3)' ) 'Integrator error code :', IERR
#if defined( MODEL_GEOS ) || defined( MODEL_WRF )
             IF ( Input_Opt%KppStop ) THEN
                CALL ERROR_STOP(ERRMSG, 'INTEGRATE_KPP')
             ! Revert to start values
             ELSE
                C = localC
             ENDIF
             IF ( ASSOCIATED(State_Diag%KppError) ) THEN
                State_Diag%KppError(I,J,L) = State_Diag%KppError(I,J,L) + 1.0
             ENDIF
#else
             ! Set a flag to break out of loop gracefully
             ! NOTE: You can set a GDB breakpoint here to examine the error
             !$OMP CRITICAL
             Failed2x = .TRUE.

             ! Print concentrations at trouble box KPP error
             PRINT*, REPEAT( '###', 79 )
             PRINT*, '### KPP DEBUG OUTPUT!'
             PRINT*, '### Species concentrations at problem box ', I, J, L
             PRINT*, REPEAT( '###', 79 )
             DO N = 1, NSPEC
                PRINT*, '### ', C(N), TRIM( ADJUSTL( SPC_NAMES(N) ) )
             ENDDO

             ! Print rate constants at trouble box KPP error
             PRINT*, REPEAT( '###', 79 )
             PRINT*, '### KPP DEBUG OUTPUT!'
             PRINT*, '### Reaction rates at problem box ', I, J, L
             PRINT*, REPEAT( '###', 79 )
             DO N = 1, NREACT
                PRINT*, RCONST(N), TRIM( ADJUSTL( EQN_NAMES(N) ) )
             ENDDO
             !$OMP END CRITICAL
#endif
          ENDIF

       ENDIF


       !=====================================================================
       ! Continue upon successful return...
       !=====================================================================

       ! Revert Alkalinity (only when using sulfur chemistry in KPP)
       IF ( .not. State_Chm%Do_SulfateMod_SeaSalt ) THEN
          CALL fullchem_ConvertEquivToAlk()
       ENDIF

       ! Save Hnew (the last predicted but not taken step) from the 3rd slot
       ! of RSTATE into State_Chm so that it can be written to the restart
       ! file.  For simulations that are broken into multiple stages,
       ! Hstart will be initialized to the value of Hnew from the restart
       ! file at startup (see above).
       State_Chm%KPPHvalue(I,J,L) = RSTATE(Nhnew)

       ! Save cpu time spent for bulk of KPP-related routines for History archival
       ! (hplin, 11/8/21)
       IF ( State_Diag%Archive_KppTime ) THEN
          call cpu_time(TimeEnd)
          State_Diag%KppTime(I,J,L) = TimeEnd - TimeStart
       ENDIF

       !=====================================================================
       ! Check we have no negative values and copy the concentrations
       ! calculated from the C array back into State_Chm%Species%Conc
       !=====================================================================

       ! Loop over KPP species
       DO N = 1, NSPEC

          ! GEOS-Chem species ID
          SpcID = State_Chm%Map_KppSpc(N)

          ! Skip if this is not a GEOS-Chem species
          IF ( SpcID <= 0 ) CYCLE

          ! Set negative concentrations to zero
          C(N) = MAX( C(N), 0.0_dp )

          ! Copy concentrations back into State_Chm%Species
          State_Chm%Species(SpcID)%Conc(I,J,L) = REAL( C(N), kind=fp )

       ENDDO

       IF ( Input_Opt%useTimers ) THEN

          ! Stop main KPP timer
          CALL Timer_End( TimerName  = "  -> KPP",                           &
                          InLoop     = .TRUE.,                               &
                          ThreadNum  = Thread,                               &
                          RC         = RC                                   )

          ! Start Prod/Loss timer
          CALL Timer_Start( TimerName = "  -> Prod/loss diags",              &
                            InLoop    = .TRUE.,                              &
                            ThreadNum = Thread,                              &
                            RC        = RC                                  )
       ENDIF

#ifdef BPCH_DIAG
#ifdef TOMAS
       !always calculate rate for TOMAS
       DO F = 1, NFAM

          ! Determine dummy species index in KPP
          KppID =  PL_Kpp_Id(F)

          !-----------------------------------------------------------------
          ! FOR TOMAS MICROPHYSICS:
          !
          ! Obtain P/L with a unit [kg S] for tracing
          ! gas-phase sulfur species production (SO2, SO4, MSA)
          ! (win, 8/4/09)
          !-----------------------------------------------------------------

          ! Calculate H2SO4 production rate [kg s-1] in each
          ! time step (win, 8/4/09)
          IF ( TRIM(FAM_NAMES(F)) == 'PSO4' ) THEN
             ! Hard-coded MW
             H2SO4_RATE(I,J,L) = C(KppID) / AVO * 98.e-3_fp * &
                                 State_Met%AIRVOL(I,J,L)    * &
                                 1.0e+6_fp / DT

            IF ( H2SO4_RATE(I,J,L) < 0.0d0) THEN
              write(*,*) "H2SO4_RATE negative in fullchem_mod.F90!!", &
                 I, J, L, "was:", H2SO4_RATE(I,J,L), "  setting to 0.0d0"
              H2SO4_RATE(I,J,L) = 0.0d0
            ENDIF
          ENDIF
       ENDDO

#endif
#endif

#ifdef MODEL_CESM
       ! Calculate H2SO4 production rate for coupling to CESM (interface to MAM4 nucleation)
       DO F = 1, NFAM

          ! Determine dummy species index in KPP
          KppID =  PL_Kpp_Id(F)

          ! Calculate H2SO4 production rate [mol mol-1] in this timestep (hplin, 1/25/23)
          IF ( TRIM(FAM_NAMES(F)) == 'PSO4' ) THEN
             ! mol/mol = molec cm-3 * g * mol(Air)-1 * kg g-1 * m-3 cm3 / (molec mol-1 * kg m-3)
             !         = mol/molAir
             State_Chm%H2SO4_PRDR(I,J,L) = C(KppID) * AIRMW * 1e-3_fp * 1.0e+6_fp / &
                                 (AVO * State_Met%AIRDEN(I,J,L))

             IF ( State_Chm%H2SO4_PRDR(I,J,L) < 0.0d0) THEN
               write(*,*) "H2SO4_PRDR negative in fullchem_mod.F90!!", &
                  I, J, L, "was:", State_Chm%H2SO4_PRDR(I,J,L), "  setting to 0.0d0"
               State_Chm%H2SO4_PRDR(I,J,L) = 0.0d0
             ENDIF
          ENDIF
       ENDDO
#endif

#ifdef MODEL_GEOS
       !--------------------------------------------------------------------
       ! Archive NOx lifetime [h]
       !--------------------------------------------------------------------
       IF ( State_Diag%Archive_NoxTau .OR. State_Diag%Archive_TropNOxTau ) THEN
          CALL Fun( V       = C(1:NVAR),                                     &
                    F       = C(NVAR+1:NSPEC),                               &
                    RCT     = RCONST,                                        &
                    Vdot    = Vloc,                                          &
                    Aout    = Aout                                          )
          NOxTau = Vloc(ind_NO) + Vloc(ind_NO2) + Vloc(ind_NO3)         &
                 + 2.*Vloc(ind_N2O5) + Vloc(ind_ClNO2) + Vloc(ind_HNO2) &
                 + Vloc(ind_HNO4)
          NOxConc = C(ind_NO) + C(ind_NO2) + C(ind_NO3) + 2.*C(ind_N2O5)         &
                  + C(ind_ClNO2) + C(ind_HNO2) + C(ind_HNO4)
          ! NOx chemical lifetime per grid cell
          IF ( State_Diag%Archive_NoxTau ) THEN
             NoxTau = ( NOxConc / (-1.0_f4*NOxTau) ) / 3600.0_f4
             IF ( NoxTau > 0.0_f4 ) THEN
                State_Diag%NOxTau(I,J,L) = min(1.0e10_f4,max(1.0e-10_f4,NOxTau))
             ELSE
                State_Diag%NOxTau(I,J,L) = max(-1.0e10_f4,min(-1.0e-10_f4,NOxTau))
             ENDIF
          ENDIF
          ! NOx chemical lifetime per trop. column
          IF ( State_Diag%Archive_TropNOxTau ) THEN
             NOx_weight = ( NOxConc )*State_Met%AIRDEN(I,J,L)*State_Met%DELP_DRY(I,J,L)
             NOx_tau_weighted = ( NOxConc / ( -1.0_f4*NOxTau*3600.0_f4 ) )*NOx_weight
             IF ( ABS(NOx_tau_weighted) < 1.0e8 ) THEN
               NOx_tau_weighted = ( NINT(NOx_tau_weighted)*1.0e6 )*1.0e-6_f4
             ELSE
                IF ( NOx_tau_weighted > 0.0 ) THEN
                   NOx_tau_weighted = 1.0e8
                ELSE
                   NOx_tau_weighted = -1.0e8
                ENDIF
             ENDIF
             IF ( State_Met%InTroposphere(I,J,L) ) THEN
               TROPv_NOx_mass(I,J) = TROPv_NOx_mass(I,J) + NOx_weight
               TROPv_NOx_tau(I,J)  = TROPv_NOx_tau(I,J) + NOx_tau_weighted
             ENDIF
          ENDIF
       ENDIF
#endif

       !====================================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Prod and loss of families or species [molec/cm3/s]
       !
       ! NOTE: KppId is the KPP ID # for each of the prod and loss
       ! diagnostic species.  This is the value used to index the
       ! KPP "C" array (in module gckpp_Global.F90).
       !====================================================================

       ! Chemical loss of species or families [molec/cm3/s]
       IF ( State_Diag%Archive_Loss ) THEN
          DO S = 1, State_Diag%Map_Loss%nSlots
             KppId = State_Diag%Map_Loss%slot2Id(S)
             State_Diag%Loss(I,J,L,S) = C(KppID) / DT
          ENDDO
       ENDIF

       ! Chemical production of species or families [molec/cm3/s]
       IF ( State_Diag%Archive_Prod ) THEN
          DO S = 1, State_Diag%Map_Prod%nSlots
             KppID = State_Diag%Map_Prod%slot2Id(S)
             State_Diag%Prod(I,J,L,S) = C(KppID) / DT
          ENDDO
       ENDIF

       ! Satellite diagnostic: Chemical loss [molec/cm3/s]
       IF ( State_Diag%Archive_SatDiagnLoss ) THEN
          DO S = 1, State_Diag%Map_SatDiagnLoss%nSlots
             KppId = State_Diag%Map_SatDiagnLoss%slot2Id(S)
             State_Diag%SatDiagnLoss(I,J,L,S) = C(KppID) / DT
          ENDDO
       ENDIF

       ! Satellite diagnostic: Chemical production [molec/cm3/s]
       IF ( State_Diag%Archive_SatDiagnProd ) THEN
          DO S = 1, State_Diag%Map_SatDiagnProd%nSlots
             KppID = State_Diag%Map_SatDiagnProd%slot2Id(S)
             State_Diag%SatDiagnProd(I,J,L,S) = C(KppID) / DT
          ENDDO
       ENDIF

       !--------------------------------------------------------------------
       ! Archive prod/loss fields for the TagCO simulation [molec/cm3/s]
       ! (In practice, we only need to do this from benchmark simulations)
       !--------------------------------------------------------------------
       IF ( State_Diag%Archive_ProdCOfromCH4     .or.                        &
            State_Diag%Archive_ProdCOfromNMVOC ) THEN

          ! Total production of CO
          PCO_TOT   = C(id_PCO) / DT

          ! Loss of CO from CH4
          LCH4      = C(id_LCH4) / DT

          ! P(CO)_CH4 is LCH4. Cap so that it is never greater
          ! than total P(CO) to prevent negative P(CO)_NMVOC.
          PCO_CH4   = MIN( LCH4, PCO_TOT )

          ! P(CO) from NMVOC is the remaining P(CO)
          PCO_NMVOC = PCO_TOT - PCO_CH4

          ! Archive P(CO) from CH4 for tagCO simulations
          IF ( State_Diag%Archive_ProdCOfromCH4 ) THEN
             State_Diag%ProdCOfromCH4(I,J,L) = PCO_CH4
          ENDIF

          ! Archive P(CO) from NMVOC for tagCO simulations
          IF ( State_Diag%Archive_ProdCOfromNMVOC ) THEN
             State_Diag%ProdCOfromNMVOC(I,J,L) = PCO_NMVOC
          ENDIF

       ENDIF

       ! Stop Prod/Loss timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( TimerName =  "  -> Prod/loss diags",               &
                          InLoop    = .TRUE.,                                &
                          ThreadNum = Thread,                                &
                          RC        = RC                                    )
       ENDIF

       !====================================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Write out OH reactivity.  The OH reactivity is defined here as the
       ! inverse of its life-time. In a crude ad-hoc approach, manually add
       ! all OH reactants (ckeller, 9/20/2017)
       !====================================================================
       IF ( State_Diag%Archive_OHreactivity           .or.                   &
            State_Diag%Archive_SatDiagnOHreactivity ) THEN

          ! Start timer
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( TimerName = "  -> OH reactivity diag",        &
                               InLoop    = .TRUE.,                           &
                               ThreadNum = Thread,                           &
                               RC        = RC                               )
          ENDIF

          ! Archive OH reactivity diagnostic
          CALL Get_OHreactivity ( C, RCONST, OHreact )
          IF ( State_Diag%Archive_OHreactivity ) THEN
             State_Diag%OHreactivity(I,J,L) = OHreact
          ENDIF
          IF ( State_Diag%Archive_SatDiagnOHreactivity ) THEN
             State_Diag%SatDiagnOHreactivity(I,J,L) = OHreact
          ENDIF

          ! Stop timer
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( TimerName = "  -> OH reactivity diag",          &
                             InLoop    = .TRUE.,                             &
                             ThreadNum = Thread,                             &
                             RC        = RC                                 )
          ENDIF
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Stop timer
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End( "  -> FlexChem loop", RC )
    ENDIF

    ! Compute sum of in-loop timers
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_Sum_Loop( "  -> Init KPP",            RC )
       CALL Timer_Sum_Loop( "  -> Het chem rates",      RC )
       CALL Timer_Sum_Loop( "  -> Photolysis rates",    RC )
       CALL Timer_Sum_Loop( "  -> KPP",                 RC )
       CALL Timer_Sum_Loop( "     RCONST",              RC )
       CALL Timer_Sum_Loop( "     Integrate 1",         RC )
       CALL Timer_Sum_Loop( "     Integrate 2",         RC )
       CALL Timer_Sum_Loop( "  -> Prod/loss diags",     RC )
       CALL Timer_Sum_Loop( "  -> OH reactivity diag",  RC )
    ENDIF

    !=======================================================================
    ! Return gracefully if integration failed 2x anywhere
    ! (as we cannot break out of a parallel DO loop!)
    !=======================================================================
    IF ( Failed2x ) THEN
       ErrMsg = 'KPP failed to converge after 2 iterations!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

#if defined( MODEL_GEOS )
    IF ( State_Diag%Archive_TropNOxTau ) THEN
       WHERE (TROPv_NOx_mass > 0.0_f4 )
          State_Diag%TropNOxTau = TROPv_NOx_tau / TROPv_NOx_mass
       END WHERE
    ENDIF
#endif

    !=======================================================================
    ! Archive OH, HO2, O1D, O3P concentrations after solver
    !=======================================================================
    IF ( Do_Diag_OH_HO2_O1D_O3P ) THEN

       ! Save OH, HO2, O1D, O3P for the ND43 diagnostic
       ! NOTE: These might not be needed for netCDF, as they will already
       ! have been archived in State_Chm%Species%Conc output.
       CALL Diag_OH_HO2_O1D_O3P( Input_Opt,  State_Chm, State_Diag, &
                                 State_Grid, State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Diag_OH_HO2_O1D_O3P"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### Do_FullChem: after OHSAVE' )
       ENDIF
    ENDIF

    !=======================================================================
    ! Archive quantities for computing OH metrics
    !=======================================================================
    CALL Diag_Metrics( Input_Opt,  State_Chm, State_Diag,                    &
                       State_Grid, State_Met, RC                            )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Diag_Mean_OH_and_CH4'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### Do_FullChem: after Diag_Metrics' )
    ENDIF

    !=======================================================================
    ! Archive concentrations after chemistry
    !=======================================================================
    IF ( State_Diag%Archive_ConcAfterChem ) THEN
       ! Point to mapping obj specific to ConcAfterChem diagnostic collection
       mapData => State_Diag%Map_ConcAfterChem

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( N, S   )
       DO S = 1, mapData%nSlots
          N = mapData%slot2id(S)
          State_Diag%ConcAfterChem(:,:,:,S) = State_Chm%Species(N)%Conc(:,:,:)
       ENDDO
       !$OMP END PARALLEL DO

       ! Free pointer
       mapData => NULL()
    ENDIF

    !=======================================================================
    ! Convert species back to original units (ewl, 8/16/16)
    !=======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm,  State_Grid, State_Met, &
                            OrigUnit,  RC )

    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'fullchem_mod.F90' )
       RETURN
    ENDIF

    !=======================================================================
    ! Apply high-altitude active nitrogen partitioning and H2SO4
    ! photolysis approximations outside the chemistry grid
    !=======================================================================
    CALL UCX_NOX( Input_Opt, State_Chm, State_Grid, State_Met )
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### CHEMDR: after UCX_NOX' )
    ENDIF

    CALL UCX_H2SO4PHOT( Input_Opt, State_Chm, State_Grid, State_Met )
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### CHEMDR: after UCX_H2SO4PHOT' )
    ENDIF

    ! Set State_Chm arrays for surface J-values used in HEMCO and
    ! saved to restart file
    IF ( RXN_O3_1 >= 0 ) THEN
       State_Chm%JOH(:,:) = ZPJ(1,RXN_O3_1,:,:)
    ENDIF
    IF ( RXN_NO2 >= 0 ) THEN
       State_Chm%JNO2(:,:) = ZPJ(1,RXN_NO2,:,:)
    ENDIF

    ! Set FIRSTCHEM = .FALSE. -- we have gone thru one chem step
    FIRSTCHEM = .FALSE.

  END SUBROUTINE Do_FullChem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PrintFirstTimeInfo
!
! !DESCRIPTION: Prints information the first time DO_FULLCHEM is called
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE PrintFirstTimeInfo( Input_Opt, State_Chm, FirstChem, Do_PhotChem )
!
! !USES:
!
    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN) :: Input_Opt     ! Input Options object
    TYPE(ChmState), INTENT(IN) :: State_Chm     ! Chemistry State object
    LOGICAL,        INTENT(IN) :: FirstChem     ! Is this the first call?
    LOGICAL,        INTENT(IN) :: Do_PhotChem   ! Is photolysis turned on?
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Exit if we are not on the first chemistry timestep
    IF ( .not. FirstChem ) RETURN

    ! Print on the root CPU only
    IF ( Input_Opt%AmIRoot ) THEN

       ! Write header
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
       WRITE( 6, 100 )
 100   FORMAT('DO_FULLCHEM : Interface between GEOS-Chem and the solver')
       WRITE( 6, '(a)' )

       ! Print where sulfur seasalt chemistry is done
       IF ( State_Chm%Do_SulfateMod_Seasalt ) THEN
          WRITE( 6, 110 )
 110      FORMAT( '* Sulfur sea salt chemistry is computed in sulfate_mod' )
       ELSE
          WRITE( 6, 120 )
 120      FORMAT( '* Sulfur sea salt chemistry is computed in KPP' )
       ENDIF

       ! Print where sulfur in-cloud chemistry is done
       IF ( State_Chm%Do_SulfateMod_Cld ) THEN
          WRITE( 6, 130 )
 130      FORMAT( '* Sulfur in-cloud chemistry is computed in sulfate_mod' )
       ELSE
          WRITE( 6, 140 )
 140      FORMAT( '* Sulfur in-cloud chemistry is computed in KPP' )
       ENDIF

       ! Print the status of photolysis: on or off
       IF ( Do_Photchem ) THEN
          WRITE( 6, 150 )
 150      FORMAT(  '* Photolysis is activated -- rates computed by FAST-JX' )
       ELSE
          WRITE( 6, 160 )
 160      FORMAT(  '* Photolysis has been deactivated for testing purposes'  )
       ENDIF

       ! Write footer
       WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    ENDIF

  END SUBROUTINE PrintFirstTimeInfo
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_sulfur_chem_rates
!
! !DESCRIPTION: Calls functions from the KPP rate-law library to compute
!  rates for sulfate chemistry reactions.  These are passed to the KPP
!  mechanism.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Sulfur_Chem_Rates( I,          J,          L,               &
                                    Input_Opt,  State_Chm,  State_Diag,      &
                                    State_Grid, State_Met,  size_res,        &
                                    RC                                      )
!
! !USES:
!
    USE ErrCode_Mod
    USE GcKpp_Global
    USE GcKpp_Parameters
    USE fullchem_SulfurChemFuncs, ONLY : fullchem_ConvertAlktoEquiv
    USE fullchem_SulfurChemFuncs, ONLY : fullchem_SulfurAqChem
    USE fullchem_SulfurChemFuncs, ONLY : fullchem_SulfurCldChem
    USE Input_Opt_Mod,            ONLY : OptInput
    USE State_Chm_Mod,            ONLY : ChmState
    USE State_Diag_Mod,           ONLY : DgnState
    USE State_Grid_Mod,           ONLY : GrdState
    USE State_Met_Mod,            ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L      ! X, Y, Z gridbox indices
    TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,        INTENT(OUT)   :: size_res     ! Should we call HetDropChem?
    INTEGER,        INTENT(OUT)   :: RC           ! Success or failure!
!
! !REMARKS:
!  The routines below are based on or meant to replace reactions computed
!  outside of KPP within sulfate_mod.
!
!  Rates are set for the following:
!  1) Cloud sulfur chemistry
!  2) (If Dust acid uptake) Dust acid uptake reactions
!  - MSL, 29-Mar-2021, 7-May-2021
!
!  Seasalt aerosol chemistry reaction rate-law functions are now contained
!  in module fullchem_RateLawFuncs.F90.
!
!  NOTE This routine defines the variables State_Chm%HSO3_aq and
!  State_Chm%SO3aq.  Therefore, we must call Set_Sulfur_Chem_Rates after
!  Set_KPP_GridBox_Values, but before fullchem_SetStateHet.  Otherwise we
!  will not be able to copy State_Chm%HSO3_aq to State_Het%HSO3_aq and
!  State_Chm%SO3_aq to State_Het%SO3_aq properly.
!    -- Mike Long, Bob Yantosca (08 Mar 2022)
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Initialize
    RC       = GC_SUCCESS
    size_res = .FALSE.

    !========================================================================
    ! Do this when KPP is handling aqueous sulfur chemistry
    !
    ! From Alexander et al., buffering capacity (or alkalinity) of
    ! sea-salt aerosols is equal to 0.07 equivalents per kg dry sea salt
    ! emitted Gurciullo et al., 1999. JGR 104(D17) 21,719-21,731.
    ! tdf; MSL
    !========================================================================
    IF ( .not. State_Chm%Do_SulfateMod_SeaSalt ) THEN

       ! Convert alkalinity from [molec/cm3] to equivalents
       CALL fullchem_ConvertAlkToEquiv()

       ! Compute reaction rates for aqueous sulfur chemistry
       ! (i.e. S(IV)->S(VI), HCl,  and HNO3)
       CALL fullchem_SulfurAqChem( I          = I,                           &
                                   J          = J,                           &
                                   L          = L,                           &
                                   Input_Opt  = Input_Opt,                   &
                                   State_Chm  = State_Chm,                   &
                                   State_Grid = State_Grid,                  &
                                   State_Met  = State_Met,                   &
                                   RC         = RC                          )
    ENDIF

    !========================================================================
    ! Do this when KPP is handling SO2 cloud chemistry ...
    !========================================================================
    IF ( .not. State_Chm%Do_SulfateMod_Cld ) THEN

       ! Compute reaction rates [1/s] for sulfate in cloud for KPP chem mech
       CALL fullchem_SulfurCldChem(  I          = I,                         &
                                     J          = J,                         &
                                     L          = L,                         &
                                     Input_Opt  = Input_Opt,                 &
                                     State_Chm  = State_Chm,                 &
                                     State_Diag = State_Diag,                &
                                     State_Grid = State_Grid,                &
                                     State_Met  = State_Met,                 &
                                     size_res   = size_res,                  &
                                     RC         = RC                        )

       ! Update HSO3- and SO3-- concentrations [molec/cm3]
       State_Chm%fupdateHOBr(I,J,L) = 1.0_fp
       State_Chm%fupdateHOCl(I,J,L) = 1.0_fp

    ENDIF

  END SUBROUTINE Set_Sulfur_Chem_Rates
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Kpp_GridBox_Values
!
! !DESCRIPTION: Populates KPP variables in the gckpp_Global.F90 module
!  for a particular (I,J,L) grid box.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Kpp_GridBox_Values( I,         J,         L,                &
                                     Input_Opt, State_Chm, State_Grid,       &
                                     State_Met, RC                          )
!
! !USES:
!
    USE ErrCode_Mod
    USE GcKpp_Global
    USE GcKpp_Parameters
    USE Input_Opt_Mod,          ONLY : OptInput
    USE PhysConstants,          ONLY : CONSVAP, RGASLATM, RSTARG
    USE Pressure_Mod,           ONLY : Get_Pcenter
    USE State_Chm_Mod,          ONLY : ChmState
    USE State_Grid_Mod,         ONLY : GrdState
    USE State_Met_Mod,          ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)  :: I, J, L
    TYPE(OptInput), INTENT(IN)  :: Input_Opt
    TYPE(ChmState), INTENT(IN)  :: State_Chm
    TYPE(GrdState), INTENT(IN)  :: State_Grid
    TYPE(MetState), INTENT(IN)  :: State_Met
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: F,       N,        NA,    KppId,    SpcId
    REAL(f8)           :: CONSEXP, VPRESH2O

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !========================================================================
    ! Set_Kpp_GridBox_Values begins here!
    !========================================================================

    ! Initialization
    RC      = GC_SUCCESS
    NA      = State_Chm%nAeroType
    ErrMsg  = ''
    ThisLoc = &
     ' -> at Set_Kpp_Gridbox_Values (in module GeosCore/fullchem_mod.F90)'

    !========================================================================
    ! Populate global variables in gckpp_Global.F90
    !========================================================================

    ! Solar quantities
    SUNCOS          = State_Met%SUNCOSmid(I,J)

    ! Pressure and density quantities
    NUMDEN          = State_Met%AIRNUMDEN(I,J,L)
    H2O             = State_Met%AVGW(I,J,L) * NUMDEN
    PRESS           = Get_Pcenter( I, J, L )

    ! Temperature quantities
    TEMP            = State_Met%T(I,J,L)
    INV_TEMP        = 1.0_dp   / TEMP
    TEMP_OVER_K300  = TEMP     / 300.0_dp
    K300_OVER_TEMP  = 300.0_dp / TEMP
    SR_TEMP         = SQRT( TEMP )
    FOUR_R_T        = 4.0_dp * CON_R    * TEMP
    FOUR_RGASLATM_T = 4.0_dp * RGASLATM * TEMP
    EIGHT_RSTARG_T  = 8.0_dp * RSTARG   * TEMP

    ! Relative humidity quantities
    CONSEXP         = 17.2693882_dp * (TEMP - 273.16_dp) / (TEMP - 35.86_dp)
    VPRESH2O        = CONSVAP * EXP( CONSEXP ) / TEMP
    RELHUM          = ( H2O / VPRESH2O ) * 100_dp

  END SUBROUTINE Set_Kpp_GridBox_Values
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diag_OH_HO2_O1D_O3P
!
! !DESCRIPTION: Archives the chemical production of OH, HO2, O1D, O3P.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diag_OH_HO2_O1D_O3P( Input_Opt,  State_Chm, State_Diag, &
                                  State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE Species_Mod,    ONLY : SpcConc
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  This routine replaces both OHSAVE and DIAGOH.  Those routines were needed
!  for SMVGEAR, when we had separate arrays for the non-advected species.
!  But now, all species are stored in State_Chm%Species%Conc, so the various
!  arrays (SAVEOH, SAVEHO2, etc.) are no longer necessary.  We can now just
!  just get values directly from State_Chm%Species%Conc.
!
!  Also note: for the netCDF diagnostics, we have removed multiplication by
!  LTOH etc arrays.  These are almost always set between 0 and 24.
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

    ! Scalars
    LOGICAL            :: Do_Diag
    INTEGER            :: I,       J,       L

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    ! Pointers
    REAL(fp),      POINTER  :: AirNumDen(:,:,:)
    TYPE(SpcConc), POINTER  :: Spc      (:    )

    !=======================================================================
    ! Diag_OH_HO2_O1D_O3P begins here!
    !=======================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Point to the array of species concentrations
    AirNumDen => State_Met%AirNumDen
    Spc       => State_Chm%Species

    ! Zero the netCDF diagnostic arrays (if activated) above the
    ! tropopause or mesopause to avoid having leftover values
    ! from previous timesteps
#ifdef MODEL_GEOS
    IF ( State_Diag%Archive_O3concAfterChem  ) THEN
       State_Diag%O3concAfterChem  = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_RO2concAfterChem ) THEN
       State_Diag%RO2concAfterChem = 0.0_f4
    ENDIF
#endif
    IF ( State_Diag%Archive_OHconcAfterChem  ) THEN
       State_Diag%OHconcAfterChem  = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_HO2concAfterChem ) THEN
       State_Diag%HO2concAfterChem = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_O1DconcAfterChem ) THEN
       State_Diag%O1DconcAfterChem = 0.0_f4
    ENDIF
    IF ( State_Diag%Archive_O3PconcAfterChem ) THEN
       State_Diag%O3PconcAfterChem = 0.0_f4
    ENDIF

!$OMP PARALLEL DO        &
!$OMP DEFAULT( SHARED )  &
!$OMP PRIVATE( I, J, L ) &
!$OMP SCHEDULE( DYNAMIC )
      DO L = 1, State_Grid%NZ
      DO J = 1, State_Grid%NY
      DO I = 1, State_Grid%NX

         ! Skip non-chemistry boxes
         IF ( .not. State_Met%InChemGrid(I,J,L) ) THEN
            ! Skip to next grid box
            CYCLE
         ENDIF

         !------------------------------------------------------------------
         ! OH concentration [molec/cm3]
         !------------------------------------------------------------------
         IF ( ok_OH ) THEN

            ! HISTORY (aka netCDF diagnostics)
            IF ( State_Diag%Archive_OHconcAfterChem ) THEN
               State_Diag%OHconcAfterChem(I,J,L) = Spc(id_OH)%Conc(I,J,L)
            ENDIF
#ifdef MODEL_GEOS
            IF ( State_Diag%Archive_O3concAfterChem ) THEN
               State_Diag%O3concAfterChem(I,J,L) = Spc(id_O3)%Conc(I,J,L)
            ENDIF
            IF ( Archive_RO2concAfterChem ) THEN
               IF ( id_A3O2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_A3O2)%Conc(I,J,L)
               IF ( id_ATO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_ATO2)%Conc(I,J,L)
               IF ( id_B3O2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_B3O2)%Conc(I,J,L)
               IF ( id_BRO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_BRO2)%Conc(I,J,L)
               IF ( id_ETO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_ETO2)%Conc(I,J,L)
               IF ( id_HO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_HO2)%Conc(I,J,L)
               IF ( id_LIMO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_LIMO2)%Conc(I,J,L)
               IF ( id_MO2 > 0  ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_MO2)%Conc(I,J,L)
               IF ( id_PIO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_PIO2)%Conc(I,J,L)
               IF ( id_PO2 > 0  ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_PO2)%Conc(I,J,L)
               IF ( id_PRN1 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_PRN1)%Conc(I,J,L)
               IF ( id_R4N1 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_R4N1)%Conc(I,J,L)
               IF ( id_R4O2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_R4O2)%Conc(I,J,L)
               IF ( id_TRO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_TRO2)%Conc(I,J,L)
               IF ( id_XRO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_XRO2)%Conc(I,J,L)
               IF ( id_IHOO1 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_IHOO1)%Conc(I,J,L)
               IF ( id_IHOO4 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_IHOO4)%Conc(I,J,L)
               IF ( id_ICHOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_ICHOO)%Conc(I,J,L)
               IF ( id_IHPOO1 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_IHPOO1)%Conc(I,J,L)
               IF ( id_IHPOO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_IHPOO2)%Conc(I,J,L)
               IF ( id_IHPOO3 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_IHPOO3)%Conc(I,J,L)
               IF ( id_IEPOXAOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_IEPOXAOO)%Conc(I,J,L)
               IF ( id_IEPOXBOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_IEPOXBOO)%Conc(I,J,L)
               IF ( id_C4HVP1 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_C4HVP1)%Conc(I,J,L)
               IF ( id_C4HVP2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_C4HVP2)%Conc(I,J,L)
               IF ( id_HPALD1OO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_HPALD1OO)%Conc(I,J,L)
               IF ( id_HPALD2OO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_HPALD2OO)%Conc(I,J,L)
               IF ( id_ISOPNOO1 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_ISOPNOO1)%Conc(I,J,L)
               IF ( id_ISOPNOO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_ISOPNOO2)%Conc(I,J,L)
               IF ( id_INO2D > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_INO2D)%Conc(I,J,L)
               IF ( id_INO2B > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_INO2B)%Conc(I,J,L)
               IF ( id_IDHNBOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_IDHNBOO)%Conc(I,J,L)
               IF ( id_IDHNDOO1 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_IDHNDOO1)%Conc(I,J,L)
               IF ( id_IDHNDOO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_IDHNDOO2)%Conc(I,J,L)
               IF ( id_IHPNBOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_IHPNBOO)%Conc(I,J,L)
               IF ( id_IHPNDOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_IHPNDOO)%Conc(I,J,L)
               IF ( id_ICNOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_ICNOO)%Conc(I,J,L)
               IF ( id_IDNOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(id_IDNOO)%Conc(I,J,L)
            ENDIF
#endif

         ENDIF

         !------------------------------------------------------------------
         ! HO2 concentration [v/v]
         !------------------------------------------------------------------
         IF ( ok_HO2 ) THEN
            IF ( State_Diag%Archive_HO2concAfterChem ) THEN
               State_Diag%HO2concAfterChem(I,J,L) = ( Spc(id_HO2)%Conc(I,J,L) &
                                                  /   AirNumDen(I,J,L)      )
            ENDIF

         ENDIF

         !---------------------------------------------------------------
         ! O1D concentration [molec/cm3]
         !---------------------------------------------------------------
         IF ( ok_O1D ) THEN
            IF ( State_Diag%Archive_O1DconcAfterChem ) THEN
               State_Diag%O1DconcAfterChem(I,J,L) = Spc(id_O1D)%Conc(I,J,L)
            ENDIF
         ENDIF


         !---------------------------------------------------------------
         ! O3P concentration [molec/cm3]
         !---------------------------------------------------------------
         IF ( ok_O3P ) THEN
            IF ( State_Diag%Archive_O3PconcAfterChem ) THEN
               State_Diag%O3PconcAfterChem(I,J,L) = Spc(id_O3P)%Conc(I,J,L)
            ENDIF
         ENDIF

      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointers
      AirNumDen => NULL()
      Spc       => NULL()

  END SUBROUTINE Diag_OH_HO2_O1D_O3P
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Diag_Metrics
!
! !DESCRIPTION: Computes mass-weighted mean OH columns (full-atmosphere and
!  trop-only) that are needed to compute the overall mean OH concentration.
!  This is used as a metric as to how reactive, or "hot" the chemistry
!  mechanism is.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Diag_Metrics( Input_Opt,  State_Chm, State_Diag,                &
                           State_Grid, State_Met, RC                        )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE PhysConstants,  ONLY : AVO
    USE PhysConstants,  ONLY : XNUMOLAIR
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm    ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!
! !REMARKS:
!  References:
!  (1) Prather, M. and C. Spivakovsky, "Tropospheric OH and
!       the lifetimes of hydrochlorofluorocarbons", JGR,
!       Vol 95, No. D11, 18723-18729, 1990.
!  (2) Lawrence, M.G, Joeckel, P, and von Kuhlmann, R., "What
!       does the global mean OH concentraton tell us?",
!       Atm. Chem. Phys, 1, 37-49, 2001.
!  (3) WMO/UNEP Scientific Assessment of Ozone Depletion: 2010
!
! !REVISION HISTORY:
!  18 Aug 2020 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(f8), PARAMETER :: M3toCM3        = 1.0e+6_f8
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL,  SAVE      :: first          = .TRUE.
    INTEGER,  SAVE      :: id_OH          = -1
    REAL(f8), SAVE      :: MCM3toKGM3_OH  = -1.0_f8

    ! Scalars
    INTEGER             :: I,           J,           L
    REAL(f8)            :: airMass_m,   airmass_kg,  airMassFull
    REAL(f8)            :: airMassTrop, Ktrop,       LossOHbyCH4
    REAL(f8)            :: LossOHbyMCF, OHconc_mcm3, OHmassWgt
    REAL(f8)            :: OHmassFull,  OHmassTrop,  volume

    ! Strings
    CHARACTER(LEN=255)  :: errMsg,      thisLoc

    !========================================================================
    ! Compute_Mean_OH_and_CH4 begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Compute_Mean_OH (in module GeosCore/diagnostics_mod.F90)'

    ! Exit if we have not turned on the Metrics collection
    IF ( .not. State_Diag%Archive_Metrics ) RETURN

    !========================================================================
    ! First-time setup
    !========================================================================
    IF ( first ) THEN

       ! Get the species ID for OH
       id_OH  = Ind_('OH')
       IF ( id_OH < 0 ) THEN
          errMsg = 'OH is not a defined species in this simulation!!!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Convert [molec OH cm-3] --> [kg OH m-3]
       MCM3toKGM3_OH  = M3toCM3                                              &
                      * ( State_Chm%SpcData(id_OH)%Info%MW_g * 1.0e-3_f8 )   &
                      / AVO


       ! Reset first-time flag
       first  = .FALSE.
    ENDIF

    !========================================================================
    ! Loop over surface boxes and compute mean OH in columns
    !========================================================================
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I,           J,           L,           airMass_kg        )&
    !$OMP PRIVATE( airMass_m,   airMassFull, airMassTrop, Ktrop             )&
    !$OMP PRIVATE( LossOHbyCH4, LossOHbyMCF, OHconc_mcm3, OHmassWgt         )&
    !$OMP PRIVATE( OHmassFull,  OHmassTrop,  volume                         )&
    !$OMP SCHEDULE( DYNAMIC, 4                                              )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !--------------------------------------------------------------------
       ! Zero column-specific quantities
       !--------------------------------------------------------------------
       airMass_kg  = 0.0_f8
       airMass_m   = 0.0_f8
       airMassFull = 0.0_f8
       airMassTrop = 0.0_f8
       Ktrop       = 0.0_f8
       LossOHbyCH4 = 0.0_f8
       LossOHbyMCF = 0.0_f8
       OHconc_mcm3 = 0.0_f8
       OHmassWgt   = 0.0_f8
       OHmassFull  = 0.0_f8
       OHmassTrop  = 0.0_f8
       volume      = 0.0_f8

       !--------------------------------------------------------------------
       ! Loop over the number of levels in the full-atmosphere column
       ! Limit the computations to boxes in the chemistry grid
       !--------------------------------------------------------------------
       DO L = 1, State_Met%ChemGridLev(I,J)

          ! Compute box volume [cm3], and air mass ([molec] and [kg])
          ! Note: air mass in [molec] is also the atmospheric burden of
          ! methyl chloroform (aka MCF, formula=CH3CCl3), since we assume
          ! a uniform mixing ratio (=1) of MCF in air.
          volume      = State_Met%AIRVOL(I,J,L)    * M3toCM3
          airMass_m   = State_Met%AIRNUMDEN(I,J,L) * volume
          airMass_kg  = airMass_m / XNUMOLAIR

          ! OH concentration [molec cm-3]
          OHconc_mcm3 = State_Chm%Species(id_OH)%Conc(I,J,L)

          ! Airmass-weighted OH [kg air * (kg OH  m-3)]
          OHmassWgt   = airMass_kg * ( OHconc_mcm3  * MCM3toKGM3_OH  )

          ! Sum the air mass, mass-weighted CH4,
          ! and mass-weighted OH in the full-atm column
          airMassFull = airMassFull + airMass_kg
          OHmassFull  = OHmassFull  + OHMassWgt

          !------------------------------------------------------------------
          ! Only do the following for tropospheric boxes ...
          !------------------------------------------------------------------
          IF ( State_Met%InTroposphere(I,J,L) ) THEN

             ! Sum the air mass, mass-weighted CH4,
             ! and mass-weighted OH in the trop-only column
             airMassTrop = airMassTrop + airMass_kg
             OHmassTrop  = OHmassTrop  + OHmassWgt

             ! Compute CH4 loss rate in troposphere
             ! Ktrop (Arrhenius parameter) has units [cm3/molec/s]
             ! OHconc has units [molec/cm3]
             ! AirMass has units [molec]
             ! Resultant units of CH4 loss rate = [molec/s]
             Ktrop = 2.45e-12_f8 * EXP( -1775.0_f8 / State_Met%T(I,J,L) )
             LossOHbyCH4 = LossOHbyCH4 + ( Ktrop * OHconc_mcm3 * airMass_m )

             ! Compute MCF loss rate in the troposphere
             ! Ktrop (Arrhenius parameter) has units [cm3/molec/s]
             ! OHconc has units [molec/cm3]
             ! AirMass has units [molec]
             ! Resultant units of MCF loss rate = [molec/s]
             Ktrop = 1.64e-12_f8 * EXP( -1520.0_f8 / State_Met%T(I,J,L) )
             LossOHbyMCF = LossOHbyMCF + ( Ktrop * OHconc_mcm3 * airMass_m )

          ENDIF
       ENDDO

       !---------------------------------------------------------------------
       ! HISTORY (aka netCDF diagnostics)
       ! Air mass [kg], full-atmosphere and trop-only column sums
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_AirMassColumnFull ) THEN
          State_Diag%AirMassColumnFull(I,J) = airMassFull
       ENDIF

       IF ( State_Diag%Archive_AirMassColumnTrop ) THEN
          State_Diag%AirMassColumnTrop(I,J) = airMassTrop
       ENDIF

       !---------------------------------------------------------------------
       ! HISTORY (aka netCDF diagnostics)
       ! Airmass-weighted mean OH [kg air * (kg OH m-3)]
       ! Full-atmosphere and trop-only column sums
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_OHwgtByAirMassColumnFull ) THEN
          State_Diag%OHwgtByAirMassColumnFull(I,J) = OHmassFull
       ENDIF

       IF ( State_Diag%Archive_OHwgtByAirMassColumnTrop ) THEN
          State_Diag%OHwgtByAirMassColumnTrop(I,J) = OHmassTrop
       ENDIF

       !-----------------------------------------------------------------
       ! HISTORY (aka netCDF diagnostics)
       !
       ! OH loss by CH4 + OH in troposphere [molec/s] and
       ! OH loss by MCF + OH in troposphere [molec/s]
       ! Full-atmosphere and trop-only column sums
       !----------------------------------------------------------------
       IF ( State_Diag%Archive_LossOHbyCH4columnTrop ) THEN
          State_Diag%LossOHbyCH4columnTrop(I,J) = LossOHbyCH4
       ENDIF

       IF ( State_Diag%Archive_LossOHbyMCFcolumnTrop ) THEN
          State_Diag%LossOHbyMCFcolumnTrop(I,J) = LossOHbyMCF
       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Diag_Metrics
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_fullchem
!
! !DESCRIPTION: Subroutine Init\_FullChem is used to allocate arrays for the
!  KPP solver.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_FullChem( Input_Opt, State_Chm, State_Diag, RC )
!
! !USES:
!
    USE aciduptake_DustChemFuncs, ONLY : aciduptake_InitDustChem
    USE ErrCode_Mod
    USE fullchem_SulfurChemFuncs, ONLY : fullchem_InitSulfurChem
    USE Gckpp_Monitor,            ONLY : Eqn_Names, Fam_Names
    USE Gckpp_Precision
    USE Gckpp_Parameters,         ONLY : nFam, nReact
    USE Gckpp_Global,             ONLY : Henry_K0, Henry_CR, MW, SR_MW
    USE Input_Opt_Mod,            ONLY : OptInput
    USE State_Chm_Mod,            ONLY : ChmState
    USE State_Chm_Mod,            ONLY : Ind_
    USE State_Diag_Mod,           ONLY : DgnState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Diagnostics State object
    TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  14 Dec 2015 - M.S. Long   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: KppId,    N
    LOGICAL            :: prtDebug

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg,   ThisLoc

    !=======================================================================
    ! Init_FullChem begins here!
    !=======================================================================

    ! Assume success
    RC       = GC_SUCCESS

    ! Do the following only if it is a full-chemistry simulation
    ! NOTE: If future specialty simulations use the KPP solver,
    ! modify the IF statement accordingly to allow initialization
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) RETURN

    !=======================================================================
    ! Initialize variables
    !=======================================================================
    ErrMsg   = ''
    ThisLoc  = ' -> at Init_FullChem (in module GeosCore/FullChem_mod.F90)'
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Debug output
    IF ( prtDebug ) THEN
       WRITE( 6, 100 )
100    FORMAT( '     - INIT_FULLCHEM: Allocating arrays' )

       WRITE( 6 ,'(a)' ) ' KPP Reaction Reference '
       DO N = 1, NREACT
          WRITE( 6, '(i8,a3,a85)' ) N,' | ',EQN_NAMES(N)
       END DO
    ENDIF

    ! Initialize species flags
    id_CH4      = Ind_( 'CH4', 'A'     ) ! CH4 advected species
    id_HO2      = Ind_( 'HO2'          )
    id_NH3      = Ind_( 'NH3'          )
    id_O3P      = Ind_( 'O'            )
    id_O1D      = Ind_( 'O1D'          )
    id_OH       = Ind_( 'OH'           )
    id_SALAAL   = Ind_( 'SALAAL'       )
    id_SALCAL   = Ind_( 'SALCAL'       )
    id_SO4      = Ind_( 'SO4'          )
    id_SALC     = Ind_( 'SALC'         )

#ifdef MODEL_GEOS
    ! ckeller
    id_O3       = Ind_( 'O3'           )
    id_A3O2     = Ind_( 'A3O2'         )
    id_ATO2     = Ind_( 'ATO2'         )
    id_BRO2     = Ind_( 'BRO2'         )
    id_ETO2     = Ind_( 'ETO2'         )
    id_LIMO2    = Ind_( 'LIMO2'        )
    id_MO2      = Ind_( 'MO2'          )
    id_PIO2     = Ind_( 'PIO2'         )
    id_PO2      = Ind_( 'PO2'          )
    id_PRN1     = Ind_( 'PRN1'         )
    id_R4N1     = Ind_( 'R4N1'         )
    id_R4O2     = Ind_( 'R4O2'         )
    id_TRO2     = Ind_( 'TRO2'         )
    id_XRO2     = Ind_( 'XRO2'         )
    id_IHOO1    = Ind_( 'IHOO1'        )
    id_IHOO4    = Ind_( 'IHOO4'        )
    id_IHCOO    = Ind_( 'IHCOO'        )
    id_IHPOO1   = Ind_( 'IHPOO1'       )
    id_IHPOO2   = Ind_( 'IHPOO2'       )
    id_IHPOO3   = Ind_( 'IHPOO3'       )
    id_IEPOXAOO = Ind_( 'IEPOXAOO'     )
    id_IEPOXBOO = Ind_( 'IEPOXBOO'     )
    id_C4HVP1   = Ind_( 'C4HVP1'       )
    id_C4HVP2   = Ind_( 'C4HVP2'       )
    id_HPALD1OO = Ind_( 'HPALD1OO'     )
    id_HPALD2OO = Ind_( 'HPALD2OO'     )
    id_ISOPNOO1 = Ind_( 'ISOPNOO1'     )
    id_ISOPNOO2 = Ind_( 'ISOPNOO2'     )
    id_INO2B    = Ind_( 'INO2B'        )
    id_INO2D    = Ind_( 'INO2D'        )
    id_IDHNBOO  = Ind_( 'IDHNBOO'      )
    id_IDHNDOO1 = Ind_( 'IDHNDOO1'     )
    id_IDHNDOO2 = Ind_( 'IDHNDOO2'     )
    id_IHPNBOO  = Ind_( 'IHPNBOO'      )
    id_IHPNDOO  = Ind_( 'IHPNDOO'      )
    id_ICNOO    = Ind_( 'ICNOO'        )
    id_IDNOO    = Ind_( 'IDNOO'        )
#endif

    ! Set flags to denote if each species is defined
    ok_HO2      = ( id_HO2 > 0         )
    ok_O1D      = ( id_O1D > 0         )
    ok_O3P      = ( id_O3P > 0         )
    ok_OH       = ( id_OH  > 0         )

    ! Should we archive OH, HO2, O1D, O3P diagnostics?
    Do_Diag_OH_HO2_O1D_O3P = (                                               &
#ifdef MODEL_GEOS
                               State_Diag%Archive_O3concAfterChem       .or. &
                               State_Diag%Archive_RO2concAfterChem      .or. &
#endif
                               State_Diag%Archive_OHconcAfterChem       .or. &
                               State_Diag%Archive_HO2concAfterChem      .or. &
                               State_Diag%Archive_O1DconcAfterChem      .or. &
                               State_Diag%Archive_O3PconcAfterChem          )

    !=======================================================================
    ! Save physical parameters from the species database into KPP arrays
    ! in gckpp_Global.F90.  These are for the hetchem routines.
    !=======================================================================
    DO KppId = 1, State_Chm%nKppSpc + State_Chm%nOmitted
       N                  = State_Chm%Map_KppSpc(KppId)
       IF ( N > 0 ) THEN
          MW(KppId)       = State_Chm%SpcData(N)%Info%MW_g
          SR_MW(KppId)    = SQRT( MW(KppId ) )
          HENRY_K0(KppId) = State_Chm%SpcData(N)%Info%Henry_K0
          HENRY_CR(KppId) = State_Chm%SpcData(N)%Info%Henry_CR
       ENDIF
    ENDDO

    !=======================================================================
    ! Allocate arrays
    !=======================================================================

    ! Initialize
    id_PCO =  -1
    id_LCH4 = -1

    !--------------------------------------------------------------------
    ! Pre-store the KPP indices for each KPP prod/loss species or family
    !--------------------------------------------------------------------
    IF ( nFam > 0 ) THEN

       ! Allocate mapping array for KPP Ids for ND65 bpch diagnostic
       ALLOCATE( PL_Kpp_Id( nFam ), STAT=RC )
       CALL GC_CheckVar( 'fullchem_mod.F90:PL_Kpp_Id', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN

       ! Loop over all KPP prod/loss species
       DO N = 1, nFam

          ! NOTE: KppId is the KPP ID # for each of the prod and loss
          ! diagnostic species.  This is the value used to index the
          ! KPP "VAR" array (in module gckpp_Global.F90).
          KppID = Ind_( TRIM ( Fam_Names(N) ), 'K' )

          ! Find the KPP Id corresponding to PCO and LCH4
          ! so that we can save output for tagged CO simulations
          IF ( TRIM( Fam_Names(N) ) == 'PCO'  ) id_PCO  = KppId
          IF ( TRIM( Fam_Names(N) ) == 'LCH4' ) id_LCH4 = KppId

          ! Exit if an invalid ID is encountered
          IF ( KppId <= 0 ) THEN
             ErrMsg = 'Invalid KPP ID for prod/loss species: '            // &
                       TRIM( Fam_Names(N) )
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! If the species ID is OK, save in ND65_Kpp_Id
          PL_Kpp_Id(N) = KppId
       ENDDO

    ENDIF

    !--------------------------------------------------------------------
    ! If we are archiving the P(CO) from CH4 and from NMVOC from a fullchem
    ! simulation for the tagCO simulation, throw an error if we cannot find
    ! the PCO or LCH4 prod/loss families in this KPP mechanism.
    !--------------------------------------------------------------------
    IF ( State_Diag%Archive_ProdCOfromCH4    .or.                            &
         State_Diag%Archive_ProdCOfromNMVOC ) THEN

       IF ( id_PCO < 1 ) THEN
          ErrMsg = 'Could not find PCO in list of KPP families!  This   ' // &
                   'is needed to archive the ProdCOfromCH4 diagnostic!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( id_LCH4 < 1 ) THEN
          ErrMsg = 'Could not find LCH4 in list of KPP families!  This '  // &
                   'is needed to archive the ProdCOfromNMVOC diagnostic!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    !--------------------------------------------------------------------
    ! Initialize sulfate chemistry code (cf Mike Long)
    !--------------------------------------------------------------------
    CALL fullchem_InitSulfurChem( RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "fullchem_InitSulfurCldChem"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !--------------------------------------------------------------------
    ! Initialize dust acid uptake code (Mike Long, Bob Yantosca)
    !--------------------------------------------------------------------
    IF ( Input_Opt%LDSTUP ) THEN
       CALL aciduptake_InitDustChem( RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "aciduptake_InitDustChem"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE Init_FullChem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_fullchem
!
! !DESCRIPTION: Subroutine Cleanup\_FullChem deallocates module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_FullChem( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  24 Aug 2016 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! Cleanup_FullChem begins here!
    !=================================================================

    ! Initialize
    RC = GC_SUCCESS

    IF ( ALLOCATED( PL_Kpp_Id ) ) THEN
       DEALLOCATE( PL_Kpp_Id, STAT=RC  )
       CALL GC_CheckVar( 'fullchem_mod.F90:PL_Kpp_Id', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( JvCountDay ) ) THEN
       DEALLOCATE( JvCountDay, STAT=RC  )
       CALL GC_CheckVar( 'fullchem_mod.F90:JvCountDay', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( JvSumDay ) ) THEN
       DEALLOCATE( JvSumDay, STAT=RC  )
       CALL GC_CheckVar( 'fullchem_mod.F90:JvCountDay', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( JvCountMon ) ) THEN
       DEALLOCATE( JvCountMon, STAT=RC  )
       CALL GC_CheckVar( 'fullchem_mod.F90:JvCountMon', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( JvSumMon ) ) THEN
       DEALLOCATE( JvSumMon, STAT=RC  )
       CALL GC_CheckVar( 'fullchem_mod.F90:JvCountMon', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Cleanup_FullChem
!EOC
END MODULE FullChem_Mod
