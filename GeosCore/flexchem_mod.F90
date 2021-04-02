!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: flexchem_mod.F90
!
! !DESCRIPTION: Module FlexChem\_Mod contines arrays and routines for the
!  FlexChem chemical solver.
!\\
!\\
! !INTERFACE:
!
MODULE FlexChem_Mod
!
! !USES:
!
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Do_FlexChem
  PUBLIC  :: Init_FlexChem
  PUBLIC  :: Cleanup_FlexChem
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Diag_OH_HO2_O1D_O3P
  PRIVATE :: Diag_Metrics
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
  INTEGER               :: id_PCO, id_LCH4, id_SALA
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
  LOGICAL               :: ok_OH, ok_HO2, ok_O1D, ok_O3P

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
! !ROUTINE: do_flexchem
!
! !DESCRIPTION: Subroutine Do\_FlexChem is the driver subroutine for
!  full chemistry with KPP.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_FlexChem( Input_Opt,  State_Chm, State_Diag,                 &
                          State_Grid, State_Met, RC                         )
!
! !USES:
!
    USE AEROSOL_MOD,          ONLY : SOILDUST, AEROSOL_CONC, RDAER
    USE CMN_FJX_MOD
    USE DUST_MOD,             ONLY : RDUST_ONLINE
    USE ErrCode_Mod
    USE ERROR_MOD
    USE FAST_JX_MOD,          ONLY : PHOTRATE_ADJ, FAST_JX
    USE GCKPP_HetRates,       ONLY : SET_HET
    USE GCKPP_Monitor,        ONLY : SPC_NAMES, FAM_NAMES
    USE GCKPP_Parameters
    USE GCKPP_Integrator,     ONLY : INTEGRATE, NHnew
    USE GCKPP_Function
    USE GCKPP_Model
    USE GCKPP_Global
    USE GCKPP_Rates,          ONLY : UPDATE_RCONST, RCONST
    USE GCKPP_Initialize,     ONLY : Init_KPP => Initialize
    USE GcKPP_Util,           ONLY : Get_OHreactivity
    USE Input_Opt_Mod,        ONLY : OptInput
    USE PhysConstants,        ONLY : AVO
    USE PRESSURE_MOD
    USE Species_Mod,          ONLY : Species
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Chm_Mod,        ONLY : Ind_
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE Strat_Chem_Mod,       ONLY : SChem_Tend
    USE TIME_MOD,             ONLY : GET_TS_CHEM
    USE TIME_MOD,             ONLY : Get_Day
    USE TIME_MOD,             ONLY : Get_Month
    USE TIME_MOD,             ONLY : Get_Year
    USE Timers_Mod
    USE UnitConv_Mod,         ONLY : Convert_Spc_Units
    USE UCX_MOD,              ONLY : CALC_STRAT_AER
    USE UCX_MOD,              ONLY : SO4_PHOTFRAC
    USE UCX_MOD,              ONLY : UCX_NOX
    USE UCX_MOD,              ONLY : UCX_H2SO4PHOT
#ifdef TOMAS
#ifdef BPCH_DIAG
    USE TOMAS_MOD,            ONLY : H2SO4_RATE
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
    LOGICAL                :: prtDebug,   IsLocNoon
    INTEGER                :: I,          J,         L,       N
    INTEGER                :: NA,         F,         SpcID,   KppID
    INTEGER                :: P,          MONTH,     YEAR,    Day
    INTEGER                :: WAVELENGTH, IERR,      S,       Thread
    REAL(fp)               :: SO4_FRAC,   T,         TIN
    REAL(fp)               :: TOUT

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
#ifdef MODEL_GEOS
    ! Special arrays only for running GEOS-Chem in the GEOS ESM
    REAL(f4)               :: GLOB_RCONST(State_Grid%NX, State_Grid%NY,      &
                                          State_Grid%NZ, NREACT             )
    REAL(f4)               :: GLOB_JVAL(State_Grid%NX,   State_Grid%NY,      &
                                        State_Grid%NZ,   JVN_               )
#endif

    ! For tagged CO saving
    REAL(fp)               :: LCH4, PCO_TOT, PCO_CH4, PCO_NMVOC

    ! Objects
    TYPE(Species), POINTER :: SpcInfo

    ! For testing purposes
    LOGICAL                :: DO_HETCHEM
    LOGICAL                :: DO_PHOTCHEM

    ! OH reactivity and KPP reaction rate diagnostics
    REAL(fp)               :: OHreact
    REAL(dp)               :: Vloc(NVAR), Aout(NREACT)

    !=======================================================================
    ! Do_FlexChem begins here!
    !=======================================================================

    ! Initialize
    RC        =  GC_SUCCESS
    ErrMsg    =  ''
    ThisLoc   =  ' -> at Do_FlexChem (in module GeosCore/flexchem_mod.F90)'
    SpcInfo   => NULL()
    prtDebug  =  ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    Day       =  Get_Day()    ! Current day
    Month     =  Get_Month()  ! Current month
    Year      =  Get_Year()   ! Current year
    Thread    =  1

    ! Turn heterogeneous chemistry and photolysis on/off for testing
    DO_HETCHEM  = .TRUE.
    DO_PHOTCHEM = .TRUE.
    IF ( FIRSTCHEM .AND. Input_Opt%amIRoot ) THEN
       IF ( .not. DO_HETCHEM ) THEN
          WRITE( 6, '(a)' ) REPEAT( '#', 32 )
          WRITE( 6, '(a)' )  ' # Do_FlexChem: Heterogeneous chemistry' // &
                             ' is turned off for testing purposes.'
          WRITE( 6, '(a)' ) REPEAT( '#', 32 )
       ENDIF
       IF ( .not. DO_PHOTCHEM ) THEN
          WRITE( 6, '(a)' ) REPEAT( '#', 32 )
          WRITE( 6, '(a)' )  ' # Do_FlexChem: Photolysis chemistry' // &
                             ' is turned off for testing purposes.'
          WRITE( 6, '(a)' ) REPEAT( '#', 32 )
       ENDIF
    ENDIF

    ! Zero diagnostic archival arrays to make sure that we don't have any
    ! leftover values from the last timestep near the top of the chemgrid
    IF (State_Diag%Archive_Loss           ) State_Diag%Loss           = 0.0_f4
    IF (State_Diag%Archive_Prod           ) State_Diag%Prod           = 0.0_f4
    IF (State_Diag%Archive_Jval           ) State_Diag%Jval           = 0.0_f4
    IF (State_Diag%Archive_JNoon          ) State_Diag%JNoon          = 0.0_f4
    IF (State_Diag%Archive_ProdCOfromCH4  ) State_Diag%ProdCOfromCH4  = 0.0_f4
    IF (State_Diag%Archive_ProdCOfromNMVOC) State_Diag%ProdCOfromNMVOC= 0.0_f4
    IF (State_Diag%Archive_OHreactivity   ) State_Diag%OHreactivity   = 0.0_f4
    IF (State_Diag%Archive_RxnRate        ) State_Diag%RxnRate        = 0.0_f4
    IF (State_Diag%Archive_KppDiags) THEN
       IF (State_Diag%Archive_KppIntCounts) State_Diag%KppIntCounts   = 0.0_f4
       IF (State_Diag%Archive_KppJacCounts) State_Diag%KppJacCounts   = 0.0_f4
       IF (State_Diag%Archive_KppTotSteps ) State_Diag%KppTotSteps    = 0.0_f4
       IF (State_Diag%Archive_KppAccSteps ) State_Diag%KppAccSteps    = 0.0_f4
       IF (State_Diag%Archive_KppRejSteps ) State_Diag%KppRejSteps    = 0.0_f4
       IF (State_Diag%Archive_KppLuDecomps) State_Diag%KppLuDecomps   = 0.0_f4
       IF (State_Diag%Archive_KppSubsts   ) State_Diag%KppSubsts      = 0.0_f4
       IF (State_Diag%Archive_KppSmDecomps) State_Diag%KppSmDecomps   = 0.0_f4
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

#ifdef MODEL_GEOS
    GLOB_RCONST = 0.0_f4
    GLOB_JVAL   = 0.0_f4
#endif

    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End  ( "=> FlexChem",     RC ) ! started in Do_Chemistry
       CALL Timer_Start( "=> Aerosol chem", RC )
    ENDIF

    !=======================================================================
    ! Get concentrations of aerosols in [kg/m3]
    ! for FAST-JX and optical depth diagnostics
    !=======================================================================
    IF ( Input_Opt%LSULF .or. Input_Opt%LCARB .or. &
         Input_Opt%LDUST .or. Input_Opt%LSSALT ) THEN

       ! Special handling for UCX
       IF ( Input_Opt%LUCX ) THEN

          ! Calculate stratospheric aerosol properties (SDE 04/18/13)
          CALL CALC_STRAT_AER( Input_Opt, State_Chm, State_Grid, &
                               State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Calc_Strat_Aer"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Debug output
          IF ( prtDebug ) THEN
             CALL DEBUG_MSG( '### Do_FlexChem: after CALC_PSC' )
          ENDIF

       ENDIF

       ! Compute aerosol concentrations
       CALL AEROSOL_CONC( Input_Opt,  State_Chm, State_Diag, &
                          State_Grid, State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "AEROSOL_CONC"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    !=======================================================================
    ! Call RDAER -- computes aerosol optical depths for FAST-JX
    !=======================================================================
    WAVELENGTH = 0
    CALL RDAER( Input_Opt, State_Chm, State_Diag, State_Grid, State_Met, RC,  &
                MONTH,     YEAR,       WAVELENGTH )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "RDAER"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !### Debug
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### Do_FlexChem: after RDAER' )
    ENDIF

    !=======================================================================
    ! If LDUST is turned on, then we have online dust aerosol in
    ! GEOS-CHEM...so just pass SOILDUST to RDUST_ONLINE in order to
    ! compute aerosol optical depth for FAST-JX, etc.
    !
    ! If LDUST is turned off, then we do not have online dust aerosol
    ! in GEOS-CHEM...so read monthly-mean dust files from disk.
    ! (rjp, tdf, bmy, 4/1/04)
    !=======================================================================
    IF ( Input_Opt%LDUST ) THEN
       CALL RDUST_ONLINE( Input_Opt, State_Chm,  State_Diag, &
                          State_Grid, State_Met, SOILDUST,   WAVELENGTH, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "RDUST_ONLINE"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    !### Debug
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### Do_FlexChem: after RDUST' )
    ENDIF

    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End  ( "=> Aerosol chem", RC )
       CALL Timer_Start( "=> FlexChem",   RC )
    ENDIF

    !=======================================================================
    ! Zero out certain species:
    !    - isoprene oxidation counter species (dkh, bmy, 6/1/06)
    !    - isoprene-NO3 oxidation counter species (hotp, 6/1/10)
    !    - if SOA or SOA_SVPOA, aromatic oxidation counter species
    !      (dkh, 10/06/06)
    !    - if SOA_SVPOA, LNRO2H and LNRO2N for NAP (hotp 6/25/09
    !=======================================================================
    DO N = 1, State_Chm%nSpecies

       ! Get info about this species from the species database
       SpcInfo => State_Chm%SpcData(N)%Info

       ! isoprene oxidation counter species
       IF ( TRIM( SpcInfo%Name ) == 'LISOPOH' .or. &
            TRIM( SpcInfo%Name ) == 'LISOPNO3' ) THEN
          State_Chm%Species(:,:,:,N) = 0e+0_fp
       ENDIF

       ! aromatic oxidation counter species
       IF ( Input_Opt%LSOA .or. Input_Opt%LSVPOA ) THEN
          SELECT CASE ( TRIM( SpcInfo%Name ) )
             CASE ( 'LBRO2H', 'LBRO2N', 'LTRO2H', 'LTRO2N', &
                    'LXRO2H', 'LXRO2N', 'LNRO2H', 'LNRO2N' )
                State_Chm%Species(:,:,:,N) = 0e+0_fp
          END SELECT
       ENDIF

       ! Temporary fix for CO2
       ! CO2 is a dead species and needs to be set to zero to
       ! match the old SMVGEAR code (mps, 6/14/16)
       IF ( TRIM( SpcInfo%Name ) == 'CO2' ) THEN
          State_Chm%Species(:,:,:,N) = 0e+0_fp
       ENDIF

       ! Free pointer
       SpcInfo => NULL()

    ENDDO

    !=======================================================================
    ! Archive initial species mass for stratospheric tendency
    !=======================================================================
    IF ( Input_Opt%LUCX  ) THEN

       ! Loop over advected species
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( N      )
       DO NA = 1, State_Chm%nAdvect

          ! Get the species ID from the advected species ID
          N = State_Chm%Map_Advect(NA)

          ! Archive initial mass [kg]
          Before(:,:,:,N) = State_Chm%Species(:,:,:,N)

       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !======================================================================
    ! Convert species to [molec/cm3] (ewl, 8/16/16)
    !======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'molec/cm3', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'flexchem_mod.F90')
       RETURN
    ENDIF

    !=======================================================================
    ! Call photolysis routine to compute J-Values
    !=======================================================================
    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End  ( "=> FlexChem",     RC )
       CALL Timer_Start( "=> FAST-JX photolysis", RC )
    ENDIF

    ! Do Photolysis
    CALL FAST_JX( WAVELENGTH, Input_Opt,  State_Chm, &
                  State_Diag, State_Grid, State_Met, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "FAST_JX"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !### Debug
    IF ( prtDebug ) THEN
       CALL DEBUG_MSG( '### Do_FlexChem: after FAST_JX' )
    ENDIF

    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End  ( "=> FAST-JX photolysis", RC )
       CALL Timer_Start( "=> FlexChem",     RC ) ! ended in Do_Chemistry
    ENDIF

#if defined( MODEL_GEOS ) || defined( MODEL_WRF )
    ! Init diagnostics
    IF ( ASSOCIATED(State_Diag%KppError) ) THEN
       State_Diag%KppError(:,:,:) = 0.0
    ENDIF
#endif

    !=======================================================================
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
    !=======================================================================

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
    IF ( Input_Opt%LUCX  ) THEN
       ! UCX-based mechanisms
       !RTOL      = 2e-2_dp
       !RTOL      = 1e-2_dp
       RTOL      = 0.5e-2_dp
    ELSE
       ! Non-UCX mechanisms
       RTOL      = 1e-2_dp
    ENDIF

    !%%%%% SOLVER OPTIONS %%%%%

    ! Zero all slots of ICNTRL
    ICNTRL    = 0

    ! 0 - non-autonomous, 1 - autonomous
    ICNTRL(1) = 1

    ! 0 - vector tolerances, 1 - scalars
    ICNTRL(2) = 0

    ! Select Integrator
    ! ICNTRL(3)  -> selection of a particular method.
    ! For Rosenbrock, options are:
    ! = 0 :  default method is Rodas3
    ! = 1 :  method is  Ros2
    ! = 2 :  method is  Ros3
    ! = 3 :  method is  Ros4
    ! = 4 :  method is  Rodas3
    ! = 5:   method is  Rodas4
    ICNTRL(3) = 4

    ! 0 - adjoint, 1 - no adjoint
    ICNTRL(7) = 1

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

    !-----------------------------------------------------------------------
    ! MAIN LOOP: Compute reaction rates and call chemical solver
    !
    ! Variables not listed here are held THREADPRIVATE in gckpp_Global.F90
    ! !$OMP COLLAPSE(3) vectorizes the loop and !$OMP DYNAMIC(24) sends
    ! 24 boxes at a time to each core... then when that core is finished,
    ! it gets a nother chunk of 24 boxes.  This should lead to better
    ! load balancing, and will spread the sunrise/sunset boxes across
    ! more cores.
    !-----------------------------------------------------------------------
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I,        J,        L,       N                           )&
    !$OMP PRIVATE( SO4_FRAC, IERR,     RCNTRL,  ISTATUS,   RSTATE           )&
    !$OMP PRIVATE( SpcID,    KppID,    F,       P,         Vloc             )&
    !$OMP PRIVATE( Aout,     Thread,   RC,      S,         LCH4             )&
    !$OMP PRIVATE( OHreact,  PCO_TOT,  PCO_CH4, PCO_NMVOC                   )&
    !$OMP COLLAPSE( 3                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 24                                             )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !=====================================================================
       ! Initialize private loop  variables for each (I,J,L)
       ! Other private variables will be assigned in Set_Kpp_GridBox_Values
       !=====================================================================
       HET       = 0.0_dp                   ! Het chem array
       IERR      = 0                        ! KPP success or failure flag
       ISTATUS   = 0.0_dp                   ! Rosenbrock output
       PHOTOL    = 0.0_dp                   ! Photolysis array for KPP
       RCNTRL    = 0.0_fp                   ! Rosenbrock input
       RSTATE    = 0.0_dp                   ! Rosenbrock output
       SO4_FRAC  = 0.0_fp                   ! Frac of SO4 avail for photolysis
       P         = 0                        ! GEOS-Chem photolyis species ID
       LCH4      = 0.0_fp                   ! P/L diag: Methane loss rate
       PCO_TOT   = 0.0_fp                   ! P/L diag: Total P(CO)
       PCO_CH4   = 0.0_fp                   ! P/L diag: P(CO) from CH4
       PCO_NMVOC = 0.0_fp                   ! P/L diag: P(CO) from NMVOC
#ifdef MODEL_CLASSIC
#ifndef NO_OMP
       Thread    = OMP_GET_THREAD_NUM() + 1 ! OpenMP thread number
#endif
#endif

       !====================================================================
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
       !====================================================================
       IF ( State_Met%SUNCOSmid(I,J) > -0.1391731e+0_fp ) THEN

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "  -> Photolysis rates", RC, &
                               InLoop=.TRUE., ThreadNum=Thread )
          ENDIF

          ! Get the fraction of H2SO4 that is available for photolysis
          ! (this is only valid for UCX-enabled mechanisms)
          IF ( Input_Opt%LUCX ) THEN
             SO4_FRAC = SO4_PHOTFRAC( I, J, L )
          ELSE
             SO4_FRAC = 0.0_fp
          ENDIF

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

#ifdef MODEL_GEOS
             ! Archive in local array
             GLOB_JVAL(I,J,L,N) = PHOTOL(N)
#endif

             !--------------------------------------------------------------
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
             !
             !    NOTE: Depending on the simulation, some GEOS-Chem
             !    species might not map to a of the FAST-JX species
             !    (e.g. SOA species will not be present in a tropchem run).
             !
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
             !
             ! To match the legacy bpch diagnostic, we archive the sum of
             ! photolysis rates for a given GEOS-Chem species over all of
             ! the reaction branches.
             !--------------------------------------------------------------

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

             ELSE IF ( P == State_Chm%nPhotol+2 ) THEN

                ! J(O3_O3P).  This used to be stored as the nPhotol+2nd
                ! diagnostic in Jval, but needed to be broken off
                ! to facilitate cleaner diagnostic indexing (bmy, 6/3/20)
                IF ( State_Diag%Archive_JvalO3O3P ) THEN
                   State_Diag%JvalO3O3P(I,J,L) =                             &
                   State_Diag%JvalO3O3P(I,J,L) + PHOTOL(N)
                ENDIF

             ENDIF
          ENDDO

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "  -> Photolysis rates", RC, &
                             InLoop=.TRUE., ThreadNum=Thread )
          ENDIF
       ENDIF

       !====================================================================
       ! Test if we need to do the chemistry for box (I,J,L),
       ! otherwise move onto the next box.
       !====================================================================

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
       ! Intialize KPP solver arrays: CFACTOR, VAR, FIX, etc.
       !=====================================================================
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start( "  -> Init KPP", RC,                             &
                            InLoop=.TRUE.,   ThreadNum=Thread               )
       ENDIF

       ! Initialize KPP for this grid box
       CALL Init_KPP()

       ! Copy values at each gridbox into variables in gckpp_Global.F90
       CALL Set_Kpp_GridBox_Values( Input_Opt, State_Chm, State_Grid,        &
                                    State_Met, I,         J,                 &
                                    L,         RC                           )

       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End  ( "  -> Init KPP", RC,                             &
                            InLoop=.TRUE.,   ThreadNum=Thread               )
       ENDIF

       !====================================================================
       ! Get rates for heterogeneous chemistry
       !====================================================================
       IF ( DO_HETCHEM ) THEN
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "  -> Het chem rates", RC, &
                               InLoop=.TRUE., ThreadNum=Thread )
          ENDIF

          ! Compute heterogeneous rates
          CALL SET_HET( I, J, L, Input_Opt, State_Chm, State_Met )

          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "  -> Het chem rates", RC, &
                             InLoop=.TRUE., ThreadNum=Thread )
          ENDIF
       ENDIF

       !====================================================================
       ! Start KPP main timer and prepare arrays
       !====================================================================
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start( "  -> KPP", RC, InLoop=.TRUE., ThreadNum=Thread )
       ENDIF

       ! Zero out dummy species index in KPP after call to SET_HET
       DO F = 1, NFAM
          KppID = PL_Kpp_Id(F)
          IF ( KppID > 0 ) C(KppID) = 0.0_dp
       ENDDO

       ! Set VAR and FIX arrays
       ! This has to be done after the zeroing above
       VAR(1:NVAR) = C(1:NVAR)
       FIX         = C(NVAR+1:NSPEC)

       !====================================================================
       ! Update reaction rates
       !====================================================================
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start( "     RCONST", RC, InLoop=.TRUE., ThreadNum=Thread )
       ENDIF

       ! Update the array of rate constants
       CALL Update_RCONST( )

       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( "     RCONST", RC, InLoop=.TRUE., ThreadNum=Thread )
       ENDIF

       !--------------------------------------------------------------------
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Archive KPP reaction rates [s-1]
       ! See gckpp_Monitor.F90 for a list of chemical reactions
       !--------------------------------------------------------------------
       IF ( State_Diag%Archive_RxnRate ) THEN
          CALL Fun( VAR, FIX, RCONST, Vloc, Aout=Aout )
          DO S = 1, State_Diag%Map_RxnRate%nSlots
             N = State_Diag%Map_RxnRate%slot2Id(S)
             State_Diag%RxnRate(I,J,L,S) = Aout(N)
          ENDDO
       ENDIF

       !=====================================================================
       ! Set options for the KPP Integrator (M. J. Evans)
       !
       ! NOTE: Because RCNTRL(3) is set to an array value that
       ! depends on (I,J,L), we must declare RCNTRL as PRIVATE
       ! within the OpenMP parallel loop and define it just
       ! before the call to to Integrate. (bmy, 3/24/16)
       !=====================================================================

       ! Zero all slots of RCNTRL
       RCNTRL    = 0.0_fp

       ! Starting value for integration time step
       RCNTRL(3) = State_Chm%KPPHvalue(I,J,L)

       !=====================================================================
       ! Integrate the box forwards
       !=====================================================================
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start( "     Integrate 1", RC,                          &
                            InLoop=.TRUE.,      ThreadNum=Thread            )
       ENDIF

       ! Call the KPP integrator
       CALL Integrate( TIN,    TOUT,    ICNTRL,      &
                       RCNTRL, ISTATUS, RSTATE, IERR )

       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( "     Integrate 1", RC, InLoop=.TRUE., ThreadNum=Thread )
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

       !------------------------------------------------------------------
       ! HISTORY: Archive KPP solver diagnostics
       !------------------------------------------------------------------
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
          IF ( State_Diag%Archive_KppTotSteps ) THEN
             State_Diag%KppAccSteps(I,J,L) = ISTATUS(4)
          ENDIF

          ! # of rejected internal timesteps
          IF ( State_Diag%Archive_KppTotSteps ) THEN
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
       ENDIF

       !------------------------------------------------------------------
       ! Try another time if it failed
       !------------------------------------------------------------------
       IF ( IERR < 0 ) THEN

          ! Reset first time step and start concentrations
          ! Retry the integration with non-optimized
          ! settings
          RCNTRL(3)  = 0e+0_fp
          CALL Init_KPP( )
          VAR = C(1:NVAR)
          FIX = C(NVAR+1:NSPEC)
          CALL Update_RCONST( )
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "     Integrate 2", RC, InLoop=.TRUE., ThreadNum=Thread )
          ENDIF
          CALL Integrate( TIN,    TOUT,    ICNTRL,                           &
                          RCNTRL, ISTATUS, RSTATE, IERR                     )
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "     Integrate 2", RC, InLoop=.TRUE., ThreadNum=Thread )
          ENDIF
          !------------------------------------------------------------------
          ! HISTORY: Archive KPP solver diagnostics
          ! This time, add to the existing value
          !------------------------------------------------------------------
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
             IF ( State_Diag%Archive_KppTotSteps ) THEN
                State_Diag%KppAccSteps(I,J,L) =                              &
                State_Diag%KppAccSteps(I,J,L) + ISTATUS(4)
             ENDIF

             ! # of rejected internal timesteps
             IF ( State_Diag%Archive_KppTotSteps ) THEN
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
             IF ( State_Diag%Archive_KppSmDecomps ) THEN
                State_Diag%KppSmDecomps(I,J,L) =                             &
                State_Diag%KppSmDecomps(I,J,L) + ISTATUS(8)
             ENDIF
          ENDIF

          !------------------------------------------------------------------
          ! Exit upon the second failure
          !------------------------------------------------------------------
          IF ( IERR < 0 ) THEN
             WRITE(6,*) '## INTEGRATE FAILED TWICE !!! '
             WRITE(ERRMSG,'(a,i3)') 'Integrator error code :',IERR
#if defined( MODEL_GEOS ) || defined( MODEL_WRF )
             IF ( Input_Opt%KppStop ) THEN
                CALL ERROR_STOP(ERRMSG, 'INTEGRATE_KPP')
             ! Revert to start values
             ELSE
                VAR = C(1:NVAR)
                FIX = C(NVAR+1:NSPEC)
             ENDIF
             IF ( ASSOCIATED(State_Diag%KppError) ) THEN
                State_Diag%KppError(I,J,L) = State_Diag%KppError(I,J,L) + 1.0
             ENDIF
#else
              CALL ERROR_STOP(ERRMSG, 'INTEGRATE_KPP')
#endif
          ENDIF

       ENDIF

       !--------------------------------------------------------------------
       ! Continue upon successful return...
       !--------------------------------------------------------------------

       ! Copy VAR and FIX back into C (mps, 2/24/16)
       C(1:NVAR)       = VAR(:)
       C(NVAR+1:NSPEC) = FIX(:)

       ! Save for next integration time step
       State_Chm%KPPHvalue(I,J,L) = RSTATE(Nhnew)

#ifdef MODEL_GEOS
       ! Save rate constants in global array (not used)
       GLOB_RCONST(I,J,L,:) = RCONST(:)
#endif

       !====================================================================
       ! Check we have no negative values and copy the concentrations
       ! calculated from the C array back into State_Chm%Species
       !====================================================================

       ! Loop over KPP species
       DO N = 1, NSPEC

          ! GEOS-Chem species ID
          SpcID = State_Chm%Map_KppSpc(N)

          ! Skip if this is not a GEOS-Chem species
          IF ( SpcID .eq. 0 ) CYCLE

          ! Set negative concentrations to zero
          C(N) = MAX( C(N), 0.0E0_dp )

          ! Copy concentrations back into State_Chm%Species
          State_Chm%Species(I,J,L,SpcID) = REAL( C(N), kind=fp )

       ENDDO

       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End( "  -> KPP", RC, InLoop=.TRUE., ThreadNum=Thread )
          CALL Timer_Start( "  -> Prod/loss diags", RC, &
                            InLoop=.TRUE., ThreadNum=Thread )
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
             H2SO4_RATE(I,J,L) = VAR(KppID) / AVO * 98.e-3_fp * &
                                 State_Met%AIRVOL(I,J,L)  * &
                                 1e+6_fp / DT

            IF ( H2SO4_RATE(I,J,L) < 0.0d0) THEN
              write(*,*) "H2SO4_RATE negative in flexchem_mod.F90!!", &
                 I, J, L, "was:", H2SO4_RATE(I,J,L), "  setting to 0.0d0"
              H2SO4_RATE(I,J,L) = 0.0d0
            ENDIF
          ENDIF
       ENDDO

#endif
#endif

       !====================================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Prod and loss of families or species [molec/cm3/s]
       !
       ! NOTE: KppId is the KPP ID # for each of the prod and loss
       ! diagnostic species.  This is the value used to index the
       ! KPP "VAR" array (in module gckpp_Global.F90).
       !====================================================================

       ! Chemical loss of species or families [molec/cm3/s]
       IF ( State_Diag%Archive_Loss ) THEN
          DO S = 1, State_Diag%Map_Loss%nSlots
             KppId = State_Diag%Map_Loss%slot2Id(S)
             State_Diag%Loss(I,J,L,S) = VAR(KppID) / DT
          ENDDO
       ENDIF

       ! Chemical production of species or families [molec/cm3/s]
       IF ( State_Diag%Archive_Prod ) THEN
          DO S = 1, State_Diag%Map_Prod%nSlots
             KppID = State_Diag%Map_Prod%slot2Id(S)
             State_Diag%Prod(I,J,L,S) = VAR(KppID) / DT
          ENDDO
       ENDIF

       !--------------------------------------------------------------------
       ! Archive prod/loss fields for the TagCO simulation [molec/cm3/s]
       ! (In practice, we only need to do this from benchmark simulations)
       !--------------------------------------------------------------------
       IF ( State_Diag%Archive_ProdCOfromCH4     .or.                        &
            State_Diag%Archive_ProdCOfromNMVOC ) THEN

          ! Total production of CO
          PCO_TOT   = VAR(id_PCO) / DT

          ! Loss of CO from CH4
          LCH4      = VAR(id_LCH4) / DT

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

       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End  ( "  -> Prod/loss diags", RC, &
                            InLoop=.TRUE., ThreadNum=Thread )
       ENDIF

       !====================================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Write out OH reactivity.  The OH reactivity is defined here as the
       ! inverse of its life-time. In a crude ad-hoc approach, manually add
       ! all OH reactants (ckeller, 9/20/2017)
       !====================================================================
       IF ( State_Diag%Archive_OHreactivity ) THEN
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_Start( "  -> OH reactivity diag", RC, &
                               InLoop=.TRUE., ThreadNum=Thread )
          ENDIF
          CALL Get_OHreactivity ( C, RCONST, OHreact )
          State_Diag%OHreactivity(I,J,L) = OHreact
          IF ( Input_Opt%useTimers ) THEN
             CALL Timer_End( "  -> OH reactivity diag", RC, &
                             InLoop=.TRUE., ThreadNum=Thread )
          ENDIF
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    IF ( Input_Opt%useTimers ) THEN
       CALL Timer_End( "  -> FlexChem loop", RC )
    ENDIF

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
    ! Archive OH, HO2, O1D, O3P concentrations after FlexChem solver
    !=======================================================================
    IF ( Do_Diag_OH_HO2_O1D_O3P ) THEN

       ! Save OH, HO2, O1D, O3P for the ND43 diagnostic
       ! NOTE: These might not be needed for netCDF, as they will already
       ! have been archived in State_Chm%Species output.
       CALL Diag_OH_HO2_O1D_O3P( Input_Opt,  State_Chm, State_Diag, &
                                 State_Grid, State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Diag_OH_HO2_O1D_O3P"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### Do_FlexChem: after OHSAVE' )
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
       CALL DEBUG_MSG( '### Do_FlexChem: after Diag_Metrics' )
    ENDIF

    !=======================================================================
    ! Convert species back to original units (ewl, 8/16/16)
    !=======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm,  State_Grid, State_Met, &
                            OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error!'
       CALL GC_Error( ErrMsg, RC, 'flexchem_mod.F90' )
       RETURN
    ENDIF

    !=======================================================================
    ! Special handling for UCX chemistry
    !=======================================================================
    IF ( Input_Opt%LUCX ) THEN

       !--------------------------------------------------------------------
       ! If using stratospheric chemistry, applying high-altitude
       ! active nitrogen partitioning and H2SO4 photolysis
       ! approximations  outside the chemgrid
       !--------------------------------------------------------------------
       CALL UCX_NOX( Input_Opt, State_Chm, State_Grid, State_Met )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMDR: after UCX_NOX' )
       ENDIF

       CALL UCX_H2SO4PHOT( Input_Opt, State_Chm, State_Grid, State_Met )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMDR: after UCX_H2SO4PHOT' )
       ENDIF

       !--------------------------------------------------------------------
       ! Compute stratospheric chemical tendency for UCX simulations
       !--------------------------------------------------------------------
       IF ( Input_Opt%LSCHEM ) THEN

          ! Loop over advected species
          !$OMP PARALLEL DO               &
          !$OMP DEFAULT( SHARED         ) &
          !$OMP PRIVATE( I, J, L, N, NA )
          DO NA = 1, State_Chm%nAdvect

             ! Get the species ID from the advected species ID
             N = State_Chm%Map_Advect(NA)

             ! Loop over grid boxes
             DO L = 1, State_Grid%NZ
             DO J = 1, State_Grid%NY
             DO I = 1, State_Grid%NX

                ! Aggregate stratospheric chemical tendency [kg box-1]
                ! for tropchem simulations
                IF ( State_Met%InStratosphere(I,J,L) ) THEN
                   SChem_Tend(I,J,L,N) = SChem_Tend(I,J,L,N) + &
                        ( State_Chm%Species(I,J,L,N) - Before(I,J,L,N) )
                ENDIF

             ENDDO
             ENDDO
             ENDDO
          ENDDO
          !$OMP END PARALLEL DO
       ENDIF
    ENDIF

#ifdef MODEL_GEOS
    ! Archive all needed reaction rates in state_diag
    IF ( Input_Opt%NN_RxnRconst > 0 ) THEN
       DO N = 1, Input_Opt%NN_RxnRconst
          State_Diag%RxnRconst(:,:,:,N) = GLOB_RCONST(:,:,:,Input_Opt%RxnRconst_IDs(N))
       ENDDO
    ENDIF
    IF ( Input_Opt%NN_Jvals > 0 ) THEN
       DO N = 1, Input_Opt%NN_Jvals
          State_Diag%JvalIndiv(:,:,:,N) = GLOB_JVAL(:,:,:,Input_Opt%Jval_IDs(N))
       ENDDO
    ENDIF
#endif

    ! Set FIRSTCHEM = .FALSE. -- we have gone thru one chem step
    FIRSTCHEM = .FALSE.

  END SUBROUTINE Do_FlexChem
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
  SUBROUTINE Set_Kpp_GridBox_Values( Input_Opt, State_Chm, State_Grid,       &
                                     State_Met, I,         J,                &
                                     L,         RC                          )
!
! !USES:
!
    USE ErrCode_Mod
    USE GcKpp_Global
    USE GcKpp_Parameters
    USE Input_Opt_Mod,   ONLY : OptInput
    USE PhysConstants,   ONLY : CONSVAP
    USE Pressure_Mod,    ONLY : Get_Pcenter
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Met_Mod,   ONLY : MetState

!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt
    TYPE(ChmState), INTENT(IN)  :: State_Chm
    TYPE(GrdState), INTENT(IN)  :: State_Grid
    TYPE(MetState), INTENT(IN)  :: State_Met
    INTEGER,        INTENT(IN)  :: I, J, L
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC
!
! !RETURN VALUE:
!
!
! !REMARKS:
!
! !REVISION HISTORY:
!  06 Jan 2015 - R. Yantosca - Initial version
!  See the subsequent Git history with the gitk browser!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER  :: F,       N,       NA,    KppId,    SpcId
    REAL(f8) :: CONSEXP, VPRESH2O

    !========================================================================
    ! Set_Kpp_GridBox_Values begins here!
    !========================================================================

    ! Initialization
    RC          = GC_SUCCESS                        ! Success or failure?
    nAeroType   = State_Chm%nAeroType               ! # of aerosol types
    NA          = nAeroType                         ! Local copy

    !========================================================================
    ! Copy various quantities at each grid box into gckpp_Global variables
    !========================================================================
    AClAREA     = State_Chm%AClArea(I,J,L)          ! SSA+SNA area [cm2/cm3]
    AClRADI     = State_Chm%AClRadi(I,J,L)          ! SSA+SNA radius [cm]
    AClVOL      = AClArea * AClRadi / 3e+0_fp       ! SSA+SNA volume [cm3/cm3]
    AWATER(:)   = State_Chm%IsorropAeroH2O(I,J,L,:) ! Aerosol water
    GAMMA_HO2   = Input_Opt%GAMMA_HO2
    H_PLUS      = State_Chm%IsorropHplus(I,J,L,1)   ! Proton activity [1]
    MHSO4       = State_Chm%IsorropBisulfate(I,J,L) ! Bisulvate conc [M]
    MNO3        = State_Chm%IsorropNitrate(I,J,L,1) ! Nitrate conc [M]
    MSO4        = State_Chm%IsorropSulfate(I,J,L)   ! Sulfate conc [M]
    NUMDEN      = State_Met%AIRNUMDEN(I,J,L)        ! Air density [molec/cm3]
    H2O         = State_Met%AVGW(I,J,L) * NUMDEN    ! H2O conc [molec/cm3]
    OMOC_POA    = State_Chm%OMOC_POA(I,J)           ! OM:OC ratio, POA [1]
    OMOC_OPOA   = State_Chm%OMOC_OPOA(I,J)          ! OM:OC ratio, OPOA [1]
    PRESS       = Get_Pcenter( I, J, L )            ! Pressure [hPa]
    QICE        = State_Met%QI(I,J,L)               ! Ice MR [kg/kg dry air]
    QLIQ        = State_Met%QL(I,J,L)               ! Water MR [kg/kg dry air]
    SPC_SALA    = State_Chm%Species(I,J,L,id_SALA)  ! Seasalt conc [molec/cm3]
    SUNCOS      = State_Met%SUNCOSmid(I,J)          ! cos(solar zenith angle)
    TEMP        = State_Met%T(I,J,L)                ! Temperature [K]
    VAIR        = State_Met%AIRVOL(I,J,L)*1.0e6_f8  ! Volume of air (cm3)
    XTEMP       = SQRT( TEMP )                      ! Temp**0.5 [K]
    XAREA(1:NA) = State_Chm%AeroArea(I,J,L,1:NA)    ! Aer area [cm2/cm3]
    XDENA       = NUMDEN                            ! Air dens [molec/cm3]
    XRADI(1:NA) = State_Chm%AeroRadi(I,J,L,1:NA)    ! Aer eff radius [cm]
    XVOL(1:NA)  = XAREA(1:NA) * XRADI(1:NA) / 3e+0_fp     ! Aer vol [cm3/cm3]
    XH2O(1:NA)  = State_Chm%AeroH2O(I,J,L,1:NA) * 1e-6_fp ! Aer water [cm3/cm3]

    !========================================================================
    ! Copy species concentrations into gckpp_Global variables
    !========================================================================
    DO N = 1, NSPEC
       SpcID = State_Chm%Map_KppSpc(N)
       IF ( SpcId > 0 ) THEN
          C(N) = State_Chm%Species(I,J,L,SpcID)
       ELSE
          C(N) = 0.0_f8
       ENDIF
    ENDDO

    !========================================================================
    ! Copy quantities for UCX into gckpp_Global variables
    !========================================================================
    IF ( Input_Opt%LUCX ) THEN

       !---------------------------
       ! If UCX is turned on ...
       !---------------------------

       ! ... copy sticking coefficients for PSC reactions on SLA
       ! ... to the proper gckpp_Global variable
       KHETI_SLA  = State_Chm%KHETI_SLA(I,J,L,:)

       ! ... check if we are in the stratosphere
       STRATBOX   = State_Met%InStratosphere(I,J,L)

       ! ... check if there are solid PSCs at this grid box
       PSCBOX     = ( ( Input_Opt%LPSCCHEM                ) .and.            &
                      ( State_Chm%STATE_PSC(I,J,L) >= 2.0 ) .and. STRATBOX  )

       ! ... check if there is surface NAT at this grid box
       NATSURFACE = ( PSCBOX .and. ( C(ind_NIT) > 0.0_f8 )                  )

    ELSE

       !---------------------------
       ! If UCX is turned off ...
       !---------------------------

       ! ... set H2O concentration from the meteorology
       C(ind_H2O) = H2O

       ! ... zero out UCX-related quantities
       KHETI_SLA  = 0.0_f8
       STRATBOX   = .FALSE.
       PSCBOX     = .FALSE.
       NATSURFACE = .FALSE.

    ENDIF

    !========================================================================
    ! Calculate RH [%] and copy into gckpp_Global variables
    ! Not clear why this calc is slightly different than State_Met%RH
    !========================================================================
    CONSEXP  = 17.2693882_f8 * ( TEMP - 273.16_f8 ) / ( TEMP - 35.86_f8 )
    VPRESH2O = CONSVAP * EXP( CONSEXP ) / TEMP
    RELHUM   = ( H2O / VPRESH2O ) * 100_f8

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
!  But now, all species are stored in State_Chm%SPECIES, so the various
!  arrays (SAVEOH, SAVEHO2, etc.) are no longer necessary.  We can now just
!  just get values directly from State_Chm%SPECIES.
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
    REAL(fp), POINTER  :: AirNumDen(:,:,:  )
    REAL(fp), POINTER  :: Spc      (:,:,:,:)

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
               State_Diag%OHconcAfterChem(I,J,L) = Spc(I,J,L,id_OH)
            ENDIF
#ifdef MODEL_GEOS
            IF ( State_Diag%Archive_O3concAfterChem ) THEN
               State_Diag%O3concAfterChem(I,J,L) = Spc(I,J,L,id_O3)
            ENDIF
            IF ( Archive_RO2concAfterChem ) THEN
               IF ( id_A3O2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) +  Spc(I,J,L,id_A3O2)
               IF ( id_ATO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_ATO2)
               IF ( id_B3O2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_B3O2)
               IF ( id_BRO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_BRO2)
               IF ( id_ETO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_ETO2)
               IF ( id_HO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_HO2)
               IF ( id_LIMO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_LIMO2)
               IF ( id_MO2 > 0  ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_MO2)
               IF ( id_PIO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_PIO2)
               IF ( id_PO2 > 0  ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_PO2)
               IF ( id_PRN1 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_PRN1)
               IF ( id_R4N1 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_R4N1)
               IF ( id_R4O2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_R4O2)
               IF ( id_TRO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_TRO2)
               IF ( id_XRO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_XRO2)
               IF ( id_IHOO1 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_IHOO1)
               IF ( id_IHOO4 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_IHOO4)
               IF ( id_ICHOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_ICHOO)
               IF ( id_IHPOO1 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_IHPOO1)
               IF ( id_IHPOO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_IHPOO2)
               IF ( id_IHPOO3 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_IHPOO3)
               IF ( id_IEPOXAOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_IEPOXAOO)
               IF ( id_IEPOXBOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_IEPOXBOO)
               IF ( id_C4HVP1 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_C4HVP1)
               IF ( id_C4HVP2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_C4HVP2)
               IF ( id_HPALD1OO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_HPALD1OO)
               IF ( id_HPALD2OO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_HPALD2OO)
               IF ( id_ISOPNOO1 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_ISOPNOO1)
               IF ( id_ISOPNOO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_ISOPNOO2)
               IF ( id_INO2D > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_INO2D)
               IF ( id_INO2B > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_INO2B)
               IF ( id_IDHNBOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_IDHNBOO)
               IF ( id_IDHNDOO1 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_IDHNDOO1)
               IF ( id_IDHNDOO2 > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_IDHNDOO2)
               IF ( id_IHPNBOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_IHPNBOO)
               IF ( id_IHPNDOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_IHPNDOO)
               IF ( id_ICNOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_ICNOO)
               IF ( id_IDNOO > 0 ) &
                  State_Diag%RO2concAfterChem(I,J,L) = &
                     State_Diag%RO2concAfterChem(I,J,L) + Spc(I,J,L,id_IDNOO)
            ENDIF
#endif

         ENDIF

         !------------------------------------------------------------------
         ! HO2 concentration [v/v]
         !------------------------------------------------------------------
         IF ( ok_HO2 ) THEN
            IF ( State_Diag%Archive_HO2concAfterChem ) THEN
               State_Diag%HO2concAfterChem(I,J,L) = ( Spc(I,J,L,id_HO2)      &
                                                  /   AirNumDen(I,J,L)      )
            ENDIF

         ENDIF

         IF ( Input_Opt%LUCX ) THEN

            !---------------------------------------------------------------
            ! O1D concentration [molec/cm3]
            !---------------------------------------------------------------
            IF ( ok_O1D ) THEN
               IF ( State_Diag%Archive_O1DconcAfterChem ) THEN
                  State_Diag%O1DconcAfterChem(I,J,L) = Spc(I,J,L,id_O1D)
               ENDIF

            ENDIF


            !---------------------------------------------------------------
            ! O3P concentration [molec/cm3]
            !---------------------------------------------------------------
            IF ( ok_O3P ) THEN
               IF ( State_Diag%Archive_O3PconcAfterChem ) THEN
                  State_Diag%O3PconcAfterChem(I,J,L) = Spc(I,J,L,id_O3P)
               ENDIF

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
          OHconc_mcm3 = State_Chm%Species(I,J,L,id_OH)

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
! !IROUTINE: init_flexchem
!
! !DESCRIPTION: Subroutine Init\_FlexChem is used to allocate arrays for the
!  KPP solver.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_FlexChem( Input_Opt, State_Chm, State_Diag, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Gckpp_Monitor,    ONLY : Eqn_Names, Fam_Names
    USE Gckpp_Parameters, ONLY : nFam, nReact
    USE Input_Opt_Mod,    ONLY : OptInput
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Chm_Mod,    ONLY : Ind_
    USE State_Diag_Mod,   ONLY : DgnState
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
    ! Init_FlexChem begins here!
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
    ThisLoc  = ' -> at Init_FlexChem (in module GeosCore/flexchem_mod.F90)'
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Debug output
    IF ( prtDebug ) THEN
       WRITE( 6, 100 )
100    FORMAT( '     - INIT_FLEXCHEM: Allocating arrays for FLEX_CHEMISTRY' )

       WRITE( 6 ,'(a)' ) ' KPP Reaction Reference '
       DO N = 1, NREACT
          WRITE( 6, '(i8,a3,a85)' ) N,' | ',EQN_NAMES(N)
       END DO
    ENDIF

    ! Initialize species flags
    id_CH4                   = Ind_( 'CH4', 'A'     ) ! CH4 advected species
    id_HO2                   = Ind_( 'HO2'          )
    id_O3P                   = Ind_( 'O'            )
    id_O1D                   = Ind_( 'O1D'          )
    id_OH                    = Ind_( 'OH'           )
    id_SALA                  = Ind_( 'SALA'         )

#ifdef MODEL_GEOS
    ! ckeller
    id_O3                    = Ind_( 'O3'           )
    id_A3O2                  = Ind_( 'A3O2'         )
    id_ATO2                  = Ind_( 'ATO2'         )
    id_BRO2                  = Ind_( 'BRO2'         )
    id_ETO2                  = Ind_( 'ETO2'         )
    id_LIMO2                 = Ind_( 'LIMO2'        )
    id_MO2                   = Ind_( 'MO2'          )
    id_PIO2                  = Ind_( 'PIO2'         )
    id_PO2                   = Ind_( 'PO2'          )
    id_PRN1                  = Ind_( 'PRN1'         )
    id_R4N1                  = Ind_( 'R4N1'         )
    id_R4O2                  = Ind_( 'R4O2'         )
    id_TRO2                  = Ind_( 'TRO2'         )
    id_XRO2                  = Ind_( 'XRO2'         )
    id_IHOO1                 = Ind_( 'IHOO1'        )
    id_IHOO4                 = Ind_( 'IHOO4'        )
    id_IHCOO                 = Ind_( 'IHCOO'        )
    id_IHPOO1                = Ind_( 'IHPOO1'       )
    id_IHPOO2                = Ind_( 'IHPOO2'       )
    id_IHPOO3                = Ind_( 'IHPOO3'       )
    id_IEPOXAOO              = Ind_( 'IEPOXAOO'     )
    id_IEPOXBOO              = Ind_( 'IEPOXBOO'     )
    id_C4HVP1                = Ind_( 'C4HVP1'       )
    id_C4HVP2                = Ind_( 'C4HVP2'       )
    id_HPALD1OO              = Ind_( 'HPALD1OO'     )
    id_HPALD2OO              = Ind_( 'HPALD2OO'     )
    id_ISOPNOO1              = Ind_( 'ISOPNOO1'     )
    id_ISOPNOO2              = Ind_( 'ISOPNOO2'     )
    id_INO2B                 = Ind_( 'INO2B'        )
    id_INO2D                 = Ind_( 'INO2D'        )
    id_IDHNBOO               = Ind_( 'IDHNBOO'      )
    id_IDHNDOO1              = Ind_( 'IDHNDOO1'     )
    id_IDHNDOO2              = Ind_( 'IDHNDOO2'     )
    id_IHPNBOO               = Ind_( 'IHPNBOO'      )
    id_IHPNDOO               = Ind_( 'IHPNDOO'      )
    id_ICNOO                 = Ind_( 'ICNOO'        )
    id_IDNOO                 = Ind_( 'IDNOO'        )
#endif

    ! Set flags to denote if each species is defined
    ok_HO2                   = ( id_HO2 > 0         )
    ok_O1D                   = ( id_O1D > 0         )
    ok_O3P                   = ( id_O3P > 0         )
    ok_OH                    = ( id_OH  > 0         )

    ! Throw an error if certain diagnostics for UCX are turned on,
    ! but the UCX mechanism is not used in this fullchem simulation
    ! NOTE: Maybe eventually move this error check to state_diag_mod.F90
    IF ( .not. Input_Opt%LUCX ) THEN

       ! O1D diagnostic is only used w/ UCX
       IF ( State_Diag%Archive_O1DconcAfterChem ) THEN
          ErrMsg = 'The "O1DconcAfterChem" diagnostic is turned on ' //      &
                   'but the UCX mechanism is not being used!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! O3P diagnostic is only used w/ UCX
       IF ( State_Diag%Archive_O3PconcAfterChem ) THEN
          ErrMsg = 'The "O3PconcAfterChem" diagnostic is turned on ' //      &
                   'but the UCX mechanism is not being used!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

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
       CALL GC_CheckVar( 'flexchem_mod.F90:PL_Kpp_Id', 0, RC )
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

  END SUBROUTINE Init_FlexChem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_flexchem
!
! !DESCRIPTION: Subroutine Cleanup\_FlexChem deallocate module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_FlexChem( RC )
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
    ! Cleanup_FlexChem begins here!
    !=================================================================

    ! Initialize
    RC = GC_SUCCESS

    IF ( ALLOCATED( PL_Kpp_Id ) ) THEN
       DEALLOCATE( PL_Kpp_Id, STAT=RC  )
       CALL GC_CheckVar( 'flexchem_mod.F90:PL_Kpp_Id', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( JvCountDay ) ) THEN
       DEALLOCATE( JvCountDay, STAT=RC  )
       CALL GC_CheckVar( 'flexchem_mod.F90:JvCountDay', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( JvSumDay ) ) THEN
       DEALLOCATE( JvSumDay, STAT=RC  )
       CALL GC_CheckVar( 'flexchem_mod.F90:JvCountDay', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( JvCountMon ) ) THEN
       DEALLOCATE( JvCountMon, STAT=RC  )
       CALL GC_CheckVar( 'flexchem_mod.F90:JvCountMon', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

    IF ( ALLOCATED( JvSumMon ) ) THEN
       DEALLOCATE( JvSumMon, STAT=RC  )
       CALL GC_CheckVar( 'flexchem_mod.F90:JvCountMon', 2, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
    ENDIF

  END SUBROUTINE Cleanup_FlexChem
!EOC
END MODULE FlexChem_Mod
