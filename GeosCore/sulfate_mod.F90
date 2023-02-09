!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: sulfate_mod.F90
!
! !DESCRIPTION: Module SULFATE\_MOD contains arrays and routines for performing
!  either a coupled chemistry/aerosol run or an offline sulfate aerosol
!  simulation. Original code taken from Mian Chin's GOCART model and modified
!  accordingly. (rjp, bdf, bmy, 6/22/00, 8/26/10)
!\\
!\\
! !INTERFACE:
!
MODULE SULFATE_MOD
!
! !USES:
!
  USE HCO_ERROR_MOD    ! For HEMCO error handling
  USE PhysConstants    ! Physical constants
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp, f4, f8)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
#ifdef APM
  PRIVATE :: WET_SETTLINGBIN
#endif
  PUBLIC :: CHEMSULFATE
  PUBLIC :: CLEANUP_SULFATE
  PUBLIC :: INIT_SULFATE
#ifdef TOMAS
  PUBLIC :: EMISSSULFATETOMAS
#endif

!
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Andreae, M.O. & P. Merlet, "Emission of trace gases and aerosols from
!        biomass burning", Global Biogeochem. Cycles, 15, 955-966, 2001.
!  (2 ) Nightingale et al [2000a], J. Geophys. Res, 14, 373-387
!  (3 ) Nightingale et al [2000b], Geophys. Res. Lett, 27, 2117-2120
!  (4 ) Wanninkhof, R., "Relation between wind speed and gas exchange over
!        the ocean", J. Geophys. Res, 97, 7373-7382, 1992.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  !========================================================================
  ! MODULE PARAMETERS:
  !
  ! XNUMOL_OH  : Molecules OH  per kg OH          [molec/kg]
  ! XNUMOL_O3  : Molecules O3  per kg O3          [molec/kg]
  ! XNUMOL_NO3 : Molecules NO3 per kg NO3         [molec/kg]
  ! TCVV_S     : Ratio: Molwt air / Molwt S       [unitless]
  !=======================================================================
  REAL(fp),  PARAMETER   :: XNUMOL_OH   = AVO / 17e-3_fp ! hard-coded MW
  REAL(fp),  PARAMETER   :: XNUMOL_O3   = AVO / 48e-3_fp ! hard-coded MW
  REAL(fp),  PARAMETER   :: XNUMOL_NO3  = AVO / 62e-3_fp ! hard-coded MW
  REAL(fp),  PARAMETER   :: XNUMOL_H2O2 = AVO / 34e-3_fp ! hard-coded MW
  REAL(fp),  PARAMETER   :: TCVV_S      = AIRMW / 32e+0_fp ! hard-coded MW
  REAL(fp),  PARAMETER   :: TCVV_N      = AIRMW / 14e+0_fp ! hard-coded MW
  REAL(fp),  PARAMETER   :: SMALLNUM    = 1e-20_fp
  REAL(fp),  PARAMETER   :: CM3PERM3    = 1.e6_fp

#ifdef TOMAS
  !---------------------------------------------------------------
  ! For TOMAS microphysics: Add parameter for scaling anthro SO2
  !---------------------------------------------------------------
  REAL(fp),  PARAMETER   :: scaleanthso2 = 1.0e+0_fp
#endif
!
! !PRIVATE TYPES:
!
  !========================================================================
  ! MODULE VARIABLES:
  !
  ! DMSo       : DMS oceanic emissions            [v/v/timestep]
  ! DRYSO4s    : Pointer to SO4s  in DEPVEL array [unitless]
  ! DRYNITs    : Pointer to NITs  in DEPVEL array [unitless]
  !
  !%%% NOTE: THESE ARE NOW OBTAINED VIA HEMCO (bmy, 5/22/15) %%%%%%%%%%%
  !% ENH3_an    : NH3 anthropogenic emissions      [kg NH3/box/s]
  !% ENH3_bb    : NH3 biomass emissions            [kg NH3/box/s]
  !% ENH3_bf    : NH3 biofuel emissions            [kg NH3/box/s]
  !% ENH3_na    : NH73 natural source emissions    [kg NH3/box/s]
  !% ESO2_ac    : SO2 aircraft emissions           [kg SO2/box/s]
  !% ESO2_an    : SO2 anthropogenic emissions      [kg SO2/box/s]
  !% ESO2_ev    : SO2 eruptive volcanic em.        [kg SO2/box/s]
  !% ESO2_nv    : SO2 non-eruptive volcanic em.    [kg SO2/box/s]
  !% ESO2_bb    : SO2 biomass burning emissions    [kg SO2/box/s]
  !% ESO2_bf    : SO2 biofuel burning emissions    [kg SO2/box/s]
  !% ESO2_sh    : SO2 ship emissions               [kg SO2/box/s]
  !% ESO4_an    : SO4 anthropogenic emissions      [kg SO4/box/s]
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! PMSA_DMS   : P(MSA) from DMS                  [v/v/timestep]
  ! PSO2_DMS   : P(SO2) from DMS                  [v/v/timestep]
  ! PSO4_SO2   : P(SO4) from SO2                  [v/v/timestep]
  ! PHMS_SO2   : P(HMS) from SO2                  [v/v/timestep] (jmm, 6/15/18)
  ! PSO2_HMS   : P(SO2) from HMS                  [v/v/timestep] (jmm, 6/15/18)
  ! PSO4_HMS   : P(SO4) from HMS & radical chem   [v/v/timestep] (jmm, 6/26/18)
  ! SSTEMP     : Sea surface temperatures         [K]
  ! Eev        : SO2 em. from eruptive volcanoes  [kg SO2/box/s]
  ! Env        : SO2 em. from non-erup volcanoes  [kg SO2/box/s]
  ! TCOSZ      : Sum of cos(SZA) for offline run  [unitless]
  ! TTDAY      : Total daylight length at (I,J)   [minutes]
  ! SMALLNUM   : Small number - prevent underflow [unitless]
  ! COSZM      : Array for MAX(cos(SZA)) at (I,J) [unitless]
  ! LVOLC      : Number of volcanic levels (20)   [unitless]
  !========================================================================

  ! Time variable
  INTEGER                :: ELAPSED_SEC

  ! Allocatable arrays
  REAL(fp),  ALLOCATABLE :: DMSo(:,:)
  REAL(fp),  ALLOCATABLE :: PMSA_DMS(:,:,:)
  REAL(fp),  ALLOCATABLE :: PSO2_DMS(:,:,:)
  REAL(fp),  ALLOCATABLE :: PSO4_SO2(:,:,:)
  REAL(fp),  ALLOCATABLE :: PSO4_SS(:,:,:)
  REAL(fp),  ALLOCATABLE :: PHMS_SO2(:,:,:)    ! jmm 06/13/2018
  REAL(fp),  ALLOCATABLE :: PSO2_HMS(:,:,:)    ! jmm 06/13/2018
  REAL(fp),  ALLOCATABLE :: PSO4_HMS(:,:,:)    ! jmm 06/26/2018
  REAL(fp),  ALLOCATABLE :: PNITs(:,:,:)
  REAL(fp),  ALLOCATABLE :: PNIT_dust(:,:,:,:) ! tdf
  REAL(fp),  ALLOCATABLE :: PSO4_dust(:,:,:,:) ! tdf
  REAL(f4),  ALLOCATABLE :: SOx_SCALE(:,:)
  REAL(fp),  ALLOCATABLE :: SSTEMP(:,:)
  REAL(fp),  ALLOCATABLE :: TCOSZ(:,:)
  REAL(fp),  ALLOCATABLE :: TTDAY(:,:)
  REAL(fp),  ALLOCATABLE :: COSZM(:,:)
  REAL(fp),  ALLOCATABLE :: GLOBAL_OH(:,:,:)
  REAL(fp),  ALLOCATABLE :: GLOBAL_HNO3(:,:,:)
  REAL(fp),  ALLOCATABLE :: GLOBAL_HCl(:,:,:)
  REAL(fp),  ALLOCATABLE :: GLOBAL_HCOOH(:,:,:)
  REAL(fp),  ALLOCATABLE :: GLOBAL_ACTA(:,:,:)
  REAL(fp),  ALLOCATABLE :: PNIT(:,:,:) ! xnw
  REAL(fp),  ALLOCATABLE :: PACL(:,:,:) ! xnw
  REAL(fp),  ALLOCATABLE :: PCCL(:,:,:) ! xnw

#ifdef APM
  REAL(fp),  ALLOCATABLE :: PSO4_SO2APM(:,:,:)
  REAL(fp),  ALLOCATABLE :: PSO4_SO2SEA(:,:,:)
#endif
#ifdef TOMAS
  !---------------------------------------------------------------
  ! For TOMAS microphysics: Define PSO4_SO2aq array
  !---------------------------------------------------------------
  REAL(fp),  ALLOCATABLE :: PSO4_SO2AQ(:,:,:)
  REAL(fp),  ALLOCATABLE :: SO4_ANTH(:,:,:)
#endif

  ! These are pointers to fields in the HEMCO data structure.
  ! Declare these with REAL(fp), aka REAL*4. (bmy, 3/4/15)
  REAL(f4), POINTER      :: OH(:,:,:)       => NULL()
  REAL(f4), POINTER      :: NDENS_SALA(:,:) => NULL()
  REAL(f4), POINTER      :: NDENS_SALC(:,:) => NULL()

  ! Emission timestep (imported from HEMCO)
  REAL(fp)               :: TS_EMIS

  ! Species ID flags
  INTEGER                :: id_AS,     id_AHS,    id_AW1
  INTEGER                :: id_DAL1,   id_DAL2,   id_DAL3
  INTEGER                :: id_DAL4,   id_DMS,    id_DST1
  INTEGER                :: id_DST2,   id_DST3,   id_DST4
  INTEGER                :: id_H2O2,   id_HNO3,   id_LET
  INTEGER                :: id_MSA,    id_NH3,    id_NH4
  INTEGER                :: id_NH4aq,  id_NIT,    id_NITd1
  INTEGER                :: id_NITd2,  id_NITd3,  id_NITd4
  INTEGER                :: id_NITs,   id_NK1,    id_NK5
  INTEGER                :: id_NK8,    id_NK10,   id_NK20
  INTEGER                :: id_NO3,    id_O3,     id_OH
  INTEGER                :: id_SALA,   id_SALC,   id_SF1
  INTEGER                :: id_SO2,    id_SO4,    id_SO4aq
  INTEGER                :: id_SO4d1,  id_SO4d2,  id_SO4d3
  INTEGER                :: id_SO4d4,  id_SO4s,   id_pFe
  INTEGER                :: id_SALACL, id_HCL,    id_SALCCL
  INTEGER                :: id_SALAAL, id_SALCAL
  INTEGER                :: id_HOBr,   id_SO4H1,  id_SO4H2
  INTEGER                :: id_HOCl,   id_SO4H3,  id_SO4H4
  INTEGER                :: id_HCOOH,  id_ACTA,   id_PSO4
  INTEGER                :: id_HMS,    id_CH2O


  ! Species drydep ID flags
  INTEGER                :: DRYSO4s,   DRYNITs,   DRYSO4d1
  INTEGER                :: DRYSO4d2,  DRYSO4d3,  DRYSO4d4
  INTEGER                :: DRYNITd1,  DRYNITd2,  DRYNITd3
  INTEGER                :: DRYNITd4

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemsulfate
!
! !DESCRIPTION: Subroutine CHEMSULFATE is the interface between the GEOS-CHEM
!  main program and the sulfate chemistry routines.  The user has the option of
!  running a coupled chemistry-aerosols simulation or an offline aerosol
!  simulation. (rjp, bdf, bmy, 5/31/00, 3/16/06)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEMSULFATE( Input_Opt,  State_Chm, State_Diag, State_Grid, &
                          State_Met,  FullRun,   RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD
    USE HCO_State_GC_Mod,   ONLY : HcoState
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : SpcConc
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_MONTH
    USE TIME_MOD,           ONLY : GET_TS_CHEM
    USE TIME_MOD,           ONLY : GET_ELAPSED_SEC
    USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
    USE UCX_MOD,            ONLY : SETTLE_STRAT_AER
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
#ifdef APM
    USE HCO_STATE_MOD,      ONLY : HCO_GetHcoID
    USE APM_DRIV_MOD,       ONLY : EMITNH3,EMITSO2
#endif
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: FullRun     ! Modify species conc?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE            :: FIRSTCHEM = .TRUE.
    INTEGER, SAVE            :: LASTMONTH = -99

    ! Scalars
    LOGICAL                  :: LGRAVSTRAT
    LOGICAL                  :: LDSTUP
    LOGICAL                  :: prtDebug
    INTEGER                  :: I, J, L, N, MONTH
    REAL(fp)                 :: DTCHEM
    CHARACTER(LEN=63)        :: OrigUnit

    ! Strings
    CHARACTER(LEN=255)       :: ErrMsg, ThisLoc

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)
#ifdef APM
    INTEGER           :: IDNH3,IDSO2
    REAL(fp)          :: A_M2
#endif

    !=================================================================
    ! CHEMSULFATE begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at CHEMSULFATE (in module GeosCore/sulfate_mod.F90)'

      ! Copy fields from INPUT_OPT to local variables for use below
    LGRAVSTRAT           = Input_Opt%LGRAVSTRAT
    LDSTUP               = Input_Opt%LDSTUP

    ! Initialize pointers
    Spc                  => State_Chm%Species  ! Chemistry species [kg]

    ! Should we print debug output?
    prtDebug             = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Get current month
    MONTH                = GET_MONTH()

    ! If it's an offline simulation ...
    IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       ! Evaluate offline global OH from HEMCO
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'GLOBAL_OH', GLOBAL_OH, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get data for GLOBAL_OH from HEMCO!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Evaluate offline global HNO3 from HEMCO
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'GLOBAL_HNO3', GLOBAL_HNO3, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get data for GLOBAL_HNO3 from HEMCO!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Evaluate offline global HCl from HEMCO
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'GLOBAL_HCl', GLOBAL_HCl, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get data for GLOBAL_HCl from HEMCO!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! And compute time scaling arrays for offline OH, NO3
       CALL OHNO3TIME( State_Grid )

    ENDIF

    ! Store NTIME in a shadow variable
    ELAPSED_SEC = GET_ELAPSED_SEC()

    ! DTCHEM is the chemistry timestep in seconds
    DTCHEM = GET_TS_CHEM()

    ! TS_EMIS is the emission timestep (in seconds). This is a module
    ! variable, hence define only on first call.
    IF ( FIRSTCHEM ) THEN
#if defined( MODEL_CESM )
       ! Do not use HEMCO state in CESM
       TS_EMIS = REAL( Input_Opt%TS_EMIS, fp )
#else

       IF ( .NOT. ASSOCIATED(HcoState) ) THEN
          ErrMsg = 'Cannot get HEMCO state variable "HCOState"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          IF ( RC /= GC_SUCCESS ) RETURN
       ENDIF
       TS_EMIS = HcoState%TS_EMIS
#endif
    ENDIF

#ifdef APM
    IDNH3 = HCO_GetHcoID( 'NH3', HcoState )
    !$OMP PARALLEL DO              &
    !$OMP DEFAULT( SHARED )        &
    !$OMP PRIVATE( L, J, I, A_M2 ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Grid box surface area [m2]
       A_M2 = HcoState%Grid%AREA_M2%Val( I, J )

       DO L = 1, State_Grid%NZ
          ! Get emissions [kg/m2/s] and convert to [kg/box-sec]
          EMITNH3(I,J,L) = HcoState%Spc(IDNH3)%Emis%Val(I,J,L)*A_M2
       ENDDO

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    IDSO2 = HCO_GetHcoID( 'SO2', HcoState )
    !$OMP PARALLEL DO              &
    !$OMP DEFAULT( SHARED )        &
    !$OMP PRIVATE( L, J, I, A_M2 ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Grid box surface area [m2]
       A_M2 = HcoState%Grid%AREA_M2%Val( I, J )

       DO L = 1, State_Grid%NZ
          ! Get emissions [kg/m2/s] and convert to [kg/box-sec]
          EMITSO2(I,J,L) = HcoState%Spc(IDSO2)%Emis%Val(I,J,L)*A_M2
       ENDDO

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO
#endif

    ! Initialize module arrays
    PSO2_DMS   = 0e+0_fp
    PMSA_DMS   = 0e+0_fp
    PSO4_SO2   = 0e+0_fp
    PHMS_SO2   = 0e+0_fp     ! jmm 06/13/2018
    PSO2_HMS   = 0e+0_fp     ! jmm 06/13/2018
    PSO4_HMS   = 0e+0_fp     ! jmm 06/26/2018
    PSO4_SS    = 0e+0_fp
    PNITs      = 0e+0_fp
    PSO4_dust  = 0e+0_fp     ! tdf 04/17/08
    PNIT_dust  = 0e+0_fp     ! tdf 04/17/08
    PNIT       = 0e+0_fp
    PACL       = 0e+0_fp
    PCCL       = 0e+0_fp
#ifdef APM
    PSO4_SO2APM = 0d0
    PSO4_SO2SEA = 0d0
#endif
#ifdef TOMAS
    PSO4_SO2AQ = 0e+0_fp     ! For TOMAS microphysics
#endif

    !=========================================================================
    ! Call individual chemistry routines for sulfate/aerosol speccies
    !=========================================================================

    ! Perform all routines only when doing a "full" run
    IF ( FullRun ) THEN

       !---------------------------------------------------------------------
       ! FullRun = T: Do all sulfate chemistry
       !---------------------------------------------------------------------

       ! SO4s [kg] gravitational settling
       IF ( id_SO4s > 0 ) THEN
          CALL GRAV_SETTLING( Input_Opt, State_Chm, State_Diag, State_Grid, &
                              State_Met, id_SO4s,   RC )
          IF ( prtDebug ) THEN
             CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, SO4S' )
          ENDIF
       ENDIF

       ! NITs [kg] gravitational settling
       IF ( id_NITs > 0 ) THEN
          CALL GRAV_SETTLING( Input_Opt, State_Chm, State_Diag, State_Grid, &
                              State_Met, id_NITs,   RC )
          IF ( prtDebug ) THEN
             CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, NITS' )
          ENDIF
       ENDIF

       !----------------------------------------------------------------
       ! These species are only used for the aciduptake simulations
       !----------------------------------------------------------------
       IF ( LDSTUP ) THEN

          ! SO4d1 [kg] gravitational settling
          IF ( id_SO4d1 > 0 ) THEN
             CALL GRAV_SETTLING( Input_Opt, State_Chm, State_Diag, State_Grid, &
                                 State_Met, id_SO4d1,  RC )
             IF ( prtDebug ) THEN
                CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, SO4d1')
             ENDIF
          ENDIF

          ! SO4d2 [kg] gravitational settling
          IF ( id_SO4d2 > 0 ) THEN
             CALL GRAV_SETTLING( Input_Opt, State_Chm, State_Diag, State_Grid, &
                                 State_Met, id_SO4d2,  RC )
             IF ( prtDebug ) THEN
                CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, SO4d2')
             ENDIF
          ENDIF

          ! SO4d3 [kg] gravitational settling
          IF ( id_SO4d3 > 0 ) THEN
             CALL GRAV_SETTLING( Input_Opt, State_Chm, State_Diag, State_Grid, &
                                 State_Met, id_SO4d3,  RC )
             IF ( prtDebug ) THEN
                CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, SO4d3')
             ENDIF
          ENDIF

          ! SO4d4 [kg] gravitational settling
          IF ( id_SO4d4 > 0 ) THEN
             CALL GRAV_SETTLING( Input_Opt, State_Chm, State_Diag, State_Grid, &
                                 State_Met, id_SO4d4,  RC )
             IF ( prtDebug ) THEN
                CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, SO4d4')
             ENDIF
          ENDIF

          ! NITd1 [kg] gravitational settling
          IF ( id_NITd1 > 0 ) THEN
             CALL GRAV_SETTLING( Input_Opt, State_Chm, State_Diag, State_Grid, &
                                 State_Met, id_NITd1,  RC )
             IF ( prtDebug ) THEN
                CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, NITd1')
             ENDIF
          ENDIF

          ! NITd2 [kg] gravitational settling
          IF ( id_NITd2 > 0 ) THEN
             CALL GRAV_SETTLING( Input_Opt, State_Chm, State_Diag, State_Grid, &
                                 State_Met, id_NITd2,  RC )
             IF ( prtDebug ) THEN
                CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, NITd2')
             ENDIF
          ENDIF

          ! NITd3 [kg] gravitational settling
          IF ( id_NITd3 > 0 ) THEN
             CALL GRAV_SETTLING( Input_Opt, State_Chm, State_Diag, State_Grid, &
                                 State_Met, id_NITd3,  RC )
             IF ( prtDebug ) THEN
                CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, NITd3')
             ENDIF
          ENDIF

          ! NITd4 [kg] gravitational settling
          IF ( id_NITd4 > 0 ) THEN
             CALL GRAV_SETTLING( Input_Opt, State_Chm, State_Diag, State_Grid, &
                                 State_Met, id_NITd4,  RC )
             IF ( prtDebug ) THEN
                CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, NITd4')
             ENDIF
          ENDIF
       ENDIF

       ! Stratospheric aerosol gravitational settling
       IF ( LGRAVSTRAT ) THEN
          CALL SETTLE_STRAT_AER( Input_Opt, State_Chm, State_Grid, &
                                 State_Met, RC )
          IF ( prtDebug ) THEN
             CALL DEBUG_MSG( '### CHEMSULFATE: GRAV_SET, STRAT' )
          ENDIF
       ENDIF

       ! Convert species to [v/v dry]
       CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                               'v/v dry', RC, OrigUnit=OrigUnit )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error('Unit conversion error', RC, &
                        'Start of CHEM_SULFATE in sulfate_mod.F90')
          RETURN
       ENDIF
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMSULFATE: a CONVERT UNITS' )
       ENDIF

       ! For offline runs only ...
       IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

          !--------------------------------
          ! DMS chemistry (offline only)
          !--------------------------------
          CALL CHEM_DMS( Input_Opt, State_Chm, State_Diag, State_Grid, &
                         State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_DMS"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Debug info
          IF ( prtDebug ) THEN
             CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_DMS' )
          ENDIF

          !--------------------------------
          ! H2O2 (offline only)
          !--------------------------------
          CALL CHEM_H2O2( Input_Opt, State_Chm, State_Diag, State_Grid, &
                          State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "Chem_H2O2"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          IF ( prtDebug ) THEN
             CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_H2O2' )
          ENDIF

       ENDIF

       !-----------------------
       ! SO2 chemistry
       !-----------------------
       CALL CHEM_SO2( Input_Opt, State_Chm, State_Diag, State_Grid, &
                      State_Met, .TRUE.,    RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Chem_SO2"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Debug info
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_SO2' )
       ENDIF

       !-----------------------
       ! SO4 chemistry
       !-----------------------
       CALL CHEM_SO4( Input_Opt,  State_Chm, State_Diag, State_Grid, &
                      State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Chem_SO4"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_SO4' )
       ENDIF

#ifdef TOMAS
       !-----------------------------------------------------------------
       ! For TOMAS microphysics:
       !
       ! SO4 from aqueous chemistry of SO2 (in-cloud oxidation)
       !
       ! SO4 produced via aqueous chemistry is distributed onto 30-bin
       ! aerosol by TOMAS subroutine AQOXID.   NOTE: This may be moved
       ! to tomas_mod.F90 in the future, but for now it still needs to get
       ! the PSO4_SO2AQ value while CHEMSULFATE is called
       !-----------------------------------------------------------------
       CALL CHEM_SO4_AQ( Input_Opt, State_Chm, State_Grid, State_Met, RC )
       IF ( prtDebug ) CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_SO4_AQ' )
#endif

       ! MSA
       CALL CHEM_MSA( Input_Opt, State_Chm, State_Grid, State_Met, RC )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_MSA' )
       ENDIF

       ! Sulfur Nitrate.
       ! CHEM_NIT includes a source term from sea salt aerosols, so keep
       ! here.
       CALL CHEM_NIT( Input_Opt, State_Chm, State_Grid, State_Met, RC )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_NIT' )
       ENDIF

       ! Calculate the HCl uptake by alkalinity, xnw
       CALL CHEM_CL( Input_Opt, State_Met, State_Chm, State_Grid, RC )
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_CL' )
       ENDIF

    ELSE

       !---------------------------------------------------------------------
       ! FullRun = F: Just set up Cloud pH & related parameters, and exit
       !---------------------------------------------------------------------

       ! Convert species to [v/v dry]
       CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met,  &
                               'v/v dry', RC, OrigUnit=OrigUnit )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL GC_Error('Unit conversion error', RC, &
                        'Start of CHEM_SULFATE in sulfate_mod.F90')
          RETURN
       ENDIF
       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMSULFATE: a CONVERT UNITS' )
       ENDIF

       ! Call the SO2 routine to get cloud pH parameters
       CALL CHEM_SO2( Input_Opt, State_Chm, State_Diag, State_Grid,          &
                      State_Met, .FALSE.,   RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Chem_SO2"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       IF ( prtDebug ) THEN
          CALL DEBUG_MSG( '### CHEMSULFATE: a CHEM_SO2 false' )
       ENDIF

    ENDIF ! FullRun

    ! Convert species units back to original unit
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met,     &
                            OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'End of CHEM_SULFATE in sulfate_mod.F90')
       RETURN
    ENDIF

    ! Free pointer
    Spc => NULL()

    ! We have already gone thru one chemistry iteration
    FIRSTCHEM = .FALSE.

  END SUBROUTINE CHEMSULFATE
!EOC
#ifdef TOMAS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emisssulfatetomas
!
! !DESCRIPTION: Subroutine EMISSSULFATETOMAS connects HEMCO bulk emissions to
! the TOMAS tracers. Only use this for TOMAS sims. This should be quite similar
! to the TOMAS relevant parts of 'emisssulfate' in v9 (Jkodros 6/2/15)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMISSSULFATETOMAS( Input_Opt, State_Chm, State_Grid, &
                                State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : SpcConc
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TOMAS_MOD,          ONLY : IBINS,      ICOMP,   IDIAG
    USE TOMAS_MOD,          ONLY : NH4BULKTOBIN
    USE TOMAS_MOD,          ONLY : SRTNH4
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State objectt
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    ! Fields for TOMAS simulation
    REAL*8           :: BINMASS(State_Grid%NX,State_Grid%NY,State_Grid%NZ, &
                                IBINS*ICOMP)
    INTEGER          :: TID, I, J, L, M
    INTEGER          :: ii=53, jj=29, ll=1
    REAL(fp)         :: NH4_CONC
    CHARACTER(LEN=63):: OrigUnit

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    ! Arrays
    REAL(fp)         :: tempnh4(ibins)
    REAL(fp)         :: MK_TEMP2(IBINS)

    !=================================================================
    ! EMISSSULFATETOMAS begins here!
    !=================================================================

    ! Convert species to [kg] for TOMAS. This will be removed once
    ! TOMAS uses mixing ratio instead of mass as tracer units (ewl, 9/11/15)
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'kg', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'Start of EMISSSULFATETOMAS in sulfate_mod.F90')
       RETURN
    ENDIF

    ! Point to chemical species array [kg]
    Spc => State_Chm%Species

    IF (id_SF1 > 0 .and. id_NK1 > 0 ) THEN

       ! Get NH4 and aerosol water into the same array
       DO M = 1, IBINS*(ICOMP-IDIAG)
          BINMASS(:,:,:,M) = Spc(id_SF1+M-1)%Conc(:,:,:)
       ENDDO

       IF ( SRTNH4 > 0 ) THEN
          TID = IBINS*(ICOMP-IDIAG) + 1

          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( I, J, L, M, TEMPNH4, MK_TEMP2, NH4_CONC ) &
          !$OMP SCHEDULE( DYNAMIC )
          DO L=1,State_Grid%NZ
          DO J=1,State_Grid%NY
          DO I=1,State_Grid%NX

             ! Change pointer to a variable to avoid array temporary
             ! (bmy, 7/7/17)
             DO M = 1, IBINS
                MK_TEMP2(M) = Spc(id_SF1+M-1)%Conc(I,J,L)
             ENDDO
             NH4_CONC = Spc(id_NH4)%Conc(I,J,L)
             CALL NH4BULKTOBIN( MK_TEMP2, NH4_CONC, TEMPNH4 )

             BINMASS(I,J,L,TID:TID+IBINS-1) = TEMPNH4(1:IBINS)

          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ENDIF

       TID = IBINS*(ICOMP-1) +1
       DO M = 1, IBINS
          BINMASS(:,:,:,TID+M-1) = Spc(id_AW1+M-1)%Conc(:,:,:)
       ENDDO

       !IF ( id_SF1 > 0 ) THEN
       CALL SRCSF30( Input_Opt, State_Grid, State_Met, &
                     State_Chm, BINMASS(:,:,:,:), RC )

       ! Return the aerosol mass after emission subroutine to Spc
       ! excluding the NH4 aerosol and aerosol water (win, 9/27/08)
       DO M = 1, IBINS*(ICOMP-IDIAG)
          Spc(id_SF1+M-1)%Conc(:,:,:) = BINMASS(:,:,:,M)
       ENDDO
    ENDIF

    ! Free pointer
    NULLIFY( Spc )

    ! Convert species back to original units (ewl, 9/11/15)
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'End of EMISSSULFATETOMAS in sulfate_mod.F90')
       RETURN
    ENDIF

  END SUBROUTINE EMISSSULFATETOMAS
!EOC
#endif
!-----------------------------------------------------------------------------
!                  Jack Kodros re-writing this
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE: srcsf30
!
! !DESCRIPTION: Subroutine SRCSF30 (Jkodros 6/2/15)
!\\
!\\
! !INTERFACE:
!
#ifdef TOMAS
  SUBROUTINE SRCSF30( Input_Opt, State_Grid, State_Met, State_Chm, TC2, RC )
!
! !USES:
!
#ifdef BPCH_DIAG
    USE CMN_DIAG_MOD             ! ND13 (for now)
    USE DIAG_MOD,             ONLY : AD59_SULF,     AD59_NUMB
#endif
    USE ERROR_MOD,            ONLY : ERROR_STOP,  IT_IS_NAN
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Species_Mod,          ONLY : SpcConc
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE State_Chm_Mod,        ONLY : ChmState
    USE TOMAS_MOD,            ONLY : IBINS, AVGMASS, ICOMP
    USE TOMAS_MOD,            ONLY : Xk
    USE TOMAS_MOD,            ONLY : SUBGRIDCOAG, MNFIX
    USE TOMAS_MOD,            ONLY : SRTSO4, SRTNH4,  DEBUGPRINT

    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_GetDiagn
    USE HCO_State_GC_Mod,     ONLY : HcoState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
    REAL(fp), INTENT(INOUT) :: TC2(State_Grid%NX,State_Grid%NY,State_Grid%NZ, &
                                   IBINS*ICOMP)
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    INTEGER                :: I, J, K, L, DOW_LT, NTOP, C, Bi
    REAL*8                 :: SO4(State_Grid%NZ)
    REAL*8                 :: DTSRCE
    REAL*8                 :: EFRAC(State_Grid%NZ)
    REAL*8                 :: TSO4,       FEMIS
    REAL*8                 :: AREA_CM2
    REAL*8                 :: SO4an(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL*8                 :: SO4bf(State_Grid%NX,State_Grid%NY)
    REAL*8                 :: SO4anbf(State_Grid%NX,State_Grid%NY,2)
    REAL*8 BFRAC(IBINS)     ! Mass fraction emitted to each bin

#if  defined( TOMAS12) || defined( TOMAS15)
    DATA BFRAC/ &
# if  defined( TOMAS15)
         0.0d0     , 0.0d0     , 0.0d0, &
# endif
         4.3760E-02, 6.2140E-02, 3.6990E-02, 1.8270E-02, &
         4.2720E-02, 1.1251E-01, 1.9552E-01, 2.2060E-01, &
         1.6158E-01, 7.6810E-02, 2.8884E-02, 2.0027E-04/

#else

    DATA BFRAC/ &
# if  defined( TOMAS40)
         0.000d00 , 0.000d00 , 0.000d00 , 0.000d00 , 0.000d00 , &
         0.000d00 , 0.000d00 , 0.000d00 , 0.000d00 , 0.000d00 , &
# endif
         1.728d-02, 2.648d-02, 3.190d-02, 3.024d-02, 2.277d-02, &
         1.422d-02, 9.029d-03, 9.241d-03, 1.531d-02, 2.741d-02, &
         4.529d-02, 6.722d-02, 8.932d-02, 1.062d-01, 1.130d-01, &
         1.076d-01, 9.168d-02, 6.990d-02, 4.769d-02, 2.912d-02, &
         1.591d-02, 7.776d-03, 3.401d-03, 1.331d-03, 4.664d-04, &
         1.462d-04, 4.100d-05, 1.029d-05, 2.311d-06, 4.645d-07/
#endif

    REAL(fp)               :: NDISTINIT(IBINS)
    REAL(fp)               :: NDISTFINAL(IBINS)
    REAL(fp)               :: MADDFINAL(IBINS)
    REAL(fp)               :: NDIST(IBINS),  MDIST(IBINS,ICOMP)
    REAL(fp)               :: NDIST2(IBINS), MDIST2(IBINS,ICOMP)
    REAL*4                 :: TSCALE, BOXVOL, TEMP, PRES

    !REAL(fp)              :: N0(State_Grid%NZ,IBINS)
    !REAL(fp)              :: M0(State_Grid%NZ,IBINS)
    REAL(fp)               :: Ndiag(IBINS), Mdiag(IBINS)
    REAL(fp)               :: MADDTOTAL !optimization variable for diag
    !REAL(fp)              :: Avginit(IBINS), Avgfinal(IBINS)
    !REAL(fp)              :: Avginner(IBINS)

    REAL(fp)               :: AREA(State_Grid%NX, State_Grid%NY)
    !REAL(fp)              :: AREA3D(State_Grid%NX, State_Grid%NY,2)

    ! Pointers
    REAL(f4),      POINTER :: Ptr2D(:,:  )
    REAL(f4),      POINTER :: Ptr3D(:,:,:)
    TYPE(SpcConc), POINTER :: TC1(:)

    INTEGER                :: N_TRACERS

    LOGICAL                :: ERRORSWITCH, SGCOAG = .TRUE.
    INTEGER                :: FLAG, ERR
    logical                :: pdbug !(temporary) win, 10/24/07
    !integer               :: ii, jj, ll
    !data ii, jj, ll / 61, 1, 7 /
    INTEGER                :: ii=53, jj=29, ll=1

    ! Ratio of molecular weights: S/SO4
    REAL*8,  PARAMETER     :: S_SO4 = 32d0 / 96d0

    ! debugging
    real*8   dummy

    ! For fields from Input_Opt
    LOGICAL :: prtDebug, LNLPBL
    LOGICAL :: jkdbg=.true.

    ! Strings
    CHARACTER(LEN= 63)       :: DgnName
    CHARACTER(LEN=255)       :: MSG
    CHARACTER(LEN=255)       :: LOC='srcsf30 (sulfate_mod.F90)'

    !=================================================================
    ! SRCSF30 begins here!
    !=================================================================

    ! Free pointers
    Ptr2D    => NULL()
    Ptr3D    => NULL()

    ! COpy values from Input_Opt
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    LNLPBL   = Input_Opt%LNLPBL

    ! Import emissions from HEMCO (through HEMCO state)
    IF ( .NOT. ASSOCIATED(HcoState) ) THEN
       CALL ERROR_STOP ( 'HcoState not defined!', LOC )
    ENDIF

    ! Emission timestep [seconds]
    DTSRCE = HcoState%TS_EMIS

    ! Grid box aarea
    AREA = HcoState%Grid%AREA_M2%Val(:,:)
    !AREA3D(:,:,1) = AREA(:,:)
    !AREA3D(:,:,2) = AREA(:,:)

    ! Define subgrid coagulation timescale (win, 10/28/08)
    IF ( TRIM(State_Grid%GridRes) == '4.0x5.0' ) THEN
       TSCALE = 10.*3600.  ! 10 hours
    ELSE IF ( TRIM(State_Grid%GridRes) == '2.0x2.5' ) THEN
       TSCALE = 5.*3600.
    ELSE IF ( TRIM(State_Grid%GridRes) == '0.5x0.625' ) THEN
       TSCALE = 1.*3600.
    ELSE IF ( TRIM(State_Grid%GridRes) == '0.25x0.3125' ) THEN
       TSCALE = 0.5*3600.
    ENDIF

    ! Point to species array
    TC1 => State_Chm%Species

    !================================================================
    ! READ IN HEMCO EMISSIONS
    !================================================================
    DgnName = 'SO4_ANTH'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., ERR, Ptr3D=Ptr3D )
    IF ( .NOT. ASSOCIATED(Ptr3D) ) THEN
       CALL HCO_WARNING('Not found: '//TRIM(DgnName),ERR,THISLOC=LOC)
    ELSE
       SO4_ANTH = Ptr3D(:,:,:)
    ENDIF
    Ptr3D => NULL()

    ! convert to kg/box/sec
    DO L = 1, State_Grid%NZ
       SO4an(:,:,L) = 0.0d0
       SO4an(:,:,L) = SO4_ANTH(:,:,L) * AREA(:,:)
    END DO

    ! NOTE: Biofuels are now lumped into anthro,
    ! so set SO4bf to zero (bmy, 10/1/19)
    SO4bf = 0.0_fp

    !=================================================================
    ! Compute SO4 emissions
    !=================================================================

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, NTOP, SO4, TSO4, L, FEMIS, EFRAC, K )      &
    !$OMP PRIVATE( NDISTINIT, NDIST, MDIST, NDISTFINAL, MADDFINAL ) &
    !$OMP PRIVATE( Ndiag, Mdiag)                                    &
    !$OMP PRIVATE( MADDTOTAL, NDIST2, MDIST2, C , ERRORSWITCH)      &
    !$OMP PRIVATE( BOXVOL, TEMP, PRES, pdbug )                      &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !initialize diagnostics
       Ndiag(:) = 0.0D0
       Mdiag(:) = 0.0D0

       ! Top level of boundary layer at (I,J)
       NTOP = CEILING( State_Met%PBL_TOP_L(I,J) )

       ! Zero SO4 array at all levels
       DO L = 1, State_Grid%NZ
          SO4(L) = 0.0
       ENDDO

       ! Compute total anthro SO4 (surface + 100m) plus biofuel SO4
       TSO4 = 0.d0
       TSO4 = SUM( SO4an(I,J,:) ) + SO4bf(I,J)
       IF ( TSO4 <  0d0 ) THEN
          WRITE(*,*) ' Negative Sulfate emis from hemco at IJ=', I, J
       ENDIF
       IF ( TSO4 == 0d0 ) CYCLE

       !=============================================================
       ! First calculate emission distribution vertically within PBL
       !=============================================================
       ! EFRAC(30) = fraction of total emission splitted for each
       !             level until reaching PBL top.
       EFRAC = 0d0

       !==============================================================
       ! Partition the total anthro SO4 emissions thru the entire
       ! boundary layer (if PBL top is higher than level 2)
       !==============================================================
       ! Add option for non-local PBL (Lin, 03/31/09)
       IF (.NOT. LNLPBL) THEN
          !IF ( NTOP > 2 ) THEN

          ! Loop thru boundary layer
          DO L = 1, State_Grid%NZ

             ! Fraction of PBL spanned by grid box (I,J,L) [unitless]
             EFRAC(L)  = State_Met%F_OF_PBL(I,J,L)

          ENDDO

          !ELSE
          !   EFRAC(1) = ( SO4an(I,J,1) + SO4bf(I,J) ) / TSO4
          !   EFRAC(2) = SO4an(I,J,2) / TSO4
          !ENDIF

          IF ( ABS( SUM( EFRAC(:)) - 1.d0 ) > 1.D-5 ) THEN
             PRINT*, '### ERROR in SRCSF30!'
             PRINT*, '### I, J : ', I, J
             print*, 'EFRAC',EFRAC(:)
             PRINT*, '### SUM(EFRAC) : ', SUM( EFRAC(:) )
             PRINT*, '### This should exactly 1.00'
             CALL ERROR_STOP( 'Check SO4 redistribution', &
                              'SRCSF30 (sulfate_mod.F90)' )
          ENDIF

       ELSE
          ! stop the program for now since I don't totally implement
          ! the subgrid coagulation option w/ Lin's new PBL scheme (win, 1/25/10)
          print *,'If the program stops here, that means you are ',    &
                  'running TOMAS simulation with the new PBL scheme ', &
                  'implemented since GEOS-Chem v.8-02-01.',            &
                  '-----> Try not using the non-local PBL option'
          CALL ERROR_STOP( 'Code does not support new PBL scheme', &
                           'SRCSF30 (sulfate_mod.F90)')

       ENDIF  ! .not. LNLPBL

       !=============================================================
       ! Add the size-resolved SO4 emission to tracer array
       ! Having the options to do sub-grid coagulation or simply
       ! emit.
       ! Sub-grid coagulation reduces the number being emitted
       ! and modifies the mass size distribution of existing particle
       ! as well as the size distribution being emitted.
       ! (win, 10/4/07)
       !=============================================================
       IF ( SGCOAG ) THEN

! ewl: TC1 = Spc(:,:,:,id_NK1:id_NK1+IBINS-1)

          !save number and mass before emission
          !DO M = 1, IBINS
          !   N0(:,M) = TC1(id_NK1+M-1)%Conc(I,J,:)
          !ENDDO
          !M0(:,:) = TC2(I,J,:,1:IBINS)

          DO L = 1, State_Grid%NZ
             !only really need to loop L=1,NTOP
             SO4(L) = TSO4 * EFRAC(L) * DTSRCE
             IF ( SO4(L) == 0.d0 ) CYCLE
             DO K = 1, IBINS
                !set number of sulfate particles emitted
                !as emitted mass * fraction in this bin / avg mass per particle
                ! for this bin
                NDISTINIT(K) = SO4(L) * BFRAC(K) / AVGMASS(K)
                !sfarina - sqrt is expensive.
                !NDISTINIT(K) = SO4(L) * BFRAC(K) / ( SQRT( XK(K)*XK(K+1) ) )
                !set existing number of particles
                NDIST(K) = TC1(id_NK1+K-1)%Conc(I,J,L)
                !sfarina - what are the chances aerosol water and ammonium
                ! are properly equilibrated?
                DO C = 1, ICOMP
                   !set existing mass of each component
                   MDIST(K,C) = TC2(I,J,L,K+(C-1)*IBINS)
                   IF( IT_IS_NAN( MDIST(K,C) ) ) THEN
                      PRINT *,'+++++++ Found NaN in SRCSF30  +++++++'
                      PRINT *,'Location (I,J,L):',I,J,L,'Bin',K,'comp',C
                      CALL  ERROR_STOP('SRCSF30 SGCCOAG','sulfate_mod.F90')
                   ENDIF
                ENDDO
                !initialize emitted sulfate number and mass returned
                ! from subgridcoag
                NDISTFINAL(K) = 0.0D0
                MADDFINAL(K)  = 0.0D0
             ENDDO
             !sfarina subgridcoag does its own mnfix. this call might be
             ! unnecessary?
             CALL MNFIX( NDIST, MDIST, ERRORSWITCH )
             IF( ERRORSWITCH ) PRINT *,'SRCSF30: MNFIX found error ', &
                                       'before SUBGRIDCOAG at ',I,J,L
             ERRORSWITCH = .FALSE.

             !!debug
             !DO K = 1, IBINS
             !   ! Overwrite number and mass before emission for diagnostic
             !   ! just in case there was any change by MNFIX (win, 10/27/08)
             !   N0(L,K) = NDIST(K)
             !   M0(L,K) = MDIST(K,SRTSO4)
             !   Avginit(K) = SUM(MDIST(K,:)) / NDIST(K)
             !ENDDO

             BOXVOL  = State_Met%AIRVOL(I,J,L) * 1.e6 !convert from m3 -> cm3
             TEMP    = State_Met%T(I,J,L)
             PRES    = State_Met%PMID(i,j,l)*100.0 ! in Pa

             pdbug=.false.

             CALL SUBGRIDCOAG( NDISTINIT, NDIST, MDIST, BOXVOL,TEMP, &
                               PRES, TSCALE, NDISTFINAL, MADDFINAL,pdbug)
             DO K = 1, IBINS
                !add number from emissions
                NDIST2(K) = NDIST(K) + NDISTFINAL(K)

                !use this number for the diag
                Ndiag(K)  = Ndiag(K) + NDISTFINAL(K)

                !sfarina - An example to illustrate what's happening here:
                !assuming mass doubling
                !avgmass = .01, .02, .04 (AVGMASS(K) = sqrt(XK(K)*XK(K+1)))
                !N0 = 100, 50, 25
                !M0 = 1, 1, 1
                !emitted SO4 particles: 10, 5, 1 (for a total mass of 0.24)
                !but with subgrid coagulation, we lose 2 particles from bin 1
                !onto particles in each of bins 2 and 3, and 1 from bin 2 to
                ! bin 3, giving us a final distribution of
                !Ndistfinal = 6, 4, 1 (for a total mass of
                !             6*.01 + 4*.02 + 1*.04 = 0.18)
                !but that doesn't conserve mass... those particles are now a
                !little heavier than they were before subgrid coag, so you have
                ! to add that additional mass (maddfinal)
                !(6*.01)+ (4*.02 + 2*.01) + (1*0.04 + 1*.02 + 2*.01) = 0.24
                MADDTOTAL   = NDISTFINAL(K) * AVGMASS(K) + MADDFINAL(K)

                !copy mass from all species
                DO C = 1, ICOMP
                   MDIST2(K,C) = MDIST(K,C)
                ENDDO

                !add mass from emissions as explained above
                MDIST2(K,SRTSO4) = MDIST2(K,SRTSO4) + MADDTOTAL

                !save this for the diag
                Mdiag(K)   = Mdiag(K) + MADDTOTAL

                !sanity check
                if(NDISTFINAL(K) < 0d0) then
                   CALL  ERROR_STOP('negative number emis','sulfate_mod.F90')
                endif
                if(MADDTOTAL < 0d0) then
                   CALL  ERROR_STOP('negative mass emis','sulfate_mod.F90')
                endif
             ENDDO

             !debug - avg particle mass after emission but before mnfix
             !DO K = 1, IBINS
             !   Avginner(K) = SUM(MDIST2(K,:)) / NDIST2(K)
             !ENDDO

             ERRORSWITCH = .FALSE.
             CALL MNFIX( NDIST2, MDIST2, ERRORSWITCH )

             IF( ERRORSWITCH ) PRINT *,'SRCSF30: MNFIX found error ', &
                                       'after SUBGRIDCOAG at ',I,J,L

             DO K = 1, IBINS
                TC1(id_NK1+K-1)%Conc(I,J,L) = NDIST2(K)
                DO C=1,ICOMP
                   TC2(I,J,L,K+(C-1)*IBINS) = MDIST2(K,C)
                ENDDO
             ENDDO

             !debug - avg particle mass final
             !DO K = 1, IBINS
             !   Avgfinal(K) = SUM(MDIST2(K,:)) / NDIST2(K)
             !ENDDO
             !
             !DO K = 1, IBINS
             !!sfarina debug
             !if(TC1(id_NK1+K-1)%Conc(I,J,L)-N0(L,K) < 0d0) then
             ! write(*,*) '"Negative NK emis" details:'
             ! write(*,*) 'NTOP       ', NTOP
             ! write(*,*) 'S_SO4:     ', S_SO4
             ! write(*,*) 'TSO4:      ', TSO4
             ! write(*,*) 'EFRAC(L):  ', EFRAC(L)
             ! DO Bi=1,IBINS
             !  write(*,*) 'Bin        ',Bi
             !  write(*,*) 'n0, TC1    ', N0(l,bi),  TC1(id_NK1+Bi-1)%Conc(I,J,L)
             !  write(*,*) 'ndist1,2   ', NDIST(Bi), NDIST2(Bi)
             !  write(*,*) 'ndistfinal ', NDISTFINAL(Bi)
             !  write(*,*) 'MADDFINAL  ', MADDFINAL(Bi)
             !  write(*,*) 'Avginit    ', Avginit(Bi)
             !  write(*,*) 'Avginner   ', Avginner(Bi)
             !  write(*,*) 'Avgfinal   ', Avgfinal(Bi)
             !  write(*,*) 'M0(so4)    ', M0(l,bi)
             !  DO C=1,ICOMP
             !  write(*,*) 'Component  ', C
             !  write(*,*) 'mdist      ', MDIST(Bi, C)
             !  write(*,*) 'mdist2     ', MDIST2(Bi, C)
             !  write(*,*) 'TC2        ', TC2(i,j,l,(C-1)*IBINS+Bi)
             !  END DO !c
             ! END DO !bi
             !end if
             !
             !ENDDO

          ENDDO ! L loop

          !==============================================================
          ! ND59 Diagnostic: Size-resolved primary sulfate emission in
          !                 [kg S/box/timestep] and the corresponding
          !                  number emission [no./box/timestep]
          !==============================================================
#ifdef BPCH_DIAG
          IF ( ND59 > 0 ) THEN
             !print*, 'JACK IN ND59 SULFATE'
             DO K = 1, IBINS
                !if(TC2(I,J,L,K)-M0(L,K) < 0d0)
                !  print *,'Negative SF emis ',TC2(I,J,L,K)-M0(L,K), &
                !     'at',I,J,L,K
                !if(TC1(id_NK1+K-1)%Conc(I,J,L)-N0(L,K) < 0d0) then
                !   print *,'Negative NK emis ',TC1(id_NK1+K-1)%Conc(I,J,L)-N0(L,K), &
                !     'at',I,J,L,K
                !   print *,'tc1, N0 ',TC1(id_NK1+K-1)%Conc(I,J,L),N0(L,K)
                !end if

                !sfarina - I have studied this extensively and determined that
                !negative NK emis as defined here is not an accurate statement.
                !The particle number in a given bin IS reduced by this
                !subroutine, but it is not reduced by emission. Particle
                !number is reduced by mnfix.
                !if the bin is just about to boil over, that added sulfate mass
                !will trigger a big particle shift in mnfix and it will look
                !like a 'negative number emission' event as defined by this
                !inequality

                !sfarina - changing the definition of this diagnostic to ignore
                ! changes to the distribution by mnfix
                !AD59_SULF(I,J,1,K) = AD59_SULF(I,J,1,K) + &
                !                      (TC2(I,J,L,K)-M0(L,K))*S_SO4
                !AD59_NUMB(I,J,1,K) = AD59_NUMB(I,J,1,K) +
                !                      TC1(id_NK1+K-1)%Conc(I,J,L)-N0(L,K) &
                AD59_SULF(I,J,1,K) = AD59_SULF(I,J,1,K) + Mdiag(K)
                AD59_NUMB(I,J,1,K) = AD59_NUMB(I,J,1,K) + Ndiag(K)
             ENDDO
          ENDIF
#endif

       ELSE
          ! Distributing primary emission without sub-grid coagulation
          !=============================================================
          ! Add SO4 emissions to tracer array
          ! For SF: Convert from [kg SO4/box/s] -> [kg SO4/box/timestep]
          ! For NK: Convert from [kg SO4/box/s] -> [No.   /box/timestep]
          !=============================================================
          DO L = 1, State_Grid%NZ
             SO4(L) = TSO4 * EFRAC(L)
             DO K = 1, IBINS
                TC1(id_NK1+K-1)%Conc(I,J,L) = TC1(id_NK1+K-1)%Conc(I,J,L) + &
                     ( SO4(L) * DTSRCE * BFRAC(K) / AVGMASS(K) )
               TC2(I,J,L,K) = TC2(I,J,L,K) + &
                     ( SO4(L) * DTSRCE * BFRAC(K)               )
             ENDDO
          ENDDO

#ifdef BPCH_DIAG
          !==============================================================
          ! ND59 Diagnostic: Size-resolved primary sulfate emission in
          !                 [kg S/box/timestep] and the corresponding
          !                  number emission [no./box/timestep]
          !==============================================================
          IF ( ND59 > 0 ) THEN
             SO4anbf(:,:,1) = SO4an(:,:,1) + SO4bf(:,:)
             SO4anbf(:,:,2) = SO4an(:,:,2)

             DO L = 1, 2
             DO K = 1, IBINS
                AD59_SULF(I,J,L,K) = AD59_SULF(I,J,L,K) + &
                     ( SO4anbf(I,J,L) * BFRAC(K) * S_SO4 * DTSRCE        )
                AD59_NUMB(I,J,L,K) = AD59_NUMB(I,J,L,K) + &
                     ( SO4anbf(I,J,L) * BFRAC(K) / AVGMASS(K) * DTSRCE   )
             ENDDO
             ENDDO

          ENDIF
#endif

       ENDIF !SGCOAG

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    NULLIFY(TC1)

    IF ( prtDebug ) print *,'   ### Finish SRCSF30'

  END SUBROUTINE SRCSF30
#endif
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: grav_settling
!
! !DESCRIPTION: Subroutine GRAV\_SETTLING performs gravitational settling of
!  sulfate and nitrate in coarse sea salt (SO4S and NITS).
!  (bec, rjp, bmy, 4/20/04, 7/20/04, 10/25/05)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GRAV_SETTLING( Input_Opt, State_Chm, State_Diag, State_Grid, &
                            State_Met, N, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE Species_Mod,        ONLY : Species
    USE TIME_MOD,           ONLY : GET_ELAPSED_SEC
    USE TIME_MOD,           ONLY : GET_TS_CHEM
#ifdef TOMAS
#ifdef BPCH_DIAG
    USE CMN_DIAG_MOD
    USE DIAG_MOD,           ONLY : AD44
#endif
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    INTEGER,        INTENT(IN)    :: N           ! Species index
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !REMARKS:
!  N=1 is SO4S; N=2 is NITS
!                                                                             .
!  tdf Include Coarse Mode DUST size bins
!  N=3 is SO4d2; N=4  is NIT_d1
!  N=5 is SO4d3; N=6  is NIT_d2
!  N=7 is SO4d4; N=8  is NIT_d3
!  N=9 is SO4d4; N=10 is NIT_d4
!  tdf Treat these coated DUSTs as DRY for now
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(fp),  PARAMETER   :: C1     =  0.7674e+0_fp
    REAL(fp),  PARAMETER   :: C2     =  3.079e+0_fp
    REAL(fp),  PARAMETER   :: C3     =  2.573e-11_fp
    REAL(fp),  PARAMETER   :: C4     = -1.424e+0_fp
    REAL(fp),  PARAMETER   :: DEN_SS = 2200.0e+0_fp ! [kg/m3] sea-salt density

    ! Parameters for polynomial coefficients to derive seawater
    ! density. From Tang et al. (1997) (bec, jaegle, 5/11/11)
    REAL(fp),  PARAMETER   :: A1   =  7.93e-3_fp
    REAL(fp),  PARAMETER   :: A2   = -4.28e-5_fp
    REAL(fp),  PARAMETER   :: A3   =  2.52e-6_fp
    REAL(fp),  PARAMETER   :: A4   = -2.35e-8_fp
    REAL(f8),  PARAMETER   :: EPSI = 1.0e-4_f8
!
! !LOCAL VARIABLES:
!
    !tdf from dry_settling
    ! P    Pressure in Kpa 1 mb = 100 pa = 0.1 kPa
    ! Dp   Diameter of aerosol [um]
    ! PDp  Pressure * DP
    ! TEMP Temperature (K)
    ! Slip Slip correction factor
    ! Visc Viscosity of air (Pa s)
    ! VTS  Settling velocity of particle (m/s)
    LOGICAL                :: IS_UPTAKE_SPC
    INTEGER                :: I,         J,      L
    INTEGER                :: DryDep_ID, S
    REAL(fp)               :: DTCHEM
    REAL(fp)               :: DELZ,      DELZ1,  REFF
    REAL(fp)               :: P,         DP,     PDP,      TEMP
    REAL(fp)               :: CONST,     SLIP,   VISC,     FAC1
    REAL(fp)               :: FAC2,      FLUX,   AREA_CM2, RHB
    REAL(fp)               :: RUM,       RWET,   RATIO_R
    REAL(fp)               :: TOT1,      TOT2
    REAL(fp)               :: DEN
    REAL(fp)               :: MW_g
    REAL(f8)               :: RHO1,      WTP,    RHO

    ! Arrays
    REAL(fp)               :: SALA_REDGE_um(2)
    REAL(fp)               :: SALC_REDGE_um(2)
    REAL(fp)               :: VTS(State_Grid%NZ)
    REAL(fp)               :: TC0(State_Grid%NZ)

    ! Pointers
    REAL(fp),      POINTER :: TC(:,:,:)
    TYPE(Species), POINTER :: ThisSpc

    !=================================================================
    ! GRAV_SETTLING begins here!
    !=================================================================

    ! Initialize
    RC       =  GC_SUCCESS

    ! Return if tracers are undefined
    IF ( Input_Opt%LDSTUP ) THEN
       IF ( id_SO4d1 < 0 .and. id_NITd1 < 0 ) RETURN ! tdf
    ENDIF

    ! Copy fields from INPUT_OPT to local variables for use below
    SALA_REDGE_um =  Input_Opt%SALA_REDGE_um
    SALC_REDGE_um =  Input_Opt%SALC_REDGE_um

    ! Chemistry timestep [s]
    DTCHEM        =  GET_TS_CHEM()

    ! Look up this species in the species database
    ThisSpc       => State_Chm%SpcData(N)%Info

    ! Point to the species concentration array
    TC            => State_Chm%Species(N)%Conc(:,:,:)

    ! Set a logical to denote that the species is one of the dust
    ! uptake species, i.e. SO4d{1-4}, NITd{1-4} (bmy, 3/17/17)
    IS_UPTAKE_SPC =  ( ( N /= id_SO4s ) .and. ( N /= id_NITs ) )

    ! Drydep species index
    DryDep_Id     =  ThisSpc%DryDepId

    ! Molecular weight [g], aerosol radius [m], and density [kg/m3]
    MW_g          =  ThisSpc%MW_g
    REFF          =  ThisSpc%Radius
    DEN           =  ThisSpc%Density

    ! Sea salt radius [cm]
    ! The Gerber formula for hygroscopic growth uses the radius in
    ! micrometers instead of centimeters. This fix is implemented by using
    ! RUM instead of RCM (jaegle 5/5/11)
    RUM           =  REFF * 1e+6_fp

    ! Exponential factors
    ! replace with radius in microns (jaegle 5/5/11)
    FAC1          =  C1 * ( RUM**C2 )
    FAC2          =  C3 * ( RUM**C4 )

    !$OMP PARALLEL DO                                                       &
    !$OMP DEFAULT( SHARED                                                 ) &
    !$OMP PRIVATE( I,       J,     L,    VTS,  P,        TEMP, RHB,  RWET ) &
    !$OMP PRIVATE( RATIO_R, RHO,   DP,   PDP,  CONST,    SLIP, VISC, TC0  ) &
    !$OMP PRIVATE( DELZ,    DELZ1, TOT1, TOT2, AREA_CM2, FLUX             ) &
    !$OMP PRIVATE( RHO1,    WTP,   S                                      ) &
    !$OMP SCHEDULE( DYNAMIC                                               )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize
       DO L = 1, State_Grid%NZ
          VTS(L) = 0e+0_fp
       ENDDO

       ! Loop over levels
       DO L = 1, State_Grid%NZ

          ! Pressure at center of the level [kPa]
          ! Use moist air pressure for mean free path (ewl, 3/2/15)
          P       = State_Met%PMID(I,J,L) * 0.1e+0_fp

          ! Temperature [K]
          TEMP    = State_Met%T(I,J,L)

          ! Cap RH at 0.99
          RHB    = MIN( 0.99e+0_fp, State_Met%RH(I,J,L) * 1e-2_fp )

          ! Safety check (phs, 5/1/08)
          RHB     = MAX( TINY(RHB), RHB           )

          ! Aerosol growth with relative humidity in radius [m]
          ! (Gerber, 1985)
          ! Several bug fixes to the Gerber formulation: a log10 (instead of
          ! ln) should be used and the dry radius should be expressed in
          ! micrometers (instead of cm) also add more significant digits to
          ! the exponent (jaegle 5/5/11)
          !RWET    = 1d-6*(FAC1/(FAC2-LOG10(RHB))+RUM**3.e+0_fp)**0.33333e+0_fp

          ! Use equation 5 in Lewis and Schwartz (2006) [m] for sea salt
          ! growth (jaegle 5/11/11)
          RWET = REFF * (4.e+0_fp / 3.7e+0_fp) * &
               ( (2.e+0_fp - RHB)/(1.e+0_fp - RHB) )**(1.e+0_fp/3.e+0_fp)

          ! Ratio dry over wet radii at the cubic power
          RATIO_R = ( REFF / RWET )**3.e+0_fp

          ! Density of the wet aerosol (kg/m3)
          !RHO     = RATIO_R * DEN + ( 1.e+0_fp - RATIO_R ) * 1000.e+0_fp

          ! Above density calculation is chemically unsound because it
          ! ignores chemical solvation.
          ! Iteratively solve Tang et al., 1997 equation 5 to calculate
          ! density of wet aerosol (kg/m3)
          ! (bec, jaegle 5/11/11)
          RATIO_R = ( REFF / RWET )
          ! Assume an initial density of 1000 kg/m3
          RHO  = 1000.e+0_f8
          RHO1 = 0.e+0_f8 !initialize (bec, 6/21/10)
          DO WHILE ( ABS( RHO1-RHO ) .gt. EPSI )
             ! First calculate weight percent of aerosol (kg_RH=0.8/kg_wet)
             WTP    = 100.e+0_f8 * DEN/RHO * RATIO_R**3.e+0_f8
             ! Then calculate density of wet aerosol using equation 5
             ! in Tang et al., 1997 [kg/m3]
             RHO1   = ( 0.9971e+0_f8 + (A1 * WTP) &
                      + (A2 * WTP**2.e+0_f8) &
                      + (A3 * WTP**3.e+0_f8) &
                      + (A4 * WTP**4.e+0_f8) ) * 1000.e+0_f8
             ! Now calculate new weight percent using above density
             ! calculation
             WTP    = 100.e+0_f8 * DEN/RHO1 * RATIO_R**3.e+0_f8
             ! Now recalculate new wet density [kg/m3]
             RHO    = ( 0.9971e+0_f8 + (A1 * WTP) &
                      + (A2 * WTP**2.e+0_f8) &
                      + (A3 * WTP**3.e+0_f8) &
                      + (A4 * WTP**4.e+0_f8) ) * 1000.e+0_f8
          ENDDO

          ! Dp = particle diameter [um]
          ! Use dry radius for dust uptake species
          ! Use wet radius for SO4s, NITs (tdf, bmy, 3/17/17)
          IF ( IS_UPTAKE_SPC ) THEN
             DP = 2.e+0_fp * REFF * 1.e+6_fp  ! SO4d*, NITd*
          ELSE
             DP = 2.e+0_fp * RWET * 1.e+6_fp  ! SO4s,  NITs
          ENDIF

          ! PdP = P * dP [hPa * um]
          PDp = P * Dp

          ! Constant
          ! Use dry radius for dust uptake species
          ! Use wet radius for SO4s, NITs (tdf, bmy, 3/17/17)
          IF ( IS_UPTAKE_SPC ) THEN
             CONST = 2.e+0_fp * DEN * REFF**2 * g0 / 9.e+0_fp ! SO4d*, NITd*
          ELSE
             CONST = 2.e+0_fp * RHO * RWET**2 * g0 / 9.e+0_fp ! SO4s,  NITs
          ENDIF

          !===========================================================
          ! NOTE: Slip correction factor calculations following
          ! Seinfeld, pp464 which is thought to be more accurate
          ! but more computation required. (rjp, 1/24/02)
          !
          ! # air molecule number density
          ! num = P * 1d3 * 6.023d23 / (8.314 * Temp)
          !
          ! # gas mean free path
          ! lamda = 1.d6/( 1.41421 * num * 3.141592 * (3.7d-10)**2 )
          !
          ! # Slip correction
          ! Slip = 1. + 2. * lamda * (1.257 + 0.4 * exp( -1.1 * Dp
          !     &     / (2. * lamda))) / Dp
          !
          ! NOTE: Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
          ! which produces slip correction factore with small error
          ! compared to the above with less computation.
          !===========================================================

          ! Slip correction factor (as function of P*dp)
          Slip = 1.e+0_fp+(15.60e+0_fp + 7.0e+0_fp * &
               EXP(-0.059e+0_fp * PDp)) / PDp

          !=====================================================
          ! NOTE, Eq) 3.22 pp 50 in Hinds (Aerosol Technology)
          ! which produce slip correction factor with small
          ! error compared to the above with less computation.
          ! tdf
          !=====================================================

          ! Viscosity [Pa*s] of air as a function of temperature
          VISC = 1.458e-6_fp * (Temp)**(1.5e+0_fp) / &
                 ( Temp + 110.4e+0_fp )

          ! Settling velocity [m/s]
          VTS(L) = CONST * Slip / VISC
       ENDDO

       ! Method is to solve bidiagonal matrix which is
       ! implicit and first order accurate in z (rjp, 1/24/02)

       ! Save initial tracer concentration in column
       DO L = 1, State_Grid%NZ
          TC0(L) = TC(I,J,L)
       ENDDO

       ! We know the boundary condition at the model top
       L    = State_Grid%MaxChemLev
       DELZ = State_Met%BXHEIGHT(I,J,L)

       TC(I,J,L) = TC(I,J,L) / ( 1.e+0_fp + DTCHEM * VTS(L) / DELZ )

       DO L = State_Grid%MaxChemLev-1, 1, -1
          DELZ  = State_Met%BXHEIGHT(I,J,L)
          DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
          TC(I,J,L) = 1.e+0_fp / ( 1.e+0_fp + DTCHEM * VTS(L) / DELZ ) &
                      * ( TC(I,J,L) + DTCHEM * VTS(L+1) / DELZ1 &
                      *  TC(I,J,L+1) )
       ENDDO

       !==============================================================
       ! DIAGNOSTIC: Drydep flux [molec/cm2/s]
       ! (specifically sea salt loss diagnostics)
       !==============================================================
#ifdef TOMAS
#ifdef BPCH_DIAG
       !-----------------------------------------------------------
       ! ND44 DIAGNOSTIC (bpch)
       ! Dry deposition flux loss [molec/cm2/s]
       !
       ! NOTE: Bpch diagnostics are being phased out.
       !-----------------------------------------------------------
       IF ( ND44 > 0 ) THEN

          ! Initialize
          TOT1 = 0e+0_fp
          TOT2 = 0e+0_fp

          ! Compute column totals of TCO(:) and TC(I,J,:,N)
          DO L = 1, State_Grid%NZ
             TOT1 = TOT1 + TC0(L)
             TOT2 = TOT2 + TC(I,J,L)
          ENDDO

          ! Surface area [cm2]
          AREA_CM2 = State_Grid%Area_M2(I,J) * 1e+4_fp

          ! Convert sea salt/dust flux from [kg/s] to [molec/cm2/s]
          FLUX     = ( TOT1 - TOT2 ) / DTCHEM
          FLUX     = FLUX * AVO / ( MW_g * 1.e-3_fp ) / AREA_CM2

          ! Store in global AD44 array for bpch diagnostic output
          AD44(I,J,DryDep_Id,1) = AD44(I,J,DryDep_Id,1) + FLUX

       ENDIF
#endif
#endif

       !-----------------------------------------------------------
       ! HISTORY (aka netCDF diagnostics)
       ! Dry deposition flux loss [molec/cm2/s]
       !
       ! NOTE: Eventually think about converting this
       ! diagnostic to more standard units [kg/m2/s]
       !-----------------------------------------------------------
       IF ( State_Diag%Archive_DryDepChm .OR. &
            State_Diag%Archive_DryDep        ) THEN

          ! Initialize
          TOT1 = 0e+0_fp
          TOT2 = 0e+0_fp

          ! Compute column totals of TCO(:) and TC(I,J,:,N)
          DO L = 1, State_Grid%NZ
             TOT1 = TOT1 + TC0(L)
             TOT2 = TOT2 + TC(I,J,L)
          ENDDO

          ! Surface area [cm2]
          AREA_CM2 = State_Grid%Area_M2(I,J) * 1e+4_fp

          ! Convert sea salt/dust flux from [kg/s] to [molec/cm2/s]
          FLUX     = ( TOT1 - TOT2 ) / DTCHEM
          FLUX     = FLUX * AVO / ( MW_g * 1.e-3_fp ) / AREA_CM2

          ! Drydep flux in chemistry only
          S = State_Diag%Map_DryDepChm%id2slot(DryDep_Id)
          IF ( S > 0 ) THEN
             State_Diag%DryDepChm(I,J,DryDep_Id) = FLUX
          ENDIF
       ENDIF

    ENDDO ! I
    ENDDO ! J
    !$OMP END PARALLEL DO

    ! Free pointers
    ThisSpc => NULL()
    TC      => NULL()

  END SUBROUTINE GRAV_SETTLING
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_dms
!
! !DESCRIPTION: Subroutine CHEM\_DMS is the DMS chemistry subroutine from Mian
!  Chin's GOCART model, modified for use with the GEOS-CHEM model.
!  (rjp, bdf, bmy, 5/31/00, 10/15/09)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_DMS( Input_Opt,  State_Chm, State_Diag, &
                       State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE Input_Opt_Mod,    ONLY : OptInput
    USE Species_Mod,      ONLY : SpcConc
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Diag_Mod,   ONLY : DgnState
    USE State_Grid_Mod,   ONLY : GrdState
    USE State_Met_Mod,    ONLY : MetState
    USE TIME_MOD,         ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Reaction List (by Mian Chin, chin@rondo.gsfc.nasa.gov)
!  ============================================================================
!                                                                             .
!  R1:    DMS + OH  -> a*SO2 + b*MSA                OH addition channel
!         k1 = { 1.7e-42*exp(7810/T)*[O2] / (1+5.5e-31*exp(7460/T)*[O2] }
!         a = 0.75, b = 0.25
!                                                                             .
!  R2:    DMS + OH  ->   SO2 + ...                  OH abstraction channel
!         k2 = 1.2e-11*exp(-260/T)
!                                                                             .
!         DMS_OH = DMS0 * exp(-(r1+r2)* NDT1)
!         where DMS0 is the DMS concentration at the beginning,
!         r1 = k1*[OH], r2 = k2*[OH].
!                                                                             .
!  R3:    DMS + NO3 ->   SO2 + ...
!         k3 = 1.9e-13*exp(500/T)
!                                                                             .
!         DMS = DMS_OH * exp(-r3*NDT1)
!         where r3 = k3*[NO3].
!                                                                             .
!  R4:    DMS + X   ->   SO2 + ...
!         assume to be at the rate of DMS+OH and DMS+NO3 combined.
!                                                                             .
!  The production of SO2 and MSA here, PSO2_DMS and PMSA_DMS, are saved
!  for use in CHEM_SO2 and CHEM_MSA subroutines as a source term.  They
!  are in unit of [v/v/timestep].
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: FX = 1.0e+0_fp
    REAL(fp), PARAMETER :: A  = 0.75e+0_fp
    REAL(fp), PARAMETER :: B  = 0.25e+0_fp

    ! From D4: only 0.8 efficiency, also some goes to DMSO and lost.
    ! So we assume 0.75 efficiency for DMS addtion channel to form
    ! products.
    REAL(fp), PARAMETER :: EFF = 1e+0_fp
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL             :: IS_FULLCHEM
    INTEGER             :: I,    J,    L
    REAL(fp)            :: TK,   O2,   RK1,    RK2,    RK3,   F
    REAL(fp)            :: DMS,  DMS0, DMS_OH, DTCHEM, XOH,   XN3
    REAL(fp)            :: XX,   OH,   OH0,    XNO3,   XNO30, LOH
    REAL(fp)            :: LNO3, BOXVL

    ! Strings
    CHARACTER(LEN=255)  :: ErrMsg, ThisLoc

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    ! Arrays
    REAL(fp)            :: GLOBAL_NO3(State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    !=================================================================
    ! CHEM_DMS begins here!
    !=================================================================
    IF ( id_DMS < 0 ) RETURN

    ! Assume success
    RC          = GC_SUCCESS

    ! Set location for error messages
    ThisLoc  = ' -> at CHEM_DMS (in module GeosCore/sulfate_mod.F90)'

    ! Copy fields from INPUT_OPT to local variables for use below
    IS_FULLCHEM = Input_Opt%ITS_A_FULLCHEM_SIM

    ! Point to chemical species array [v/v dry]
    Spc         => State_Chm%Species

    ! DTCHEM is the chemistry timestep in seconds
    DTCHEM      = GET_TS_CHEM()

    ! Factor to convert AIRDEN from kgair/m3 to molecules/cm3:
    f           = 1000.e+0_fp / AIRMW * AVO * 1.e-6_fp

    ! Evaluate offline global NO3 from HEMCO
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'GLOBAL_NO3', GLOBAL_NO3, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get data for GLOBAL_NO3 from HEMCO!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Do the chemistry over all chemically-active grid boxes!
    !=================================================================
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, TK, O2, DMS0,OH, XNO3, RK1, RK2, BOXVL ) &
    !$OMP PRIVATE( RK3, DMS_OH, DMS, OH0, XNO30, XOH, XN3, XX, LOH, LNO3 ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Skip non-chemistry boxes
       IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

       ! Temperature [K]
       TK     = State_Met%T(I,J,L)

       ! Get O2 [molec/cm3], DMS [v/v], OH [molec/cm3]
       O2     = State_Met%AIRDEN(I,J,L) * f * 0.21e+0_fp
       DMS0   = Spc(id_DMS)%Conc(I,J,L)
       OH     = GET_OH(  I, J, L, Input_Opt, State_Chm, State_Met )

       ! Get NO3 [molec/cm3]
       !==============================================================
       ! Offline simulation: Read monthly mean GEOS-CHEM NO3 fields
       ! If at nighttime, use the monthly mean NO3 concentration from
       ! the NO3 array.  If during the daytime, set the NO3 concentration
       ! to zero.  We don't have to relax to the monthly mean
       ! concentration every 3 hours (as for HNO3) since NO3 has a
       ! very short lifetime. (rjp, bmy, 12/16/02)
       !==============================================================
       IF ( State_Met%SUNCOS(I,J) > 0e+0_fp ) THEN
          ! NO3 goes to zero during the day
          XNO3 = 0e+0_fp
       ELSE
          ! At night: Get global offline NO3 [v/v] and convert to [molec/cm3]
          XNO3 = GLOBAL_NO3(I,J,L) * State_Met%AIRDEN(I,J,L) * 1.0e-3_fp &
                                   * AVO / AIRMW
       ENDIF
       ! Make sure NO3 is not negative
       XNO3  = MAX( XNO3, 0e+0_fp )

       !==============================================================
       ! (1) DMS + OH:  RK1 - addition channel
       !                RK2 - abstraction channel
       !==============================================================
       RK1 = 0.e+0_fp
       RK2 = 0.e+0_fp
       RK3 = 0.e+0_fp

       IF ( OH > 0.e+0_fp ) THEN
          RK1 = ( 1.7e-42_fp * EXP( 7810.e+0_fp / TK ) * O2 ) / &
                ( 1.e+0_fp + 5.5e-31_fp * EXP( 7460.e+0_fp / TK ) * O2 ) * OH

          ! Update reaction rate to match JPL06 and full chem
          ! (jaf, bmy, 10/15/09)
          RK2 = 1.1e-11_fp * EXP( -240.e+0_fp / TK ) * OH
       ENDIF

       !==============================================================
       ! (2) DMS + NO3 (only happens at night):
       !==============================================================
       IF ( State_Met%SUNCOS(I,J) <= 0e+0_fp ) THEN
          RK3 = 1.9e-13_fp * EXP( 500.e+0_fp / TK ) * XNO3
       ENDIF

       !==============================================================
       ! Update DMS concentrations after reaction with OH and NO3,
       ! and also account for DMS + X assuming at a rate as
       ! (DMS+OH)*Fx in the day and (DMS+NO3)*Fx at night:
       !
       ! DMS_OH :  DMS concentration after reaction with OH
       ! DMS    :  DMS concentration after reaction with NO3
       !           (min(DMS) = 1.0E-32)
       !
       ! NOTE: If we are doing a coupled fullchem/aerosol run, then
       ! also modify OH and NO3 concentrations after rxn w/ DMS.
       !==============================================================
       DMS_OH = DMS0   * EXP( -( RK1 + RK2 ) * Fx * DTCHEM )
       DMS    = DMS_OH * EXP( -( RK3       ) * Fx * DTCHEM )
       IF ( DMS < SMALLNUM ) DMS = 0e+0_fp

       ! Archive initial OH and NO3 for diagnostics
       OH0    = OH
       XNO30  = XNO3

       IF ( IS_FULLCHEM ) THEN

          ! Update OH after rxn w/ DMS (coupled runs only)
          OH    = OH0 - ( ( DMS0 - DMS_OH ) * State_Met%AIRDEN(I,J,L) * f )
          IF ( OH < SMALLNUM ) OH = 0e+0_fp

          ! Update NO3 after rxn w/ DMS (coupled runs only)
          XNO3  = XNO30 - ( ( DMS_OH - DMS ) * State_Met%AIRDEN(I,J,L) * f )
          IF ( XNO3 < SMALLNUM ) XNO3 = 0e+0_fp

       ENDIF

       ! Save DMS back to the tracer array
       Spc(id_DMS)%Conc(I,J,L) = DMS

       !==============================================================
       ! Save SO2 and MSA production from DMS oxidation
       ! in [mixing ratio/timestep]:
       !
       ! SO2 is formed in DMS+OH addition (0.85) and abstraction
       ! (1.0) channels as well as DMS + NO3 reaction.  We also
       ! assume that SO2 yield from DMS + X is 1.0.
       !
       ! MSA is formed in DMS + OH addition (0.15) channel.
       !==============================================================
       IF ( ( RK1 + RK2 ) == 0.e+0_fp ) THEN
          PMSA_DMS(I,J,L) = 0.e+0_fp
       ELSE
          PMSA_DMS(I,J,L) = ( DMS0 - DMS_OH ) * &
                              B*RK1 / ( ( RK1 + RK2 ) * Fx ) * EFF
       ENDIF

       PSO2_DMS(I,J,L) =  DMS0 - DMS - PMSA_DMS(I,J,L)

       !==============================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Production and loss diagnostics
       !==============================================================

       ! P(SO2) from DMS+OH, DMS+NO3, and DMS+X
       XOH  = ( DMS0   - DMS_OH - PMSA_DMS(I,J,L) ) / &
                  Fx   * State_Met%AD(I,J,L) / TCVV_S
       XN3  = ( DMS_OH - DMS ) / Fx * State_Met%AD(I,J,L) / TCVV_S
       XX   = ( ( DMS0 - DMS ) * State_Met%AD(I,J,L) / TCVV_S ) &
                - XOH  - XN3

       ! Grid box volume [cm3]
       BOXVL = State_Met%AIRVOL(I,J,L) * 1e+6_fp

       ! Convert L(OH) and L(NO3) from [molec/cm3] to [kg/timestep]
       LOH  = ( OH0   - OH  ) * BOXVL / XNUMOL_OH
       LNO3 = ( XNO30 - XNO3) * BOXVL / XNUMOL_NO3

       ! Store P(SO2) from DMS + OH [kg S/s]
       IF ( State_Diag%Archive_ProdSO2fromDMSandOH ) THEN
          State_Diag%ProdSO2fromDMSandOH(I,J,L) = XOH / DTCHEM
       ENDIF

       ! Store P(SO2) from DMS + NO3 [kg S/s]
       IF ( State_Diag%Archive_ProdSO2fromDMSandNO3 ) THEN
          State_Diag%ProdSO2fromDMSandNO3(I,J,L) = XN3 / DTCHEM
       ENDIF

       ! Store P(SO2) from DMS + NO3 [kg S/s]
       IF ( State_Diag%Archive_ProdSO2fromDMS ) THEN
          State_Diag%ProdSO2fromDMS(I,J,L) = &
               ( PSO2_DMS(I,J,L) * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
       ENDIF

       ! Store P(MSA) from DMS [kg S/s]
       IF ( State_Diag%Archive_ProdMSAfromDMS ) THEN
          State_Diag%ProdMSAfromDMS(I,J,L) = &
               ( PMSA_DMS(I,J,L) * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
       ENDIF

       !==============================================================
       ! For a coupled fullchem/aerosol run, save OH [molec/cm3]
       ! and NO3 [molec/cm3] back into State_Chm%Species
       !==============================================================
       IF ( IS_FULLCHEM ) THEN
          Spc(id_OH  )%Conc(I,J,L) = OH
          Spc(id_NO3 )%Conc(I,J,L) = XNO3
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE CHEM_DMS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_h2o2
!
! !DESCRIPTION: Subroutine CHEM\_H2O2 is the H2O2 chemistry subroutine for
!  offline sulfate simulations.  For coupled runs, H2O2 chemistry is already
!  computed by the SMVGEAR module. (rjp, bmy, 11/26/02, 10/25/05)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_H2O2( Input_Opt, State_Chm, State_Diag, State_Grid, &
                        State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE Input_Opt_Mod,    ONLY : OptInput
    USE Species_Mod,      ONLY : SpcConc
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Diag_Mod,   ONLY : DgnState
    USE State_Grid_Mod,   ONLY : GrdState
    USE State_Met_Mod,    ONLY : MetState
    USE TIME_MOD,         ONLY : GET_MONTH
    USE TIME_MOD,         ONLY : GET_TS_CHEM
    USE TIME_MOD,         ONLY : ITS_A_NEW_MONTH
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOC
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: A = 2.9e-12_fp
!
! !LOCAL VARIABLES:
!
    ! SAVEd Scalars
    LOGICAL            :: FIRST     = .TRUE.
    INTEGER, SAVE      :: LASTMONTH = -99

    ! Scalars
    INTEGER            :: I,     J,    L
    REAL(fp)           :: DT,    Koh,  DH2O2, M,    F ,   XTAU
    REAL(fp)           :: H2O20, H2O2, ALPHA, FREQ, PHOTJ

    ! Strings
    CHARACTER(LEN=255) :: FILENAME, ErrMsg, ThisLoc

    ! Arrays
    REAL(fp)           :: PH2O2m(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(fp)           :: JH2O2(State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    ! Pointers
    TYPE(SpcConc), POINTER  :: Spc(:)

    !=================================================================
    ! CHEM_H2O2 begins here!
    !=================================================================
    IF ( id_H2O2 < 0 ) RETURN

    ! Assume success
    RC        = GC_SUCCESS

    ! Set location for error messages
    ThisLoc  = ' -> at CHEM_H2O2 (in module GeosCore/sulfate_mod.F90)'

    ! Point to chemical species array [v/v dry]
    Spc       => State_Chm%Species

    ! Chemistry timestep [s]
    DT        = GET_TS_CHEM()

    ! Factor to convert AIRDEN from kgair/m3 to molecules/cm3:
    F         = 1000.e+0_fp / AIRMW * AVO * 1.e-6_fp

    ! Evaluate offline fields from HEMCO for P(H2O2)
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'PH2O2', PH2O2m, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get data for PH2O2 from HEMCO!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Evaluate offline fields from HEMCO for J(H2O2)
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'JH2O2', JH2O2, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get data for JH2O2 from HEMCO!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Loop over tropopsheric grid boxes and do chemistry
    !=================================================================
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, M, H2O20, KOH, FREQ, ALPHA, DH2O2, H2O2 ) &
    !$OMP PRIVATE( PHOTJ )  &
    !$OMP SCHEDULE( DYNAMIC )
    DO L  = 1, State_Grid%NZ
    DO J  = 1, State_Grid%NY
    DO I  = 1, State_Grid%NX

       ! Initialize for safety's sake
       FREQ = 0e+0_fp

       ! Skip non-chemistry boxes
       IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

       ! Density of air [molec/cm3]
       M     = State_Met%AIRDEN(I,J,L) * f

       ! Initial H2O2 [v/v]
       H2O20 = Spc(id_H2O2)%Conc(I,J,L)

       ! Loss frequenty due to OH oxidation [s-1]
       KOH   = A * EXP( -160.e+0_fp / State_Met%T(I,J,L) ) * &
               GET_OH( I, J, L, Input_Opt, State_Chm, State_Met )

       ! Now do all dry deposition in mixing_mod.F90 (ckeller, 3/5/15)
       FREQ = 0.e+0_fp

       ! Impose a diurnal variation of jH2O2 by multiplying COS of
       ! solar zenith angle normalized by maximum solar zenith angle
       ! because the archived JH2O2 is for local noon time
       IF ( COSZM(I,J) > 0.e+0_fp ) THEN
          PHOTJ = JH2O2(I,J,L) * State_Met%SUNCOS(I,J) / COSZM(I,J)
          PHOTJ = MAX( PHOTJ, 0e+0_fp )
       ELSE
          PHOTJ = 0e+0_fp
       ENDIF

       ! Compute loss fraction from OH, photolysis, drydep [unitless].
       ALPHA = 1.e+0_fp + ( KOH + PHOTJ + FREQ ) * DT

       ! Delta H2O2 [v/v]
       ! PH2O2m is in kg/m3 (from HEMCO), convert to molec/cm3/s (mps,9/18/14)
       DH2O2 = ( PH2O2m(I,J,L) / TS_EMIS * XNUMOL_H2O2 / CM3PERM3 ) &
               * DT / ( ALPHA * M )

       ! Final H2O2 [v/v]
       H2O2  = ( H2O20 / ALPHA + DH2O2 )
       IF ( H2O2 < SMALLNUM ) H2O2 = 0e+0_fp

       ! Store final H2O2 in Spc
       Spc(id_H2O2)%Conc(I,J,L) = H2O2

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE CHEM_H2O2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_so2
!
! !DESCRIPTION: Subroutine CHEM\_SO2 is the SO2 chemistry subroutine.
!  (rjp, bmy, 11/26/02, 8/26/10)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_SO2( Input_Opt, State_Chm,  State_Diag, State_Grid, &
                       State_Met, FullRun,    RC )
!
! !USES:
!
    USE CMN_SIZE_Mod,         ONLY : NDSTBIN
    USE DUST_MOD,             ONLY : GET_DUST_ALK      ! tdf 04/08/08
    USE ErrCode_Mod
    USE ERROR_MOD,            ONLY : IS_SAFE_EXP
    USE ERROR_MOD,            ONLY : SAFE_DIV
    USE Input_Opt_Mod,        ONLY : OptInput
    USE PRESSURE_MOD,         ONLY : GET_PCENTER
    USE Species_Mod,          ONLY : SpcConc
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE TIME_MOD,             ONLY : GET_TS_CHEM, GET_MONTH
    USE TIME_MOD,             ONLY : ITS_A_NEW_MONTH
    USE HCO_State_GC_Mod,     ONLY : HcoState
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld, HCO_GC_GetDiagn, LoadHcoValEmis
#ifdef APM
    USE APM_DRIV_MOD,         ONLY : PSO4GAS
    USE APM_DRIV_MOD,         ONLY : XO3
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    LOGICAL,        INTENT(IN)    :: FullRun     ! Modify species conc?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC          ! Success or failure?
!
! !REMARKS:
!  Reaction List (by Rokjin Park)
!  ============================================================================
!  (1 ) SO2 production:
!       (a) DMS + OH, DMS + NO3 (saved in CHEM_DMS)
!       (b) HMS -> SO2 + HCHO (aq)
!                                                                             .
!  (2 ) SO2 loss:
!       (a) SO2 + OH  -> SO4
!       (b) SO2       -> drydep
!       (c) SO2 + H2O2 or O3 (aq) -> SO4
!       (d) SO2 + HCHO (aq)-> HMS
!       (d) SO2 + HMS -> 2 SO4
!                                                                             .
!  (3 ) SO2 = SO2_0 * exp(-bt) +  PSO2_DMS/bt * [1-exp(-bt)]
!                                                                             .
!       where b is the sum of the reaction rate of SO2 + OH and the dry
!       deposition rate of SO2, PSO2_DMS is SO2 production from DMS in
!       MixingRatio/timestep.
!                                                                             .
!  If there is cloud in the gridbox (fraction = fc), then the aqueous
!  phase chemistry also takes place in cloud. The amount of SO2 oxidized
!  by H2O2 in cloud is limited by the available H2O2; the rest may be
!  oxidized due to additional chemistry, e.g, reaction with O3 or O2
!  (catalyzed by trace metal).
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(fp),  PARAMETER  :: HPLUS_45  = 3.16227766016837953e-5_fp  !pH = 4.5
    REAL(fp),  PARAMETER  :: HPLUS_50  = 1.0e-5_fp  !pH = 5.0
    REAL(fp),  PARAMETER  :: MINDAT    = 1.e-20_fp
    REAL(fp),  PARAMETER  :: TNA_CONV  = 31.6_fp * 0.359_fp / 23.0_fp
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL               :: DO_SEASALT_CHEM
    LOGICAL               :: IS_OFFLINE
    LOGICAL               :: IS_FULLCHEM
    LOGICAL               :: LDSTUP
    INTEGER               :: I,      J,       L
    INTEGER               :: II,     NSTEP
    INTEGER               :: BULK,   SIZE_RES
    INTEGER               :: IBIN
    REAL(fp)              :: K0,     Ki,      KK,     M,    L1
    REAL(fp)              :: L2,     L3,      Ld,     F,    Fc
    REAL(fp)              :: RK,     RKT,     DTCHEM, DT_T, TK
    REAL(fp)              :: F1,     RK1,     RK2,    RK3,  SO20
    REAL(fp)              :: SO2_cd, H2O20,   O3,     L2S,  L3S
    REAL(fp)              :: LWC,    KaqH2O2, KaqO3,  PATM
    REAL(fp)              :: ALK,    ALK1,    ALK2,    SO2_ss
    REAL(fp)              :: AlkA,   AlkC
    REAL(fp)              :: Kt1,    Kt2
    REAL(fp)              :: PSO4E,  PSO4F,   Kt1N,    Kt2N
    REAL(fp)              :: XX, Kt1L, Kt2L
    REAL(fp)              :: HPLUS,  SO4nss, TNH3,   TNO3,  GNO3, ANIT
    REAL(fp)              :: LSTOT,  ALKdst, ALKss,  ALKds, NH3, CL, TNA
    REAL(fp)              :: SSCvv,  aSO4,   SO2_sr, SR,    TANIT
    REAL(fp)              :: TFA,  TAA,   TDCA    ! (jmm, 12/03/2018)
    REAL(fp)              :: HCHO0,  HMSc,    KaqHCHO, KaqHMS    ! (jmm, 06/07/2018)
    REAL(fp)              :: L7,     L7S,    L7_b, L7S_b, HMS0   ! (jmm, 06/07/2018)
    REAL(fp)              :: L8,     L8S,    OH0, KaqHMS2       ! (jmm, 06/26/2018)
    REAL(fp)              :: PSO4d_tot, PNITd_tot
    REAL(fp)              :: SO2_gas,   PH2SO4d_tot
    REAL(fp)              :: H2SO4_cd,  H2SO4_gas

    ! (qjc, 04/10/16)
    REAL(fp)              :: L5,L5S,SRo3,SRhobr
    REAL(fp)              :: L5_1,L5S_1,L3_1,L3S_1,KaqO3_1
    REAL(fp)              :: HSO3aq, SO3aq
    REAL(fp)              :: SO4H1_vv, SO4H2_vv, LSTOT0
    REAL(fp)              :: SO2_ss0, rSIV, fupdateHOBr_0
    REAL(fp)              :: HCO3, HCHOBr, KO3, KHOBr, f_srhobr, HOBr0
    REAL(fp)              :: TMP, LSTOT_HMS                ! (jmm, 06/15/18)

    REAL(fp)              :: KaqO2, L4, L4S, MnII, FeIII
    REAL(fp)              :: DUST,  Mn_ant,  Mn_nat
    REAL(fp)              :: Mn_tot, Mn_d,    Fe_d
    REAL(fp)              :: Fe_ant, Fe_nat,  Fe_tot
    REAL(fp)              :: Fe_d_ant, Fe_d_nat
    REAL(fp)              :: IONIC


    REAL(fp)              :: L6,L6S,SRhocl,L6_1,L6S_1      !XW
    REAL(fp)              :: SO4H3_vv, SO4H4_vv            !XW
    REAL(fp)              :: fupdateHOCl_0  !XW
    REAL(fp)              :: HCHOCl, KHOCl, f_srhocl, HOCl0 !XW
    REAL(fp)              :: one_m_KRATE

    ! Arrays
    ! tdf 04/07/08
    REAL(fp)              :: ALK_d   (NDSTBIN)
    REAL(fp)              :: ALKA_d  (NDSTBIN)
    REAL(fp)              :: PSO4_d  (NDSTBIN)
    REAL(fp)              :: PNIT_d  (NDSTBIN)
    REAL(fp)              :: PH2SO4_d(NDSTBIN)
    REAL(fp)              :: KTN     (NDSTBIN)
    REAL(fp)              :: KTS     (NDSTBIN)
    REAL(fp)              :: KTH     (NDSTBIN)
    !tdf KTH now contains the fraction of uptake of H2SO4 on to each of the
    ! dust size bins, based on a size- and area-weighted formulism
    ! (GET_DUST_ALK)
    REAL(fp)              :: O3m(State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)
    !REAL(fp), POINTER     :: SSAlk(:,:,:,:)
    REAL(fp), POINTER     :: H2O2s(:,:,:)
    REAL(fp), POINTER     :: SO2s(:,:,:)
    REAL(f4), POINTER     :: Ptr2D(:,:) => NULL()

    ! For HEMCO update
    LOGICAL, SAVE         :: FIRST = .TRUE.

    CHARACTER(LEN=255)    :: ErrMsg, ThisLoc

    !=================================================================
    ! CHEM_SO2 begins here!
    !=================================================================
    IF ( id_H2O2 < 0 .or. id_SO2 < 0  ) RETURN

    ! Assume success
    RC          = GC_SUCCESS
    ErrMsg      = ''
    ThisLoc     = ' -> at CHEM_SO2 (in module GeosCore/sulfate_mod.F90)'

    ! Copy fields from INPUT_OPT to local variables for use below
    IS_FULLCHEM          =  Input_Opt%ITS_A_FULLCHEM_SIM
    IS_OFFLINE           =  Input_Opt%ITS_AN_AEROSOL_SIM
    LDSTUP               =  Input_Opt%LDSTUP
    DTCHEM               =  GET_TS_CHEM()
    Spc                  => State_Chm%Species
    H2O2s                => State_Chm%H2O2AfterChem
    SO2s                 => State_Chm%SO2AfterChem
    State_Chm%isCloud    =  0.0_fp
!    State_Chm%pHcloud    =  0.0_fp
    State_Chm%pHcloud    =  4.5_fp
    State_Chm%QLxpHcloud =  0.0_fp
#ifdef LUO_WETDEP
    State_Chm%pHrain     =  5.6_fp
    State_Chm%QQpHrain   =  0.0_fp
    State_Chm%QQrain     =  0.0_fp
#endif

    ! Set a flag for when to call the SeaSalt_Chem routine
    DO_SEASALT_CHEM =                                                        &
         ( State_Chm%Do_SulfateMod_SeaSalt                           ) .or.  &
         ( .not. FullRun .and. .not. State_Chm%Do_SulfateMod_SeaSalt )

    ! Factor to convert AIRDEN from [kg air/m3] to [molec air/cm3]
    F        = 1000.e+0_fp / AIRMW * AVO * 1.e-6_fp

    ! On first call, get pointers to HEMCO diagnostics arrays.
    ! These are the sea salt aerosol number densities for the fine
    ! and coarse mode, respectively. Values are in # / surface grid
    ! box. These values are needed in the GET_ALK call below.
    ! If the diagnostics are not being found, e.g. because the
    ! sea salt emissions extension is turned off, the passed
    ! pointers NDENS_SALA and NDENS_SALC will stay nullified.
    ! Values of zero will be used in this case! (ckeller, 01/12/2015)
    !IF ( FIRST ) THEN

       ! Sea salt density, fine mode
#if !defined( MODEL_CESM )
      CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, 'SEASALT_DENS_FINE', &
                       StopIfNotFound=.FALSE., RC=RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
        ErrMsg = 'Cannot get HEMCO field SEASALT_DENS_COARSE!'
        CALL GC_Error( ErrMsg, RC, ThisLoc )
        RETURN
      ENDIF

      IF ( ASSOCIATED( Ptr2D ) ) THEN
        ALLOCATE( NDENS_SALA( State_Grid%NX, State_Grid%NY ), STAT=RC )
        NDENS_SALA(:,:) = Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()

      ! Sea salt density, coarse mode
      CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, 'SEASALT_DENS_COARSE', &
                        StopIfNotFound=.FALSE., RC=RC, Ptr2D=Ptr2D )

      ! Trap potential errors
      IF ( RC /= HCO_SUCCESS ) THEN
         ErrMsg = 'Cannot get HEMCO field SEASALT_DENS_COARSE!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      IF ( ASSOCIATED( Ptr2D ) ) THEN
        ALLOCATE( NDENS_SALC( State_Grid%NX, State_Grid%NY ), STAT=RC )
        NDENS_SALC(:,:) = Ptr2D(:,:)
      ENDIF
      Ptr2D => NULL()
#endif

    !IF ( FIRST ) THEN
       ! Adjust first flag
    !   FIRST = .FALSE.
    !ENDIF

    ! If offline aerosol simulation, evaluate fields from HEMCO
    IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'GLOBAL_O3', O3m, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get data for GLOBAL_O3 from HEMCO!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'GLOBAL_HCOOH', GLOBAL_HCOOH, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get data for GLOBAL_HCOOH from HEMCO!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'GLOBAL_ACTA', GLOBAL_ACTA, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Cannot get data for GLOBAL_ACTA from HEMCO!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Load emissions into buffer first for ALK1, ALK2
    ! BEFORE entering loop (hplin, 9/27/20)
    CALL LoadHcoValEmis ( Input_Opt, State_Grid, id_SALA )
    CALL LoadHcoValEmis ( Input_Opt, State_Grid, id_SALC, AltBuffer=.true. )

    ! Loop over chemistry grid boxes
    ! NOTE: Bob Yantosca verified that these !$OMP PRIVATE statements
    ! are correct (12/11/20).  Make sure you add variables to the !$OMP
    ! PRIVATE declaration if they are (1) Scalar variables; (2) Pointers
    ! to other variables; (3) Arrays that have less than (I,J,L) scope.
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I,        J,             L,         SO20,     H2O20      )&
    !$OMP PRIVATE( O3,       PATM,          TK,        K0,       M          )&
    !$OMP PRIVATE( KK,       F1,            RK1,       RK2,      RK         )&
    !$OMP PRIVATE( RKT,      SO2_cd,        L1,        Ld,       L2         )&
    !$OMP PRIVATE( L2S,      L3,            L3S,       FC,       LWC        )&
    !$OMP PRIVATE( KaqH2O2,  KaqO3,         ALK,       ALK1,     ALK2       )&
    !$OMP PRIVATE( Kt1,      Kt2,           SO2_ss,    Kt1N,     Kt2N       )&
    !$OMP PRIVATE( PSO4E,    PSO4F,         XX,        Kt1L,     Kt2L       )&
    !$OMP PRIVATE( TFA,      TAA,           TDCA,      HPLUS,    SO4nss     )&
    !$OMP PRIVATE( TNH3,     TNO3,          CL,        GNO3,     ANIT       )&
    !$OMP PRIVATE( LSTOT,    ALKdst,        ALKds,     ALKss,    NH3        )&
    !$OMP PRIVATE( SSCvv,    aSO4,          SO2_sr,    SR,       TANIT      )&
    !$OMP PRIVATE( BULK,     SIZE_RES,      RC,        AlkA,     AlkC       )&
    !$OMP PRIVATE( ALK_d,    KTS,           KTN,       PSO4_d,   PH2SO4_d   )&
    !$OMP PRIVATE( PNIT_d,   SO2_gas,       KTH,       H2SO4_cd, H2SO4_gas  )&
    !$OMP PRIVATE( Ki,       PH2SO4d_tot,   PSO4d_tot, IBIN,     PNITd_tot  )&
    !$OMP PRIVATE( ALKA_d,   L5,            L5S,       SRo3,     SRhobr     )&
    !$OMP PRIVATE( L3_1,     L3S_1,         KaqO3_1,   L5_1,     L5S_1      )&
    !$OMP PRIVATE( HSO3aq,   SO3aq,         SO4H1_vv,  SO4H2_vv, LSTOT0     )&
    !$OMP PRIVATE( SO2_ss0,  rSIV,          L6S_1,     HCO3,     HCHOBr     )&
    !$OMP PRIVATE( KO3,      KHOBr,         f_srhobr,  HOBr0,    TMP        )&
    !$OMP PRIVATE( L4,       L4S,           DUST,      Mn_ant,   Mn_nat     )&
    !$OMP PRIVATE( Mn_tot,   Fe_ant,        Fe_nat,    Fe_tot,   Fe_d       )&
    !$OMP PRIVATE( Mn_d,     FeIII,         MnII,      Fe_d_ant, Fe_d_nat   )&
    !$OMP PRIVATE( HCHOCl,   KHOCl,         f_srhocl,  HOCl0,    L6         )&
    !$OMP PRIVATE( L6S,      fupdateHOBr_0, SRhocl,    L6_1,     SO4H3_vv   )&
    !$OMP PRIVATE( SO4H4_vv, fupdateHOCl_0, KaqO2,     TNA,      one_m_KRATE)&
    !$OMP PRIVATE( HCHO0,    HMSc,          HMS0,      OH0,      KaqHCHO    )&
    !$OMP PRIVATE( KaqHMS,   KaqHMS2,       L7,        L7S,      L7_b       )&
    !$OMP PRIVATE( L7S_b,    L8,            L8S,       LSTOT_HMS            )&
    !$OMP SCHEDULE( DYNAMIC, 1                                              )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize for safety's sake
       Ld       = 0.0_fp
       LSTOT0   = 0.0_fp
       LSTOT    = 0.0_fp
       LSTOT_HMS = 0.0_fp

       ! Skip non-chemistry boxes
       IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

#ifdef LUO_WETDEP
       ! Save the value of 1.0 - KRATE in a variable,
       ! as this only has to be computed once per grid box
       one_m_KRATE = 1.0_fp - State_Chm%KRATE(I,J,L)
#endif

       ! Initialize [v/v]
       SO20   = Spc(id_SO2)%Conc(I,J,L)
       H2O20  = Spc(id_H2O2)%Conc(I,J,L)

       ! These species are only needed for fullchem simulations
       HOBr0  = 0.0_fp
       HOCl0  = 0.0_fp
       HCHO0  = 0.0_fp
       HMS0   = 0.0_fp
       OH0    = 0.0_fp
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
          HOBr0  = Spc(id_HOBr)%Conc(I,J,L)
          HOCl0  = Spc(id_HOCl)%Conc(I,J,L)
          HCHO0  = Spc(id_CH2O)%Conc(I,J,L)
          HMS0   = Spc(id_HMS)%Conc(I,J,L)
          OH0    = Spc(id_OH)%Conc(I,J,L)
       ENDIF

       ! Calculate O3, defined only in the chemistry grid
       IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
          ! Get O3 from State_Chm%Species%Conc [v/v]
          O3 = 0.0_fp
          IF ( State_Met%InChemGrid(I,J,L) ) THEN
             O3 = State_Chm%Species(id_O3)%Conc(I,J,L)
          ENDIF
       ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
          ! Get offline mean O3 [v/v] for this gridbox and month
          O3 = 0.0_fp
          IF ( L <= State_Grid%MaxChemLev ) THEN
             O3 = O3m(I,J,L)
          ENDIF
       ENDIF

#ifdef APM
       XO3(I,J,L)= O3
#endif

       ! PATM  : Atmospheric pressure in atm
       ! Now use dry air partial pressure (ewl, 4/28/15)
       PATM = State_Met%PMID_DRY( I, J, L ) / ( ATM * 1.e-2_fp )

       ! TK : Temperature [K]
       TK = State_Met%T(I,J,L)

       ! Updated to match JPL 2006 + full chem (jaf, 10/14/09)
       K0 = 3.3e-31_fp * ( 300.e+0_fp / TK )**4.3e+0_fp
       Ki = 1.6e-12_fp

       IF ( IS_OFFLINE ) THEN

          ! Gas phase SO4 production is done here in offline run only
          M   = State_Met%AIRDEN(I,J,L) * F
          KK  = K0 * M / Ki
          F1  = ( 1.e+0_fp + ( LOG10( KK ) )**2 )**( -1 )
          RK1 = ( K0 * M / ( 1.e+0_fp + KK ) ) * 0.6e+0_fp**F1 * &
                  GET_OH( I, J, L, Input_Opt, State_Chm, State_Met )

       ELSE

          ! For online runs, SMVGEAR deals w/ this computation,
          ! so we can simply set RK1 = 0 (rjp, bmy, 3/23/03)
          M   = 0.e+0_fp
          KK  = 0.e+0_fp
          F1  = 0.e+0_fp
          RK1 = 0.e+0_fp

       ENDIF

       ! Now do all dry deposition in mixing_mod.F90 (ckeller, 3/5/15)
       RK2  = 0.e+0_fp

       ! RK: total reaction rate [1/s]
       RK     = ( RK1 + RK2 )

       ! RKT: RK * DTCHEM [unitless] (bmy, 6/1/00)
       RKT    =  RK * DTCHEM

       !==============================================================
       ! Update SO2 conc. after gas phase chemistry and deposition
       !==============================================================
       IF ( RK > 0.e+0_fp ) THEN
          SO2_cd = ( SO20  * EXP( -RKT ) ) + &
                   ( PSO2_DMS(I,J,L) * ( 1.e+0_fp - EXP( -RKT ) ) / RKT )

          L1     = ( SO20 - SO2_cd + PSO2_DMS(I,J,L) ) * RK1/RK

          Ld     = ( SO20 - SO2_cd + PSO2_DMS(I,J,L) ) * RK2/RK

       ELSE
          SO2_cd = SO20
          L1     = 0.e+0_fp
       ENDIF

       ! Isolate H2SO4 for reaction with dust    tdf 3/6/2K9
       IF ( LDSTUP ) THEN
          H2SO4_cd = 0.0_fp

          ! Safety check: only proceed if the Prod diagnostic is archived,
          ! or else this will result in a segmentation fault (bmy, 22 Mar 2022)
          IF ( State_Diag%Archive_Prod .and. id_PSO4 > 0 ) THEN

             ! Compute gas phase SO4 production again, as in offline case
             ! RK1: SO2 + OH(g) [s-1]  (rjp, bmy, 3/23/03)
             M    = State_Met%AIRDEN(I,J,L) * F

             ! Convert State_Diag%Prod from [molec/cm3/s] to [v/v/timestep].
             ! Update by Shixian Zhai added by Bob Yantosca (22 Mar 2022)
             ! See https://github.com/geoschem/geos-chem/discussions/874
             KK       = State_Diag%Prod(I,J,L,id_PSO4)
             H2SO4_cd = KK / M * DTCHEM        
          ENDIF
          
          !tdf Reset these constants to zero to avoid any problems below
          M   = 0.0_fp
          KK  = 0.0_fp
          F1  = 0.0_fp
          RK1 = 0.0_fp
       ENDIF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@ NOTE: The computation of sea salt alkalinity should eventually
!@@@@@ be abstracted out of sulfate_mod.F90 so that the KPP/fullchem/fullchem_*
!@@@@@ modules can make use of it.
!@@@@@   -- Bob Yantosca (09 Sep 2021)
       !==============================================================
       ! Update SO2 conc. after seasalt chemistry (bec, 12/7/04)
       !==============================================================

       ! Get alkalinity of accum (ALK1) and coarse (ALK2) [kg]
       CALL GET_ALK( I,         J,           L,          ALK1,               &
                     ALK2,      Kt1,         Kt2,        Kt1N,               &
                     Kt2N,      Kt1L,        Kt2L,       Input_Opt,          &
                     State_Chm, State_Grid,  State_Met,  RC                 )

       ! Total alkalinity [kg]
       ALK = ALK1 + ALK2

       ! Compute seasalt reaction rates here if:
       ! (1 ) there is alkalinity,
       ! (2 ) there is SO2 present, and
       ! (3 ) O3 is in excess. AND
       ! (4a) we are computing sulfate chemistry in sulfate_mod.F90, or
       ! (4b) we are computing sulfate chemistry in KPP and FULLRUN = F.
       IF ( ( DO_SEASALT_CHEM ) .and. ( ALK    > MINDAT )  .and.             &
            ( SO2_cd > MINDAT ) .and. ( SO2_cd < O3     ) ) THEN


          ! Compute oxidation of SO2 -> SO4 and condensation of
          ! HNO3 -> nitrate within the seasalt aerosol
          CALL SEASALT_CHEM( I,          J,          L,         ALK1,        &
                             ALK2,       SO2_cd,     Kt1,       Kt2,         &
                             Kt1N,       Kt2N,       Kt1L,      Kt2L,        &
                             SO2_ss,     PSO4E,      PSO4F,     AlkA,        &
                             AlkC,       Input_Opt,  State_Met, State_Chm,   &
                             State_Diag, FullRun,    RC                     )

       ELSE

          ! Otherwise set equal to zero
          SO2_ss       = SO2_cd
          PSO4E        = 0.0_fp
          PSO4F        = 0.0_fp
          PNITS(I,J,L) = 0.0_fp
          AlkA         = 0.0_fp
          AlkC         = 0.0_fp
          PNIT(I,J,L)  = 0.0_fp
          PACL(I,J,L)  = 0.0_fp
          PCCL(I,J,L)  = 0.0_fp

       ENDIF

       ! If we are not using KPP to compute seasalt reaction rates,
       ! then update sea salt alkalinity [v/v] in FullRun (XW 12/8/17)
       ! This will make sure that SALAAL and SALCAL have the proper
       ! values befoe sulfate chemistry is computed below.
       IF ( FullRun .and. State_Chm%Do_SulfateMod_SeaSalt ) Then
          Spc(id_SALAAL)%Conc(I,J,L) = AlkA
          Spc(id_SALCAL)%Conc(I,J,L) = AlkC
       ENDIF
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

       IF ( LDSTUP .and. FullRun ) THEN

          !==============================================================
          ! %%% NOTE: THIS IS ONLY DONE FOR ACID UPTAKE SIMULATIONS %%%
          !
          ! Update SO2 conc. after DUST chemistry (tdf, 04/07/08)
          !==============================================================

          ! Get dust alkalinity ALK_d (NDSTBIN) [v/v], Uptake rates for
          ! sulfate, KTS(NDSTBIN), and nitrate, KTN(NDSTBIN) on dust [s-1]
          CALL GET_DUST_ALK( I, J, L, ALK_d, KTS, KTN, KTH, &
                             Input_Opt, State_Met, State_Chm )

          ! Total alkalinity [kg]
          ALK = 0.0e+0_fp

          DO IBIN = 1, NDSTBIN
             ALK = ALK + ALK_d (IBIN)
          END DO

          ! If (1) there is alkalinity, (2) there is SO2 present, and
          ! (3) O3 is in excess, then compute dust SO2 chemistry
          IF  ( ( ALK    > MINDAT )  .AND. &
                ( SO2_cd > MINDAT )  .AND. &
                ( SO2_cd < O3     ) ) THEN

             ! Compute oxidation of SO2 -> SO4 and condensation of
             ! HNO3 -> nitrate within the dust aerosol

             !tdf Call DUST_CHEM using updated SO2_ss after sea salt chemistry
             CALL DUST_CHEM( I,         J,         L,         ALK_d,         &
                             SO2_ss,    H2SO4_cd,  KTS,       KTN,           &
                             KTH,       SO2_gas,   H2SO4_gas, PSO4_d,        &
                             PH2SO4_d,  PNIT_d,    ALKA_d,    Input_Opt,     &
                             State_Met, State_Chm, RC                       )

             ! tdf "SO2_ss" is SO2 mixing ratio remaining after interaction
             ! with dust
             SO2_ss = SO2_gas

             ! tdf "H2SO4_cd" is H2SO4 remaining after interaction with dust
             H2SO4_cd = H2SO4_gas

          ELSE

             ! Otherwise set equal to zero
             SO2_ss  = SO2_ss
             DO IBIN = 1, NDSTBIN
                PSO4_d    (IBIN)   = 0.e+0_fp
                PH2SO4_d  (IBIN)   = 0.e+0_fp
                PNIT_d    (IBIN)   = 0.e+0_fp
             END DO

          ENDIF

       ELSE

          ! Otherwise set equal to zero
          SO2_ss  = SO2_ss
          DO IBIN = 1, NDSTBIN
             PSO4_d    (IBIN)   = 0.e+0_fp
             PH2SO4_d  (IBIN)   = 0.e+0_fp
             PNIT_d    (IBIN)   = 0.e+0_fp
          END DO

       ENDIF    !tdf end of if (LDSTUP) condition

       !==============================================================
       ! Update SO2 concentration after cloud chemistry
       ! SO2 chemical loss rate = SO4 production rate [v/v/timestep]
       !==============================================================

       ! Get cloud fraction from met fields
       FC      = State_Met%CLDF(I,J,L)

       ! Get liquid water content [m3 H2O/m3 air] within cloud from met flds
       ! Units: [kg H2O/kg air] * [kg air/m3 air] * [m3 H2O/1e3 kg H2O]
#ifdef LUO_WETDEP
       ! Luo et al wetdep scheme
       LWC = State_Met%QL(I,J,L)                                             &
           * State_Met%AIRDEN(I,J,L)                                         &
           * 1e-3_fp                                                         &
           + MAX( 0.0_fp, State_Chm%QQ3D(I,J,L) * DTCHEM )

       IF( LWC > 0.d0 ) THEN
          FC = FC                                                            &
             * LWC                                                           &
             / ( LWC + MAX( 0.0_fp, State_Met%QI(I,J,L) )                    &
                     *         State_Met%AIRDEN(I,J,L)                       &
                     *         1e-3_fp                    )
       ELSE
          LWC = 0.d0
       ENDIF
#else
       ! Default scheme
       LWC = State_Met%QL(I,J,L) * State_Met%AIRDEN(I,J,L) * 1e-3_fp
#endif

       LWC  = MAX( 0.0e+0_fp, LWC )

       ! LWC is a grid-box averaged quantity. To improve the representation
       ! of sulfate chemistry, we divide LWC by the cloud fraction and
       ! compute sulfate chemistry based on the LWC within the cloud.  We
       ! get the appropriate grid-box averaged mass of SO2 and sulfate by
       ! multiplying these quantities by FC AFTER computing the aqueous
       ! sulfur chemistry within the cloud. (lzh, jaf, bmy, 5/27/11)
       LWC     = SAFE_DIV( LWC, FC, 0e+0_fp )

       ! Zero variables
       KaqH2O2 = 0.e+0_fp
       KaqO3   = 0.e+0_fp
       KaqO3_1 = 0.e+0_fp !(qjc, 04/10/16)
       L2      = 0.e+0_fp
       L3      = 0.e+0_fp
       L3_1    = 0.e+0_fp !(qjc, 04/10/16)
       L5      = 0.e+0_fp !(qjc, 04/10/16)
       L5_1    = 0.e+0_fp !(qjc, 04/10/16)
       L6      = 0.e+0_fp !XW
       L6_1    = 0.e+0_fp !XW
       L2S     = 0.e+0_fp
       L3S     = 0.e+0_fp
       L3S_1   = 0.e+0_fp !(qjc, 04/10/16)
       L5S     = 0.e+0_fp !(qjc, 04/10/16)
       L5S_1   = 0.e+0_fp !(qjc, 04/10/16)
       KaqO2   = 0.e+0_fp
       L4      = 0.e+0_fp
       L4S     = 0.e+0_fp
       L6S     = 0.e+0_fp !XW
       L6S_1   = 0.e+0_fp !XW
       KaqHCHO = 0.e+0_fp !(jmm, 06/07/18)
       KaqHMS  = 0.e+0_fp !(jmm, 06/07/18)
       KaqHMS2 = 0.e+0_fp !(jmm, 06/26/18)
       L7      = 0.e+0_fp !(jmm, 06/13/18)
       L7_b    = 0.e+0_fp !(jmm, 06/13/18)
       L7S     = 0.e+0_fp !(jmm, 06/13/18)
       L7S_b   = 0.e+0_fp !(jmm, 06/13/18)
       L8      = 0.e+0_fp !(jmm, 06/26/18)
       L8S     = 0.e+0_fp !(jmm, 06/26/18)

       ! If (1) there is cloud, (2) there is SO2 present, (3) T > -15 C, and
       ! (4) liquid water content (LWC) is present (but not small enough to
       ! make divisions blow up), then compute sulfate production in cloud.
       IF ( ( State_Chm%Do_SulfateMod_Cld )                            .and. &
            ( FC     > 1.e-4_fp           )                            .and. &
            ( SO2_ss > MINDAT             )                            .and. &
#ifdef LUO_WETDEP
            ( TK     > 237.0_fp           )                            .and. &
#else
            ( TK     > 258.0_fp           )                            .and. &
#endif
            ( LWC    > 1.0e-20_fp         ) ) THEN

          !===========================================================
          ! NOTE...Sulfate production from aquatic reactions of SO2
          ! with H2O2 & O3 is computed here and followings are
          ! approximations or method used for analytical (integral)
          ! solution of these computations.
          !
          ! 1) with H2O2(aq)
          !      [HSO3-] + [H+] + [H2O2(aq)] => [SO4=]     (rxn)
          !      d[SO4=]/dt = k[H+][HSO3-][H2O2(aq)] (M/s) (rate)
          !
          ! we can rewrite k[H+][HSO3-] as K1 pSO2 hSO2,
          ! where pSO2 is equilibrium vapor pressure of SO2(g)
          ! in atm, and hSO2 is henry's law constant for SO2
          !
          ! Therefore, rate can be written as
          !
          !       k * K1 * pSO2 * hSO2 * pH2O2 * hH2O2,
          !
          ! where pH2O2 is equilibrium vapor pressure of H2O2(g),
          ! and hH2O2 is henry's law constant for H2O2. Detailed
          ! values are given in AQCHEM_SO2 routine.
          !
          ! Let us define a fraction of gas phase of A species
          ! in equilibrium with aqueous phase as
          !
          !        xA  = 1/(1+f),
          !
          ! where  f   = hA * R * T * LWC,
          !        hA  = Henry's constant,
          !        R   = gas constant,
          !        T   = temperature in kelvin,
          !        LWC = liquid water content [m3/m3]
          !
          ! As a result, the rate would become:
          !
          !    d[SO4=]
          !    ------- = k K1 hSO2 hH2O2 xSO2 xH2O2 P P [SO2][H2O2]
          !      dt
          !      ^       ^                            ^   ^    ^
          !      |       |____________________________|   |    |
          !
          !   mole/l/s               mole/l/s            v/v  v/v
          !
          !
          ! And we multiply rate by (LWC * R * T / P) in order to
          ! convert unit from mole/l/s to v/v/s
          !
          ! Finally we come to
          !
          !    d[SO4=]
          !    ------- = KaqH2O2 [SO2][H2O2],
          !      dt
          !
          ! where
          !
          !   KaqH2O2 = k K1 hSO2 hH2O2 xSO2 xH2O2 P LWC R T,
          !
          ! this new rate corresponds to a typical second order
          ! reaction of which analytical (integral) solution is
          !
          !   X  = A0 B0 ( exp[(A0-B0) Ka t] - 1 )
          !      / ( A0 exp[(A0-B0) Ka t] - B0 )
          !
          ! inserting variables into solution then we get
          ! [SO4=] =  [SO2][H2O2](exp[([SO2]-[H2O2]) KaqH2O2 t] - 1 )
          !        / ( [SO2] exp[([SO2]-[H2O2]) KaqH2O2 t] - [H2O2] )
          !
          ! Note...Exactly same method can be applied to O3 reaction
          ! in aqueous phase with different rate constants.
          !===========================================================

	  ! Get concentrations for cloud pH calculation (bec, 12/23/11)

	  ! Get sulfate concentration and convert from [v/v] to
          ! [moles/liter]
	  ! Use a cloud scavenging ratio of 0.7
#ifdef LUO_WETDEP
          SO4nss = Spc(id_SO4)%Conc(I,J,L)    * State_Met%AIRDEN(I,J,L) *         &
                   one_m_KRATE          / ( AIRMW * LWC )
#else
          SO4nss = ( Spc(id_SO4)%Conc(I,J,L)  * State_Met%AIRDEN(I,J,L) *         &
                     0.7_fp             / ( AIRMW * LWC )           )       &
                 + ( Spc(id_SO4s)%Conc(I,J,L) * State_Met%AIRDEN(I,J,L) /         &
                     ( AIRMW * LWC )                                )
#endif

          ! Get HMS cloud concentration and convert from [v/v] to
          ! [moles/liter] (jmm, 06/13/2018)
          ! Use a cloud scavenging ratio of 0.7
          ! assume nonvolatile like sulfate for realistic cloud pH
          ! NOTE: Only needed for fullchem sims, otherwise it's zero
          HMSc = 0.0_fp
          IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
             HMSc  =  Spc(id_HMS)%Conc(I,J,L) * State_Met%AIRDEN(I,J,L) * &
                                0.7e+0_fp / ( AIRMW * LWC )
          ENDIF

          ! Get total ammonia (NH3 + NH4+) concentration [v/v]
          ! Use a cloud scavenging ratio of 0.7 for NH4+
#ifdef LUO_WETDEP
          TNH3 = ( Spc(id_NH4)%Conc(I,J,L) * one_m_KRATE ) + Spc(id_NH3)%Conc(I,J,L)
#else
          TNH3 = ( Spc(id_NH4)%Conc(I,J,L) * 0.7_fp      ) + Spc(id_NH3)%Conc(I,J,L)
#endif

          ! Get total chloride (SALACL + HCL) concentration [v/v]
          ! Use a cloud scavenging ratio of 0.7
#ifdef LUO_WETDEP
          CL = ( Spc(id_SALACL)%Conc(I,J,L) * one_m_KRATE ) + Spc(id_SALCCL)%Conc(I,J,L)
#else
          CL = ( Spc(id_SALACL)%Conc(I,J,L) * 0.7_fp )      + Spc(id_SALCCL)%Conc(I,J,L)
#endif
          IF ( id_HCl > 0 ) THEN
             CL = CL + Spc(id_HCL)%Conc(I,J,L)
          ELSE
             CL = CL + GLOBAL_HCL(I,J,L)
          ENDIF

          ! Get total formic acid concentration [v/v]
          ! jmm (12/3/18)
          ! no cloud scavenging because gases?
          IF ( id_HCOOH > 0 ) THEN
             TFA = Spc(id_HCOOH)%Conc(I,J,L)
          ELSE
             TFA = GLOBAL_HCOOH(I,J,L)
          ENDIF

          ! Get total acetic acid concentration [v/v]
          ! jmm (12/3/18)
          ! no cloud scavenging b/c gases?
          IF ( id_ACTA > 0 ) THEN
             TAA = Spc(id_ACTA)%Conc(I,J,L)
          ELSE
             TAA = GLOBAL_ACTA(I,J,L)
          ENDIF

          ! Get total sea salt NVC concentration expressed as NA+ equivalents
          ! and convert from [v/v] to [moles/liter]
          ! NVC is calculated to balance initial Cl- + alkalinity in
          ! seas salt. Note that we should not consider SO4ss here.
          ! Use a cloud scavenging ratio of 0.7 for fine aerosols
#ifdef LUO_WETDEP
          TNA = ( Spc(id_SALA)%Conc(I,J,L) * State_Met%AIRDEN(I,J,L)   *           &
                  TNA_CONV           * one_m_KRATE               /           &
                  ( AIRMW * LWC )                                  )         &
              + ( Spc(id_SALC)%Conc(I,J,L) * State_Met%AIRDEN(I,J,L)   *           &
                  TNA_CONV           / ( AIRMW * LWC )             )
#else
          TNA = ( Spc(id_SALA)%Conc(I,J,L) * State_Met%AIRDEN(I,J,L)   *           &
                  TNA_CONV           * 0.7_fp                    /           &
                  ( AIRMW * LWC )                                  )         &
              + ( Spc(id_SALC)%Conc(I,J,L) * State_Met%AIRDEN(I,J,L)   *           &
                  TNA_CONV           / ( AIRMW * LWC )             )
#endif
          ! Get total dust cation concentration [mol/L]
          ! Use a cloud scavenging ratio of 1 for dust
          ! to be consistent for how it was calculated for
          ! metal catalyzed SO2 oxidation
          ! Use asumption of dust being 3% soluble Ca2+ and
          ! 0.6% soluble Mg2+ by mass (Fairlie et al., 2010)
          !
          ! Dust treated at non-volatile cation and charge applied in
          ! pH calculation
          !
          ! Move dust calculation from SO2 Metal catalzyed oxidation
          ! up here becasue needed for cloud pH
          ! jmm (12/3/18)
          !
          ! Get dust concentrations [v/v -> ng/m3]
#ifdef LUO_WETDEP
          DUST = ( Spc(id_DST1)%Conc(I,J,L)*one_m_KRATE + Spc(id_DST2)%Conc(I,J,L) +     &
                   Spc(id_DST3)%Conc(I,J,L)             + Spc(id_DST4)%Conc(I,J,L)   )   &
               * 1.e+12_fp * State_Met%AD(I,J,L)                             &
               / ( AIRMW   / State_Chm%SpcData(id_DST1)%Info%MW_g        )   &
               / State_Met%AIRVOL(I,J,L)
#else
          DUST = ( Spc(id_DST1)%Conc(I,J,L)*0.7_fp + Spc(id_DST2)%Conc(I,J,L) +          &
                   Spc(id_DST3)%Conc(I,J,L)        + Spc(id_DST4)%Conc(I,J,L)   )        &
               * 1.e+12_fp * State_Met%AD(I,J,L)                             &
               / ( AIRMW   / State_Chm%SpcData(id_DST1)%Info%MW_g )          &
               / State_Met%AIRVOL(I,J,L)
#endif

          ! Conversion from dust mass to Ca2+ and Mg2+ mol:
          !     0.071*(1/40.08)+0.011*(1/24.31) = 2.22e-3
          !     (Engelbrecht et al., 2016)
          !     1e-12_fp from m3->L & ng->g
          TDCA     = DUST * 2.22e-15_fp / LWC

          ! Get total nitrate (HNO3 + NIT) concentrations [v/v]
          ! Use a cloud scavenging ratio of 0.7 for NIT
          IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
#ifdef LUO_WETDEP
             TNO3 = ( Spc(id_HNO3)%Conc(I,J,L)                                     &
                  +   Spc(id_NIT )%Conc(I,J,L)                                     &
                  +   Spc(id_NITs)%Conc(I,J,L) ) * one_m_KRATE
#else
             TNO3 = Spc(id_HNO3)%Conc(I,J,L)                                       &
                  + ( Spc(id_NIT)%Conc(I,J,L) * 0.7_fp )                           &
                  + Spc(id_NITs)%Conc(I,J,L)
#endif
             GNO3 = Spc(id_HNO3)%Conc(I,J,L) !For Fahey & Pandis decision algorithm
          ELSE IF ( IS_OFFLINE ) THEN
             TANIT = Spc(id_NIT)%Conc(I,J,L) !aerosol nitrate [v/v]
             GNO3  = GLOBAL_HNO3(I,J,L) - TANIT ! gas-phase nitric acid [v/v]
             ANIT  = TANIT * 0.7e+0_fp ! aerosol nitrate in the cloud drops [v/v]
             TNO3  = GNO3 + ANIT   ! total nitrate for cloud pH calculations
          ENDIF

          ! Calculate cloud pH
          CALL GET_HPLUS( SO4nss, HMSc, TNH3, TNO3, SO2_ss, CL,  TNA,        &
                          TDCA,   TFA,  TAA,  TK,   PATM,   LWC, HPLUS_45,   &
                          HPLUS                                             )

          ! Store the cloud pH quantities
          State_Chm%isCloud(I,J,L)    =  1.0_fp
          State_Chm%pHCloud(I,J,L)    = -1.0_fp * log10(HPLUS)
          State_Chm%QLxpHCloud(I,J,L) = State_Chm%pHCloud(I,J,L)             &
                                      * State_Met%QL(I,J,L)

          IF ( Input_Opt%LMETALCATSO2 ) THEN

             !--------------------------------------------------------
             ! Metal catalyzed oxidation of SO2 pathway
             !--------------------------------------------------------

             ! Get dust concentrations [v/v -> ng/m3]
#ifdef TOMAS
             ! TOMAS uses its own dust tracers and does not
             ! carry DST1-4.  Set DUST to zero here. (mps, 2/2/18)
             DUST = 0e+0_fp
#else
             DUST = ( Spc(id_DST1)%Conc(I,J,L)*0.7 + Spc(id_DST2)%Conc(I,J,L) + &
                      Spc(id_DST3)%Conc(I,J,L) + Spc(id_DST4)%Conc(I,J,L) ) * &
                      1.e+12_fp * State_Met%AD(I,J,L) &
                      / ( AIRMW / State_Chm%SpcData(id_DST1)%Info%MW_g ) &
                      / State_Met%AIRVOL(I,J,L)
#endif

             ! Calculate Fe and Mn natural [ng m-3]
             ! Assume that Fe is 3.5% of total dust mass based on
             ! Taylor and McLennan [1985]
             Fe_nat = DUST * 35e-3_fp
             ! and Mn is 50 times less than Fe based on Desbouefs et al.[2005]
             Mn_nat = Fe_nat / 50e+0_fp

             ! Anthropogenic Fe concentrations [v/v -> ng/m3]
             IF ( id_pFe > 0 ) THEN
                Fe_ant = Spc(id_pFe)%Conc(I,J,L) * &
                         1.e+12_fp * State_Met%AD(I,J,L) &
                         / ( AIRMW / State_Chm%SpcData(id_pFe)%Info%MW_g ) &
                         / State_Met%AIRVOL(I,J,L)
             ELSE
                Fe_ant = 0e+0_fp
             ENDIF

             ! Calculate Mn anthropogenic [ng m-3]
             ! assume anthropogenic Mn is 1/30 times anthropogenic Fe
             Mn_ant = Fe_ant / 10e+0_fp

             ! Calculate total Mn and Fe [ng m-3]
             Mn_tot = Mn_ant + Mn_nat
             Fe_tot = Fe_ant + Fe_nat

             ! Convert Mn and Fe [ng m-3] to [mole l-1]

             ! Assume that 50% of Mn is dissolved [Spokes et al., 1994]
             ! Hardcoded MW for Mn
             IF ( LWC > 0e+0_fp ) THEN
                ! Units: ng/m3 * (g/ng) / (g/mol) / (m3 H2O / m3 air) * (m3/L)
                Mn_d = Mn_tot * 1e-9_fp / 54.94e+0_fp / LWC * 1e-3_fp
                Mn_d = Mn_d * 0.5e+0_fp
             ELSE
                Mn_d = 0e+0_fp
             ENDIF

             ! Solubility of Fe is 10% for anthropogenic, and 1% for dust
             IF ( LWC > 0e+0_fp ) THEN
                Fe_d_ant = Fe_ant * 1e-9_fp / &
                           State_Chm%SpcData(id_pFe)%Info%MW_g / &
                           LWC * 1e-3_fp
                Fe_d_nat = Fe_nat * 1e-9_fp / &
                           State_Chm%SpcData(id_pFe)%Info%MW_g / &
                           LWC * 1e-3_fp
                Fe_d     = Fe_d_ant * 0.1e+0_fp + &
                           Fe_d_nat * 0.01e+0_fp
             ELSE
                Fe_d     = 0e+0_fp
             ENDIF

             ! Impose a dependence of Fe speciation on sunlight
             IF ( Input_Opt%ITS_A_FULLCHEM_SIM .or. &
                  Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
                IF ( State_Met%SUNCOS(I,J) > 0e+0_fp ) THEN
                   ! Assume 10% of dissolved Fe is in Fe(III)
                   !oxidation state during the daytime
                   FeIII = Fe_d * 0.1e+0_fp
                ELSE
                   ! Assume 90% of dissolved Fe is in Fe(III)
                   ! oxidation state during the nighttime
                   FeIII = Fe_d * 0.9e+0_fp
                ENDIF

             ENDIF

             ! Assume that dissolved Mn is in Mn(II) oxidation state all of
             ! the time
             MnII = Mn_d

          ELSE

             ! Set Fe and Mn concentrations to zero
             FeIII = 0.0e+0_fp
             MnII  = 0.0e+0_fp

          ENDIF

          IONIC = 0.5e+0_fp * ( 4.0e+0_fp * SO4nss +  &
               TNH3*State_Met%AIRDEN(I,J,L)/(28.97e+0_fp*LWC) + &
               TNO3*State_Met%AIRDEN(I,J,L)/(28.97e+0_fp*LWC) )

          ! Compute aqueous rxn rates for SO2
          CALL AQCHEM_SO2( LWC,     TK,    PATM,    SO2_ss, H2O20, &
                           O3,      HPLUS, MnII,    FeIII, IONIC, &
                           KaqH2O2, KaqO3, KaqO3_1, KaqO2, &
                           HSO3aq,  SO3aq, HCHO0,   KaqHCHO, &
                           KaqHMS, KaqHMS2 )

          !----------------------------------------------------------
          ! Compute loss by H2O2.  Prevent floating-point exception
          ! by not allowing the exponential to go to infinity if
          ! the argument is too large.  (win, bmy, 1/4/09)
          !----------------------------------------------------------

          ! Argument of the exponential
          XX  = ( SO2_ss - H2O20 ) * KaqH2O2 * DTCHEM

          ! Test if EXP(XX) can be computed w/o numerical exception
          IF ( IS_SAFE_EXP( XX ) .and. ABS( XX ) > 0e+0_fp ) THEN

             ! Aqueous phase SO2 loss rate w/ H2O2 [v/v/timestep]
             L2  = EXP( XX )

             ! Loss by H2O2
             L2S = SO2_ss * H2O20 * ( L2 - 1.e+0_fp ) / &
                   ( (SO2_ss * L2) - H2O20 )
          ELSE

             ! NOTE from Jintai Lin (4/28/10):
             ! However, in the case of a negative XX, L2S should be
             ! approximated as SO2_ss, instead of H2O20. In other words,
             ! L2S = SO2_ss * H2O20 * ( L2 - 1.D0 ) / ( (SO2_ss*L2) - H2O20 )
             ! reaches different limits when XX reaches positive infinity
             ! and negative infinity.
             IF ( XX > 0.e+0_fp ) THEN
                L2S = H2O20
             ELSE IF ( XX < 0.e+0_fp) THEN
                L2S = SO2_ss
             ELSE
                !(qjc, 04/10/16) different solution when SO2_ss = H2O20
                L2S = SO2_ss - 1/(KaqH2O2*DTCHEM+1/SO2_ss)
             ENDIF

          ENDIF

          !----------------------------------------------------------
          ! Compute loss by O3.  Prevent floating-point exception
          ! by not allowing the exponential to go to infinity if
          ! the argument is too large. (win, bmy, 1/4/09)
          !----------------------------------------------------------

          ! Argument of the exponential
          XX = ( SO2_ss - O3 ) * KaqO3 * DTCHEM

          ! Test if EXP(XX) can be computed w/o numerical exception
          IF ( IS_SAFE_EXP( XX ) .and. ABS( XX ) > 0e+0_fp ) THEN

             ! Aqueous phase SO2 loss rate w/ O3 [v/v/timestep]
             L3  = EXP( XX )

             ! Loss by O3
             L3S = SO2_ss * O3 * (L3 - 1.e+0_fp)/((SO2_ss * L3) - O3)

          ELSE

             ! Follow the same logic for L3S as described in
             ! Jintai Lin's note above (bmy, 4/28/10)
             IF ( XX > 0.e+0_fp ) THEN
                L3S = O3
             ELSE IF ( XX < 0.e+0_fp) THEN
                L3S = SO2_ss
             ELSE
                !(qjc, 04/10/16) different solution when SO2_ss = O3
                L3S = SO2_ss - 1/(KaqO3*DTCHEM+1/SO2_ss)
             ENDIF
          ENDIF

          !(qjc, 04/10/16)
          !----------------------------------------------------------
          ! Compute loss by O3, but SO3-- only.  Prevent floating-point
          ! exception by not allowing the exponential to go to infinity
          ! if the argument is too large.
          !----------------------------------------------------------

          ! Argument of the exponential
          XX = ( SO2_ss - O3 ) * KaqO3_1 * DTCHEM

          ! Test if EXP(XX) can be computed w/o numerical exception
          IF ( IS_SAFE_EXP( XX ) .and. ABS( XX ) > 0e+0_fp ) THEN

             ! Aqueous phase SO2 loss rate w/ O3 [v/v/timestep]
             L3_1  = EXP( XX )

             ! Loss by O3
             L3S_1 = SO2_ss * O3*(L3_1 - 1.e+0_fp)/((SO2_ss*L3_1)-O3)

          ELSE

             ! Follow the same logic for L3S_1 as described in
             ! Jintai Lin's note above
             IF ( XX > 0.e+0_fp ) THEN
                L3S_1 = O3
             ELSE IF ( XX < 0.e+0_fp) THEN
                L3S_1 = SO2_ss
             ELSE
                !(qjc, 04/10/16) different solution when SO2_ss = O3
                L3S_1 = SO2_ss - 1/(KaqO3_1*DTCHEM+1/SO2_ss)
             ENDIF
          ENDIF

          IF ( Input_Opt%LMETALCATSO2 ) THEN

             !--------------------------------------------------------
             ! Metal catalyzed oxidation of SO2 pathway
             !
             ! Compute loss by O2.  I did not do what Jintai Lin did
             ! above for the other aqueous-phase reactions because it
             ! doesn't make sense for this reaction (bec, 7/7/15)
             !----------------------------------------------------------

             ! Argument of the exponential
             XX = -KaqO2 * DTCHEM

             IF ( IS_SAFE_EXP( XX ) ) THEN
                ! Aqueous phase SO2 loss rate w/ O2 [v/v/timestep]
                L4  = EXP( XX )

                ! Loss by O2
                L4S = SO2_ss * (1.e+0_fp - L4)
             ELSE
                L4S = SO2_ss
             ENDIF

          ELSE

             ! Set loss by O2 to zero
             L4S = 0.0e+0_fp

          ENDIF

          !----------------------------------------------------------
          ! Compute loss by HOBr.  Prevent floating-point exception
          ! by not allowing the exponential to go to infinity if
          ! the argument is too large. !qjc (04/05/16)
          ! Add loss by HOCl, XW (11/08/18)
          !----------------------------------------------------------

          ! Get SO4H (sulfate produced via HOBr) from the Spc array [v/v dry]
          SO4H1_vv = 0.e+0_fp
          SO4H2_vv = 0.e+0_fp
          SO4H3_vv = 0.e+0_fp
          SO4H4_vv = 0.e+0_fp
          IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
             SO4H1_vv = Spc(id_SO4H1)%Conc(I,J,L)
             SO4H2_vv = Spc(id_SO4H2)%Conc(I,J,L)
             SO4H3_vv = Spc(id_SO4H3)%Conc(I,J,L)
             SO4H4_vv = Spc(id_SO4H4)%Conc(I,J,L)
          ENDIF

          L5S   = SO4H1_vv + SO4H2_vv
          L5S_1 = SO4H2_vv
          L6S   = SO4H3_vv + SO4H4_vv
          L6S_1 = SO4H3_vv

          ! make sure sulfate produced is less than SO2 available
          ! (qjc, 06/20/16)
          IF (L5S > SO2_ss) THEN
             fupdateHOBr_0 = SO2_ss/L5S
             L5S   = SO2_ss
             L5S_1 = SO2_ss * L5S_1/L5S
          ELSE
             L5S = L5S
             L5S_1 = L5S_1
             fupdateHOBr_0 = 1.e+0_fp
          ENDIF

          IF (L6S > SO2_ss) THEN
             fupdateHOCl_0 = SO2_ss/L6S
             L6S   = SO2_ss
             L6S_1 = SO2_ss * L6S_1/L6S
          ELSE
             L6S = L6S
             L6S_1 = L6S_1
             fupdateHOCl_0 = 1.e+0_fp
          ENDIF

          !----------------------------------------------------------
          ! Compute loss by HCHO.  Prevent floating-point exception
          ! by not allowing the exponential to go to infinity if
          ! the argument is too large.  (jmm, 06/13/18)
          !----------------------------------------------------------
          IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
             ! Argument of the exponential
             XX  = ( SO2_ss - HCHO0 ) * KaqHCHO * DTCHEM

             ! Test if EXP(XX) can be computed w/o numerical exception
             IF ( IS_SAFE_EXP( XX ) .and. ABS( XX ) > 0e+0_fp ) THEN

                ! Aqueous phase SO2 loss rate w/ HCHO [v/v/timestep]
                L7  = EXP( XX )

                ! Loss by HCHO
                L7S = SO2_ss * HCHO0 * ( L7 - 1.e+0_fp ) / &
                     ( (SO2_ss * L7) - HCHO0 )
             ELSE

                ! NOTE from Jintai Lin (4/28/10):
                ! However, in the case of a negative XX, L7S should be
                ! approximated as SO2_ss, instead of HCHO0. In other words,
                ! L7S = SO2_ss * HCHO0 * ( L7 - 1.D0 ) / ( (SO2_ss*L7) - HCHO0 )
                ! reaches different limits when XX reaches positive infinity
                ! and negative infinity.
                IF ( XX > 0.e+0_fp ) THEN
                   L7S = HCHO0
                ELSE IF ( XX < 0.e+0_fp) THEN
                   L7S = SO2_ss
                ELSE
                   !(qjc, 04/10/16) different solution when SO2_ss = H2O20
                   L7S = SO2_ss - 1/(KaqHCHO*DTCHEM+1/SO2_ss)
                ENDIF

             ENDIF

             IF ( L7S > 1.0e+15_fp .or. L7S < 0 ) THEN
                RC = GC_FAILURE
                PRINT *,'Loc:',I,J,L
                PRINT *,'L7S:',L7S
                PRINT *,'HCHO :',HCHO0
                PRINT *,'SO2_ss',SO2_ss
                PRINT *,'in exp:',XX
                PRINT *,'L7:',L7
                PRINT *,'KaqHCHO',KaqHCHO
             ENDIF


             !----------------------------------------------------------
             ! Compute aqueous SO2 from HMS decomposition. Prevent floating-point
             ! exception by not allowing the exponential to go to infinity if
             ! the argument is too large. This reaction also produces HCHO.
             ! (jmm, 06/13/18)
             !
             ! This is treated as a pseudo first order reaction as it depends
             ! on [HMS] and [OH-]. The [OH-] is incorporated into the rate
             ! constant (KaqHMS) with HPLUS determined seperately. Updates to
             ! [OH-] are incorporated in the next call of GET_HPLUS, which
             ! includes [HMS] in the pH calculation
             !
             ! Also note that [HMS] is in mol/L, and the conversion for the
             ! rate from [mol/L/s] to [v/v/s] is done here
             !
             ! Followed proceedure as done with oxidation by O2 (metal catalyzed)
             !----------------------------------------------------------

             ! Argument of the exponential
             XX  = -KaqHMS * DTCHEM

             ! Test if EXP(XX) can be computed w/o numerical exception
             IF ( IS_SAFE_EXP( XX ) ) THEN

                ! Aqueous phase SO2 production rate w/ HMS [v/v/timestep]
                L7_b  = EXP( XX )

                ! Production by HMS
                L7S_b =  HMSc * ( 1.e+0_fp -  L7_b )
             ELSE
                L7S_b = HMSc
             ENDIF

             ! Convert HMS produced from [mol/L] to [v/v] (jmm, 06/15/18)
             L7S_b = L7S_b *LWC * 0.08205e+0_fp * TK / PATM



             IF ( L7S_b > 1.0e+15_fp .or. L7S_b < 0 ) THEN
                RC = GC_FAILURE
                PRINT *,'Loc: ',I,J,L
                PRINT *,'L7S_b: ',L7S_b
                PRINT *,'HMSc: ',HMSc
                PRINT *,'in exp: ',XX
                PRINT *,'L7_b: ',L7_b
                PRINT *,'KaqHMS: ',KaqHMS
                PRINT *,'TK: ',TK
                PRINT *,'LWC: ',LWC
                PRINT *,'PATM :',PATM
                !                  CALL ERROR_STOP( 'L6s_b >1e15 or <0,  point 3', LOC )
             ENDIF

             L7S = 0.e+0_fp
             L7S_b = 0.e+0_fp
             L8S = 0.e+0_fp

             !----------------------------------------------------------
             ! Compute HMS loss by aqueous OH. Prevent floating-point exception
             ! by not allowing the exponential to go to infinity if
             ! the argument is too large.  (jmm, 06/27/18)
             !
             ! HMS oxidation by aqueous OH is converted to v/v/s in the KaqHMS2
             ! term. See aqchem_so2 subroutine for reaction details.
             !
             ! This is a simplified version of a radical chain.
             ! L8S consumes 1 HMS, 1 OH, and 1 SO2 and produces 2 SO4--
             !----------------------------------------------------------

             ! Convert KaqHMS2 from [m^6 kg^-2 s^-1] to [v/v/s]
             KaqHMS2 = KaqHMS2 * State_Met%AIRDEN(I,J,L) * &
                            State_Met%AIRDEN(I,J,L)
             ! Argument of the exponential
             XX  = ( HMS0 - OH0 ) * KaqHMS2 * DTCHEM

             ! Test if EXP(XX) can be computed w/o numerical exception
             IF ( IS_SAFE_EXP( XX ) .and. ABS( XX ) > 0e+0_fp ) THEN

                ! Aqueous phase HMS loss rate w/ OH(aq) [v/v/timestep]
                L8  = EXP( XX )

                ! Loss by OH(aq)
                L8S = HMS0 * OH0 * ( L8 - 1.e+0_fp ) / &
                               ( (HMS0 * L8) - OH0 )
             ELSE

                ! NOTE from Jintai Lin (4/28/10):
                ! However, in the case of a negative XX, L8S should be
                ! approximated as , instead of HCHO0. In other words,
                ! L8S = HMS0 * OH0 * ( L7 - 1.D0 ) / ( (HMS0*L8) - OH0 )
                ! reaches different limits when XX reaches positive infinity
                ! and negative infinity.
                IF ( XX > 0.e+0_fp ) THEN
                   L8S = OH0
                ELSE IF ( XX < 0.e+0_fp) THEN
                   L8S = HMS0
                ELSE
                   !(qjc, 04/10/16) different solution when HMS0 = OH0
                   L8S = HMS0 - 1/(KaqHMS2*DTCHEM+1/HMS0)
                ENDIF

             ENDIF

          ELSE
             L7S = 0.e+0_fp
             L7S_b = 0.e+0_fp
             L8S = 0.e+0_fp
          ENDIF


          L2S   =  L2S   * FC
          L3S   =  L3S   * FC
          L3S_1 =  L3S_1 * FC   !(qjc, 06/20/16)
          L4S   =  L4S   * FC
          L5S   =  L5S          !(qjc, 06/20/16), do not multiply by FC if it
                                ! is not divided by FC at the begining
          L5S_1 =  L5S_1        !(qjc, 06/20/16)
          L6S   =  L6S
          L6S_1 =  L6S_1

          L7S   =  L7S   * FC   ! Note: consumes SO2 but for HMS (jmm, 06/13/18)
          L7S_b =  L7S_b * FC   ! Note: releases SO2 from HMS (jmm, 06/13/18)
          L8S   =  L8S   * FC   ! Notes: releases 2 sulfate and 1 HCHO for 1 HMS and 1 SO2 (jmm, 06/29/18)


          ! make sure HMS produced is less than HCHO available
          ! (jmm, 06/28/18)
          LSTOT_HMS = L7S - L7S_b - L8S
          IF (LSTOT_HMS > HCHO0) THEN
             L7S = HCHO0 + L7S_b + L8S
          ELSE
             L7S = L7S
          ENDIF

          ! make sure HMS decomposed to HCHO and SO2 is less than HMS available
          ! (jmm, 06/28/18)
          IF (-LSTOT_HMS > HMS0) THEN
             L7S_b = ( HMS0 + L7S ) * L7S_b / ( L7S_b + L8S )
             L7S   = ( HMS0 + L7S ) * L8S   / ( L7S_b + L8S )
          ELSE
             L7S_b = L7S_b
             L8S   = L8S
          ENDIF

          LSTOT0 = L2S + L3S + L4S + L5S + L6S + L7S - L7S_b + L8S ! (qjc, 11/04/16), (jmm, 06/13/18)

          ! make sure sulfate produced is less than SO2 available
          ! (qjc, 06/20/16)
          IF (LSTOT0 > SO2_ss) THEN

             L2S   = SO2_ss * L2S   / LSTOT0
             L3S   = SO2_ss * L3S   / LSTOT0
             L4S   = SO2_SS * L4S   / LSTOT0
             L5S   = SO2_ss * L5S   / LSTOT0
             L3S_1 = SO2_ss * L3S_1 / LSTOT0
             L5S_1 = SO2_ss * L5S_1 / LSTOT0
             L6S_1 = SO2_ss * L6S_1 / LSTOT0
             L7S   = SO2_SS * L7S   / LSTOT0 ! (jmm, 06/13/18)
             L8S   = SO2_SS * L8S   / LSTOT0 ! (jmm, 06/29/18)


             ! This is the ratio used to calculate the actual removal of
             ! HOBr by SO2 for use in gckpp_HetRates.F90 (qjc, 06/20/16)
             State_Chm%fupdateHOBr(I,J,L) = fupdateHOBr_0 * SO2_ss/ LSTOT0
             State_Chm%fupdateHOCl(I,J,L) = fupdateHOCl_0 * SO2_ss/ LSTOT0

          ELSE

             L2S   = L2S
             L3S   = L3S
             L4S   = L4S
             L5S   = L5S
             L6S   = L6S
             L3S_1 = L3S_1
             L5S_1 = L5S_1
             L6S_1 = L6S_1
             L7S   = L7S  ! (jmm, 06/13/18)
             L8S   = L8S  ! (jmm, 06/29/18)


             ! This is the ratio used to calculate the actual removal of
             ! HOBr by SO2 (qjc, 06/20/16)
             State_Chm%fupdateHOBr(I,J,L) = fupdateHOBr_0
             State_Chm%fupdateHOCl(I,J,L) = fupdateHOCl_0

          ENDIF

          ! Decide whether or not it is necessary to use heterogeneous cloud
          ! pH calculations based on the Fahey and Pandis, 2001 decision
          ! algorithm (bec, 12/23/11)

          ! Add up total seasalt and dust and convert to ug/m3
          ! Note that it is better to use dust and sea-salt alkalinity
          ! tracers if these are being transported (bec, 12/23/11)

#ifdef TOMAS
          !%%%%%%%%%%%%%%%%% BUG FIX FOR TOMAS %%%%%%%%%%%%%%%%%%%%%%%
          ! NOTE: TOMAS uses its own dust tracers and does not
          ! carry ALKdst.  Set ALKdst to zero here. (bmy, 1/28/14)
          ALKdst = 0e+0_fp
#else

          ! For other simulations, Sum up the contributions from
          ! DST1 thru DST4 tracers into ALKdst. (bmy, 1/28/14)
          ALKdst = ( Spc(id_DST1)%Conc(I,J,L) + Spc(id_DST2)%Conc(I,J,L) +            &
                     Spc(id_DST3)%Conc(I,J,L) + Spc(id_DST4)%Conc(I,J,L) ) *          &
                     1.e+9_fp * State_Met%AD(I,J,L)                       &
                     / ( AIRMW / State_Chm%SpcData(id_DST1)%Info%MW_g ) &
                     / State_Met%AIRVOL(I,J,L)
#endif

          ALKss  = ( Spc(id_SALA)%Conc(I,J,L) + Spc(id_SALC)%Conc(I,J,L) ) *        &
                     1.e+9_fp * State_Met%AD(I,J,L)                       &
                     / ( AIRMW / State_Chm%SpcData(id_SALA)%Info%MW_g ) &
                     / State_Met%AIRVOL(I,J,L)

          ALKds = ALKdst + ALKss

          ! Get NH3 concentrations (v/v)
          NH3 = Spc(id_NH3)%Conc(I,J,L)

          ! Initialize
          BULK = 0
          SIZE_RES = 0

          ! Fahey and Seinfeld decision algorithm
          IF ( H2O20 > SO2_ss + 1e-9_fp ) THEN
             BULK = 1
          ELSEIF( LWC < 0.1e-6_fp ) THEN !10^-6 coversion from g/m3 --> m3/m3
             SIZE_RES = 1
          ELSEIF( gno3 > NH3 ) THEN
             IF ( SO2_ss >= 5.e-9_fp          .and. &
                  H2O20  >= SO20   )                &
                  BULK    = 1
             IF ( LWC    >= 0.3e-6_fp         .and. &
                  SO2_ss >= 3.e-9_fp          .and. &
                  H2O20  >= SO2_ss )                &
                  BULK    = 1
             IF ( ALKds  >= 5.e+0_fp          .and. &
                  LWC    >= 0.5e-6_fp         .and. &
                  H2O20  >= SO2_ss )                &
                  BULK    = 1
             IF ( LWC    >= 0.1e-6_fp         .and. &
                  gno3   <= (NH3 + 2.e-9_fp) )      &
                  BULK    = 1
          ELSEIF( LWC    >= 0.5e-6_fp ) THEN
             IF ( H2O20  >= (0.9e+0_fp * SO2_ss) )  &
                  BULK    = 1
             IF ( NH3    <= 1.e-9_fp          .and. &
                  ALKds  >= 5.e+0_fp          .and. &
                  SO2_ss <= 10.e-9_fp )             &
                  BULK    = 1
          ELSEIF( LWC    >= 0.3e-6_fp ) THEN
             IF ( NH3    >= (gno3 + 5.e-9_fp) .and. &
                  SO2_ss <= 10.e-9_fp )             &
                  BULK    = 1
             IF ( gno3   <= 1.e-9_fp          .and. &
                  NH3    >= (gno3 + 2.e-9_fp) )     &
                  BULK    = 1
             IF ( gno3   <= 7.e-9_fp          .and. &
                  NH3    >= (gno3 + 3.e-9_fp) )     &
                  BULK    = 1
             IF ( ALKds  >= 3.e+0_fp          .and. &
                  NH3 <= 10e-9_fp             .and. &
                  SO2_ss <= 5e-9_fp )               &
                  BULK    = 1
             IF ( ALKds  >= 5.e+0_fp          .and. &
                  NH3    <= 10.e-9_fp         .and. &
                  SO2_ss <= 5.e-9_fp )              &
                  BULK    = 1
             IF ( SO2_ss >= 1.5e-9_fp         .and. &
                  H2O20  >= SO2_ss )                &
                  BULK    = 1
             IF ( NH3    <= 12.e-9_fp         .and. &
                  ALKds  >=10.e+0_fp )              &
                  BULK    = 1
             IF ( NH3    <= 1.e-9_fp          .and. &
                  ALKds  >= 4.e+0_fp          .and. &
                  SO2_ss <= 10.e-9_fp )             &
                  BULK    = 1
             IF ( NH3    <= 5.e-9_fp          .and. &
                  ALKds  >= 6.e+0_fp          .and. &
                  SO2_ss <= 10.e-9_fp )             &
                  BULK    = 1
             IF ( NH3    <= 7.e-9_fp          .and. &
                  ALKds   >-8.e+0_fp          .and. &
                  SO2_ss <= 10.e-9_fp )             &
                  BULK    = 1
          ELSEIF( LWC    >= 0.1e-6_fp ) THEN
             IF ( NH3    <= 1.e-9_fp          .and. &
                  ALKds  >= 5.e+0_fp   )            &
                  BULK    = 1
             IF ( NH3    <= 5.e-9_fp          .and. &
                  ALKds  >= 10.e+0_fp  )            &
                  BULK    = 1
             IF ( gno3   <= 1.e-9_fp          .and. &
                  NH3    >= (gno3 + 2.e-9_fp) .and. &
                  SO2_ss <= 7.e-9_fp )              &
                  BULK    = 1
             IF ( gno3   <= 1.e-9_fp          .and. &
                  NH3    >= (gno3 + 2.e-9_fp) .and. &
                  ALKds  >= 2.e+0_fp )              &
                  BULK = 1
             IF ( gno3   <= 3.e-9_fp          .and. &
                  NH3    >= (gno3 + 4.e-9_fp) )     &
                  BULK    = 1
             IF ( gno3   <= 7.e-9_fp          .and. &
                  NH3    >= (gno3 + 3.e-9_fp) .and. &
                  SO2_ss <= 5.e-9_fp )              &
                  BULK    = 1
             IF ( gno3   <= 7.e-9_fp          .and. &
                  NH3    >= (gno3 + 3.e-9_fp) .and. &
                  ALKds  >= 4.e+0_fp          .and. &
                  SO2_ss <= 9.e-9_fp  )             &
                  BULK    = 1
             IF ( ALKds  >= 3.e+0_fp          .and. &
                  NH3    <= 3.e-9_fp          .and. &
                  SO2_ss <= 4.e-9_fp )              &
                  BULK    = 1
             IF ( ALKds  >= 5.e+0_fp          .and. &
                  SO2_ss <= 5.e-9_fp          .and. &
                  NH3    <= 7.e-9_fp )              &
                  BULK    = 1
             IF ( NH3    >= (gno3 + 2.e-9_fp) .and. &
                  SO2_ss <= 5.e-9_fp )              &
                  BULK    = 1
             IF ( NH3    >= (gno3 + 4.e-9_fp) .and. &
                  SO2_ss <= 10.e-9_fp )             &
                  BULK    = 1
             IF ( ALKds  >= 2.e+0_fp          .and. &
                  NH3    <= 10.e-9_fp         .and. &
                  H2O20  >= SO2_ss )                &
                  BULK    = 1
             IF ( NH3    <= 1.e-9_fp          .and. &
                  SO2_ss >= 3.e-9_fp          .and. &
                  H2O20  >= SO2_ss )                &
                  BULK    = 1
          ELSE
             SIZE_RES = 1
          ENDIF

          ! Decide whether or not to perform sulfate production rate
          ! enhancement due to cloud drop heterogenity in pH over the oceans
          ! (bec, 12/23/11)
          IF ( SIZE_RES == 1 .AND. State_Met%IsWater(I,J) .AND. &
               TK > 268.15 ) THEN

             ! Get total in-cloud sulfate production based on bulk cloud pH
             ! calculations for use in HET_DROP_CHEM
             ! added in SO4 from HMS (jmm, 06/29/18)
             LSTOT = (L2S + L3S + L4S + L5S + L6S + L8S + L8S) / FC !(qjc, 06/20/16)

             ! Get coarse-mode sea-salt concentration for use in
             ! HET_DROP_CHEM [v/v]
             ! Note that it is better to use coarse sea salt alkalinity
             ! tracer if it is being transported (bec, 12/23/11)
             SSCvv  = Spc(id_SALC)%Conc(I,J,L)

             ! Get sulfate concentrations for use in HET_DROP_CHEM [v/v]
             aSO4  =  Spc(id_SO4)%Conc(I,J,L)

             ! This is to make sure HET_DROP_CHEM does not compute more
             ! sulfate then there is SO2
             ! Add L5S (qjc, 11/04/16)
             ! Add L7S, L7S_b, and L8S (jmm, 06/29/18)
             SO2_sr = MAX( SO2_ss - ( L2S + L3S + L4S + L5S + L6S + &
                  L7S - L7S_b + L8S), MINDAT)

             CALL HET_DROP_CHEM( I,    J,   L,      LSTOT, SSCvv, &
                                 aSO4, NH3, SO2_sr, H2O20, GNO3,  SR, &
                                 Input_Opt, State_Met, State_Chm )

             ! Henry's Law constant of O3 and HOBr (M atm-1)
             HCO3   = 1.13e-2_fp * EXP( 8.51e+0_fp * &
                      ( 298.15e+0_fp / TK - 1.e+0_fp ) )
             HCHOBr = 1.3e+3_fp
             HCHOCl = 6.6e+2_fp

             ! Rate coefficient (M-1 s-1)
             KO3   = 7.32e+14_fp * EXP( -4.03e+3_fp / TK ) ! for O3+SO3
             KHOBr = 5.0e+9_fp                             ! for HOBr+SO3
             KHOCl = 7.6e+8_fp                             ! for HOCl+SO3
             ! Make sure we don't divide by zero (ckeller, 1/25/18)
             TMP = KHOBr*HOBr0*HCHOBr + KO3*O3*HCO3 + KHOCl*HOCl0*HCHOCl
             IF ( TMP > 0.0_fp ) THEN
                  f_srhobr = KHOBr*HOBr0*HCHOBr / TMP
                  f_srhocl = KHOCl*HOCl0*HCHOCl / TMP

                  SRhobr   = SR * f_srhobr
                  SRhocl = SR * f_srhocl
                  SRo3   = SR * (1-f_srhobr-f_srhocl)
             ELSE
                SR     = 0.e+0_fp
                SRhobr = 0.e+0_fp
                SRhocl = 0.e+0_fp
                SRo3   = 0.e+0_fp
                so2_sr = 0.e+0_fp
             ENDIF
          ELSE
             SR     = 0.0_fp
             SRhobr = 0.0_fp
             SRhocl = 0.e+0_fp
             SRo3   = 0.0_fp
             so2_sr = 0.0_fp
             LSTOT  = 0.0_fp
          ENDIF

          ! We have used the in-cloud LWC to compute the sulfate
          ! aqueous chemistry.  We get the appropriate grid-box averaged
          ! mass of SO2 and sulfate by multiplying the reaction rates
          ! L2S and L3s by the cloud fraction after the aqueous chemistry
          ! has been done.  (lzh, jaf, bmy, 5/27/11)
          SR     =  SR     * FC
          SRo3   =  SRo3   * FC !(qjc, 04/10/16)
          SRhobr =  SRhobr * FC !(qjc, 04/10/16)
          SRhocl =  SRhocl * FC

          ! Store initial SO2_ss (qjc, 06/20/16)
          SO2_ss0 = SO2_ss

          ! Make sure SO2_ss and H2O20 are in the proper range
          ! Add L5S (qjc, 11/04/16)
          ! Add L7S and L7S_b (jmm, 06/13/18)
          ! Add loss and production of HCHO (jmm, 06/13/18)
          ! Add loss and production of HMS (jmm, 06/15/18)
          SO2_ss = MAX( SO2_ss - ( L2S+L3S+L4S+L5S+L6S+ &
               L7S-L7S_b+L8S+SR ), MINDAT )
          H2O20  = MAX( H2O20  - L2S,                    MINDAT )
          HCHO0  = MAX( HCHO0  -  L7S + L7S_b + L8S, MINDAT )
          HMS0  = MAX( HMS0  +  L7S - L7S_b - L8S, MINDAT )
          O3  = MAX( O3 - L3S, MINDAT )
          OH0 = MAX( OH0 - L8S, MINDAT )


          ! Factor to calculate effective SO3 and HSO3 in aqchem_so2 to
          ! be used in HOBr+HSO3/SO3 (qjc, 06/20/16)
          rSIV = SO2_sr/SO2_ss0
          rSIV = MAX(rSIV, 0.e+0_fp)
          rSIV = MIN(rSIV, 1.e+0_fp)

          ! Store HSO3aq, SO3aq for use in gckpp_HetRates.F90
          State_Chm%HSO3_AQ(I,J,L) = HSO3aq * ( 1.0_fp + rSIV ) / 2.0_fp
          State_Chm%SO3_AQ(I,J,L)  = SO3aq  * ( 1.0_fp + rSIV ) / 2.0_fp

          ! Update SO2 level, save SO2[ppv], H2O2[ppv] for WETDEP
          SO2s( I,J,L) = SO2_ss
          H2O2s(I,J,L) = H2O20

          ! SO2 chemical loss rate  = SO4 production rate [v/v/timestep]
          ! Add HMS (jmm, 06/13/18)
          PSO4_SO2(I,J,L) = LSTOT + PSO4E
          PSO4_ss (I,J,L) = PSO4F
          PHMS_SO2(I,J,L) = L7S
          PSO2_HMS(I,J,L) = L7S_b
          PSO4_HMS(I,J,L) = L8S + L8S


       ELSE

          ! Otherwise, don't do aqueous chemistry, and
          ! save the original concentrations into SO2 and H2O2
          H2O2s(I,J,L) = MAX( H2O20,  1.0e-32_fp )
          SO2s(I,J,L ) = MAX( SO2_ss, 1.0e-32_fp )
          L2S          = 0.e+0_fp
          L3S          = 0.e+0_fp
          L3S_1        = 0.e+0_fp !(qjc, 11/04/16)
          L4S          = 0.e+0_fp
          L5S          = 0.e+0_fp !(qjc, 11/04/16)
          L5S_1        = 0.e+0_fp !(qjc, 11/04/16)
          L6S          = 0.e+0_fp
          L6S_1        = 0.e+0_fp
          SR           = 0.e+0_fp
          SRhobr       = 0.e+0_fp
          SRhocl       = 0.e+0_fp
          SRo3         = 0.e+0_fp
          HPLUS        = 0.e+0_fp
          L7S          = 0.e+0_fp !(jmm, 06/13/18)
          L7S_b        = 0.e+0_fp !(jmm, 06/13/18)
          L8S          = 0.e+0_fp !(jmm, 06/13/18)

          ! Store HSO3aq, SO3aq for use in gckpp_HetRates.F90
          ! Avoid divide-by-zero errors
          State_Chm%HSO3_AQ(I,J,L)     = 1.0e-32_fp
          State_Chm%SO3_AQ(I,J,L)      = 1.0e-32_fp

          ! This is the ratio used to calculate the actual removal of
          ! HOBr by SO2 for use in gckpp_HetRates.F90 (qjc, 06/20/16)
          State_Chm%fupdateHOBr(I,J,L) = 0.e+0_fp
          State_Chm%fupdateHOCl(I,J,L) = 0.e+0_fp

          ! SO2 chemical loss rate  = SO4 production rate [v/v/timestep]
          ! Add HMS (jmm, 06/13/18)
          PSO4_SO2(I,J,L) = PSO4E
          PSO4_ss (I,J,L) = PSO4F
          PHMS_SO2(I,J,L) = 0.0e+0_fp
          PSO2_HMS(I,J,L) = 0.0e+0_fp
          PSO4_HMS(I,J,L) = 0.0e+0_fp


       ENDIF

       ! Store updated SO2, H2O2 back to the tracer arrays
       ! Add HCHO, OH, O3, and HMS (jmm, 06/13/18)
       If (FullRun) Then
          Spc(id_SO2)%Conc(I,J,L)  = SO2s( I,J,L)
          Spc(id_H2O2)%Conc(I,J,L) = H2O2s(I,J,L)

          IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
             Spc(id_OH)%Conc(I,J,L)   = OH0
             Spc(id_O3)%Conc(I,J,L)   = O3
             Spc(id_CH2O)%Conc(I,J,L) = HCHO0
             Spc(id_HMS)%Conc(I,J,L)  = HMS0
          ENDIF

          ! Set SO4H1 and SO4H2 to zero at end of each timestep
          IF ( id_SO4H1 > 0 ) Spc(id_SO4H1)%Conc(I,J,L) = 0.0e+0_fp
          IF ( id_SO4H2 > 0 ) Spc(id_SO4H2)%Conc(I,J,L) = 0.0e+0_fp
          IF ( id_SO4H3 > 0 ) Spc(id_SO4H3)%Conc(I,J,L) = 0.0e+0_fp
          IF ( id_SO4H4 > 0 ) Spc(id_SO4H4)%Conc(I,J,L) = 0.0e+0_fp
       End If

       ! SO2 chemical loss rate  = SO4 production rate [v/v/timestep]
       ! Add L5S (qjc, 11/04/16)
       ! Add L8S (jmm, 06/29/18)
       PSO4_SO2(I,J,L) = L1 + L2S + L3S + L4S + L5S + PSO4E + SR + &
            L6S + L8S + L8S

       ! Production of sulfate on sea salt
       PSO4_ss (I,J,L) = PSO4F

#ifdef APM
       PSO4_SO2APM(I,J,L) = L2S + L3S + L4S + L5S + SR
       PSO4_SO2APM(I,J,L) = MAX(0.d0,PSO4_SO2APM(I,J,L))
       PSO4_SO2SEA(I,J,L) = PSO4E + PSO4F
       PSO4_SO2SEA(I,J,L) = MAX(0.d0,PSO4_SO2SEA(I,J,L))
#endif
#ifdef TOMAS
       PSO4_SO2AQ(I,J,L) = L2S + L3S + SR ! For TOMAS microphysics
#endif

       ! tdf Production of sulfate and nitrate on dust
       IF ( LDSTUP ) THEN

          ! NB Fine dust mass excluded from PSO4_SO2 - kept separately
          ! (tdf 07/24/08)
          ! tdf PNIT_d, PH2SO4_d, and PSO4_d computed in DUST_CHEM

          DO IBIN = 1, NDSTBIN
             ! included P(SO4) due to uptake of H2SO4(g)        !tdf 3/2/2K9
             PSO4_dust(I,J,L,IBIN) = PSO4_d(IBIN) + PH2SO4_d(IBIN)
             PNIT_dust(I,J,L,IBIN) = PNIT_d(IBIN)
          END DO

          ! tdf Subtract from PSO4_SO2 that which is now diverted to dust
          DO IBIN = 1, NDSTBIN
             PSO4_SO2(I,J,L) =  PSO4_SO2(I,J,L) - PH2SO4_d(IBIN)
          END DO

       ENDIF    !tdf end of if (LDSTUP) condition

       !=================================================================
       ! HISTORY (aka netCDF diagnostics)
       !=================================================================
       IF ( FullRun ) THEN

          ! P(SO4) from gas-phase oxidation [kg S/s]
          IF ( State_Diag%Archive_ProdSO4fromGasPhase ) THEN
             State_Diag%ProdSO4fromGasPhase(I,J,L) = &
                  ( L1  * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
          ENDIF

          ! P(SO4) from aqueous-phase oxidation with H2O2 [kg S/s]
          IF ( State_Diag%Archive_ProdSO4fromH2O2inCloud ) THEN
             State_Diag%ProdSO4fromH2O2inCloud(I,J,L) = &
                  ( L2S * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
          ENDIF

          ! P(SO4) from aqueous-phase oxidation with O3 [kg S/s]
          IF ( State_Diag%Archive_ProdSO4fromO3InCloud ) THEN
             State_Diag%ProdSO4fromO3InCloud(I,J,L) = &
                  ( ( L3S + SRo3 ) * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
          ENDIF

          ! P(SO4) from aqueous-phase oxidation with HOBr [kg S/s]
          IF ( State_Diag%Archive_ProdSO4fromHOBrInCloud ) THEN
             State_Diag%ProdSO4fromHOBrInCloud(I,J,L) = &
                  ( ( L5S + SRhobr ) * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
          ENDIF

          ! P(SO4) from aqueous-phase oxidation with O2 metal-catalyzed
          ! [kg S/s]
          IF ( State_Diag%Archive_ProdSO4fromO2InCloudMetal ) THEN
             State_Diag%ProdSO4fromO2InCloudMetal(I,J,L) = &
                  ( L4S * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
          ENDIF

          ! P(SO4) from O3 in sea salt aerosol [kg S/s]
          IF ( State_Diag%Archive_ProdSO4fromO3inSeaSalt ) THEN
             State_Diag%ProdSO4fromO3inSeaSalt(I,J,L) = &
                  ( ( PSO4E + PSO4F ) * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
          ENDIF

          ! P(SO4) by SRo3 [kg S/s]
          IF ( State_Diag%Archive_ProdSO4fromSRO3 ) THEN
             State_Diag%ProdSO4fromSRO3(I,J,L) = &
                  ( SRo3 * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
          ENDIF

          ! P(SO4) by SRhobr [kg S/s]
          IF ( State_Diag%Archive_ProdSO4fromSRHOBr ) THEN
             State_Diag%ProdSO4fromSRHOBr(I,J,L) = &
                  ( SRhobr * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
          ENDIF

          ! P(SO4) by o3s [kg S/s]
          IF ( State_Diag%Archive_ProdSO4fromO3s ) THEN
             State_Diag%ProdSO4fromO3s(I,J,L) = &
                  ( L3S_1 * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
          ENDIF

          ! P(HMS) by SO2 and HCHO [kg S/s]
          ! jmm (06/13/18)
          IF ( State_Diag%Archive_ProdHMSfromSO2andHCHOinCloud ) THEN
             State_Diag%ProdHMSfromSO2andHCHOinCloud(I,J,L) = &
                  ( L7S * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
          ENDIF

          ! P(SO2 and HCHO) from HMS decomp [kg S/s]
          ! jmm (06/13/18)
          IF ( State_Diag%Archive_ProdSO2andHCHOfromHMSinCloud ) THEN
             State_Diag%ProdSO2andHCHOfromHMSinCloud(I,J,L) = &
                  ( L7S_b * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
          ENDIF

          ! P(SO4) from aqueous-phase oxidation of HMS [kg S/s]
          ! jmm (07/06/18)
          IF ( State_Diag%Archive_ProdSO4fromHMSinCloud ) THEN
             State_Diag%ProdSO4fromHMSinCloud(I,J,L) = &
                  ( 2.e+0_fp * L8S * State_Met%AD(I,J,L) / &
                  TCVV_S ) / DTCHEM
          ENDIF

          !-----------------------------------------------------------
          ! Diagnostics for acid uptake on dust aerosol simulations
          !-----------------------------------------------------------
          IF ( LDSTUP ) THEN

             ! Zero
             PSO4d_tot   = 0.e+0_fp
             PH2SO4d_tot = 0.e+0_fp
             PNITd_tot   = 0.e+0_fp

             DO IBIN = 1, NDSTBIN
                PSO4d_tot   = PSO4d_tot   + PSO4_d(IBIN)
                PNITd_tot   = PNITd_tot   + PNIT_d(IBIN)
                PH2SO4d_tot = PH2SO4d_tot + PH2SO4_d(IBIN)
             END DO

             ! P(SO4) from O3 oxidation on dust aerosols [kg S/s]
             IF ( State_Diag%Archive_ProdSO4fromOxidationOnDust ) THEN
                State_Diag%ProdSO4fromOxidationOnDust(I,J,L) = &
                     ( PSO4d_tot * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
             ENDIF

             ! P(NIT) from HNO3 uptake on dust [kg N/s]
             IF ( State_Diag%Archive_ProdNITfromHNO3uptakeOnDust ) THEN
                State_Diag%ProdNITfromHNO3uptakeOnDust(I,J,L) = &
                     ( PNITd_tot * State_Met%AD(I,J,L) / TCVV_N ) / DTCHEM
             ENDIF

             ! P(SO4) from uptake of H2SO4 on dust aerosols [kg S/s]
             IF ( State_Diag%Archive_ProdSO4fromUptakeOfH2SO4g ) THEN
                State_Diag%ProdSO4fromUptakeOfH2SO4g(I,J,L) = &
                     ( PH2SO4d_tot * State_Met%AD(I,J,L) / TCVV_S )/DTCHEM
             ENDIF
          ENDIF
       ENDIF
#ifdef LUO_WETDEP
       ! Luo et al 2020 wtdep
       IF( SUM( State_Chm%QQ3D    (I,J,L:State_Grid%NZ) *                    &
                State_Met%BXHEIGHT(I,J,L:State_Grid%NZ)   ) > 1.D-30 ) THEN

         State_Chm%pHRain(I,J,L) =                                           &
              SUM( State_Chm%pHCloud (I,J,L:State_Grid%NZ)    *              &
                   State_Chm%QQ3D    (I,J,L:State_Grid%NZ)    *              &
                   State_Met%BXHEIGHT(I,J,L:State_Grid%NZ) )  /              &
              SUM( State_Chm%QQ3D    (I,J,L:State_Grid%NZ)    *              &
                   State_Met%BXHEIGHT(I,J,L:State_Grid%NZ) )

         State_Chm%QQpHRain(I,J,L) =                                         &
              SUM( State_Chm%pHCloud (I,J,L:State_Grid%NZ)    *              &
                   State_Chm%QQ3D    (I,J,L:State_Grid%NZ)    *              &
                   State_Met%BXHEIGHT(I,J,L:State_Grid%NZ) )

         State_Chm%QQRain(I,J,L) =                                           &
              SUM( State_Chm%QQ3D    (I,J,L:State_Grid%NZ)    *              &
                   State_Met%BXHEIGHT(I,J,L:State_Grid%NZ) )
       ELSE
          ! Default wetdep scheme
         State_Chm%pHrain(I,J,L)   = 5.6D0
         State_Chm%QQpHrain(I,J,L) = 0.D0
         State_Chm%QQrain(I,J,L)   = 0.D0
       ENDIF
#endif

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Deallocate if allocated
    IF ( ASSOCIATED( NDENS_SALA ) ) DEALLOCATE ( NDENS_SALA )
    IF ( ASSOCIATED( NDENS_SALC ) ) DEALLOCATE ( NDENS_SALC )

    ! Free pointers
    Spc        => NULL()
    !SSAlk      => NULL()
    H2O2s      => NULL()
    SO2s       => NULL()
    NDENS_SALA => NULL()
    NDENS_SALC => NULL()

  END SUBROUTINE CHEM_SO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: seasalt_chem
!
! !DESCRIPTION: Subroutine SEASALT\_CHEM computes SO4 formed from S(IV) + O3 on
!  seasalt aerosols as a function of seasalt alkalinity. (bec, bmy, 4/13/05,
!  10/7/08)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SEASALT_CHEM( I,          J,         L,         &
                           ALK1,       ALK2,      SO2_cd,    &
                           Kt1,        Kt2,       Kt1N,      &
                           Kt2N,       Kt1L,      Kt2L,      &
                           SO2_ss,     PSO4E,                &
                           PSO4F,      AlkA,      AlkC,      &
                           Input_Opt,  State_Met, State_Chm, &
                           State_Diag, FullRun,   RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : SpcConc
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM
    USE TIME_MOD,           ONLY : GET_ELAPSED_SEC
    USE TIME_MOD,           ONLY : GET_MONTH
    USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)   :: I, J, L    ! Grid box indices
    REAL(fp),       INTENT(IN)   :: SO2_cd     ! SO2 mixing ratio [v/v] after
                                               !  gas phase chemistry and
                                               !  dry deposition
    REAL(fp),       INTENT(IN)   :: Kt1, Kt2   ! Rate constant [s-1] for
                                               !  sulfate formation on sea
                                               !  salt aerosols from GET_ALK
                                               !  (1=fine; 2=coarse)
    REAL(fp),       INTENT(IN)   :: Kt1N, Kt2N
    REAL(fp),       INTENT(IN)   :: Kt1L, Kt2L
    REAL(fp),       INTENT(IN)   :: ALK1, ALK2 ! Alkalinity [kg] from
                                               !  seasalt_mod
    TYPE(MetState), INTENT(IN)   :: State_Met  ! Meteorology State object
    TYPE(OptInput), INTENT(IN)   :: Input_Opt  ! Input Options object
    LOGICAL,        INTENT(IN)   :: FullRun    ! Modify species conc?
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT)   :: SO2_ss    ! SO2 mixing ratio [v/v]
                                               !  after sea salt chemistry
    REAL(fp),       INTENT(OUT)   :: PSO4E     ! SO4E (sulfate produced by
                                               !  S(IV)+O3 on fine seasalt)
                                               !  mixing ratio [v/v]
    REAL(fp),       INTENT(OUT)   :: PSO4F     ! SO4F (sulfate produced by
                                               !  S(IV)+O3 on coarse seasalt)
    REAL(fp),       INTENT(OUT)   :: AlkA      ! Modified SSA alkalinity [v/v]
    REAL(fp),       INTENT(OUT)   :: AlkC      ! Modified SSA alkalinity [v/v]
    INTEGER,        INTENT(OUT)   :: RC        ! Success or failure?
!
! !REMARKS:
!  Chemical reactions:
!  ============================================================================
!  (R1) SO2 + O3 + ALK => SO4 + O2
!       Modeled after Chamedies and Stelson, 1992?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER :: MINDAT    = 1.0e-20_fp
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL(fp)            :: SO2_chem,    DTCHEM
    REAL(fp)            :: EQ_1_C,      EQ_2_C
    REAL(fp)            :: SO4E,        SO2_new,    SO4F
    REAL(fp)            :: SO2_eq,      N_FLUX_A,   N_FLUX_C
    REAL(fp)            :: END_ALK,     L5A,        L5C
    REAL(fp)            :: EQ1,         EQ2,        TITR_SO2
    REAL(fp)            :: TITR_HNO3,   NIT_vv,     NITs_vv
    REAL(fp)            :: NIT0,        NITS0
    REAL(fp)            :: F_SO2,       FALK_A_SO2, FALK_C_SO2
    REAL(fp)            :: EQ_BEG,      F_SO2_A,    F_SO2_C
    REAL(fp)            :: TOTAL_ACID_FLUX
    REAL(fp)            :: HNO3_EQ,     TOT_FLUX_A, TOT_FLUX_C
    REAL(fp)            :: FALK_A_HNO3, HNO3_vv
    REAL(fp)            :: FALK_C_HNO3, F_HNO3_A,   F_HNO3_C
    REAL(fp)            :: EQ_1_N,      EQ_2_N,     F_HNO3
    REAL(fp)            :: HNO3_SSA,    HNO3_SSC,   N_FLUX
    REAL(fp)            :: HNO3_EQ_C,   L6A,        L6C
    REAL(fp)            :: C_FLUX_A,    C_FLUX_C,   C_FLUX
    REAL(fp)            :: HNO3_ss,     HNO3_kg
    REAL(fp)            :: MW_SAL1,     MW_SAL2  ! for salinity/alkalinity
    REAL(fp)            :: L_FLUX_A,    L_FLUX_C,   TITR_HCl
    REAL(fp)            :: ACL_vv,      CCL_vv,     ACL0
    REAL(fp)            :: CCL0,        HCl_eq,     FALK_A_HCl
    REAL(fp)            :: FALK_C_HCl,  HCl_vv,     F_HCl_A
    REAL(fp)            :: F_HCl_C,     EQ_1_L,     EQ_2_L
    REAL(fp)            :: F_HCl,       L_FLUX,     HCl_SSA
    REAL(fp)            :: HCl_ss,      HCl_kg,     HCl_SSC
    REAL(fp)            :: L7A,         L7C,        HCl_EQ_C

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)
    REAL(fp), POINTER   :: AD(:,:,:)
    REAL(fp), POINTER   :: AIRVOL(:,:,:)

    !=================================================================
    ! SEASALT_CHEM begins here!
    !=================================================================

    ! Assume success
    RC      = GC_SUCCESS

    ! Initialize pointers
    Spc     => State_Chm%Species
    AD      => State_Met%AD
    AIRVOL  => State_Met%AIRVOL

    ! Uncomment if transporting salinity/alkalinity as needed
    MW_SAL1 =  State_Chm%SpcData(id_SALAAL)%Info%MW_g
    MW_SAL2 =  State_Chm%SpcData(id_SALCAL)%Info%MW_g

    ! DTCHEM is the chemistry timestep in seconds
    DTCHEM  = GET_TS_CHEM()

    ! Convert SO2 [v/v] to [eq/gridbox]
    ! Remove species molecular weights from equation (bmy, 2/10/17)
    SO2_eq  = ( ( 2.0_fp * SO2_cd * AD(I,J,L) ) / AIRMW  ) * 1000.0_fp
    SO2_eq  = MAX( SO2_eq, MINDAT )

    ! Get the HNO3 and HCl concentration [v/v], either from the species
    ! array (fullchem sims) or from HEMCO (aerosol-only sims)
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
       HNO3_vv = Spc(id_HNO3)%Conc(I,J,L)
       HCl_vv  = Spc(id_HCL)%Conc(I,J,L)
    ELSE
       HNO3_vv = GLOBAL_HNO3(I,J,L)
       HCl_vv  = GLOBAL_HCl(I,J,L)
    ENDIF

    ! Convert HNO3 and HCl [v/v] to [equivalents]
    HNO3_eq = ( ( HNO3_vv * AD(I,J,L) ) / AIRMW ) * 1000.0_fp
    HCl_eq  = ( ( HCl_vv  * AD(I,J,L) ) / AIRMW ) * 1000.0_fp

    !-----------
    ! SO2
    !-----------

    ! Available flux of SO2 to accum sea salt aerosols [v/v/timestep]
    L5A      = EXP( -Kt1 * DTCHEM )
    F_SO2_A  = SO2_cd * ( 1.e+0_fp - L5A )
    F_SO2_A  = MAX( F_SO2_A, 1.e-32_fp )

    ! Convert to [eq/timestep]
    ! Remove species molecular weight from equation (bmy, 2/10/17)
    C_FLUX_A = ( 2.0_fp * F_SO2_A * AD(I,J,L) / AIRMW ) * 1000.0_fp

    ! Available flux of SO2 to coarse sea salt aerosols [v/v/timestep]
    L5C      = EXP( - Kt2 * DTCHEM )
    F_SO2_C  = SO2_cd * ( 1.e+0_fp - L5C )
    F_SO2_C  = MAX( F_SO2_C, 1.0e-32_fp )

    ! Convert to [eq/timestep]
    ! Remove species molecular weight from equation (bmy, 2/10/17)
    C_FLUX_C = ( 2.0_fp * F_SO2_C * AD(I,J,L) / AIRMW ) * 1000.0_fp

    ! Total flux of SO2 [v/v/timestep]
    F_SO2    = F_SO2_A + F_SO2_C

    ! Total flux of SO2 [eq/timestep]
    C_FLUX   = C_FLUX_A + C_FLUX_C

    !-----------
    ! HNO3
    !-----------

    ! Available flux of HNO3 to accum sea salt aerosols [v/v/timestep]
    L6A = EXP( - Kt1N * DTCHEM )
    F_HNO3_A = HNO3_vv * ( 1.e+0_fp - L6A )
    F_HNO3_A = MAX( F_HNO3_A, 1.0e-32_fp )

    ! Convert to [eq/timestep]
    ! Remove species molecular weight from equation (bmy, 2/10/17)
    N_FLUX_A = ( F_HNO3_A * AD(I,J,L) / AIRMW ) * 1000.0_fp

    ! Available flux of HNO3 to coarse sea salt aerosols [v/v/timestep]
    L6C = EXP( - Kt2N * DTCHEM )
    F_HNO3_C = HNO3_vv * ( 1.e+0_fp - L6C )
    F_HNO3_C = MAX( F_HNO3_C, 1.0e-32_fp )

    ! convert to [eq/timestep]
    ! Remove species molecular weight from equation (bmy, 2/10/17)
    N_FLUX_C = ( F_HNO3_C * AD(I,J,L) / AIRMW ) * 1000.0_fp

    ! Total flux of HNO3
    F_HNO3 = F_HNO3_A + F_HNO3_C ![v/v/timestep]
    N_FLUX = N_FLUX_A + N_FLUX_C ![eq/timestep]

    !-----------
    ! HCl
    !-----------

    ! Available flux of HCl to accum sea salt aerosols [v/v/timestep]
    L7A = EXP( - Kt1L * DTCHEM )
    F_HCl_A = HCl_vv * ( 1.e+0_fp - L7A )
    F_HCl_A = MAX( F_HCl_A, 1.0e-32_fp )

    ! Convert to [eq/timestep]
    ! Remove species molecular weight from equation (bmy, 2/10/17)
    L_FLUX_A = ( F_HCl_A * AD(I,J,L) / AIRMW ) * 1000.0_fp

    ! Available flux of HCl to coarse sea salt aerosols
    ! [v/v/timestep]
    L7C = EXP( - Kt2L * DTCHEM )
    F_HCl_C = HCl_vv * ( 1.e+0_fp - L7C )
    F_HCl_C = MAX( F_HCl_C, 1.0e-32_fp )

    ! convert to [eq/timestep]
    ! Remove species molecular weight from equation (bmy, 2/10/17)
    L_FLUX_C = ( F_HCl_C * AD(I,J,L) / AIRMW ) * 1000.0_fp

    ! Total flux of HCl
    F_HCl  = F_HCl_A  + F_HCl_C ![v/v/timestep]
    L_FLUX = L_FLUX_A + L_FLUX_C ![eq/timestep]

    !-----------
    ! Acid
    !-----------

    ! Total acid flux to accum sea-salt aerosols [eq/box/timestep]
    TOT_FLUX_A = C_FLUX_A + N_FLUX_A + L_FLUX_A
    TOT_FLUX_A = MAX( TOT_FLUX_A, MINDAT )

    ! Total acid flux to coarse sea-salt aerosols [eq/box/timestep]
    TOT_FLUX_C = C_FLUX_C + N_FLUX_C + L_FLUX_C
    TOT_FLUX_C = MAX( TOT_FLUX_C, MINDAT )

    ! Total  acid flux to sea salt aerosols
    TOTAL_ACID_FLUX = TOT_FLUX_A + TOT_FLUX_C

    ! Total available alkalinity [eq]
    !----------------------------------------------------------------------
    ! From Alexander et al., buffering capacity (or alkalinity) of sea-salt
    ! aerosols is equal to 0.07 equivalents per kg dry sea salt emitted
    ! Gurciullo et al., 1999. JGR 104(D17) 21,719-21,731.
    ! tdf
    !----------------------------------------------------------------------
    EQ1 = ALK1 * 0.07e+0_fp
    EQ2 = ALK2 * 0.07e+0_fp

    !----------------------------------------------------------------------
    ! NOTE: This was a sensitivity simulation, keep for future reference
    !       cf Alexander et al 2005 (bec, bmy, 4/13/05)
    !! Total available alkalinity [eq] doubled for Sievering run
    !EQ1 = ALK1 * 0.14e+0_fp
    !EQ2 = ALK2 * 0.14e+0_fp
    !----------------------------------------------------------------------

    IF ( TOT_FLUX_A > EQ1 ) THEN

       ! Fraction of alkalinity available for each acid
       FALK_A_SO2  = C_FLUX_A / TOT_FLUX_A
       FALK_A_HNO3 = N_FLUX_A / TOT_FLUX_A
       FALK_A_HCl  = L_FLUX_A / TOT_FLUX_A
       FALK_A_SO2  = MAX( FALK_A_SO2, MINDAT )
       FALK_A_HNO3 = MAX( FALK_A_HNO3, MINDAT )
       FALK_A_HCl  = MAX( FALK_A_HCl, MINDAT )

    ELSE

       FALK_A_SO2  = 1.0e+0_fp
       FALK_A_HNO3 = 1.0e+0_fp
       FALK_A_HCl = 1.0e+0_fp

    ENDIF

    IF ( TOT_FLUX_C > EQ2 ) THEN

       ! Fraction of flkalinity available for each acid
       FALK_C_SO2  = C_FLUX_C/TOT_FLUX_C
       FALK_C_HNO3 = N_FLUX_C/TOT_FLUX_C
       FALK_C_HCl  = L_FLUX_C/TOT_FLUX_C
       FALK_C_SO2  = MAX( FALK_C_SO2, MINDAT )
       FALK_C_HNO3 = MAX( FALK_C_HNO3, MINDAT )
       FALK_C_HCl  = MAX( FALK_C_HCl, MINDAT )

    ELSE

       FALK_C_SO2  = 1.0e+0_fp
       FALK_C_HNO3 = 1.0e+0_fp
       FALK_C_HCl  = 1.0e+0_fp

    ENDIF

    ! Alkalinity available for S(IV) --> S(VI)
    EQ_1_C       = EQ1 * FALK_A_SO2
    EQ_1_C       = MAX( EQ_1_C, MINDAT )
    EQ_1_N       = EQ1 * FALK_A_HNO3
    EQ_1_N       = MAX( EQ_1_N, MINDAT )
    EQ_1_L       = EQ1 * FALK_A_HCl
    EQ_1_L       = MAX( EQ_1_L, MINDAT )

    EQ_2_C       = EQ2 * FALK_C_SO2
    EQ_2_C       = MAX( EQ_2_C, MINDAT )
    EQ_2_N       = EQ2 * FALK_C_HNO3
    EQ_2_N       = MAX( EQ_2_N, MINDAT )
    EQ_2_L       = EQ2 * FALK_C_HCl
    EQ_2_L       = MAX( EQ_2_L, MINDAT )

    !-----------------
    ! Fine Seasalt
    !-----------------

    ! don't produce more SO4 than available ALK or SO2
    SO4E         = MIN( C_FLUX_A, EQ_1_C, SO2_eq )
    SO4E         = MAX( SO4E, MINDAT )

    ! Update SO2 concentration [eq/box]
    SO2_new      = SO2_eq - SO4E
    SO2_new      = MAX( SO2_new, MINDAT )

    !-----------------
    ! Coarse Seasalt
    !-----------------
    IF ( SO2_new > MINDAT ) THEN

       ! don't produce more SO4 than available ALK or SO2
       SO4F      = MIN( C_FLUX_C, SO2_new, EQ_2_C )
       SO4F      = MAX( SO4F, MINDAT )

       !Update SO2 concentration [eq]
       SO2_chem  = SO2_new - SO4F
       SO2_chem  = MAX( SO2_chem, MINDAT )
    ELSE
       SO4F      = MINDAT
       SO2_chem  = MINDAT
    ENDIF

    ! Alkalinity titrated by S(IV) --> S(VI) [eq]
    TITR_SO2     = SO4E + SO4F

    ! Modified SO2 [eq] converted back to [v/v]
    ! Remove species molecular weights from equation (bmy, 2/10/17)
    SO2_ss       = ( SO2_chem * AIRMW / AD(I,J,L) ) / 2000.0_fp
    SO2_ss       = MAX( SO2_ss, MINDAT )

    ! SO4E produced converted from [eq/timestep] to [v/v/timestep]
    ! Remove species molecular weights from equation (bmy, 2/10/17)
    PSO4E        = ( SO4E * AIRMW / AD(I,J,L) ) / 2000.0_fp

    ! SO4F produced converted from [eq/timestep] to [v/v/timestep]
    ! Remove species molecular weights from equation (bmy, 2/10/17)
    !
    ! NOTE: This new equation will correct the prior 3X overestimate
    ! caused  by switching the MW of SO4S from 96 to 31.4 (bmy, 2/10/17)
    PSO4F        = ( SO4F * AIRMW / AD(I,J,L) ) / 2000.0_fp

    ! Alkalinity titrated by HNO3
    HNO3_SSA     = MIN(N_FLUX_A, HNO3_EQ, EQ_1_N)
    HNO3_SSA     = MAX(HNO3_SSA, MINDAT)
    HNO3_EQ_C    = HNO3_EQ - HNO3_SSA
    HNO3_EQ_C    = MAX(HNO3_EQ_C, MINDAT)
    HNO3_SSC     = MIN(N_FLUX_C, HNO3_EQ_C, EQ_2_N)
    HNO3_SSC     = MAX(HNO3_SSC, MINDAT)
    TITR_HNO3    = HNO3_SSA + HNO3_SSC

    ! Alkalinity titrated by HCl
    HCl_SSA     = MIN(L_FLUX_A, HCl_EQ, EQ_1_L)
    HCl_SSA     = MAX(HCl_SSA, MINDAT)
    HCl_EQ_C    = HCl_EQ - HCl_SSA
    HCl_EQ_C    = MAX(HCl_EQ_C, MINDAT)
    HCl_SSC     = MIN(L_FLUX_C, HCl_EQ_C, EQ_2_L)
    HCl_SSC     = MAX(HCl_SSC, MINDAT)
    TITR_HCl    = HCl_SSA + HCl_SSC

    ! HNO3 lost [eq/timestep] converted back to [v/v/timestep]
    IF ( id_HNO3 > 0 ) THEN

       ! Fullchem simulations:  Get HNO3 from the species array
       ! Remove species molecular weights from equation (bmy, 2/10/17)
       !HNO3_ss = ( HNO3_SSC * AIRMW / AD(I,J,L) ) / 1000.0_fp
       ! Should remvoe HNO3 by both SALA and SALC (xnw, 12/8/17)
       HNO3_ss = ( TITR_HNO3 * AIRMW / AD(I,J,L) ) / 1000.0_fp

       If (FullRun) Then
          ! Store back into the species array
          Spc(id_HNO3)%Conc(I,J,L) = MAX( HNO3_vv - HNO3_ss, MINDAT )
       End If

    ELSE

       ! Aerosol-only simulations: Use TITR_HNO3
       ! Remove species molecular weight from equation (bmy, 2/10/17)
       HNO3_ss = ( TITR_HNO3 * AIRMW / AD(I,J,L) ) / 1000.0_fp

    ENDIF

    ! HCl lost [eq/timestep] converted back to [v/v/timestep]
    IF ( FullRun .and. id_HCl > 0 ) THEN
       HCl_ss = ( TITR_HCl * AIRMW / AD(I,J,L) ) / 1000.0_fp
       Spc(id_HCl)%Conc(I,J,L) = MAX( HCl_vv - HCl_ss, MINDAT )
    ENDIF

    !=================================================================
    ! HISTORY (aka netCDF diagnostics)
    !
    ! Loss of HNO3 on sea salt [kg/s]
    !=================================================================
    IF ( FullRun .and. State_Diag%Archive_LossHNO3onSeaSalt ) THEN
       State_Diag%LossHNO3onSeaSalt(I,J,L) = &
            ( HNO3_ss * State_Met%AD(I,J,L) / TCVV_N ) / DTCHEM
    ENDIF

    ! NITS produced converted from [eq/timestep] to [v/v/timestep]
    ! Remove species molecular weight from equation (bmy, 2/10/17)
    !
    ! NOTE: This new equation will correct the prior 2X overestimate
    ! caused  by switching the MW of NITs from 62 to 31.4 (bmy, 2/10/17)
    PNITs(I,J,L) = ( HNO3_SSC * AIRMW / AD(I,J,L) ) / 1000.0_fp

    !Add NIT and Cl productions, xnw 12/8/17
    PNIT(I,J,L) = ( HNO3_SSA * AIRMW / AD(I,J,L) ) / 1000.0_fp
    PACL(I,J,L) = ( HCl_SSA * AIRMW / AD(I,J,L) ) / 1000.0_fp
    PCCL(I,J,L) = ( HCl_SSC * AIRMW / AD(I,J,L) ) / 1000.0_fp

    ! Modified accum alkalinity
    ALKA         = EQ1 - (SO4E + HNO3_SSA + HCl_SSA)
    !ALKA         = MAX( ALKA, MINDAT )

    !------------------------------------------------------------------------
    ! Uncomment this if you want to transport alkalinity (bec, bmy, 4/13/05)
    ![eq] --> [kg] --> [v/v] use this only if transporting alkalinity
    ALKA = (ALKA * ( AIRMW / MW_SAL1) ) / ( 7.0d-2 * AD(I,J,L) )
    !ALKA = MAX( ALKAvv, MINDAT )
    IF (ALKA .LE. MINDAT) ALKA = 0.e+0_fp
    !------------------------------------------------------------------------

    ! Modified accum alkalinity
    ALKC         = EQ2 - (SO4F + HNO3_SSC + HCl_SSC)
    !ALKC         = MAX( ALKC, MINDAT )

    !------------------------------------------------------------------------
    ! Uncomment this if you want to transport alkalinity (bec, bmy, 4/13/05)
    !! [eq] --> [kg] --> [v/v] use this only if transporting alkalinity
    ALKC = (ALKC * ( AIRMW / MW_SAL2 ))/(7.0d-2 * AD(I,J,L))
    !ALKC = MAX(ALKCvv, MINDAT)
    IF (ALKC .LE. MINDAT) ALKC = 0.e+0_fp
    !------------------------------------------------------------------------

    ! Free pointers
    Spc    => NULL()
    AD     => NULL()
    AIRVOL => NULL()

  END SUBROUTINE SEASALT_CHEM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dust_chem
!
! !DESCRIPTION: Subroutine DUST\_CHEM computes SO4 formed from S(IV) + O3 on
!  dust aerosols as a function of dust alkalinity  (tdf 3/28/2K8)
!  Based on routine SEASALT\_CHEM (bec, bmy, 4/13/05, 10/25/05)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE DUST_CHEM( I,         J,         L,         &
                        ALK,       SO2_cd,    H2SO4_cd,  &
                        KTS,       KTN,       KTH,       &
                        SO2_gas,   H2SO4_gas, PSO4d,     &
                        PH2SO4d,   PNITd,     ALKA,      &
                        Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
    USE CMN_SIZE_Mod,       ONLY : NDSTBIN
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : SpcConc
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM,        GET_ELAPSED_SEC
    USE TIME_MOD,           ONLY : GET_ELAPSED_SEC,    GET_MONTH
    USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L       ! Grid box indices
    REAL(fp),       INTENT(IN)    :: SO2_cd        ! SO2 mixing ratio after
                                                   !  gas phase chemistry and
                                                   !  dry deposition [v/v]
    REAL(fp),       INTENT(IN)    :: H2SO4_cd      ! H2SO4 mixing ratio after
                                                   !  gas phase chemistry and
                                                   !  dry deposition [v/v]
    REAL(fp),       INTENT(IN)    :: ALK(NDSTBIN)  ! Dust Alkalinity [v/v]
    REAL(fp),       INTENT(IN)    :: KTS(NDSTBIN)  ! Rate constant for uptake
                                                   !  of SO2 on dust [s-1]
    REAL(fp),       INTENT(IN)    :: KTN(NDSTBIN)  ! Rate constant for uptake
                                                   !  of HNO3 on dust [s-1]
    REAL(fp),       INTENT(IN)    :: KTH(NDSTBIN)  ! Size- and area-weighted
                                                   !  FRACTION for uptake of
                                                   !  H2SO4 on dust
    TYPE(MetState), INTENT(IN)    :: State_Met     ! Meteorology State object
    TYPE(OptInput), INTENT(IN)    :: Input_Opt     ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm     ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT)   :: SO2_gas         ! SO2 mixing ratio after
                                                     !  dust chem [v/v]
    REAL(fp),       INTENT(OUT)   :: PSO4d(NDSTBIN)  ! Sulfate produced by
                                                     !  S(IV)+O3 on dust in
                                                     !  each size bin
    REAL(fp),       INTENT(OUT)   :: H2SO4_gas       ! H2SO4 mixing ratio
                                                     !  after dust chem [v/v]
    REAL(fp),       INTENT(OUT)   :: PNITd (NDSTBIN) ! Nitrate produced by
                                                     !  HNO3 uptake on dust
                                                     !  in each size bin
    REAL(fp),       INTENT(OUT)   :: PH2SO4d(NDSTBIN)! Sulfate produced by
                                                     !  uptake of H2SO4 on
                                                     !  dust in each size bin
    REAL(fp),       INTENT(OUT)   :: ALKA(NDSTBIN)   ! Dust Alkalinity after
                                                     !  dust chemistry [v/v]
    INTEGER,        INTENT(OUT)   :: RC              ! Success or failure?
!
! !REMARKS:
!  Chemical reactions:
!  ============================================================================
!  (R1) SO2 + O3 + CaCO3 => CaSO4 + O2 + CO2
!                                                                             .
!  (R2) 2(HNO3) + CaCO3 => Ca(NO3)2 + CO2 + H2O
!                                                                             .
! Added sulfate production due to H2SO4 adsorption tdf 2/13/2K9
!  (R3) H2SO4  + CaCO3 => CaSO4 + H2O + CO2
!
! !REVISION HISTORY:
!  28 Mar 2008 - T.D. Fairlie- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER  :: MINDAT = 1.0e-20_fp
!
! !LOCAL VARIABLES:
!
    LOGICAL              :: IT_IS_A_FULLCHEM_SIM
    INTEGER              :: IBIN
    REAL(fp)             :: EQ1
    REAL(fp)             :: HNO3_gas
    REAL(fp)             :: T_SO2, T_HNO3, T_H2SO4, KT1
    REAL(fp)             :: RH2
    REAL(fp)             :: SO2_chem,    DTCHEM
    REAL(fp)             :: SO2_eq,      SO2_new,    SO4d
    REAL(fp)             :: H2SO4_eq,    H2SO4_new
    REAL(fp)             :: HNO3_eq,     HNO3_vv,    HNO3_new
    REAL(fp)             :: HNO3_kg,     HNO3d
    REAL(fp)             :: F_SO2_A,     F_SO2_T,    F_SO2
    REAL(fp)             :: S_FLUX_A,    S_FLUX_T,   S_FLUX(NDSTBIN)
    REAL(fp)             :: F_H2SO4_A,   F_H2SO4_T,  F_H2SO4
    REAL(fp)             :: H_FLUX_A,    H_FLUX_T,   H_FLUX(NDSTBIN)
    REAL(fp)             :: F_HNO3_A,    F_HNO3_T,   F_HNO3
    REAL(fp)             :: N_FLUX_A,    N_FLUX_T,   N_FLUX(NDSTBIN)
    REAL(fp)             :: T_FLUX_A,    TOT_FLUX(NDSTBIN)
    REAL(fp)             :: FALK_SO2,    FALK_HNO3,  FALK_H2SO4
    REAL(fp)             :: ALK_EQ_S (NDSTBIN), ALK_EQ_N (NDSTBIN)
    REAL(fp)             :: ALK_EQ_H (NDSTBIN), TITR_H2SO4(NDSTBIN)
    REAL(fp)             :: TITR_SO2 (NDSTBIN), TITR_HNO3(NDSTBIN)
    REAL(fp)             :: ALK1_vv,     ALK1_kg,    ALK1_eq
    REAL(fp)             :: ALKA_vv,     ALKA_kg,    ALKA_eq
    REAL(fp)             :: END_ALK,     L5A,        L6A,      L7A
    REAL(fp)             :: EQ_BEG,      MW_SO2,     MW_SO4
    REAL(fp)             :: MW_NIT,      MW_HNO3

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)
    REAL(fp), POINTER    :: AD(:,:,:)
    REAL(fp), POINTER    :: AIRVOL(:,:,:)

    !=================================================================
    ! DUST_CHEM begins here!
    !=================================================================

    ! Assume success
    RC                   = GC_SUCCESS

    IT_IS_A_FULLCHEM_SIM = Input_Opt%ITS_A_FULLCHEM_SIM

    ! Set pointers
    Spc                 => State_Chm%Species   ! Chemical species [kg]
    AD                  => State_Met%AD
    AIRVOL              => State_Met%AIRVOL

    ! Set molecular weights locally
    MW_SO2  = State_Chm%SpcData(id_SO2)%Info%MW_g
    MW_SO4  = State_Chm%SpcData(id_SO4)%Info%MW_g
    MW_NIT  = State_Chm%SpcData(id_NIT)%Info%MW_g
    MW_HNO3 = State_Chm%SpcData(id_HNO3)%Info%MW_g

    ! DTCHEM is the chemistry timestep in seconds
    DTCHEM = GET_TS_CHEM()

    ! Convert SO2 [v/v] to  [eq/gridbox]

    !tdf Equivalence defined as moles of a substance * its valence
    !tdf Note 0.064D0 is Mw of SO2 in kg
    !    Hence, SO2_eq =  2 * moles(SO2) / gridbox
    SO2_eq = ( 2.e+0_fp * SO2_cd * AD(I,J,L) ) / &
             ( ( AIRMW / MW_SO2  ) * 0.064e+0_fp )
    SO2_eq = MAX( SO2_eq, MINDAT )

    !tdf 2/13/2K9
    ! Convert H2SO4 [v/v] to  [eq/gridbox]
    ! Note: H2SO4_cd [v/v] provided by H2SO4 production * DTCHEM

    !tdf Equivalence defined as moles of a substance * its valence
    !tdf Note 0.098D0 is Mw of H2SO4 in kg
    !    Hence, H2SO4_eq =  2 * moles(H2SO4) / gridbox

    H2SO4_eq = ( 2.e+0_fp * H2SO4_cd * AD(I,J,L) ) &
               / ( AIRMW / 98.e+0_fp ) / 98.e-3_fp
    H2SO4_eq = MAX( H2SO4_eq, MINDAT )
    !tdf 2/13/2K9

    ! HNO3 mixing ratio
    IF ( IT_IS_A_FULLCHEM_SIM ) THEN

       ! Convert HNO3 [v/v] to [equivalents]

       !tdf Note 28.97/63.0 = Mw(air)/Mw(HNO3)
       !    Hence, HNO3_eq =  1. * moles(HNO3) / gridbox

       HNO3_vv = Spc(id_HNO3)%Conc(I,J,L)
       HNO3_eq = HNO3_vv * AD(I,J,L) / ( AIRMW / 63.e+0_fp ) / 63.e-3_fp

    ELSE

       ! HNO3 is in v/v (from HEMCO)
       HNO3_vv = GLOBAL_HNO3(I,J,L)

       ! Convert HNO3 [v/v] to [equivalents]
       HNO3_eq = HNO3_vv * AD(I,J,L) / ( AIRMW / 63e+0_fp ) / 63.e-3_fp

    ENDIF

    !--------------------------------------------------------
    ! Compute Available SO2 fluxes to dust, S_FLUX (NDSTBIN)
    !--------------------------------------------------------

    S_FLUX_T = 0.0e+0_fp
    F_SO2_T  = 0.0e+0_fp
    RH2      = State_Met%RH(I,J,L) * 1.0e-2_fp

    DO IBIN = 1, NDSTBIN

       KT1      = 0.e+0_fp

       ! Choose a threshold of 18% RH for SO2 uptake flux    tdf 4/22/08
       IF ( RH2 >= 0.18e+0_fp ) THEN
          KT1      = KTS(IBIN)
       ENDIF

       ! Available flux of SO2 to dust aerosols [v/v/timestep]
       L5A      = EXP( -KT1 * DTCHEM )
       F_SO2_A  = SO2_cd * ( 1.e+0_fp - L5A )
       F_SO2_A  = MAX( F_SO2_A, 1.e-32_fp )

       ! Available flux of SO2 converted to [eq/timestep]
       S_FLUX_A = 2.e+0_fp * F_SO2_A * AD(I,J,L) / &
                  ( AIRMW / MW_SO2 ) / 0.064e+0_fp
       S_FLUX (IBIN) = S_FLUX_A

       ! Total available flux of SO2 [v/v/timestep]
       F_SO2_T  = F_SO2_T + F_SO2_A

       ! Total available flux of SO2 [eq/timestep]
       S_FLUX_T = S_FLUX_T + S_FLUX_A

    END DO

    !--------------------------------------------------------
    ! Compute Available H2SO4 fluxes to dust, H_FLUX (NDSTBIN)
    ! tdf 2/13/2K9
    !--------------------------------------------------------

    H_FLUX_T  = 0.0e+0_fp
    F_H2SO4_T = 0.0e+0_fp

    DO IBIN = 1, NDSTBIN

       ! Supplied uptake rates, KTH, for H2SO4 uptake tdf 2/13/2K9
       ! Now KTH is a fraction, so the flux is H2SO4_cd * KTH(IBIN)
       !tdf 08/20/09
       KT1      = KTH(IBIN)

       ! Available flux of H2SO4 to dust aerosols [v/v/timestep]
       !tdf L7A      = EXP( -KT1 * DTCHEM )
       !tdf F_H2SO4_A  = H2SO4_cd * ( 1.e+0_fp - L7A )
       F_H2SO4_A  = H2SO4_cd * KT1

       F_H2SO4_A  = MAX( F_H2SO4_A, 1.e-32_fp )

       ! Available flux of H2SO4 converted to [eq/timestep]
       H_FLUX_A = 2.e+0_fp * F_H2SO4_A * AD(I,J,L) &
                  / ( AIRMW / 98.e+0_fp ) / 0.098e+0_fp
       H_FLUX (IBIN) = H_FLUX_A

       ! Total available flux of H2SO4 [v/v/timestep]
       F_H2SO4_T  = F_H2SO4_T + F_H2SO4_A

       ! Total available flux of H2SO4 [eq/timestep]
       H_FLUX_T = H_FLUX_T + H_FLUX_A

    END DO

    !--------------------------------------------------------
    ! Compute Available HNO3 fluxes to dust, N_FLUX (NDSTBIN)
    !--------------------------------------------------------

    F_HNO3_T = 0.0e+0_fp
    N_FLUX_T = 0.0e+0_fp

    DO IBIN = 1, NDSTBIN

       ! Available flux of HNO3 to dust aerosols [v/v/timestep]
       L6A = EXP( - KTN(IBIN) * DTCHEM )
       F_HNO3_A = HNO3_vv * ( 1.e+0_fp - L6A )
       F_HNO3_A = MAX( F_HNO3_A, 1.0e-32_fp )

       ! Available flux of HNO3 converted to [eq/timestep]
       N_FLUX_A = F_HNO3_A * AD(I,J,L) / ( AIRMW / 63.e+0_fp ) / 0.063e+0_fp
       N_FLUX (IBIN) = N_FLUX_A

       !tdf 3/28/2K8
       ! Accumulate Total available flux of HNO3
       F_HNO3_T = F_HNO3_T + F_HNO3_A ![v/v/timestep]
       N_FLUX_T = N_FLUX_T + N_FLUX_A ![eq/timestep]

    END DO

    !------------------------------------------
    ! Compute Total Available Acid Flux to dust
    !------------------------------------------

    DO IBIN = 1, NDSTBIN

       ! Total acid flux to DUST aerosols [eq/box/timestep] by size bin
       !tdf T_FLUX_A  = S_FLUX (IBIN) + N_FLUX (IBIN)
       !tdf Include sulfuric acid flux                       tdf 2/13/2K9
       T_FLUX_A  = S_FLUX (IBIN) + N_FLUX (IBIN) + H_FLUX (IBIN)
       T_FLUX_A  = MAX( T_FLUX_A, MINDAT )

       ! Total acid flux to DUST aerosols
       TOT_FLUX (IBIN) = T_FLUX_A

    END DO

    !-------------------------------------
    ! Find Total Available Alkalinity [eq]
    !-------------------------------------

    DO IBIN = 1, NDSTBIN

       ALK1_vv = ALK (IBIN)

       !tdf 04/08/08
       ! Convert dust alkalinity from vv to eq., using Mw(Ca) for
       ! Mw (alkalinity). Recall, Equvalents = moles * valency
       ! In this case, we have taken the valency of alkalinity as 2.
       ! Units of ALK1_eq (below) work out to be moles * 2.
       ALK1_eq     = 2.e+0_fp * ALK1_vv * AD(I,J,L) &
                     / ( AIRMW / 40.e+0_fp ) / 40.e-3_fp

       ! total acid flux available; exclude flux from H2SO4, since it is
       ! not limited by dust alkalinity               ! tdf 3/02/2K9
       T_FLUX_A = S_FLUX (IBIN) + N_FLUX (IBIN)
       T_FLUX_A = MAX ( T_FLUX_A, MINDAT )

       ! if the total acid flux available exceeds the available alkalinity
       ! then compute the fraction of the available alkalinity for each acid
       IF ( T_FLUX_A > ALK1_eq ) THEN

          S_FLUX_A  = S_FLUX (IBIN)
          N_FLUX_A  = N_FLUX (IBIN)

          ! Fraction of alkalinity available for each acid
          FALK_SO2  = S_FLUX_A / T_FLUX_A
          FALK_SO2  = MAX( FALK_SO2, MINDAT )
          FALK_HNO3 = N_FLUX_A / T_FLUX_A
          FALK_HNO3 = MAX( FALK_HNO3, MINDAT )

       ELSE

          ! Fraction of alkalinity available for each acid
          FALK_SO2  = 1.0e+0_fp
          FALK_HNO3 = 1.0e+0_fp

       ENDIF

       !tdf Add sulfuric acid flux (not limited by alkalinity)    tdf 2/13/2K9
       FALK_H2SO4 = 1.0e+0_fp

       ! Alkalinity available for S(IV) --> S(VI)
       EQ1             = ALK1_eq * FALK_SO2
       EQ1             = MAX( EQ1, MINDAT )
       ALK_EQ_S (IBIN) = EQ1

       ! Alkalinity available for HNO3 update    tdf 04/07/08
       EQ1             = ALK1_eq * FALK_HNO3
       EQ1             = MAX( EQ1, MINDAT )
       ALK_EQ_N (IBIN) = EQ1

       ! H2SO4 not limited by dust alkalinity     tdf 3/02/2K9

    END DO

    !-------------------
    ! Sulfate production
    !-------------------

    SO2_new       = SO2_eq

    ! Don't produce more SO4 than available ALK or SO2
    DO IBIN = 1, NDSTBIN

       S_FLUX_A = S_FLUX (IBIN)
       EQ1      = ALK_EQ_S (IBIN)
       SO4d     = MIN( S_FLUX_A, EQ1, SO2_new )
       SO4d     = MAX( SO4d, MINDAT )

       ! Update SO2 concentration [eq/box]
       SO2_new    = SO2_new - SO4d
       SO2_new    = MAX( SO2_new, MINDAT )

       ! Alkalinity titrated by S(IV) --> S(VI) [eq]
       TITR_SO2 (IBIN) =  SO4d

       !SO4d produced converted from [eq/timestep] to [v/v/timestep]
       PSO4d (IBIN) = SO4d * 0.096e+0_fp * ( AIRMW / MW_SO4 ) / &
                      AD(I,J,L) / 2.0e+0_fp

       ! tdf
       !if (I .eq. 1 .and. J .eq. 63 .and. L .eq. 6) then
       !   print *,' CHEM_SO4: SO4 production, SO2'
       !   write (6,30) IBIN, KTS(IBIN)
       !   print *,' IBIN,EQ1,S_FLUX_A,SO2_new,SO4d,PSO4d(IBIN)'
       !   write (6,31) IBIN,EQ1,S_FLUX_A,SO2_new,SO4d,PSO4d(IBIN)
       !30 format (' IBIN, KTS(IBIN) ',I4,E12.3)
       !31 format (' ',I4,5E12.3)
       !endif

    END DO

    !Modified SO2 [eq] converted back to [v/v]
    SO2_gas       = SO2_new * 0.064e+0_fp * ( AIRMW / MW_SO2 ) / &
                    AD(I,J,L) / 2.0e+0_fp
    SO2_gas       = MAX( SO2_gas, MINDAT )

    !------------------------------------------------
    ! Additional sulfate production from H2SO4 uptake
    !------------------------------------------------

    H2SO4_new   = H2SO4_eq

    ! Don't produce more SO4 than available H2SO4
    ! Uptake not limited by alkalinity
    DO IBIN = 1, NDSTBIN

       H_FLUX_A = H_FLUX (IBIN)

       ! H2SO4 uptake not limited by dust alkalinity, EQ1
       SO4d     = MIN( H_FLUX_A, H2SO4_new )
       SO4d     = MAX( SO4d, MINDAT )

       ! Update H2SO4 concentration [eq/box]
       H2SO4_new    = H2SO4_new - SO4d
       H2SO4_new    = MAX( H2SO4_new, MINDAT )

       ! Alkalinity titrated by H2SO4 uptake [eq]
       TITR_H2SO4 (IBIN) =  SO4d

       !SO4d produced converted from [eq/timestep] to [v/v/timestep]
       PH2SO4d (IBIN) = SO4d * 0.096e+0_fp * ( AIRMW / MW_SO4 ) / &
                        AD(I,J,L) / 2.0e+0_fp

       ! tdf
       !if (I .eq. 1 .and. J .eq. 63 .and. L .eq. 6) then
       !   print *,' CHEM_SO4: SO4 production, H2SO4'
       !   write (6,40) IBIN, KTH(IBIN)
       !   print *,' IBIN,H_FLUX_A,H2SO4_new,SO4d,PH2SO4d(IBIN)'
       !   write (6,41) IBIN,H_FLUX_A,H2SO4_new,SO4d,PH2SO4d(IBIN)
       !40 format (' IBIN, KTH(IBIN) ',I4,E12.3)
       !41 format (' ',I4,4E12.3)
       !endif

    END DO

    !Modified H2SO4 [eq] converted back to [v/v]
    H2SO4_gas       = H2SO4_new * 0.098e+0_fp * ( AIRMW &
                      / 98.e+0_fp ) / AD(I,J,L) / 2.0e+0_fp ! Hard-coded MW
    H2SO4_gas       = MAX( H2SO4_gas, MINDAT )

    !-------------------
    ! Nitrate production
    !-------------------

    HNO3_new    = HNO3_eq

    ! Alkalinity titrated by HNO3 in dust
    DO IBIN = 1, NDSTBIN

       N_FLUX_A = N_FLUX (IBIN)
       EQ1      = ALK_EQ_N (IBIN)
       HNO3d    = MIN( N_FLUX_A, EQ1, HNO3_new )
       HNO3d    = MAX( HNO3d, MINDAT )

       ! Update HNO3 concentration [eq/box]
       HNO3_new = HNO3_new - HNO3d
       HNO3_new = MAX( HNO3_new, MINDAT )

       ! Alkalinity titrated by HNO3 [eq]
       TITR_HNO3 (IBIN) = HNO3d

       ! NIT produced converted from [eq/timestep] to [v/v/timestep]
       PNITd (IBIN) = HNO3d * 0.063e+0_fp * ( AIRMW / MW_NIT ) / AD(I,J,L)

       ! tdf
       !if (I .eq. 1 .and. J .eq. 63 .and. L .eq. 6) then
       !   print *,' CHEM_SO4: NIT production, HNO3'
       !   write (6,50) IBIN, KTN(IBIN)
       !   print *,' IBIN,EQ1,N_FLUX_A,HNO3_new,HNO3d,PNITd(IBIN)'
       !   write (6,51) IBIN,EQ1,N_FLUX_A,HNO3_new,HNO3d,PNITd(IBIN)
       !50 format (' IBIN, KTN(IBIN) ',I4,E12.3)
       !51 format (' ',I4,5E12.3)
       !endif

    END DO

    !Modified HNO3 [eq/timestep] converted back to [v/v/timestep]
    HNO3_gas      = HNO3_new * 0.063e+0_fp * ( AIRMW / MW_HNO3 ) / AD(I,J,L)
    HNO3_gas      = MAX( HNO3_gas, MINDAT )

    ! HNO3 [v/v]
    IF ( id_HNO3 > 0 ) THEN
       Spc(id_HNO3)%Conc(I,J,L) = MAX( HNO3_gas, MINDAT )
    ENDIF

    DO IBIN = 1, NDSTBIN

       ALK1_vv     = ALK (IBIN)
       ALK1_eq     = 2.e+0_fp * ALK1_vv * AD(I,J,L) &
                     / ( AIRMW / 40.e+0_fp ) / 40.e-3_fp ! Hard-coded MW
       T_SO2       = TITR_SO2 (IBIN)
       T_HNO3      = TITR_HNO3 (IBIN)

       !tdf Include alkalinity titrated by sulfuric acid flux 2/13/2K9
       T_H2SO4     = TITR_H2SO4 (IBIN)

       ! tdf
       ! if (I .eq. 1 .and. J .eq. 63 .and. L .eq. 6) then
       !    print *,' CHEM_DUST: Titrate Alkalinity'
       !    print *,' IBIN,  ALK1_eq,  T_SO2,   T_HNO3,   T_H2SO4'
       !    write (6,61) IBIN,ALK1_eq,T_SO2,T_HNO3,T_H2SO4
       !61 format (' ',I4,4E12.3)
       !endif

       ! Titrate DUST alkalinity  [eq]
       ALKA_eq     = ALK1_eq - ( T_SO2 + T_HNO3 + T_H2SO4 )
       ALKA_eq     = MAX( ALKA_eq, MINDAT )

       ! Note:  Although we don't let the alkalinity go negative,
       ! the uptake of H2SO4 can continue when the alkalinity is
       ! fully titrated.                            ! tdf 3/02/2K9

       ! Return remaining DUST Alkalinity [v/v]
       ALKA_vv     = ALKA_eq / AD(I,J,L) * ( AIRMW / 40.e+0_fp ) * &
                     40e-3_fp / 2.e+0_fp ! Hard-coded MW
       ALKA (IBIN) = ALKA_vv

    END DO

    ! Update dust alkalinity
    ! NB Hardwired for 4 size bins                    tdf 04/08/08

    Spc(id_DAL1)%Conc(I,J,L) = MAX( ALKA(1), MINDAT )
    Spc(id_DAL2)%Conc(I,J,L) = MAX( ALKA(2), MINDAT )
    Spc(id_DAL3)%Conc(I,J,L) = MAX( ALKA(3), MINDAT )
    Spc(id_DAL4)%Conc(I,J,L) = MAX( ALKA(4), MINDAT )

    ! Free pointers
    Spc    => NULL()
    AD     => NULL()
    AIRVOL => NULL()

  END SUBROUTINE DUST_CHEM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_hplus
!
! !DESCRIPTION: Subroutine GET\_HPLUS computes H+ concentrations in cloud
!  liquid water for pH dependent cloud chemistry. (bec, 4/11/11)
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE GET_HPLUS( SO4nss, HMsc, TNH3, TNO3, SO2, CL, TNA, TDCA, TFA, &
                          TAA,  T, PRES, LWC,  iHPLUS, HPLUS )
!
! !USES:
!
    USE ERROR_MOD,       ONLY : IT_IS_NAN, GEOS_CHEM_STOP
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN)    :: SO4nss ! Total nss sulfate mixing ratio [M]
    REAL(fp),  INTENT(IN)    :: HMSc ! Total HMS mixing ratio [M]
    REAL(fp),  INTENT(IN)    :: TNO3   ! Total nitrate (gas+particulate) mixing
                                       ! ratio [v/v]
    REAL(fp),  INTENT(IN)    :: TNH3   ! NH3 mixing ratio [v/v]
    REAL(fp),  INTENT(IN)    :: SO2    ! SO2 mixing ratio [v/v]
    REAL(fp),  INTENT(IN)    :: CL     ! Total chloride (gas+particulate) mixing
    REAL(fp),  INTENT(IN)    :: TNA    ! Sodium (particulate) [v/v]
    REAL(fp),  INTENT(IN)    :: TDCA   ! Total Ca2+ and Mg2+ mixing ratio [M] ! jmm 12/3/18
    REAL(fp),  INTENT(IN)    :: TAA    ! Acetic acid mixing ratio [v/v] ! jmm 12/3/18
    REAL(fp),  INTENT(IN)    :: TFA    ! Formic acid mixing ratio [v/v] ! jmm 12/3/18
    REAL(fp),  INTENT(IN)    :: T      ! Temperature [K]
    REAL(fp),  INTENT(IN)    :: PRES   ! Dry air partial ressure [atm]
    REAL(fp),  INTENT(IN)    :: LWC    ! Cloud liquid water content [m3/m3]
    REAL(fp),  INTENT(IN)    :: iHPLUS ! Initial [H+] [M]
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),  INTENT(OUT)   :: HPLUS  ! Calculated [H+] [M]
! !REMARKS:
!  Calculation:
!  ============================================================================
!  Solve the following electroneutrality equation:
!  [H+] = 2[SO4--] + [Cl-] + [OH-] + [HCO3-] + 2[CO3--] + [HSO3-] + 2[SO3--] +!
!          [NO3-] + [HCOO-] + [CH3COO-] +[HMS] - [Na] - 2[Ca] - [NH4]
!  Uses Newton's method to solve the equation:
!     x_1 = x_0 -f(x_0)/f'(x_0)
!     iterate until converge
!
!  Let concentrations of [HCO3], [CO3], [HSO3], [SO3], [NO3] and [NH4] evolve
!  according to Henry's law equilibrium.
!
!  To add new species:
!    - Add species not affected by HPLUS to the "D' term
!    - Add species that disassociate once using the kHNO3 and dkHNO3
!    functions
!      as a template
!    - Add species that disassociate twice using the kSO21 and dkSO21
!    functions
!      as a template for the single charged ion and kSO22 and dkSO22
!      functions for
!      the double charged ion

!  Assume [S(VI)] = [SO4]nss (this applies for pH > 3)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Water dissociation constants
    REAL(fp),  PARAMETER   :: Kw   = 1.0e-14_fp
    REAL(fp),  PARAMETER   :: DhrKw = -6710.e+0_fp
    REAL(fp),  PARAMETER   :: MINVAL = 0.01
!
! !LOCAL VARIABLES:
!
    REAL(fp)               :: D, Kw_T, ipH, newpH, nHPLUS
    REAL(fp)               :: fHCO3, fCO3
    REAL(fp)               :: fHSO3, fSO3
    REAL(fp)               :: fHNO3, fNH4, fHCl
    REAL(fp)               :: dHCO3, dCO3
    REAL(fp)               :: dHSO3, dSO3
    REAL(fp)               :: dHNO3, dNH4, dHCl
    REAL(fp)               :: fAA, fFA, dAA, dFA
    REAL(fp)               :: f, df, nnHPLUS, fCa, dCa
    INTEGER                :: count

    !=================================================================
    ! GET_HPLUS begins here!
    !=================================================================

    ! Initial pH guess
    ipH = -log10(iHPLUS)

    ! Non-volatile aerosol concentration [M]
    ! For now sulfate is the only non-volatile species
#ifdef LUO_WETDEP
    D = ( 1.5_fp * SO4nss ) - TNA - ( 2.0_fp * TDCA ) + HMSc
#else
    D = ( 2.0_fp * SO4nss ) - TNA - ( 2.0_fp * TDCA ) + HMSc
#endif

    ! Temperature dependent water equilibrium constant
    Kw_T = Kw*exp(DhrKw*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    ! Initialize
    newpH   = 0.0
    COUNT = 0

    DO WHILE ( ABS(ipH-newpH) .gt. MINVAL )

       COUNT = COUNT+1

       IF ( COUNT .EQ. 1 ) THEN
          ipH = ipH
       ELSE
          ipH = newpH
       ENDIF

       nHPLUS = 10.e+0_fp**(-ipH)

       ! Get f(x) terms
       fHCO3  = kCO21 ( PRES, T, LWC, nHPLUS )

       fCO3 = kCO22 ( PRES, T, LWC, nHPLUS )

       fHSO3  = kSO21 ( PRES, T, LWC, nHPLUS, SO2 )

       fSO3 = kSO22 ( PRES, T, LWC, nHPLUS, SO2 )

       fHNO3 = kHNO3 ( PRES, T, LWC, nHPLUS, TNO3 )

       fNH4  = kNH3  ( PRES, T, LWC, nHPLUS, TNH3, Kw_T )

       ! include HCl in cloud pH calculations, xnw 10/17/17
       fHCl  = kHCl  ( PRES, T, LWC, nHPLUS, CL  )

       fFA   = kFA   ( PRES, T, LWC, nHPLUS, TFA ) ! jmm 12/3/18

       fAA   = kAA   ( PRES, T, LWC, nHPLUS, TAA ) ! jmm 12/3/18

       ! Get f'(x) terms
       dHCO3  = dkCO21 ( PRES, T, LWC, nHPLUS )

       dCO3 = dkCO22 ( PRES, T, LWC, nHPLUS )

       dHSO3  = dkSO21 ( PRES, T, LWC, nHPLUS, SO2 )

       dSO3 = dkSO22 ( PRES, T, LWC, nHPLUS, SO2 )

       dHNO3 = dkHNO3 ( PRES, T, LWC, nHPLUS, TNO3 )

       dNH4  = dkNH3  ( PRES, T, LWC, nHPLUS, TNH3, Kw_T )

       dHCl = dkHCl ( PRES, T, LWC, nHPLUS, CL )

       dFA   = dkFA   ( PRES, T, LWC, nHPLUS, TFA ) ! jmm 12/3/18

       dAA   = dkAA   ( PRES, T, LWC, nHPLUS, TAA ) ! jmm 12/3/18
       ! Calculate [Ca2+] in equilibrium with CaCO3(s)
       CALL CaCO3_PRECIP ( PRES, T, nHPLUS, fCa, dCa )

       ! if [Ca2+] in equilibrium with CacO3(s) is greater than total [Ca2+]
       ! then all Ca is dissolved else [Ca2+] varies with [H+]
       IF ( fCa .ge. TDCA ) THEN
          ! Non-volatile aerosol concentration [M]
#ifdef LUO_WETDEP
          D = (1.5e+0_fp*SO4nss) - (TNA+2.e+0_fp*TDCA)
#else
          D = (2.e+0_fp*SO4nss) - (TNA+2.e+0_fp*TDCA)
#endif

          ! Define f(x)
          f = D - nHPLUS + Kw/nHPLUS + fHCO3 + 2.e+0_fp * &
               fCO3 + fHSO3 + 2.e+0_fp * fSO3 + fHNO3 - fNH4 + &
               fHCl + fFA + fAA

          ! Define f'(x)
          df = - 1.d0 - Kw/nHPLUS/nHPLUS + dHCO3 + 2.e+0_fp * &
               dCO3 + dHSO3 + 2.e+0_fp * dSO3 + dHNO3 - dNH4 + &
               dHCl + dFA + dAA

       ELSE
          ! Non-volatile aerosol concentration [M]
#ifdef LUO_WETDEP
          D = (1.5e+0_fp * SO4nss) - TNA
#else
          D = (2.e+0_fp * SO4nss) - TNA
#endif

          ! Define f(x)
          f = D - nHPLUS + Kw/nHPLUS + fHCO3 + 2.e+0_fp * fCO3 + &
               fHSO3 + 2.e+0_fp * fSO3 + fHNO3 - fNH4 + &
               fHCl + fFA + fAA - 2.e+0_fp * fCa
          ! Define f'(x)
          df = - 1.d0 - Kw/nHPLUS/nHPLUS + dHCO3 + 2.e+0_fp * dCO3 + &
               dHSO3 + 2.e+0_fp * dSO3 + dHNO3 - dNH4 + &
               dHCl + dFA + dAA - 2.e+0_fp * dCa
       ENDIF

       ! Apply Newton's method
       nnHPLUS = nHPLUS - f/df

       ! Set minimum [H+] = 1.d-14 (pH = 14)
       nnHPLUS = MAX(nnHPLUS,1.0e-14_fp)

       ! Set maximum [H+] = 1.d-1 (pH = 1)
       nnHPLUS = MIN(nnHPLUS,1.0e-1_fp)

       ! If solution does not converge after 50 iterations
       ! average last 2 pH calculations
       IF (count > 50) THEN
          newpH = ((-log10(nnHPLUS)) + (-log10(nHPLUS))) / 2.0e+0_fp

          IF (IT_IS_NAN( newpH )) THEN
             write(6,*) 'newpH = ', newpH
             write(6,*) 'nnHPLUS = ', nnHPLUS
             write(6,*) 'nHPLUS = ', nHPLUS
             CALL GEOS_CHEM_STOP
          ENDIF

          EXIT
       ELSE
          newpH = -log10(nnHPLUS)

          IF (IT_IS_NAN( newpH )) THEN
             write(6,*) 'newpH = ', newpH
             write(6,*) 'nnHPLUS = ', nnHPLUS
             CALL GEOS_CHEM_STOP
          ENDIF

       ENDIF

    ENDDO

    HPLUS = 10.0e+0_fp**(-newpH)

  END SUBROUTINE GET_HPLUS

!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kCO21
!
! !DESCRIPTION: Function kCO21
!\\
!\\
! !INTERFACE:
!
  FUNCTION kCO21 ( P, T, LWC, HPLUS ) RESULT ( KCO2p )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KCO2p, KCO2p2
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! CO2 dissociation constants
    REAL(fp),  PARAMETER  :: Kc1 = 4.3e-7
    REAL(fp),  PARAMETER  :: Kc2 =4.68e-11
    REAL(fp),  PARAMETER  :: DhrKc1 = -1000.
    REAL(fp),  PARAMETER  :: DhrKc2 = -1760.
    REAL(fp),  PARAMETER  :: Hco2 = 3.4e-2
    REAL(fp),  PARAMETER  :: Dhco2 = 2.44e+3_fp
    ! CO2 concentration [v/v]
    REAL(fp),  PARAMETER  :: CO2 = 390.0e-6_fp
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: Hco2_T, Kc1_T, Kc2_T
    REAL(fp)              :: Hco2eff, xCO2, pCO2

    !=================================================================
    ! kCO21 begins here!
    !=================================================================

    !CO2 dissolution constants
    Hco2_T  = Hco2*exp(Dhco2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Kc1_T   = Kc1*exp(DhrKc1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Kc2_T   = Kc2*exp(DhrKc2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    !CO2 dissolution
    Hco2eff = Hco2_T*(1.e+0_fp+(Kc1_T/HPLUS)+((Kc1_T*Kc2_T)/(HPLUS*HPLUS)))
    xCO2    = 1.e+0_fp / ( 1.e+0_fp + ( Hco2eff * 0.08205e+0_fp * T * LWC ) )
    pCO2    = CO2 * P * xCO2

    KCO2p  = Hco2_T / HPLUS * Kc1_T * pCO2

  END FUNCTION kCO21
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkCO21
!
! !DESCRIPTION: Function dkCO21
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkCO21 ( P, T, LWC, HPLUS ) RESULT ( KCO2p )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KCO2p, KCO2p2
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  22 Mar 2017 - M. Sulprizio- Dhco2 value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (V. Shah)
!  15 Feb 2019 - J. Moch     - updated function to make output
!  derivative of [HCO3-]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! CO2 dissociation constants
      REAL(fp),  PARAMETER  :: Kc1 = 4.3e-7
      REAL(fp),  PARAMETER  :: Kc2 =4.68e-11
      REAL(fp),  PARAMETER  :: DhrKc1 = -1000.
      REAL(fp),  PARAMETER  :: DhrKc2 = -1760.
      REAL(fp),  PARAMETER  :: Hco2 = 3.4e-2
      REAL(fp),  PARAMETER  :: Dhco2 = 2.44e+3_fp
      ! CO2 concentration [v/v]
      REAL(fp),  PARAMETER  :: CO2 = 390.0e-6_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hco2_T, Kc1_T, Kc2_T

! !REMARKS:

      !=================================================================
      ! dkCO21 begins here!
      !=================================================================

      !CO2 dissolution constants
      Hco2_T = Hco2*exp(Dhco2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc1_T = Kc1*exp(DhrKc1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc2_T = Kc2*exp(DhrKc2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      !CO2 dissolution

      KCO2p  = Kc1_T * Hco2_T * CO2 * P * ( Kc1_T * Kc2_T * Hco2_T * &
          0.08205e+0_fp * T * LWC - Hco2_T * 0.08205e+0_fp * T *     &
          LWC * HPLUS * HPLUS - HPLUS * HPLUS) / (Kc1_T * Kc2_T *    &
          Hco2_T * 0.08205e+0_fp * T * LWC + Kc1_T * Hco2_T *        &
          0.08205e+0_fp * T * LWC * HPLUS + Hco2_T * 0.08205e+0_fp * &
          T * LWC * HPLUS * HPLUS + HPLUS * HPLUS)**2

      END FUNCTION dkCO21
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kCO22
!
! !DESCRIPTION: Function kCO22
!\\
!\\
! !INTERFACE:
!
  FUNCTION kCO22 ( P, T, LWC, HPLUS ) RESULT ( KCO2p2 )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KCO2p, KCO2p2
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! CO2 dissociation constants
    REAL(fp),  PARAMETER  :: Kc1 = 4.3e-7
    REAL(fp),  PARAMETER  :: Kc2 =4.68e-11
    REAL(fp),  PARAMETER  :: DhrKc1 = -1000.
    REAL(fp),  PARAMETER  :: DhrKc2 = -1760.
    REAL(fp),  PARAMETER  :: Hco2 = 3.4e-2
    REAL(fp),  PARAMETER  :: Dhco2 = 2.44e+3_fp
    ! CO2 concentration [v/v]
    REAL(fp),  PARAMETER  :: CO2 = 390.0e-6_fp
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: Hco2_T, Kc1_T, Kc2_T
    REAL(fp)              :: Hco2eff, xCO2, pCO2

    !=================================================================
    ! kCO22 begins here!
    !=================================================================

    !CO2 dissolution constants
    Hco2_T  = Hco2*exp(Dhco2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Kc1_T   = Kc1*exp(DhrKc1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Kc2_T   = Kc2*exp(DhrKc2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    !CO2 dissolution
    Hco2eff = Hco2_T*(1.e+0_fp+(Kc1_T/HPLUS)+((Kc1_T*Kc2_T)/(HPLUS*HPLUS)))
    xCO2    = 1.e+0_fp / ( 1.e+0_fp  + ( Hco2eff * 0.08205e+0_fp * T * LWC ) )
    pCO2    = CO2 * P * xCO2

    KCO2p2 = Kc1_T * Kc2_T * Hco2_T / HPLUS / HPLUS * pCO2

  END FUNCTION kCO22
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkCO22
!
! !DESCRIPTION: Function dkCO22
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkCO22 ( P, T, LWC, HPLUS ) RESULT ( KCO2p2 )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KCO2p, KCO2p2
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  22 Mar 2017 - M. Sulprizio- Dhco2 value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (V. Shah)
!  15 Feb 2019 - J. Moch     - updated function to make output deriviate
!  of [CO3--]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! CO2 dissociation constants
      REAL(fp),  PARAMETER  :: Kc1 = 4.3e-7
      REAL(fp),  PARAMETER  :: Kc2 =4.68e-11
      REAL(fp),  PARAMETER  :: DhrKc1 = -1000.
      REAL(fp),  PARAMETER  :: DhrKc2 = -1760.
      REAL(fp),  PARAMETER  :: Hco2 = 3.4e-2
      REAL(fp),  PARAMETER  :: Dhco2 = 2.44e+3_fp
      ! CO2 concentration [v/v]
      REAL(fp),  PARAMETER  :: CO2 = 390.0e-6_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hco2_T, Kc1_T, Kc2_T

      !=================================================================
      ! dkCO22 begins here!
      !=================================================================

      !CO2 dissolution constants
      Hco2_T = Hco2*exp(Dhco2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc1_T = Kc1*exp(DhrKc1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc2_T = Kc2*exp(DhrKc2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      !CO2 dissolution

      KCO2p2 = -1.e+0_fp * Kc1_T * Kc2_T * Hco2_T * CO2 * P * ( Kc1_T * &
           Hco2_T * 0.08205e+0_fp * T * LWC + 2.0e+0_fp * Hco2_T *      &
           0.08205e+0_fp * T * LWC * HPLUS + 2.0e+0_fp * HPLUS ) /      &
           ( Kc1_T * Kc2_T * Hco2_T * 0.08205e+0_fp * T * LWC +         &
           Kc1_T * Hco2_T * 0.08205e+0_fp * T * LWC * HPLUS +           &
           Hco2_T *0.08205e+0_fp * T * LWC * HPLUS * HPLUS +            &
           HPLUS * HPLUS )**2

      END FUNCTION dkCO22
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kSO21
!
! !DESCRIPTION: Function kSO21
!\\
!\\
! !INTERFACE:
!
  FUNCTION kSO21 ( P, T, LWC, HPLUS, SO2 ) RESULT ( KSO2p )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, SO2
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KSO2p, KSO2p2
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! SO2 dissociation constants
    REAL(fp),  PARAMETER  :: Ks1 = 1.3e-2
    REAL(fp),  PARAMETER  :: Ks2 = 6.6e-8
    REAL(fp),  PARAMETER  :: Hso2 = 1.23
    REAL(fp),  PARAMETER  :: Dhso2 = 3.14e+3_fp
    REAL(fp),  PARAMETER  :: DhrKso21 = 1960.
    REAL(fp),  PARAMETER  :: DhrKso22 = 1500.
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: Hso2_T, Ks1_T, Ks2_T
    REAL(fp)              :: Hso2eff, xSO2, pSO2

    !=================================================================
    ! kSO21 begins here!
    !=================================================================

    ! SO2 dissolution constants
    Hso2_T  = Hso2*exp(Dhso2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Ks1_T   = Ks1*exp(DhrKso21*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Ks2_T   = Ks2*exp(DhrKso22*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    ! SO2 dissolution
    Hso2eff = Hso2_T*(1.e+0_fp+(Ks1_T/HPLUS)+((Ks1_T*Ks2_T)/(HPLUS*HPLUS)))
    xSO2    = 1.e+0_fp / ( 1.e+0_fp  + ( Hso2eff * 0.08205e+0_fp * T * LWC ) )
    pSO2    = SO2 * P * xSO2

    KSO2p   = Hso2_T * Ks1_T * pSO2 / HPLUS

  END FUNCTION kSO21
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkSO21
!
! !DESCRIPTION: Function dkSO21
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkSO21 ( P, T, LWC, HPLUS, SO2 ) RESULT ( KSO2p )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, SO2
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KSO2p, KSO2p2
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  22 Mar 2017 - M. Sulprizio- Dhso2 value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (V. Shah)
!  15 Feb 2019 - J. Moch     - updated function to make output
!  derivative of [HSO3-]

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! SO2 dissociation constants
      REAL(fp),  PARAMETER  :: Ks1 = 1.3e-2
      REAL(fp),  PARAMETER  :: Ks2 = 6.6e-8
      REAL(fp),  PARAMETER  :: Hso2 = 1.23
      REAL(fp),  PARAMETER  :: Dhso2 = 3.14e+3_fp
      REAL(fp),  PARAMETER  :: DhrKso21 = 1960.
      REAL(fp),  PARAMETER  :: DhrKso22 = 1500.
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hso2_T, Ks1_T, Ks2_T

      !=================================================================
      ! dkSO21 begins here!
      !=================================================================



      ! SO2 dissolution constants
      Hso2_T = Hso2*exp(Dhso2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Ks1_T = Ks1*exp(DhrKso21*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Ks2_T = Ks2*exp(DhrKso22*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))


      KSO2p  = Ks1_T * Hso2_T * SO2 * P * ( Ks1_T * Ks2_T * Hso2_T *  &
           0.08205e+0_fp * T * LWC - Hso2_T * 0.08205e+0_fp * T *     &
           LWC * HPLUS * HPLUS - HPLUS * HPLUS) / (Ks1_T * Ks2_T *    &
           Hso2_T * 0.08205e+0_fp * T * LWC + Ks1_T * Hso2_T *        &
           0.08205e+0_fp * T * LWC * HPLUS + Hso2_T * 0.08205e+0_fp * &
           T * LWC * HPLUS * HPLUS + HPLUS * HPLUS)**2

      END FUNCTION dkSO21
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kSO22
!
! !DESCRIPTION: Function kSO22
!\\
!\\
! !INTERFACE:
!
  FUNCTION kSO22 ( P, T, LWC, HPLUS, SO2 ) RESULT ( KSO2p2 )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, SO2
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KSO2p, KSO2p2
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! SO2 dissociation constants
    REAL(fp),  PARAMETER  :: Ks1 = 1.3e-2
    REAL(fp),  PARAMETER  :: Ks2 = 6.6e-8
    REAL(fp),  PARAMETER  :: Hso2 = 1.23
    REAL(fp),  PARAMETER  :: Dhso2 = 3.14e+3_fp
    REAL(fp),  PARAMETER  :: DhrKso21 = 1960.
    REAL(fp),  PARAMETER  :: DhrKso22 = 1500.
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: Hso2_T, Ks1_T, Ks2_T
    REAL(fp)              :: Hso2eff, xSO2, pSO2

    !=================================================================
    ! kSO22 begins here!
    !=================================================================

    ! SO2 dissolution constants
    Hso2_T  = Hso2*exp(Dhso2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Ks1_T   = Ks1 *exp(DhrKso21*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Ks2_T   = Ks2 *exp(DhrKso22*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    !SO2 dissolution
    Hso2eff = Hso2_T*(1.e+0_fp+(Ks1_T/HPLUS)+((Ks1_T*Ks2_T)/(HPLUS*HPLUS)))
    xSO2    = 1.e+0_fp / ( 1.e+0_fp + ( Hso2eff * 0.08205e+0_fp * T * LWC ) )
    pSO2    = SO2 * P * xSO2

    KSO2p2 = Ks1_T * Ks2_T * Hso2_T / HPLUS / HPLUS * pSO2

  END FUNCTION kSO22
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkSO22
!
! !DESCRIPTION: Function dkSO22
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkSO22 ( P, T, LWC, HPLUS, SO2 ) RESULT ( KSO2p2 )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, SO2
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KSO2p, KSO2p2
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  22 Mar 2017 - M. Sulprizio- Dhso2 value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (V. Shah)
!  15 Feb 2019 - J. Moch     - updated function to make output
!  derivative [SO3--]

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! SO2 dissociation constants
      REAL(fp),  PARAMETER  :: Ks1 = 1.3e-2
      REAL(fp),  PARAMETER  :: Ks2 = 6.6e-8
      REAL(fp),  PARAMETER  :: Hso2 = 1.23
      REAL(fp),  PARAMETER  :: Dhso2 = 3.14e+3_fp
      REAL(fp),  PARAMETER  :: DhrKso21 = 1960.
      REAL(fp),  PARAMETER  :: DhrKso22 = 1500.
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hso2_T, Ks1_T, Ks2_T

      !=================================================================
      ! dkSO22 begins here!
      !=================================================================
      ! SO2 dissolution constants
      Hso2_T = Hso2*exp(Dhso2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Ks1_T  = Ks1 *exp(DhrKso21*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Ks2_T  = Ks2 *exp(DhrKso22*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      KSO2p2 = -1.e+0_fp * Ks1_T * Ks2_T * Hso2_T * SO2 * P * ( Ks1_T * &
           Hso2_T * 0.08205e+0_fp * T * LWC + 2.0e+0_fp * Hso2_T *      &
           0.08205e+0_fp * T * LWC * HPLUS + 2.0e+0_fp * HPLUS ) /      &
           ( Ks1_T * Ks2_T * Hso2_T * 0.08205e+0_fp * T * LWC +         &
           Ks1_T * Hso2_T * 0.08205e+0_fp * T * LWC * HPLUS +           &
           Hso2_T *0.08205e+0_fp * T * LWC * HPLUS * HPLUS +            &
           HPLUS * HPLUS )**2

      END FUNCTION dkSO22
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kHNO3
!
! !DESCRIPTION: Function kNO3
!\\
!\\
! !INTERFACE:
!
  FUNCTION kHNO3 ( P, T, LWC, HPLUS, HNO3 ) RESULT ( KHNO3p )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, HNO3
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KHNO3p
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! HNO3 dissociation constants
    REAL(fp),  PARAMETER  :: Kn1 = 15.4
    REAL(fp),  PARAMETER  :: Hhno3 = 2.1e5
    REAL(fp),  PARAMETER  :: Dhhno3 = 0.
    REAL(fp),  PARAMETER  :: DhrKn1 = 8700.
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: Hhno3_T, Kn1_T
    REAL(fp)              :: Hhno3eff, xHNO3, pHNO3

    !=================================================================
    ! kHNO3 begins here!
    !=================================================================

    ! HNO3 dissolution constants
    Hhno3_T  = Hhno3*exp(Dhhno3*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Kn1_T    = Kn1*exp(DhrKn1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    ! HNO3 dissolution
    ! The original Hhno3eff expression is valid for 298K (Seinfeld and Pandis
    ! 2006, pp 299-301), and Kn1 has a strong temperature dependence. The
    ! fix follows Eq. 7.59 of Seinfeld and Pandis (2006, pp 301).
    !Hhno3eff = 3.2e6/HPLUS
    Hhno3eff = Hhno3_T*(1.0e+0_fp+(Kn1_T/HPLUS))
    xHNO3    = 1.e+0_fp / ( 1.e+0_fp + ( Hhno3eff * 0.08205e+0_fp * T * LWC ) )
    pHNO3    = HNO3 * P * xHNO3

    kHNO3p = Hhno3_T * Kn1_T * pHNO3 / HPLUS

  END FUNCTION kHNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkHNO3
!
! !DESCRIPTION: Function dkNO3
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkHNO3 ( P, T, LWC, HPLUS, HNO3 ) RESULT ( KHNO3p )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, HNO3
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KHNO3p
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  22 Mar 2017 - M. Sulprizio- Add fix for Hhno3eff from V. Shah
!  15 Feb 2019 - J. Moch     - updated function to make output
!  derivative of [HNO3-]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! HNO3 dissociation constants
      REAL(fp),  PARAMETER  :: Kn1 = 15.4
      REAL(fp),  PARAMETER  :: Hhno3 = 2.1e5
      REAL(fp),  PARAMETER  :: Dhhno3 = 0.
      REAL(fp),  PARAMETER  :: DhrKn1 = 8700.
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hhno3_T, Kn1_T

      !=================================================================
      ! dkHNO3 begins here!
      !=================================================================

      ! HNO3 dissolution constants
      Hhno3_T = Hhno3*exp(Dhhno3*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kn1_T = Kn1*exp(DhrKn1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      ! HNO3 dissolution
      ! The original Hhno3eff expression is valid for 298K (Seinfeld and
      ! Pandis
      ! 2006, pp 299-301), and Kn1 has a strong temperature dependence.
      ! The
      ! fix follows Eq. 7.59 of Seinfeld and Pandis (2006, pp 301).

      kHNO3p = -1.0e+0_fp * Kn1_T * Hhno3_T * HNO3 * P * &
          ( 1.0e+0_fp + Hhno3_T * 0.08205e+0_fp * T * LWC ) / &
          ( Kn1_T * Hhno3_T * 0.08205e+0_fp * T * LWC + &
          Hhno3_T * 0.08205e+0_fp * T * LWC * HPLUS + &
          HPLUS )**2

      END FUNCTION dkHNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kHCl
!
! !DESCRIPTION: Function kHCl
!\\
!\\
! !INTERFACE:
!
  FUNCTION kHCl ( P, T, LWC, HPLUS, Cl ) RESULT ( KHClp )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, Cl
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KHClp
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! HNO3 dissociation constants
    REAL(fp),  PARAMETER  :: Kcl = 1.74e+6_fp
    REAL(fp),  PARAMETER  :: Hcl = 1.5e+3_fp
    REAL(fp),  PARAMETER  :: Dhcl = 2.3e+3_fp
    REAL(fp),  PARAMETER  :: DhrKcl = 6900.e+0_fp
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: Hcl_T, Kcl_T
    REAL(fp)              :: Hcleff, xCl, pHCl

    !=================================================================
    ! kHCl begins here!
    !=================================================================

    ! HCl dissolution constants
    HCl_T  = Hcl*exp(Dhcl*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))
    Kcl_T  = Kcl*exp(DhrKcl*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))

    !HCl dissolution
    Hcleff = Hcl_T*(1.0e+0_fp+(Kcl_T/HPLUS))
    xCl    = 1.0e+0_fp / ( 1.0e+0_fp + ( Hcleff * 0.08205e+0_fp * T * LWC ) )
    pHCl   = Cl * P * xCl

    kHClp  = Hcl_T * Kcl_T * pHCl / HPLUS

  END FUNCTION kHCl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkHCl
!
! !DESCRIPTION: Function dkHCl
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkHCl ( P, T, LWC, HPLUS, Cl ) RESULT ( KHClp )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, Cl
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KHClp
!
! !REVISION HISTORY:
!  03 Apr 2019 - X. Wang    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! HCl dissociation constants
      REAL(fp),  PARAMETER  :: Kcl = 1.74e+6_fp
      REAL(fp),  PARAMETER  :: Hcl = 1.5e+3_fp
      REAL(fp),  PARAMETER  :: Dhcl = 2.3e+3_fp
      REAL(fp),  PARAMETER  :: DhrKcl = 6900.e+0_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hcl_T, Kcl_T

      !=================================================================
      ! dkHCl begins here!
      !=================================================================

      ! HCl dissolution constants
      Hcl_T = Hcl*exp(Dhcl*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))
      Kcl_T = Kcl*exp(DhrKcl*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))

      ! HCl dissolution
      ! The fix follows Eq. 7.59 of Seinfeld and Pandis (2006, pp 301).

      kHClp = -1.0e+0_fp * Kcl_T * Hcl_T * Cl * P *          &
           ( 1.0e+0_fp + Hcl_T * 0.08205e+0_fp * T * LWC ) / &
           ( Kcl_T * Hcl_T * 0.08205e+0_fp * T * LWC +       &
           Hcl_T * 0.08205e+0_fp * T * LWC * HPLUS +         &
           HPLUS )**2

      END FUNCTION dkHCl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kNH3
!
! !DESCRIPTION: Function kNH3
!\\
!\\
! !INTERFACE:
!
  FUNCTION kNH3 ( P, T, LWC, HPLUS, NH3, Kw ) RESULT ( KNH3p )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, NH3, Kw
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KNH3p
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! NH3 dissociation contants
    REAL(fp),  PARAMETER  :: Ka1    =  1.7e-5_fp
    REAL(fp),  PARAMETER  :: Dhnh3  =  4200.0_fp
#ifdef LUO_WETDEP
    REAL(fp),  PARAMETER  :: Hnh3   =  59.8_fp
    REAL(fp),  PARAMETER  :: DhrKa1 = -4325.0_fp
#else
    REAL(fp),  PARAMETER  :: Hnh3   =  60.0_fp
    REAL(fp),  PARAMETER  :: DhrKa1 = -450.0_fp
#endif

    ! Variables
    REAL(fp)              :: Hnh3_T,  Ka1_T
    REAL(fp)              :: Hnh3eff, xNH3, pNH3

    !=================================================================
    ! kNH3 begins here!
    !=================================================================

    !NH3 dissolution constants
    Hnh3_T  = Hnh3*exp(Dhnh3*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Ka1_T   = Ka1*exp(DhrKa1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    !NH3 dissolution
    Hnh3eff = Hnh3_T*(1.e+0_fp+((Ka1_T* HPLUS) / Kw))
    xNH3    = 1.e+0_fp / ( 1.e+0_fp + ( Hnh3eff * 0.08205e+0_fp * T * LWC ) )
    pNH3    = NH3 * P * xNH3

    KNH3p   = HPLUS * Hnh3_T * Ka1_T * pNH3 / Kw

  END FUNCTION kNH3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkNH3
!
! !DESCRIPTION: Function dkNH3
!\\
!\\
! !INTERFACE:
!
  FUNCTION dkNH3 ( P, T, LWC, HPLUS, NH3, Kw ) RESULT ( KNH3p )
!

! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, NH3, Kw
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KNH3p
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  22 Mar 2017 - M. Sulprizio- Dhnh3 value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (V. Shah)
!  15 Feb 2019 - J. Moch     - updated function to make output
!  derivative of [NH4+]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! NH3 dissociation contants
    REAL(fp),  PARAMETER  :: Ka1    =  1.7e-5_fp
    REAL(fp),  PARAMETER  :: Dhnh3  =  4200.0_fp
#ifdef LUO_WETDEP
    REAL(fp),  PARAMETER  :: Hnh3   =  59.8_fp
    REAL(fp),  PARAMETER  :: DhrKa1 = -4325.0_fp
#else
    REAL(fp),  PARAMETER  :: Hnh3   =  60.0_fp
    REAL(fp),  PARAMETER  :: DhrKa1 = -450.0_fp
#endif

    ! Variables
    REAL(fp)              :: Hnh3_T, Ka1_T

    !=================================================================
    ! dkNH3 begins here!
    !=================================================================

    !NH3 dissolution constants
    Hnh3_T = Hnh3*exp(Dhnh3*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Ka1_T = Ka1*exp(DhrKa1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    !NH3 dissolutionnyn

    KNH3p = Ka1_T * Hnh3_T * NH3 * Kw * P * ( 1.0e+0_fp +    &
         Hnh3_T * 0.08205e+0_fp * T * LWC ) /                &
         ( Hnh3_T * 0.08205e+0_fp * T * LWC * ( Kw + Ka1_T * &
         HPLUS ) + Kw)**2

  END FUNCTION dkNH3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kFA
!
! !DESCRIPTION: Function kFA
!\\
!\\
! !INTERFACE:
!
  FUNCTION kFA ( P, T, LWC, HPLUS, FA ) RESULT ( kFAp )
!

! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, FA
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KFAp
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  17 Oct 2017 - M. Sulprizio- Dhck value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (Qianjie Chen)
!  03 Dec 2018 - J. Moch     - Modified for formic acid (HCOOH). Values
!                              taken from Sienfeld and Pandis. Made it
!                              to output is [FA]
!  01 May 2020 - V. Shah     - Use correct equilibrium constants
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! HCOOH dissociation constants
    REAL(fp),  PARAMETER  :: Kformate = 1.8e-4_fp ! equib const
    REAL(fp),  PARAMETER  :: Hfa      = 8800.0_fp ! henry const
    REAL(fp),  PARAMETER  :: Dhfa     = 6100.0_fp ! henry temp
    REAL(fp),  PARAMETER  :: DhrKfa   = 151.0_fp  ! equib temp
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: Hfa_T, Kfa_T
    REAL(fp)              :: Hfaeff, xFA, pFA

    !=================================================================
    ! kFA begins here!
    !=================================================================

    ! Formic acid dissolution constants
    HFA_T = Hfa*exp(Dhfa*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))
    Kfa_T = Kformate*exp(DhrKfa*((1.0e+0_fp/T) &
         - (1.0e+0_fp/298.0e+0_fp)))

    !HCOOH  dissolution
    Hfaeff = Hfa_T*(1.0e+0_fp+(Kfa_T/HPLUS))
    xFA = 1.0e+0_fp / ( 1.0e+0_fp &
         + ( Hfaeff * 0.08205e+0_fp * T * LWC ) )
    pFA = FA * P * xFA

    kFAp = Hfa_T * Kfa_T * pFA / HPLUS

  END FUNCTION kFA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkFA
!
! !DESCRIPTION: Function dkFA
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkFA ( P, T, LWC, HPLUS, FA ) RESULT ( kFAp )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, FA
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KFAp
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  17 Oct 2017 - M. Sulprizio- Dhck value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (Qianjie Chen)
!  03 Dec 2018 - J. Moch     - Modified for formic acid (HCOOH). Values
!  taken from
!                              Sienfeld and Pandis. Made it to output is
!                              [FA]
!  01 May 2020 - V. Shah     - Use correct equilibrium constants
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! HCOOH dissociation constants
      REAL(fp),  PARAMETER  :: Kformate = 1.8e-4_fp ! equib const
      REAL(fp),  PARAMETER  :: Hfa = 8800e+0_fp ! henry const
      REAL(fp),  PARAMETER  :: Dhfa = 6100e+0_fp ! henry temp
      REAL(fp),  PARAMETER  :: DhrKfa = 151.e+0_fp ! equib temp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hfa_T, Kfa_T

      !=================================================================
      ! dkFA begins here!
      !=================================================================

      ! Formic acid dissolution constants
      HFA_T = Hfa*exp(Dhfa*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))
      Kfa_T = Kformate*exp(DhrKfa*((1.0e+0_fp/T) &
           - (1.0e+0_fp/298.0e+0_fp)))

      !HCOOH  dissolution

      kFAp = -1.0e+0_fp * Kfa_T * HFA_T * FA * P *           &
           ( 1.0e+0_fp + HFA_T * 0.08205e+0_fp * T * LWC ) / &
           ( Kfa_T * HFA_T * 0.08205e+0_fp * T * LWC +       &
           HFA_T * 0.08205e+0_fp * T * LWC * HPLUS +         &
           HPLUS )**2

      END FUNCTION dkFA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kAA
!
! !DESCRIPTION: Function kAA
!\\
!\\
! !INTERFACE:
!
      FUNCTION kAA ( P, T, LWC, HPLUS, AA ) RESULT ( kAAp )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, AA
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KAAp
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  17 Oct 2017 - M. Sulprizio- Dhck value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (Qianjie Chen)
!  03 Dec 2018 - J. Moch     - Modified for acetic acid (CH3COOH).
!  Values taken from
!                              Sienfeld and Pandis, value of [HCOOH]
!  01 May 2020 - V. Shah     - Use correct equilibrium constants
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! CH3HCOOH dissociation constants
      REAL(fp),  PARAMETER  :: Kacetate = 1.75e-5_fp
      REAL(fp),  PARAMETER  :: Haa      = 4100.0_fp
      REAL(fp),  PARAMETER  :: Dhaa     = 6200.0_fp
      REAL(fp),  PARAMETER  :: DhrKaa   = 50.0_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Haa_T, Kaa_T
      REAL(fp)              :: Haaeff, xAA, pAA
      !=================================================================
      ! kAA begins here!
      !=================================================================

      ! Formic acid dissolution constants
      HAA_T = Haa*exp(Dhaa*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))
      Kaa_T = Kacetate*exp(DhrKaa*((1.0e+0_fp/T) &
          - (1.0e+0_fp/298.0e+0_fp)))

      !HCOOH  dissolution
      Haaeff = Haa_T*(1.0e+0_fp+(Kaa_T/HPLUS))
      xAA = 1.0e+0_fp / ( 1.0e+0_fp &
         + ( Haaeff * 0.08205e+0_fp * T * LWC ) )
      pAA = AA * P * xAA

      kAAp = Haa_T * Kaa_T * pAA / HPLUS

      END FUNCTION kAA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkAA
!
! !DESCRIPTION: Function kdAA
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkAA ( P, T, LWC, HPLUS, AA ) RESULT ( kAAp )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, AA
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KAAp
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  17 Oct 2017 - M. Sulprizio- Dhck value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (Qianjie Chen)
!  03 Dec 2018 - J. Moch     - Modified for acetic acid (CH3COOH).
!  Values taken from
!                              Sienfeld and Pandis. Output is
!                              derivative.
!  01 May 2020 - V. Shah     - Use correct equilibrium constants
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! HCOOH dissociation constants
      REAL(fp),  PARAMETER  :: Kacetate = 1.75e-5_fp
      REAL(fp),  PARAMETER  :: Haa      = 4100.0_fp
      REAL(fp),  PARAMETER  :: Dhaa     = 6200.0_fp
      REAL(fp),  PARAMETER  :: DhrKaa   = 50.0_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Haa_T, Kaa_T
      !=================================================================
      ! kAA begins here!
      !=================================================================

      ! Formic acid dissolution constants
      HAA_T = Haa*exp(Dhaa*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))
      Kaa_T = Kacetate*exp(DhrKaa*((1.0e+0_fp/T) &
          - (1.0e+0_fp/298.0e+0_fp)))

      !HCOOH  dissolution
      kAAp =  -1.0e+0_fp * Kaa_T * HAA_T * AA * P * &
          ( 1.0e+0_fp + HAA_T * 0.08205e+0_fp * T * LWC ) / &
          ( Kaa_T * HAA_T * 0.08205e+0_fp * T * LWC + &
          HAA_T * 0.08205e+0_fp * T * LWC * HPLUS + &
          HPLUS )**2


      END FUNCTION dkAA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CaCO3_PRECIP
!
! !DESCRIPTION: Subroutine CaCO3 to calculate [Ca++] in equilibrium with
! CaCO3(s) (dust particles) depending on [H+]
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CaCO3_PRECIP ( P,  T, HPLUS, fCa, dCa )
!
! !INPUT PARAMETERS:
!
      REAL(fp),        INTENT(IN) :: T, P, HPLUS
!
! !OUTPUT PARAMETERS:
!
      REAL(fp),  INTENT(OUT):: fCa, dCa ! [Ca2+] and d([Ca2+])/d[H+]
!
! !REVISION HISTORY:
!  25 Dec 2019 - V. Shah - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! !DEFINED PARAMETERS:
!
      REAL(fp),  PARAMETER  :: Kc1 = 4.3e-7_fp
      REAL(fp),  PARAMETER  :: Kc2 = 4.68e-11_fp
      REAL(fp),  PARAMETER  :: DhrKc1 = -1000.
      REAL(fp),  PARAMETER  :: DhrKc2 = -1760.
      REAL(fp),  PARAMETER  :: Hco2 = 3.4e-2_fp
      REAL(fp),  PARAMETER  :: Dhco2 = 2.44e+3_fp
      ! CO2 concentration [v/v]
      REAL(fp),  PARAMETER  :: CO2 = 390.0e-6_fp
      REAL(fp),  PARAMETER  :: Ksp = 3.3e-9_fp
      REAL(fp),  PARAMETER  :: DHrKsp = -1200e+0_fp

! !LOCAL VARIABLES:
      REAL(fp)              :: HCO2_T, Kc1_T, Kc2_T, Ksp_T

! !REMARKS:

      !=================================================================
      ! CaCO3_PRECIP begins here!
      !=================================================================
      !Temperature adjusted eq. constants
      !CO2 dissolution constants
      Hco2_T = Hco2*exp(Dhco2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc1_T = Kc1*exp(DhrKc1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc2_T = Kc2*exp(DhrKc2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      ! CaCO3 eq constants
      Ksp_T = Ksp*exp(DhrKsp*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      !Ca concentrations [M]
      fCa = Ksp_T * HPLUS * HPLUS / (Kc1_T * Kc2_T * Hco2_T * CO2 * P)
      !derivative d[Ca2+]/dH+
      dCa  = 2e+0_fp * Ksp_T * HPLUS / (Kc1_T * Kc2_T * Hco2_T * CO2 * P)

      END SUBROUTINE CaCO3_PRECIP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aqchem_so2
!
! !DESCRIPTION: Subroutine AQCHEM\_SO2 computes the reaction rates for aqueous
! SO2 chemistry. (rjp, bmy, 10/31/02, 12/12/02)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AQCHEM_SO2( LWC,     T,     P,       SO2,   H2O2, &
                         O3,      Hplus, MnII,    FeIII, IONIC, &
                         KaqH2O2, KaqO3, KaqO3_1, KaqO2, &
                         HSO3aq,  SO3aq, HCHO, KaqHCHO,  &
                         KaqHMS, KaqHMS2 )
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  :: LWC     ! Liq water content [m3/m3]=1.E-6*L [g/m3]
    REAL(fp), INTENT(IN)  :: T       ! Temperature [K]
    REAL(fp), INTENT(IN)  :: P       ! Dry air partial pressure [atm]
    REAL(fp), INTENT(IN)  :: SO2     ! SO2  mixing ratio [v/v]
    REAL(fp), INTENT(IN)  :: H2O2    ! H2O2 mixing ratio [v/v]
    REAL(fp), INTENT(IN)  :: O3      ! O3   mixing ratio [v/v]
    REAL(fp), INTENT(IN)  :: HPLUS   ! Concentration of H+ ion (i.e. pH) [v/v]
    REAL(fp), INTENT(IN)  :: MnII    ! Concentration of MnII [mole/l]
    REAL(fp), INTENT(IN)  :: FeIII   ! Concentration of FeIII [mole/l]
    REAL(fp), INTENT(IN)  :: IONIC   ! Ionic strength [mole/l]?
    REAL(fp), INTENT(IN)  :: HCHO    ! HCHO   mixing ratio [v/v] (jmm, 06/13/18)

!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: KaqH2O2 ! Reaction rate for H2O2
    REAL(fp), INTENT(OUT) :: KaqO3   ! Reaction rate for O3
    REAL(fp), INTENT(OUT) :: KaqO3_1 ! only the SO3-- oxidation, (qjc, 04/10/16)
    REAL(fp), INTENT(OUT) :: KaqO2   ! Reaction rate for O2 (metal cat)
    REAL(fp), INTENT(OUT) :: KaqHCHO ! Reaction rate for SO2 and HCHO (jmm, 06/13/18)
    REAL(fp), INTENT(OUT) :: KaqHMS  ! Reaction rate for HMS and OH- (jmm, 06/13/18)
    REAL(fp), INTENT(OUT) :: KaqHMS2  ! Reaction rate for HMS and OH(aq) (jmm, 06/28/18)
    REAL(fp), INTENT(OUT) :: HSO3aq  ! Cloud bisulfite [mol/l] (qjc, 06/10/16)
    REAL(fp), INTENT(OUT) :: SO3aq   ! Cloud sulfite   [mol/l] (qjc, 06/10/16)
!
! !REMARKS:
!  Chemical Reactions:
!  ============================================================================
!  (R1) HSO3- + H2O2(aq) + H+ => SO4-- + 2H+ + H2O [Jacob, 1986]
!                                                                             .
!      d[S(VI)]/dt = k[H+][H2O2(aq)][HSO3-]/(1 + K[H+])
!      [Seinfeld and Pandis, 1998, page 366]
!                                                                             .
!  (R2) SO2(aq) + O3(aq) =>
!       HSO3-   + O3(aq) =>
!       SO3--   + O3(aq) =>
!       [Jacob, 1986; Jacobson, 1999]
!                                                                             .
!       d[S(VI)]/dt = (k0[SO2(aq)] + k1[HSO3-] + K2[SO3--])[O3(aq)]
!       [Seinfeld and Pandis, 1998, page 363]
!                                                                             .
!  (R3) HSO3-   + HCHO(aq) => HMS
!       SO3--   + HCHO(aq) => HMS + OH-
!       [Moch et al., 2018; Olson and Hoffman, 1986]
!                                                                             .
!       d[S(HMS)]/dt = (k1[HSO3-] + k2[SO3--])[HCHO(aq)]
!       [Seinfeld and Pandis, 2016, 309]
!
!  (R4) HMS + OH- => HCHO(aq) + SO3--
!       [Moch et al., 2018; Deister et al., 1986]
!        (note treated as 1st order in contrast to other reactions here)
!
!  (R5) HMS + OH(aq) =(SO2,HO2,O2)=> HCHO + 2SO4-- + O2 + 3H+ + 2H2O
!       [Jacob et al, 1986, Olson and Fessenden, 1992;
!        Seinfeld and Pandis, 2016, Table 7A.7]
!          Net reaction (R5):
!           HMS + OH(aq) =(O2)=> SO5- + HCHO + H2O
!           HO2 <=> H+ + O2-
!           SO5- + O2- =(H2O)=> HSO5- + OH- + O2
!           SO2(aq) <=> HSO3- + H+
!           H+ + OH- <=> H2O
!           HSO5- + HSO3- => 2SO4-- + 2H+
!
!  Reaction rates can be given as
!       Ra     = k [H2O2(ag)] [S(IV)]  [mole/liter*s]  OR
!       Krate  = Ra LWC R T / P        [1/s]
!                                                                             .
!  Where:
!       LWC = Liquid water content(g/m3)*10-6 [m3(water)/m3(gas)]
!       R   = 0.08205  (atm L / mol-K), Universal gas const.
!       T   = Temperature (K)
!       P   = Pressure (atm)
!                                                                             .
!  Procedure:
!  ============================================================================
!  (a ) Given [SO2] which is assumed to be total SO2 (gas+liquid) in
!        equilibrium between gas and liquid phase.
!                                                                             .
!  (b ) We can compute SO2(g) using Henry's law
!          P(so2(g)) = Xg * [SO2]
!          Xg = 1/(1 + Faq), Fraction of SO2 in gas
!       where:
!          Faq   = Kheff * R * T * LWC,
!          KHeff = Effective Henry's constant
!                                                                             .
!  (c ) Then Calculate Aquous phase, S[IV] concentrations
!        S[IV] = Kheff * P(so2(g) in atm) [M]
!                                                                             .
!  (d ) The exact same procedure is applied to calculate H2O2(aq) and HCHO(aq)
!
! !REVISION HISTORY:
!  (1 ) Updated by Rokjin Park (rjp, bmy, 12/12/02)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER   :: R = 0.08205e+0_fp
    REAL(fp), PARAMETER   :: dOH = 1.0e-19_fp ! [M cm^3 molec^-1] (jmm, 06/28/18)
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: KH2O2,  RA,     KS1, KS2,    HCSO2
    REAL(fp)              :: FHCSO2, XSO2G,  SIV, HSO3,   XSO2AQ
    REAL(fp)              :: XHSO3,  XSO3,   KH1, HCH2O2, FHCH2O2
    REAL(fp)              :: XH2O2G, H2O2aq, KO0, KO1,    KO2
    REAL(fp)              :: HCO3,   XO3g,   O3aq, XHCHOg, HCHCHO ! (jmm, 06/13/18)
    REAL(fp)              :: FHCHCHO,KHCHO1, KHCHO2,  KHMS  ! (jmm, 06/13/18)
    REAL(fp)              :: KW1,    KHC1,   KHMS2          ! (jmm, 06/15/18)
    REAL(fp)              :: Eff_mn, Eff_fe !jys



    !=================================================================
    ! AQCHEM_SO2 begins here!
    !
    ! Aqueous reaction rate
    ! HSO3- + H2O2 + H+ => SO4-- + 2H+ + H2O [Jacob, 1986]
    !=================================================================

    ! [Jacob, 1986]
    KH2O2 = 6.31e+14_fp * EXP( -4.76e+3_fp / T )

    !! [Jacobson, 1999]
    !KH2O2 = 7.45e+0_fp7 * EXP( -15.96e+0_fp * ( (298.15/T) - 1.) ) / &
    !        ( 1.e+0_fp + 13.e+0_fp * Hplus)

    !=================================================================
    ! Equilibrium reaction of SO2-H2O
    !    SO2 + H2O = SO2(aq)        (s0)
    !    SO2(ag)   = HSO3- + H+     (s1)
    !    HSO3-     = SO3-- + H+     (s2)
    !
    ! Reaction constant for Aqueous chemistry -- No big difference
    ! between Jacob and Jacobson, choose one of them.
    !
    ! Reaction rate dependent on Temperature is given
    !   H = A exp ( B (T./T - 1) )
    !
    ! For equilibrium reactions of SO2:
    !            As1      Bs1   As2      Bs2
    !  Seinfeld  1.30d-2  7.02  6.60d-8  3.76   [1998]
    !  Jacob     1.30d-2  6.75  6.31d-8  5.05   [1986]
    !  Jacobson  1.71d-2  7.04  5.99d-8  3.74   [1996]
    !=================================================================
    Ks1    = 1.30e-2_fp * EXP( 6.75e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp ) )
    Ks2    = 6.31e-8_fp * EXP( 5.05e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp ) )

    ! SIV Fraction
    XSO2aq = 1.e+0_fp/(1.e+0_fp + Ks1/Hplus + Ks1*Ks2/(Hplus*Hplus))
    XHSO3  = 1.e+0_fp/(1.e+0_fp + Hplus/Ks1 + Ks2/Hplus)
    XSO3   = 1.e+0_fp/(1.e+0_fp + Hplus/Ks2 + Hplus*Hplus/(Ks1*Ks2))

    ! Henry's constant [mol/l-atm] and Effective Henry's constant for SO2
    HCSO2  = 1.22e+0_fp * EXP( 10.55e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp) )
    FHCSO2 = HCSO2 * (1.e+0_fp + (Ks1/Hplus) + (Ks1*Ks2 / (Hplus*Hplus)))

    XSO2g  = 1.e+0_fp / ( 1.e+0_fp + ( FHCSO2 * R * T * LWC ) )
    SIV    = FHCSO2 * XSO2g * SO2 * P
    !HSO3   = Ks1 * HCSO2 * XSO2g * SO2 * P

    ! Effective HSO3aq for HOBr+HSO3
    HSO3aq = SIV * XHSO3           ! unit: M (qjc, 06/10/16)

    ! Effective SO3aq for HOBr+SO3
    SO3aq  = SIV * XSO3            ! unit: M (qjc, 06/10/16)

    !=================================================================
    ! H2O2 equilibrium reaction
    ! H2O2 + H2O = H2O2.H2O
    ! H2O2.H2O   = HO2- + H+   1)
    !
    ! Reaction rate dependent on Temperature is given
    !   H = A exp ( B (T./T - 1) )
    !
    ! For equilibrium reactions of SO2
    !            Ah1       Bh1
    !  Jacob     1.58E-12  -12.49  [1986]
    !  Jacobson  2.20E-12  -12.52  [1996]
    !=================================================================
    Kh1 = 2.20e-12_fp * EXP( -12.52e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp ) )

    ! Henry's constant [mol/l-atm] and Effective Henry's constant for H2O2
    ! [Seinfeld and Pandis, 1998]
    ! HCH2O2  = 7.45D4 * EXP( 24.48e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp) )

    ! [Jacobson,1999]
    HCH2O2  = 7.45e+4_fp * EXP( 22.21e+0_fp * (298.15e+0_fp / T - 1.e+0_fp) )
    FHCH2O2 = HCH2O2 * (1.e+0_fp + (Kh1 / Hplus))

    XH2O2g  = 1.e+0_fp / ( 1.e+0_fp + ( FHCH2O2 * R * T * LWC ) )
    !H2O2aq  = FHCH2O2 * XH2O2g * H2O2 * P

    ! Conversion rate from SO2 to SO4 via reaction with H2O2
    KaqH2O2  = kh2o2 * Ks1 * FHCH2O2 * HCSO2 * XH2O2g * XSO2g &
               * P * LWC * R * T            ! [v/v/s]

    !=================================================================
    !  Aqueous reactions of SO2 with O3
    !  SO2(aq) + O3 =>                       (0)
    !  HSO3-   + O3 => SO4-- + H+ + O2       (1)
    !  SO3--   + O3 => SO4-- + O2            (2)
    !
    ! NOTE
    ! [Jacob, 1986]
    !    KO1  = 3.49E12 * EXP( -4.83E3 / T )
    !    KO2  = 7.32E14 * EXP( -4.03E3 / T )
    !
    ! [Jacobson, 1999]
    !    KO0  = 2.40E+4
    !    KO1  = 3.70E+5 * EXP( -18.56 * ((298.15/T) - 1.))
    !    KO2  = 1.50E+9 * EXP( -17.72 * ((298.15/T) - 1.))
    !
    ! Rate constants from Jacobson is larger than those of Jacob
    ! and results in faster conversion from S(IV) to S(VI)
    ! We choose Jacob 1) 2) and Jacobson 0) here
    !=================================================================
    KO0 = 2.40e+4_fp
    KO1 = 3.49e+12_fp * EXP( -4.83e+3_fp / T )
    KO2 = 7.32e+14_fp * EXP( -4.03e+3_fp / T )

    !=================================================================
    ! H2O2 equilibrium reaction
    ! O3 + H2O = O3.H2O
    !  HCO3  = 1.13E-2 * EXP( 8.51 * (298.15/T -1.) ), S & P
    !  HCO3  = 1.13E-2 * EXP( 7.72 * (298.15/T -1.) ), Jacobson
    !=================================================================

    ! Calculate Henry's Law constant for atmospheric temperature
    HCO3  = 1.13e-2_fp * EXP( 8.51e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp ) )

    XO3g  = 1.e+0_fp / ( 1.e+0_fp + ( HCO3 * R * T * LWC ) )
    !O3aq  = HCO3 * XO3g * O3 * P

    ! Conversion rate from SO2 to SO4 via reaction with O3
    KaqO3 = (KO0*XSO2AQ + KO1*XHSO3 + KO2*XSO3) * FHCSO2 * XSO2g &
            * P * HCO3 * XO3g * LWC * R * T   ! [v/v/s]

    !(qjc, 04/10/16)
    KaqO3_1 = KO2*XSO3 * FHCSO2 * XSO2g &
              * P * HCO3 * XO3g * LWC * R * T   ! [v/v/s]

    ! ===================================================================
    ! Metal (Fe, Mn) catalyzed O2 oxidation (bec, 7/12/04)
    ! R = d[S(VI)]/dt = 750*[Mn(II)]*[S(IV)] + 2600*[Fe(III)]*[S(IV)] +
    !               1.d10*[Mn(II)]*[Fe(III)]*[S(IV)]
    ! from Seinfeld and Pandis, 1998 pg. 371
    ! S(IV) = HFCSO2 * XSO2*P*[SO2]
    ! R = KaqO2*[SO2] (v/v/s)
    ! KaqO2 = FHCSO2 * XSO2g * P *
    !        ((750*[Mn(II)])+(2600[Fe(III)])+(1.d10*[Mn(II)]*[Fe(III)]))
    ! in units of [M/s]
    ! KaqO2 = FHCSO2 * XSO2g * P *
    !        ((750*[Mn(II)])+(2600[Fe(III)])+(1.d10*[Mn(II)]*[Fe(III)])) *
    !        LWC * R * T/P
    ! in units of [v/v/s]
    ! ===================================================================

    ! Conversion rate from SO2 to SO4 via reaction with O2 (met cat)
    ! Commented out becasue using ionic strength pH modifiers version
    !KaqO2 = FHCSO2 * XSO2g * ( (750e+0_fp * MnII ) + &
    !        ( 2600e+0_fp * FeIII ) + (1e+10_fp * MnII * FeIII ) ) * &
    !        LWC * R * T   ! [s-1]

    ! Conversion rate from SO2 to SO4 via reaction with O2 (met cat)
    ! added by shaojy16  10/13/2017
    ! takes into account pH and ionic strength
    Eff_mn = 10.0**(-4.0*(SQRT(IONIC)/(1.0+SQRT(IONIC))))
    Eff_fe = 10.0**(-2.0*(SQRT(IONIC)/(1.0+SQRT(IONIC))))

    IF ( Hplus > 10.0**(-4.2) ) THEN
       KaqO2 = FHCSO2 * XSO2g * &
            (3.72e+7_fp*Hplus**(-0.74)* &
            (MnII*FeIII*Eff_fe*Eff_mn)) * &
            LWC * R * T   ! [s-1]
    ELSE
       KaqO2 = FHCSO2 * XSO2g * &
            (2.51e+13_fp*Hplus**(0.67) * &
            (MnII*FeIII*Eff_fe*Eff_mn)) * &
            LWC * R * T   ! [s-1]
    ENDIF

    !=================================================================
    !  Aqueous reactions of SO2 with HCHO
    !     HSO3-   + HCHO(aq) => HMS + OH-           (1)
    !     SO3--   + HCHO(aq) => HMS                 (2)
    !
    !     NOTE:
    !     [Boyce and Hoffman, 1984]
    !        KHCHO1  = 7.9E2 * EXP( -16.435 * ((298.15/T) - 1.))
    !        KHCHO2  = 2.5E7 * EXP( -6.037 * ((298.15/T) - 1.))
    !
    !
    !  Aqueous reaction of HMS with OH-
    !    HMS + OH- => HCHO(aq) + SO3--             (3)
    !
    !     NOTE: unclear where B (E/R) value in Seinfeld and Pandis from,
    !     but close to Deister. Using Seinfeld and Pandis value for now
    !     [Deister et al., 1986]
    !        KHMS    = 3.6E3 * EXP( -22.027 * ((298.15/T) - 1.))
    !     [Seinfeld and Pandis, 2016; Munger et al., 1986]
    !        KHMS    = 3.6E3 * EXP( -15.09 * ((298.15/T) - 1.))
    !
    !
    !  Aqueous reaction of HMS with OH(aq)
    !    HMS + OH(aq) =(SO2,O2,HO2)=> 2SO4-- + HCHO + O2 + 3H+ + 2H2O  (4)
    !
    !    NOTE: O2, SO2, and HO2 particpate in the stoichiometry but not kinetics.
    !          Assume steady state for sulfur radicals and the following reaction chain:
    !            HMS + OH(aq) =(O2)=> SO5- + HCHO + H2O [Olsen and Fessenden, 1992]
    !            HO2 <=> H+ + O2-                       [Jacob, 1986]
    !            SO5- + O2- =(H2O)=> HSO5- + OH- + O2   [Jacob, 1986]
    !            SO2(aq) <=> HSO3- + H+
    !            H+ + OH- <=> H2O
    !            HSO5- + HSO3- => 2SO4-- + 2H+          [Jacob, 1986]
    !       Instead of assuming Henry's law for OH, use the parameter from
    !       Jacob et al, 2005 that relates gas phase OH to aqueous phase OH
    !       accounting for the HO2(aq)/O2- cylcing in cloud droplets:
    !        dOH = 1E-19 [M cm^3 molec^-1]
    !     [Olson and Fessenden, 1992]
    !        KHMS2    = 6.2E8 * EXP( -5.03 * ((298.15/T) -1.))
    !
    !
    ! (jmm, 06/28/18)
    !=================================================================
    KHCHO1 = 7.9e+2_fp * EXP( -16.44e+0_fp &
         * ( 298.15e+0_fp / T - 1.e+0_fp ) )
    KHCHO2 = 2.5e+7_fp * EXP( -6.04e+0_fp &
         * ( 298.15e+0_fp / T - 1.e+0_fp ) )
    KHMS = 3.6e+3_fp * EXP( -15.09e+0_fp &
         * ( 298.15e+0_fp / T - 1.e+0_fp ) )
    KHMS2 = 6.2e+8_fp * EXP( -5.03e+0_fp &
         * ( 298.15e+0_fp / T - 1.e+0_fp ) )

    !=================================================================
    ! HCHO equilibrium reaction
    ! HCHO(aq) + H2O   = HCH(OH)2
    !
    ! Reaction rate dependent on Temperature is given
    !   H = A exp ( B (T./T - 1) )
    !
    ! For equilibrium reactions of HCHO
    !                             Ah1       Bh1
    !  Sienfeld and Pandis      2.53E3    13.48  [2016]
    !
    ! (jmm, 06/15/18)
    !=================================================================
    Khc1 = 2.53e+3_fp * EXP( 13.48e+0_fp &
         * ( 298.15e+0_fp / T - 1.e+0_fp ) )

    !=================================================================
    ! H2O equilibrium reaction
    !  H2O   = H+ + OH-
    !
    ! Reaction rate dependent on Temperature is given
    !   H = A exp ( B (T./T - 1) )
    !
    ! For equilibrium reactions of HCHO
    !                             Ah1       Bh1
    !  Sienfeld and Pandis       1E-14     -22.51  [2016]
    !
    ! (jmm, 06/15/18)
    !=================================================================
    Kw1 = 1e-14_fp * EXP( -22.51e+0_fp &
         *  ( 298.15e+0_fp / T - 1.e+0_fp ) )

    ! Henry's constant [mol/l-atm] and Effective Henry's constant for HCHO
    ! [Seinfeld and Pandis, 2016]
    ! HCHCHO  = 2.5 * EXP( 21.6e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp) )
    ! (jmm, -6/15/18)
    HCHCHO  = 2.5e+0_fp * EXP( 21.6e+0_fp &
         *  (298.15e+0_fp / T - 1.e+0_fp) )
    FHCHCHO = HCHCHO * (1.e+0_fp + Khc1 )

    XHCHOg  = 1.e+0_fp / ( 1.e+0_fp + ( FHCHCHO * R * T * LWC ) )


    ! Conversion rate from SO2 to HMS via reaction with HCHO
    ! (jmm, 06/15/18)
    KaqHCHO = (KHCHO1*XHSO3 + KHCHO2*XSO3) * FHCSO2 * XSO2G &
         * P * HCHCHO * XHCHOg * LWC * R * T    ! [v/v/s]

    ! Conversion rate from HMS to SO2 via reaction with OH-
    ! (jmm, 06/15/18)
    KaqHMS = KHMS * ( Kw1 / Hplus )            ! [mol/L/s]

    ! Conversion rate from HMS to SO42- & HCHO via reaction with OH(aq)
    ! (jmm, 06/28/18)
    KaqHMS2 = KHMS2 * dOH / AIRMW / AIRMW * 7.e-4_fp * AVO &
         * T * R / P                           ! [m^6 kg^-2 s^-1]

  END SUBROUTINE AQCHEM_SO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: het_drop_chem
!
! !DESCRIPTION: Subroutine HET\_DROP\_CHEM estimates the in-cloud sulfate
!  production rate in heterogeneous cloud droplets based on the Yuen et al.,
!  1996 parameterization. (bec, 6/16/11)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE HET_DROP_CHEM( I,    J,    L,      LSTOT, SSCvv, &
                            aSO4, GNH3, SO2_sr, H2O20, GNO3,  &
                            SR,   Input_Opt,    State_Met, State_Chm )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)  :: I, J, L
    REAL(fp),       INTENT(IN)  :: LSTOT
    REAL(fp),       INTENT(IN)  :: SSCvv
    REAL(fp),       INTENT(IN)  :: aSO4
    REAL(fp),       INTENT(IN)  :: GNH3
    REAL(fp),       INTENT(IN)  :: SO2_sr
    REAL(fp),       INTENT(IN)  :: H2O20
    REAL(fp),       INTENT(IN)  :: GNO3
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT) :: SR          ! Sulfate production rate
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER   :: SS_DEN = 2200.e+0_fp !dry sea-salt density [kg/m3]

    ! sigma of the size distribution for sea-salt (Jaegle et al., 2011)
    REAL(fp), PARAMETER   :: SIG_S = 1.8e+0_fp

    ! geometric dry mean diameters [m] for computing lognormal size distribution
    REAL(fp), PARAMETER   :: RG_S = 0.4e-6_fp !(Jaegle et a., 2011)
    REAL(fp), PARAMETER   :: RG_D2 = 1.5e-6_fp!(Ginoux et al., 2001)
    REAL(fp), PARAMETER   :: RG_D3 = 2.5e-6_fp
    REAL(fp), PARAMETER   :: RG_D4 = 4.e-6_fp
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: alpha_NH3, alpha_SO2, alpha_H2O2
    REAL(fp)              :: alpha_HNO3, alpha_B, alpha_CN
    REAL(fp)              :: alpha_W, alpha_SO4, sum_gas, H
    REAL(fp)              :: NDss, NDd2, NDd3, NDd4, CN, DEN, REFF, W
    REAL(fp)              :: DTCHEM, APV, DSVI
    REAL(fp)              :: B, NH3, SO2, H2O2, HNO3, SO4
    REAL(fp)              :: CNss, CNd2, CNd3, CNd4
    REAL(fp)              :: MW_SO4, MW_SALC

    ! Pointers
    REAL(fp), POINTER     :: AD(:,:,:)
    REAL(fp), POINTER     :: AIRDEN(:,:,:)
    REAL(fp), POINTER     :: AIRVOL(:,:,:)
    REAL(fp), POINTER     :: OMEGA(:,:,:)
    REAL(fp), POINTER     :: U(:,:,:)
    REAL(fp), POINTER     :: V(:,:,:)

    !=================================================================
    ! HET_DROP_CHEM begins here!
    !=================================================================

    ! Initialize pointers
    AD     => State_Met%AD
    AIRDEN => State_Met%AIRDEN
    AIRVOL => State_Met%AIRVOL
    OMEGA  => State_Met%OMEGA
    U      => State_Met%U
    V      => State_Met%V

    ! Convert gas phase concentrations from [v/v] to [pptv]
    NH3  = GNH3  * 1.0e+12_fp
    SO2  = SO2_sr  * 1.0e+12_fp
    H2O2 = H2O20 * 1.0e12_fp
    HNO3 = GNO3 * 1.0e12_fp

    ! Set molecular weight local variables
    MW_SO4 = State_Chm%SpcData(id_SO4)%Info%MW_g
    MW_SALC = State_Chm%SpcData(id_SALC)%Info%MW_g

    ! Convert sulfate aerosol concentrations from [v/v] to [ug/m3]
    SO4 = ( aSO4 * AD(I,J,L) * 1.0e+9_fp ) / &
          ( ( AIRMW / MW_SO4 ) * AIRVOL(I,J,L) )

    ! Convert in cloud sulfate production rate from [v/v/timestep] to
    ! [ug/m3/timestep]
    B  = ( LSTOT * AD(I,J,L) * 1.0e+9_fp ) / &
         ( ( AIRMW / MW_SO4 ) * AIRVOL(I,J,L) )

    ! Convert coarse-mode aerosol concentrations from [v/v] to [#/cm3]
    ! based on equation in Hofmann, Science, 1990.
    ! First convert from [v/v] to [kg/m3 air]
    CNss = SSCvv * AD(I,J,L) / ( ( AIRMW / MW_SALC ) * AIRVOL(I,J,L) )

    ! Now convert from [kg/m3 air] to [#/cm3 air]
    ! Sea-salt
    NDss = ( (3.e+0_fp/4.e+0_fp) * CNss ) / (PI * SS_DEN * RG_S**3.e+0_fp * &
           exp( (9.e+0_fp/2.e+0_fp) * (LOG(SIG_S)) ** 2.e+0_fp ) ) * 1.e-6_fp

    ! Total coarse mode number concentration [#/cm3]
    CN = NDss ! sea-salt

    ! Determine regression coefficients based on the local SO2 concentration
    IF ( SO2 <= 200.0e+0_fp ) THEN
       alpha_B    = 0.5318e+0_fp
       alpha_NH3  = -1.67e-7_fp
       alpha_SO2  = 2.59e-6_fp
       alpha_H2O2 = -1.77e-7_fp
       alpha_HNO3 = -1.72e-7_fp
       alpha_W    = 1.22e-6_fp
       alpha_CN   = 4.58e-6_fp
       alpha_SO4  = -1.00e-5_fp
    ELSE IF ( SO2 > 200.0e+0_fp .AND. SO2 <= 500.0e+0_fp ) THEN
       alpha_B    = 0.5591e+0_fp
       alpha_NH3  = 3.62e-6_fp
       alpha_SO2  = 1.66e-6_fp
       alpha_H2O2 = 1.06e-7_fp
       alpha_HNO3 = -5.45e-7_fp
       alpha_W    = -5.79e-7_fp
       alpha_CN   = 1.63e-5_fp
       alpha_SO4  = -7.40e-6_fp
    ELSE IF ( SO2 > 500.0e+0_fp .AND. SO2 < 1000.0e+0_fp ) THEN
       alpha_B    = 1.1547e+0_fp
       alpha_NH3  = -4.28e-8_fp
       alpha_SO2  = -1.23e-7_fp
       alpha_H2O2 = -9.05e-7_fp
       alpha_HNO3 = 1.73e-7_fp
       alpha_W    = 7.22e-6_fp
       alpha_CN   = 2.44e-5_fp
       alpha_SO4  = 3.25e-5_fp
    ELSE IF ( SO2 >= 1000.0e+0_fp ) THEN
       alpha_B    = 1.1795e+0_fp
       alpha_NH3  = 2.57e-7_fp
       alpha_SO2  = -5.54e-7_fp
       alpha_H2O2 = -1.08e-6_fp
       alpha_HNO3 = 1.95e-6_fp
       alpha_W    = 6.14e-6_fp
       alpha_CN   = 1.64e-5_fp
       alpha_SO4  = 2.48e-6_fp
    ENDIF

    ! Updraft velocity over the oceans [cm/s]
    ! 500 cm/s is too high. Get W from the met field. (qjc, 04/10/16)
    !W = 500e+0_fp
    W = -OMEGA(I,J,L) / ( AIRDEN(I,J,L) * g0 ) * 100e+0_fp

    ! Compute H (integration time interval * air parcel velocity) [m]
    ! DTCHEM is the chemistry timestep in seconds
    DTCHEM = GET_TS_CHEM()

    ! Compute air parcel velocity [m/s]
    !APV = SQRT( (U(I,J,L) * U(I,J,L)) + (V(I,J,L) * V(I,J,L)) )
    !(qjc, 04/10/16)
    APV = SQRT( (U(I,J,L) * U(I,J,L)) + (V(I,J,L) * V(I,J,L)) + &
                 W * W *1.e-4_fp )

    H   = DTCHEM * APV          ![m]

    sum_gas = (alpha_NH3 * NH3) + (alpha_SO2 * SO2) + &
              (alpha_H2O2 * H2O2) + (alpha_HNO3 * HNO3)

    DSVI = ( alpha_B * B ) + &
           ( ( (alpha_CN * CN) + (alpha_W * W) + (alpha_SO4 * SO4) + &
                sum_gas ) * H )

    ! Only calculate SR when air parcel rises, in consistence with
    ! Yuen et al. (1996) (qjc, 04/10/16)
    IF ( W > 0e+0_fp ) THEN

       ! additional sulfate production that can be attributed to
       ! ozone [ug/m3/timestep]
       SR = DSVI - B

       ! Convert SR from [ug/m3/timestep] to [v/v/timestep]
       SR = SR * ( AIRMW / MW_SO4 ) * 1.e-9_fp / AIRDEN(I,J,L)

       ! Don't allow SR to be negative
       SR = MAX( SR, 0.e+0_fp )

       ! Don't produce more SO4 than SO2 available after AQCHEM_SO2
       SR = MIN( SR, SO2_sr )

    ELSE
       SR = 0.e+0_fp
    ENDIF

    ! Free pointers
    AD     => NULL()
    AIRDEN => NULL()
    AIRVOL => NULL()
    OMEGA  => NULL()
    U      => NULL()
    V      => NULL()

  END SUBROUTINE HET_DROP_CHEM
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_so4
!
! !DESCRIPTION: Subroutine CHEM\_SO4 is the SO4 chemistry subroutine from Mian
!  Chin's GOCART model, modified for the GEOS-CHEM model.  Now also modified to
!  account for production of crystalline and aqueous sulfur tracers.
!  (rjp, bdf, cas, bmy, 5/31/00, 5/23/06)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_SO4( Input_Opt, State_Chm, State_Diag, State_Grid, &
                       State_Met, RC )
!
! !USES:
!
    USE CMN_SIZE_Mod,   ONLY : NDSTBIN
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE Species_Mod,    ONLY : SpcConc
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
#ifdef APM
    USE APM_INIT_MOD,   ONLY : APMIDS
    USE APM_INIT_MOD,   ONLY : NGCOND,   NSO4,      NSEA
    USE APM_INIT_MOD,   ONLY : NCTSO4,NCTBC,NCTOC,NCTDST,NCTSEA
    USE APM_INIT_MOD,   ONLY : IFEMITBCOCS
    USE APM_DRIV_MOD,   ONLY : PSO4GAS
    USE APM_DRIV_MOD,   ONLY : FCLOUD
#endif
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  The only production is from SO2 oxidation (save in CHEM_SO2).  Dry
!  deposition is now handled in mixing_mod.F90, so we can must add the
!  production from SO2 into the SO4 tracers.
!                                                                             .
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL           :: LDSTUP
    INTEGER           :: I,   J,    L,    N
    INTEGER           :: IBIN  ! tdf
    REAL(fp)          :: SO4, SO40, SO4s, SO40s
    REAL(fp)          :: PSO4d, SO40_dust ! tdf 04/07/08

    ! Arrays
    REAL(fp)          :: SO4d (NDSTBIN)  ! tdf 04/07/08
    REAL(fp)          :: SO40d(NDSTBIN)  ! tdf 04/07/08

    ! Following are index arrays to hold pointers to STT
    ! tdf 04/07/08
    INTEGER           :: IDTRC(NDSTBIN)

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)
#ifdef APM
    REAL(fp), POINTER :: PSO4_SO2APM2(:,:,:)
#endif

#ifdef APM
    REAL*8            :: MASS0, MASS, PMASS
    REAL*8            :: RKTs, E_RKTs, DTCHEM
#endif

    !=================================================================
    ! CHEM_SO4 begins here!
    !=================================================================

    ! Return if tracers are not defined
    IF ( id_SO4 < 0 .or. id_SO4s < 0 ) RETURN

    ! Assume success
    RC        = GC_SUCCESS

    ! Copy fields from INPUT_OPT to local variables for use below
    LDSTUP    = Input_Opt%LDSTUP

#ifdef APM
    !------------------------------------------
    ! Call APM size-resolved drydep algorithm
    !------------------------------------------
    CALL WET_SETTLINGBIN( Input_Opt, State_Chm, State_Diag, State_Grid, &
                          State_Met, RC )

    ! Point to PSO4_SO2APM2 now moved to State_Met
    PSO4_SO2APM2 => State_Met%PSO4_SO2APM2
#endif

    ! Point to chemical species array [kg]
    Spc => State_Chm%Species

    ! Loop over chemistry grid boxes
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, SO4, SO4s, SO40, SO40s ) &
    !$OMP PRIVATE( SO4d, SO40d, SO40_dust          ) &
    !$OMP PRIVATE( IBIN, PSO4d, IDTRC              ) &
#ifdef APM
    !$OMP PRIVATE( N, MASS0, MASS, PMASS, RKTs, E_RKTs ) &
#endif
    !$OMP SCHEDULE( DYNAMIC )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Skip non-chemistry boxes
       IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

       ! Initialize for safety's sake
       SO4  = 0e+0_fp
       SO4s = 0e+0_fp
       SO4d = 0e+0_fp  ! tdf 04/07/08

       !==============================================================
       ! Initial concentrations before chemistry
       !==============================================================

       ! SO4 [v/v]
       SO40  = Spc(id_SO4)%Conc(I,J,L)

       ! SO4 within coarse seasalt aerosol [v/v]
       SO40s = Spc(id_SO4s)%Conc(I,J,L)

       IF ( LDSTUP ) THEN
          ! Initial Sulfate w/in dust bins [v/v]   ! tdf 04/07/08
          SO40d(1) = Spc(id_SO4d1)%Conc(I,J,L)
          SO40d(2) = Spc(id_SO4d2)%Conc(I,J,L)
          SO40d(3) = Spc(id_SO4d3)%Conc(I,J,L)
          SO40d(4) = Spc(id_SO4d4)%Conc(I,J,L)
       ENDIF

       !==============================================================
       ! SO4 chemistry
       !==============================================================

       ! SO4 production from SO2 and HMS [v/v/timestep]
       SO4 = SO40 + PSO4_SO2(I,J,L)

       !==============================================================
       ! SO4s (SO4 w/in seasalt aerosol) chemistry:
       !==============================================================

       ! SO4 production from SO2 [v/v/timestep]
       SO4s = SO40s + PSO4_ss(I,J,L)

       !tdf
       IF ( LDSTUP ) THEN

          !==============================================================
          ! SO4d (SO4 w/in dust aerosol) chemistry:     tdf 04/07/08
          !==============================================================

          IDTRC(1) = id_SO4d1
          IDTRC(2) = id_SO4d2
          IDTRC(3) = id_SO4d3
          IDTRC(4) = id_SO4d4

          ! tdf Loop over size bins
          DO IBIN = 1, NDSTBIN

             ! Initial amount of sulfate on dust size bin IBIN
             SO40_dust = SO40d(IBIN)

             ! Production of sulfate on dust [v/v/timestep]
             PSO4d = PSO4_dust(I,J,L,IBIN)

             ! SO4 production from SO2 [v/v/timestep]
             SO4d(IBIN) = SO40_dust + PSO4d

          ENDDO
       ENDIF    !tdf end of if ( LDSTUP) condition

       !==============================================================
       ! Final concentrations after chemistry
       !==============================================================

       ! Error check
       IF ( SO4  < SMALLNUM ) SO4  = 0e+0_fp
       IF ( SO4s < SMALLNUM ) SO4s = 0e+0_fp

       ! Final concentrations [v/v]
       Spc(id_SO4)%Conc(I,J,L)  = SO4
       Spc(id_SO4s)%Conc(I,J,L) = SO4s

!APM_GanLuo+
#ifdef APM
       IF(NSO4>=1)THEN
          DO N=1,NSO4
             ! Updated SO4 (gas phase) [v/v]
             Spc(APMIDS%id_SO4BIN1+N-1)%Conc(I,J,L) = &
                  Spc(APMIDS%id_SO4BIN1+N-1)%Conc(I,J,L) + &
                  (PSO4_SO2APM(I,J,L)+PSO4_SO2APM2(I,J,L)*(AIRMW/96.D0)/ &
                  (g0_100*State_Met%DELP_DRY(I,J,L)))* &
                  FCLOUD(I,J,L,N)
          ENDDO
       ENDIF

       IF(NCTBC>=1)THEN
          DO N=1,1
             Spc(APMIDS%id_CTBC+N-1)%Conc(I,J,L) = &
                  Spc(APMIDS%id_CTBC+N-1)%Conc(I,J,L) + &
                  (PSO4_SO2APM(I,J,L)+PSO4_SO2APM2(I,J,L)*(AIRMW/96.D0)/ &
                  (g0_100*State_Met%DELP_DRY(I,J,L)))* &
                  FCLOUD(I,J,L,(NSO4+N))
          ENDDO
       ENDIF

       IF(NCTOC>=1)THEN
          DO N=1,1
             Spc(APMIDS%id_CTOC+N-1)%Conc(I,J,L) = &
                  Spc(APMIDS%id_CTOC+N-1)%Conc(I,J,L) + &
                  (PSO4_SO2APM(I,J,L)+PSO4_SO2APM2(I,J,L)*(AIRMW/96.D0)/ &
                  (g0_100*State_Met%DELP_DRY(I,J,L)))* &
                  FCLOUD(I,J,L,(NSO4+N))
          ENDDO
       ENDIF

       IF(NCTDST>=1)THEN
          DO N=1,1
             Spc(APMIDS%id_CTDST+N-1)%Conc(I,J,L) = &
                  Spc(APMIDS%id_CTDST+N-1)%Conc(I,J,L) + &
                  (PSO4_SO2APM(I,J,L)+PSO4_SO2APM2(I,J,L)*(AIRMW/96.D0)/ &
                  (g0_100*State_Met%DELP_DRY(I,J,L)))* &
                  FCLOUD(I,J,L,(NSO4+3))
          ENDDO
       ENDIF

       IF(NCTSEA>=1)THEN
          DO N=1,1
             Spc(APMIDS%id_CTSEA+N-1)%Conc(I,J,L) = &
                  Spc(APMIDS%id_CTSEA+N-1)%Conc(I,J,L) + &
                  (PSO4_SO2APM(I,J,L)+PSO4_SO2APM2(I,J,L)*(AIRMW/96.D0)/ &
                  (g0_100*State_Met%DELP_DRY(I,J,L)))* &
                  FCLOUD(I,J,L,(NSO4+4)) + PSO4_SO2SEA(I,J,L)
          ENDDO
       ENDIF
#endif

       !tdf
       IF ( LDSTUP ) THEN
          DO IBIN = 1, NDSTBIN
             SO40_dust = SO4d(IBIN)
             IF ( SO40_dust < SMALLNUM ) SO40_dust = 0e+0_fp
             Spc(IDTRC(IBIN))%Conc(I,J,L) = SO40_dust ! dust sulfate
          ENDDO
       ENDIF    !tdf end of if ( LDSTUP) condition

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

#ifdef APM
    PSO4_SO2APM2 = 0.D0
#endif

    ! Free pointers
    Spc => NULL()

  END SUBROUTINE CHEM_SO4
!EOC
#ifdef TOMAS
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_so4_aq
!
! !DESCRIPTION: Subroutine CHEM\_SO4\_AQ takes the SO4 produced via aqueous
!  chemistry of SO2 and distribute onto the size-resolved aerosol number and
!  sulfate mass as a part of the TOMAS aerosol microphysics module
!  (win, 1/25/10)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_SO4_AQ( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD
    USE Input_Opt_Mod,      ONLY : OptINput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TOMAS_MOD,          ONLY : AQOXID, GETACTBIN
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  NOTE: This subroutine is ignored unless we compile for TOMAS microphysics.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER           :: I, J, L
    INTEGER           :: k, binact1, binact2
    INTEGER           :: KMIN
    REAL(fp)          :: SO4OXID
    CHARACTER(LEN=63) :: OrigUnit

    !=================================================================
    ! CHEM_SO4_AQ begins here!
    !=================================================================

    ! Assume success
    RC  = GC_SUCCESS

    ! Convert species from to [kg]
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'kg', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'Start of CHEM_SO4_AQ in sulfate_mod.F90')
       RETURN
    ENDIF

    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED )  &
    !$OMP PRIVATE( I, J, L ) &
    !$OMP PRIVATE( KMIN, SO4OXID, BINACT1, BINACT2 ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Skip non-chemistry boxes
       IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

       SO4OXID = PSO4_SO2AQ(I,J,L) * State_Met%AD(I,J,L) &
                 / ( AIRMW / State_Chm%SpcData(id_SO4)%Info%MW_g )

       IF ( SO4OXID > 0e+0_fp ) THEN
          ! JKodros (6/2/15 - Set activating bin based on which TOMAS bin
          !length being used)
#if defined( TOMAS12 )
          CALL GETACTBIN( I, J, L, id_NK5, .TRUE. , BINACT1, State_Chm, RC )

          CALL GETACTBIN( I, J, L, id_NK5, .FALSE., BINACT2, State_Chm, RC )
#elif defined( TOMAS15 )
          CALL GETACTBIN( I, J, L, id_NK8, .TRUE. , BINACT1, State_Chm, RC )

          CALL GETACTBIN( I, J, L, id_NK8, .FALSE., BINACT2, State_Chm, RC )
#elif defined( TOMAS30 )
          CALL GETACTBIN( I, J, L, id_NK10, .TRUE. , BINACT1, State_Chm, RC )

          CALL GETACTBIN( I, J, L, id_NK10, .FALSE., BINACT2, State_Chm, RC )
#else
          CALL GETACTBIN( I, J, L, id_NK20, .TRUE. , BINACT1, State_Chm, RC )

          CALL GETACTBIN( I, J, L, id_NK20, .FALSE., BINACT2, State_Chm, RC )
#endif

          KMIN = ( BINACT1 + BINACT2 )/ 2.

          CALL AQOXID( SO4OXID, KMIN, I, J, L, Input_Opt, &
                       State_Chm, State_Grid, State_Met, RC )
       ENDIF
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Convert species back to original units
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, &
                     'End of CHEM_SO4_AQ in sulfate_mod.F90')
       RETURN
    ENDIF

  END SUBROUTINE CHEM_SO4_AQ
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_msa
!
! !DESCRIPTION: Subroutine CHEM\_MSA is the SO4 chemistry subroutine from Mian
!  Chin's GOCART model, modified for the GEOS-CHEM model. (rjp, bdf, bmy,
!  5/31/00, 10/25/05)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_MSA( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : SpcConc
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  The only production is from DMS oxidation (saved in CHEM_DMS).
!  Dry deposition is now treaded in mixing_mod.F90.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER           :: I,    J,    L
    REAL(fp)          :: MSA0, MSA

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    !=================================================================
    ! CHEM_MSA begins here!
    !=================================================================
    IF ( id_MSA < 0 ) RETURN

    ! Assume success
    RC  = GC_SUCCESS

    ! Point to chemical species array [v/v/ dry]
    Spc => State_Chm%Species

    ! Loop over chemistry grid boxes
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, MSA0, MSA ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Skip non-chemistry boxes
       IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

       ! Initial MSA [v/v]
       MSA0 = Spc(id_MSA)%Conc(I,J,L)

       ! MSA production from DMS [v/v/timestep]
       MSA = MSA0 + PMSA_DMS(I,J,L)

       ! Final MSA [v/v]
       IF ( MSA < SMALLNUM ) MSA = 0e+0_fp
       Spc(id_MSA)%Conc(I,J,L) = MSA

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE CHEM_MSA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_nit
!
! !DESCRIPTION: Subroutine CHEM\_NIT removes SULFUR NITRATES (NIT) from the
!  surface via dry deposition. (rjp, bdf, bmy, 1/2/02, 5/23/06)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_NIT( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
    USE CMN_SIZE_Mod,       ONLY : NDSTBIN
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : SpcConc
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?

!
! !REMARKS:
!  Dry deposition is now applied in mixing_mod.F90.  Therefore we can
!                                                                             .
! !REMARKS:
!  Reaction List:
!  ============================================================================
!  (1 ) NIT = NIT_0 * EXP( -dt )  where d = dry deposition rate [s-1]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER           :: I, IBIN, J, L, N
    REAL(fp)          :: PNITd

    ! Arrays
    INTEGER           :: IDTRC(NDSTBIN)
    REAL(fp)          :: NITd (NDSTBIN)
    REAL(fp)          :: NIT0d(NDSTBIN)

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    !=================================================================
    ! CHEM_NIT begins here!
    !=================================================================

    ! Return if tracers are not defined
    IF ( id_NIT < 0 .OR.id_NITs < 0 ) RETURN

    ! If we are doing a dust uptake simulation,
    ! return if dust uptake tracers are undefined
    IF ( Input_Opt%LDSTUP ) THEN
       IF ( id_NITd1 < 0 .or. id_NITd2 < 0   .or. &
            id_NITd3 < 0 .or. id_NITd4 < 0 ) RETURN

    ENDIF

    ! Assume success
    RC       = GC_SUCCESS

    ! Point to chemical species array [v/v dry]
    Spc      => State_Chm%Species

    ! Assign pointers to Spc arrays for loop over dust size bins,
    ! since this can be done outside the parallel DO loop.
    ! These will be set to a missing value (-1) if undefined.
    IDTRC(1) = id_NITd1
    IDTRC(2) = id_NITd2
    IDTRC(3) = id_NITd3
    IDTRC(4) = id_NITd4

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, NITd, NIT0d, IBIN, PNITd ) &
    !$OMP SCHEDULE( DYNAMIC, 1 )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Skip non-chemistry boxes
       IF ( .not. State_Met%InChemGrid(I,J,L) ) CYCLE

       !==============================================================
       ! NITs chemistry
       !==============================================================

       ! Add NIT prod from HNO3 uptake on fine sea-salt [v/v/timestep]
       ! to the NITs concentration in the Spc array [v/v dry]
       ! Should do for both fine and coarse mode (xnw 12/8/17)
       Spc(id_NIT)%Conc(I,J,L)  = Spc(id_NIT)%Conc(I,J,L) + PNIT(I,J,L)
       Spc(id_NITs)%Conc(I,J,L) = Spc(id_NITs)%Conc(I,J,L) + PNITs(I,J,L)

       !==============================================================
       ! NITd chemistry (cf. Duncan Fairlie)
       !==============================================================
       IF ( Input_Opt%LDSTUP ) THEN

          ! Initialize variables
          NITd     = 0.0_fp

          ! Initial NITRATE w/in dust bins [v/v]
          NIT0d(1) = Spc(id_NITd1)%Conc(I,J,L)
          NIT0d(2) = Spc(id_NITd2)%Conc(I,J,L)
          NIT0d(3) = Spc(id_NITd3)%Conc(I,J,L)
          NIT0d(4) = Spc(id_NITd4)%Conc(I,J,L)

          ! Loop over size bins
          DO IBIN = 1, NDSTBIN

             ! Production of nitrate on dust [v/v/timestep]
             PNITd                  = PNIT_dust(I,J,L,IBIN)

             ! NIT prod from HNO3 uptake on dust [v/v/timestep]
             NITd(IBIN)             = NIT0d(IBIN) + PNITd

             ! Store final concentration in Spc [v/v]
             Spc(IDTRC(IBIN))%Conc(I,J,L) = NITd(IBIN)

          ENDDO

       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE CHEM_NIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_cl
!
! !DESCRIPTION: Subroutine CHEM\_CL cacluate sea-salt chloride from
!  uptake of HCl by alkalinity. (xnw, 12/8/17)
!\\
!\\
! !INTERFACE:

      SUBROUTINE CHEM_CL( Input_Opt, &
           State_Met, State_Chm, State_Grid, RC )

! !USES:

#ifdef BPCH_DIAG
        USE CMN_DIAG_MOD
#endif

      USE CMN_SIZE_MOD
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE Species_Mod,        ONLY : SpcConc
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE State_Grid_Mod,     ONLY : GrdState
! !INPUT PARAMETERS:

      !LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      TYPE(GrdState), INTENT(IN)    :: State_Grid

! !INPUT/OUTPUT PARAMETERS:

      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!-----------------------------------------------------------------------------
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER           :: I, J, L

      ! Pointers
      TYPE(SpcConc), POINTER :: Spc(:)

      !=================================================================
      ! CHEM_NIT begins here!
      !=================================================================

      ! Return if tracers are not defined
      IF ( id_SALACL < 0 .OR.id_SALCCL < 0 ) RETURN

      ! Assume success
      RC       = GC_SUCCESS

      ! Point to chemical species array [v/v dry]
      Spc      => State_Chm%Species

      !$OMP PARALLEL DO          &
      !$OMP DEFAULT( SHARED )    &
      !$OMP PRIVATE( I, J, L)    &
      !$OMP SCHEDULE( DYNAMIC, 1 )
      DO L = 1, State_Grid%MaxChemLev
      DO J = 1, State_Grid%NY
      DO I = 1, State_Grid%NX

         ! Skip non-chemistry boxes
!         IF ( ITS_IN_THE_NOCHEMGRID( I, J, L, State_Met ) ) CYCLE

         !==============================================================
         ! Cl chemistry
         !==============================================================

         ! Add Cl prod from HCl uptake on fine and coarse sea-salt [v/v/timestep]
         Spc(id_SALACL)%Conc(I,J,L)  = Spc(id_SALACL)%Conc(I,J,L) + PACL(I,J,L)
         Spc(id_SALCCL)%Conc(I,J,L)  = Spc(id_SALCCL)%Conc(I,J,L) + PCCL(I,J,L)

      ENDDO
      ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      ! Free pointer
      Spc => NULL()

      END SUBROUTINE CHEM_CL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_oh
!
! !DESCRIPTION: Function GET\_OH returns OH from State\_Chm%Species (for
!  coupled runs) or monthly mean OH (for offline runs).  Imposes a diurnal
!  variation on OH for offline simulations. (bmy, 12/16/02, 7/20/04)
!\\
!\\
! !INTERFACE:
!
  FUNCTION GET_OH( I, J, L, Input_Opt, State_Chm, State_Met ) &
       RESULT( OH_MOLEC_CM3 )
!
! !USES:
!
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: I, J, L     ! Lon, lat, level indices
    TYPE(OptInput), INTENT(IN) :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
!
! !RETURN VALUE:
!
    REAL(fp)                   :: OH_MOLEC_CM3 ! OH conc [molec/cm3]
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! GET_OH begins here!
    !=================================================================

    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       !---------------------
       ! Coupled simulation
       !---------------------

       ! OH is defined only in the chemistry grid
       IF ( State_Met%InChemGrid(I,J,L) ) THEN

          ! Get OH from State_Chm%Species [v/v] converted to [molec/cm3]
          OH_MOLEC_CM3 = State_Chm%Species(id_OH)%Conc(I,J,L) * &
                         State_Met%AIRNUMDEN(I,J,L)
       ELSE
          OH_MOLEC_CM3 = 0e+0_fp
       ENDIF

    ELSE IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       !---------------------
       ! Offline simulation
       !---------------------

       ! Test for sunlight...
       IF ( State_Met%SUNCOS(I,J) > 0e+0_fp .and. TCOSZ(I,J) > 0e+0_fp ) THEN

          ! OH from HEMCO is in mol/mol, convert to molec/cm3
          OH_MOLEC_CM3 = GLOBAL_OH(I,J,L) * State_Met%AIRNUMDEN(I,J,L)

          ! Impose a diurnal variation on OH during the day
          OH_MOLEC_CM3 = OH_MOLEC_CM3 * &
                         ( State_Met%SUNCOS(I,J) / TCOSZ(I,J) ) * &
                         ( 86400e+0_fp           / GET_TS_CHEM() )

          ! Make sure OH is not negative
          OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0e+0_fp )

       ELSE

          ! At night, OH goes to zero
          OH_MOLEC_CM3 = 0e+0_fp

       ENDIF

    ENDIF

  END FUNCTION GET_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ohno3time
!
! !DESCRIPTION: Subroutine OHNO3TIME computes the sum of cosine of the solar
!  zenith angle over a 24 hour day, as well as the total length of daylight.
!  This is needed to scale the offline OH and NO3 concentrations.
!  (rjp, bmy, 12/16/02, 3/30/04)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE OHNO3TIME( State_Grid )
!
! !USES:
!
    USE TIME_MOD,       ONLY : GET_NHMSb,   GET_ELAPSED_SEC
    USE TIME_MOD,       ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE       :: FIRST = .TRUE.
    INTEGER             :: I, J, L, N, NT, NDYSTEP
    REAL(fp)            :: A0, A1, A2, A3, B1, B2, B3
    REAL(fp)            :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
    REAL(fp)            :: SUNTMP(State_Grid%NX,State_Grid%NY)

    !=================================================================
    ! OHNO3TIME begins here!
    !=================================================================

    !  Solar declination angle (low precision formula, good enough for us):
    A0 = 0.006918
    A1 = 0.399912
    A2 = 0.006758
    A3 = 0.002697
    B1 = 0.070257
    B2 = 0.000907
    B3 = 0.000148
    R  = 2.* PI * float( GET_DAY_OF_YEAR() - 1 ) / 365.

    DEC = A0 - A1*cos(  R) + B1*sin(  R) &
             - A2*cos(2*R) + B2*sin(2*R) &
             - A3*cos(3*R) + B3*sin(3*R)

    LHR0 = int(float( GET_NHMSb() )/10000.)

    ! Only do the following at the start of a new day
    IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN

       ! Zero arrays
       TTDAY(:,:) = 0e+0_fp
       TCOSZ(:,:) = 0e+0_fp
       COSZM(:,:) = 0e+0_fp

       ! NDYSTEP is # of chemistry time steps in this day
       NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 3600 / GET_TS_CHEM()

       ! NT is the elapsed time [s] since the beginning of the run
       NT = GET_ELAPSED_SEC()

       ! Loop forward through NDYSTEP "fake" timesteps for this day
       DO N = 1, NDYSTEP

          ! Zero SUNTMP array
          SUNTMP = 0e+0_fp

          ! Loop over surface grid boxes
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Grid box latitude center [radians]
             YMID_R = State_Grid%YMid_R(I,J)

             TIMLOC = real(LHR0) + real(NT)/3600.0 + State_Grid%XMid(I,J)/15.0

             DO WHILE (TIMLOC .lt. 0)
                TIMLOC = TIMLOC + 24.0
             ENDDO

             DO WHILE (TIMLOC .gt. 24.0)
                TIMLOC = TIMLOC - 24.0
             ENDDO

             AHR = abs(TIMLOC - 12.) * 15.0 * PI_180

             !===========================================================
             ! The cosine of the solar zenith angle (SZA) is given by:
             !
             !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)
             !
             ! where LAT = the latitude angle,
             !       DEC = the solar declination angle,
             !       AHR = the hour angle, all in radians.
             !
             ! If SUNCOS < 0, then the sun is below the horizon, and
             ! therefore does not contribute to any solar heating.
             !===========================================================

             ! Compute Cos(SZA)
             SUNTMP(I,J) = sin(YMID_R) * sin(DEC) + &
                           cos(YMID_R) * cos(DEC) * cos(AHR)

             ! TCOSZ is the sum of SUNTMP at location (I,J)
             ! Do not include negative values of SUNTMP
             TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(I,J), 0e+0_fp )

             ! COSZM is the peak value of SUMTMP during a day at (I,J)
             ! (rjp, bmy, 3/30/04)
             COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(I,J) )

             ! TTDAY is the total daylight time at location (I,J)
             IF ( SUNTMP(I,J) > 0e+0_fp ) THEN
                TTDAY(I,J) = TTDAY(I,J) + DBLE( GET_TS_CHEM() ) * 60e+0_fp
             ENDIF
          ENDDO
          ENDDO

          !### Debug
          !PRINT*, '### IN OHNO3TIME'
          !PRINT*, '### N       : ', N
          !PRINT*, '### NDYSTEP : ', NDYSTEP
          !PRINT*, '### NT      : ', NT
          !PRINT*, '### JDAY    : ', JDAY
          !PRINT*, '### RLAT    : ', RLAT
          !PRINT*, '### XMID    : ', XMID
          !PRINT*, '### SUNTMP  : ', SUNTMP
          !PRINT*, '### TCOSZ   : ', MINVAL( TCOSZ ), MAXVAL( TCOSZ )
          !PRINT*, '### TTDAY   : ', MINVAL( TCOSZ ), MAXVAL( TCOSZ )

          ! Increment elapsed time [sec]
          NT = NT + GET_TS_CHEM()
       ENDDO

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

  END SUBROUTINE OHNO3TIME
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_alk
!
! !DESCRIPTION: Subroutine GET\_ALK returns the seasalt alkalinity emitted at
!  each timestep to sulfate\_mod.F90 for chemistry on seasalt aerosols.
!  (bec, 12/7/04, 11/23/09)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GET_ALK( I, J, L, ALK1, ALK2, Kt1, Kt2, Kt1N, Kt2N, Kt1L, Kt2L, &
                      Input_Opt, State_Chm, State_Grid, State_Met, RC)
!
! !USES:
!
    USE CMN_SIZE_MOD,         ONLY : NDUST
    USE ErrCode_Mod
    USE HCO_Utilities_GC_Mod, ONLY : GetHcoValEmis
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L     ! Lon-lat-alt indices
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT)   :: ALK1, ALK2  ! [kg]
    REAL(fp),       INTENT(OUT)   :: Kt1, Kt2, Kt1N, Kt2N, Kt1L, Kt2L ! [s-1]
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(INOUT) :: RC
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp)            :: N1, N2, Kt
    REAL(fp)            :: HGF, ALK
    REAL(fp)            :: RAD1, RAD2, RAD3
    REAL(fp)            :: term1a, term2a, term3a
    REAL(fp)            :: term1b, term2b, term3b
    REAL(fp)            :: term1aN, term2aN, term3aN
    REAL(fp)            :: term1bN, term2bN, term3bN
    REAL(fp)            :: const1, const2, const1N, const2N
    REAL(fp)            :: const1L, const2L
    REAL(fp)            :: a1, a2, b1, b2, a1N, a2N, b1N, b2N
    REAL(fp)            :: SLA, SLC !SeaSalt concentration [kg]
    REAL(fp), PARAMETER :: MINDAT = 1.e-20_fp
    INTEGER             :: IRH
    REAL(fp), PARAMETER :: gamma_SO2 = 0.11e+0_fp !from Worsnop et al.(1989)
    REAL(fp), PARAMETER :: gamma_HNO3 = 0.5e+0_fp !from JPL [2015]
    ! Mass accommodation coefficient alpha for HCl on liquid water surfaces, from JPL,2015
    REAL(fp), PARAMETER :: gamma_HCl = 0.07e+0_fp
    REAL(fp)            :: Dg    !gas phase diffusion coeff. [cm2/s]
    REAL(fp), PARAMETER :: v = 3.0e+4_fp  !cm/s
    REAL(fp)            :: SA1, SA2 ! SSA surface area [cm2/cm3]
    REAL(fp)            :: R1, R2 ! SSA radius [cm]

    ! HEMCO update
    REAL(fp)            :: FEMIS
    !INTEGER             :: NTOP
    !LOGICAL             :: FOUND

    ! Pointers
    REAL(fp), POINTER   :: ERADIUS(:,:,:,:)
    REAL(fp), POINTER   :: TAREA(:,:,:,:)

    !=================================================================
    ! GET_ALK begins here!
    !=================================================================

    ! Zero variables
    KT1   = 0.e+0_fp
    KT2   = 0.e+0_fp
    KT1N  = 0.e+0_fp
    KT2N  = 0.e+0_fp
    KT1L  = 0.e+0_fp
    KT2L  = 0.e+0_fp
    N1    = 0.e+0_fp
    N2    = 0.e+0_fp
    ALK1  = 0.e+0_fp
    ALK2  = 0.e+0_fp

    ! Initialize pointers
    ERADIUS => State_Chm%AeroRadi
    TAREA   => State_Chm%AeroArea

    !-----------------------------------------------------------------------
    ! Get alkalinity from HEMCO. This is just the current emissions,
    ! converted from kg/m2/s to kg. In the original seasalt code, the
    ! alkalinity was set to the total surface flux for every layer within
    ! the PBL, and to zero above. This seems unrealistic, and in the code
    ! below the alkalinity (and number density) are scaled by the fraction
    ! of PBL. This approach ensures that the total number density and
    ! alkalinity is preserved.
    !                                                 (ckeller, 10/31/2014)
    !-----------------------------------------------------------------------
    ! [kg] use this when not transporting alk
    !ALK1  = ALK_EMIS(I,J,L,1)
    !ALK2  = ALK_EMIS(I,J,L,2)

    ! Layer in which the PBL top occurs
    !NTOP = CEILING( State_Met%PBL_TOP_L(I,J) )

    ! Do the following only if we are within the PBL
    !IF ( L <= NTOP ) THEN

       ! Fraction of the PBL spanned by box (I,J,L) [unitless]
       FEMIS = State_Met%F_OF_PBL(I,J,L)

       ! Uncomment the following line to reproduce pre-HEMCO code.
       ! FEMIS = 1.0d0

       ! Get ALK1 and ALK2 surface emissions from HEMCO. These are in
       ! kg/m2/s.
       !CALL GetHcoValEmis( Input_Opt, State_Grid, id_SALA, I, J, 1, FOUND, ALK1 )
       !CALL GetHcoValEmis( Input_Opt, State_Grid, id_SALC, I, J, 1, FOUND, ALK2, AltBuffer=.true. )

       ! kg/m2/s --> kg. Weight by fraction of PBL
       ALK1 = MAX(ALK1,0.0e+0_fp) * State_Grid%Area_M2(I,J) * TS_EMIS * FEMIS
       ALK2 = MAX(ALK2,0.0e+0_fp) * State_Grid%Area_M2(I,J) * TS_EMIS * FEMIS

       ! Get number densities in # / cm3. Weight by fraction of PBL
       !IF ( ASSOCIATED(NDENS_SALA) ) THEN
       !   N1 = NDENS_SALA(I,J) / State_Met%AIRVOL(I,J,L) &
       !      * 1.e-6_fp * FEMIS
       !ENDIF
       !IF ( ASSOCIATED(NDENS_SALC) ) THEN
       !   N2 = NDENS_SALC(I,J) / State_Met%AIRVOL(I,J,L) &
       !      * 1.e-6_fp * FEMIS
       !ENDIF

    !ENDIF

    !-----------------------------------------------------------------------
    ! NOTE: If you want to transport alkalinity then uncomment this section
    ! (bec, bmy, 4/13/05)
    !
    !! alkalinity [v/v] to [kg] use this when transporting alk
    !! or using Liao et al [2004] assumption of a continuous supply of
    ! alkalinity based on Laskin et al. [2003]
    !ALK1 = Spc(id_SALA)%Conc(I,J,L) * State_Met%AD(I,J,L)/
    !  & ( AIRMW / State_Chm%SpcData(id_SALA)%Info%MW_g )
    !ALK2 = Spc(id_SALC)%Conc(I,J,L) * State_Met%AD(I,J,L)/
    !  & ( AIRMW / State_Chm%SpcData(id_SALC)%Info%MW_g )
    !-----------------------------------------------------------------------

    ! Conversion from [m-3] --> [cm-3]
    !N1 = N_DENS(I,J,L,1) * 1.d-6
    !N2 = N_DENS(I,J,L,2) * 1.d-6
       !Read Alkalinity from Alkalinity tracers [v/v] to [kg], xnw 12/8/17
    ALK1 = State_Chm%Species(id_SALAAL)%Conc(I,J,L) * State_Met%AD(I,J,L)/ &
         ( AIRMW / State_Chm%SpcData(id_SALAAL)%Info%MW_g )
    ALK2 = State_Chm%Species(id_SALCAL)%Conc(I,J,L) * State_Met%AD(I,J,L)/ &
      ( AIRMW / State_Chm%SpcData(id_SALCAL)%Info%MW_g )
    !Seasalt mass, [v/v] to [kg]
    !SLA = State_Chm%Species(id_SALA)%Conc(I,J,L) * State_Met%AD(I,J,L)/
    !     & ( AIRMW / State_Chm%SpcData(id_SALA)%Info%MW_g )
    !SLC = State_Chm%Species(id_SALC)%Conc(I,J,L) * State_Met%AD(I,J,L)/
    !     & ( AIRMW / State_Chm%SpcData(id_SALC)%Info%MW_g )


    ALK = ALK1 + ALK2

    ! If there is any alkalinity ...
    IF ( ALK > MINDAT ) THEN

       SA1= TAREA(I,J,L,4+NDUST) !in cm2/cm3
       SA2= TAREA(I,J,L,5+NDUST) !in cm2/cm3

       R1 = ERADIUS(I,J,L,4+NDUST) !in cm
       R2 = ERADIUS(I,J,L,5+NDUST) !in cm

       ! set humidity index IRH as a percent
       !IRH = State_Met%RH(I,J,L)
       !IRH = MAX(  1, IRH )
       !IRH = MIN( 99, IRH )

       ! Hygroscopic growth factor for sea-salt from Chin et al. (2002)
       ! Updated (bec, bmy, 11/23/09)
       !IF ( IRH < 100 ) HGF = 4.8e+0_fp
       !IF ( IRH < 99  ) HGF = 2.9e+0_fp
       !IF ( IRH < 95  ) HGF = 2.4e+0_fp
       !IF ( IRH < 90  ) HGF = 2.0e+0_fp
       !IF ( IRH < 80  ) HGF = 1.8e+0_fp
       !IF ( IRH < 70  ) HGF = 1.6e+0_fp
       !IF ( IRH < 50  ) HGF = 1.0e+0_fp

       ! radius of sea-salt aerosol size bins [cm] accounting for
       ! hygroscopic growth
       !RAD1 = Input_Opt%SALA_REDGE_um(1) * HGF * 1.e-4_fp
       !RAD2 = Input_Opt%SALA_REDGE_um(2) * HGF * 1.e-4_fp
       !RAD3 = Input_Opt%SALC_REDGE_um(2) * HGF * 1.e-4_fp

       !----------------------------------
       ! SO2 uptake onto fine particles
       !----------------------------------
        DG = 9.45E+17_fp/State_Met%AIRNUMDEN(I,J,L) * &
             SQRT(State_Met%T(I,J,L)) * SQRT(3.472E-2_fp +1.E+0_fp/64.e+0_fp)

       ! calculate gas-to-particle rate constant for uptake of
       ! SO2 onto fine sea-salt aerosols [Jacob, 2000] analytical solution
       CONST1 = 4.e+0_fp/(V*GAMMA_SO2)
       !A1     = (RAD1/DG)+CONST1
       !B1     = (RAD2/DG)+CONST1
       !TERM1A = ((B1**2)/2.0e+0_fp) - ((A1**2)/2.0e+0_fp)
       !TERM2A = 2.e+0_fp*CONST1*(B1-A1)
       !TERM3A = (CONST1**2)*LOG(B1/A1)
       !KT1    = 4.e+0_fp*PI*N1*(DG**3)*(TERM1A - TERM2A + TERM3A)
       ! now calculate rate of uptake as kt = [SSA]/( r/Dg + 4/(v*gamma) )
       KT1    = SA1/((R1/DG) + CONST1)

       !----------------------------------
       ! SO2 uptake onto coarse particles
       !----------------------------------

       ! calculate gas-to-particle rate constant for uptake of
       ! SO2 onto coarse sea-salt aerosols [Jacob, 2000] analytical solution
       CONST2 = 4.e+0_fp/(V*GAMMA_SO2)
       !A2     = (RAD2/DG)+CONST2
       !B2     = (RAD3/DG)+CONST2
       !TERM1B = ((B2**2)/2.0e+0_fp) - ((A2**2)/2.0e+0_fp)
       !TERM2B = 2.e+0_fp*CONST2*(B2-A2)
       !TERM3B = (CONST2**2)*LOG(B2/A2)
       !KT2    = 4.e+0_fp*PI*N2*(DG**3)*(TERM1B - TERM2B + TERM3B)
       KT2    = SA2/((R2/DG) + CONST2)
       KT     = KT1 + KT2

       !----------------------------------
       ! HNO3 uptake onto fine particles
       !----------------------------------
       DG = 9.45E+17_fp/State_Met%AIRNUMDEN(I,J,L) * &
            SQRT(State_Met%T(I,J,L)) * SQRT(3.472E-2_fp +1.E+0_fp/63.e+0_fp)

       ! calculate gas-to-particle rate constant for uptake of
       ! HNO3 onto fine sea-salt aerosols [Jacob, 2000] analytical solution
       CONST1N = 4.e+0_fp/(V*GAMMA_HNO3)
       !A1N     = (RAD1/DG)+CONST1N
       !B1N     = (RAD2/DG)+CONST1N
       !TERM1AN = ((B1N**2)/2.0e+0_fp) - ((A1N**2)/2.0e+0_fp)
       !TERM2AN = 2.e+0_fp*CONST1N*(B1N-A1N)
       !TERM3AN = (CONST1N**2)*LOG(B1N/A1N)
       !KT1N    = 4.e+0_fp*PI*N1*(DG**3)*(TERM1AN - TERM2AN + TERM3AN)
       KT1N    = SA1/((R1/DG) + CONST1N)

       !----------------------------------
       ! HNO3 uptake onto coarse particles
       !----------------------------------

       ! calculate gas-to-particle rate constant for uptake of
       ! HNO3 onto coarse sea-salt aerosols [Jacob, 2000] analytical solution
       CONST2N = 4.e+0_fp/(V*GAMMA_HNO3)
       !A2N     = (RAD2/DG)+CONST2N
       !B2N     = (RAD3/DG)+CONST2N
       !TERM1BN = ((B2N**2)/2.0e+0_fp) - ((A2N**2)/2.0e+0_fp)
       !TERM2BN = 2.e+0_fp*CONST2N*(B2N-A2N)
       !TERM3BN = (CONST2N**2)*LOG(B2N/A2N)
       !KT2N    = 4.e+0_fp*PI*N2*(DG**3)*(TERM1BN - TERM2BN + TERM3BN)
       KT2N    = SA2/((R2/DG) + CONST2N)
       !----------------------------------
       ! HCl uptake onto fine particles
       !----------------------------------
       DG = 9.45E+17_fp/State_Met%AIRNUMDEN(I,J,L) * &
            SQRT(State_Met%T(I,J,L))*SQRT(3.472E-2_fp +1.E+0_fp/36.45e+0_fp)
       CONST1L = 4.e+0_fp/(V*GAMMA_HCl)
       KT1L    = SA1/((R1/DG) + CONST1L)
       !----------------------------------
       ! HCl uptake onto coarse particles
       !----------------------------------
       CONST2L = 4.e+0_fp/(V*GAMMA_HCl)
       KT2L    = SA2/((R2/DG) + CONST2L)

    ELSE

       ! If no alkalinity, set everything to zero
       ALK1 = 0.e+0_fp
       ALK2 = 0.e+0_fp
       KT1  = 0.e+0_fp
       KT1N = 0.e+0_fp
       KT2  = 0.e+0_fp
       KT2N = 0.e+0_fp

    ENDIF

    ! Free pointers
    NULLIFY(  ERADIUS, TAREA )

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE GET_ALK
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_sulfate
!
! !DESCRIPTION: Subroutine INIT\_SULFATE initializes and zeros all allocatable
!  arrays declared in "sulfate\_mod.F90" (bmy, 6/2/00, 10/15/09)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_SULFATE( Input_Opt, State_Chm, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE CMN_SIZE_Mod,   ONLY : NDSTBIN
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(IN)  :: State_Diag  ! Diagnostics State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  02 Jun 2000 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, J, N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! INIT_SULFATE begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Init_Sulfate (in module GeosCore/sulfate_mod.F90)'

    ! Exit immediately if this is a dry-run simulation
    IF ( Input_Opt%DryRun ) RETURN

    !=================================================================
    ! Error check
    !=================================================================
    IF ( ( .not. Input_Opt%ITS_A_FULLCHEM_SIM )   .and. &
         ( .not. Input_Opt%ITS_AN_AEROSOL_SIM ) ) THEN
       ErrMsg = 'Invalid simulation for sulfate_mod!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! Allocate arrays
    !=================================================================
    ALLOCATE( SSTEMP( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:SSTEMP', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SSTEMP = 0e+0_fp

    ALLOCATE( DMSo( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:DMSo', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    DMSo = 0e+0_fp

    ALLOCATE( PMSA_DMS( State_Grid%NX, State_Grid%NY, State_Grid%MaxChemLev ), &
              STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:PMSA_DMS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PMSA_DMS = 0e+0_fp

    ALLOCATE( PSO2_DMS( State_Grid%NX, State_Grid%NY, State_Grid%MaxChemLev ), &
              STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:PSO2_DMS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PSO2_DMS = 0e+0_fp

    ALLOCATE( PSO4_SO2( State_Grid%NX, State_Grid%NY, State_Grid%MaxChemLev ), &
              STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:PSO4_SO2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PSO4_SO2 = 0e+0_fp

    ALLOCATE( PHMS_SO2( State_Grid%NX, State_Grid%NY, State_Grid%MaxChemLev ), &
              STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:PHMS_SO2', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PHMS_SO2 = 0e+0_fp

    ALLOCATE( PSO2_HMS( State_Grid%NX, State_Grid%NY, State_Grid%MaxChemLev ), &
              STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:PSO2_HMS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PSO2_HMS = 0e+0_fp

    ALLOCATE( PSO4_HMS( State_Grid%NX, State_Grid%NY, State_Grid%MaxChemLev ), &
              STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:PSO4_HMS', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PSO4_HMS = 0e+0_fp

    ALLOCATE( PSO4_ss( State_Grid%NX, State_Grid%NY, State_Grid%MaxChemLev ), &
              STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:PSO4_ss', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PSO4_ss = 0e+0_fp

    ALLOCATE( PNITs( State_Grid%NX, State_Grid%NY, State_Grid%MaxChemLev ), &
              STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:PNITs', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PNITs = 0e+0_fp
    !xnw
    ALLOCATE( PNIT( State_Grid%NX, State_Grid%NY, State_Grid%MaxChemLev), STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F:PNIT', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PNIT = 0e+0_fp
    !xnw
    ALLOCATE( PACL( State_Grid%NX, State_Grid%NY, State_Grid%MaxChemLev ), STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F:PACL', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PACL = 0e+0_fp
    !xnw
    ALLOCATE( PCCL( State_Grid%NX, State_Grid%NY, State_Grid%MaxChemLev ), STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F:PCCL', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PCCL = 0e+0_fp

    !tdf
    ALLOCATE( PSO4_dust( State_Grid%NX, State_Grid%NY, State_Grid%MaxChemLev, &
                         NDSTBIN ), STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:PSO4_dust', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PSO4_dust = 0e+0_fp

    !tdf
    ALLOCATE( PNIT_dust( State_Grid%NX, State_Grid%NY, State_Grid%MaxChemLev, &
                         NDSTBIN ), STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:PNIT_dust', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PNIT_dust = 0e+0_fp

    ALLOCATE( SOx_SCALE( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:SOx_SCALE', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SOx_SCALE = 0e+0_fp

#ifdef APM
    ! Allocate for APM microphysics
    ALLOCATE( PSO4_SO2APM( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:PSO4_SO2APM', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PSO4_SO2APM = 0e+0_fp
    ALLOCATE( PSO4_SO2SEA( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
                           STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:PSO4_SO2SEA', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PSO4_SO2SEA = 0e+0_fp
#endif

#ifdef TOMAS
    ! Allocate for TOMAS microphysics
    ALLOCATE( PSO4_SO2AQ( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:PSO4_SO2aq', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    PSO4_SO2AQ = 0e+0_fp

    ALLOCATE( SO4_ANTH( State_Grid%NX, State_Grid%NY, State_Grid%NZ ), &
              STAT=RC )
    CALL GC_CheckVar( 'sulfate_mod.F90:SO4_ANTH', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SO4_ANTH = 0e+0_fp
#endif

    !=================================================================
    ! Only initialize the following for offline aerosol simulations
    !=================================================================
    IF ( Input_Opt%ITS_AN_AEROSOL_SIM ) THEN

       ALLOCATE( TCOSZ( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'sulfate_mod.F90:TCOSZ', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       TCOSZ = 0e+0_fp

       ALLOCATE( TTDAY( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'sulfate_mod.F90:TTDAY', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       TTDAY = 0e+0_fp

       ALLOCATE( COSZM( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'sulfate_mod.F90:COSZM', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       COSZM = 0e+0_fp

       ALLOCATE(GLOBAL_OH(State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=RC)
       CALL GC_CheckVar( 'sulfate_mod.F90:GLOBAL_OH', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GLOBAL_OH = 0e+0_fp

       ALLOCATE(GLOBAL_HNO3(State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=RC)
       CALL GC_CheckVar( 'sulfate_mod.F90:GLOBAL_HNO3', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GLOBAL_HNO3 = 0e+0_fp

       ALLOCATE(GLOBAL_HCl(State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=RC)
       CALL GC_CheckVar( 'sulfate_mod.F90:GLOBAL_HCl', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GLOBAL_HCl = 0e+0_fp

       ALLOCATE(GLOBAL_HCOOH(State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=RC)
       CALL GC_CheckVar( 'sulfate_mod.F90:GLOBAL_HCOOH', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GLOBAL_HCOOH = 0e+0_fp

       ALLOCATE(GLOBAL_ACTA(State_Grid%NX,State_Grid%NY,State_Grid%NZ), STAT=RC)
       CALL GC_CheckVar( 'sulfate_mod.F90:GLOBAL_ACTA', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       GLOBAL_ACTA = 0e+0_fp
    ENDIF

    !================================================================
    ! Find drydep species
    !=================================================================

    ! Define flags for species ID's
    id_AS    = Ind_('AS'       )
    id_AHS   = Ind_('AHS'      )
    id_AW1   = Ind_('AW1'      )
    id_DAL1  = Ind_('DSTAL1'   )
    id_DAL2  = Ind_('DSTAL2'   )
    id_DAL3  = Ind_('DSTAL3'   )
    id_DAL4  = Ind_('DSTAL4'   )
    id_DMS   = Ind_('DMS'      )
    id_DST1  = Ind_('DST1'     )
    id_DST2  = Ind_('DST2'     )
    id_DST3  = Ind_('DST3'     )
    id_DST4  = Ind_('DST4'     )
    id_H2O2  = Ind_('H2O2'     )
    id_HNO3  = Ind_('HNO3'     )
    id_HOBr  = Ind_('HOBr'     )
    id_HOCl  = Ind_('HOCl'     )
    id_LET   = Ind_('LET'      )
    id_MSA   = Ind_('MSA'      )
    id_NH3   = Ind_('NH3'      )
    id_NH4   = Ind_('NH4'      )
    id_NH4aq = Ind_('NH4aq'    )
    id_NIT   = Ind_('NIT'      )
    id_NITd1 = Ind_('NITD1'    )
    id_NITd2 = Ind_('NITD2'    )
    id_NITd3 = Ind_('NITD3'    )
    id_NITd4 = Ind_('NITD4'    )
    id_NITs  = Ind_('NITs'     )
    id_NK1   = Ind_('NK1'      )
    id_NK5   = Ind_('NK5'      )
    id_NK8   = Ind_('NK8'      )
    id_NK10  = Ind_('NK10'     )
    id_NK20  = Ind_('NK20'     )
    id_NO3   = Ind_('NO3'      )
    id_O3    = Ind_('O3'       )
    id_OH    = Ind_('OH'       )
    id_pFe   = Ind_('pFe'      )
    id_PSO4  = Ind_('PSO4'     )
    id_SALA  = Ind_('SALA'     )
    id_SALC  = Ind_('SALC'     )
    id_SF1   = Ind_('SF1'      )
    id_SO2   = Ind_('SO2'      )
    id_SO4   = Ind_('SO4'      )
    id_SO4aq = Ind_('SO4aq'    )
    id_SO4d1 = Ind_('SO4D1'    )
    id_SO4d2 = Ind_('SO4D2'    )
    id_SO4d3 = Ind_('SO4D3'    )
    id_SO4d4 = Ind_('SO4D4'    )
    id_SO4s  = Ind_('SO4s'     )
    id_SALACL= Ind_('SALACL'   )
    id_HCL   = Ind_('HCL'     )
    !id_NH4s  = Ind_('NH4s'     )
    id_SALCCL= Ind_('SALCCL'   )
    id_SALAAL= Ind_('SALAAL'   )
    id_SALCAL= Ind_('SALCAL'   )
    id_SO4H1 = Ind_('SO4H1'    )
    id_SO4H2 = Ind_('SO4H2'    )
    id_SO4H3 = Ind_('SO4H3'    )
    id_SO4H4 = Ind_('SO4H4'    )
    id_HCOOH = Ind_('HCOOH'    ) ! (jmm, 12/3/18)
    id_ACTA  = Ind_('ACTA'     ) ! (jmm, 12/3/18)
    id_CH2O  = Ind_('CH2O'     ) ! (jmm, 06/15/18)
    id_HMS   = Ind_('HMS'      ) ! (jmm, 06/15/18)

    ! Define flags for species drydep indices
    DRYNITs  = Ind_('NITs', 'D')
    DRYNITd1 = Ind_('NITD1','D')
    DRYNITd2 = Ind_('NITD2','D')
    DRYNITd3 = Ind_('NITD3','D')
    DRYNITd4 = Ind_('NITD4','D')
    DRYSO4s  = Ind_('SO4s', 'D')
    DRYSO4d1 = Ind_('SO4d1','D')
    DRYSO4d2 = Ind_('SO4d2','D')
    DRYSO4d3 = Ind_('SO4d3','D')
    DRYSO4d4 = Ind_('SO4d4','D')

    ! Error check the dust uptake species
    IF ( Input_Opt%LDSTUP ) THEN
       IF ( id_DAL1  < 0 .or. id_DAL2  < 0 .or. id_DAL3  < 0   .or. &
            id_DAL4  < 0 .or. id_NITd1 < 0 .or. id_NITd2 < 0   .or. &
            id_NITd3 < 0 .or. id_NITd4 < 0 .or. id_SO4d1 < 0   .or. &
            id_SO4d2 < 0 .or. id_SO4d3 < 0 .or. id_SO4d4 < 0 ) THEN
          ErrMsg =' Dust uptake tracers are undefined!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE INIT_SULFATE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_sulfate
!
! !DESCRIPTION: Subroutine CLEANUP\_SULFATE deallocates all previously
!  allocated arrays for sulfate emissions -- call at the end of the run
!  (bmy, 6/1/00, 10/15/09)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_SULFATE()
!
! !REVISION HISTORY:
!  01 Jun 2000 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    IF ( ALLOCATED( DMSo        ) ) DEALLOCATE( DMSo        )
    IF ( ALLOCATED( PMSA_DMS    ) ) DEALLOCATE( PMSA_DMS    )
    IF ( ALLOCATED( PNITs       ) ) DEALLOCATE( PNITs       )
    IF ( ALLOCATED( PNIT        ) ) DEALLOCATE( PNIT        )
    IF ( ALLOCATED( PACL        ) ) DEALLOCATE( PACL        )
    IF ( ALLOCATED( PCCL        ) ) DEALLOCATE( PCCL        )
    IF ( ALLOCATED( PSO2_DMS    ) ) DEALLOCATE( PSO2_DMS    )
    IF ( ALLOCATED( PSO4_SO2    ) ) DEALLOCATE( PSO4_SO2    )
    IF ( ALLOCATED( PHMS_SO2    ) ) DEALLOCATE( PHMS_SO2    )
    IF ( ALLOCATED( PSO2_HMS    ) ) DEALLOCATE( PSO2_HMS    )
    IF ( ALLOCATED( PSO4_HMS    ) ) DEALLOCATE( PSO4_HMS    )
#ifdef APM
    IF ( ALLOCATED( PSO4_SO2APM ) ) DEALLOCATE( PSO4_SO2APM )
    IF ( ALLOCATED( PSO4_SO2SEA ) ) DEALLOCATE( PSO4_SO2SEA )
#endif
#ifdef TOMAS
    IF ( ALLOCATED( PSO4_SO2AQ  ) ) DEALLOCATE( PSO4_SO2AQ  )
    IF ( ALLOCATED( SO4_ANTH    ))  DEALLOCATE( SO4_ANTH    )
#endif
    IF ( ALLOCATED( PSO4_ss     ) ) DEALLOCATE( PSO4_ss     )
    IF ( ALLOCATED( PSO4_dust   ) ) DEALLOCATE( PSO4_dust   )
    IF ( ALLOCATED( PNIT_dust   ) ) DEALLOCATE( PNIT_dust   )
    IF ( ALLOCATED( SOx_SCALE   ) ) DEALLOCATE( SOx_SCALE   )
    IF ( ALLOCATED( SSTEMP      ) ) DEALLOCATE( SSTEMP      )
    IF ( ALLOCATED( TCOSZ       ) ) DEALLOCATE( TCOSZ       )
    IF ( ALLOCATED( TTDAY       ) ) DEALLOCATE( TTDAY       )
    IF ( ALLOCATED( COSZM       ) ) DEALLOCATE( COSZM       )
    IF ( ALLOCATED( GLOBAL_OH   ) ) DEALLOCATE( GLOBAL_OH   )
    IF ( ALLOCATED( GLOBAL_HNO3 ) ) DEALLOCATE( GLOBAL_HNO3 )
    IF ( ALLOCATED( GLOBAL_HCl  ) ) DEALLOCATE( GLOBAL_HCl  )
    IF ( ALLOCATED( GLOBAL_HCOOH) ) DEALLOCATE( GLOBAL_HCOOH)
    IF ( ALLOCATED( GLOBAL_ACTA ) ) DEALLOCATE( GLOBAL_ACTA )

    ! Free pointers
    IF ( ASSOCIATED( NDENS_SALA ) ) NDENS_SALA => NULL()
    IF ( ASSOCIATED( NDENS_SALC ) ) NDENS_SALC => NULL()

  END SUBROUTINE CLEANUP_SULFATE
!EOC
#ifdef APM
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dry_settlingbin
!
! !DESCRIPTION: Subroutine WET\_SETTLINGBIN computes the dry settling of
!  aerosol tracers. Modified for APM simulation. (G. Luo)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE WET_SETTLINGBIN( Input_Opt, State_Chm, State_Diag, State_Grid, &
                              State_Met, RC  )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE Species_Mod,    ONLY : SpcConc
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Met_Mod,  ONLY : MetState
    USE PRESSURE_MOD,   ONLY : GET_PCENTER
    USE TIME_MOD,       ONLY : GET_TS_CHEM
    USE PhysConstants
    USE APM_INIT_MOD,   ONLY : APMIDS
    USE APM_INIT_MOD,   ONLY : NCTSO4,NSO4
    USE APM_INIT_MOD,   ONLY : RDRY
    USE APM_DRIV_MOD,   ONLY : GFTOT3D,DENWET3D
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Local variables
    INTEGER               :: I, J, L, N, K
    INTEGER               :: IDTEMP1, IDTEMP2
    REAL*8                :: DT_SETTL, DELZ,  DELZ1
    REAL*8                :: REFF,     DEN,   CONST
    REAL*8                :: NUM,      LAMDA, FLUX
    REAL*8                :: AREA_CM2, TC0(State_Grid%NZ)
    REAL*8                :: TOT1,     TOT2

    ! Pressure in Kpa 1 mb = 100 pa = 0.1 kPa
    REAL*8                :: P

    ! Diameter of aerosol [um]
    REAL*8                :: Dp

    ! Pressure * DP
    REAL*8                :: PDp

    ! Temperature (K)
    REAL*8                :: TEMP

    ! Slip correction factor
    REAL*8                :: Slip

    ! Viscosity of air (Pa s)
    REAL*8                :: Visc

    ! Settling velocity of particle (m/s)
    REAL*8                :: VTS(State_Grid%NZ)
    REAL*8                :: MASS(State_Grid%NZ)
    REAL*8                :: OLD(State_Grid%NZ,NCTSO4+2)

    ! Make a pointer to the tracer array
    TYPE(SpcConc), POINTER :: Spc(:)

    !=================================================================
    ! WET_SETTLINGBIN begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS

    ! Point to Spc
    Spc => State_Chm%Species

    ! Aerosol settling timestep [s]
    DT_SETTL = GET_TS_CHEM()

    IDTEMP1 = APMIDS%id_SO4BIN1
    IDTEMP2 = APMIDS%id_SO4BIN1+NSO4-1

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, N, K, DEN, REFF, DP )       &
    !$OMP PRIVATE( CONST, VTS, TEMP, P, PDP, SLIP )     &
    !$OMP PRIVATE( MASS, OLD, VISC, TC0, DELZ, DELZ1  ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       DO L = 1, State_Grid%NZ
          MASS(L) = 0.d8
          DO N = IDTEMP1, IDTEMP2
             MASS(L) = MASS(L) + Spc(N)%Conc(I,J,L)
          ENDDO
          DO K = 1, NCTSO4
             OLD(L,K) = Spc(APMIDS%id_CTSO4+K-1)%Conc(I,J,L)
             Spc(APMIDS%id_CTSO4+K-1)%Conc(I,J,L) = 0.D0
          ENDDO
          OLD(L,NCTSO4+1) = Spc(APMIDS%id_NIT)%Conc(I,J,L)
          OLD(L,NCTSO4+2) = Spc(APMIDS%id_NH4)%Conc(I,J,L)
          Spc(APMIDS%id_NIT)%Conc(I,J,L) = 0.D0
          Spc(APMIDS%id_NH4)%Conc(I,J,L) = 0.D0
       ENDDO

       ! Loop over aerosol bins
       DO N = 1, NSO4

          DO L = 1, State_Grid%NZ

             TC0(L) = Spc(APMIDS%id_SO4BIN1+N-1)%Conc(I,J,L)

             IF(TC0(L)>1.D-30)THEN
                ! Initialize
                DEN   = DENWET3D(I,J,L,1)*1.d3
                REFF  = RDRY(N)*GFTOT3D(I,J,L,1)
                DP    = 2D0 * REFF * 1.D6 ! Dp [um] = particle diameter
                CONST = 2D0 * DEN * REFF**2 * g0 / 9D0

                ! Get P [kPa], T [K], and P*DP
                P    = GET_PCENTER(I,J,L) * 0.1d0
                TEMP = State_Met%T(I,J,L)
                PDP  = P * DP

                ! Slip correction factor as function of (P*dp)
                SLIP = 1d0 + ( 15.60d0 + 7.0d0 * EXP(-0.059d0*PDP) ) / PDP

                ! Viscosity [Pa s] of air as a function of temp (K)
                VISC = 1.458d-6 * (TEMP)**(1.5d0) / ( TEMP + 110.4d0 )

                ! Settling velocity [m/s]
                VTS(L) = CONST * SLIP / VISC
             ELSE
                VTS(L) = 0.D0
             ENDIF

          ENDDO

          ! Method is to solve bidiagonal matrix
          ! which is implicit and first order accurate in Z
          L    = State_Grid%NZ
          IF(MASS(L)>1.D-30)THEN
             DELZ = State_Met%BXHEIGHT(I,J,L)
             Spc(APMIDS%id_SO4BIN1+N-1)%Conc(I,J,L) =                &
                  Spc(APMIDS%id_SO4BIN1+N-1)%Conc(I,J,L) /           &
                  ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )

             DO K = 1, NCTSO4
                Spc(APMIDS%id_CTSO4+K-1)%Conc(I,J,L) =               &
                     Spc(APMIDS%id_CTSO4+K-1)%Conc(I,J,L)+           &
                     OLD(L,K)*TC0(L)/MASS(L) /                   &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )
             ENDDO
             Spc(APMIDS%id_NIT)%Conc(I,J,L) =                          &
                  Spc(APMIDS%id_NIT)%Conc(I,J,L)+                      &
                  OLD(L,NCTSO4+1)*TC0(L)/MASS(L) /               &
                  ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )
             Spc(APMIDS%id_NH4)%Conc(I,J,L) =                          &
                  Spc(APMIDS%id_NH4)%Conc(I,J,L)+                      &
                  OLD(L,NCTSO4+2)*TC0(L)/MASS(L) /               &
                  ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )
          ENDIF

          DO L = State_Grid%NZ-1, 1, -1
             IF((MASS(L)*MASS(L+1))>1.D-30)THEN
                DELZ  = State_Met%BXHEIGHT(I,J,L)
                DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
                Spc(APMIDS%id_SO4BIN1+N-1)%Conc(I,J,L) = 1.e+0_fp /  &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )     &
                     * (Spc(APMIDS%id_SO4BIN1+N-1)%Conc(I,J,L)       &
                     + DT_SETTL * VTS(L+1) / DELZ1  * TC0(L+1) )

                DO K = 1, NCTSO4
                   Spc(APMIDS%id_CTSO4+K-1)%Conc(I,J,L) =            &
                        Spc(APMIDS%id_CTSO4+K-1)%Conc(I,J,L)+        &
                        1.e+0_fp /                               &
                        ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )  &
                        * (OLD(L,K)*TC0(L)/MASS(L)               &
                        + DT_SETTL * VTS(L+1) / DELZ1            &
                        * OLD(L+1,K)*TC0(L+1)/MASS(L+1) )
                ENDDO
                Spc(APMIDS%id_NIT)%Conc(I,J,L) =                       &
                     Spc(APMIDS%id_NIT)%Conc(I,J,L)+                   &
                     1.e+0_fp /                                  &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )     &
                     * (OLD(L,NCTSO4+1)*TC0(L)/MASS(L)           &
                     + DT_SETTL * VTS(L+1) / DELZ1               &
                     * OLD(L+1,NCTSO4+1)*TC0(L+1)/MASS(L+1) )
                Spc(APMIDS%id_NH4)%Conc(I,J,L) =                       &
                     Spc(APMIDS%id_NH4)%Conc(I,J,L)+                   &
                     1.e+0_fp /                                  &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )     &
                     * (OLD(L,NCTSO4+2)*TC0(L)/MASS(L)           &
                     + DT_SETTL * VTS(L+1) / DELZ1               &
                     * OLD(L+1,NCTSO4+2)*TC0(L+1)/MASS(L+1) )

             ELSE IF(MASS(L)>1.D-30)THEN
                DELZ  = State_Met%BXHEIGHT(I,J,L)
                DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
                Spc(APMIDS%id_SO4BIN1+N-1)%Conc(I,J,L) = 1.e+0_fp /  &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )     &
                     * (Spc(APMIDS%id_SO4BIN1+N-1)%Conc(I,J,L)       &
                     + DT_SETTL * VTS(L+1) / DELZ1  * TC0(L+1) )

                DO K = 1, NCTSO4
                   Spc(APMIDS%id_CTSO4+K-1)%Conc(I,J,L) =            &
                        Spc(APMIDS%id_CTSO4+K-1)%Conc(I,J,L)+        &
                        1.e+0_fp /                               &
                        ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )  &
                        * OLD(L,K)*TC0(L)/MASS(L)
                ENDDO
                Spc(APMIDS%id_NIT)%Conc(I,J,L) =                       &
                     Spc(APMIDS%id_NIT)%Conc(I,J,L)+                   &
                     1.e+0_fp /                                  &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )     &
                     * OLD(L,NCTSO4+1)*TC0(L)/MASS(L)
                Spc(APMIDS%id_NH4)%Conc(I,J,L) =                       &
                     Spc(APMIDS%id_NH4)%Conc(I,J,L)+                   &
                     1.e+0_fp /                                  &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )     &
                     * OLD(L,NCTSO4+2)*TC0(L)/MASS(L)

             ELSE IF(MASS(L+1)>1.D-30)THEN
                DELZ  = State_Met%BXHEIGHT(I,J,L)
                DELZ1 = State_Met%BXHEIGHT(I,J,L+1)
                Spc(APMIDS%id_SO4BIN1+N-1)%Conc(I,J,L) = 1.e+0_fp /  &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )     &
                     * (Spc(APMIDS%id_SO4BIN1+N-1)%Conc(I,J,L)       &
                     + DT_SETTL * VTS(L+1) / DELZ1 * TC0(L+1) )

                DO K = 1, NCTSO4
                   Spc(APMIDS%id_CTSO4+K-1)%Conc(I,J,L) =            &
                        Spc(APMIDS%id_CTSO4+K-1)%Conc(I,J,L)+        &
                        1.e+0_fp /                               &
                        ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )  &
                        * DT_SETTL * VTS(L+1) / DELZ1            &
                        * OLD(L+1,K)*TC0(L+1)/MASS(L+1)
                ENDDO
                Spc(APMIDS%id_NIT)%Conc(I,J,L) =                       &
                     Spc(APMIDS%id_NIT)%Conc(I,J,L)+                   &
                     1.e+0_fp /                                  &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )     &
                     * DT_SETTL * VTS(L+1) / DELZ1               &
                     * OLD(L+1,NCTSO4+1)*TC0(L+1)/MASS(L+1)
                Spc(APMIDS%id_NH4)%Conc(I,J,L) =                       &
                     Spc(APMIDS%id_NH4)%Conc(I,J,L)+                   &
                     1.e+0_fp /                                  &
                     ( 1.e+0_fp + DT_SETTL * VTS(L) / DELZ )     &
                     * DT_SETTL * VTS(L+1) / DELZ1               &
                     * OLD(L+1,NCTSO4+2)*TC0(L+1)/MASS(L+1)

             ENDIF

          ENDDO

       ENDDO

       DO L = 1, State_Grid%NZ
          DO K = 1, NCTSO4
             Spc(APMIDS%id_CTSO4+K-1)%Conc(I,J,L) = &
                  MAX(1.d-30,Spc(APMIDS%id_CTSO4+K-1)%Conc(I,J,L))
          ENDDO
          Spc(APMIDS%id_NIT)%Conc(I,J,L) = MAX(1.d-30,Spc(APMIDS%id_NIT)%Conc(I,J,L))
          Spc(APMIDS%id_NH4)%Conc(I,J,L) = MAX(1.d-30,Spc(APMIDS%id_NH4)%Conc(I,J,L))
       ENDDO

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Clear the pointer
    NULLIFY( Spc )

  END SUBROUTINE WET_SETTLINGBIN
#endif
END MODULE SULFATE_MOD
