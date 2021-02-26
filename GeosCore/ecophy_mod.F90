!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ecophy_mod.F90
!
! !DESCRIPTION: Module ECOPHY\_MOD contains variables and routines for the
!  GEOS-Chem ecophysiology scheme. It modifies the bulk canopy stomatal
!  resistance RIX in the dry deposition scheme.
!\\
!\\
! !INTERFACE:
!
      MODULE ECOPHY_MOD
!
! !USES:
!
      USE PhysConstants, ONLY: RSTARG, AIRMW    ! Physical constants
      USE PRECISION_MOD                         ! For GEOS-Chem Precision (fp)
      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: DO_ECOPHY
      PUBLIC :: INIT_ECOPHY
      PUBLIC :: CLEANUP_ECOPHY
!
! !REVISION HISTORY:
!
! !REMARKS:
!  References:
!  ============================================================================
!  (1 ) Clark, D. B., et al., "The Joint UK Land Environment Simulator (JULES), 
!        model description – Part 2: Carbon fluxes and vegetation dynamics", 
!        Geosci. Model Dev., 4, 701–722, 
!        https://doi.org/10.5194/gmd-4-701-2011, 2011.
!  (2 ) Best, M. J., et al. "The Joint UK Land Environment Simulator (JULES), 
!        model description – Part 1: Energy and water fluxes", 
!        Geosci. Model Dev., 4, 677–699, 
!        https://doi.org/10.5194/gmd-4-677-2011, 2011.
!  (3 ) P.M Cox, C Huntingford, R.J Harding, "A canopy conductance and 
!        photosynthesis model for use in a GCM land surface scheme",
!        Journal of Hydrology, Volumes 212–213, Pages 79-94, ISSN 0022-1694,
!        https://doi.org/10.1016/S0022-1694(98)00203-0, 1998.
!  (4 ) Raoult, N. M., "Land-surface parameter optimisation using data 
!        assimilation techniques: the adJULES system V1.0",
!        Geosci. Model Dev., 9, 2833–2852, 
!        https://doi.org/10.5194/gmd-9-2833-2016, 2016.
!  (5 ) Pacifico, F., et al., "Photosynthesis-based biogenic isoprene emission 
!        scheme in JULES",
!        Atmos. Chem. Phys., 11, 4371–4389, 2011
!        https://doi.org/10.5194/acp-11-4371-2011, 2011
!  ============================================================================
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE DATA MEMBERS:
!
      !---------------------------------------------------------------------------------------
      ! NUMPFT        : Total number of PFTs                              []
      ! MAX_ITER      : Maximum number of iterations                      []
      ! IS_C3_PLANT   : IS_C3_PLANT = 1 for C3 plants, else 0             []
      ! ALPHA         : Quantum efficiency of photosynthesis              [mol CO2 mol^-1 PAR]
      ! V_CMAX25      : V_CMAX at 25 deg C                                [mol CO2 m^-2 s^-1]
      ! T_UPP         : PFT-specific parameter for V_CMAX                 [deg C]
      ! T_LOW         : PFT-specific parameter for V_CMAX                 [deg C]
      ! F_DARKRESP    : Dark respiration coefficient                      []
      ! D_STAR        : PFT-specific parameter for closure eq.            [kg H2O / kg air]
      ! f0            : PFT-specific parameter for closure eq.            []
      ! G_LEAF_MIN    : PFT-specific min leaf conductance for closure eq. [m s^-1]
      ! K_EXTINCT     : Light extinction coefficient                      []
      ! PARAM_A_LOW   : PFT-specific parameter for low sensitivity   
      !                 ozone damage scheme                               [m^2 s nmol^-1]
      ! PARAM_A_HIGH  : PFT-specific parameter for high sensitivity 
      !                 ozone damage scheme                               [m^2 s nmol^-1]
      ! FLUXO3_CRIT   : PFT-specific threshold for ozone uptake           [nmol m^-2 s^-1]
      ! THRESHOLD     : Threshold of relative error                       []
      ! CO2_O2_RATIO  : Ratio of diffusivity of CO2 compared to H2O       []
      ! Leaf_density  : Specific leaf density                             [kg C m^-2 leaf]
      ! IEF           : Isoprene Emission Factors under standard cond.    [ug C g^-1 (dry weight) hr^-1]
      ! TEMP_st       : Standard conditions (temperature)                 [K]
      ! PAR_st        : Standard conditions (PAR)                         []
      !---------------------------------------------------------------------------------------
      INTEGER,  PARAMETER   :: NUMPFT                 = 5 
      INTEGER,  PARAMETER   :: MAX_ITER               = 1
      INTEGER,  PARAMETER   :: IS_C3_PLANT   (NUMPFT) = (/ 1,1,1,0,1 /)
      ! REAL(fp), PARAMETER   :: ALPHA         (NUMPFT) = (/ 0.08, 0.08, 0.12, 0.06, 0.08 /)
      ! REAL(fp), PARAMETER   :: V_CMAX25      (NUMPFT) = (/ 0.046, 0.033, 0.073, 0.060, 0.060 /) * &
      !                                                   (/ 0.0008, 0.0008, 0.0008, 0.0004, 0.0008 /)
      ! REAL(fp), PARAMETER   :: T_UPP         (NUMPFT) = (/ 36.0, 26.0, 36.0, 45.0, 36.0 /)
      ! REAL(fp), PARAMETER   :: T_LOW         (NUMPFT) = (/ 0.0, -10.0, 0.0,  13.0, 0.0  /)
      REAL(fp), PARAMETER   :: F_DARKRESP    (NUMPFT) = (/ .015, .015, .015, .025, .015 /)
      ! REAL(fp), PARAMETER   :: D_STAR        (NUMPFT) = (/ 0.09, 0.06, 0.1,  .075, 0.1  /)
      ! REAL(fp), PARAMETER   :: f0            (NUMPFT) = (/ .875, .875, 0.9,  0.8,  0.9  /)
      REAL(fp), PARAMETER   :: G_LEAF_MIN    (NUMPFT) = 1.0e-4                                ! gives RS ~ 9999
      REAL(fp), PARAMETER   :: K_EXTINCT     (NUMPFT) = 0.5
      REAL(fp), PARAMETER   :: PARAM_A_LOW   (NUMPFT) = (/ 0.04, 0.02, 0.25, 0.13, 0.03 /)    ! low sensitivity
      REAL(fp), PARAMETER   :: PARAM_A_HIGH  (NUMPFT) = (/ 0.15, 0.075, 1.40, 0.735, 0.10 /)  ! high sensistivity
      REAL(fp), PARAMETER   :: FLUXO3_CRIT   (NUMPFT) = (/ 1.6,  1.6,  5.0,  5.0,  1.6  /)
      REAL(fp), PARAMETER   :: THRESHOLD              = 1.0e-3
      ! Second set of optimized parameters (from Raoult et al. 2016)
      REAL(fp), PARAMETER   :: ALPHA         (NUMPFT) = (/ 0.131, 0.096, 0.179, 0.118, 0.102 /)
      REAL(fp), PARAMETER   :: V_CMAX25      (NUMPFT) = (/ 0.061, 0.065, 0.070, 0.051, 0.041 /) * &
                                                        (/ 0.0008, 0.0008, 0.0008, 0.0004, 0.0008 /)
      REAL(fp), PARAMETER   :: T_UPP         (NUMPFT) = (/ 38.578, 34.721, 36.242, 44.897, 35.385 /)
      REAL(fp), PARAMETER   :: T_LOW         (NUMPFT) = (/ 1.203, -8.698, -1.985, 11.37, -5.208 /)
      REAL(fp), PARAMETER   :: D_STAR        (NUMPFT) = (/ 0.048, 0.036, 0.086, 0.046, 0.077 /)
      REAL(fp), PARAMETER   :: f0            (NUMPFT) = (/ 0.765, 0.737, 0.817, 0.765, 0.782 /)  
      REAL(fp), PARAMETER   :: CO2_O2_RATIO = 1.6

      ! Parameters for photosynthesis-dependent isoprene emission
      REAL(fp), PARAMETER   :: Leaf_density  (NUMPFT) = (/ 0.0375, 0.1000, 0.0250, 0.0500, 0.0500 /) 
      REAL(fp), PARAMETER   :: IEF           (NUMPFT) = (/ 35.0, 12.0, 16.0, 8.0,  20.0 /)
      REAL(fp), PARAMETER   :: Dry_fraction           = 0.4
      REAL(fp), PARAMETER   :: TEMPK_st               = 303.15e+0_fp

!
! MODULE VARIABLES:
!
      !=========================================================================================
      ! SATU          : Soil moisture at saturation point                 [m^3 water / m^3 soil]
      ! CRIT          : Soil moisture at critical point                   [m^3 water / m^3 soil]
      ! WILT          : Soil moisture at wilting point                    [m^3 water / m^3 soil]
      ! O3dmg_opt     : Control switch for ozone damage scheme            [high/low/off]
      ! A_NET_st      : Leaf photosynthesis rate under standard cond.     [mol CO2 m^-2 leaf s^-1]
      ! RESP_st       : Leaf respiration rate under standard cond.        [mol CO2 m^-2 leaf s^-1]
      ! CO2_IN_st     : Leaf CO2 partial pressure under standard cond.    [Pa]
      !=========================================================================================
      REAL(fp)              :: SATU
      REAL(fp)              :: CRIT
      REAL(fp)              :: WILT
      CHARACTER(len=4)      :: O3dmg_opt 
      REAL(fp)              :: A_NET_st      (NUMPFT)
      REAL(fp)              :: RESP_st       (NUMPFT)
      REAL(fp)              :: CO2_IN_st     (NUMPFT)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_ecophy
!
! !DESCRIPTION: Subroutine DO\_ECOPHY is the interface between dry
!  deposition module and the ecophysiology module. It computes the
!  bulk canopy stomatal resistance r_s according to meterological inputs,
!  plant functional types, and soil types.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_ECOPHY ( Input_Opt, State_Met,             &
                             State_Chm, State_Diag, RC,        &
                             I, J,      LDT, PFT,   RA, RB_O3, &
                             PRESSU,    Isop_from_Ecophy,      &
                             RS,        SumLAI_PFT, IUSE_PFT   )
!
! !USES:
!
      USE ErrCode_Mod
      USE Input_Opt_Mod,      ONLY : OptInput
      USE Species_Mod,        ONLY : Species
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      USE State_Diag_Mod,     ONLY : DgnState
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)    :: I           ! longitude index
      INTEGER,        INTENT(IN)    :: J           ! latitude index
      INTEGER,        INTENT(IN)    :: PFT         ! PFT index
      INTEGER,        INTENT(IN)    :: LDT         ! Land type index (for archiving)
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
      REAL(fp),       INTENT(IN)    :: RA          ! Aerodynamic resistance 
      REAL(fp),       INTENT(IN)    :: RB_O3       ! Boundary layer resistance
      REAL(fp),       INTENT(IN)    :: PRESSU      ! Surface Pressure (Pa)
      REAL(fp),       INTENT(IN)    :: SumLAI_PFT  ! leaf area of the PFT
      INTEGER,        INTENT(IN)    :: IUSE_PFT    ! fraction of grid box 
                                                   ! occupied by the PFT
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
      REAL(fp),       INTENT(OUT)   :: RS          ! Bulk canopy stomatal resistance
      REAL(fp),       INTENT(OUT)   :: Isop_from_Ecophy  ! Isoprene emission
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      !  Note: Read subroutine DO_PHOTOSYNTHESIS for descriptions of variables.
      REAL(fp)   :: TEMPK
      REAL(fp)   :: QV2M
      REAL(fp)   :: APAR
      ! REAL(fp)   :: PRESSURE
      REAL(fp)   :: CO2
      REAL(fp)   :: O2
      REAL(fp)   :: O3
      REAL(fp)   :: SOIL_WETNESS
      REAL(fp)   :: G_CAN_OUT
      REAL(fp)   :: G_LEAF_OUT
      REAL(fp)   :: CO2_IN
      REAL(fp)   :: A_CAN_OUT
      REAL(fp)   :: A_NET_OUT
      REAL(fp)   :: RESP_CAN_OUT
      REAL(fp)   :: RESP_OUT
      REAL(fp)   :: FLUXO3_CAN
      REAL(fp)   :: FLUXO3
      REAL(fp)   :: FACTOR_O3
      REAL(fp)   :: BETA
      REAL(fp)   :: V_CMAX
      REAL(fp)   :: RATE_LIGHT
      REAL(fp)   :: RATE_RUBISCO
      REAL(fp)   :: RATE_PRODUCT
      REAL(fp)   :: A_GROSS
      REAL(fp)   :: LAI
      REAL(fp)   :: ISOP_EMIS
      INTEGER    :: IOLSON
      INTEGER    :: IUSE

      !=================================================================
      ! DO_ECOPHY begins here!
      !=================================================================
      ! Assume success
      RC = GC_SUCCESS

      ! Initialize
      ErrMsg  = ''
      ThisLoc = &
      ' -> at Do_ECOPHY (in module GeosCore/ecophysiology.F90)'

      ! get inputs for the module
      CALL GET_ECOPHY_INPUTS( State_Met,    State_Chm, Input_Opt,&
                              I, J, LDT,                         &
                              TEMPK,        QV2M,                &
                              APAR,         CO2,                 &
                              O2,           LAI,       O3,       &
                              SOIL_WETNESS, IUSE                 &
                              )

      ! simulate plant processes
      CALL DO_PHOTOSYNTHESIS( LAI, PFT,     PRESSU,      APAR,         &
                              CO2, O2,      TEMPK,       G_CAN_OUT,    &
                              G_LEAF_OUT,   CO2_IN,      A_CAN_OUT,    &
                              A_NET_OUT,    RESP_CAN_OUT, RESP_OUT,    &
                              FLUXO3_CAN,   FLUXO3,      FACTOR_O3,    &
                              BETA,         V_CMAX,      RATE_LIGHT,   &
                              RATE_RUBISCO, RATE_PRODUCT, A_GROSS, RC, &
                              QV2M=QV2M,    RA=RA,       RB_O3=RB_O3,  &
                              O3=O3,        SOIL_WETNESS=SOIL_WETNESS  &
                              )

      IF ( Input_Opt%LIsop_from_Ecophy ) THEN
         ! calculate isoprene emission from photosynthesis 
         CALL DO_ISOP_EMIS( A_CAN_OUT, RESP_CAN_OUT,  &
                            TEMPK, CO2_IN,            &
                            LAI, PFT,                 &
                            ISOP_EMIS                 )
      ENDIF

      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in call to "GET_ECOPHY_INPUTS!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Output bulk stomatal resistance to dry deposition module
      RS = 1.0 / G_CAN_OUT

      ! Output isoprene emission rate to dry deposition module
      ! Convert from kg C m-2 land s-1 to kg C m-2 grid s-1 
      ! IUSE / 1000 is the the fraction of the grid for that land type
      Isop_from_Ecophy = ISOP_EMIS * DBLE( IUSE ) / 1.e+3_fp

      ! Send output to State_Diag
      CALL Ecophy_Diagn( I, J,         LDT,        PFT,          &
                         IUSE,         LAI,        SumLAI_PFT,   &
                         IUSE_PFT,     RA,         RB_O3,        &
                         G_CAN_OUT,    A_CAN_OUT,  RESP_CAN_OUT, &
                         CO2_IN,       FLUXO3_CAN,               &
                         FACTOR_O3,    BETA,       V_CMAX,       &
                         RATE_LIGHT,   RATE_RUBISCO,             &
                         RATE_PRODUCT, A_GROSS,    ISOP_EMIS,    &
                         State_Met,    State_Diag, RC            &
                       )
      ! Trap potential errors
      IF ( RC /= GC_SUCCESS ) THEN
         ErrMsg = 'Error encountered in call to "GET_ECOPHY_INPUTS!'
         CALL GC_Error( ErrMsg, RC, ThisLoc )
         RETURN
      ENDIF

      ! Return to calling program
      END SUBROUTINE DO_ECOPHY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_photosynthesis
!
! !DESCRIPTION: Subroutine DO\_PHOTOSYNTHESIS is the main driver of this
!  module.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_PHOTOSYNTHESIS( LAI, PFT,     PRESSURE,    APAR,         &
                                    CO2, O2,      TEMPK,       G_CAN_OUT,    &
                                    G_LEAF_OUT,   CO2_IN,      A_CAN_OUT,    &
                                    A_NET_OUT,    RESP_CAN_OUT, RESP_OUT,    &
                                    FLUXO3_CAN,   FLUXO3,      FACTOR_O3,    &
                                    BETA,         V_CMAX,      RATE_LIGHT,   &
                                    RATE_RUBISCO, RATE_PRODUCT, A_GROSS, RC, &
                                    V_CMAX_IN,    QV2M,        DEFICIT_Q_IN, &
                                    RA, RB_O3,    O3,          FACTOR_O3_IN, &
                                    SOIL_WETNESS, BETA_IN                    &
                                    )
!
! !USES:
!
      USE ErrCode_Mod
!
!INPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! LAI           : Leaf area index for the PFT                       [m^2 m^-2]
      ! PFT           : Index for PFT                                     []
      ! APAR          : Absorbed PAR                                      [mol photon m^-2 s^-1]
      ! PRESSURE      : Atmospheric Pressure in canopy layer              [Pa]
      ! CO2           : Ambient CO2 mole fraction                         [mol/mol air]
      ! O2            : Ambient O2 mole fraction                          [mol/mol air]
      ! TEMPK         : Leaf temperature in Kelvin                        [K]
      ! QV2M          : Specific humidity in canopy layer                 [kg H2O / kg air]
      ! DEFICIT_Q_IN  : Alternative input of specific humidity deficit    [kg H2O / kg air]
      ! V_CMAX_IN     : Alternative input of V_CMAX                       [mol CO2 m^-2 s^-1]
      ! RA            : Aerodynamic resistance                            [s m^-1]
      ! RB_O3         : Boundary layer resistance                         [s m^-1]
      ! O3            : Ozone mole fraction in canopy layer               [mol/mol air]
      ! FACTOR_O3_IN  : Alternative input of ozone damage factor          []
      ! SOIL_WETNESS  : Fraction of moisture in soil pores                []
      ! BETA_IN       : Alternative input of soil moisture stress factor  []
      !---------------------------------------------------------------------------------------
      REAL(fp), INTENT(IN)             :: LAI
      INTEGER,  INTENT(IN)             :: PFT
      REAL(fp), INTENT(IN)             :: PRESSURE
      REAL(fp), INTENT(IN)             :: APAR
      REAL(fp), INTENT(IN)             :: CO2
      REAL(fp), INTENT(IN)             :: O2
      REAL(fp), INTENT(IN)             :: TEMPK

      ! If V_CMAX_IN is provided, it will be used and the module will not 
      ! calculate V_CMAX according to leaf temperature.
      REAL(fp), INTENT(IN), OPTIONAL   :: V_CMAX_IN

      ! Parameters for calculating specific humidity deficit
      ! If DEFICIT_Q_IN is provided, it will be used and the module will not 
      ! calculate DEFICIT_Q according to QV2M and leaf temperature.
      REAL(fp), INTENT(IN), OPTIONAL   :: QV2M
      REAL(fp), INTENT(IN), OPTIONAL   :: DEFICIT_Q_IN

      ! Parameters for ozone damage scheme
      ! If FACTOR_O3_IN is provided, it will be used as the ozone damage factor
      ! and the module will not calculate FACTOR_O3 according to 
      ! RA, RB and O3 from other modules.
      REAL(fp), INTENT(IN), OPTIONAL   :: RA
      REAL(fp), INTENT(IN), OPTIONAL   :: RB_O3
      REAL(fp), INTENT(IN), OPTIONAL   :: O3
      REAL(fp), INTENT(IN), OPTIONAL   :: FACTOR_O3_IN

      ! Parameters for soil moisture stress factor
      ! If BETA_IN is provided, it will be used as the soil moisture stress factor
      ! and the module will not calculate BETA according to SOIL_WETNESS from 
      ! input meteorology data.
      REAL(fp), INTENT(IN), OPTIONAL   :: SOIL_WETNESS
      REAL(fp), INTENT(IN), OPTIONAL   :: BETA_IN
!
!OUTPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! G_CAN_OUT     : Canopy conductance for H2O (output)               [m s^-1]
      ! G_LEAF_OUT    : Leaf level stomatal conductance for H2O (output)  [m s^-1]
      ! CO2_IN        : Leaf internal partial pressure of CO2             [Pa]
      ! A_CAN_OUT     : Canopy net photosynthetic rate (output)           [mol CO2 m^-2 s^-1]
      ! A_NET_OUT     : Leaf level net photosynthetic rate (output)       [mol CO2 m^-2 s^-1]
      ! RESP_CAN_OUT  : Canopy dark respiration (output)                  [mol CO2 m^-2 s^-1]
      ! RESP_OUT      : Leaf level dark respiration (output)              [mol CO2 m^-2 s^-1]
      ! FLUXO3_CAN    : Canopy ozone uptake                               [nmol m^-2 s^-1]
      ! FLUXO3        : Stomatal ozone uptake                             [nmol m^-2 s^-1]
      ! FACTOR_O3     : Ozone damage factor                               []
      ! BETA          : Soil moisture stress factor                       []
      ! V_CMAX        : Max Rubisco carboxylation rate                    [mol CO2 m^-2 s^-1]
      ! RATE_LIGHT    : Light-limited rate                                [mol CO2 m^-2 s^-1]
      ! RATE_PRODUCT  : Product-limited rate                              [mol CO2 m^-2 s^-1]
      ! RATE_RUBISCO  : Rubisco-limited rate                              [mol CO2 m^-2 s^-1]
      ! A_GROSS       : Leaf level gross photosynthesis                   [mol CO2 m^-2 s^-1]
      !---------------------------------------------------------------------------------------
      REAL(fp), INTENT(OUT) :: G_CAN_OUT
      REAL(fp), INTENT(OUT) :: G_LEAF_OUT
      REAL(fp), INTENT(OUT) :: CO2_IN
      REAL(fp), INTENT(OUT) :: A_CAN_OUT
      REAL(fp), INTENT(OUT) :: A_NET_OUT
      REAL(fp), INTENT(OUT) :: RESP_CAN_OUT
      REAL(fp), INTENT(OUT) :: RESP_OUT
      REAL(fp), INTENT(OUT) :: FLUXO3_CAN
      REAL(fp), INTENT(OUT) :: FLUXO3
      REAL(fp), INTENT(OUT) :: FACTOR_O3
      REAL(fp), INTENT(OUT) :: BETA
      REAL(fp), INTENT(OUT) :: V_CMAX
      REAL(fp), INTENT(OUT) :: RATE_LIGHT
      REAL(fp), INTENT(OUT) :: RATE_RUBISCO
      REAL(fp), INTENT(OUT) :: RATE_PRODUCT
      REAL(fp), INTENT(OUT) :: A_GROSS
      INTEGER,  INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
!LOCAL VARIABLES:
!
      !---------------------------------------------------------------------------------------
      ! TEMPC         : Leaf temperature in degree Celsius                [deg C]
      ! SPHU_SAT      : Saturation specific humidity in canopy layer      [kg H2O / kg air]
      ! DEFICIT_Q     : Specific humidity deficit at leaf surface         [kg H2O / kg air]
      ! CO2_AMBIENT   : Ambient CO2 partial pressure                      [Pa]
      ! O3_CONC       : Ozone concentration in canopy layer               [nmol m^-3]
      ! BIGLEAFSCALE  : Scaling with Big-leaf approach                    []
      ! G_CAN         : Canopy conductance for H2O                        [m s^-1]
      ! G_LEAF        : Leaf level stomatal conductance for H2O           [m s^-1]
      ! G_LEAF_PREV   : G_LEAF in previous iteration                      [m s^-1]
      ! CO2_IN_PREV   : CO2_IN in previous iteration                      [Pa]
      ! A_NET         : Leaf level net photosynthetic rate                [mol CO2 m^-2 s^-1]
      ! RESP          : Leaf level dark respiration                       [mol CO2 m^-2 s^-1]
      ! A_NET_PREV    : A_NET in previous iteration                       [mol CO2 m^-2 s^-1]
      ! V_CMAX        : Max Rubisco carboxylation rate                    [mol CO2 m^-2 s^-1]
      ! CO2_GAMMA     : CO2 Compensation point                            [Pa]
      ! RATE_LIGHT    : Light-limited rate                                [mol CO2 m^-2 s^-1]
      ! RATE_PRODUCT  : Product-limited rate                              [mol CO2 m^-2 s^-1]
      ! RATE_RUBISCO  : Rubisco-limited rate                              [mol CO2 m^-2 s^-1]
      ! A_GROSS       : Leaf level gross photosynthesis                   [mol CO2 m^-2 s^-1]
      ! TAU           : Rubisco specificity for CO2 to O2                 []
      ! DENOM         : Denominator for calculating V_CMAX                []
      ! ITER          : No. of iterations in calculating conductance      []
      ! ERR1          : Relative change for G_LEAF between iterations     []
      ! ERR2          : Relative change for CO2_IN between iterations     []
      ! ERR3          : Relative change for A_NET between iterations      []
      ! DELTA         : Maximum of the 3 relative changes                 []
      ! RAB           : Aerodynamic and boundary layer resistance         [s m^-1]
      !---------------------------------------------------------------------------------------
      REAL(fp)    :: TEMPC
      REAL(fp)    :: SPHU_SAT
      REAL(fp)    :: DEFICIT_Q
      REAL(fp)    :: CO2_AMBIENT
      REAL(fp)    :: O3_CONC
      REAL(fp)    :: BIGLEAFSCALE
      REAL(fp)    :: G_CAN
      REAL(fp)    :: G_LEAF
      ! REAL(fp)    :: G_LEAF_PREV
      ! REAL(fp)    :: CO2_IN_PREV
      REAL(fp)    :: A_NET
      REAL(fp)    :: RESP
      ! REAL(fp)    :: A_NET_PREV
      ! REAL(fp)    :: V_CMAX
      REAL(fp)    :: CO2_GAMMA
      ! REAL(fp)    :: RATE_LIGHT
      ! REAL(fp)    :: RATE_PRODUCT
      ! REAL(fp)    :: RATE_RUBISCO
      ! REAL(fp)    :: A_GROSS
      REAL(fp)    :: TAU
      REAL(fp)    :: DENOM
      INTEGER     :: ITER
      ! REAL(fp)    :: ERR1
      ! REAL(fp)    :: ERR2
      ! REAL(fp)    :: ERR3
      REAL(fp)    :: DELTA
      REAL(fp)    :: RAB
      ! Strings
      CHARACTER(LEN=255) :: ErrMsg, ThisLoc

      !=================================================================
      ! DO_PHOTOSYNTHESIS begins here!
      !=================================================================
      ! Initialize output variables
      G_CAN_OUT      = 0.e+0_fp
      A_CAN_OUT      = 0.e+0_fp
      RESP_CAN_OUT   = 0.e+0_fp
      G_LEAF_OUT     = 0.e+0_fp
      CO2_IN         = 0.e+0_fp
      A_NET_OUT      = 0.e+0_fp
      RESP_OUT       = 0.e+0_fp
      FLUXO3_CAN     = 0.e+0_fp
      FLUXO3         = 0.e+0_fp
      FACTOR_O3      = 1.e+0_fp
      BETA           = 1.e+0_fp
      V_CMAX         = 0.e+0_fp
      RATE_LIGHT     = 0.e+0_fp
      RATE_RUBISCO   = 0.e+0_fp
      RATE_PRODUCT   = 0.e+0_fp
      A_GROSS        = 0.e+0_fp

      ! Initialize local variables
      TEMPC          = 0.e+0_fp
      SPHU_SAT       = 0.e+0_fp
      DEFICIT_Q      = 0.e+0_fp
      CO2_AMBIENT    = 0.e+0_fp
      O3_CONC        = 0.e+0_fp
      BIGLEAFSCALE   = 0.e+0_fp
      G_CAN          = 0.e+0_fp
      G_LEAF         = 0.e+0_fp
      A_NET          = 0.e+0_fp
      RESP           = 0.e+0_fp 
      CO2_GAMMA      = 0.e+0_fp
      TAU            = 0.e+0_fp
      DENOM          = 0.e+0_fp
      DELTA          = 0.e+0_fp

      ! Initialize
      RC         = GC_SUCCESS
      ErrMsg     = ''
      ThisLoc    = ' -> at DO_PHOTOSYNTHESIS (in GeosCore/ecophy_mod.F90)'

      ! Calculations start.
      TEMPC          = TEMPK - 273.15e+0_fp
      ! Calculate V_CMAX if V_CMAX_IN is not provided
      IF ( PRESENT( V_CMAX_IN ) ) THEN
         V_CMAX = V_CMAX_IN 
      ELSE 
         ! Calculate V_CMAX and respiration which depends on V_CMAX only
         DENOM          = ( 1.e+0_fp + EXP( 0.3e+0_fp*( TEMPC - T_UPP(PFT) ) ) ) &
                        * ( 1.e+0_fp + EXP( 0.3e+0_fp*( T_LOW(PFT) - TEMPC ) ) )
         V_CMAX         = V_CMAX25(PFT) * FACTOR_Q10( 2.e+0_fp, TEMPC ) / DENOM
      END IF
      RESP           = F_DARKRESP(PFT) * V_CMAX

      ! Calculate CO2 compensation point
      TAU            = 2.6e+3_fp * FACTOR_Q10( 0.57e+0_fp, TEMPC )    ! CO2/O2 Specificity Ratio
      CO2_GAMMA      = IS_C3_PLANT(PFT) / ( 2.e+0_fp * TAU ) &
                     * PRESSURE * O2

      ! Calculate canopy scaling factor
      BIGLEAFSCALE   = MAX( ( 1.e+0_fp - EXP( - K_EXTINCT(PFT) * LAI ) ) / K_EXTINCT(PFT), 0.e+0_fp )

      ! Calculate CO2 partial pressure in ambient air
      CO2_AMBIENT    = PRESSURE * CO2

      ! Calculate inputs for ozone damage scheme if FACTOR_O3_IN is not provided
      IF ( .NOT. PRESENT( FACTOR_O3_IN ) ) THEN
         IF ( .NOT. PRESENT( RA    ) ) RETURN
         IF ( .NOT. PRESENT( RB_O3 ) ) RETURN
         IF ( .NOT. PRESENT( O3    ) ) RETURN
         ! Calculate O3 molar concentration in canopy layer
         O3_CONC        = O3 * PRESSURE / RSTARG / TEMPK * 1.e+9_fp
         ! Calculate aerodynamic and boundary layer resistance for ozone damage 
         RAB            = RA + RB_O3
      END IF

      ! Calculate soil moisture stress factor if BETA_IN is not provided.
      IF ( PRESENT( BETA_IN ) ) THEN
         BETA = BETA_IN
      ELSE
         IF ( .NOT. PRESENT( SOIL_WETNESS ) ) RETURN
         ! To modify net photosynthesis rate by soil moisture stress later
         ! Not needed to be inside the loop
         CALL MOIST_STRESS( SOIL_WETNESS, BETA )
      END IF

      ! Iterate to find a self-consistent set of photosynthesis,
      ! stomatal conductance and leaf internal CO2 concentration
      ! Initial guess: G_LEAF = 0 and other initializations
      ITER           = 1
      G_LEAF         = 0.e+0_fp
      G_CAN          = 0.e+0_fp
      ! CO2_IN_PREV    = 0.e+0_fp
      ! A_NET_PREV     = 0.e+0_fp
      ! G_LEAF_PREV    = 0.e+0_fp
      ! ERR1           = 1.e+0_fp
      ! ERR2           = 1.e+0_fp
      ! ERR3           = 1.e+0_fp
      DELTA          = 1.e+0_fp
      DO WHILE ( DELTA >= THRESHOLD .AND. ITER <= MAX_ITER )
         ! Step 1: Calculate internal CO2 partial pressure from 
         !         the closure condition by Jacobs (1994)
         SPHU_SAT    = 0.622e+0_fp * E_SAT( TEMPC ) / PRESSURE
         G_CAN       = G_LEAF * BIGLEAFSCALE
         ! Calculate specific humidity deficit if DEFICIT_Q_IN is not provided
         ! Remark: Consider moving this out of the (potential) loop. Currently 
         ! the total iteration of the loop is set to 1 only.
         IF ( PRESENT( DEFICIT_Q_IN ) ) THEN
            DEFICIT_Q = DEFICIT_Q_IN
         ELSE 
            IF ( .NOT. PRESENT( QV2M ) ) RETURN
            DEFICIT_Q   = MAX( SPHU_SAT - QV2M, 0.e+0_fp )
         END IF
         CO2_IN      = CO2_GAMMA + f0(PFT)*( 1 - DEFICIT_Q / D_STAR(PFT) ) &
                     * ( CO2_AMBIENT - CO2_GAMMA )
         IF ( BETA <= 0.e+0_fp .OR. DEFICIT_Q >= D_STAR(PFT) & 
                               .OR. APAR <= 0.e+0_fp ) THEN
            ! Close stomata if the above conditions are satisfied
            A_NET_OUT   = - RESP * BETA
            RESP_OUT    = RESP
            G_LEAF_OUT  = G_LEAF_MIN(PFT)
            EXIT
         ELSE
            ! Step 2: Photosynthesis model by Collatz et al. (1991) and
            !         Collatz et al. (1992)
            CALL PHOTOSYNTHESIS_LIMITS( CO2_IN,       CO2_GAMMA,    &
                                        O2, APAR,     PRESSURE,     &
                                        TEMPC,        V_CMAX,       &
                                        PFT,          RATE_LIGHT,   &
                                        RATE_PRODUCT, RATE_RUBISCO  )
            CALL SOLVE_COLIMIT( RATE_LIGHT,   RATE_PRODUCT, &
                                RATE_RUBISCO, A_GROSS      )
            ! Calculate net photosynthesis
            A_NET    = ( A_GROSS - RESP ) * BETA

            ! Step 3: Calculate leaf-level stomatal conductance by 
            !         considering diffusive CO2 flux thru open stomata
            CALL LEAF_CONDUCTANCE( A_NET, CO2_AMBIENT, CO2_IN,  &
                                   TEMPK, G_LEAF                )
            ! Close stomata if net photosynthesis <= 0 or
            ! stomatal conductance is too small
            IF ( A_NET <= 0.e+0_fp .OR. G_LEAF <= G_LEAF_MIN(PFT) ) THEN
               A_NET_OUT   = - RESP * BETA
               RESP_OUT    = RESP
               G_LEAF_OUT  = G_LEAF_MIN(PFT)
               EXIT
            END IF

            ! Calculate ozone damage if FACTOR_O3_IN is not provided
            IF ( PRESENT ( FACTOR_O3_IN ) ) THEN 
               FACTOR_O3   = FACTOR_O3_IN 
               A_NET_OUT   = FACTOR_O3 * A_NET
               RESP_OUT    = FACTOR_O3 * RESP
               G_LEAF_OUT  = FACTOR_O3 * G_LEAF
                  ! Close stomata if net photosynthesis <= 0 or
                  ! stomatal conductance is too small
                  IF ( A_NET_OUT <= 0.e+0_fp .OR. G_LEAF_OUT <= G_LEAF_MIN(PFT) ) THEN
                     A_NET_OUT   = - RESP * BETA
                     RESP_OUT    = RESP
                     G_LEAF_OUT  = G_LEAF_MIN(PFT)
                     EXIT
                  END IF
            ELSE 
               ! Apply ozone damage scheme by Sitch et al. (2007)
               SELECT CASE( O3dmg_opt )
               CASE( 'LOW', 'HIGH' )
                  CALL OZONE_DAMAGE ( O3_CONC,   RAB,         &
                                      G_LEAF,    PFT,         &
                                      FLUXO3,    FACTOR_O3,   &
                                      O3dmg_opt, RC           )
                  A_NET_OUT   = FACTOR_O3 * A_NET
                  RESP_OUT    = FACTOR_O3 * RESP
                  G_LEAF_OUT  = FACTOR_O3 * G_LEAF
                     ! Close stomata if net photosynthesis <= 0 or
                     ! stomatal conductance is too small
                     IF ( A_NET_OUT <= 0.e+0_fp .OR. G_LEAF_OUT <= G_LEAF_MIN(PFT) ) THEN
                        A_NET_OUT   = - RESP * BETA
                        RESP_OUT    = RESP
                        G_LEAF_OUT  = G_LEAF_MIN(PFT)
                        EXIT
                     END IF
               CASE( 'OFF' )
                  A_NET_OUT   = A_NET
                  RESP_OUT    = RESP
                  G_LEAF_OUT  = G_LEAF
               CASE DEFAULT 
                  ErrMsg = 'No ozone damage option chosen. Please ' // &
                           'choose from HIGH, LOW or OFF.'
                  CALL GC_Error( ErrMsg, RC, ThisLoc )
                  RETURN               
               END SELECT   ! O3 damage
            END IF   ! FACTOR_O3_IN is present 
         END IF      ! Open or closed stomata
         ! IF ( ITER >= 2 ) THEN
         ! ! calculate error from step 2 onwards
         !    ERR1  = REL_ERR( G_LEAF, G_LEAF_PREV )
         !    ERR2  = REL_ERR( CO2_IN, CO2_IN_PREV )
         !    ERR3  = REL_ERR( A_NET,  A_NET_PREV  )
         !    DELTA = ABS( MAX( ERR1, ERR2, ERR3 ) )
         ! END IF
         ! CO2_IN_PREV  = CO2_IN
         ! A_NET_PREV   = A_NET
         ! G_LEAF_PREV  = G_LEAF
         ITER = ITER + 1
      END DO  ! Do while loop

      ! Canopy scaling: scale leaf-level net photosynthesis, 
      !                 bulk stomatal conductance, repiration and 
      !                 stomatal ozone flux to canopy-level
      A_CAN_OUT      = BIGLEAFSCALE * A_NET_OUT
      G_CAN_OUT      = BIGLEAFSCALE * G_LEAF_OUT
      RESP_CAN_OUT   = BIGLEAFSCALE * RESP_OUT
      FLUXO3_CAN     = BIGLEAFSCALE * FLUXO3

      ! Return to calling program
      END SUBROUTINE DO_PHOTOSYNTHESIS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_isop_emis
!
! !DESCRIPTION: Subroutine DO\_ISOP\_EMIS calculates the photosyntheis-
!  dependent isoprene emission based on Pacifico et al. (2011) 
!  Atm. Chem. Phys.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_ISOP_EMIS( A_CAN_OUT, RESP_CAN_OUT,  &
                               TEMPK, CO2_IN,            &
                               LAI, PFT,                 &
                               ISOP_EMIS                 )
!
! !INPUT PARAMETERS: 
!
      !---------------------------------------------------------------------------------------
      ! A_CAN_OUT     : Canopy net photosynthetic rate      [mol CO2 m^-2 s^-1]
      ! RESP_CAN_OUT  : Canopy dark respiration             [mol CO2 m^-2 s^-1]
      ! TEMPK         : Leaf temperature in Kelvin          [K]
      ! CO2_IN        : Leaf partial pressure of CO2        [Pa]
      ! LAI           : Leaf area index for the PFT         [m^2 m^-2]
      ! PFT           : Index for PFT                       []
      !---------------------------------------------------------------------------------------
      REAL(fp), INTENT(IN)    :: A_CAN_OUT
      REAL(fp), INTENT(IN)    :: RESP_CAN_OUT
      REAL(fp), INTENT(IN)    :: TEMPK
      REAL(fp), INTENT(IN)    :: CO2_IN
      REAL(fp), INTENT(IN)    :: LAI
      INTEGER,  INTENT(IN)    :: PFT
!
! !OUTPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! ISOP_EMIS     : Rate of isoprene emission           [kg C m-2 land s-1]
      !---------------------------------------------------------------------------------------
      REAL(fp), INTENT(OUT)    :: ISOP_EMIS
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
!
! !LOCAL VARIABLES:
!
      !---------------------------------------------------------------------------------------
      ! FACTOR_CO2    : CO2 factor for online isoprene emission           [] 
      ! FACTOR_TEMP   : Temparture factor for online isoprene emission    [] 
      ! Dry_weight    : Dry weight of leaf per leaf area                  [g dry weight m-2 leaf]
      !---------------------------------------------------------------------------------------
      REAL(fp)    :: FACTOR_CO2
      REAL(fp)    :: FACTOR_TEMP
      REAL(fp)    :: Dry_weight
      
      !=================================================================
      ! DO_ISOP_EMIS begins here!
      !=================================================================
      ! Initialize local variables
      FACTOR_CO2     = 0.e+0_fp
      FACTOR_TEMP    = 0.e+0_fp

      ! Calculate CO2 and temperature factors for photosynthesis-dependent 
      ! isoprene Emission
      FACTOR_CO2     = CO2_IN_st(PFT) / CO2_IN 
      FACTOR_TEMP    = MIN( 2.3e+0_fp , EXP( 0.1e+0_fp * (TEMPK - TEMPK_st) ) )
      
      ! Calculate dry weight of leaf per leaf area (g dry weight m-2 leaf)
      Dry_weight     = Leaf_density(PFT) * 1.e+3_fp / Dry_fraction 

      ! Calculate isoprene emission (kg C m-2 land s-1)
      ISOP_EMIS      = IEF(PFT) * Dry_weight * 1.e-9_fp / 3.6e+3_fp  &
                     * ( A_CAN_OUT + RESP_CAN_OUT )                  &
                     / ( A_NET_st(PFT) + RESP_st (PFT) )             &
                     * FACTOR_CO2 * FACTOR_TEMP   
      ISOP_EMIS      = MAX( 0.e+0_fp, ISOP_EMIS )

      ! Return to calling program
      END SUBROUTINE DO_ISOP_EMIS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: photosynthesis_limits
!
! !DESCRIPTION: Subroutine PHOTOSYNTHESIS\_LIMITES determines potential 
!  leaf-level photosynthesis, according to C3 and C4 photosynthesis
!  model from Collatz et al. (1991) and Collatz et al. (1992)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE PHOTOSYNTHESIS_LIMITS( CO2_IN,       CO2_GAMMA,    &
                                        O2, APAR,     PRESSURE,     &
                                        TEMPC,        V_CMAX,       &
                                        PFT,          RATE_LIGHT,   &
                                        RATE_PRODUCT, RATE_RUBISCO  )
!
! !INPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! CO2_IN        : Leaf internal partial pressure of CO2              [Pa]
      ! CO2_GAMMA     : CO2 Compensation point                             [Pa]
      ! O2            : Ambient O2 mole fraction                           [mol/mol air]
      ! APAR          : Absorbed photosynthetically active radiation (PAR) [mol photon m^-2 s^-1]
      ! PRESSURE      : Surface air pressure                               [Pa]
      ! TEMPC         : Temperature                                        [deg C]
      ! V_CMAX        : Maximum rate of carboxylation of Rubisco           [mol CO2 m^-2 s^-1]
      ! PFT           : Index for PFT                                      []
      !---------------------------------------------------------------------------------------
      REAL(fp), INTENT(IN)  :: CO2_IN
      REAL(fp), INTENT(IN)  :: CO2_GAMMA
      REAL(fp), INTENT(IN)  :: O2
      REAL(fp), INTENT(IN)  :: APAR
      REAL(fp), INTENT(IN)  :: PRESSURE
      REAL(fp), INTENT(IN)  :: TEMPC
      REAL(fp), INTENT(IN)  :: V_CMAX
      INTEGER,  INTENT(IN)  :: PFT
!
! !OUTPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! RATE_LIGHT         : Light-limited rate                            [mol CO2 m^-2 s^-1]
      ! RATE_PRODUCT       : Product-limited rate                          [mol CO2 m^-2 s^-1]
      ! RATE_RUBISCO       : Rubisco-limited rate                          [mol CO2 m^-2 s^-1]
      !---------------------------------------------------------------------------------------
      REAL(fp), INTENT(OUT) :: RATE_LIGHT
      REAL(fp), INTENT(OUT) :: RATE_PRODUCT
      REAL(fp), INTENT(OUT) :: RATE_RUBISCO
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Michaelis-Menten parameters for CO2 and O2 respectively
      REAL(fp)              :: K_C
      REAL(fp)              :: K_O

      !=================================================================
      ! PHOTOSYNTHESIS_LIMITS begins here!
      !=================================================================
      K_C               = 3.e+1_fp * FACTOR_Q10( 2.1e+0_fp, TEMPC )
      K_O               = 3.e+4_fp * FACTOR_Q10( 1.2e+0_fp, TEMPC )
      ! For C4 plants
      IF ( IS_C3_PLANT(PFT) == 0 ) THEN
         RATE_RUBISCO   = V_CMAX
         RATE_LIGHT     = ALPHA(PFT) * APAR
         RATE_PRODUCT   = 2.e+4_fp * V_CMAX * CO2_IN / PRESSURE
      ELSE    ! For C3 plants
         RATE_RUBISCO   = V_CMAX * ( CO2_IN - CO2_GAMMA )  &
                        / ( CO2_IN + K_C * ( 1.e+0_fp + PRESSURE * O2 / K_O ) )
         RATE_LIGHT     = ALPHA(PFT) * APAR           &
                        * ( CO2_IN - CO2_GAMMA )      &
                        / ( CO2_IN + CO2_GAMMA * 2.e+0_fp )
         RATE_RUBISCO   = MAX( RATE_RUBISCO, 0.e+0_fp )
         RATE_LIGHT     = MAX( RATE_LIGHT, 0.e+0_fp )
         RATE_PRODUCT   = 0.5e+0_fp * V_CMAX
      END IF

      ! Return to calling program
      END SUBROUTINE PHOTOSYNTHESIS_LIMITS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: solve_colimit
!
! !DESCRIPTION: Subroutine SOLVE\_COLIMIT determines gross leaf-level 
!  photosynthesis from the three unstressed photosynthesis rates, 
!  according to the co-limiting regime in C3 and C4 photosynthesis
!  model from Collatz et al. (1991) and Collatz et al. (1992)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SOLVE_COLIMIT( RATE_LIGHT,   RATE_PRODUCT, &
                                RATE_RUBISCO, A_GROSS       )
!
! !INPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! RATE_LIGHT         : Light-limited rate [mol CO2 m^-2 s^-1]
      ! RATE_PRODUCT       : Product-limited rate [mol CO2 m^-2 s^-1]
      ! RATE_RUBISCO       : Rubisco-limited rate [mol CO2 m^-2 s^-1]
      !---------------------------------------------------------------------------------------
      REAL(fp), INTENT(IN)  :: RATE_LIGHT
      REAL(fp), INTENT(IN)  :: RATE_PRODUCT
      REAL(fp), INTENT(IN)  :: RATE_RUBISCO
!
! !OUTPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! A_GROSS            : Gross rate of photosynthesis [mol CO2 m^-2 s^-1]
      !---------------------------------------------------------------------------------------
      REAL(fp), INTENT(OUT) :: A_GROSS
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp), PARAMETER   :: BETA1 = 0.83e+0_fp
      REAL(fp), PARAMETER   :: BETA2 = 0.93e+0_fp
      REAL(fp)              :: TEMP
      REAL(fp)              :: B
      REAL(fp)              :: C

      !=================================================================
      ! SOLVE_COLIMIT begins here!
      !=================================================================
      ! 1st quadratic
      B = -( RATE_RUBISCO + RATE_LIGHT ) / BETA1
      C = RATE_RUBISCO * RATE_LIGHT / BETA1
      ! Note that C > 0, SQRT( B^2 - 4*C ) < ABS(B)
      ! Take smaller root
      TEMP = 0.5e+0_fp * ( - B - SQRT( B * B - 4.e+0_fp * C ) )

      ! 2nd quadratic
      B = - ( TEMP + RATE_PRODUCT ) / BETA2
      C = TEMP * RATE_PRODUCT / BETA2
      ! Note that C > 0, SQRT( B^2 - 4*C ) < ABS(B)
      ! Take smaller root
      A_GROSS  = 0.5e+0_fp * ( - B - SQRT( B * B - 4.e+0_fp * C ) )

      ! Return to calling program
      END SUBROUTINE SOLVE_COLIMIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: moist_stress
!
! !DESCRIPTION: Subroutine MOIST\_STRESS determines a soil moisutre 
!  stress factor that reduces photosynthesis rate and stomatal 
!  conductance. It is currently a linear function of soil moisture.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE MOIST_STRESS( SOIL_WETNESS, BETA )
!
! !INPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! SOIL_WETNESS    : Volumetric mean moisture concentration in root zone divided by porosity
      !---------------------------------------------------------------------------------------
      REAL(fp), INTENT(IN)  :: SOIL_WETNESS
!
! !OUTPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! BETA            : Moisture stress factor
      !---------------------------------------------------------------------------------------
      REAL(fp), INTENT(OUT) :: BETA
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES: (Declared as module variables)
!
      !---------------------------------------------------------------------------------------
      ! SATU            : Volumetric soil moisture at saturation (= porosity)
      ! CRIT            : Volumetric soil moisture at critical point (above which
      !                   plants are no longer stressed by soil moisture)
      ! WILT            : Volumetric soil moisture at wilting point (below which
      !                   photosynthesis is stopped by limited soil moisture)
      !---------------------------------------------------------------------------------------

      !=================================================================
      ! MOIST_STRESS begins here!
      !=================================================================
      BETA =  ( SOIL_WETNESS*SATU - WILT ) / ( CRIT - WILT )
      BETA = MIN( MAX( 0.e+0_fp, BETA ), 1.e+0_fp ) 

      ! Return to calling program
      END SUBROUTINE MOIST_STRESS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: leaf_conductance
!
! !DESCRIPTION: Subroutine LEAF_CONDUCTANCE determines leaf-level stomatal 
!  conductance by considering the diffusive CO2 exchange through stomatal 
!  exchange.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LEAF_CONDUCTANCE( A_NET, CO2_AMBIENT, CO2_IN,  &
                                   TEMPK, G_LEAF                )
!
! !INPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! A_NET             : Net photosynthesis                      [mol CO2 m^-2 s^-1]
      ! CO2_AMBIENT       : Ambient CO2 partial pressure            [Pa]
      ! CO2_IN            : Leaf internal CO2 partial pressure      [Pa]
      ! TEMPK             : Leaf temperature                        [K]
      !---------------------------------------------------------------------------------------
      REAL(fp), INTENT(IN)    :: A_NET
      REAL(fp), INTENT(IN)    :: CO2_AMBIENT
      REAL(fp), INTENT(IN)    :: CO2_IN
      REAL(fp), INTENT(IN)    :: TEMPK
!
! !OUTPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! G_LEAF            : Leaf conductance for H2O                [m s^-1]
      !---------------------------------------------------------------------------------------
      REAL(fp), INTENT(OUT)   :: G_LEAF
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC     

      !=================================================================
      ! LEAF_CONDUCTANCE begins here!
      !=================================================================
      G_LEAF = CO2_O2_RATIO * RSTARG * TEMPK  &
             * A_NET / ( CO2_AMBIENT - CO2_IN )

      ! Return to calling program
      END SUBROUTINE LEAF_CONDUCTANCE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ozone_damage
!
! !DESCRIPTION: Subroutine MOIST\_STRESS determines an ozone damage 
!  factor that reduces photosynthesis rate and stomatal conductance, 
!  based on Sitch et al. (2007).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OZONE_DAMAGE ( O3_CONC,   RAB,       &
                                G_LEAF,    PFT,       &
                                FLUXO3,    FACTOR_O3, &
                                O3dmg_opt, RC         )
!
! !USES:
!
      USE ErrCode_Mod
!
! !INPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! O3_CONC         : Ozone concentration in canopy layer                     [nmol m^-3]
      ! RAB             : Aerodynamic and boundary resistance                     [s m^-1]
      ! G_LEAF          : Leaf conductance for H2O in the absence of O3 effects   [m s^-1]
      ! PFT             : Index for PFT                                           []
      ! O3dmg_opt       : Control switch for ozone damage scheme                  [HIGH/LOW/OFF]
      !---------------------------------------------------------------------------------------
      REAL(fp),         INTENT(IN)  :: O3_CONC
      REAL(fp),         INTENT(IN)  :: RAB
      REAL(fp),         INTENT(IN)  :: G_LEAF
      INTEGER,          INTENT(IN)  :: PFT
      CHARACTER(LEN=4), INTENT(IN)  :: O3dmg_opt
!
! !OUTPUT PARAMETERS:
!  
      !---------------------------------------------------------------------------------------
      ! FLUXO3          : Leaf uptake of O3                                       [nmol m^-2 s^-1]
      ! FACTOR_O3       : Ozone damage factor                                     []
      ! RC              : Return Code                                             []
      !---------------------------------------------------------------------------------------
      REAL(fp),         INTENT(OUT) :: FLUXO3
      REAL(fp),         INTENT(OUT) :: FACTOR_O3
      INTEGER,          INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!    
      REAL(fp) :: B       ! Second coefficient in quadratic equation
      REAL(fp) :: C       ! Constant in quadratic equation
      REAL(fp) :: TEMP1
      REAL(fp) :: TEMP2
      REAL(fp) :: F
      REAL(fp) :: PARAM_A
      ! Strings
      CHARACTER(LEN=255) :: ErrMsg,    ThisLoc

      !=================================================================
      ! OZONE_DAMAGE begins here!
      !=================================================================
      ! Initialize
      RC         = GC_SUCCESS
      ErrMsg     = ''
      ThisLoc    = ' -> at OZONE_DAMAGE (in GeosCore/ecophy_mod.F90)'
      ! Choose ozone damage parameters based on sensistivity
      SELECT CASE( O3dmg_opt )
         CASE( 'HIGH' )
            PARAM_A = PARAM_A_HIGH( PFT )
         CASE( 'LOW' )
            PARAM_A = PARAM_A_LOW( PFT )
         CASE DEFAULT
            PARAM_A = 0
            ErrMsg = 'Any O3dmg_opt other than "HIGH" or "LOW" should ' // &
                     'not have reached this subroutine.'
            CALL GC_Error( ErrMsg, RC, ThisLoc )
            RETURN
      END SELECT
      TEMP1       = 1.e+0_fp + PARAM_A * FLUXO3_CRIT(PFT)
      TEMP2       = 1.61e+0_fp / G_LEAF

      ! Solve for ozone damage factor
      ! Calculate coefficients for quadratic equation F^2 + B*F + C = 0
      IF ( ABS(RAB) < EPSILON(1.e+0_fp) ) THEN
         ! RAB        = 0.0
         FACTOR_O3 = TEMP1 / ( 1.e+0_fp + PARAM_A * O3_CONC / TEMP2 )
      ELSE
         B         = TEMP2 / RAB - TEMP1 + PARAM_A * O3_CONC / RAB
         C         = - TEMP1 * TEMP2 / RAB
         ! Note that C < 0, SQRT( B^2 - 4*C ) > ABS(B)
         ! Take positive root
         F         = 0.5e+0_fp * ( SQRT( B * B - 4.e+0_fp * C ) - B )
         FACTOR_O3 = MIN( MAX( F, 0.e+0_fp ), 1.e+0_fp )
      END IF

      ! Calculate stomatal ozone flux
      FLUXO3      = O3_CONC / ( RAB + TEMP2 / FACTOR_O3 )

      ! Return to calling program
      END SUBROUTINE OZONE_DAMAGE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: factor_q10
!
! !DESCRIPTION: Function FACTOR\_Q10 determines the q10 temperature 
!  sensitivity factor that scales the maximum rubisco carboxylation rate
!  at 25 degree Celsius to other temperatures.
!\\
!\\
! !INTERFACE:
!
      FUNCTION FACTOR_Q10( Q10, TEMPC ) RESULT( FACTOR )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN)    :: Q10         ! Coefficient
      REAL(fp), INTENT(IN)    :: TEMPC       ! Temperature [deg C]
!
! !RETURN VALUE:
!   
      REAL(fp)                :: FACTOR
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
      

      !=================================================================
      ! FACTOR_Q10 begins here!
      !=================================================================
      FACTOR = Q10**( 0.1e+0_fp * ( TEMPC - 25.e+0_fp ) )

      ! Return to calling program
      END FUNCTION FACTOR_Q10
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rel_err
!
! !DESCRIPTION: Function REL\_ERR calculates the relative error of a 
!  quantity between two iterations.
!\\
!\\
! !INTERFACE:
!
      FUNCTION REL_ERR( ITEM, ITEM_PREV ) RESULT( ERROR )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN)    :: ITEM
      REAL(fp), INTENT(IN)    :: ITEM_PREV
!
! !RETURN VALUE:
!
      REAL(fp)                :: ERROR
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! REL_ERR begins here!
      !=================================================================
      ERROR = ABS( ITEM - ITEM_PREV ) / ABS( ITEM_PREV )

      !Return to calling program
      END FUNCTION REL_ERR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: e_sat
!
! !DESCRIPTION: Function E\_SAT calculates the saturation vapour 
!  pressure (in Pa) using the empirical formula by Lowe and Ficke (1974)
!\\
!\\
! !INTERFACE:
!
      FUNCTION E_SAT( TEMPC ) RESULT ( Esat )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN)  :: TEMPC
!
! !RETURN VALUE:
!
      REAL(fp)              :: Esat
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp), PARAMETER   :: a0 = 6.107799961e+0_fp
      REAL(fp), PARAMETER   :: a1 = 4.436518521e-1_fp
      REAL(fp), PARAMETER   :: a2 = 1.428945805e-2_fp
      REAL(fp), PARAMETER   :: a3 = 2.650648471e-4_fp
      REAL(fp), PARAMETER   :: a4 = 3.031240396e-6_fp
      REAL(fp), PARAMETER   :: a5 = 2.034080948e-8_fp
      REAL(fp), PARAMETER   :: a6 = 6.136820929e-11_fp

      !=================================================================
      ! E_SAT begins here!
      !=================================================================
      Esat = 1.e+2_fp * ( a0 + TEMPC * ( a1 + TEMPC * ( a2 + TEMPC   &
                      * ( a3 + TEMPC * ( a4 + TEMPC * ( a5 + TEMPC * a6 ) ) ) ) ) )

      ! Return to calling program
      END FUNCTION E_SAT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_ecophy_inputs
!
! !DESCRIPTION: Subroutine GET\_ECOPHY_INPUTS get inputs
!  from Met and Chem State objects and Input Options.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_ECOPHY_INPUTS( State_Met,    State_Chm, Input_Opt,&
                                    I, J, LDT,                         &
                                    TEMPK,        QV2M,                &
                                    PAR_ABSORBED, CO2,                 &
                                    O2,           LAI,       O3,       &
                                    SOIL_WETNESS, IUSE                 &
                                    )
!
! !USES:
!
      USE State_Chm_Mod,      ONLY : ChmState, Ind_
      USE State_Met_Mod,      ONLY : MetState
      USE Input_Opt_Mod,      ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! State_Met     : Meteorology State Object
      ! State_Chm     : Chemistry State Object
      ! Input_Opt.    : Input Options Object
      ! I             : Current lon index
      ! J             : Current lat index
      ! LDT           : Land type index
      !---------------------------------------------------------------------------------------
      Type(MetState), INTENT(IN)  :: State_Met
      Type(ChmState), INTENT(IN)  :: State_Chm
      TYPE(OptInput), INTENT(IN)  :: Input_Opt 
      INTEGER,        INTENT(IN)  :: I
      INTEGER,        INTENT(IN)  :: J
      INTEGER,        INTENT(IN)  :: LDT
!
! !OUTPUT PARAMETERS:
!
      !---------------------------------------------------------------------------------------
      ! TEMPK         : Temperature in Kelvin                             [K]
      ! QV2M          : Specific humidity in canopy layer                 [kg H2O / kg air]
      ! PAR_ABSORBED  : Absorbed PAR                                      [mol photon m^-2 s^-1]
      ! PRESSURE      : Atmospheric Pressure in canopy layer              [Pa]
      ! CO2           : Ambient CO2 mole fraction                         [kg / kg dry]
      ! O2            : Ambient O2 mole fraction                          [kg / kg dry]
      ! O3            : Ozone mole fraction in canopy layer               [kg / kg dry]
      ! LAI           : Leaf area index for the PFT                       [m^2 m^-2]
      ! SOIL_WETNESS  : Fraction of moisture in soil pores                []
      ! IUSE          : Fraction(*1000) of grid occupied by land type LDT []
      !---------------------------------------------------------------------------------------
      REAL(fp), INTENT(OUT) :: TEMPK
      REAL(fp), INTENT(OUT) :: QV2M
      REAL(fp), INTENT(OUT) :: PAR_ABSORBED
      ! REAL(fp), INTENT(OUT) :: PRESSURE
      REAL(fp), INTENT(OUT) :: CO2
      REAL(fp), INTENT(OUT) :: O2
      REAL(fp), INTENT(OUT) :: O3
      REAL(fp), INTENT(OUT) :: LAI
      REAL(fp), INTENT(OUT) :: SOIL_WETNESS
      INTEGER,  INTENT(OUT) :: IUSE
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: PARDR, PARDF, ALBD
      ! INTEGER :: I, J
      INTEGER  :: id_O2, id_O3

      !=================================================================
      ! get_ecophy_inputs begins here!
      !=================================================================
      ! Find tracer indices with function the Ind_() function
      id_O2     = IND_( 'O2'  )
      id_O3     = IND_( 'O3'  )

      ! ! Loop over surface grid boxes
      ! DO J = 1, JJPAR
      ! DO I = 1, IIPAR
      ! Surface temperature [K]
      TEMPK         = State_Met%TS( I,J )
      ! Specific humidity [kg/kg]
      QV2M          = State_Met%QV2M( I,J )
      ! Photosynthetically active radiation absorbed [W m^-2]
      PARDR         = State_Met%PARDR( I,J )
      PARDF         = State_Met%PARDF( I,J )
      ALBD          = State_Met%ALBD ( I,J )
      PAR_ABSORBED  = ( 1 - ALBD ) * ( PARDR + PARDF ) * 4.6e-6_fp
      ! ! Pressure [Pa]
      ! PRESSURE      = State_Met%SLP( I,J ) * 1.e+2_fp
      ! CO2 mole fraction [mol/mol]
      ! CO2           = State_Chm%Species( I,J,1,id_CO2 ) * AIRMW &
      !               / State_Chm%SpcData( id_CO2 )%Info%MW_g
      CO2           = Input_Opt%Ecophy_CO2 * 1.e-6_fp
      ! O2 mole fraction [mol/mol]
      O2            = State_Chm%Species( I,J,1,id_O2  ) * AIRMW &
                    / State_Chm%SpcData( id_O2  )%Info%MW_g
      ! O3 mole fraction [mol/mol]
      O3            = State_Chm%Species( I,J,1,id_O3  ) * AIRMW &
                    / State_Chm%SpcData( id_O3  )%Info%MW_g
      ! LAI [m^2 m^-2]
      LAI           = State_Met%XLAI( I,J,LDT )
      ! Root zone soil wetness
      SOIL_WETNESS  = State_Met%GWETROOT( I,J )
      ! Soil moisture at saturation
      SATU          = State_Met%THETA_SATU( I,J )
      ! Soil moisture at critical point
      CRIT          = State_Met%THETA_CRIT( I,J )
      ! Soil moisture at wilting point
      WILT          = State_Met%THETA_WILT( I,J )
      ! NOTE: Leave possibility to call soil matric potentials instead
      !       of soil moistures
      ! Fraction(*1000) of grid occupied by land type LDT 
      ! IUSE=500 => 50% of the grid 
      IUSE          = State_Met%IUSE( I,J,LDT )
       ! END DO
       ! END DO

      ! Return to calling program
      END SUBROUTINE GET_ECOPHY_INPUTS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ecophy_diagn
!
! !DESCRIPTION: Subroutine ECOPHY\_DIAGN deals with saving ecophysiology
!  module outputs (per land type) to State_Diag arrays (per PFT)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Ecophy_Diagn( I, J,         LDT,        PFT,          &
                               IUSE,         LAI,        SumLAI_PFT,   &
                               IUSE_PFT,     RA,         RB_O3,        &
                               G_CAN_OUT,    A_CAN_OUT,  RESP_CAN_OUT, &
                               CO2_IN,       FLUXO3_CAN,               &
                               FACTOR_O3,    BETA,       V_CMAX,       &
                               RATE_LIGHT,   RATE_RUBISCO,             &
                               RATE_PRODUCT, A_GROSS,    ISOP_EMIS,    &
                               State_Met,    State_Diag, RC            &
                             )
!
! !USES:
!
      USE ErrCode_Mod
      USE State_Met_Mod,  ONLY : MetState
      USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)    :: I, J, LDT, PFT
      INTEGER,        INTENT(IN)    :: IUSE
      REAL(fp),       INTENT(IN)    :: LAI
      REAL(fp),       INTENT(IN)    :: SumLAI_PFT
      INTEGER,        INTENT(IN)    :: IUSE_PFT
      REAL(fp),       INTENT(IN)    :: RA
      REAL(fp),       INTENT(IN)    :: RB_O3
      REAL(fp),       INTENT(IN)    :: G_CAN_OUT
      REAL(fp),       INTENT(IN)    :: A_CAN_OUT
      REAL(fp),       INTENT(IN)    :: RESP_CAN_OUT
      REAL(fp),       INTENT(IN)    :: CO2_IN
      REAL(fp),       INTENT(IN)    :: FLUXO3_CAN 
      REAL(fp),       INTENT(IN)    :: FACTOR_O3
      REAL(fp),       INTENT(IN)    :: BETA
      REAL(fp),       INTENT(IN)    :: V_CMAX
      REAL(fp),       INTENT(IN)    :: RATE_LIGHT
      REAL(fp),       INTENT(IN)    :: RATE_RUBISCO
      REAL(fp),       INTENT(IN)    :: RATE_PRODUCT
      REAL(fp),       INTENT(IN)    :: A_GROSS
      REAL(fp),       INTENT(IN)    :: ISOP_EMIS
      Type(MetState), INTENT(IN)    :: State_Met   ! Met State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                :: IOLSON
      REAL(fp)               :: Tmp
      ! Strings
      CHARACTER(LEN=255)     :: Msg, ErrMsg, ThisLoc

      !=================================================================
      ! Ecophy_Diagn begins here!
      !=================================================================
      ! Initialize
      RC        = GC_SUCCESS
      ErrMsg    = ''
      ThisLoc   = ' -> at Ecophy_Diagn (in module GeosCore/ecophysiology.F90)'

      ! send to diagnostics outputs
      IF ( State_Diag%Archive_EcophyPFTFrac      .AND. SumLAI_PFT /= 0 ) THEN
         State_Diag%EcophyPFTFrac      ( I,J,PFT ) = DBLE( IUSE_PFT ) 
      END IF
      IF ( State_Diag%Archive_EcophyLAI          .AND. IUSE_PFT /= 0 ) THEN
         State_Diag%EcophyLAI          ( I,J,PFT ) = SumLAI_PFT / DBLE( IUSE_PFT )
      END IF
      IF ( State_Diag%Archive_EcophyRa           .AND. SumLAI_PFT /= 0 ) THEN
         State_Diag%EcophyRa     ( I,J ) = RA
      END IF
      IF ( State_Diag%Archive_EcophyRbO3        .AND. SumLAI_PFT /= 0 ) THEN
         State_Diag%EcophyRbO3  ( I,J ) = RB_O3
      END IF
      !-----------------------------------------------------------------
      ! Canopy-level diagnostics: weight by land area occupied by the PFT
      !-----------------------------------------------------------------
      ! G_CAN_OUT is archived in drydep_mod.F instead of here.
      IF ( State_Diag%Archive_EcophyACan    .AND. SumLAI_PFT /= 0 ) THEN
         Tmp = State_Diag%EcophyACan    ( I,J,PFT )
         State_Diag%EcophyACan    ( I,J,PFT ) = Tmp    &
            + A_CAN_OUT    * DBLE( IUSE ) / DBLE( IUSE_PFT )
      END IF
      IF ( State_Diag%Archive_EcophyRespCan .AND. SumLAI_PFT /= 0 ) THEN
         Tmp = State_Diag%EcophyRespCan ( I,J,PFT )
         State_Diag%EcophyRespCan ( I,J,PFT ) = Tmp    &
            + RESP_CAN_OUT * DBLE( IUSE ) / DBLE( IUSE_PFT )
      END IF
      IF ( State_Diag%Archive_EcophyFluxO3Can   .AND. SumLAI_PFT /= 0 ) THEN
         Tmp = State_Diag%EcophyFluxO3Can   ( I,J,PFT )
         State_Diag%EcophyFluxO3Can   ( I,J,PFT ) = Tmp    &
            + FLUXO3_CAN   * DBLE( IUSE ) / DBLE( IUSE_PFT )
      END IF
      IF ( State_Diag%Archive_EcophyIsopEmisPFT    .AND. SumLAI_PFT /= 0 ) THEN
         Tmp = State_Diag%EcophyIsopEmisPFT    ( I,J,PFT )
         State_Diag%EcophyIsopEmisPFT    ( I,J,PFT ) = Tmp    &
            + ISOP_EMIS    * DBLE( IUSE ) / DBLE( IUSE_PFT )
      END IF
      !-----------------------------------------------------------------
      ! Leaf-level diagnostics: weight by leaf area
      !-----------------------------------------------------------------
      IF ( State_Diag%Archive_EcophyCO2Leaf       .AND. SumLAI_PFT /= 0 ) THEN
         Tmp = State_Diag%EcophyCO2Leaf       ( I,J,PFT )
         State_Diag%EcophyCO2Leaf       ( I,J,PFT ) = Tmp    &
            + CO2_IN       * DBLE( IUSE ) * LAI / SumLAI_PFT
      END IF
      IF ( State_Diag%Archive_EcophyO3DmgFac    .AND. SumLAI_PFT /= 0 ) THEN
         Tmp = State_Diag%EcophyO3DmgFac    ( I,J,PFT )
         State_Diag%EcophyO3DmgFac    ( I,J,PFT ) = Tmp    &
            + FACTOR_O3    * DBLE( IUSE ) * LAI / SumLAI_PFT
      END IF
      IF ( State_Diag%Archive_EcophySoilStress   .AND. SumLAI_PFT /= 0 ) THEN
         Tmp = State_Diag%EcophySoilStress   ( I,J,PFT )
         State_Diag%EcophySoilStress   ( I,J,PFT ) = Tmp    &
            + BETA         * DBLE( IUSE ) * LAI / SumLAI_PFT
      END IF
      IF ( State_Diag%Archive_EcophyVCMax       .AND. SumLAI_PFT /= 0 ) THEN
         Tmp = State_Diag%EcophyVCMax       ( I,J,PFT )
         State_Diag%EcophyVCMax       ( I,J,PFT ) = Tmp    &
            + V_CMAX       * DBLE( IUSE ) * LAI / SumLAI_PFT
      END IF
      IF ( State_Diag%Archive_EcophyLightLmtRate   .AND. SumLAI_PFT /= 0 ) THEN
         Tmp = State_Diag%EcophyLightLmtRate   ( I,J,PFT )
         State_Diag%EcophyLightLmtRate   ( I,J,PFT ) = Tmp    &
            + RATE_LIGHT   * DBLE( IUSE ) * LAI / SumLAI_PFT
      END IF
      IF ( State_Diag%Archive_EcophyRubisLmtRate .AND. SumLAI_PFT /= 0 ) THEN
         Tmp = State_Diag%EcophyRubisLmtRate ( I,J,PFT )
         State_Diag%EcophyRubisLmtRate ( I,J,PFT ) = Tmp    &
            + RATE_RUBISCO * DBLE( IUSE ) * LAI / SumLAI_PFT
      END IF
      IF ( State_Diag%Archive_EcophyProdLmtRate .AND. SumLAI_PFT /= 0 ) THEN
         Tmp = State_Diag%EcophyProdLmtRate ( I,J,PFT )
         State_Diag%EcophyProdLmtRate ( I,J,PFT ) = Tmp    &
            + RATE_PRODUCT * DBLE( IUSE ) * LAI / SumLAI_PFT
      END IF
      IF ( State_Diag%Archive_EcophyAGrossLeaf      .AND. SumLAI_PFT /= 0 ) THEN
         Tmp = State_Diag%EcophyAGrossLeaf      ( I,J,PFT )
         State_Diag%EcophyAGrossLeaf      ( I,J,PFT ) = Tmp    &
            + A_GROSS      * DBLE( IUSE ) * LAI / SumLAI_PFT
      END IF

      ! Return to calling program 
      END SUBROUTINE Ecophy_Diagn
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_ecophy
!
! !DESCRIPTION: Subroutine INIT\_ECOPHY initializes certain variables for the
!  GEOS-CHEM ecophysiology subroutines and calculates the photosynthesis under 
!  standard conditions.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_ECOPHY( Input_Opt, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE Charpak_Mod,    ONLY : To_UpperCase
      USE Input_Opt_Mod,  ONLY : OptInput
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: PFT
      REAL(fp) :: LAI            = 1.e+0_fp
      ! Standard conditions 
      REAL(fp) :: PRESSURE       = 1.01325e+5_fp
      REAL(fp) :: APAR           = 1000.e-6_fp
      REAL(fp) :: CO2            = 370.e-6_fp
      REAL(fp) :: O2             = 0.209e+0_fp
      REAL(fp) :: TEMPK          = 303.15e+0_fp 
      REAL(fp) :: DEFICIT_Q_IN   = 0.e+0_fp
      REAL(fp) :: FACTOR_O3_IN   = 1.e+0_fp
      REAL(fp) :: BETA_IN        = 1.e+0_fp
      REAL(fp) :: G_CAN_OUT
      REAL(fp) :: G_LEAF_OUT
      REAL(fp) :: CO2_IN
      REAL(fp) :: A_CAN_OUT
      REAL(fp) :: A_NET_OUT
      REAL(fp) :: RESP_CAN_OUT
      REAL(fp) :: RESP_OUT
      REAL(fp) :: FLUXO3_CAN
      REAL(fp) :: FLUXO3
      REAL(fp) :: FACTOR_O3
      REAL(fp) :: BETA
      REAL(fp) :: V_CMAX
      REAL(fp) :: RATE_LIGHT
      REAL(fp) :: RATE_RUBISCO
      REAL(fp) :: RATE_PRODUCT
      REAL(fp) :: A_GROSS

      ! Strings
      CHARACTER(LEN=255)     :: Msg, ErrMsg, ThisLoc

      !=================================================================
      ! INIT_ECOPHY begins here!
      !=================================================================

      ! Initialize
      RC        = GC_SUCCESS
      O3dmg_opt = To_Uppercase( TRIM( Input_Opt%O3dmg_opt ) )
      ErrMsg    = ''
      ThisLoc   = ' -> at Init_Ecophy (in module GeosCore/ecophysiology.F90)'

      IF ( Input_Opt%LIsop_from_Ecophy ) THEN
         DO PFT = 1, NUMPFT
            ! simulate plant processes under standard conditions
            CALL DO_PHOTOSYNTHESIS( LAI, PFT,     PRESSURE,     APAR,         &
                                    CO2, O2,      TEMPK,        G_CAN_OUT,    &
                                    G_LEAF_OUT,   CO2_IN,       A_CAN_OUT,    &
                                    A_NET_OUT,    RESP_CAN_OUT, RESP_OUT,     &
                                    FLUXO3_CAN,   FLUXO3,       FACTOR_O3,    &
                                    BETA,         V_CMAX,       RATE_LIGHT,   &
                                    RATE_RUBISCO, RATE_PRODUCT, A_GROSS, RC,  &
                                    DEFICIT_Q_IN=DEFICIT_Q_IN,                &
                                    FACTOR_O3_IN=FACTOR_O3_IN,                &
                                    BETA_IN=BETA_IN                           &
                                    )
            A_NET_st    (PFT)  = A_NET_OUT
            RESP_st     (PFT)  = RESP_OUT
            CO2_IN_st   (PFT)  = CO2_IN
         END DO
         WRITE( 6, *   ) 'Ecophysiology variables under standard conditions'
         WRITE( 6, *   ) 'A_NET_st :   Canopy Photosynthesis (mol CO2 m^-2 s^-1) '
         WRITE( 6, *   ) 'RESP_st  :   Canopy Respiration    (mol CO2 m^-2 s^-1) '
         WRITE( 6, *   ) 'CO2_IN_st:   Leaf CO2 Pressure     (Pa)                '
         WRITE( 6, *   ) 'A_NET_st,      RESP_st,       CO2_IN_st                '
         DO PFT = 1,5
            WRITE( 6, 100 ) A_NET_st(PFT), RESP_st (PFT), CO2_IN_st(PFT)
         END DO 
 100     FORMAT( ES15.4E4, 3X, ES15.4E4, 3X, ES15.4E4 )
      END IF

      END SUBROUTINE INIT_ECOPHY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CLEANUP\_ECOPHY
!
! !DESCRIPTION: Subroutine CLEANUP_ECOPHY deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_ECOPHY()
!
! !REVISION HISTORY:
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !=================================================================
      ! CLEANUP_ECOPHY begins here!
      !=================================================================

      ! Return to calling program
      END SUBROUTINE CLEANUP_ECOPHY
!EOC
      END MODULE ECOPHY_MOD
