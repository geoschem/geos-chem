  !==========================================================================
  ! Parameters
  !==========================================================================

  ! Universal gas consatant [bar/(mol/kg)/K]  (Source: NIST, 2014)
  ! NOTE: Make sure this is consistent w/ the value in physconsts.F90!
  REAL(dp), PARAMETER :: CON_R = 0.083144598_dp

  !==========================================================================
  ! Molecular weights and Henry's law constants
  !==========================================================================

  ! Species molecular weight [g/mole]
  REAL(dp) :: MW(NSPEC)

  ! Square root of species molecular weight [g/mole]
  REAL(dp) :: SR_MW(NSPEC)

  ! Henry's law constants
  REAL(dp) :: HENRY_K0(NSPEC)  ! Solubility constant [M/atm]
  REAL(dp) :: HENRY_CR(NSPEC)  ! Volatility constant [K    ]

  !==========================================================================
  ! Meterological variables
  !==========================================================================

  ! Inverse of temperature [1/K]
  REAL(dp) :: INV_TEMP
  !$OMP THREADPRIVATE( INV_TEMP )

  ! Square root of Temperature [K**0.5]
  REAL(dp) :: SR_TEMP
  !$OMP THREADPRIVATE( SR_TEMP )

  ! TEMP_OVER_K300 - Temperature [K] / 300.0
  REAL(dp) :: TEMP_OVER_K300
  !$OMP THREADPRIVATE( TEMP_OVER_K300 )

  ! K300_OVER_TEMP - 300.0 / Temperature [K]
  REAL(dp) :: K300_OVER_TEMP
  !$OMP THREADPRIVATE( K300_OVER_TEMP )

  ! NUMDEN - Air number density [#/cm3]
  REAL(dp) :: NUMDEN
  !$OMP THREADPRIVATE( NUMDEN )

  ! Pressure and relative humidity
  REAL(dp) :: PRESS
  !$OMP THREADPRIVATE( PRESS )

  REAL(dp) :: RELHUM
  !$OMP THREADPRIVATE( RELHUM )

  ! Cosine of solar zenith angle
  REAL(dp) :: SUNCOS
  !$OMP THREADPRIVATE( SUNCOS )

  !=========================================================================
  ! Rate array variables (used by more than one mechanism)
  !=========================================================================

  ! Photolysis rates (increase size if necessary)
  ! Used by fullchem and Hg mechanisms
  REAL(dp) :: PHOTOL(200)
  !$OMP THREADPRIVATE( PHOTOL )


  !=========================================================================
  ! Variables specific to the "fullchem" mechanism
  !==========================================================================

  ! 4.0 * CON_R * TEMP
  REAL(dp) :: FOUR_R_T
  !$OMP THREADPRIVATE( FOUR_R_T )

  ! 4.0 * RGASLATM * TEMP
  REAL(dp) :: FOUR_RGASLATM_T
  !$OMP THREADPRIVATE( FOUR_RGASLATM_T )

  ! 8.0 * RSTARG * TEMP
  REAL(dp) :: EIGHT_RSTARG_T
  !$OMP THREADPRIVATE( EIGHT_RSTARG_T )

  ! H2O concentration
  REAL(dp) :: H2O
  !$OMP THREADPRIVATE( H2O )

  ! Array for aqueous sulfur chemistry rates (1/s)
  REAL(dp) :: K_MT(6)
  !$OMP THREADPRIVATE( K_MT )

  ! Array for sulfur chemistry rates in cloud (1/s)
  REAL(dp) :: K_CLD(6)
  !$OMP THREADPRIVATE( K_CLD )

  ! Liquid water conversion factor
  ! - Size = number of aerosol types (nAeroType)
  REAL(dp) :: CVFAC(14)
  !$OMP THREADPRIVATE( CVFAC )

  ! Proton activity [unitless] and H+ concentration [M]
  ! (assumed equivalent - for now):
  REAL(dp) :: H_PLUS
  !$OMP THREADPRIVATE( H_PLUS )

  TYPE, PUBLIC :: HetState
     REAL(dp) :: AVO            ! Avogadro's constant [molec/mol]
     LOGICAL  :: debugBox       ! Are we in a debugging box?
     LOGICAL  :: natSurface     ! Is there NAT in this box?
     LOGICAL  :: pscBox         ! Are there polar strat clouds?
     LOGICAL  :: stratBox       ! Are we in the stratosphere
     INTEGER  :: NAEROTYPE      ! Number of aerosol types
     REAL(dp) :: aClArea        ! Fine SSA+SNA aerosol area [cm2/cm3]
     REAL(dp) :: aClRadi        ! Fine SSA+SNA aerosol radius [cm]
     REAL(dp) :: aClVol         ! Fine SSA+SNA aerosol volume [cm3/cm3]
     REAL(dp) :: aIce           ! Ice surface area [cm2/cm3]
     REAL(dp) :: aLiq           ! Liquid surface area [cm2/cm3]
     REAL(dp) :: aWater(2)      ! ISORROPIA aerosol water (fine & coarse)
     REAL(dp) :: Br_conc_CldA   ! Br- in cloud (A=fine SSA) [mol/kg H2O]
     REAL(dp) :: Br_conc_CldC   ! Br- in cloud (C=coarse SSA) [mol/kg H2O]
     REAL(dp) :: Br_conc_CldG   ! Br- in cloud (G=gas-phase) [mol/kg H2O]
     REAL(dp) :: Br_conc_Cld    ! Br- total in cloud = A + C + G
     REAL(dp) :: Br_conc_SSA    ! Br- in fine SSA [mol/kg H2O]
     REAL(dp) :: Br_conc_SSC    ! Br- in coarse SSA  [mol/kg H2O]
     REAL(dp) :: Br_over_Cl_Cld ! Br/Cl ratio in cloud
     REAL(dp) :: Br_over_Cl_SSA ! Br/Cl ratio in fine sea salt
     REAL(dp) :: Br_over_Cl_SSC ! Br/Cl ratio in coarse sea salt
     REAL(dp) :: Cl_conc_CldA   ! Cl- in cloud (A=fine SSA) [mol/kg H2O]
     REAL(dp) :: Cl_conc_CldC   ! Cl- in cloud (C=coarse SSA) [mol/kg H2O]
     REAL(dp) :: Cl_conc_CldG   ! Cl- in cloud (G=gas-phase) [mol/kg H2O]
     REAL(dp) :: Cl_conc_Cld    ! Cl- total = A + C + G
     REAL(dp) :: Cl_conc_SSA    ! Cl- in fine   sea salt [mol/kg H2O]
     REAL(dp) :: Cl_conc_SSC    ! Cl- in coarse sea salt [mol/kg H2O]
     REAL(dp) :: cldFr          ! Cloud fraction
     REAL(dp) :: clearFr        ! Clear sky fraction
     REAL(dp) :: frac_Br_CldA   ! Frac of Br- in fine SSA(in cloud)
     REAL(dp) :: frac_Br_CldC   ! Frac of Br- in coarse SSA (in cloud)
     REAL(dp) :: frac_Br_CldG   ! Frac of Br- in gas phase
     REAL(dp) :: frac_Cl_CldA   ! Frac of Cl- in fine SSA (in cloud)
     REAL(dp) :: frac_Cl_CldC   ! Frac of Cl- in coarse SSA (in cloud)
     REAL(dp) :: frac_Cl_CldG   ! Frac of Cl- in gas phase
     REAL(dp) :: frac_HSO3_aq   ! HSO3_aq / ( HSO3_aq + SO3_aq )
     REAL(dp) :: frac_SALACL    ! Frac of SALACL / total fine sea salt
     REAL(dp) :: frac_SO3_aq    ! SO3_aq  / ( HSO3_aq + SO3_aq )
     REAL(dp) :: fupdateHOBr    ! Correction factor - HOBr
     REAL(dp) :: fupdateHOCl    ! Correction factor - HOCl
     REAL(dp) :: gamma_HO2      ! Uptake probability for HO2  [1]
     REAL(dp) :: H2O            ! H2O concentration
     REAL(dp) :: HBr_theta      ! HBr theta for uptake on ice
     REAL(dp) :: HCl_theta      ! HCl theta for uptake on ice
     REAL(dp) :: H_conc_ICl     ! Ice phase pH, Cl-
     REAL(dp) :: H_conc_LCl     ! Liquid phase pH, Cl-
     REAL(dp) :: H_conc_SSA     ! Liquid phase pH, fine sea salt
     REAL(dp) :: H_conc_SSC     ! Liquid phase pH, coarse sea salt
     REAL(dp) :: H_conc_Sul     ! Liquid phase pH, sulfate
     REAL(dp) :: HNO3_theta     ! HNO3 theta for uptake on ice
     REAL(dp) :: HSO3_aq        ! Cloud bisulfite (aka HSO3-) [mol/L]
     REAL(dp) :: SO3_aq         ! Cloud sulfite (aka SO3--) [mol/L]
     REAL(dp) :: HSO3m          ! HSO3-/SO2 ratio [1]
     REAL(dp) :: SO3mm          ! SO3--/SO2 ratio [1]
     REAL(dp) :: TSO3_aq        ! Cloud sulfite + bisulfite [mol/L]
     REAL(dp) :: H_plus         ! Proton activity [1] and H+ conc [M]
     REAL(dp) :: KHETI_SLA(11)  ! Probs for PSC uptk rxns on SLA [1]
     REAL(dp) :: HSO4_molal     ! Bisulfate concentration [M]
     REAL(dp) :: NO3_molal      ! Nitrate concentration [M]
     REAL(dp) :: SO4_molal      ! Sulfate concentration [M]
     REAL(dp) :: NIT_conc_SSA   ! Cl- in fine sea salt [mol/kg H2O]
     REAL(dp) :: NIT_conc_SSC   ! Cl- in coarse sea salt [mol/kg H2O]
     REAL(dp) :: PI             ! PI constant
     REAL(dp) :: pHCloud        ! Cloud PH
     REAL(dp) :: pHSSA(2)       ! Sea salt pH (1=fine, 2=coarse)
     REAL(dp) :: OMOC_POA       ! Org matter/orgc carbon in POA [1]
     REAL(dp) :: OMOC_OPOA      ! Org matter/org carbon in POA [1]
     REAL(dp) :: qIce           ! Ice mixing ratio [kg/kg]
     REAL(dp) :: qLIq           ! Water mixing ratio [kg/kg]
     REAL(dp) :: rIce           ! Ice radius
     REAL(dp) :: rLiq           ! Liquid radius
     REAL(dp) :: ssAlk(2)       ! Sea salt alk'nty (1=fine, 2=coarse)
     LOGICAL  :: SSA_is_Alk     ! Is fine sea-salt alkaline?
     LOGICAL  :: SSA_is_Acid    ! Is fine sea-salt acidic?
     LOGICAL  :: SSC_is_Alk     ! Is coarse sea-salt alkaline?
     LOGICAL  :: SSC_is_Acid    ! Is coarse sea-salt acidic?
     REAL(dp) :: vAir           ! Volume of air [cm3]
     REAL(dp) :: vIce           ! Ice volume [cm3]
     REAL(dp) :: vLiq           ! Liquid volume [cm3]
     REAL(dp) :: wetArea(14)    ! Aerosol specific wet sfc area [cm3/cm3]
     REAL(dp) :: xArea(14)      ! Aerosol specific sfc area  [cm3/cm3]
     REAL(dp) :: xH2O(14)       ! Aerosol water content [cm3/cm3]
     REAL(dp) :: xRadi(14)      ! Aerosol effective radius [cm]
     REAL(dp) :: xVol(14)       ! Aerosol specific volume [cm3/cm3]
     !
     ! Additional fields only used by the Hg simulation
     !
     LOGICAL  :: cloudBox       ! Are we in a box with cloud?
     REAL(dp) :: fracOrgAer     ! Fraction of organic aerosol [1]
     REAL(dp) :: fracInorgAer   ! Fraction of organic aerosol [1]
  END TYPE HetState
  TYPE(HetState), TARGET, PUBLIC :: State_Het
  !$OMP THREADPRIVATE( State_Het )

  !==========================================================================
  ! Variables specific to the carbon mechanism
  !==========================================================================

  ! Array for strat CH4/CO/CO2 rates
  REAL(dp) :: K_STRAT(5)
  !$OMP THREADPRIVATE( K_STRAT )

  ! Array for trop CH4/CO/CO2 rates
  REAL(dp) :: K_TROP(5)
  !$OMP THREADPRIVATE( K_TROP )

  REAL(dp) :: TROP
  !$OMP THREADPRIVATE( TROP )
