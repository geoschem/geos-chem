!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gckpp_HetRates
!
! !DESCRIPTION: FlexChem module for heterogeneous chemistry, via KPP.
!\\
!\\
! !INTERFACE:

MODULE GcKpp_HetRates
!
! !USES:
!
  USE CMN_FJX_MOD,      ONLY : NDUST
  USE CMN_FJX_MOD,      ONLY : NAER
  USE Error_Mod,        ONLY : ERROR_STOP
  USE Error_Mod,        ONLY : GEOS_CHEM_STOP
  USE Error_Mod,        ONLY : IS_SAFE_DIV, SAFE_DIV
  USE GcKpp_Global
  USE GcKpp_Precision
  USE GcKpp_Parameters
  USE State_Chm_Mod,    ONLY : ChmState
  USE State_Chm_Mod,    ONLY : Ind_
  USE State_Met_Mod,    ONLY : MetState
  USE Input_Opt_Mod,    ONLY : OptInput
  USE PhysConstants,    ONLY : AVO, RGASLATM, CONSVAP, RSTARG, PI
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Cld_Params
  PUBLIC :: Halide_Conc
  PUBLIC :: Set_Het
!
! !REMARKS:
!  There are NAEROTYPE=14 aerosol types used in het chem rate computations:
!  (1 ) Mineral dust (reff = 0.151 um)   (8 ) Tropospheric sulfate
!  (2 ) Mineral dust (reff = 0.253 um)   (9 ) Black Carbon
!  (3 ) Mineral dust (reff = 0.402 um)   (10) Organic Carbon
!  (4 ) Mineral dust (reff = 0.818 um)   (11) Fine (accum-mode) sea salt
!  (5 ) Mineral dust (reff = 1.491 um)   (12) Coarse sea salt
!  (6 ) Mineral dust (reff = 2.417 um)   (13) Stratospheric sulfate
!  (7 ) Mineral dust (reff = 3.721 um)   (14) Irregular ice cloud
!
! !REFERENCES:
!  Eastham et al., Development and evaluation of the unified tropospheric-
!    stratospheric chemistry extension (UCX) for the global chemistry-transport
!    model GEOS-Chem, Atmos. Env., doi:10.1016/j.atmosenv.2014.02.001, 2014.
!  Fisher et al, Organic nitrate chemistry and its implications for nitrogen
!    budgets in an isoprene- and monoterpene-rich atmosphere: constraints from
!    aircraft (SEAC4RS) and ground-based (SOAS) observations in the Southeast
!    US. Atmos. Chem. Phys., 16, 2961-2990, 2016.
!  Holmes, C.D., Bertram, T. H., Confer, K. L., Ronan, A. C., Wirks, C. K.,
!    Graham, K. A., Shah, V. (2019) The role of clouds in the tropospheric
!    NOx cycle: a new modeling approach for cloud chemistry and its global
!    implications, Geophys. Res. Lett. 46, 4980-4990,
!    https://doi.org/10.1029/2019GL081990
!  Marais et al., Aqueous-phase mechanism for secondary organic aerosol
!    formation from isoprene: application to the southeast United States and
!    co-benefit of SO2 emission controls, Atmos. Chem. Phys., 16, 1603-1618,
!    doi:10.5194/acp-16-1603-2016, 2016.
!  Parrella et al, Tropospheric bromine chemistry: implications for present and
!    pre-industrial ozone and mercury, Atmos. Chem. Phys., 12, 6,723-6,740,
!    doi:10.5194/acp-12-6723-2012, 2012.
!  Schmidt, J., et al., “Modelling the observed tropospheric BrO background:
!    Importance of multiphase chemistry & implications for ozone, OH, &
!    mercury”, J Geophys. Res-Atmos., 121, 024229,
!    https://doi.org/10.1002/2015JD024229, 2016.
!  Sherwen, T., et al., Global impacts of tropospheric halogens (Cl, Br, I) on
!    oxidants and composition in GEOS-Chem, Atmos. Chem. Phys., 16, 12239-12271,
!    https://doi.org/10.5194/acp-16-12239-2016, 2016.
!
! !REVISION HISTORY:
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
! !IROUTINE: Set_Het
!
! !DESCRIPTION: Main heterogenous chemistry driver routine.  Sets up the
!  vector of heterogeneous chemistry rates for the KPP chemistry solver.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_HET( I, J, L, Input_Opt, State_Chm, State_Met )
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L    ! Lon, lat, level indices
    TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
!EOP
!-----------------------------------------------------------------------------
!BOC

!LOCAL VARIABLES:

    ! SAVEd scalars
    INTEGER,           SAVE :: id_SALA = -1
    LOGICAL,           SAVE :: FIRST   = .TRUE.

    ! Scalars
    LOGICAL                 :: SAFEDIV
    INTEGER                 :: IND
    REAL(dp)                :: kITemp, kIITemp, rate, gamma

    ! Pointers and objects
    TYPE(HetState), POINTER :: H
!
! !DEFINED PARAMETERS:
!
    LOGICAL, PARAMETER      :: YES = .TRUE.
    LOGICAL, PARAMETER      :: NO  = .FALSE.

    !====================================================================
    ! SET_HET begins here!
    !====================================================================

    ! Initialize
    H         => State_Het  ! Point to State_Het in gckpp_Global.F90

    !--------------------------------------------------------------------
    ! Calculate and pass het rates to the KPP rate array
    !--------------------------------------------------------------------

    ! Zero the HET array
    HET = 0.0_dp

!    !========================================================================
!    ! Br/Cl heterogeneous chemistry
!    !========================================================================
!
!    !------------------------------------------------------------------------
!    ! BrNO3 hydrolysis (update: XW 2019-06-08)
!    !------------------------------------------------------------------------
!
!    HET(ind_BrNO3, 1) = BrNO3_Hydrolysis( SR_MW(ind_NO3), H )
!    ! Reaction probabilities on liquid and ice
!    gammaLiq    = 0.0021_dp * TEMP - 0.561_dp           ! Deiber et al, 2004
!    gammaIce    = 5.3e-4_dp * EXP( 1100.0_dp / TEMP )
!
!    ! BrNO3 + H2O hydrolysis in clear sky
!    rate        = HetBrNO3( SR_MW(ind_BrNO3), gammaLiq, rate )
!
!    ! Add BrNO3 + H2O hydrolysis in tropospheric cloud
!    rate = rate + CloudHet( CldFr,  Aliq,             Aice,    rLiq,         &
!                            rIce,   SR_MW(ind_BrNO3), gammaLiq, gammaIce,    &
!                            1.0_dp, 1.0_dp                                  )
!
!    ! BrNO3 + H2O overall rate
!    HET(ind_BrNO3, 1) = kIIR1Ltd( C(ind_BrNO3), C(ind_H2O), rate            )
!
!    !-----------------------------------------------------------------------
!    ! ClNO3 hydrolysis (update: XW 2019-06-08)
!    !-----------------------------------------------------------------------
!
!    ! ClNO3 + H2O hydrolysis rate in clear sky
!    rate = HetClNO3( ClearFr, SR_MW(ind_ClNO3), clConc_SALA, brConc_SALA    )
!
!    ! Add ClNO3 + H2O hydrolysis in tropospheric cloud
!    rate = rate + HetClNO3_TCld( rLiq,        rIce,        ALiq,             &
!                                 AIce,        VAir,        CldFr,            &
!                                 clConc_CldA, clConc_CldC, clConc_Cldg,      &
!                                 brConc_CldA, brConc_CldC, brConc_Cldg,      &
!                                 hno3_th,     hcl_th,      hbr_th,      7   )
!
!    ! ClNO3 + H2O Overall rate
!    HET(ind_ClNO3, 1) = kIIR1Ltd( C(ind_ClNO3), C(ind_H2O), rate            )
!!!^^^^ Identical w/r/t ref up to here ^^^^

!
!    !------------------------------------------------------------------------
!    ! HOBr + HBr (update: XW 2019-06-08)
!    !------------------------------------------------------------------------
!    rate = HETHOBr_HBr( H%HOBr%MW_g, 0.0_dp, Input_Opt                      )
!
!    rate = rate + HETHOBr_TCld( rLiq,         rIce,        ALiq,             &
!                                AIce,         VAir,        CldFr,            &
!                                hConc_Sul,    hConc_LCl,   hConc_ICl,        &
!                                clConc_CldA,  clConc_CldC, clConc_Cldg,      &
!                                brConc_CldA,  brConc_CldC, brConc_Cldg,      &
!                                HSO3conc_Cld, SO3conc_Cld, hno3_th,          &
!                                hcl_th,       hbr_th,      6                )
!
!    HET(ind_HOBr,  1) = kIIR1Ltd( C(ind_HOBr), C(ind_HBr), rate             )
!
!    !------------------------------------------------------------------------
!    ! HOBr + HCl (update: XW 2019-06-08)
!    !------------------------------------------------------------------------
!    kITemp = HETHOBr_HCl( H%HOBr%MW_g, 0.0_dp, Input_Opt                    )
!
!    kITemp = kITemp                                                          &
!           + HETHOBr_TCld( rLiq,        rIce,                   &
!                           ALiq,        AIce,        VAir,                   &
!                           CldFr,       hConc_Sul,              &
!                           hConc_LCl,   hConc_ICl,   clConc_CldA,            &
!                           clConc_CldC, clConc_Cldg, brConc_CldA,            &
!                           brConc_CldC, brConc_Cldg, HSO3conc_Cld,           &
!                           SO3conc_Cld, hno3_th,     hcl_th,                 &
!                           hbr_th,      5,           H                      )
!
!    HET(ind_HOBr,  2) = kIIR1Ltd( C(ind_HOBr), C(ind_HCl), kITemp           )
!
!    !------------------------------------------------------------------------
!    ! HOBr + BrSalA/C (update: XW 2019-06-08)
!    !------------------------------------------------------------------------
!
!    ! First consider reaction in troposphere cloud
!    kITemp = HETHOBr_TCld( rLiq,        rIce,                   &
!                           ALiq,        AIce,        VAir,                   &
!                           CldFr,        hConc_Sul,              &
!                           hConc_LCl,   hConc_ICl,   clConc_CldA,            &
!                           clConc_CldC, clConc_Cldg, brConc_CldA,            &
!                           brConc_CldC, brConc_Cldg, HSO3conc_Cld,           &
!                           SO3conc_Cld, hno3_th,     hcl_th,                 &
!                           hbr_th,      3,           H                      )
!
!    ! Then calculate reaction on aerosols
!    kITemp = kITemp +                                                        &
!             HETHOBr_SS( NUMDEN,       AClRADI,     (1-CldFr)*AClAREA,        &
!                         SSAlk(1),    Temp,        hConc_SSA,                &
!                         clConc_SALA, brConc_SALA, 2,                        &
!                         H                                                  )
!
!    HET(ind_HOBr,  5) = kIIR1Ltd( C(ind_HOBr), C(ind_BrSALA), kITemp        )
!
!    ! First consider reaction in troposphere cloud
!    kITemp = HETHOBr_TCld( rLiq,        rIce,                   &
!                           ALiq,        AIce,        VAir,                   &
!                           CldFr,       hConc_Sul,              &
!                           hConc_LCl,   hConc_ICl,   clConc_CldA,            &
!                           clConc_CldC, clConc_Cldg, brConc_CldA,            &
!                           brConc_CldC, brConc_Cldg, HSO3conc_Cld,           &
!                           SO3conc_Cld, hno3_th,     hcl_th,                 &
!                           hbr_th,      4,           H                      )
!
!    ! Then calculate reaction on aerosols
!    kITemp = kITemp                                                          &
!           + HETHOBr_SS( NUMDEN,       xRadi(12),   (1-CldFr)*xArea(12),      &
!                         SSAlk(2),    Temp,        hConc_SSC,                &
!                         clConc_SALC, brConc_SALC, 2,                        &
!                         H                                                  )
!
!    HET(ind_HOBr,  6) = kIIR1Ltd( C(ind_HOBr),  C(ind_BrSALC), kITemp       )
!
!    !------------------------------------------------------------------------
!    ! HOBr + Cl-(p) (update: XW 2019-06-08)
!    !------------------------------------------------------------------------
!    ! First consider reaction in troposphere cloud
!    kITemp = HETHOBr_TCld( rLiq,        rIce,                   &
!                           ALiq,        AIce,        VAir,                   &
!                           CldFr,       hConc_Sul,              &
!                           hConc_LCl,   hConc_ICl,   clConc_CldA,            &
!                           clConc_CldC, clConc_Cldg, brConc_CldA,            &
!                           brConc_CldC, brConc_Cldg, HSO3conc_Cld,           &
!                           SO3conc_Cld, hno3_th,     hcl_th,                 &
!                           hbr_th,      1,           H                      )
!
!    ! Then calculate reaction on aerosols
!    kITemp = kITemp                                                          &
!           + HETHOBr_SS( NUMDEN, AClRADI,     (1-CldFr)*AClAREA,        &
!                         SSAlk(1),    Temp,        hConc_SSA,                &
!                         clConc_SALA, brConc_SALA, 1,                        &
!                         H                                                  )
!
!    HET(ind_HOBr,  3) = kIIR1Ltd( C(ind_HOBr), C(ind_SALACL), kITemp        )
!
!    ! First consider reaction in troposphere cloud
!    kITemp = HETHOBr_TCld( rLiq,        rIce,                   &
!                           ALiq,        AIce,        VAir,                   &
!                           CldFr,       hConc_Sul,              &
!                           hConc_LCl,   hConc_ICl,   clConc_CldA,            &
!                           clConc_CldC, clConc_Cldg, brConc_CldA,            &
!                           brConc_CldC, brConc_Cldg, HSO3conc_Cld,           &
!                           SO3conc_Cld, hno3_th,     hcl_th,                 &
!                           hbr_th,      2,           H                      )
!
!    ! Then calculate reaction on aerosols
!    kITemp = kITemp                                                          &
!           + HETHOBr_SS( xRadi(12), ClearFr*xArea(12), SSAlk(2),             &
!                         hConc_SSC, clConc_SALC,       brConc_SALC, 1       )
!
!    HET(ind_HOBr,  4) = kIIR1Ltd( C(ind_HOBr), C(ind_SALCCL), kITemp        )
!
!    !------------------------------------------------------------------------
!    ! HOBr + HSO3-(aq) (update: XW 2019-06-08)
!    !------------------------------------------------------------------------
!
!    ! This reaction is first order, so no kII calculation is required
!    kITemp = HETHOBr_TCld( rLiq,        rIce,                   &
!                           ALiq,        AIce,        VAir,                   &
!                           CldFr,       hConc_Sul,              &
!                           hConc_LCl,   hConc_ICl,   clConc_CldA,            &
!                           clConc_CldC, clConc_Cldg, brConc_CldA,            &
!                           brConc_CldC, brConc_Cldg, HSO3conc_Cld,           &
!                           SO3conc_Cld, hno3_th,     hcl_th,                 &
!                           hbr_th,      7,           H                      )
!
!    ! Make sure sulfate produced is less than SO2 available (qjc, 06/20/16)
!    HET(ind_HOBr,  7) = kITemp * fupdateHOBr
!
!    !------------------------------------------------------------------------
!    ! HOBr + SO3--(aq) (update: XW 2019-06-08)
!    !------------------------------------------------------------------------
!    ! This reaction is first order, so no kII calculation is required
!    kITemp = HETHOBr_TCld( rLiq,        rIce,                   &
!                           ALiq,        AIce,        VAir,                   &
!                           CldFr,       hConc_Sul,              &
!                           hConc_LCl,   hConc_ICl,   clConc_CldA,            &
!                           clConc_CldC, clConc_Cldg, brConc_CldA,            &
!                           brConc_CldC, brConc_Cldg, HSO3conc_Cld,           &
!                           SO3conc_Cld, hno3_th,     hcl_th,                 &
!                           hbr_th,      8,           H                      )
!
!    ! Make sure sulfate produced is less than SO2 available (qjc, 06/20/16)
!    HET(ind_HOBr,  8) = kITemp * fupdateHOBr
!
!    !-----------------------------------------------------------------------
!    ! ClNO3 + BrSALA/C (update: XW 2019-06-08)
!    !-----------------------------------------------------------------------
!    ! First consider reaction in troposphere cloud
!    kITemp = HETClNO3_TCld( rLiq,        rIce,                 &
!                            ALiq,        AIce,        VAir,                  &
!                            CldFr,       clConc_CldA,           &
!                            clConc_CldC, clConc_Cldg, brConc_CldA,           &
!                            brConc_CldC, brConc_Cldg, hno3_th,               &
!                            hcl_th,      hbr_th,      3,                     &
!                            H                                               )
!
!    ! Then calculate reaction on aerosols
!    kITemp = kITemp                                                          &
!           + HETClNO3_SS( AClRADI,     ClearFr*AClAREA,  SSAlk(1),           &
!                          clConc_SALA, brConc_SALA,      2,        1        )
!
!    HET(ind_ClNO3, 4) = kIIR1Ltd( C(ind_ClNO3), C(ind_BrSALA), kITemp       )
!
!    ! First consider reaction in troposphere cloud
!    kITemp = HETClNO3_TCld( rLiq,        rIce,                  &
!                            ALiq,        AIce,        VAir,                  &
!                            CldFr,       clConc_CldA,           &
!                            clConc_CldC, clConc_Cldg, brConc_CldA,           &
!                            brConc_CldC, brConc_Cldg, hno3_th,               &
!                            hcl_th,      hbr_th,      4,                     &
!                            H                                               )
!
!    ! Then calculate reaction on aerosols
!    kITemp = kITemp                                                          &
!           + HETClNO3_SS( xRadi(12),   ClearFr*xArea(12), SSAlk(2),          &
!                          clConc_SALC, brConc_SALC,       2,         2      )
!
!    HET(ind_ClNO3, 5) = kIIR1Ltd( C(ind_ClNO3), C(ind_BrSALC), kITemp       )
!
!    !-----------------------------------------------------------------------
!    ! ClNO3 + Cl-(p) (update: XW 2019-06-08)
!    !-----------------------------------------------------------------------
!    ! First consider reaction in troposphere cloud
!    kITemp = HETClNO3_TCld( rLiq,        rIce,                  &
!                            ALiq,        AIce,        VAir,                  &
!                            CldFr,       clConc_CldA,           &
!                            clConc_CldC, clConc_Cldg, brConc_CldA,           &
!                            brConc_CldC, brConc_Cldg, hno3_th,               &
!                            hcl_th,      hbr_th,      1,                     &
!                            H                                               )
!
!    ! Then calculate reaction on aerosols
!    kITemp = kITemp                                                          &
!           + HETClNO3_SS( AClRADI,     ClearFr*AClAREA, SSAlk(1),            &
!                          clConc_SALA, brConc_SALA,     1,        1         )
!
!    HET(ind_ClNO3, 6) = kIIR1Ltd( C(ind_ClNO3), C(ind_SALACL), kITemp       )
!
!    ! First consider reaction in troposphere cloud
!    kITemp = HETClNO3_TCld( rLiq,        rIce,                  &
!                            ALiq,        AIce,        VAir,                  &
!                            CldFr,       clConc_CldA,           &
!                            clConc_CldC, clConc_Cldg, brConc_CldA,           &
!                            brConc_CldC, brConc_Cldg, hno3_th,               &
!                            hcl_th,      hbr_th,      2,                     &
!                            H                                               )
!
!    ! Then calculate reaction on aerosols out of cloud
!    kITemp = kITemp                                                          &
!           + HETClNO3_SS( xRadi(12),   ClearFr*xArea(12), SSAlk(2),          &
!                          clConc_SALC, brConc_SALC,       1,        2       )
!
!    HET(ind_ClNO3, 7) = kIIR1Ltd( C(ind_ClNO3), C(ind_SALCCL), kITemp       )
!
!    !----------------------------------------------------------------------
!    ! ClNO3 + HCl (update: XW 2019-06-08)
!    !----------------------------------------------------------------------
!    kITemp = HETClNO3_HCl( H%ClNO3%MW_g, 0.0_dp                             )
!
!    kITemp = kITemp                                                          &
!           + HETClNO3_TCld( rLiq,        rIce,                  &
!                            ALiq,        AIce,        VAir,                  &
!                            CldFr,       clConc_CldA,           &
!                            clConc_CldC, clConc_Cldg, brConc_CldA,           &
!                            brConc_CldC, brConc_Cldg, hno3_th,               &
!                            hcl_th,      hbr_th,      5,                     &
!                            H                                               )
!
!    HET(ind_ClNO3, 2) = kIIR1Ltd( C(ind_ClNO3), C(ind_HCl), kITemp          )
!
!    !----------------------------------------------------------------------
!    ! ClNO3 + HBr (update: XW 2019-06-08)
!    !----------------------------------------------------------------------
!    kITemp = HETClNO3_HBr( H%ClNO3%MW_g, 0.0_dp, Input_Opt                  )
!
!    kITemp = kITemp                                                          &
!           + HETClNO3_TCld( rLiq,        rIce,                  &
!                            ALiq,        AIce,        VAir,                  &
!                            CldFr,       clConc_CldA,           &
!                            clConc_CldC, clConc_Cldg, brConc_CldA,           &
!                            brConc_CldC, brConc_Cldg, hno3_th,               &
!                            hcl_th,      hbr_th,      6,                     &
!                            H                                               )
!
!    HET(ind_ClNO3, 3) = kIIR1Ltd( C(ind_ClNO3), C(ind_HBr), kITemp          )
!
!    !-----------------------------------------------------------------------
!    ! HOCl + HCl and HOCl + HBr (update: XW 2019-06-08)
!    !-----------------------------------------------------------------------
!    kITemp = HETHOCl_HCl( H%HOCl%MW_g, 0.0_dp, Input_opt                    )
!    kITemp = kITemp                                                          &
!           + HETHOCl_TCld( NUMDEN,       rLiq,        rIce,                   &
!                           ALiq,        AIce,        VAir,                   &
!                           Temp,        CldFr,       hConc_Sul,              &
!                           hConc_LCl,   hConc_ICl,   clConc_CldA,            &
!                           clConc_CldC, clConc_Cldg, HSO3conc_Cld,           &
!                           SO3conc_Cld, hcl_th,      3,                      &
!                           H                                                )
!
!    HET(ind_HOCl,  1) = kIIR1Ltd( C(ind_HOCl), C(ind_HCl), kITemp           )
!
!    kiTemp            = HETHOCl_HBr( H%HOCl%MW_g, 0E+0_dp, Input_Opt        )
!
!    HET(ind_HOCl,  2) = kIIR1Ltd( C(ind_HOCl), C(ind_HBr), kiTemp           )
!
!    !------------------------------------------------------------------------
!    ! HOCl + Cl-(p) (update: XW 2019-06-08)
!    !------------------------------------------------------------------------
!    ! First consider reaction in troposphere cloud
!    kITemp = HETHOCl_TCld( NUMDEN,       rLiq,        rIce,                   &
!                           ALiq,        AIce,        VAir,                   &
!                           Temp,        CldFr,       hConc_Sul,              &
!                           hConc_LCl,   hConc_ICl,   clConc_CldA,            &
!                           clConc_CldC, clConc_Cldg, HSO3conc_Cld,           &
!                           SO3conc_Cld, hcl_th,      1,                      &
!                           H                                                )
!
!    ! Then calculate reaction on aerosols
!    kITemp = kITemp                                                          &
!           + HETHOCl_SS( NUMDEN,       AClRADI, (1-CldFr)*AClAREA,            &
!                         SSAlk(1),    Temp,   hConc_SSA,                     &
!                         clConc_SALA, H                                     )
!
!    HET(ind_HOCl,  3) = kIIR1Ltd( C(ind_HOCl), C(ind_SALACL), kITemp        )
!
!    ! First consider reaction in troposphere cloud
!    kITemp = HETHOCl_TCld( NUMDEN,       rLiq,        rIce,                   &
!                           ALiq,        AIce,        VAir,                   &
!                           Temp,        CldFr,       hConc_Sul,              &
!                           hConc_LCl,   hConc_ICl,   clConc_CldA,            &
!                           clConc_CldC, clConc_Cldg, HSO3conc_Cld,           &
!                           SO3conc_Cld, hcl_th,      2,                      &
!                           H                                                )
!
!    ! Then calculate reaction on aerosols
!    kITemp = kITemp                                                          &
!           + HETHOCl_SS( NUMDEN,       xRadi(12), (1-CldFr)*xArea(12),        &
!                         SSAlk(2),    Temp,      hConc_SSC,                  &
!                         clConc_SALC, H                                     )
!
!    HET(ind_HOCl,  4) = kIIR1Ltd( C(ind_HOCl), C(ind_SALCCL), kITemp        )
!
!    !------------------------------------------------------------------------
!    ! HOCl + HSO3-/SO3--(aq)
!    !------------------------------------------------------------------------
!    ! This reaction is first order, so no kII calculation is required
!    kITemp = HETHOCl_TCld( NUMDEN,       rLiq,        rIce,                   &
!                           ALiq,        AIce,        VAir,                   &
!                           Temp,        CldFr,       hConc_Sul,              &
!                           hConc_LCl,   hConc_ICl,   clConc_CldA,            &
!                           clConc_CldC, clConc_Cldg, HSO3conc_Cld,           &
!                           SO3conc_Cld, hcl_th,      4,                      &
!                           H                                                )
!
!    ! Make sure sulfate produced is less than SO2 available
!    HET(ind_HOCl,  5) = kITemp * fupdateHOCl
!
!    ! This reaction is first order, so no kII calculation is required
!    kITemp = HETHOCl_TCld( NUMDEN,       rLiq,        rIce,                   &
!                           ALiq,        AIce,        VAir,                   &
!                           Temp,        CldFr,       hConc_Sul,              &
!                           hConc_LCl,   hConc_ICl,   clConc_CldA,            &
!                           clConc_CldC, clConc_Cldg, HSO3conc_Cld,           &
!                           SO3conc_Cld, hcl_th,      5,                      &
!                           H                                                )
!
!    ! Make sure sulfate produced is less than SO2 available
!    HET(ind_HOCl,  6) = kITemp * fupdateHOCl
!
!    !------------------------------------------------------------------------
!    ! O3 + Br- calculation (TS index: hhc12)
!    !------------------------------------------------------------------------
!    kITemp = HETO3_TCld( NUMDEN,       rLiq,        rIce,                     &
!                         ALiq,        AIce,        VAir,                     &
!                         Temp,        CldFr,       brConc_CldA,              &
!                         brConc_CldC, brConc_Cldg, C(ind_O3),                &
!                         3,           H                                     )
!
!    HET(ind_O3,    1) = kIIR1Ltd( C(ind_O3), C(ind_HBr), kITemp             )
!
!    !------------------------------------------------------------------------
!    ! O3 + BrSALA/C calculations (TMS index: hhc13/14)
!    !------------------------------------------------------------------------
!    kITemp = HETO3_TCld( NUMDEN,       rLiq,        rIce,                     &
!                         ALiq,        AIce,        VAir,                     &
!                         Temp,        CldFr,       brConc_CldA,              &
!                         brConc_CldC, brConc_Cldg, C(ind_O3),                &
!                         1,           H                                     )
!
!    kITemp = kITemp                                                          &
!           + HETO3_SS( NUMDEN,         AClRADI, (1-CldFr)*AClAREA,            &
!                       SSAlk(1),      Temp,    brConc_SALA,                  &
!                       C(ind_O3),     H                                     )
!
!    HET(ind_O3,    2) = kIIR1Ltd( C(ind_O3), C(ind_BrSALA), kITemp          )
!
!    kITemp  = HETO3_TCld( NUMDEN,       rLiq,        rIce,                    &
!                          ALiq,        AIce,        VAir,                    &
!                          Temp,        CldFr,       brConc_CldA,             &
!                          brConc_CldC, brConc_Cldg, C(ind_O3),               &
!                          2,           H                                    )
!
!    kITemp = kITemp                                                          &
!           + HETO3_SS( NUMDEN,          xRadi(12), (1-CldFr)*xArea(12),       &
!                       SSAlk(2),       Temp,      brConc_SALC,               &
!                       C(ind_O3),      H                                    )
!
!    HET(ind_O3,    3) = kIIR1Ltd( C(ind_O3), C(ind_BrSALC), kITemp          )
!
!    !------------------------------------------------------------------------
!    ! Br uptake calculation - forms BrSALA/C (TS index: hhc17/18)
!    !------------------------------------------------------------------------
!    ! First-order reactions, no calculation of kII required
!    kITemp            = HETHXUptake( NUMDEN, AClRADI, (1-CldFr)*ACLAREA,      &
!                                     Temp, 2                                )
!    HET(ind_HBr,   1) = kITemp
!
!    kITemp            = HETHXUptake( NUMDEN, xRadi(12), (1-CldFr)*xArea(12),  &
!                                     Temp, 2                                )
!    HET(ind_HBr,   2) = kITemp
!
!    !------------------------------------------------------------------------
!    ! BrNO3 + HCl in stratosphere
!    !------------------------------------------------------------------------
!    ! NOTE: the restriction of these reactions to the troposphere has been
!    ! restored - TMS (2017/04/06 )
!    HET(ind_BrNO3, 2) = kIIR1Ltd( C(ind_BrNO3), C(ind_HCl),                  &
!                                  HETBrNO3_HCl( 1.42E2_dp, 0.0_dp )         )
!
!    !------------------------------------------------------------------------
!    ! N2O5 + HCl in stratosphere
!    !------------------------------------------------------------------------
!    ! NOTE: this extension of calculation in troposphere has been removed
!    !  (TMS 17/04/10)
!    kITemp            = HETN2O5_HCl( 1.08E2_dp, 0.0_dp, Input_Opt           )
!    HET(ind_N2O5,  2) = kIIR1Ltd( C(ind_N2O5), C(ind_HCl), kITemp           )
!
!    !------------------------------------------------------------------------
!    ! Reaction of N2O5 with Cl-, XW 1/24/18)
!    !------------------------------------------------------------------------
!    HetTemp(1:3)      = HETN2O5_SS( CldFr, 1 )
!    HET(ind_N2O5,  4) = kIIR1Ltd( C(ind_N2O5), C(ind_SALACL), HetTemp(1)    )
!    State_Chm%GammaN2O5(I,J,L,2) = HetTemp(1)
!
!    HetTemp(1:3)      = HETN2O5_SS( CldFr, 2)
!    HET(ind_N2O5,  5) = kIIR1Ltd( C(ind_N2O5), C(ind_SALCCL), HetTemp(1)    )
!    State_Chm%GammaN2O5(I,J,L,3) = HetTemp(1)
!
!    !------------------------------------------------------------------------
!    ! Reaction of ClNO2 with Cl-, XW 3/12/18)
!    !------------------------------------------------------------------------
!
!    ! First consider reaction in troposphere cloud
!    kITemp = HETClNO2_TCld( NUMDEN,       rLiq,        rIce,                  &
!                            ALiq,        AIce,        VAir,                  &
!                            Temp,        CldFr,       pHCloud,               &
!                            clConc_CldA, clConc_CldC, clConc_Cldg,           &
!                            brConc_CldA, brConc_CldC, brConc_Cldg,           &
!                            1,           H                                  )
!
!    ! Then calculate reaction on aerosols out of cloud
!    kITemp = kITemp                                                          &
!           + HETClNO2( NUMDEN,       AClRADI,     (1-CldFr)*AClAREA,          &
!                       SSAlk(1),    Temp,       pHSSA(1),                    &
!                       clConc_SALA, brConc_SALA, 1,                          &
!                        H                                                   )
!
!    HET(ind_ClNO2, 1) = kIIR1Ltd( C(ind_ClNO2), C(ind_SALACL), kITemp       )
!
!    ! First consider reaction in troposphere cloud
!    kITemp = HETClNO2_TCld( NUMDEN,       rLiq,        rIce,                  &
!                            ALiq,        AIce,        VAir,                  &
!                            Temp,        CldFr,       pHCloud,               &
!                            clConc_CldA, clConc_CldC, clConc_Cldg,           &
!                            brConc_CldA, brConc_CldC, brConc_Cldg,           &
!                            2,           H                                  )
!    ! not on SALC
!    HET(ind_ClNO2, 2) = kIIR1Ltd( C(ind_ClNO2), C(ind_SALCCL), kITemp       )
!
!    !------------------------------------------------------------------------
!    ! ClNO2 + dissolved HCl in cloud
!    !-----------------------------------------------------------------------
!    kITemp = HETClNO2_TCld( NUMDEN,       rLiq,        rIce,                  &
!                            ALiq,        AIce,        VAir,                  &
!                            Temp,        CldFr,       pHCloud,               &
!                            clConc_CldA, clConc_CldC, clConc_Cldg,           &
!                            brConc_CldA, brConc_CldC, brConc_Cldg,           &
!                            3,           H                                  )
!
!    HET(ind_ClNO2, 3) = kIIR1Ltd( C(ind_ClNO2), C(ind_HCl), kITemp          )
!
!    !------------------------------------------------------------------------
!    ! Reaction of ClNO2 with sea-salt Br-, XW 8/8/18)
!    !------------------------------------------------------------------------
!
!    ! First consider reaction in troposphere cloud
!    kITemp = HETClNO2_TCld( NUMDEN,       rLiq,        rIce,                  &
!                            ALiq,        AIce,        VAir,                  &
!                            Temp,        CldFr,       pHCloud,               &
!                            clConc_CldA, clConc_CldC, clConc_Cldg,           &
!                            brConc_CldA, brConc_CldC, brConc_Cldg,           &
!                            4,           H                                  )
!
!    ! Then calculate reaction on aerosols
!    kITemp = kITemp +                                                        &
!             HETClNO2( NUMDEN,       AClRADI,     (1-CldFr)*AClAREA,          &
!                       SSAlk(1),    Temp,        pHSSA(1),                   &
!                       clConc_SALA, brConc_SALA, 2,                          &
!                       H                                                    )
!
!    HET(ind_ClNO2, 4) = kIIR1Ltd( C(ind_ClNO2), C(ind_BrSALA), kITemp       )
!
!    ! First consider reaction in troposphere cloud
!    kITemp = HETClNO2_TCld( NUMDEN,       rLiq,        rIce,                  &
!                            ALiq,        AIce,        VAir,                  &
!                            Temp,        CldFr,       pHCloud,               &
!                            clConc_CldA, clConc_CldC, clConc_Cldg,           &
!                            brConc_CldA, brConc_CldC, brConc_Cldg,           &
!                            5,           H                                  )
!
!    ! Then calculate reaction on aerosols
!    kITemp = kITemp +                                                        &
!             HETClNO2( NUMDEN,       xRadi(12),   (1-CldFr)*xArea(12),        &
!                       SSAlk(2),    Temp,        pHSSA(2),                   &
!                       clConc_SALC, brConc_SALC, 2,                          &
!                       H                                                    )
!
!    HET(ind_ClNO2, 5) = kIIR1Ltd( C(ind_ClNO2), C(ind_BrSALC), kITemp       )
!
!    !------------------------------------------------------------------------
!    ! ClNO2 + dissolved HBr in cloud
!    !------------------------------------------------------------------------
!    kITemp = HETClNO2_TCld( NUMDEN,       rLiq,        rIce,                  &
!                            ALiq,        AIce,        VAir,                  &
!                            Temp,        CldFr,       pHCloud,               &
!                            clConc_CldA, clConc_CldC, clConc_Cldg,           &
!                            brConc_CldA, brConc_CldC, brConc_Cldg,           &
!                            6,           H                                  )
!
!    HET(ind_ClNO2, 6) = kIIR1Ltd( C(ind_ClNO2), C(ind_HBr), kITemp          )

    !========================================================================
    ! Cleanup & quit
    !========================================================================
!### KPP DEBUG
9999 continue
!#include "print_het.H"

    H => NULL()

  END SUBROUTINE SET_HET
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kiir1r2l
!
! !DESCRIPTION: Determine removal rates for both species in an uptake reaction
! without assuming which reactant is limiting.
!\\
!\\
! !INTERFACE:
!
  FUNCTION kIIR1R2L( spcVec, indGasA, indGasB, kIASource, kIBSource ) &
       RESULT( kII )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: spcVec(:)
    INTEGER,  INTENT(IN) :: indGasA
    INTEGER,  INTENT(IN) :: indGasB
    REAL(dp), INTENT(IN) :: kIASource
    REAL(dp), INTENT(IN) :: kIBSource
!
! !RETURN VALUE:
!
    REAL(dp)             :: kII
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(dp) :: concGasA, concGasB, kIA, kIB
    REAL(dp) :: R_GasA, R_GasB
    LOGICAL  :: nonZeroRate

    ! Get the base concentrations
    concGasA = spcVec(indGasA)
    concGasB = spcVec(indGasB)

    ! Copy the first estimates of kI for each species
    kIA = kIASource
    kIB = kIBSource

    ! Assume for now that the reaction will not proceed
    nonZeroRate = .False.
    kII = 0.0_dp

    ! Prevent reaction if either concentration is too low
    IF ((concGasA.gt.100.0e+0_dp).and.(concGasB.gt.100.0e+0_dp).and.&
        (kIA.gt.0.0e+0_dp).and.(kIB.gt.0.0e+0_dp)) THEN
       ! Calculate the overall rate based on each reactant
       R_GasA = kIA*concGasA
       R_GasB = kIB*concGasB
       IF (R_GasA > R_GasB) THEN
          ! Limited by uptake of B
          nonZeroRate = Is_Safe_Div( R_GasB, concGasA )
          IF (nonZeroRate) THEN
             kII = kIB / concGasA
          ENDIF
       ELSE
          ! Limited by uptake of A
          nonZeroRate = Is_Safe_Div( R_GasA, concGasB )
          IF (nonZeroRate) THEN
             kII = kIA / concGasB
          ENDIF
       ENDIF
    ENDIF

    ! If no tests were passed, zero out both rates
    IF (.not.nonZeroRate) THEN
       kII = 0.0_dp
    ENDIF

  END FUNCTION kIIR1R2L
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHXUptake
!
! !DESCRIPTION: Sets the uptake rate of HCl and HBr on sea salt using Johan
!  Schmidt's updated code.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHXUptake( denAir, rAer, AAer, TK, X ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(dp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(dp), INTENT(IN) :: rAer        ! Radius of aerosol (cm)
      REAL(dp), INTENT(IN) :: AAer        ! Area of aerosol (cm2/cm3)
      REAL(dp), INTENT(IN) :: TK          ! Temperature (K)
      INTEGER,  INTENT(IN) :: X           ! 1: Cl-, 2: Br-
!
! !RETURN VALUE:
!
      REAL(dp)             :: kISum
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(dp) :: XSTKCF, ADJUSTEDRATE, XSqM
      Real(dp), Parameter :: XMolWeightHCl=36.45e+0_dp
      Real(dp), Parameter :: XSqMHCl=SQRT(XMolWeightHCl)
      Real(dp), Parameter :: XMolWeightHBr=80.91e+0_dp
      Real(dp), Parameter :: XSqMHBr=SQRT(XMolWeightHBr)

!      ! Initialize
!      kISum        = 0.0_dp
!
!      ! Select between halogens
!      IF (X.eq.1) THEN
!         ! This would never be used since HCl uptake is
!         ! handled by ISSOROPIA now, xnw 1/25/18
!         XSqM = XSqMHCl
!      ELSEIF (X.eq.2) THEN
!         XSqM = XSqMHBr
!      ENDIF
!
!      XStkCf = Gamma_HX_Uptake( rAer, denAir, X, TK )
!
!      ! Reaction rate for surface of aerosol
!      kISum = Arsl1K(AAer,rAer,denAir,XStkCf,SR_TEMP,XSqM)

    END FUNCTION HETHXUptake
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: HetClNO2_TCld
!!
!! !DESCRIPTION: Sets the rate of the multiphase reaction ClNO2 + Cl-/Br- in
!!  troposphere cloud
!!\\
!!\\
!! !INTERFACE:
!!
    FUNCTION HETClNO2_TCld( denAir,   rLiq,     rIce,     ALiq,              &
                            AIce,     VAir,     TK,       CldFr,             &
                            pH,       clConc_A, clConc_C, clConc_g,          &
                            brConc_A, brConc_C, brConc_g, X,                 &
                            H                                                &
                          ) RESULT( rate )
!
! !INPUT PARAMETERS:
!
      REAL(dp),       INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(dp),       INTENT(IN) :: rLiq        ! Radius of liquid cloud
                                                !  droplets (cm)
      REAL(dp),       INTENT(IN) :: rIce        ! Radius of ice cloud
                                                !  crystals (cm)
      REAL(dp),       INTENT(IN) :: ALiq        ! Area of liquid cloud
                                                !  droplets (cm2/cm3)
      REAL(dp),       INTENT(IN) :: AIce        ! Area of ice cloud
                                                !  crystals (cm2/cm3)
      REAL(dp),       INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(dp),       INTENT(IN) :: TK          ! Temperature (K)
      REAL(dp),       INTENT(IN) :: CldFr       ! Cloud fraction
      REAL(dp),       INTENT(IN) :: clConc_A    ! Fine Chloride
                                                !  concentration (mol/L)
      REAL(dp),       INTENT(IN) :: clConc_C    ! Coarse Chloride
                                                !  concentration (mol/L)
      REAL(dp),       INTENT(IN) :: clConc_g
      REAL(dp),       INTENT(IN) :: brConc_A    ! Fine Bromide
                                                !  concentration (mol/L)
      REAL(dp),       INTENT(IN) :: brConc_C    ! Coarse Bromide
                                                !  concentration (mol/L)
      REAL(dp),       INTENT(IN) :: brConc_g
      REAL(dp),       INTENT(IN) :: pH
      INTEGER,        INTENT(IN) :: X           ! 1=fineCl;2=coarseCl;
                                                ! 3=HCl,4=fineBr;5=coarseBr;
                                                ! 6=HBr
      TYPE(HetState), POINTER    :: H           ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(dp)                   :: rate        ! Rxn rate [1/s]
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      REAL(dp) :: X1, X2, ADJUSTEDRATE
!      REAL(dp) :: GAM_ClNO2, r_gp, clConc, r_ac, brConc, r1
!      INTEGER  :: Y
!
!      ! Initialize the return value
!      rate = 0.0_dp
!
!      ! Return if we are in the stratosphere
!      ! This skips unneccesary computations & function calls
!      IF ( StratBox ) RETURN
!
!      ! Continue initializing
!      X1 = 0.0_dp
!
!      ! Cloud halide concentration, put fine and coarse mode together
!      clConc = clConc_A + clConc_C + clConc_g
!      brConc = brConc_A + brConc_C + brConc_g
!
!      If (X == 1) THEN
!         r_ac = clConc_A / clConc
!         Y = 1
!      ELSEIF (X == 2) THEN
!         r_ac = clConc_C / clConc
!         Y = 1
!      ElSEIF (X == 3) THEN
!         r_ac = clConc_g / clConc
!         Y = 1
!      ELSEIF (X == 4) THEN
!         r_ac = brConc_A / brConc
!         Y = 2
!      ELSEIF (X == 5) THEN
!         r_ac = brConc_C / brConc
!         Y = 2
!      ELSEIF (X == 6) THEN
!         r_ac = brConc_g / brConc
!         Y = 2
!      ENDIF
!
!      CALL GAMMA_ClNO2( rLiq,   denAir, TK,        pH,   clConc,             &
!                        brConc, Y,      GAM_ClNO2, r_gp, H                  )
!      X1   = GAM_ClNO2
!      r1   = r_gp * r_ac
!      X2   = 0.0_fp
!      rate = CloudHet( CldFr,            Aliq, Aice, rLiq, rIce,             &
!                       SR_MW(ind_ClNO2), X1,   X2,   r1,   1.0_dp           )
!
    END FUNCTION HetClNO2_TCld
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetClNO2
!
! !DESCRIPTION: Sets the ClNO2 + Cl- (in fine mode) rate.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETClNO2( denAir, rAer,   AAer,   alkAer, TK,                   &
                       pH,     clConc, brConc, X,      H   ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(dp),       INTENT(IN) :: denAir   ! Density of air (#/cm3)
      REAL(dp),       INTENT(IN) :: rAer     ! Radius of aerosol (cm)
      REAL(dp),       INTENT(IN) :: AAer     ! Area of aerosol (cm2/cm3)
      REAL(dp),       INTENT(IN) :: alkAer   ! Aerosol alkalinity (?)
      REAL(dp),       INTENT(IN) :: TK       ! Temperature (K)
      REAL(dp),       INTENT(IN) :: clConc   ! Cloride concentration (mol/L)
      REAL(dp),       INTENT(IN) :: brConc   ! Bromide concentration (mol/L)
      REAL(dp),       INTENT(IN) :: pH       ! Aerosol pH
      INTEGER,        INTENT(IN) :: X        ! 1: Cl-, 2: Br-
      TYPE(HetState), POINTER    :: H        ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(dp)                   :: kISum
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(dp) :: XSQM, GAM_ClNO2, r_gp

!      ! Initialize
!      kISum     = 0.0_dp
!      GAM_ClNO2 = 0.0_dp
!      r_gp      = 0.0_dp
!      XSQM      = SQRT( H%ClNO2%MW_g )
!
!      CALL GAMMA_ClNO2( rAer,   denAir, TK,        pH,   clConc,             &
!                        brConc, X,      GAM_ClNO2, r_gp, H                  )
!
!      ! Reaction rate for surface of aerosol
!      kISum = Arsl1K( AAer, rAer, denAir, GAM_ClNO2, SR_TEMP, XSqM ) * r_gp

    END FUNCTION HETClNO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
    SUBROUTINE GAMMA_HOCl_CLD( Radius,   n_air, X,    T,                     &
                               C_Y1,     C_Y3,  C_Y4, C_Hp,                  &
                               GAM_HOCl, r_gp,  H                           )
!
! !USES:
!
      USE PhysConstants,      ONLY : Pi, RStarG
!
! !INPUT PARAMETERS:
!
      ! Radius (cm), n_air (#/cm), and X (1 for Cl and 2 for Br)
      REAL(dp),       INTENT(IN)  :: Radius   ! Radius (cm)
      REAL(dp),       INTENT(IN)  :: n_air    ! n_air (#/cm)
      INTEGER,        INTENT(IN)  :: X        ! 1=Cl-,2=HSO3-/SO3--
      REAL(dp),       INTENT(IN)  :: T        ! Temperature (K)
      REAL(dp),       INTENT(IN)  :: C_Y1     ! Concentration (mol/L)
      REAL(dp),       INTENT(IN)  :: C_Y3     ! Concentration (mol/L)
      REAL(dp),       INTENT(IN)  :: C_Y4     ! Concentration (mol/L)
      REAL(dp),       INTENT(IN)  :: C_Hp     ! Concentration (mol/L)
      TYPE(HetState), POINTER     :: H        ! Hetchem species metadata
!
! !OUTPUT PARAMETERS:
!
      ! Reactive uptake coefficient (unitless)
      REAL(dp),       INTENT(OUT) :: GAM_HOCl
      REAL(dp),       INTENT(OUT) :: r_gp
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!      ! Scalars
!      REAL(dp) :: ab,     M_X,  k1,  k2,  H_X,   C_Y2
!      REAL(dp) :: gb_tot, cavg, D_l, l_r, k_tot, H_HOCl
!
!      C_Y2     = C_Y3 + C_Y4
!
!      ! MW of HOCl (kg/mol)
!      M_X      = H%HOCl%MW_g * 1.0e-3_dp
!
!      ! Mass accommodation coefficient
!      ab       = 0.8_dp
!
!      ! thermal velocity (cm/s)
!      cavg     = SQRT( 8.0_dp * RStarG * T / ( Pi * M_X) ) *1.0e2_dp
!
!      ! Liquid phase diffusion coefficient [cm2/s] for HOCl
!      ! (Ammann et al., Atmos. Chem. Phys., 2013)
!      D_l      = 2.0e-5_dp
!
!      k1       = 1.5e+4_dp !M-1 s-1
!      k2       = 2.8e+5_dp !M-1 s-1 Liu and Abbatt, Geophys. Res. Lett., 2020
!      k_tot    = k1 * C_Hp * C_Y1                                            &
!               + k2 * C_Y2
!
!      ! Henry's law
!      H_HOCl   = H%HOCl%K0 * con_atm_bar
!      H_X      = H_HOCl * EXP( H%HOCl%CR * ( 1.0_dp/T - 1.0_dp/H%HOCl%TK ) )
!
!      l_r      = SQRT( D_l / k_tot )
!      gb_tot   = 4.0_dp * H_X * con_R * T * l_r * k_tot / cavg
!      gb_tot   = gb_tot * REACTODIFF_CORR( Radius, l_r )
!
!      ! Reactive uptake coefficient [unitless]
!      GAM_HOCl = 1.0_dp / ( 1.0_dp/ab  +  1.0_dp/gb_tot )
!
!      ! turn off HOCl+S(IV), Xuan Wang
!      !gb2 = 0.0e0_dp
!
!      IF ( X==1 ) THEN
!         r_gp = k1 * C_Hp * C_Y1  /k_tot
!      ELSEIF ( X==2 ) THEN
!         r_gp = k2 * C_Y2 / k_tot
!      ENDIF

    END SUBROUTINE GAMMA_HOCl_CLD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gamma_HOCl_AER
!
! !DESCRIPTION: Function GAMMA\_HOCl calculates reactive uptake coef.
!               for HOCL + Cl-
!\\
!\\
! !INTERFACE:
!
    FUNCTION GAMMA_HOCl_AER( Radius, n_air, T, C_H, C_X, H ) RESULT(GAM_HOCl)
!
! !USES:
!
      USE PhysConstants, ONLY : Pi, RStarG
!
! !INPUT PARAMETERS:
!
      REAL(dp),       INTENT(IN) :: Radius   ! Radius (cm)
      REAL(dp),       INTENT(IN) :: n_air    ! n_air (#/cm)
      REAL(dp),       INTENT(IN) :: T        ! Temperature (K)
      REAL(dp),       INTENT(IN) :: C_H      ! H+ concentration
      REAL(dp),       INTENT(IN) :: C_X      ! Cl- Concentration (mol/L)
      TYPE(HetState), POINTER    :: H        ! Hetchem species metadata
!
! !RETURN VALUE:
!
      ! Reactive uptake coefficient (unitless)
      REAL(dp)                   :: GAM_HOCl
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!      ! Scalars
!      REAL(dp) :: ab, M_X, k_ter, H_X, H_HOCl, cavg, D_l, gb, l_r
!
!      ! MW of HOCl (kg/mol)
!      M_X      = H%HOCl%MW_g * 1.0e-2_dp
!
!      ! Mass accommodation coefficient
!      ab       = 0.8_dp
!
!      ! thermal velocity (cm/s)
!      cavg     = SQRT( 8.0_dp * RStarG * T / ( Pi * M_X ) ) *1.0e2_dp
!
!      ! Liquid phase diffusion coefficient [cm2/s] for HOCl
!      ! (Ammann et al., Atmos. Chem. Phys., 2013)
!      D_l      = 2.0e-5_dp
!
!      k_ter    = 1.5e+4_dp                 ! M-1s-1
!      H_HOCl   = H%HOCl%K0 * con_atm_bar   ! M/bar
!      H_X      = H_HOCl*EXP( H%HOCl%CR *( 1.0_dp/T - 1.0_dp/H%HOCl%TK ) )
!
!      l_r = SQRT(D_l / (k_ter * C_H * C_X))
!      gb       = 4.0_dp * H_X * con_R * T * l_r * k_ter * C_H * C_X / cavg
!      gb       = gb * REACTODIFF_CORR( Radius, l_r)
!
!      GAM_HOCl = 1.0_dp / (1.0_dp/ab  +  1.0_dp/gb)
!
!      !GAM_HOCl = MIN(GAM_HOCl, 2.0e-4_dp)

    END FUNCTION GAMMA_HOCl_AER

!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gamma_ClNO2
!
! !DESCRIPTION: Function GAMMA\_ClNO2 calculates reactive uptake coef.
!               for ClNO2 + Cl-/Br-
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE GAMMA_ClNO2( Radius, n_air, T,         pH,   C_X1,            &
                            C_X2,   X,     GAM_ClNO2, r_gp, H               )
!
! !USES:
!
      USE PhysConstants, ONLY : Pi, RStarG
!
! !INPUT PARAMETERS:
!
      REAL(dp),       INTENT(IN)   :: Radius   ! Radius (cm)
      REAL(dp),       INTENT(IN)   :: n_air    ! n_air (#/cm)
      REAL(dp),       INTENT(IN)   :: T        ! Temperature (K)
      REAL(dp),       INTENT(IN)   :: C_X1     ! Cl-/Br- Concentration (mol/L)
      REAL(dp),       INTENT(IN)   :: C_X2     ! Cl-/Br- Concentration (mol/L)
      INTEGER,        INTENT(IN)   :: X        ! 1=Cl-,2=Br-
      REAL(dp),       INTENT(IN)   :: pH
      TYPE(HetState), POINTER      :: H        ! Hetchem species metadata
!
! !OUTPUT PARAMETERS:
!
      ! Reactive uptake coefficient (unitless)
      REAL(dp),       INTENT(OUT)  :: GAM_ClNO2
      REAL(dp),       INTENT(OUT)  :: r_gp
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!
!      ! Scalars
!      REAL(dp) :: ab, M_X, H_X, fCl, k_tot
!      REAL(dp) :: cavg, D_l, k_b1, k_b2, l_r, gb_tot
!
!      ! MW of ClNO2 (kg/mol)
!      M_X       = H%ClNO2%MW_g * 1.0e-3_dp
!
!      ! Mass accommodation coefficient
!      ab        = 0.01_dp
!
!      ! thermal velocity (cm/s)
!      cavg      = SQRT( 8.0_dp * RStarG * T / ( Pi * M_X ) ) *1.0e2_dp
!
!      ! Liquid phase diffusion coefficient [cm2/s] for ClNO2
!      ! (Ammann et al., Atmos. Chem. Phys., 2013)
!      D_l       = 1.0e-5_dp
!      H_X       = 4.5e-2_dp         !M atm-1
!      H_X       = H_X * con_atm_bar !M/bar
!
!      !ClNO2+Cl-
!      k_b1      = 1.0e+7_dp !M-1s-1
!      IF (pH >= 2) THEN
!         k_b1   = 0.0e0_dp
!      ENDIF
!
!      !ClNO2+Br-
!      k_b2      = 1.01e-1_dp / (H_X*H_X*D_l)
!
!      k_tot     = k_b1*C_X1 + k_b2*C_X2
!      l_r       = SQRT(D_l / k_tot)
!      gb_tot    = 4.0_dp * H_X * con_R * T * l_r * k_tot / cavg
!      gb_tot    = gb_tot * REACTODIFF_CORR( Radius, l_r )
!
!      GAM_ClNO2 = 1.0_dp / ( 1.0_dp/ab  +  1.0_dp/gb_tot )
!
!      IF ( X==1 ) THEN
!         r_gp = k_b1*C_X1 / k_tot
!      ELSEIF ( X==2 ) THEN
!         r_gp = k_b2*C_X2 / k_tot
!      ENDIF

    END SUBROUTINE GAMMA_ClNO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetClNO3_SS
!
! !DESCRIPTION: Sets the ClNO3 + Br- (in sea salt) rate using Johan
!  Schmidt's updated code.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HETClNO3_SS( rAer,   AAer, alkAer, clConc, brConc, X, M           &
                      ) RESULT( rate )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: rAer     ! Radius of aerosol (cm)
    REAL(dp), INTENT(IN) :: aAer     ! Area of aerosol (cm2/cm3)
    REAL(dp), INTENT(IN) :: alkAer   ! Aerosol alkalinity (?)
    REAL(dp), INTENT(IN) :: clConc   ! Cloride concentration (mol/L)
    REAL(dp), INTENT(IN) :: brConc   ! Bromide concentration (mol/L)
    INTEGER,  INTENT(IN) :: X        ! 1: Cl-, 2: Br-
    INTEGER,  INTENT(IN) :: M        ! 1: fine, 2: coarse
!
! !RETURN VALUE:
!
    REAL(dp)             :: rate     ! Rxn rate ClNO3 + Br- in sea salt
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER  :: N
    REAL(dp) :: gamma, r_gp, srMw

!    ! Initialize
!    rate = 0.0_dp
!    srMw = SR_MW(ind_ClNO3)
!
!    ! Reaction probability [1]
!    CALL Gamma_ClNO3_AER( rAer, X, M, clConc, brConc, gamma, r_gp )
!
!    ! Reaction rate for surface of aerosol
!    rate = Arsl1K( aAer, rAer, NUMDEN, gamma, SR_TEMP, srMw ) * r_gp

  END FUNCTION HETClNO3_SS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOBr_SS
!
! !DESCRIPTION: Sets the HOBr + Br- or Cl- (in sea salt) rate using Johan
!  Schmidt's updated code.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HETHOBr_SS( denAir, rAer,   AAer,   alkAer, TK,                 &
                       hConc,  clConc, brConc, X,      H                 ) &
                       RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(dp),       INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(dp),       INTENT(IN) :: rAer        ! Radius of aerosol (cm)
      REAL(dp),       INTENT(IN) :: AAer        ! Area of aerosol (cm2/cm3)
      REAL(dp),       INTENT(IN) :: alkAer      ! Aerosol alkalinity (?)
      REAL(dp),       INTENT(IN) :: TK          ! Temperature (K)
      REAL(dp),       INTENT(IN) :: hConc       ! H+ concentration (mol/L)
      REAL(dp),       INTENT(IN) :: clConc      ! Cloride concentration (mol/L)
      REAL(dp),       INTENT(IN) :: brConc      ! Bromide concentration (mol/L)
      INTEGER,        INTENT(IN) :: X           ! 1: Cl-, 2: Br-
      TYPE(HetState), POINTER    :: H           ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(dp)             :: kISum
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!      REAL(dp) :: XSTKCF, GAM_HOBr, r_gp, XSQM
!
!      ! Initialize
!      kISum = 0.0_dp
!      XSQM  = SQRT( H%HOBr%MW_g )
!
!      ! Reaction can only proceed on acidic aerosol
!      IF ( alkAer > 0.05_dp ) THEN
!         XStkCf = 0.0_dp
!         r_gp   = 0.0_dp
!      ELSE
!         CALL Gamma_HOBr_AER( rAer,   denAir,  X,       TK,   clConc,        &
!                              brConc, hConc,  GAM_HOBr, r_gp, H             )
!         XStkCf = GAM_HOBr
!      ENDIF
!
!      ! Reaction rate for surface of aerosol
!      kISum = Arsl1K( AAer, rAer, denAir, XStkCf, SR_TEMP, XSqM ) * r_gp

    END FUNCTION HETHOBr_SS
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: HetO3_TCld
!!
!! !DESCRIPTION: Sets the O3 + Br- rate in troposphere cloud
!!
!!\\
!!\\
!! !INTERFACE:
!!
    FUNCTION HETO3_TCld( denAir,   rLiq,   rIce,  ALiq,     AIce,            &
                         VAir,     TK,     CldFr, brConc_a, brConc_c,        &
                         brConc_g, O3Conc, X,     H                        ) &
                         RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(dp),       INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(dp),       INTENT(IN) :: rLiq        ! Radius of liquid cloud
                                                !  droplets (cm)
      REAL(dp),       INTENT(IN) :: rIce        ! Radius of ice cloud crystals
                                                !  (cm)
      REAL(dp),       INTENT(IN) :: ALiq        ! Area of liquid cloud droplets
                                                !  (cm2/cm3)
      REAL(dp),       INTENT(IN) :: AIce        ! Area of ice cloud crystals
                                                !  (cm2/cm3)
      REAL(dp),       INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(dp),       INTENT(IN) :: TK          ! Temperature (K)
      REAL(dp),       INTENT(IN) :: CldFr       ! Cloud fraction
      REAL(dp),       INTENT(IN) :: brConc_a    ! Bromide concentration (mol/L)
      REAL(dp),       INTENT(IN) :: brConc_c    ! Bromide concentration (mol/L)
      REAL(dp),       INTENT(IN) :: brConc_g    ! Bromide concentration (mol/L)
      REAL(dp),       INTENT(IN) :: O3Conc      ! Ozone concentration (mol/L)
      INTEGER,        INTENT(IN) :: X           ! X=1 _a; 2 _c; 3 _g
      TYPE(HetState), POINTER    :: H           ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(dp)                   :: kISum
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      ! Scalars
!      INTEGER  :: N
!      REAL(dp) :: brConc, r_ac, X1, X2
!
!      ! Initialize the return value
!      kISum  = 0.0_dp
!
!      ! Return if we are in the stratosphere
!      ! This skips unneeded computations & function calls
!      IF ( StratBox ) RETURN
!
!      ! Continue initializing
!      r_ac   = 0.0_dp
!      brConc = 0.0_dp
!      X1     = 0.0_dp
!
!      brConc = brConc_a + brConc_g + brConc_c
!
!      IF (X.EQ.1) THEN
!         r_ac = brConc_a / brConc
!      ELSE IF (X .EQ. 2) THEN
!         r_ac = brConc_c / brConc
!      ELSE
!         r_ac = brConc_g / brConc
!      ENDIF
!
!      ! Reaction on liquid clouds (tropospheric only)
!      X1    = Gamma_O3_Br( rLiq, denAir, TK, brConc, O3Conc, H )
!      X2    = 0.0_dp
!      kISum = CloudHet( CldFr,         Aliq, Aice, rLiq,  rIce,              &
!                        SR_MW(ind_O3), X1,   X2,   r_ac,  1.0_dp            )
!
    END FUNCTION HETO3_TCld
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gamma_O3_Br
!
! !DESCRIPTION: Function GAMMA\_O3\_Br calculates reactive uptake coef.
!               for bromide oxidation by O3
!\\
!\\
! !INTERFACE:

    FUNCTION GAMMA_O3_Br( Radius, n_air, T, C_Y, C_X_g, H ) RESULT( GAM )
!
! !USES:
!
      USE PhysConstants, ONLY : Pi, RStarG
!
! !INPUT PARAMETERS:
!
      REAL(dp),       INTENT(IN) :: Radius   ! Radius (cm)
      REAL(dp),       INTENT(IN) :: n_air    ! Number density of air (#/cm)
      REAL(dp),       INTENT(IN) :: T        ! Temparature (K)
      REAL(dp),       INTENT(IN) :: C_Y      ! Concentration
      REAL(dp),       INTENT(IN) :: C_X_g    ! Gas-phase concentration
      TYPE(HetState), POINTER    :: H        ! Hetchem species metadata
!
! !RETURN VALUE
!
      REAL(dp)                   :: GAM      ! Reactive uptake coeff. (1)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      REAL(dp) :: ab,     gb,  gd,       gs
      REAL(dp) :: cavg,   H_X, H_O3,     M_X
      REAL(dp) :: KLangC, k_s, C_Y_surf, Nmax
      REAL(dp) :: k_b,    D_l, l_r

!      ! Henry's law
!      H_O3     = H%O3%K0 * con_atm_bar
!      H_X      = H_O3*EXP( H%O3%CR * ( 1.0_dp/T - 1.0_dp/H%O3%TK ) )
!
!      ! Molwt of O3 (kg/mol)
!      M_X      = H%O3%MW_g * 1.0e-3_dp
!
!      ! Thermal velocity (cm/s)
!      cavg     = SQRT( 8 * RStarG * T / ( Pi * M_X ) ) *1.0e2_dp
!
!      Nmax     = 3.0e14_dp  ! #/cm2
!      KLangC   = 1.0e-13_dp !cm3
!      k_s      = 1.0e-16_dp !cm2s-1, from ks*Nmax=0.03s-1
!
!     ! [Br-(surf)] = 3.41E14 cm-2/M * [Br-(bulk)], but not gt Nmax.
!      C_Y_surf = MIN( 3.41e14_dp * C_Y, Nmax )
!      gs       = ( 4.0_dp * k_s * C_Y_surf * KLangC * Nmax )                 &
!               / ( cavg * ( 1.0_dp + KLangC * C_X_g )      )
!
!      k_b      = 6.3e8_dp *  EXP(-4.45e3_dp / T) !M-1 s-1
!      D_l      = 8.9e-6_dp !cm2 s-1.
!      l_r      = SQRT( D_l / (k_b * C_Y ) )! cm
!      gb       = 4.0_dp * H_X * con_R * T * l_r * k_b * C_Y / cavg
!      gb       = gb * REACTODIFF_CORR( Radius, l_r)
!
!      GAM      = gb + gs

      END FUNCTION GAMMA_O3_Br
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetClNO3_HBr
!
! !DESCRIPTION: Sets the ClNO3 + Br- rate using Johan Schmidt's
!  updated code.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETClNO3_HBr( A, B, Input_Opt ) RESULT( kISum )
!
! !USES:
!
      USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(dp),       INTENT(IN) :: A, B
      TYPE(OptInput), INTENT(IN) :: Input_Opt
!
! !RETURN VALUE:
!
      REAL(dp)                   :: kISum
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!      INTEGER  :: N
!      REAL(dp) :: XSTKCF, ADJUSTEDRATE
!      Real(dp) :: XSQM
!
!      ! Initialize
!      kISum        = 0.0_dp
!      ADJUSTEDRATE = 0.0_dp
!      XSTKCF       = 0.0_dp
!      XSQM         =sqrt(A)
!
!      ! Loop over aerosol types
!      DO N = 1, NAEROTYPE
!
!         ! Get the aerosol type
!         XStkCf = 0.0e+0_dp
!
!         IF ( N == 8 ) THEN
!
!            ! sulfate aerosol
!            ! This seems not to happen
!
!         ELSEIF ( Input_Opt%LUCX .and. STRATBOX ) THEN
!            IF (N.eq.13) THEN
!               XSTKCF = KHETI_SLA(5)
!            ELSEIF (N.eq.14) THEN
!               IF (NATSURFACE) THEN
!                  XSTKCF = 0.3e+0_dp ! NAT
!               ELSE
!                  XSTKCF = 0.3e+0_dp ! Ice
!               ENDIF
!            ENDIF
!
!         ENDIF
!
!         IF (XStkCf.gt.0.0e+0_dp) THEN
!            IF (N.eq.13) THEN
!               ! Calculate for stratospheric liquid aerosol
!               ! Note that XSTKCF is actually a premultiplying
!               ! factor in this case, including c-bar
!               ADJUSTEDRATE = XAREA(N) * XSTKCF
!            ELSE
!               ! Reaction rate for surface of aerosol
!               ADJUSTEDRATE=ArsL1k(XAREA(N),XRADI(N),NUMDEN,XSTKCF,SR_TEMP,XSQM)
!            ENDIF
!
!            ! Add to overall reaction rate
!            kISum = kISum + ADJUSTEDRATE
!         ENDIF
!      END DO

    END FUNCTION HETClNO3_HBr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetClNO3
!
! !DESCRIPTION: Sets the hydrolysis rate for ClNO3 using Johan Schmidt's
!  updated code.
!\\
!\\
! !INTERFACE:
!
  FUNCTION HetClNO3( clf, srMw, clConc, brConc ) RESULT( rate )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: clf      ! Clear-sky fraction
    REAL(dp), INTENT(IN) :: srMw     ! SQRT( molec weight [g/mol] )
    REAL(dp), INTENT(IN) :: clconc   ! Cl- concentrtaion, fine mode
    REAL(dp), INTENT(IN) :: brconc   ! Br- concentration, fine mode
!
! !RETURN VALUE:
!
    REAL(dp)             :: rate     ! Hydrolysis. rate for ClNO3 [1/s]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL(dp) :: gamma, r_gp

!    ! Initialize
!    gamma = 0.0_dp
!    rate  = 0.0_dp
!    r_gp  = 0.0_dp
!
!    ! Fine sea salt (aerosol type #11)
!    CALL Gamma_ClNO3_Aer( AClRADI, 3, 1, clconc, brconc, gamma, r_gp        )
!    rate = rate +                                                            &
!         ArsL1K( clf*AClAREA, AClRADI, NUMDEN, gamma, SR_TEMP, srMw ) * r_gp
!
!    ! Stratospheric liquid aerosol (aerosol type #13)
!    rate = rate + XAREA(13) * KHETI_SLA(3)
!
!    ! Irregular ice cloud (aerosol type #14)
!    gamma = 0.3_fp                                  ! Rxn prob, ice [1]
!    IF ( NatSurface ) gamma = 0.004_fp              ! Rxn prob, NAT [1]
!    rate = rate + ArsL1K(XAREA(14), XRADI(14), NUMDEN, gamma, SR_TEMP, srMw )

  END FUNCTION HetClNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOBr_HCl
!
! !DESCRIPTION: Sets the rate of the multiphase reaction HOBr + Cl- in
!  sulfate aerosols, on cloud droplets and on PSCs
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOBr_HCl( A, B, Input_Opt ) RESULT( kISum )
!
! !USES:
!
      USE Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
      ! Rate coefficients
      REAL(dp),       INTENT(IN) :: A, B
      TYPE(OptInput), INTENT(IN) :: Input_Opt
!
! !RETURN VALUE:
!
      REAL(dp)                   :: kISum
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(dp) :: XSTKCF, ADJUSTEDRATE
      Real(dp) :: XSQM
      REAL(dp) :: GAM_HOBr, r_gp

!      ! Initialize
!      kISum        = 0.0_dp
!      ADJUSTEDRATE = 0.0_dp
!      XSTKCF       = 0.0_dp
!      XSQM         = SQRT(A)
!
!      ! Loop over aerosol types
!      DO N = 1, NAEROTYPE
!
!         XSTKCF = 0e+0_dp
!         ! Get the aerosol type
!
!         IF ( Input_Opt%LUCX .and. STRATBOX) THEN
!            ! add limitation to stratosphere, xnw 1/25/18
!            IF ( N == 8 ) THEN
!               XSTKCF = 0.2e+0_dp ! Sulfate, [Hanson and Ravishankara, 1995]
!            ELSEIF (N.eq.13) THEN
!               ! SSA/STS
!               XSTKCF = KHETI_SLA(10)
!            ELSEIF (N.eq.14) THEN
!               ! Ice/NAT PSC
!               IF (NATSURFACE) THEN
!                  XSTKCF = 0.1e+0_dp ! NAT
!               ELSE
!                  XSTKCF = 0.3e+0_dp ! Ice
!               ENDIF
!            ELSE
!               XSTKCF = 0e+0_dp
!            ENDIF
!         ENDIF
!
!         IF (N.eq.13) THEN
!            ! Calculate for stratospheric liquid aerosol
!            ! Note that XSTKCF is actually a premultiplying
!            ! factor in this case, including c-bar
!            ADJUSTEDRATE = XAREA(N) * XSTKCF
!         ELSEIF (XStkCf.gt.0.0e+0_dp) THEN
!            ! Reaction rate for surface of aerosol
!            ADJUSTEDRATE=ArsL1k(XAREA(N),XRADI(N),NUMDEN,XSTKCF,SR_TEMP,XSQM)
!         ENDIF
!
!         ! Add to overall reaction rate
!         kISum = kISum + ADJUSTEDRATE
!      END DO

    END FUNCTION HETHOBr_HCl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOBr_HBr
!
! !DESCRIPTION: Sets the rate of the multiphase reaction HOBr + Br- in
!  sulfate aerosols, on cloud droplets and on PSCs
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOBr_HBr( A, B, Input_Opt ) RESULT( kISum )
!
! !USES:
!
      USE Input_Opt_Mod, ONLY : OptInput
!
! INPUT PARAMETERS:
!
      REAL(dp),       INTENT(IN) :: A, B        ! Rate coefficients
      TYPE(OptInput), INTENT(IN) :: Input_Opt
!
! !RETURN VALUE:
!
      REAL(dp)                   :: kISum
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(dp) :: XSTKCF, ADJUSTEDRATE
      Real(dp) :: XSQM
      Real(dp) :: SADen

      ! Initialize
      kISum        = 0.0_dp
      ADJUSTEDRATE = 0.0_dp
      XSTKCF       = 0.0_dp
      XSQM=SQRT(A)

!      ! Loop over aerosol types
!      DO N = 1, NAEROTYPE
!         XSTKCF = 0e+0_dp
!
!         ! Get the aerosol type
!         IF ( Input_Opt%LUCX .and. STRATBOX ) THEN
!            ! add limitation to stratosphere, xnw 1/25/18
!            IF ( N == 8 ) THEN
!               XSTKCF = 0.25e+0_dp ! Sulfate, [Abbatt, 1995]
!            ELSEIF ( N == 13 ) THEN
!               ! SSA/STS
!               XSTKCF = KHETI_SLA(6)
!            ELSEIF ( N == 14 ) THEN
!               ! Ice/NAT PSC
!               IF (NATSURFACE) THEN
!                  XSTKCF = 0.001e+0_dp
!               ELSE
!                  XSTKCF = 0.3e+0_dp
!               ENDIF
!            ELSE
!               XSTKCF = 0e+0_dp
!            ENDIF
!         ENDIF
!
!         IF (N.eq.13) THEN
!            ! Calculate for stratospheric liquid aerosol
!            ! Note that XSTKCF is actually a premultiplying
!            ! factor in this case, including c-bar
!            ADJUSTEDRATE = XAREA(N) * XSTKCF
!         ELSEIF (XStkCf.gt.0.0e+0_dp) THEN
!            ! Reaction rate for surface of aerosol
!            ADJUSTEDRATE=ArsL1k(XAREA(N),XRADI(N),NUMDEN,XSTKCF,SR_TEMP,XSQM)
!         ENDIF
!
!         ! Add to overall reaction rate
!         kISum = kISum + ADJUSTEDRATE
!      END DO

    END FUNCTION HETHOBr_HBr
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: HetClNO3_TCld
!!
!! !DESCRIPTION: Sets the rate of the multiphase reaction ClNO3 + Cl-/Br- in
!!  troposphere cloud
!!\\
!!\\
!! !INTERFACE:
!!
  FUNCTION HETClNO3_TCld( rLiq,     rIce,     ALiq,     AIce,                &
                          VAir,     CldFr,    clConc_A, clConc_C,            &
                          clConc_g, brConc_A, brConc_C, brConc_g,            &
                          hno3_th,  hcl_th,   hbr_th,   X                    &
                        ) RESULT( rate )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: rLiq      ! Radius of liquid cloud droplets (cm)
    REAL(dp), INTENT(IN) :: rIce      ! Radius of ice cloud crystals (cm)
    REAL(dp), INTENT(IN) :: ALiq      ! Area of liquid cloud droplets (cm2/cm3)
    REAL(dp), INTENT(IN) :: AIce      ! Area of ice cloud crystals (cm2/cm3)
    REAL(dp), INTENT(IN) :: VAir      ! Box volume (cm3)
    REAL(dp), INTENT(IN) :: CldFr     ! Cloud fraction
    REAL(dp), INTENT(IN) :: clConc_A  ! Fine chloride concentration (mol/L)
    REAL(dp), INTENT(IN) :: clConc_C  ! Coarse chloride concentration (mol/L)
    REAL(dp), INTENT(IN) :: clConc_g
    REAL(dp), INTENT(IN) :: brConc_A  ! Fine bromide concentration (mol/L)
    REAL(dp), INTENT(IN) :: brConc_C  ! Coarse bromide concentration (mol/L)
    REAL(dp), INTENT(IN) :: brConc_g
    REAL(dp), INTENT(IN) :: hno3_th
    REAL(dp), INTENT(IN) :: hcl_th
    REAL(dp), INTENT(IN) :: hbr_th
    INTEGER,  INTENT(IN) :: X         ! 1=fineCl; 2=coarseCl; 3=fineBr;
                                      ! 4=coarseBr; 5=HCl; 6=HBr; 7=H2O
!
! !RETURN VALUE:
!
    REAL(dp)             :: rate      ! Reaction rate [1/s]
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!    ! Scalars
!    INTEGER  :: N, Y
!    REAL(dp) :: X1, X2, gamma, r1, r2, clConc, r_ac, brConc, r_gp
!
!    ! Initialize the return value
!    rate = 0.0_dp
!
!    ! Return if we are in the stratosphere
!    ! This will prevent unnecessary computations & function calls
!    IF ( StratBox ) RETURN
!
!    ! Continue initializing
!    X1 = 0.0_dp
!    X2 = 0.0_dp
!
!    ! Cloud halide concentration, put fine and coarse mode together
!    clConc = clConc_A + clConc_C + clConc_g
!    brConc = brConc_A + brConc_C + brConc_g
!
!    SELECT CASE( X )
!       CASE( 1 )
!          Y = 1
!          r_ac = clConc_A / clConc
!       CASE( 2 )
!          Y = 1
!          r_ac = clConc_C / clConc
!       CASE( 3 )
!          Y = 2
!          r_ac = brConc_A / brConc
!       CASE( 4 )
!          Y = 2
!          r_ac = brConc_C / brConc
!       CASE( 5 )
!          Y = 1
!          r_ac = clConc_g / clConc
!       CASE( 6 )
!          Y = 2
!          r_ac = brConc_g / brConc
!       CASE DEFAULT  ! X = 7
!          Y = 3
!          r_ac = 1.0_dp
!    END SELECT
!
!    ! Reaction probability on aerosol [1]
!    CALL Gamma_ClNO3_Aer( rLiq, Y, 3, clConc, brConc, gamma, r_gp )
!    X1 = gamma
!    r1 = r_gp * r_ac
!
!    ! Reaction probability on ice [1]
!    CALL Gamma_ClNO3_Ice( Y, hno3_th, hcl_th, hbr_th, gamma, r_gp )
!    X2 = gamma
!    r2 = 0.0_dp
!    IF ( X >= 5 ) r2 = r_gp
!
!    ! Reaction rate
!    rate = CloudHet( CldFr,            Aliq, Aice, rLiq, rIce,               &
!                     SR_MW(ind_ClNO3), X1,   X2,   r1,   r2                 )
!
  END FUNCTION HetClNO3_TCld
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: HetHOBr_TCld
!!
!! !DESCRIPTION: Sets the rate of the multiphase reaction HOBr + Cl-/Br-/S(IV) in
!!  troposphere cloud
!!\\
!!\\
!! !INTERFACE:
!!
    FUNCTION HETHOBr_TCld( denAir,    rLiq,      rIce,      ALiq,            &
                           AIce,      VAir,      TK,        CldFr,           &
                           hConc_Sul, hConc_LCl, hConc_ICl, clConc_A,        &
                           clConc_C,  clConc_g,  brConc_A,  brConc_C,        &
                           brConc_g,  hso3Conc,  so3Conc,   hno3_th,         &
                           hcl_th,    hbr_th,    X,         H              ) &
                           RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(dp),       INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(dp),       INTENT(IN) :: rLiq        ! Radius of liquid cloud
                                                !  droplets (cm)
      REAL(dp),       INTENT(IN) :: rIce        ! Radius of ice cloud crystals
                                                !  (cm)
      REAL(dp),       INTENT(IN) :: ALiq        ! Area of liquid cloud droplets
                                                !  (cm2/cm3)
      REAL(dp),       INTENT(IN) :: AIce        ! Area of ice cloud crystals
                                                !  (cm2/cm3)
      REAL(dp),       INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(dp),       INTENT(IN) :: TK          ! Temperature (K)
      REAL(dp),       INTENT(IN) :: CldFr       ! Cloud fraction
      REAL(dp),       INTENT(IN) :: hConc_Sul   ! Sulfate H+ concentration
      REAL(dp),       INTENT(IN) :: hConc_LCl   ! Liquid cloud H+ concentration
      REAL(dp),       INTENT(IN) :: hConc_ICl   ! Ice cloud H+ concentration
      REAL(dp),       INTENT(IN) :: clConc_A    ! Fine Chloride concentration
                                                !  (mol/L)
      REAL(dp),       INTENT(IN) :: clConc_C    ! Coarse Chloride concentration
                                                !  (mol/L)
      REAL(dp),       INTENT(IN) :: clConc_g
      REAL(dp),       INTENT(IN) :: brConc_A    ! Fine Bromide concentration
                                                !  (mol/L)
      REAL(dp),       INTENT(IN) :: brConc_C    ! Coarse Bromide concentration
                                                ! (mol/L)
      REAL(dp),       INTENT(IN) :: brConc_g
      REAL(dp),       INTENT(IN) :: hso3Conc    ! HSO3-    concentration (mol/L)
      REAL(dp),       INTENT(IN) :: so3Conc     ! SO3--    concentration (mol/L)
      REAL(dp),       INTENT(IN) :: hno3_th
      REAL(dp),       INTENT(IN) :: hcl_th
      REAL(dp),       INTENT(IN) :: hbr_th
      INTEGER,        INTENT(IN) :: X           ! 1=fineCl;2=coarseCl;
                                                ! 3=fineBr;4=coarseBr;5=HCl;
                                                ! 6=HBr;7=HSO3;8=SO3
      TYPE(HetState), POINTER    :: H           ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(dp)                   :: kISum
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      INTEGER  :: N,        Y
!      REAL(dp) :: X1,       X2,   ADJUSTEDRATE, XSQM
!      REAL(dp) :: GAM_HOBr, r_gp, clConc,       r_ac
!      REAL(dp) :: brConc,   r1,   r2
!
!      ! Initialize return value
!      kISum        = 0.0_dp
!
!      ! Return if we are in the stratosphere
!      ! This skips unneccessary computations & function calls
!      IF ( StratBox ) RETURN
!
!      ! Continue initializing
!      ADJUSTEDRATE = 0.0_dp
!      X1           = 0.0_dp
!      X2           = 0.0_dp
!      XSqM         = SR_MW(ind_HOBr)
!
!      ! Cloud halide concentration, put fine and coarse mode together
!      clConc = clConc_A + clConc_C + clConc_g
!      brConc = brConc_A + brConc_C + brConc_g
!
!      If (X == 1) THEN
!         Y = 1
!         r_ac = clConc_A / clConc
!      ELSEIF (X == 2) THEN
!         Y = 1
!         r_ac = clConc_C / clConc
!      ElSEIF (X == 3) THEN
!         Y = 2
!         r_ac = brConc_A / brConc
!      ELSEIF (X == 4) THEN
!         Y = 2
!         r_ac = brConc_C / brConc
!      ELSEIF (X == 5) THEN
!         Y = 1
!         r_ac = clConc_g / clConc
!      ELSEIF (X == 6) THEN
!         Y = 2
!         r_ac = brConc_g / brConc
!      ELSEIF (X == 7) THEN
!         Y = 3
!         r_ac = 1.0_dp
!      ELSEIF (X == 8) THEN
!         Y = 4
!         r_ac = 1.0_dp
!      ENDIF
!
!      CALL Gamma_HOBr_CLD( rLiq,      denAir,   Y,        TK,                &
!                           clConc,    brConc,   hso3Conc, so3Conc,           &
!                           hConc_LCl, GAM_HOBr, r_gp,     H                 )
!      X1 = GAM_HOBr
!      r1 = r_gp * r_ac
!
!      CALL Gamma_HOBr_ICE( Y,       TK,      hno3_th, hcl_th,                &
!                           hbr_th, GAM_HOBr, r_gp                           )
!      X2 = GAM_HOBr
!      IF ( (X >= 5) .AND. (X<=6) ) THEN
!         r2 = r_gp
!      ELSE
!         r2 = 0.0_dp
!      ENDIF
!
!      kISum = CloudHet( CldFr,           Aliq, Aice, rLiq, rIce,             &
!                        SR_MW(ind_HOBr), X1,   X2,   r1,   r2               )
!
    END FUNCTION HetHOBr_TCld
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gamma_HX_Uptake
!
! !DESCRIPTION: Function GAMMA\_HX\_uptake calculates mass accomidation coef.
!               for uptake of HX (HCl or HBr)
!\\
!\\
! !INTERFACE:
!
      FUNCTION Gamma_HX_Uptake( Radius, n_air, X, T )  RESULT( GAM )
!
! !OUTPUT PARAMETER:
      ! Reactive uptake coefficient (unitless)
      REAL(dp)                         :: GAM
! !INPUT PARAMETERS:
!
      ! Radius (cm), n_air (#/cm), and X (1 for Cl and 2 for Br)
      REAL(dp), INTENT(IN)             :: Radius, n_air
      INTEGER, INTENT(IN)              :: X
      REAL(dp), INTENT(IN)             :: T
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!
      REAL(dp)       :: ab

!      ! 1: Cl-, 2: Br-
!      IF (X==1) THEN
!         ! This would never be used since HCl uptake is
!         ! handled by ISORROPIA now, xnw 1/25/18
!         ab = 4.4e-6_dp * EXP( 2898.0e0_dp / T ) ! ab(RT) = 0.069
!      ELSE
!         ab = 1.3e-8_dp * EXP( 4290.0e0_dp / T ) ! ab(RT) = 0.021
!      ENDIF
!
!      GAM = ab

      END FUNCTION Gamma_HX_Uptake
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GAMMA_HOBr_CLD
!
! !DESCRIPTION: Returns GAMMA for HOBr in clouds (need better description)
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE GAMMA_HOBr_CLD( Radius, n_air, X,    T,        C_Y1, C_Y2,    &
                               C_Y3,   C_Y4,  C_Hp, GAM_HOBr, r_gp, H       )
!
! !INPUT PARAMETERS:
!
      ! Radius (cm), n_air (#/cm), and X (1 for Cl and 2 for Br)
      REAL(dp),       INTENT(IN)  :: Radius   ! Radius (cm)
      REAL(dp),       INTENT(IN)  :: n_air    ! n_air (#/cm)
      INTEGER,        INTENT(IN)  :: X        ! 1=Cl-,2=Br-,3=HSO3-,4=SO3--
      REAL(dp),       INTENT(IN)  :: T        ! Temperature (K)
      REAL(dp),       INTENT(IN)  :: C_Y1     ! Concentration (mol/L)
      REAL(dp),       INTENT(IN)  :: C_Y2     ! Concentration (mol/L)
      REAL(dp),       INTENT(IN)  :: C_Y3     ! Concentration (mol/L)
      REAL(dp),       INTENT(IN)  :: C_Y4     ! Concentration (mol/L)
      REAL(dp),       INTENT(IN)  :: C_Hp     ! Concentration (mol/L)
      TYPE(HetState), POINTER     :: H        ! Hetchem species metadata
!
! !OUTPUT PARAMETERS:
!
      ! Reactive uptake coefficient (unitless)
      REAL(dp),       INTENT(OUT) :: GAM_HOBr
      REAL(dp),       INTENT(OUT) :: r_gp
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      REAL(dp) :: ab,     gd,   M_X,  cavg,  H_X
      REAL(dp) :: gb_tot, k_b1, k_b2, k_b3,  k_b4
      REAL(dp) :: D_l,    ybr2, l_r,  k_tot, C_Hp1
      REAL(dp) :: H_HOBr, C_Hp2

!      ! Reaction rate coefficient for HOBr + Cl- [M-2 s-1]
!      ! (Liu and Margerum, Environ. Sci. Tech., 2001)
!      k_b1   = 2.3e+10_dp ! (qjc, 12/28/16)
!
!      ! Reaction rate coefficient for HOBr + Br- [M-2 s-1]
!      k_b2   = 1.6e+10_dp
!
!      ! Reaction rate coefficient for HOBr + HSO3- [M-2 s-1]
!      ! (Liu and Abbatt, Geophys. Res. Let., 2020)
!      k_b3   = 2.6e+7_dp
!
!      ! Reaction rate coefficient for HOBr + HSO3-- [M-2 s-1]
!      ! (Troy and Margerum, Inorg. Chem., 1991)
!      k_b4   = 5.0e+9_dp
!
!      ! Liquid phase diffusion coefficient [cm2/s] for HOBr
!      ! (Ammann et al., Atmos. Chem. Phys., 2013)
!      D_l    = 1.4e-5_dp
!
!      ! Henry's law
!      H_HOBr = H%HOBr%K0 * con_atm_bar
!      H_X    = H_HOBr * EXP( H%HOBr%CR *( 1.0e0_dp/T - 1.0e0_dp/H%HOBr%TK ) )
!      M_X    = H%HOBr%MW_g * 1e-3_dp
!
!      ! Mass accommodation coefficient
!      ab     = 0.6e0_dp
!
!      ! Thermal velocity [cm/s]
!      cavg   = SQRT(8*RStarG*T/(pi*M_X)) *1.0e2_dp
!
!      ! Follow Roberts et al, (2014)
!      C_Hp1  = min(C_Hp, 1.0e-6)
!      C_Hp2  = min(C_Hp, 1.0e-2)
!      C_Hp1  = max(C_Hp1, 1.0e-9)
!      C_Hp2  = max(C_Hp2, 1.0e-6)
!
!      k_tot  = k_b1 * C_Y1 * C_Hp1                                            &
!             + k_b2 * C_Y2 * C_Hp2                                            &
!             + k_b3 * C_Y3                                                    &
!             + k_b4 * C_Y4
!
!      ! l_r is diffusive length scale [cm];
!      !gb is Bulk reaction coefficient [unitless]
!      l_r    = SQRT( D_l / k_tot )
!      gb_tot = 4.0e0_dp * H_X * con_R * T * l_r * k_tot / cavg
!      gb_tot = gb_tot * REACTODIFF_CORR( Radius, l_r)
!
!      ! Reactive uptake coefficient [unitless]
!      GAM_HOBr = 1.0e0_dp / (1.0e0_dp/ab  +  1.0e0_dp/gb_tot)
!
!      ybr2     = 0.41e0*LOG10(C_Y2/C_Y1)+2.25        ! yield of Br2
!      ybr2     = MIN(ybr2, 0.9e0)
!      ybr2      = MAX(ybr2, TINY(1.0e0))
!
!      IF ( X==1 ) THEN
!
!         r_gp = (k_b1 * C_Y1 * C_Hp1 + k_b2 * C_Y2 * C_Hp2) / k_tot
!
!         IF (C_Y2/C_Y1>5.e-4) THEN
!            r_gp = 0.1e0 * r_gp
!         ELSE
!            r_gp = r_gp * (1.e0 - ybr2)
!         ENDIF
!
!      ELSEIF ( X==2 ) THEN
!
!         r_gp = (k_b1 * C_Y1 * C_Hp1 + k_b2 * C_Y2 * C_Hp2) / k_tot
!
!         IF (C_Y2/C_Y1>5.e-4) THEN
!            r_gp = 0.9e0 * r_gp
!         ELSE
!            r_gp = r_gp * ybr2
!         ENDIF
!
!      ELSEIF ( X==3 ) THEN
!
!         r_gp = (k_b3 * C_Y3) / k_tot
!
!      ELSEIF ( X==4 ) THEN
!
!         r_gp = (k_b4 * C_Y4) / k_tot
!
!      ENDIF

    END SUBROUTINE GAMMA_HOBr_CLD
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gamma_HOBr_Ice
!
! !DESCRIPTION: Function GAMMA\_HOBr\_ICE calculates reactive uptake coef.
!               for HOBr + HCl/HBr in ice clouds
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE GAMMA_HOBr_ICE( X, T, hno3_th, hcl_th, hbr_th, GAM_HOBr, r_gp)

! !INPUT PARAMETERS:
!
      INTEGER,  INTENT(IN)  :: X        ! 1=HCl,2=HBr
      REAL(dp), INTENT(IN)  :: T        ! Temperature (K)
      REAL(dp), INTENT(IN)  :: hno3_th
      REAL(dp), INTENT(IN)  :: hcl_th
      REAL(dp), INTENT(IN)  :: hbr_th
!
! !OUTPUT PARAMETERS:
!
      ! Reactive uptake coefficient (unitless)
      REAL(dp), INTENT(OUT) :: GAM_HOBr
      REAL(dp), INTENT(OUT) :: r_gp
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     REAL(dp)    :: rgs, g1, g2

!     rgs  = 0.0e0_dp
!     g1   = 0.0e0_dp
!     g2   = 0.0e0_dp
!
!     !HOBr + HCl
!     rgs = 0.25
!     g1 = rgs * hcl_th
!
!     !HOBr + HBr
!     rgs = 4.8e-4 * exp(1240_dp/T)
!     g2 = rgs * hbr_th
!
!     GAM_HOBr = g1 + g2
!
!     IF ( GAM_HOBr == 0) THEN
!        r_gp=0.0_dp
!     ELSEIF ( X==1 ) THEN
!         r_gp = g1 / GAM_HOBr
!     ELSEIF ( X==2 ) THEN
!         r_gp = g2 / GAM_HOBr
!     ELSE
!         r_gp = 0.0_dp
!     ENDIF


    END SUBROUTINE GAMMA_HOBr_ICE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
    SUBROUTINE GAMMA_HOBr_AER( Radius, n_air, X,        T,    C_Y1,           &
                               C_Y2,   C_Hp,  GAM_HOBr, r_gp, H              )
!
! !INPUT PARAMETERS:
!
      ! Radius (cm), n_air (#/cm), and X (1 for Cl and 2 for Br)
      REAL(dp),       INTENT(IN)  :: Radius   ! Radius (cm)
      REAL(dp),       INTENT(IN)  :: n_air    ! n_air (#/cm)
      INTEGER,        INTENT(IN)  :: X        ! 1=Cl-,2=Br-,3=HSO3-,4=SO3--
      REAL(dp),       INTENT(IN)  :: T        ! Temperature (K)
      REAL(dp),       INTENT(IN)  :: C_Y1     ! Concentration (mol/L)
      REAL(dp),       INTENT(IN)  :: C_Y2     ! Concentration (mol/L)
      REAL(dp),       INTENT(IN)  :: C_Hp     ! Concentration (mol/L)
      TYPE(HetState), POINTER     :: H        ! Hetchem species metadata
!
! !OUTPUT PARAMETERS
!
      ! Reactive uptake coefficient (unitless)
      REAL(dp),       INTENT(OUT) :: GAM_HOBr
      REAL(dp),       INTENT(OUT) :: r_gp
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      REAL(dp) :: ab, gd, M_X
      REAL(dp) :: cavg, H_X, H_HOBr
      REAL(dp) :: k_tot, gb_tot
      REAL(dp) :: k_b1, k_b2
      REAL(dp) :: D_l, ybr2
      REAL(dp) :: l_r, C_Hp1, C_Hp2

!      ! Reaction rate coefficient for HOBr + Cl- [M-2 s-1]
!      !k_b  = 5.9e+9_dp
!      ! (Liu and Margerum, Environ. Sci. Tech., 2001)
!      k_b1  = 2.3e+10_dp ! (qjc, 12/28/16)
!      ! Reaction rate coefficient for HOBr + Br- [M-2 s-1]
!      k_b2  = 1.6e+10_dp
!
!      ! Liquid phase diffusion coefficient [cm2/s] for HOBr
!      ! (Ammann et al., Atmos. Chem. Phys., 2013)
!      D_l    = 1.4e-5_dp
!
!      ! Henry's law
!      H_HOBr = H%HOBr%K0 * con_atm_bar
!      H_X    = H_HOBr*EXP(H%HOBr%CR*(1.0e0_dp/T - 1.0e0_dp/H%HOBr%TK))
!
!      ! Molwt of HOBr (kg/mol)
!      M_X    = H%HOBr%MW_g * 1e-3_dp
!
!      ! Mass accommodation coefficient
!      ab     = 0.6e0_dp
!
!      ! Thermal velocity [cm/s]
!      cavg   = SQRT( 8 * RStarG * T / ( pi * M_X) ) * 1.0e2_dp
!
!      ! Follow Roberts et al, (2014)
!      C_Hp1  = min(C_Hp, 1.0e-6)
!      C_Hp2  = min(C_Hp, 1.0e-2)
!      C_Hp1  = max(C_Hp1, 1.0e-9)
!      C_Hp2  = max(C_Hp2, 1.0e-6)
!
!      k_tot = k_b1 * C_Y1 * C_Hp1                                            &
!            + k_b2 * C_Y2 * C_Hp2
!
!      ! l_r is diffusive length scale [cm];
!      ! gb is Bulk reaction coefficient [unitless]
!      l_r    = SQRT( D_l / k_tot )
!      gb_tot = 4.0_dp * H_X * con_R * T * l_r * k_tot / cavg
!      gb_tot = gb_tot * REACTODIFF_CORR( Radius, l_r)
!
!
!      ! Reactive uptake coefficient [unitless]
!      GAM_HOBr = 1.0_dp / ( 1.0_dp/ab  +  1.0_dp/gb_tot )
!
!      ybr2 = 0.41e0*LOG10(C_Y2/C_Y1)+2.25        ! yield of Br2
!      ybr2 = MIN(ybr2, 0.9e0)
!      ybr2 = MAX(ybr2, TINY(1.0e0))
!
!      IF ( X==1 ) THEN
!
!         IF (C_Y2/C_Y1>5.e-4) THEN
!            r_gp = 0.1e0
!         ELSE
!            r_gp = 1.e0 - ybr2
!         ENDIF
!
!      ELSEIF ( X==2 ) THEN
!
!         IF (C_Y2/C_Y1>5.e-4) THEN
!            r_gp = 0.9e0
!         ELSE
!            r_gp = ybr2
!         ENDIF
!
!      ENDIF

    END SUBROUTINE GAMMA_HOBr_AER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gamma_ClNO3_AER
!
! !DESCRIPTION: Calculates reactive uptake coefficients for ClNO3 + Cl-
!  or ClNO3 + Br-.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Gamma_ClNO3_Aer( Radius, X, M, C_Y1, C_Y2, gamma, r_gp )
!
! !USES:
!
    USE PhysConstants, ONLY : Pi, RStarG
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN)  :: Radius   ! Radius (cm)
    INTEGER,  INTENT(IN)  :: X        ! 1=Cl-,  2=Br-,     3=hydrolysis
    INTEGER,  INTENT(IN)  :: M        ! 1=fine, 2= coarse, 3=cloud
    REAL(dp), INTENT(IN)  :: C_Y1     ! Concentration (mol/L)
    REAL(dp), INTENT(IN)  :: C_Y2     ! Concentration (mol/L)
!
! !OUTPUT PARAMETER:
!
    ! Reactive uptake coefficient (unitless)
    REAL(dp), INTENT(OUT) :: gamma
    REAl(dp), INTENT(OUT) :: r_gp
!EOP
!------------------------------------------------------------------------------
!BOC
!
! ! DEFINED PARAMETERS:
!
    REAL(dp), PARAMETER :: inv_ab = 1.0_dp / 0.108_dp   ! 1 / mass accum coeff
    REAL(dp), PARAMETER :: k_0    = 1.2e+5_dp ** 2.0_dp ! H2k0
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL(dp) :: M_X, fCl, cavg,   D_l, k_2, k_tot
    REAL(dp) :: gb1, gb2, gb_tot, gb0, gbr

    ! Mol wt of ClNO3 (kg/mol)
    M_X = MW(ind_ClNO3) * 1.0e-3_dp

    !-----------------------------------------------------------------------
    ! Calculate gb1 for ClNO3 + Cl-
    !-----------------------------------------------------------------------

    ! Following [Deiber et al., 2004], gamma is not significantly different
    ! from ClNO3 + H2O (gamma = 0.0244) independent of Cl- concentration,
    ! but Cl2 rather than HOCl formed. gb2 can be calculated reversely from
    ! gamma and ab:
    gb1 = 0.0_dp

    ! hydrolysis
    gb0 = 0.032_dp

    !-----------------------------------------------------------------------
    ! Calculate gb2 for ClNO3 + Br-
    !-----------------------------------------------------------------------

    ! thermal velocity (cm/s)
    cavg = SQRT( 8.0_dp * RStarG * TEMP / ( Pi * M_X ) ) * 1.0e+2_dp

    !cm2 s-1
    D_l  = 5.0e-6_dp

    ! H*sqrt(kb)=10^6 (M/s)^½ s-1
    gb2   = 4.0_dp * con_R * TEMP * 1.0e+6_dp * SQRT( C_Y2 * D_l ) / cavg

    ! H2k2br cm2 s-1.
    k_2 = 1.0e+12_dp * C_Y2

    !-----------------------------------------------------------------------
    ! Calculate gb1 for ClNO3 + Cl-
    ! Following [Deiber et al., 2004], gamma is not significantly different
    ! from ClNO3 + H2O (gamma = 0.0244) independent of Cl- concentration,
    ! but Cl2 rather than HOCl formed. gb2 can be calculated reversely from
    ! gb1 = gb0 hydrolysis
    !-----------------------------------------------------------------------
    gb0    = 4.0_dp * con_R * TEMP * 1.2e+5_dp * SQRT( D_l ) / cavg
    k_tot  = k_0 + k_2   !H2(k0+k2Br)
    gb_tot = 4.0_dp * con_R * TEMP * SQRT( k_tot * D_l ) /cavg
    gbr    = k_2 / k_tot

    ! Reaction probability [1]
    gamma = 1.0_dp / ( inv_ab  +  1.0_dp/gb_tot )

    ! 1=fine, 2=coarse, 3=cloud
    SELECT CASE( M )
       CASE( 1 )
          fCl = C(ind_SALACL) / ( C(ind_SALACL) + C(ind_NIT) + C(ind_SO4) )
       CASE( 2 )
          fCl = 1.0_dp
       CASE DEFAULT
          fCl = 0.0_dp
    END SELECT

    ! 1=Cl-, 2=Br-, 3=hydrolysis
    SELECT CASE( X )
       CASE( 1 )
          r_gp = ( 1.0_dp - gbr ) * fCl
       CASE( 2 )
          r_gp = gbr
       CASE DEFAULT
          r_gp = ( 1.0_dp - gbr ) * ( 1.0_dp - fCl )
    END SELECT

  END SUBROUTINE Gamma_ClNO3_Aer
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gamma_ClNO3_Ice
!
! !DESCRIPTION: Calculates reactive uptake coef.
!               for ClNO3 + H2O/HCl/HBr in ice clouds
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Gamma_ClNO3_Ice( X, hno3_th, hcl_th, hbr_th, gamma, r_gp )
!
! !USES:
!
    USE PhysConstants, ONLY : Pi, RStarG
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: X          ! 1=HCl, 2=HBr, 3=hydrolysis
    REAL(dp), INTENT(IN)  :: hno3_th
    REAL(dp), INTENT(IN)  :: hcl_th
    REAL(dp), INTENT(IN)  :: hbr_th
!
! !OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(OUT) :: gamma     ! Reactive uptake coefficient [1]
    REAL(dp), INTENT(OUT) :: r_gp      ! ??
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED_PARAMETERS
!
    REAL(dp), PARAMETER :: twenty = 1.0_dp / 0.5_dp
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL(dp) :: rgs, H2Os, cavg, kks, g1, g2, g3

    ! Initialize
    gamma = 0.0_dp
    r_gp  = 0.0_dp

    ! ClNO3 + HCl
    rgs   = 0.24_dp
    g1    = rgs * hcl_th

    ! ClNO3 + HBr
    rgs   = 0.56_dp
    g2    = rgs * hbr_th

    ! ClNO3 + H2O
    ! cavg = thermal velocity (cm/s)
    cavg  = SQRT( 8.0_dp * RSTARG * TEMP / ( PI * 9.745e-2_dp ) ) * 100.0_dp
    H2Os  = 1e+15_dp - ( 3.0_dp * 2.7e+14_dp * hno3_th )
    kks   = 4.0_dp * 5.2e-17_dp * exp( 2032.0_dp / TEMP )
    g3    = 1.0_dp / ( twenty + cavg / ( kks* H2Os ) )   ! 1.0/0.5 = 20
    gamma = g1 + g2 + g3

    ! Return the proper value
    SELECT CASE( X )
       CASE( 1 )
          r_gp = g1 / gamma     ! HCl-
       CASE( 2 )
          r_gp = g2 / gamma     ! HBr-
       CASE( 3 )
          r_gp = g3 / gamma     ! Hydrolysis
       CASE DEFAULT
          r_gp = 0.0_dp
    END SELECT

  END SUBROUTINE Gamma_ClNO3_Ice
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetClNO3_HCl
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for ClNO3(g) + HCl(l,s)
! in polar stratospheric clouds and on tropospheric sulfate.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETClNO3_HCl( A,B ) RESULT( kISum )

!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(dp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(dp)             :: kISum
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(dp) :: XSTKCF, ADJUSTEDRATE
      Real(dp) :: XSQM

      ! Initialize
      kISum          = 0.0_dp
      ADJUSTEDRATE   = 0.0_dp
      XSTKCF         = 0.0_dp
      XSQM=SQRT(A)

!      ! Loop over aerosol types
!      DO N = 1, NAEROTYPE
!
!         ! Assume zero unless proven otherwise
!         XSTKCF = 0e+0_dp
!
!!         IF (N.eq.8) THEN
!!            XSTKCF = 0.1e-4_dp ! Sulfate
!!         ELSEIF ( STRATBOX ) THEN
!!            IF (N.eq.13) THEN
!	 ! restore limitation to stratosphere - TMS 17/04/10
!         IF  ( STRATBOX ) THEN
!            IF (N.eq.8) THEN
!               XSTKCF = 0.1e-4_dp ! Sulfate
!            ELSEIF (N.eq.13) THEN
!               XSTKCF = KHETI_SLA(4)
!            ELSEIF (N.eq.14) THEN
!               IF (NATSURFACE) THEN
!                  XSTKCF = 0.2e+0_dp ! NAT
!               ELSE
!                  XSTKCF = 0.3e+0_dp ! Ice
!               ENDIF
!            ENDIF
!         ENDIF
!
!         IF (STRATBOX.and.(N.eq.13)) THEN
!            ! Calculate for stratospheric liquid aerosol
!            ! Note that XSTKCF is actually a premultiplying
!            ! factor in this case, including c-bar
!            ADJUSTEDRATE = XAREA(N) * XSTKCF
!         ELSE IF (XSTKCF .GT. 0.0e+0_dp) THEN
!            ! Reaction rate for surface of aerosol
!            ADJUSTEDRATE=ArsL1k(XAREA(N),XRADI(N),NUMDEN,XSTKCF,SR_TEMP, &
!                               XSQM)
!         ENDIF
!
!         ! Add to overall reaction rate
!         kISum = kISum + ADJUSTEDRATE
!
!      END DO

    END FUNCTION HETClNO3_HCl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetBrNO3_HCl
!
! !DESCRIPTION: Set heterogenous chemistry rate for BrNO3(g) + HCl(l,s)
!  in polar stratospheric clouds and on tropospheric sulfate.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETBrNO3_HCl( A, B ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(dp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(dp)             :: kISum
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(dp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      kISum         = 0.0_dp
      ADJUSTEDRATE  = 0.0_dp
      XSTKCF        = 0.0_dp

!      ! Loop over aerosol types
!      DO N = 1, NAEROTYPE
!
!         ! Default to zero
!         XSTKCF = 0.0e+0_dp
!
!	 ! restore limitation to stratosphere - TMS 17/04/10
!	 IF ( STRATBOX ) THEN
!	    IF (N.eq.8) THEN
!               XSTKCF = 0.9_dp ! Sulfate
!            ELSEIF (N.eq.13) THEN
!               XSTKCF = KHETI_SLA(7)
!            ELSEIF (N.eq.14) THEN
!               IF (NATSURFACE) THEN
!                  XSTKCF = 0.3_dp ! NAT   !%%% what is the difference?
!               ELSE
!                  XSTKCF = 0.3_dp ! Ice
!               ENDIF
!            ENDIF
!         ENDIF
!
!         IF (XStkCf.gt.0.0e+0_dp) THEN
!            IF (N.eq.13) THEN
!               ! Calculate for stratospheric liquid aerosol
!               ! Note that XSTKCF is actually a premultiplying
!               ! factor in this case, including c-bar
!               ADJUSTEDRATE = XAREA(N) * XSTKCF
!            ELSE
!               ! Reaction rate for surface of aerosol
!               ADJUSTEDRATE=ArsL1k(XAREA(N),XRADI(N),NUMDEN,XSTKCF,SR_TEMP, &
!                                  (A**0.5_DP))
!            ENDIF
!
!            ! Add to overall reaction rate
!            kISum = kISum + ADJUSTEDRATE
!         ENDIF
!
!      END DO

    END FUNCTION HETBrNO3_HCl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOCl_HCl
!
! !DESCRIPTION: Set heterogenous chemistry rate for HOCl(g) + HCl(l,s)
!  in polar stratospheric clouds and on sulfate aerosol.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOCl_HCl( A, B, Input_Opt ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(dp),       INTENT(IN) :: A, B       ! Rate coefficients
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !RETURN VALUE:
!
      REAL(dp)                   :: kISum
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(dp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      kISum         = 0.0_dp
      ADJUSTEDRATE  = 0.0_dp
      XSTKCF        = 0.0_dp

!      ! Loop over aerosol types
!      DO N = 1, NAEROTYPE
!
!         XSTKCF        = 0.0_dp
!
!         ! For UCX-based mechanisms only consider PSC reactions in strat
!         ! restore limitation to stratosphere - TMS 17/04/10
!         IF ( Input_Opt%LUCX .and. STRATBOX) THEN
!	    IF (N.eq.8) THEN
!	       XSTKCF = 0.8e+0_dp ! Sulfate
!            ELSEIF (N.eq.13) THEN
!               XSTKCF = KHETI_SLA(8)
!            ELSEIF (N.eq.14) THEN
!               IF (NATSURFACE) THEN
!                  XSTKCF = 0.1e+0_dp ! NAT
!               ELSE
!                  XSTKCF = 0.2e+0_dp ! Ice
!               ENDIF
!            ENDIF
!         ENDIF
!
!         IF (XStkCf.gt.0.0e+0_dp) THEN
!            IF (N.eq.13) THEN
!               ! Calculate for stratospheric liquid aerosol
!               ! Note that XSTKCF is actually a premultiplying
!               ! factor in this case, including c-bar
!               ADJUSTEDRATE = XAREA(N) * XSTKCF
!            ELSE
!               ! Reaction rate for surface of aerosol
!               ADJUSTEDRATE=ArsL1k(XAREA(N),XRADI(N),NUMDEN,XSTKCF,SR_TEMP, &
!                                  (A**0.5_DP))
!            ENDIF
!
!            ! Add to overall reaction rate
!            kISum = kISum + AdjustedRate
!
!         ENDIF
!
!      END DO

    END FUNCTION HETHOCl_HCl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOCl_HBr
!
! !DESCRIPTION: Set heterogenous chemistry rate for HOCl(g) + HBr(l,s)
!  in polar stratospheric clouds and on sulfate aerosol.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOCl_HBr( A, B, Input_Opt ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(dp),       INTENT(IN) :: A, B       ! Rate coefficients
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !RETURN VALUE:
!
      REAL(dp)                   :: kISum
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(dp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      kISum         = 0.0_dp
      ADJUSTEDRATE  = 0.0_dp
      XSTKCF        = 0.0_dp

!      ! Loop over aerosol types
!      DO N = 1, H%nAeroType
!
!         XSTKCF        = 0.0_dp
!
!         ! For UCX-based mechanisms only consider PSC reactions in strat
!         ! restore limitation to stratosphere - TMS 17/04/10
!         IF ( Input_Opt%LUCX .and. STRATBOX ) THEN
!	    IF (N.eq.8) THEN
! 	       XSTKCF = 0.8e+0_dp ! Sulfate
!            ELSEIF (N.eq.13) THEN
!               XSTKCF = KHETI_SLA(9)
!            ELSEIF (N.eq.14) THEN
!               IF (NATSURFACE) THEN
!                  XSTKCF = 0.3e+0_dp ! NAT
!               ELSE
!                  XSTKCF = 0.3e+0_dp ! Ice
!               ENDIF
!            ENDIF
!         ENDIF
!
!         IF (XStkCf.gt.0.0e+0_dp) THEN
!            IF (N.eq.13) THEN
!               ! Calculate for stratospheric liquid aerosol
!               ! Note that XSTKCF is actually a premultiplying
!               ! factor in this case, including c-bar
!               ADJUSTEDRATE = XAREA(N) * XSTKCF
!            ELSE
!               ! Reaction rate for surface of aerosol
!               ADJUSTEDRATE=ArsL1k(XAREA(N),XRADI(N),NUMDEN,XSTKCF,SR_TEMP, &
!                                  (A**0.5_DP))
!            ENDIF
!
!            ! Add to overall reaction rate
!            kISum = kISum + AdjustedRate
!
!         ENDIF
!
!      END DO

    END FUNCTION HETHOCl_HBr
!!EOC
!!------------------------------------------------------------------------------
!!                  GEOS-Chem Global Chemical Transport Model                  !
!!------------------------------------------------------------------------------
!!BOP
!!
!! !IROUTINE: HetHOCl_TCld
!!
!! !DESCRIPTION: Sets the rate of the multiphase reaction HOCl + Cl-/S(IV) in
!!  troposphere cloud
!!\\
!!\\
!! !INTERFACE:
!!
    FUNCTION HETHOCl_TCld( denAir,    rLiq,      rIce,      ALiq,            &
                           AIce,      VAir,      TK,        CldFr,           &
                           hConc_Sul, hConc_LCl, hConc_ICl, clConc_A,        &
                           clConc_C,  clConc_g,  hso3Conc,  so3Conc,         &
                           hcl_th,    X,         H                         ) &
                         RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(dp),       INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(dp),       INTENT(IN) :: rLiq        ! Radius of liquid cloud
                                                !  droplets (cm)
      REAL(dp),       INTENT(IN) :: rIce        ! Radius of ice cloud
                                                !  crystals (cm)
      REAL(dp),       INTENT(IN) :: ALiq        ! Area of liquid cloud
                                                !  droplets (cm2/cm3)
      REAL(dp),       INTENT(IN) :: AIce        ! Area of ice cloud
                                                !  crystals (cm2/cm3)
      REAL(dp),       INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(dp),       INTENT(IN) :: TK          ! Temperature (K)
      REAL(dp),       INTENT(IN) :: CldFr       ! Updated to new cloud process
      REAL(dp),       INTENT(IN) :: hConc_Sul   ! Sulfate H+ concentration
      REAL(dp),       INTENT(IN) :: hConc_LCl   ! Liquid cloud H+ concentration
      REAL(dp),       INTENT(IN) :: hConc_ICl   ! Ice cloud H+ concentration
      REAL(dp),       INTENT(IN) :: clConc_A    ! Fine Chloride
                                                !  concentration (mol/L)
      REAL(dp),       INTENT(IN) :: clConc_C    ! Coarse Chloride
                                                !  concentration (mol/L)
      REAL(dp),       INTENT(IN) :: clConc_g
      REAL(dp),       INTENT(IN) :: hso3Conc    ! HSO3-  concentration (mol/L)
      REAL(dp),       INTENT(IN) :: so3Conc     ! SO3--  concentration (mol/L)
      REAL(dp),       INTENT(IN) :: hcl_th
      INTEGER,        INTENT(IN) :: X           ! 1=fineCl;2=coarseCl;
                                                ! 3=HCl;4=HSO3;5=SO3
      TYPE(HetState), POINTER    :: H           ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(dp)                   :: kISum       ! Reaction rate
!!EOP
!!------------------------------------------------------------------------------
!!BOC
!!
!! !LOCAL VARIABLES:
!!
!      INTEGER  :: N,      Y
!      REAL(dp) :: X1,     X2,   GAM_HOCl, r_gp
!      REAL(dp) :: clConc, r_ac, r1,       r2,   rgs
!
!      ! Initialize return value
!      kISum = 0.0_dp
!
!      ! Return if we are in the stratosphere
!      ! This skips unneccessary computations & function calls
!      IF ( StratBox ) RETURN
!
!      ! Continue initializing
!      X1    = 0.0_dp
!      X2    = 0.0_dp
!
!      ! Cloud halide concentration, put fine and coarse mode together
!      clConc = clConc_A + clConc_C + clConc_g
!
!      If (X == 1) THEN
!         Y = 1
!         r_ac = clConc_A / clConc
!      ELSEIF (X == 2) THEN
!         Y = 1
!         r_ac = clConc_C / clConc
!      ElSEIF (X == 3) THEN
!         Y = 1
!         r_ac = clConc_g / clConc
!      ELSEIF (X == 4) THEN
!         Y = 2
!         r_ac = hso3Conc / (hso3Conc + so3Conc)
!      ELSEIF (X == 5) THEN
!         Y = 2
!         r_ac = so3Conc / (hso3Conc + so3Conc)
!      ENDIF
!
!      CALL Gamma_HOCl_CLD( rLiq,     denAir,   Y,       TK,                  &
!                          clConc,   hso3Conc, so3Conc, hConc_LCl,           &
!                           GAM_HOCl, r_gp,     H                            )
!      X1 = GAM_HOCl
!      r1 = r_gp*r_ac
!
!      !For uptake on ice crystals
!      rgs = 0.22_dp
!      X2 = rgs * hcl_th
!      IF (X == 3) THEN
!         r2 = 1.0_dp
!      ELSE
!         r2 = 0.0_dp
!      ENDIF
!
!      kISum = CloudHet( CldFr,           Aliq, Aice, rLiq, rIce,             &
!                        SR_MW(ind_HOCl), X1,   X2,   r1,   r2               )
!
    END FUNCTION HetHOCl_TCld
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOCl_SS
!
! !DESCRIPTION: Set heterogenous chemistry rate for HOCl on Cl- aerosols
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOCl_SS( denAir, rAer,  AAer,   alkAer,                      &
                         TK,     hconc, clconc, H       ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(dp),       INTENT(IN) :: denAir    ! Density of air (#/cm3)
      REAL(dp),       INTENT(IN) :: rAer      ! Radius of aerosol (cm)
      REAL(dp),       INTENT(IN) :: AAer      ! Area of aerosol (cm2/cm3)
      REAL(dp),       INTENT(IN) :: alkAer    ! Aerosol alkalinity (eqn)
      REAL(dp),       INTENT(IN) :: TK        ! Temperature (K)
      REAL(dp),       INTENT(IN) :: hconc     ! H+ (M)
      REAL(dp),       INTENT(IN) :: clconc    ! Cl- (M)
      TYPE(HetState), POINTER    :: H         ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(dp)                   :: kISum     ! Rate for HOCl on Cl- aerosol
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(dp) :: XSTKCF, XSQM

!      ! Initialize
!      kISum = 0.0_dp
!      XSQM  = SQRT( H%HOCl%MW_g )
!
!      ! Reaction can only proceed on acidic aerosol
!      IF ( alkAer > 0.05_dp ) THEN
!         XStkCf = 0.0_dp
!      ELSE
!         XStkCf = GAMMA_HOCl_AER( rAer, denAir, TK, hconc, clconc, H )
!      ENDIF
!
!      ! Reaction rate
!      kISum = ArsL1k( AAer, rAer, denAir, XSTKCF, SR_TEMP, XSQM )

    END FUNCTION HETHOCl_SS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cld_Params
!
! !DESCRIPTION: Subroutine CLD\_PARAMS returns ice and liquid cloud
!  parameters based on State\_Met off of cloud particles.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cld_Params( I, J, L, State_Het, State_Met )
!
! !USES:
!
    USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I           ! Longitude (X) index
    INTEGER,        INTENT(IN)    :: J           ! Latitude  (Y) index
    INTEGER,        INTENT(IN)    :: L           ! Altitude  (Z) index
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    TYPE(HetState), INTENT(INOUT) :: State_Het   ! Hetchem State object
!
! !REMARKS:
!  References:
!   Heymsfield, A. J., Winker, D., Avery, M., et al. (2014). Relationships
!    between ice water content and volume extinction coefficient from in
!    situ observations for temperatures from 0° to –86°C: implications for
!    spaceborne lidar retrievals. Journal of Applied Meteorology and
!    Climatology, 53(2), 479–505. https://doi.org/10.1175/JAMC-D-13-087.1
!
!   Schmitt, C. G., & Heymsfield, A. J. (2005). Total Surface Area Estimates
!    for Individual Ice Particles and Particle Populations. Journal of Applied
!    Meteorology, 44(4), 467–474. https://doi.org/10.1175/JAM2209.1
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Cloud droplet radius in continental warm clouds [cm]
    REAL(dp), PARAMETER :: CLDR_CONT = 6.0e-4_dp

    ! Cloud droplet radius in marine warm clouds [cm]
    REAL(dp), PARAMETER :: CLDR_MARI = 10.0e-4_dp

    ! Ice cloud droplet radius [cm]
    REAL(dp), PARAMETER :: CLDR_ICE  = 38.5e-4_dp

    ! Density of H2O liquid [kg/cm3]
    REAL(dp), PARAMETER :: DENS_LIQ  = 0.001_dp

    ! Density of H2O ice [kg/cm3]
    REAL(dp), PARAMETER :: DENS_ICE  = 0.91e-3_dp
!
! !LOCAL VARIABLES:
!
    REAL(dp) :: alpha, beta

    !=======================================================================
    ! CLD_PARAMS begins here!
    !=======================================================================

    ! Exit if there is no cloud
    IF ( ( State_Met%QL(I,J,L) + State_Met%QI(I,J,L) <= 0.0_dp )  .or.     &
         ( State_Met%CLDF(I,J,L)                     <= 0.0_dp ) ) THEN
       State_Het%rLiq = CLDR_CONT
       State_Het%rIce = CLDR_ICE
       State_Het%ALiq = 0.0_dp
       State_Het%VLiq = 0.0_dp
       State_Het%AIce = 0.0_dp
       State_Het%VIce = 0.0_dp
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! In GC 12.0 and earlier, the liquid water volume was set to zero at
    ! temperatures colder than 258K and over land ice (Antarctica &
    ! Greenland). That was likely legacy code from GEOS-4, which provided
    ! no information on cloud phase. As of GC 12.0, all met data sources
    ! provide explicit liquid and ice condensate amounts, so we use those
    ! as provided. (C.D. Holmes)
    !
    ! Liquid water clouds
    !
    ! Droplets are spheres, so
    ! Surface area = 3 * Volume / Radius
    !
    ! Surface area density = Surface area / Grid volume
    !-----------------------------------------------------------------------
    IF ( State_Met%FRLAND(I,J) > State_Met%FROCEAN(I,J) ) THEN
       State_Het%rLiq = CLDR_CONT      ! Continental cloud droplet radius [cm]
    ELSE
       State_Het%rLiq = CLDR_MARI      ! Marine cloud droplet radius [cm]

    ENDIF

    ! get the volume of cloud condensate [cm3(condensate)/cm3(air)]
    ! QL is [g/g]
    State_Het%VLiq = State_Met%QL(I,J,L) * State_Met%AD(I,J,L)               &
                   / DENS_LIQ            / State_Het%vAir
    State_Het%VIce = State_Met%QI(I,J,L) * State_Met%AD(I,J,L)               &
                   / DENS_ICE            / State_Het%vAir

    State_Het%ALiq = 3.0_dp * State_Het%vLiq / State_Het%rLiq

    !-----------------------------------------------------------------------
    ! Ice water clouds
    !
    ! Surface area calculation requires information about ice crystal size
    ! and shape, which is a function of temperature. Use Heymsfield (2014)
    ! empirical relationships between temperature, effective radius,
    ! surface area and ice water content.
    !
    ! Schmitt and Heymsfield (2005) found that ice surface area is about
    ! 9 times its cross-sectional area.
    !
    ! For any shape,
    !   Cross section area = pi * (Effective Radius)^2, so
    !   Cross section area = 3 * Volume / ( 4 * Effective Radius ).
    !
    ! Thus, for ice
    !   Surface area = 9 * Cross section area
    !                = 2.25 * 3 * Volume / Effective Radius
    ! (C.D. Holmes)
    !-----------------------------------------------------------------------

    ! Heymsfield (2014) ice size parameters
    IF ( State_Met%T(I,J,L) < 202.0_dp ) THEN          ! -71 C
       alpha = 83.3_dp
       beta  = 0.0184_dp
    ELSE IF ( State_Met%T(I,J,L) < 217.0_dp ) THEN     ! -56 C
       alpha = 9.1744e+4_dp
       beta  = 0.117_dp
    ELSE
       alpha = 308.4_dp
       beta  = 0.0152_dp
    ENDIF

    ! Effective radius, cm
    State_Het%rIce = 0.5_dp * alpha                                          &
                   * EXP(beta * (State_Met%T(I,J,L) - 273.15_dp)) / 1e+4_dp

    ! Ice surface area density, cm2/cm3
    State_Het%aIce = 3.0_dp * State_Het%vIce / State_Het%rIce * 2.25_dp

    !=======================================================================
    ! Get theta for ice cloud uptake
    !=======================================================================
    CALL Get_Theta_Ice( C(ind_HNO3), C(ind_HCl), C(ind_HBr), State_Het )

  END SUBROUTINE Cld_Params
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Theta_Ice
!
! !DESCRIPTION: Subroutine GET_THETA_ICE returns theta values for
!  HNO3, HCl, and HBr for ice uptake calculations
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Theta_Ice( HNO3, HCl, HBr, H )
!
!
! !INPUT PARAMETERS:
!
    REAL(dp),      INTENT(IN)     :: HNO3  ! HNO3 conc [molec/cm3]
    REAL(dp),      INTENT(IN)     :: HCl   ! HCl  conc [molec/cm3]
    REAL(dp),      INTENT(IN)     :: HBr   ! HBr  conc [molec/cm3]
!
! !RETURN VALUE:
!
    TYPE(HetState), INTENT(INOUT) :: H     ! Hetchem State object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(dp) :: KLangC1, KLangC2, KlinC, denom

    !=================================================================
    ! GET_THETA_ICE begins here!
    !=================================================================

    ! HNO3
    KlinC   = 7.5e-5_dp * EXP( 4585.0_dp / TEMP )          ! 1/cm
    KLangC1 = KlinC / 2.7e+14_dp                           ! cm3/molec

    ! HCl
    KlinC   = 1.3e-5_dp * EXP( 4600.0_dp / TEMP )          ! 1/cm
    KLangC2 = KlinC / 3.0e+14_dp                           ! cm3/molec

    denom        = 1.0_dp + KLangC1*HNO3 + KLangC2*HCl
    H%HNO3_theta = KLangC1*HNO3 / denom
    H%HCl_theta  = KLangC2*HCl  / denom

    ! HBr
    H%HBr_theta = 4.14e-10_dp * ( HBr**0.88_dp )
    H%HBr_theta = MIN( H%HBr_theta, 1.0_dp )

  END SUBROUTINE Get_Theta_Ice
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Halide_Conc
!
! !DESCRIPTION: Initializes halide (Br-, Cl-) concentrations at each
!  grid box for use in heterogeneous chemistry routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Halide_Conc( I, J, L, H )
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I    ! Longitude (X) index
    INTEGER,        INTENT(IN)    :: J    ! Latitude  (Y) index
    INTEGER,        INTENT(IN)    :: L    ! Altitude  (Z) index
!
! !OUTPUT PARAMETERS:
!
    TYPE(HetState), INTENT(INOUT) :: H    ! Hetchem State object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(dp) :: Br_conc, Cl_conc, denom, HBr, HCl

    !=======================================================================
    ! Get halide conc's in cloud (gas-phase, fine & coarse sea salt)
    !=======================================================================

    ! Br- and Cl- grid-box concentrations
    HBr = C(ind_HBr) + ( C(ind_BrSALA) * 0.7_dp ) + C(ind_BrSALC)
    HCl = C(ind_HCl) + ( C(ind_SALACL) * 0.7_dp ) + C(ind_SALCCL)

    ! Get overall Br- and Cl- grid box concentrations in cloud
    CALL Get_Halide_CldConc( H, HBr, HCl, Br_conc, Cl_conc )

    ! Split Br- into gas-phase (G), fine sea salt (A), coarse sea salt (C)
    ! Avoid div-by-zero (all three expressions use the same denominator)
    denom = C(ind_HBr) + ( C(ind_BrSALA) * 0.7_dp ) + C(ind_BrSALC)
    IF ( denom > 0.0_dp ) THEN
       H%Br_conc_CldG = ( Br_conc * C(ind_HBr   )          ) / denom
       H%Br_conc_CldA = ( Br_conc * C(ind_BrSALA) * 0.7_dp ) / denom
       H%Br_conc_CldC = ( Br_conc * C(ind_BrSALC)          ) / denom
    ELSE
       H%Br_conc_CldG = 0.0_dp
       H%Br_conc_CldA = 0.0_dp
       H%Br_conc_CldC = 0.0_dp
    ENDIF

    ! Enforce minimum values
    H%Br_conc_CldG = MAX( H%Br_conc_CldG, 1.0e-20_dp )
    H%Br_conc_CldA = MAX( H%Br_conc_CldA, 1.0e-20_dp )
    H%Br_conc_CldC = MAX( H%Br_conc_CldC, 1.0e-20_dp )

    ! Split Cl- into gas-phase (G), fine sea salt (A), coarse sea salt (C)
    ! Avoid div-by-zero (all three expressions use the same denominator)
    denom = C(ind_HCl) + ( C(ind_SALACL) * 0.7_dp ) + C(ind_SALCCL)
    IF ( denom > 0.0_dp ) THEN
       H%Cl_conc_CldG = ( Cl_conc * C(ind_HCl   )          ) / denom
       H%Cl_conc_CldA = ( Cl_conc * C(ind_SALACL) * 0.7_dp ) / denom
       H%Cl_conc_CldC = ( Cl_conc * C(ind_SALCCL)          ) / denom
    ELSE
       H%Cl_conc_CldG = 0.0_dp
       H%Cl_conc_CldA = 0.0_dp
       H%Cl_Conc_CldC = 0.0_dp
    ENDIF

    ! Enforce minimum values
    H%Cl_Conc_CldG = MAX( H%Cl_conc_CldG, 1.0e-20_dp )
    H%Cl_Conc_CldA = MAX( H%Cl_conc_CldA, 1.0e-20_dp )
    H%Cl_Conc_CldC = MAX( H%Cl_conc_CldC, 1.0e-20_dp )

    ! Total Br- and Cl- in cloud
    H%Br_conc_Cld    = H%Br_conc_CldA + H%Br_conc_CldC + H%Br_conc_CldG
    H%Cl_conc_Cld    = H%Cl_conc_CldA + H%Cl_conc_CldC + H%Cl_conc_CldG

    ! Branching ratios for Br- in each of the CldA, CldG, CldC paths
    H%Br_branch_CldA = H%Br_conc_CldA / H%Br_conc_Cld
    H%Br_branch_CldC = H%Br_conc_CldC / H%Br_conc_Cld
    H%Br_branch_CldG = H%Br_conc_CldG / H%Br_conc_Cld

    ! Branching ratios for Br- in each of the CldA, CldG, CldC paths
    H%Cl_branch_CldA = H%Cl_conc_CldA / H%Cl_conc_Cld
    H%Cl_branch_CldC = H%Cl_conc_CldC / H%Cl_conc_Cld
    H%Cl_branch_CldG = H%Cl_conc_CldC / H%Cl_conc_Cld

    ! Ratio of Br- in gas phase in cloud / Cl- in gas-phase in cloud
    H%Br_over_Cl_Cld = H%Br_conc_Cld / H%Cl_conc_Cld

    !=======================================================================
    ! Get halide concentrations, in aerosol
    !=======================================================================

    ! Br- concentration in fine sea salt aerosol
    CALL Get_Halide_SSAConc( n_x       = C(ind_BrSALA),                      &
                             surf_area = H%AClAREA,                          &
                             r_w       = H%AClRADI,                          &
                             conc_x    = H%Br_conc_SALA                     )

    ! Br- concentration in coarse sea salt aerosol
    CALL Get_Halide_SSAConc( n_x       = C(ind_BrSALC),                      &
                             surf_area = H%XAREA(12),                        &
                             r_w       = H%XRADI(12),                        &
                             conc_x    = H%Br_conc_SALC                     )

    ! Cl- concentration in fine sea salt aerosol
    CALL Get_Halide_SSAConc( n_x       = C(ind_SALACL),                      &
                             surf_area = H%AClAREA,                          &
                             r_w       = H%AClRADI,                          &
                             conc_x    = H%Cl_conc_SALA                     )

    ! Cl- concentration in coarse sea salt aerosol
    CALL Get_Halide_SSAConc( n_x       = C(ind_SALCCL),                      &
                             surf_area = H%xArea(12),                        &
                             r_w       = H%xRadi(12),                        &
                             conc_x    = H%Cl_conc_SALC                     )

    ! NO3- concentration in fine sea salt aerosol
    CALL Get_Halide_SSAConc( n_x       = C(ind_NIT),                         &
                             surf_area = H%AClAREA,                          &
                             r_w       = H%AClRADI,                          &
                             conc_x    = H%NIT_conc_SALA                    )

    ! NO3- concentration in coarse sea salt aerosol
    CALL Get_Halide_SSAConc( n_x       = C(ind_NITs),                        &
                             surf_area = H%xArea(12),                        &
                             r_w       = H%xRadi(12),                        &
                             conc_x    = H%NIT_conc_SALC                    )

  END SUBROUTINE Halide_Conc
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Halide_CldConc
!
! !DESCRIPTION: Subroutine GET\_HALIDE\_CLDCONC returns the in-cloud
!  concentration of bromide and chloride (Br- and Cl-).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Halide_CldConc( H, HBr, HCl, br_conc, cl_conc )
!
! !USES:
!
! !INPUT PARAMETERS:
!
    TYPE(HetState), INTENT(IN)  :: H          ! Hetchem State object
    REAL(dp),       INTENT(IN)  :: HBr        ! HBr- concentration [#/cm3]
    REAL(dp),       INTENT(IN)  :: HCl        ! HCl- concentration [#/cm3]
!
! !OUTPUT PARAMETERS:
!
    REAL(dp),       INTENT(OUT) :: Br_conc    ! Br conc [M/kg water]
    REAL(dp),       INTENT(OUT) :: Cl_conc    ! Cl conc [M/kg water]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(dp) :: V_tot, dr_ratio, t2l, L2G, F_L, pH

    !=================================================================
    ! Get_Halide_CldConc begins here!
    !=================================================================

    !---------------------------------------------------------------
    ! jas, 07/30/2014 (SETUP d/r ratio for ice cloud droplets)
    ! V_liq = 4pi/3 ( r^3 - (r - r*(d/r))^3 = (r^3 - r^3*(1 - d/r)^3) = r^3 (1
    ! - (1 - d/r)^3
    ! V_tot / V_liq = 1 / (1 - (1 - d/r)^3))
    DR_RATIO = 2e-2_dp
    T2L      = 1.0_dp / ( 1.0_dp - ( 1.0_dp - DR_RATIO)**3 )
    !---------------------------------------------------------------

    ! V_tot = VLiq + (VIce / T2L) ! (cm3(liq)/cm3(air)
    V_tot = H%VLiq
    V_tot = SAFE_DIV( V_tot, H%CldFr, 0.0_dp )

    ! Exit if not in cloud
    IF ( V_tot < 1.0e-20_dp ) THEN
       Br_conc = 1.0e-20_dp
       Cl_conc = 1.0e-20_dp
       RETURN
    ENDIF

    ! Chloride (mol/L)
    CALL Compute_L2G_Local( 1.0_dp, 9000.0_dp, -6.3_dp,                      &
                            TEMP,   V_tot,      L2G,    pH                  )
    F_L = L2G / ( 1.0_dp + L2G )
    cl_conc = F_L * HCl / (V_tot * AVO * 1.0e-3_dp)

    ! Bromide (mol/L)
    CALL Compute_L2G_Local( 7.5e-1_dp, 10200.0_dp, -9.0_dp,                  &
                            TEMP,      V_tot,       L2G,   pH               )
    F_L = L2G / ( 1.0_dp + L2G )
    br_conc = F_L * HBr / ( V_tot * AVO * 1.0e-3_dp )

  END SUBROUTINE Get_Halide_CldConc
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Halide_SsaConc
!
! !DESCRIPTION: Calculates concentration of a halide in sea salt aerosol.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Get_Halide_SsaConc( n_x, surf_area, r_w, conc_x )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN)  :: n_x         ! Number density     [#/cm3  ]
    REAL(dp), INTENT(IN)  :: surf_area   ! Surface area       [cm2/cm3]
    REAL(dp), INTENT(IN)  :: r_w         ! Aerosol wet radius [cm     ]
!
! !OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(OUT) :: conc_x      ! Halide conc in seasalt [mol/L]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(dp) :: V_tot

    !==================================================================
    ! Get_Halide_SsaConc begins here!
    !==================================================================

    ! Cloud volume
    V_tot = surf_area * r_w * 0.3333333e0_dp * 1e-3_dp ! L(liq)/cm3(air)

    ! Skip if we are not in cloud
    IF ( V_tot <= 1.0e-20_dp ) THEN
       conc_x = 1.0e-20_dp
       RETURN
    ENDIF

    ! This calculation can be used for both SSA X- concentration and for
    ! those out of cloud only. For X- out of cloud only, V_tot =
    ! V_tot*(1-CF), n_x = n_x*(1-CF), so (1-CF) is canceled.
    ! xnw, 02/05/18
    conc_x = ( n_x / AVO ) / V_tot    ! mol/L
    conc_x = MAX( conc_x, 1.0e-20_dp )

  END SUBROUTINE Get_Halide_SsaConc
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute_L2G_Local
!
! !DESCRIPTION: Subroutine COMPUTE\_L2G\_LOCAL is a local copy of the
!  liquid-gas partitioning routine in GEOS-Chem's wetscav\_mod.F file.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_L2G_Local( K0, CR, pKa, TK, H2OLIQ, L2G, pH )
!
! !USES:
!
    USE Henry_Mod, ONLY : Calc_KH
    USE Henry_Mod, ONLY : Calc_Heff
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN)  :: K0     ! Henry's solubility constant [M/atm]
    REAL(dp), INTENT(IN)  :: CR     ! Henry's volatility constant [K]
    REAL(dp), INTENT(IN)  :: pKa    ! Henry's pH correction factor [1]
    REAL(dp), INTENT(IN)  :: TK     ! Temperature [K]
    REAL(dp), INTENT(IN)  :: H2OLIQ ! Liquid water content [cm3 H2O/cm3 air]
    REAL(dp), INTENT(IN)  :: pH     ! Liquid water pH
!
! !OUTPUT PARAMETERS:
!
    REAL(dp), INTENT(OUT) :: L2G    ! Cliq/Cgas ratio [1]
!
! !REMARKS:
!  The ratio Cliq / Cgas is obtained via Henry's law.  The appropriate
!  values of Kstar298 and H298_R must be supplied for each species.
!  (cf Jacob et al 2000, p. 3)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: RC
    REAL*8  :: CR_8,  H2OLIQ_8, HEFF_8, K0_8, KH_8
    REAL*8  :: L2G_8, TK_8,     pKa_8,  pH_8

    !=================================================================
    ! COMPUTE_L2G_LOCAL begins here!
    !=================================================================

    ! Cast inputs to REAL*8
    CR_8  = CR
    K0_8  = K0
    pka_8 = pKa
    pH_8  = pH
    TK_8  = TK

    ! For wetdep, we assume a pH of 4.5 for rainwater
    !pH = 4.5_dp

    ! Calculate the Henry's law constant
    CALL CALC_KH( K0_8, CR_8, TK_8, KH_8, RC )

    ! Calculate effective Henry's law constant, corrected for pH
    ! (for those species that have a defined pKa value)
    CALL CALC_HEFF( pKa_8, pH_8, KH_8, HEFF_8, RC )

    ! Use Henry's Law to get the ratio:
    ! [ mixing ratio in liquid phase / mixing ratio in gas phase ]
    L2G_8 = HEFF_8 * H2OLIQ

    ! Cast outputs to flex-precision
    L2G = L2G_8

  END SUBROUTINE Compute_L2G_Local
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cld1k_XNO3
!
! !DESCRIPTION: Function CLD1K\_XNO3 calculates the rate constant for
!  heterogeneous cycling of XNO3 off of cloud particles.
!\\
!\\
! !INTERFACE:
!
    FUNCTION CLD1K_XNO3( denAir, TK, rLiq, rIce, ALiq, AIce, &
                         MX_gmol, AlphaX )    RESULT( cld1k )

!
! !USES:
!
!
! !INPUT PARAMETERS:
!
      REAL(dp),       Intent(IN) :: DENAIR   ! Density of air [#/cm3]
      REAL(dp),       Intent(In) :: TK       ! Air temperature [K]
      REAL(dp),       Intent(In) :: rLiq     ! Radius of liquid cloud drops [cm]
      REAL(dp),       Intent(In) :: rIce     ! Radius of ice cloud crystals [cm]
      REAL(dp),       Intent(In) :: ALiq     ! Surface area (liquid) [cm2/cm3]
      REAL(dp),       Intent(In) :: AIce     ! Surface area (ice) ) [cm2/cm3]
      REAL(dp),       Intent(IN) :: MX_gmol  ! Molecular mass of XNO3 [g/mol]
      REAL(dp),       Intent(IN) :: AlphaX   ! XNO3 accomodation coef [unitless]
!
! !RETURN VALUE:
!
      REAL(dp)                   :: cld1k    ! Rate constant for heterogeneous
                                             ! cycling of BrNO3 off of cloud
                                             ! particles
!EOP
!------------------------------------------------------------------------------
!BOC
!!
!! !LOCAL VARIABLES:
!!
!      REAL(dp)             :: SQM        ! Square root of molec. weight [g/mol]
!      REAL(dp)             :: STK        ! Square root of temperature   [K]
!      REAL(dp)             :: DFKG       ! Gas diffusion coefficient    [cm2/s]
!
!      !=================================================================
!      ! CLD1K_XNO3 begins here!
!      !=================================================================
!
!      ! Quick test - is there any cloud?
!      IF ((ALiq.le.0.0e+0_dp).and.(AIce.le.0.0e+0_dp)) THEN
!         cld1k = 0.0e+0_dp
!         Return
!      ENDIF
!
!      ! ------------------------------------------------------------
!      !   Calculate the 1st order rate constant for XNO3 hydrolysis.
!      !
!      !   (a) calculate the gas phase diffusion coefficient;
!      !
!      !   (b) calculate the hydrolysis rxn rate.
!      ! ------------------------------------------------------------
!      SQM = sqrt(MX_gmol)    ! square root of molar mass [g/mole]
!      STK = sqrt(TK) ! square root of temperature [K]
!
!      ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
!      DFKG  = 9.45E+17_dp/DENAIR * STK * SQRT(3.472E-2_dp + 1.E+0_dp/(SQM*SQM))
!
!      ! Compute ArsL1k according to the formula listed above
!      ! Sum contribution from ice and liquid clouds
!      cld1k = ALiq / ( rLiq/DFKG + 2.749064E-4 * SQM/(ALPHAX*STK) )
!      cld1k = AIce / ( rIce/DFKG + 2.749064E-4 * SQM/(ALPHAX*STK) ) + cld1k

    END FUNCTION CLD1K_XNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FcrO2HO2
!
! !DESCRIPTION: !fgap, based on saunder 2003 k14
!\\
!\\
! !INTERFACE:
!
    FUNCTION FCRO2HO2( XCARBN ) RESULT( FC_RO2HO2 )
!
! !INPUT PARAMETERS:
!
      REAL(dp), INTENT(IN) :: XCARBN
!
! !RETURN VALUE:
!
      REAL(dp)             :: FC_RO2HO2
!EOP
!------------------------------------------------------------------------------
!BOC

      FC_RO2HO2 = 1E+0_dp - EXP( -0.245E+0_dp * XCARBN )

    END FUNCTION FCRO2HO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FyHORO
!
! !DESCRIPTION: \subsection*{Overview}
!  Function FYHORO returns returns the branching ratio between
!  HOC2H4O oxidation and dissociation:
!  (1) HOC2H4 --O2--> HO2 + GLYC
!  (2) HOC2H4 ------> HO2 + 2CH2O
!
!\subsection*{References}
!  \begin{enumerate}
!  \item Orlando et al., 1998: \emph{Laboratory and theoretical study of the
!         oxyradicals in the OH- and Cl-initiated oxidation of ethene},
!        \underline{J. Phys. Chem. A}, \textbf{102}, 8116-8123.
!  \item Orlando et al., 2003: \emph{The atmospheric chemistry of alkoxy
!         radicals}, \underline{Chem. Rev.}, \textbf{103}, 4657-4689.
!  \end{enumerate}
!
! !INTERFACE:
!
    FUNCTION FYHORO( ZDNUM, TT ) RESULT( FY_HORO )
!
! !INPUT PARAMETERS:
!
      REAL(dp), INTENT(IN) :: ZDNUM   ! Air density   [molec/cm3 ]
      REAL(dp), INTENT(IN) :: TT      ! Temperature   [K         ]
!
! !RETURN VALUE:
!
      REAL(dp)             :: FY_HORO
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(dp) :: K1, K2, O2DNUM

      !=================================================================
      ! FYHORO begins here!
      !=================================================================
      O2DNUM  = ZDNUM * 0.21E+0_dp
      K1      = 6.0E-14_dp * EXP(-550.E+0_dp/TT) * O2DNUM
      K2      = 9.5E+13_dp * EXP(-5988.E+0_dp/TT)

      FY_HORO = K1 / (K1 + K2)

    END FUNCTION FYHORO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: FyrNO3
!
! !DESCRIPTION: Function FYRNO3 returns organic nitrate yields
!  YN = RKA/(RKA+RKB) from RO2+NO reactions as a function of the number
!  N of carbon atoms.
!\\
!\\
! !INTERFACE:
!
  FUNCTION FYRNO3( XCARBN, ZDNUM, TT ) RESULT( FYR_NO3 )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: XCARBN   ! Number of C atoms in RO2
    REAL(dp), INTENT(IN) :: ZDNUM    ! Air density   [molec/cm3 ]
    REAL(dp), INTENT(IN) :: TT       ! Temperature   [K         ]
!
! !RETURN VALUE:
!
    REAL(dp)             :: FYR_NO3
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(dp) :: YYYN, XXYN,  AAA,  RARB, ZZYN
    REAL(dp) :: XF,   ALPHA, Y300, BETA, XMINF, XM0

    ! Initialize static variables
    DATA   Y300,ALPHA,BETA,XM0,XMINF,XF/.826,1.94E-22,.97,0.,8.1,.411/

    !=================================================================
    ! FYRNO3 begins here!
    !=================================================================
    XXYN    = ALPHA*EXP(BETA*XCARBN)*ZDNUM*((300./TT)**XM0)
    YYYN    = Y300*((300./TT)**XMINF)
    AAA     = LOG10(XXYN/YYYN)
    ZZYN    = 1./(1.+ AAA*AAA )
    RARB    = (XXYN/(1.+ (XXYN/YYYN)))*(XF**ZZYN)
    FYR_NO3 = RARB/(1. + RARB)

  END FUNCTION FYRNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: ArsL1k
!
! !DESCRIPTION: Calculates the 1st-order loss rate of species on wet
!  aerosol surface
!\\
!\\
! !INTERFACE:
!
  FUNCTION ArsL1k( area, radius, denAir, gamma, srTk, srMw ) RESULT( rate )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: area   ! Area of wet aerosol / vol of air [cm2/cm3]
    REAL(dp), INTENT(IN) :: radius ! Radius of wet aerosol (Rd) [cm]
    REAL(dp), INTENT(IN) :: denAir ! Density of air [#/cm3]
    REAL(dp), INTENT(IN) :: gamma  ! Reaction probability gamma [1]
    REAL(dp), INTENT(IN) :: srTk   ! Square root of temperature [K]
    REAL(dp), INTENT(IN) :: srMw   ! Square root of mol wt [g/mole]
!
! !RETURN VALUE:
!
    REAL(dp)             :: rate   ! Reaction rate [1/s]
!
! !REMARKS:
!  The 1st-order loss rate on wet aerosol (Dentener's Thesis, p. 14)
!  is computed as:
!                                                                             .
!     ArsL1k [1/s] = area / [ radius/dfkg + 4./(gamma * xmms) ]
!                                                                             .
!  where XMMS = Mean molecular speed [cm/s] = sqrt(8R*TK/pi/M) for Maxwell
!                                                                             .
!     DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(dp) :: dfkg

    !========================================================================
    ! ArsL1k begins here!
    !========================================================================

    ! If gamma or radius is very small, set rate to zero and return
    IF ( gamma < 1.0e-30_dp .or. radius < 1.0e-30_dp ) THEN
       !rate = 1.0e-30_dp
       rate = 0.0_dp
       RETURN
    ENDIF

    ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
    dfkg = ( 9.45E+17_dp / denAir ) * srTk *                                 &
           SQRT( 3.472E-2_dp + 1.0_dp / ( srMw * srMw ) )

    ! Compute ArsL1k according to the formula listed above
    rate = area / ( (radius / dfkg) + 2.749064E-4_dp * srMw / (gamma * srTk) )

  END FUNCTION ArsL1k
!EOC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%              NEW RATE-LAW FUNCTIONS BEGIN HERE                    %%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !=========================================================================
  ! Rate-law functions for bromine species
  !=========================================================================


!EOC
END MODULE GcKpp_HetRates
