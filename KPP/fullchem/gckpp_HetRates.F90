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
!    US. Atmos. Chem. Phys., 16, 2961-1.02990, 2016.
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
    ! Initialization
    !=======================================================================

    ! Set Br and Cl fields of State_Het to zero
    H%Br_conc_CldG   = 0.0_dp
    H%Br_conc_CldA   = 0.0_dp
    H%Br_conc_CldC   = 0.0_dp
    H%Br_over_Cl_Cld = 0.0_dp
    H%Br_over_Cl_SSA = 0.0_dp
    H%Br_over_Cl_SSC = 0.0_dp
    H%Cl_conc_CldG   = 0.0_dp
    H%Cl_conc_CldA   = 0.0_dp
    H%Cl_conc_CldC   = 0.0_dp
    H%frac_Br_CldA   = 0.0_dp
    H%frac_Br_CldC   = 0.0_dp
    H%frac_Br_CldG   = 0.0_dp
    H%frac_Cl_CldA   = 0.0_dp
    H%frac_Cl_CldC   = 0.0_dp
    H%frac_Cl_CldG   = 0.0_dp

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
    ENDIF

    ! Split Cl- into gas-phase (G), fine sea salt (A), coarse sea salt (C)
    ! Avoid div-by-zero (all three expressions use the same denominator)
    denom = C(ind_HCl) + ( C(ind_SALACL) * 0.7_dp ) + C(ind_SALCCL)
    IF ( denom > 0.0_dp ) THEN
       H%Cl_conc_CldG = ( Cl_conc * C(ind_HCl   )          ) / denom
       H%Cl_conc_CldA = ( Cl_conc * C(ind_SALACL) * 0.7_dp ) / denom
       H%Cl_conc_CldC = ( Cl_conc * C(ind_SALCCL)          ) / denom
    ENDIF

    ! Total Br- and Cl- in cloud
    H%Br_conc_Cld = H%Br_conc_CldA + H%Br_conc_CldC + H%Br_conc_CldG
    H%Cl_conc_Cld = H%Cl_conc_CldA + H%Cl_conc_CldC + H%Cl_conc_CldG

    ! Fractions of Br- in each of the CldA, CldG, CldC paths
    IF ( H%Br_Conc_Cld > 0.0_dp ) THEN
       H%frac_Br_CldA = H%Br_conc_CldA / H%Br_conc_Cld
       H%frac_Br_CldC = H%Br_conc_CldC / H%Br_conc_Cld
       H%frac_Br_CldG = H%Br_conc_CldG / H%Br_conc_Cld
    ENDIF

    ! Branching ratios for Br- in each of the CldA, CldG, CldC paths
    IF ( H%Cl_Conc_Cld > 0.0_dp ) THEN
       H%frac_Cl_CldA = H%Cl_conc_CldA / H%Cl_conc_Cld
       H%frac_Cl_CldC = H%Cl_conc_CldC / H%Cl_conc_Cld
       H%frac_Cl_CldG = H%Cl_conc_CldG / H%Cl_conc_Cld
    ENDIF

    !=======================================================================
    ! Get halide concentrations, in aerosol
    !=======================================================================

    ! Br- concentration in fine sea salt aerosol
    CALL Get_Halide_SSAConc( n_x       = C(ind_BrSALA),                      &
                             surf_area = H%aClArea,                          &
                             r_w       = H%aClRadi,                          &
                             conc_x    = H%Br_conc_SSA                      )

    ! Br- concentration in coarse sea salt aerosol
    CALL Get_Halide_SSAConc( n_x       = C(ind_BrSALC),                      &
                             surf_area = H%xArea(12),                        &
                             r_w       = H%xRadi(12),                        &
                             conc_x    = H%Br_conc_SSC                      )

    ! Cl- concentration in fine sea salt aerosol
    CALL Get_Halide_SSAConc( n_x       = C(ind_SALACL),                      &
                             surf_area = H%aClArea,                          &
                             r_w       = H%aClRadi,                          &
                             conc_x    = H%Cl_conc_SSA                      )

    ! Cl- concentration in coarse sea salt aerosol
    CALL Get_Halide_SSAConc( n_x       = C(ind_SALCCL),                      &
                             surf_area = H%xArea(12),                        &
                             r_w       = H%xRadi(12),                        &
                             conc_x    = H%Cl_conc_SSC                      )

    ! NO3- concentration in fine sea salt aerosol
    CALL Get_Halide_SSAConc( n_x       = C(ind_NIT),                         &
                             surf_area = H%aClArea,                          &
                             r_w       = H%aClRadi,                          &
                             conc_x    = H%NIT_conc_SSA                     )

    ! NO3- concentration in coarse sea salt aerosol
    CALL Get_Halide_SSAConc( n_x       = C(ind_NITs),                        &
                             surf_area = H%xArea(12),                        &
                             r_w       = H%xRadi(12),                        &
                             conc_x    = H%NIT_conc_SSC                     )

    !=======================================================================
    ! Ratios of Br- to Cl-
    !=======================================================================
    IF ( H%Cl_conc_Cld > 0.0_dp ) THEN
       H%Br_over_Cl_Cld = H%Br_conc_Cld / H%Cl_conc_Cld   ! in gas, in-cloud
    ENDIF

    IF ( H%Cl_conc_SSA > 0.0_dp ) THEN
       H%Br_over_Cl_SSA = H%Br_conc_SSA / H%Cl_conc_SSA  ! in fine sea salt
    ENDIF

    IF ( H%Cl_conc_SSC > 0.0_dp ) THEN
       H%Br_over_Cl_SSC = H%Br_conc_SSC / H%Cl_conc_SSC  ! in coarse sea salt
    ENDIF

    !=======================================================================
    ! Fraction of SALACL in total fine sea salt
    !=======================================================================
    H%frac_SALACL = C(ind_SALACL) / ( C(ind_SALACL) + C(ind_NIT) + C(ind_SO4) )

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
       Br_conc = 0.0_dp
       Cl_conc = 0.0_dp
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
    V_tot = ( surf_area * r_w / 3.0_dp ) * 1e-3_dp ! L(liq)/cm3(air)

    ! Skip if we are not in cloud
    IF ( V_tot <= 1.0e-20_dp ) THEN
       conc_x = 0.0_dp
       RETURN
    ENDIF

    ! This calculation can be used for both SSA X- concentration and for
    ! those out of cloud only. For X- out of cloud only, V_tot =
    ! V_tot*(1-CF), n_x = n_x*(1-CF), so (1-CF) is canceled.
    ! xnw, 02/05/18
    conc_x = ( n_x / AVO ) / V_tot    ! mol/L
    conc_x = MAX( conc_x, 0.0_dp )

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
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kiir1ltd
!
! !DESCRIPTION: Determine removal rates for both species in an uptake reaction.
!\\
!\\
! !INTERFACE:
!
    FUNCTION kIIR1Ltd( concGas, concEduct, kISource, minLife ) RESULT( kII )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(dp), INTENT(IN)           :: concGas, concEduct
      REAL(dp), INTENT(IN)           :: kISource
      REAL(dp), INTENT(IN), OPTIONAL :: minLife
!
! !RETURN VALUE:
!
      REAL(dp)                       :: kII
!
! !REMARKS:
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(dp) :: kIGas, kIEduct
      REAL(dp) :: lifeA, lifeB, kIMult

      ! Copy kI as calculated assuming no limitation
      kIGas = kISource
      kIEduct = 0.0e+0_dp
      kII = 0.0e+0_dp

      IF (concEduct.lt.100.0e+0_dp) THEN
         kIGas = 0.0e+0_dp
         kIEduct = 0.0e+0_dp
         kII = 0.0e+0_dp
      ELSE
         ! Safe division here is probably overkill - may remove this
         IF (Is_Safe_Div(concGas*kIGas,concEduct)) THEN
            kIEduct = kIGas*concGas/concEduct
            kII = kIGas/concEduct
         ELSE
            kIGas = 0.0e+0_dp
            kIEduct = 0.0e+0_dp
            kII = 0.0e+0_dp
         ENDIF
      ENDIF

      ! Enforce a minimum lifetime?
      IF (PRESENT(minLife)) THEN
         IF ((kIGas.gt.0.0e+0_dp).and.(minLife.gt.0.0e+0_dp)) THEN
            ! Calculate lifetime of each reactant against removal
            lifeA = Safe_Div(1.0e+0_dp,kIGas,0.0e+0_dp)
            lifeB = Safe_Div(1.0e+0_dp,kIEduct,0.0e+0_dp)
            ! Check if either lifetime is "too short"
            IF ((lifeA.lt.lifeB).and.(lifeA.lt.minLife)) THEN
               IF (Is_Safe_Div(concGas*kIGas,concEduct)) THEN
                  kIGas = 1.0e+0_dp/minLife
                  kII = kIGas/concEduct
               ELSE
                  kIGas = 0.0e+0_dp
                  kII = 0.0e+0_dp
               ENDIF
            ELSEIF (lifeB.lt.minLife) THEN
               IF (Is_Safe_Div(concEduct*kIEduct,concGas)) THEN
                  kIEduct = 1.0e+0_dp/minLife
                  kII = kIEduct/concGas
               ELSE
                  kIEduct = 0.0e+0_dp
                  kII = 0.0e+0_dp
               ENDIF
            ENDIF
         ENDIF
      ENDIF

    END FUNCTION kIIR1Ltd
!EOC
END MODULE GcKpp_HetRates
