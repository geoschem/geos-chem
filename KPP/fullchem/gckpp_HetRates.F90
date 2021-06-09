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

MODULE GCKPP_HETRATES
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
  USE State_Chm_Mod,    ONLY : HetState
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
  PUBLIC  :: SET_HET
!
! !PRIVATE MEMBER FUNCTIONS:
!
  ! These functions are used for all mechanisms
  PRIVATE :: HetNO3
  PRIVATE :: HetNO2
  PRIVATE :: HetHO2
  PRIVATE :: HetGLYX
  PRIVATE :: HetMGLY
  PRIVATE :: HetIEPOX
  PRIVATE :: HetIMAE
  PRIVATE :: HetLVOC
  PRIVATE :: HetISOPND
  PRIVATE :: HetISOPNB
  PRIVATE :: HetMACRN
  PRIVATE :: HetMVKN
  PRIVATE :: HetR4N2
  PRIVATE :: HetISN1OG
  PRIVATE :: HetDHDN
  PRIVATE :: HetMONITS
  PRIVATE :: HetMONITU
  PRIVATE :: HetHONIT
  PRIVATE :: HetIONITA
  PRIVATE :: HetMONITA
  PRIVATE :: HetHBr
  PRIVATE :: HetN2O5
  PRIVATE :: N2O5
  PRIVATE :: HO2
  PRIVATE :: CloudHet

  ! New iodine heterogeneous chemistry
  PRIVATE :: HETIUptake
  PRIVATE :: HETIXCycleSSA      !XW modified

  ! Halogen fuctions, XW
  PRIVATE :: HETBrNO3
  PRIVATE :: HETClNO3
  PRIVATE :: HETHOBr_HBr
  PRIVATE :: HETHOBr_HCl
  PRIVATE :: HETClNO3_HBr
  PRIVATE :: HETO3_TCld
  PRIVATE :: HETHOBr_SS
  PRIVATE :: HETClNO3_SS
  PRIVATE :: HETO3_SS
  PRIVATE :: HETHXUptake
  PRIVATE :: HETN2O5_SS
  PRIVATE :: HETHOCl_TCld
  PRIVATE :: HETHOCl_SS
  PRIVATE :: HETHOBr_TCld
  PRIVATE :: HETClNO3_TCld
  PRIVATE :: HETClNO2_TCld
  PRIVATE :: HETOH
  PRIVATE :: HETClNO2
  PRIVATE :: HetNO3_Cl

  ! Halogen gamma calculation and other subrontines, XW
  PRIVATE :: Gamma_O3_Br
  PRIVATE :: Gamma_HX_Uptake
  PRIVATE :: Coth
  PRIVATE :: ReactoDiff_Corr
  PRIVATE :: Gamma_HOBr_CLD
  PRIVATE :: Gamma_HOBr_AER
  PRIVATE :: Gamma_ClNO3_AER
  PRIVATE :: Gamma_ClNO2
  PRIVATE :: Gamma_HOCl_AER
  PRIVATE :: Gamma_HOCl_CLD
  PRIVATE :: Gamma_NO3

  ! These are strat-only reactions
  PRIVATE :: HETClNO3_HCl
  PRIVATE :: HETHOCl_HBr
  PRIVATE :: HETHOCl_HCl
  PRIVATE :: HETBrNO3_HCl
  PRIVATE :: HETN2O5_HCl

  ! These are subfunctions to calculate rates on/in clouds and SSA
  PRIVATE :: CLD_PARAMS
  PRIVATE :: GET_HALIDE_CLDConc
  Private :: Get_Halide_SSAConc
  PRIVATE :: COMPUTE_L2G_LOCAL
  PRIVATE :: CLD1K_XNO3
  PRIVATE :: EPOXUPTK
  PRIVATE :: FCRO2HO2
  PRIVATE :: FYHORO
  PRIVATE :: FYRNO3
  PRIVATE :: ARSL1K
  PRIVATE :: kIIR1Ltd
  PRIVATE :: kIIR1R2L

  ! These are subfunctions related to ice uptake of bromine/chlorine
  PRIVATE :: GET_THETA_ICE
  PRIVATE :: GAMMA_ClNO3_ICE
  PRIVATE :: GAMMA_HOBr_ICE

!
! !PRIVATE DATA MEMBERS:
!
  ! Scalars
  REAL(fp) :: HSO3conc_Cld, SO3conc_Cld, fupdateHOBr, fupdateHOCl
  REAL(fp) :: nitConc_SALA, nitConc_SALC

!$OMP THREADPRIVATE( HSO3conc_Cld, SO3conc_Cld, fupdateHOBr       )
!$OMP THREADPRIVATE( fupdateHOCl, nitConc_SALA, nitConc_SALC      )

! !REMARKS:
!  Need
!  - TOTAREA (previously used for archiving N2O5 hydrolysis in the planeflight
!             diagnostic only)
!  - Air NUM. DENSITY
!  - TEMPERATURE
!  - Aerosol Surface Area
!  - Aerosol Type
!  - Gamma (XSTKCF; sticking factor)
!  - ARR
!  - Species num density (mcl cm-3)
!  - Continental PBL or no?
!  - In stratosphere or no?
!  - Reaction index (e.g. NK1HBr, NK2HBr)
!
! !REFERENCES:
!  Eastham et al., Development and evaluation of the unified tropospheric-
!    stratospheric chemistry extension (UCX) for the global chemistry-transport
!    model GEOS-Chem, Atmos. Env., doi:10.1016/j.atmosenv.2014.02.001, 2014.
!  Fisher et al, Organic nitrate chemistry and its implications for nitrogen
!    budgets in an isoprene- and monoterpene-rich atmosphere: constraints from
!    aircraft (SEAC4RS) and ground-based (SOAS) observations in the Southeast
!    US. Atmos. Chem. Phys., 16, 2961-2990, 2016.
!  Holmes, C.D., Bertram, T. H., Confer, K. L., Ronan, A. C., Wirks, C. K., Graham, K. A.,
!    Shah, V. (2019) The role of clouds in the tropospheric NOx cycle: a new modeling
!    approach for cloud chemistry and its global implications, Geophys. Res. Lett. 46,
!    4980-4990, https://doi.org/10.1029/2019GL081990
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
!
! !REMARKS:
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
      LOGICAL                :: SAFEDIV
      INTEGER                :: IND
      REAL(fp)               :: conc1, conc2, denom

      ! SAVEd scalars
      INTEGER,          SAVE :: id_SALA = -1
      LOGICAL,          SAVE :: FIRST   = .TRUE.

      ! New treatment for educt removal
      REAL(fp)               :: kITemp, kIITemp
      REAL(fp)               :: HetTemp(3)  !new temp array for gamma results

      ! Cloud parameters
      REAL(fp)               :: rLiq, ALiq, VLiq, CLDFr
      REAL(fp)               :: rIce, AIce, VIce

      ! New bromine/chlorine chemistry
      REAL(fp)               :: hConc_Sul
      REAL(fp)               :: hConc_LCl
      REAL(fp)               :: hConc_ICl
      REAL(fp)               :: hConc_SSA
      REAL(fp)               :: hConc_SSC
      REAL(fp)               :: brConc_SALA, clConc_SALA
      REAL(fp)               :: brConc_SALC, clConc_SALC
      REAL(fp)               :: brConc_CldA, clConc_CldA
      REAL(fp)               :: brConc_CldC, clConc_CldC
      REAL(fp)               :: brConc_Cldg, clConc_Cldg
      REAL(fp)               :: brConc_Cld, clConc_Cld
      REAL(fp)               :: pHCloud, pHSSA(2)
      REAL(fp)               :: SSAlk(2)
      REAL(fp)               :: hno3_th, hcl_th, hbr_th

      ! Cloud heterogeneous chemistry
      REAL(fp)               :: gammaLiq, gammaIce

      ! Pointers and objects
      TYPE(HetState), POINTER :: H

      !====================================================================
      ! SET_HET begins here!
      !====================================================================

      ! Zero scalars and arrays
      ! NOTE: we now get several quantities from gckpp_Global
      HetTemp       = 0.0_fp

      ! Point to the HetInfo objec
      H             => State_Chm%HetInfo

      !--------------------------------------------------------------------
      !  Calculate parameters for cloud halogen chemistry
      !  under the new scheme (SDE 2016-12-21)
      !--------------------------------------------------------------------

      ! Get cloud physical parameters
      CALL Cld_Params( I,     J,     L,    State_Met, rLiq, ALiq,          &
                       VLiq,  rIce,  AIce, VIce,      CLDFr          )

      ! Retrieve cloud pH and alkalinity
      pHCloud  = State_Chm%pHCloud(I,J,L)
      pHSSA(:) = State_Chm%IsorropAeropH(I,J,L,:)

      ! Since we have Alk tracers, State_Chm%SSAlk may be removed in future
      ! Here still use State_Chm%SSAlk instead of Alk tracers,
      ! to skip unit conversion. xnw 1/24/18
      SSAlk(1:2) = State_Chm%SSAlk(I,J,L,1:2)

      ! Estimate liquid phase pH (H+ concentration)
      hConc_Sul = 10.0**(-1.0_fp*pHSSA(1))
      hConc_LCl = 10.0**(-1.0_fp*pHCloud)
      hConc_ICl = 10.0**(-4.5_fp)
      hConc_SSA = hConc_Sul
      hConc_SSC = 10.0**(-5.0_fp)

      ! Get halide concentrations, in cloud
      conc1 = C(ind_HBr) + ( C(ind_BrSALA) * 0.7_fp ) + C(ind_BrSALC)
      conc2 = C(ind_HCl) + ( C(ind_SALACL) * 0.7_fp ) + C(ind_SALCCL)

      CALL Get_Halide_CldConc( conc1,       conc2,      VLiq,  VIce,         &
                               VAir,        Temp,       CLDFr, pHCloud,      &
                               brConc_Cldg, clConc_Cldg                     )

      brConc_Cld  = brConc_Cldg

      ! Avoid div-by-zero (all three expressions use the same denominator)
      denom = C(ind_HBr) + ( C(ind_BrSALA) * 0.7_fp ) + C(ind_BrSALC)
      IF ( denom > 0.0_fp ) THEN
         brConc_Cldg = ( brConc_Cld * C(ind_HBr   )          ) / denom
         brConc_CldA = ( brConc_Cld * C(ind_BrSALA) * 0.7_fp ) / denom
         brConc_CldC = ( brConc_Cld * C(ind_BrSALC)          ) / denom
      ELSE
         brConc_Cldg = 0.0_fp
         brConc_CldA = 0.0_fp
         brConc_CldC = 0.0_fp
      ENDIF

      ! Enforce minimum values
      brConc_Cldg = MAX( brConc_Cldg, 1.0e-20_fp )
      brConc_CldA = MAX( brConc_CldA, 1.0e-20_fp )
      brConc_CldC = MAX( brConc_CldC, 1.0e-20_fp )

      clConc_Cld  = clConc_Cldg

      ! Avoid div-by-zero (all three expressions use the same denominator)
      denom = C(ind_HCl) + ( C(ind_SALACL) * 0.7_fp ) + C(ind_SALCCL)
      IF ( denom > 0.0_fp ) THEN
         clConc_Cldg = ( clConc_Cld * C(ind_HCl   )          ) / denom
         clConc_CldA = ( clConc_Cld * C(ind_SALACL) * 0.7_fp ) / denom
         clConc_CldC = ( clConc_Cld * C(ind_SALCCL)          ) / denom
      ELSE
         clConc_Cldg = 0.0_fp
         clConc_CldA = 0.0_fp
         clConc_CldC = 0.0_fp
      ENDIF

      ! Enforce minimum values
      clConc_Cldg = MAX( clConc_Cldg, 1.0e-20_fp )
      clConc_CldA = MAX( clConc_CldA, 1.0e-20_fp )
      clConc_CldC = MAX( clConc_CldC, 1.0e-20_fp )

      ! Get halide concentrations, in aerosol
      ! Note that here Br/ClSALA = Br/Cl- in (SALA + sulfate)
      ! Get the concentration of Br-
      CALL Get_Halide_SSAConc(C(ind_BrSALA),AClAREA,  AClRADI,  brConc_SALA)
      CALL Get_Halide_SSAConc(C(ind_BrSALC),xArea(12),xRadi(12),brConc_SALC)

      ! Get the concentration of Cl-
      CALL Get_Halide_SSAConc(C(ind_SALACL),AClAREA,  AClRADI,  clConc_SALA)
      CALL Get_Halide_SSAConc(C(ind_SALCCL),xArea(12),xRadi(12),clConc_SALC)

      ! Get the concentration of NO3-
      CALL Get_Halide_SSAConc(C(ind_NIT),   AClAREA,  AClRADI,  nitConc_SALA)
      CALL Get_Halide_SSAConc(C(ind_NITs),  xArea(12),xRadi(12),nitConc_SALC)

      ! Get theta for ice cloud uptake
      CALL Get_Theta_Ice( C(ind_HNO3), C(ind_HCl), C(ind_HBr), Temp,         &
                          hno3_th,     hcl_th,     hbr_th                   )

      !--------------------------------------------------------------------
      !  Get parameters for HOBr + S(IV)
      !--------------------------------------------------------------------

      ! Cloud bisulfite (HSO3-) concentration [mol/l] from sulfate_mod.F
      HSO3conc_Cld = State_Chm%HSO3_AQ(I,J,L)

      ! Cloud sulfite (SO3--) concentration [mol/l] from sulfate_mod.F
      SO3conc_Cld  = State_Chm%SO3_AQ(I,J,L)

      ! Avoid div-by-zero issues in GAMMA_HOBr_X
      !IF ( HSO3conc_Cld <= 0.0_fp) HSO3conc_Cld = 1e-20_fp
      !IF (  SO3conc_Cld  <= 0.0_fp)  SO3conc_Cld = 1e-20_fp

      ! Correction factor for HOBr removal by SO2 [unitless]
      fupdateHOBr  = State_Chm%fupdateHOBr(I,J,L)
      fupdateHOCl  = State_Chm%fupdateHOCl(I,J,L)

      !--------------------------------------------------------------------
      ! Calculate and pass het rates to the KPP rate array
      !--------------------------------------------------------------------

      ! Zero the HET array
      HET = 0.0_dp

      ! Calculate genuine first-order uptake reactions first
      HET(ind_HO2,    1) = HetHO2(     H%HO2%MW_g,    2E-1_fp      )
      HET(ind_NO2,    1) = HetNO2(     H%NO2%MW_g,    1E-4_fp      )
      HET(ind_NO3,    1) = HetNO3(     H%NO3%MW_g,    1E-1_fp      )
      HET(ind_GLYX,   1) = HetGLYX(    H%GLYX%MW_g,   1E-1_fp      )
      HET(ind_MGLY,   1) = HetMGLY(    H%MGLY%MW_g,   1E-1_fp      )
      HET(ind_IEPOXA, 1) = HetIEPOX(   H%IEPOXA%MW_g, 1E-1_fp      )
      HET(ind_IEPOXB, 1) = HetIEPOX(   H%IEPOXB%MW_g, 1E-1_fp      )
      HET(ind_IEPOXD, 1) = HetIEPOX(   H%IEPOXD%MW_g, 1E-1_fp      )
      HET(ind_HMML,   1) = HetIMAE(    H%HMML%MW_g,   1E-1_fp      )
      HET(ind_PYAC,   1) = HetMGLY(    H%PYAC%MW_g,   1E-1_fp      )
      HET(ind_ICHE,   1) = HetIEPOX(   H%ICHE%MW_g,   1E-1_fp      )
      HET(ind_LVOC,   1) = HetLVOC(    H%LVOC%MW_g,   1E+0_fp      )
      HET(ind_IHN1,   1) = HetISOPND(  H%IHN1%MW_g,   5E-3_fp      )
      HET(ind_IHN2,   1) = HetISOPNB(  H%IHN2%MW_g,   5E-2_fp      )
      HET(ind_IHN3,   1) = HetISOPNB(  H%IHN3%MW_g,   5E-3_fp      )
      HET(ind_IHN4,   1) = HetISOPND(  H%IHN4%MW_g,   5E-3_fp      )
      HET(ind_INPB,   1) = HetISOPNB(  H%INPB%MW_g,   5E-3_fp      )
      HET(ind_INPD,   1) = HetISOPND(  H%INPD%MW_g,   5E-3_fp      )
      HET(ind_MCRHN,  1) = HetMACRN(   H%MCRHN%MW_g,  5E-3_fp      )
      HET(ind_MCRHNB, 1) = HetMACRN(   H%MCRHNB%MW_g, 5E-3_fp      )
      HET(ind_MVKN,   1) = HetMVKN(    H%MVKN%MW_g,   5E-3_fp      )
      HET(ind_R4N2,   1) = HetR4N2(    H%R4N2%MW_g,   5E-3_fp      )
      HET(ind_IDN,    1) = HetDHDN(    H%IDN%MW_g,    5E-3_fp      )
      HET(ind_ITHN,   1) = HetDHDN(    H%ITHN%MW_g,   5E-3_fp      )
      HET(ind_ITCN,   1) = HetDHDN(    H%ITCN%MW_g,   5E-3_fp      )
      HET(ind_MONITS, 1) = HetMONITS(  H%MONITS%MW_g, 1E-2_fp      )
      HET(ind_MONITU, 1) = HetMONITU(  H%MONITU%MW_g, 1E-2_fp      )
      HET(ind_HONIT,  1) = HetHONIT(   H%HONIT%MW_g,  1E-2_fp      )
      HET(ind_IONITA, 1) = HetIONITA(  H%IONITA%MW_g, 1E-1_fp      )
      HET(ind_MONITA, 1) = HetMONITA(  H%MONITA%MW_g, 1E-1_fp      )

      !======================================================================
      ! NOy uptake in clouds
      !======================================================================
      HET(ind_NO2,    1) = HET(ind_NO2, 1)                                   &
                         + CloudHet( CldFr, Aliq, Aice,                      &
                                     rLiq,  rIce, Temp,   XDenA,            &
                                     H%NO2%MW_g,  1e-8_fp, 0.0_fp           )
      HET(ind_NO3,    1) = HET(ind_NO3, 1)                                   &
                         + CloudHet( CldFr, Aliq, Aice,                      &
                                     rLiq,  rIce, Temp,    XDenA,           &
                                     H%NO3%MW_g,  0.002_fp, 0.001_fp        )

      ! Reactive uptake coefficient for N2O5 on liquid water cloud
      ! Value is 0.03 at 298 K (JPL, Burkholder et al., 2015)
      ! For temperature dependence, JPL recommends the same as
      ! sulfuric acid aerosol at zero percent H2SO4, which is 0.019 at 298 K.
      ! Then apply constant scale factor (0.03/0.019)
      gammaLiq = ( 0.03_fp / 0.019_fp ) * &
         exp( -25.5265_fp + 9283.76_fp / Temp - 851801.0_fp / Temp**2 )

      HET(ind_N2O5,   3) = CloudHet( CldFr, Aliq, Aice,                      &
                                     rLiq,  rIce, Temp,    XDenA,           &
                                     H%N2O5%MW_g, gammaLiq, 0.02_fp         )

      ! Now calculate reaction rates where the educt can be consumed.
      ! kIIR1Ltd: Assume that the first reactant is limiting. Assume that the
      ! second reactant is "abundant" and calculate the overall rate based on
      ! the uptake rate of the first reactant only.

      HetTemp(1:2)                 = HETN2O5( 1.08E2_fp, 1E-1_fp, CldFr )
      kITemp                       = HetTemp(1)
      HET(ind_N2O5,  1)            = kIIR1Ltd(C(ind_N2O5), C(ind_H2O), kITemp)
      State_Chm%GammaN2O5(I,J,L,1) = HetTemp(2)

      !======================================================================
      ! Br/Cl heterogeneous chemistry
      !======================================================================

      !----------------------------------------------------------------------
      ! Cl- enhanced NO3 hydrolysis (XW 2018-03-16)
      !----------------------------------------------------------------------
      kITemp          = HetNO3_Cl( XDenA, AClRADI,     (1-CldFr)*AClAREA,    &
                                   Temp, clConc_SALA, 1,                    &
                                   H                                        )

      HET(ind_NO3, 2) = kITemp

      kITemp          = HetNO3_Cl( XDenA, XRADI(12),  (1-CldFr)*XAREA(12),   &
                                   Temp, clConc_SALC, 2,                    &
                                   H                                        )

      HET(ind_NO3, 3) = kITemp

      !----------------------------------------------------------------------
      ! ClNO3 and BrNO3 hydrolysis (update: XW 2019-06-08)
      !----------------------------------------------------------------------
      gammaLiq = MAX( 0.0021_fp*Temp - 0.561_fp, 1e-30_fp )
      gammaIce = 5.3e-4_fp * exp(1100.0_fp / Temp)
      kITemp            = HETBrNO3( XDenA, Temp, CldFr, H )
      kITemp = kITemp   + CloudHet( CldFr, Aliq,  Aice,                      &
                                    rLiq,  rIce,  Temp, XDenA,               &
                                    H%BrNO3%MW_g, gammaLiq, gammaIce        )

      HET(ind_BrNO3, 1) = kIIR1Ltd( C(ind_BrNO3), C(ind_H2O), kITemp        )
      kITemp            = HETClNO3( XDenA,       Temp, clConc_SALA,         &
                                    brConc_SALA, CldFr, H                   )

      kITemp            = kITemp                                             &
                        + HETClNO3_TCld(                                     &
                                    XDenA,       rLiq,  rIce,                &
                                    ALiq,        AIce,  VAir,                &
                                    Temp,        CldFr, clConc_CldA,         &
                                    clConc_CldC, clConc_Cldg, brConc_CldA,   &
                                    brConc_CldC, brConc_Cldg, hno3_th,       &
                                    hcl_th,      hbr_th,      7,             &
                                    H                                       )

      HET(ind_ClNO3, 1) = kIIR1Ltd( C(ind_ClNO3), C(ind_H2O), kITemp        )

      !----------------------------------------------------------------------
      ! HOBr + HBr (update: XW 2019-06-08)
      !----------------------------------------------------------------------
      kITemp = HETHOBr_HBr( H%HOBr%MW_g, 0.0_fp, Input_Opt                  )

      kITemp = kITemp                                                        &
             + HETHOBr_TCld( XDenA,       rLiq,        rIce,                 &
                             ALiq,        AIce,        VAir,                 &
                             Temp,        CldFr,       hConc_Sul,            &
                             hConc_LCl,   hConc_ICl,   clConc_CldA,          &
                             clConc_CldC, clConc_Cldg, brConc_CldA,          &
                             brConc_CldC, brConc_Cldg, HSO3conc_Cld,         &
                             SO3conc_Cld, hno3_th,     hcl_th,               &
                             hbr_th,      6,           H                    )

      HET(ind_HOBr,  1) = kIIR1Ltd( C(ind_HOBr), C(ind_HBr), kITemp         )

      !----------------------------------------------------------------------
      ! HOBr + HCl (update: XW 2019-06-08)
      !----------------------------------------------------------------------
      kITemp = HETHOBr_HCl( H%HOBr%MW_g, 0.0_fp, Input_Opt                  )

      kITemp = kITemp                                                        &
             + HETHOBr_TCld( XDenA,       rLiq,        rIce,                 &
                             ALiq,        AIce,        VAir,                 &
                             Temp,        CldFr,       hConc_Sul,            &
                             hConc_LCl,   hConc_ICl,   clConc_CldA,          &
                             clConc_CldC, clConc_Cldg, brConc_CldA,          &
                             brConc_CldC, brConc_Cldg, HSO3conc_Cld,         &
                             SO3conc_Cld, hno3_th,     hcl_th,               &
                             hbr_th,      5,           H                    )

      HET(ind_HOBr,  2) = kIIR1Ltd( C(ind_HOBr), C(ind_HCl), kITemp         )

      !----------------------------------------------------------------------
      ! HOBr + BrSalA/C (update: XW 2019-06-08)
      !----------------------------------------------------------------------

      ! First consider reaction in troposphere cloud
      kITemp = HETHOBr_TCld( XDenA,       rLiq,        rIce,                 &
                             ALiq,        AIce,        VAir,                 &
                             Temp,       CldFr,        hConc_Sul,            &
                             hConc_LCl,   hConc_ICl,   clConc_CldA,          &
                             clConc_CldC, clConc_Cldg, brConc_CldA,          &
                             brConc_CldC, brConc_Cldg, HSO3conc_Cld,         &
                             SO3conc_Cld, hno3_th,     hcl_th,               &
                             hbr_th,      3,           H                    )

      ! Then calculate reaction on aerosols
      kITemp = kITemp +                                                      &
               HETHOBr_SS( XDenA,       AClRADI,     (1-CldFr)*AClAREA,      &
                           SSAlk(1),    Temp,        hConc_SSA,              &
                           clConc_SALA, brConc_SALA, 2,                      &
                           H                                                )

      HET(ind_HOBr,  5) = kIIR1Ltd( C(ind_HOBr), C(ind_BrSALA), kITemp      )

      ! First consider reaction in troposphere cloud
      kITemp = HETHOBr_TCld( XDenA,       rLiq,        rIce,                 &
                             ALiq,        AIce,        VAir,                 &
                             Temp,        CldFr,       hConc_Sul,            &
                             hConc_LCl,   hConc_ICl,   clConc_CldA,          &
                             clConc_CldC, clConc_Cldg, brConc_CldA,          &
                             brConc_CldC, brConc_Cldg, HSO3conc_Cld,         &
                             SO3conc_Cld, hno3_th,     hcl_th,               &
                             hbr_th,      4,           H                    )

      ! Then calculate reaction on aerosols
      kITemp = kITemp                                                        &
             + HETHOBr_SS( XDenA,       xRadi(12),   (1-CldFr)*xArea(12),    &
                           SSAlk(2),    Temp,        hConc_SSC,              &
                           clConc_SALC, brConc_SALC, 2,                      &
                           H                                                )

      HET(ind_HOBr,  6) = kIIR1Ltd( C(ind_HOBr),  C(ind_BrSALC), kITemp     )

      !----------------------------------------------------------------------
      ! HOBr + Cl-(p) (update: XW 2019-06-08)
      !----------------------------------------------------------------------
      ! First consider reaction in troposphere cloud
      kITemp = HETHOBr_TCld( XDenA,       rLiq,        rIce,                 &
                             ALiq,        AIce,        VAir,                 &
                             Temp,        CldFr,       hConc_Sul,            &
                             hConc_LCl,   hConc_ICl,   clConc_CldA,          &
                             clConc_CldC, clConc_Cldg, brConc_CldA,          &
                             brConc_CldC, brConc_Cldg, HSO3conc_Cld,         &
                             SO3conc_Cld, hno3_th,     hcl_th,               &
                             hbr_th,      1,           H                    )

      ! Then calculate reaction on aerosols
      kITemp = kITemp                                                        &
             + HETHOBr_SS( XDenA,       AClRADI,     (1-CldFr)*AClAREA,      &
                           SSAlk(1),    Temp,        hConc_SSA,              &
                           clConc_SALA, brConc_SALA, 1,                      &
                           H                                                )

      HET(ind_HOBr,  3) = kIIR1Ltd( C(ind_HOBr), C(ind_SALACL), kITemp      )

      ! First consider reaction in troposphere cloud
      kITemp = HETHOBr_TCld( XDenA,       rLiq,        rIce,                 &
                             ALiq,        AIce,        VAir,                 &
                             Temp,        CldFr,       hConc_Sul,            &
                             hConc_LCl,   hConc_ICl,   clConc_CldA,          &
                             clConc_CldC, clConc_Cldg, brConc_CldA,          &
                             brConc_CldC, brConc_Cldg, HSO3conc_Cld,         &
                             SO3conc_Cld, hno3_th,     hcl_th,               &
                             hbr_th,      2,           H                    )

      ! Then calculate reaction on aerosols
      kITemp = kITemp                                                        &
             + HETHOBr_SS( XDenA,       xRadi(12),   (1-CldFr)*xArea(12),    &
                           SSAlk(2),    Temp,       hConc_SSC,               &
                           clConc_SALC, brConc_SALC, 1,                      &
                           H                                                )

      HET(ind_HOBr,  4) = kIIR1Ltd( C(ind_HOBr), C(ind_SALCCL), kITemp      )

      !----------------------------------------------------------------------
      ! HOBr + HSO3-(aq) (update: XW 2019-06-08)
      !----------------------------------------------------------------------

      ! This reaction is first order, so no kII calculation is required
      kITemp = HETHOBr_TCld( XDenA,       rLiq,        rIce,                 &
                             ALiq,        AIce,        VAir,                 &
                             Temp,        CldFr,       hConc_Sul,            &
                             hConc_LCl,   hConc_ICl,   clConc_CldA,          &
                             clConc_CldC, clConc_Cldg, brConc_CldA,          &
                             brConc_CldC, brConc_Cldg, HSO3conc_Cld,         &
                             SO3conc_Cld, hno3_th,     hcl_th,               &
                             hbr_th,      7,           H                    )

      ! Make sure sulfate produced is less than SO2 available (qjc, 06/20/16)
      HET(ind_HOBr,  7) = kITemp * fupdateHOBr

      !----------------------------------------------------------------------
      ! HOBr + SO3--(aq) (update: XW 2019-06-08)
      !----------------------------------------------------------------------
      ! This reaction is first order, so no kII calculation is required
      kITemp = HETHOBr_TCld( XDenA,       rLiq,        rIce,                 &
                             ALiq,        AIce,        VAir,                 &
                             Temp,        CldFr,       hConc_Sul,            &
                             hConc_LCl,   hConc_ICl,   clConc_CldA,          &
                             clConc_CldC, clConc_Cldg, brConc_CldA,          &
                             brConc_CldC, brConc_Cldg, HSO3conc_Cld,         &
                             SO3conc_Cld, hno3_th,     hcl_th,               &
                             hbr_th,      8,           H                    )

      ! Make sure sulfate produced is less than SO2 available (qjc, 06/20/16)
      HET(ind_HOBr,  8) = kITemp * fupdateHOBr

      !---------------------------------------------------------------------
      ! ClNO3 + BrSALA/C (update: XW 2019-06-08)
      !---------------------------------------------------------------------
      ! First consider reaction in troposphere cloud
      kITemp = HETClNO3_TCld( XDenA,       rLiq,        rIce,                &
                              ALiq,        AIce,        VAir,                &
                              Temp,        CldFr,       clConc_CldA,         &
                              clConc_CldC, clConc_Cldg, brConc_CldA,         &
                              brConc_CldC, brConc_Cldg, hno3_th,             &
                              hcl_th,      hbr_th,      3,                   &
                              H                                             )

      ! Then calculate reaction on aerosols
      kITemp = kITemp                                                        &
             + HETClNO3_SS( XDenA,       AClRADI, (1-CldFr)*AClAREA,         &
                            SSAlk(1),    Temp,    clConc_SALA,               &
                            brConc_SALA, 2,       1,                         &
                            H                                               )

      HET(ind_ClNO3, 4) = kIIR1Ltd( C(ind_ClNO3), C(ind_BrSALA), kITemp     )

      ! First consider reaction in troposphere cloud
      kITemp = HETClNO3_TCld( XDenA,       rLiq,        rIce,                &
                              ALiq,        AIce,        VAir,                &
                              Temp,        CldFr,       clConc_CldA,         &
                              clConc_CldC, clConc_Cldg, brConc_CldA,         &
                              brConc_CldC, brConc_Cldg, hno3_th,             &
                              hcl_th,      hbr_th,      4,                   &
                              H                                             )

      ! Then calculate reaction on aerosols
      kITemp = kITemp                                                        &
             + HETClNO3_SS( XDenA,       xRadi(12), (1-CldFr)*xArea(12),     &
                            SSAlk(2),    Temp,      clConc_SALC,             &
                            brConc_SALC, 2,         2,                       &
                            H                                               )

      HET(ind_ClNO3, 5) = kIIR1Ltd( C(ind_ClNO3), C(ind_BrSALC), kITemp     )

      !---------------------------------------------------------------------
      ! ClNO3 + Cl-(p) (update: XW 2019-06-08)
      !---------------------------------------------------------------------
      ! First consider reaction in troposphere cloud
      kITemp = HETClNO3_TCld( XDenA,       rLiq,        rIce,                &
                              ALiq,        AIce,        VAir,                &
                              Temp,        CldFr,       clConc_CldA,         &
                              clConc_CldC, clConc_Cldg, brConc_CldA,         &
                              brConc_CldC, brConc_Cldg, hno3_th,             &
                              hcl_th,      hbr_th,      1,                   &
                              H                                             )

      ! Then calculate reaction on aerosols
      kITemp = kITemp                                                        &
             + HETClNO3_SS( XDenA,       AClRADI, (1-CldFr)*AClAREA,         &
                            SSAlk(1),    Temp,    clConc_SALA,               &
                            brConc_SALA, 1,       1,                         &
                            H                                               )

      HET(ind_ClNO3, 6) = kIIR1Ltd( C(ind_ClNO3), C(ind_SALACL), kITemp     )

      ! First consider reaction in troposphere cloud
      kITemp = HETClNO3_TCld( XDenA,       rLiq,        rIce,                &
                              ALiq,        AIce,        VAir,                &
                              Temp,        CldFr,       clConc_CldA,         &
                              clConc_CldC, clConc_Cldg, brConc_CldA,         &
                              brConc_CldC, brConc_Cldg, hno3_th,             &
                              hcl_th,      hbr_th,      2,                   &
                              H                                             )

      ! Then calculate reaction on aerosols out of cloud
      kITemp = kITemp                                                        &
             + HETClNO3_SS( XDenA,       xRadi(12), (1-CldFr)*xArea(12),     &
                            SSAlk(2),    Temp,      clConc_SALC,             &
                            brConc_SALC, 1,         2,                       &
                            H                                               )

      HET(ind_ClNO3, 7) = kIIR1Ltd( C(ind_ClNO3), C(ind_SALCCL), kITemp     )

      !--------------------------------------------------------------------
      ! ClNO3 + HCl (update: XW 2019-06-08)
      !--------------------------------------------------------------------
      kITemp = HETClNO3_HCl( H%ClNO3%MW_g, 0.0_fp                           )

      kITemp = kITemp                                                        &
             + HETClNO3_TCld( XDenA,       rLiq,        rIce,                &
                              ALiq,        AIce,        VAir,                &
                              Temp,        CldFr,       clConc_CldA,         &
                              clConc_CldC, clConc_Cldg, brConc_CldA,         &
                              brConc_CldC, brConc_Cldg, hno3_th,             &
                              hcl_th,      hbr_th,      5,                   &
                              H                                             )

      HET(ind_ClNO3, 2) = kIIR1Ltd( C(ind_ClNO3), C(ind_HCl), kITemp        )

      !--------------------------------------------------------------------
      ! ClNO3 + HBr (update: XW 2019-06-08)
      !--------------------------------------------------------------------
      kITemp = HETClNO3_HBr( H%ClNO3%MW_g, 0.0_fp, Input_Opt                )

      kITemp = kITemp                                                        &
             + HETClNO3_TCld( XDenA,       rLiq,        rIce,                &
                              ALiq,        AIce,        VAir,                &
                              Temp,        CldFr,       clConc_CldA,         &
                              clConc_CldC, clConc_Cldg, brConc_CldA,         &
                              brConc_CldC, brConc_Cldg, hno3_th,             &
                              hcl_th,      hbr_th,      6,                   &
                              H                                             )

      HET(ind_ClNO3, 3) = kIIR1Ltd( C(ind_ClNO3), C(ind_HBr), kITemp        )

      !----------------------------------------------------------------------
      ! HOCl + HCl and HOCl + HBr (update: XW 2019-06-08)
      !----------------------------------------------------------------------
      kITemp = HETHOCl_HCl( H%HOCl%MW_g, 0.0_fp, Input_opt                  )
      kITemp = kITemp                                                        &
             + HETHOCl_TCld( XDenA,       rLiq,        rIce,                 &
                             ALiq,        AIce,        VAir,                 &
                             Temp,        CldFr,       hConc_Sul,            &
                             hConc_LCl,   hConc_ICl,   clConc_CldA,          &
                             clConc_CldC, clConc_Cldg, HSO3conc_Cld,         &
                             SO3conc_Cld, hcl_th,      3,                    &
                             H                                              )

      HET(ind_HOCl,  1) = kIIR1Ltd( C(ind_HOCl), C(ind_HCl), kITemp         )

      kiTemp            = HETHOCl_HBr( H%HOCl%MW_g, 0E+0_fp, Input_Opt      )

      HET(ind_HOCl,  2) = kIIR1Ltd( C(ind_HOCl), C(ind_HBr), kiTemp         )

      !----------------------------------------------------------------
      ! HOCl + Cl-(p) (update: XW 2019-06-08)
      !----------------------------------------------------------------
      ! First consider reaction in troposphere cloud
      kITemp = HETHOCl_TCld( XDenA,       rLiq,        rIce,                 &
                             ALiq,        AIce,        VAir,                 &
                             Temp,        CldFr,       hConc_Sul,            &
                             hConc_LCl,   hConc_ICl,   clConc_CldA,          &
                             clConc_CldC, clConc_Cldg, HSO3conc_Cld,         &
                             SO3conc_Cld, hcl_th,      1,                    &
                             H                                              )

      ! Then calculate reaction on aerosols
      kITemp = kITemp                                                        &
             + HETHOCl_SS( XDenA,       AClRADI, (1-CldFr)*AClAREA,          &
                           SSAlk(1),    Temp,   hConc_SSA,                   &
                           clConc_SALA, H                                   )

      HET(ind_HOCl,  3) = kIIR1Ltd( C(ind_HOCl), C(ind_SALACL), kITemp      )

      ! First consider reaction in troposphere cloud
      kITemp = HETHOCl_TCld( XDenA,       rLiq,        rIce,                 &
                             ALiq,        AIce,        VAir,                 &
                             Temp,        CldFr,       hConc_Sul,            &
                             hConc_LCl,   hConc_ICl,   clConc_CldA,          &
                             clConc_CldC, clConc_Cldg, HSO3conc_Cld,         &
                             SO3conc_Cld, hcl_th,      2,                    &
                             H                                              )

      ! Then calculate reaction on aerosols
      kITemp = kITemp                                                        &
             + HETHOCl_SS( XDenA,       xRadi(12), (1-CldFr)*xArea(12),      &
                           SSAlk(2),    Temp,      hConc_SSC,                &
                           clConc_SALC, H                                   )

      HET(ind_HOCl,  4) = kIIR1Ltd( C(ind_HOCl), C(ind_SALCCL), kITemp      )

      !----------------------------------------------------------------
      ! HOCl + HSO3-/SO3--(aq)
      !----------------------------------------------------------------
      ! This reaction is first order, so no kII calculation is required
      kITemp = HETHOCl_TCld( XDenA,       rLiq,        rIce,                 &
                             ALiq,        AIce,        VAir,                 &
                             Temp,        CldFr,       hConc_Sul,            &
                             hConc_LCl,   hConc_ICl,   clConc_CldA,          &
                             clConc_CldC, clConc_Cldg, HSO3conc_Cld,         &
                             SO3conc_Cld, hcl_th,      4,                    &
                             H                                              )

      ! Make sure sulfate produced is less than SO2 available
      HET(ind_HOCl,  5) = kITemp * fupdateHOCl

      ! This reaction is first order, so no kII calculation is required
      kITemp = HETHOCl_TCld( XDenA,       rLiq,        rIce,                 &
                             ALiq,        AIce,        VAir,                 &
                             Temp,        CldFr,       hConc_Sul,            &
                             hConc_LCl,   hConc_ICl,   clConc_CldA,          &
                             clConc_CldC, clConc_Cldg, HSO3conc_Cld,         &
                             SO3conc_Cld, hcl_th,      5,                    &
                             H                                              )

      ! Make sure sulfate produced is less than SO2 available
      HET(ind_HOCl,  6) = kITemp * fupdateHOCl

      !----------------------------------------------------------------
      ! O3 + Br- calculation (TS index: hhc12)
      !----------------------------------------------------------------
      kITemp = HETO3_TCld( XDenA,       rLiq,        rIce,                   &
                           ALiq,        AIce,        VAir,                   &
                           Temp,        CldFr,       brConc_CldA,            &
                           brConc_CldC, brConc_Cldg, C(ind_O3),              &
                           3,           H                                   )

      HET(ind_O3,    1) = kIIR1Ltd( C(ind_O3), C(ind_HBr), kITemp           )

      !----------------------------------------------------------------
      ! O3 + BrSALA/C calculations (TMS index: hhc13/14)
      !----------------------------------------------------------------
      kITemp = HETO3_TCld( XDenA,       rLiq,        rIce,                   &
                           ALiq,        AIce,        VAir,                   &
                           Temp,        CldFr,       brConc_CldA,            &
                           brConc_CldC, brConc_Cldg, C(ind_O3),              &
                           1,           H                                   )

      kITemp = kITemp                                                        &
             + HETO3_SS( XDenA,         AClRADI, (1-CldFr)*AClAREA,          &
                         SSAlk(1),      Temp,    brConc_SALA,                &
                         C(ind_O3),     H                                   )

      HET(ind_O3,    2) = kIIR1Ltd( C(ind_O3), C(ind_BrSALA), kITemp        )

      kITemp  = HETO3_TCld( XDenA,       rLiq,        rIce,                  &
                            ALiq,        AIce,        VAir,                  &
                            Temp,        CldFr,       brConc_CldA,           &
                            brConc_CldC, brConc_Cldg, C(ind_O3),             &
                            2,           H                                  )

      kITemp = kITemp                                                        &
             + HETO3_SS( XDenA,          xRadi(12), (1-CldFr)*xArea(12),     &
                         SSAlk(2),       Temp,      brConc_SALC,             &
                         C(ind_O3),      H                                  )

      HET(ind_O3,    3) = kIIR1Ltd( C(ind_O3), C(ind_BrSALC), kITemp        )

      !----------------------------------------------------------------
      ! Br uptake calculation - forms BrSALA/C (TS index: hhc17/18)
      !----------------------------------------------------------------
      ! First-order reactions, no calculation of kII required
      kITemp            = HETHXUptake( XDenA, AClRADI, (1-CldFr)*ACLAREA,    &
                                       Temp, 2                              )
      HET(ind_HBr,   1) = kITemp

      kITemp            = HETHXUptake( XDenA, xRadi(12), (1-CldFr)*xArea(12),&
                                       Temp, 2                              )
      HET(ind_HBr,   2) = kITemp

      !----------------------------------------------------------------
      ! BrNO3 + HCl in stratosphere
      !----------------------------------------------------------------
      ! NOTE: the restriction of these reactions to the troposphere has been
      ! restored - TMS (2017/04/06 )
      HET(ind_BrNO3, 2) = kIIR1Ltd( C(ind_BrNO3), C(ind_HCl),                &
                                    HETBrNO3_HCl( 1.42E2_fp, 0.0_fp,         &
                                    Input_Opt )                             )

      !----------------------------------------------------------------
      ! N2O5 + HCl in stratosphere
      !----------------------------------------------------------------
      ! NOTE: this extension of calculation in troposphere has been removed
      !  (TMS 17/04/10)
      kITemp            = HETN2O5_HCl( 1.08E2_fp, 0.0_fp, Input_Opt         )
      HET(ind_N2O5,  2) = kIIR1Ltd( C(ind_N2O5), C(ind_HCl), kITemp         )

      !----------------------------------------------------------------
      ! Reaction of N2O5 with Cl-, XW 1/24/18)
      !----------------------------------------------------------------
      HetTemp(1:3)      = HETN2O5_SS( CldFr, 1 )
      HET(ind_N2O5,  4) = kIIR1Ltd( C(ind_N2O5), C(ind_SALACL), HetTemp(1)  )
      State_Chm%GammaN2O5(I,J,L,2) = HetTemp(1)

      HetTemp(1:3)      = HETN2O5_SS( CldFr, 2)
      HET(ind_N2O5,  5) = kIIR1Ltd( C(ind_N2O5), C(ind_SALCCL), HetTemp(1)  )
      State_Chm%GammaN2O5(I,J,L,3) = HetTemp(1)

      !----------------------------------------------------------------
      ! Reaction of OH with Cl-, XW 3/12/18)
      !----------------------------------------------------------------
      kITemp          = HETOH( XDenA, AClRADI, AClAREA, Temp, clConc_SALA   )
      HET(ind_OH,  1) = kIIR1Ltd( C(ind_OH), C(ind_SALACL), kITemp          )

      kITemp          = HETOH( XDenA, xRadi(12), xArea(12),                  &
                               Temp, clConc_SALC                            )
      HET(ind_OH,  2) = kIIR1Ltd( C(ind_OH), C(ind_SALCCL), kITemp          )

      !----------------------------------------------------------------
      ! Reaction of ClNO2 with Cl-, XW 3/12/18)
      !----------------------------------------------------------------

      ! First consider reaction in troposphere cloud
      kITemp = HETClNO2_TCld( XDenA,       rLiq,        rIce,                &
                              ALiq,        AIce,        VAir,                &
                              Temp,        CldFr,       pHCloud,             &
                              clConc_CldA, clConc_CldC, clConc_Cldg,         &
                              brConc_CldA, brConc_CldC, brConc_Cldg,         &
                              1,           H                                )

      ! Then calculate reaction on aerosols out of cloud
      kITemp = kITemp                                                        &
             + HETClNO2( XDenA,       AClRADI,     (1-CldFr)*AClAREA,        &
                         SSAlk(1),    Temp,       pHSSA(1),                  &
                         clConc_SALA, brConc_SALA, 1,                        &
                          H                                                 )

      HET(ind_ClNO2, 1) = kIIR1Ltd( C(ind_ClNO2), C(ind_SALACL), kITemp     )

      ! First consider reaction in troposphere cloud
      kITemp = HETClNO2_TCld( XDenA,       rLiq,        rIce,                &
                              ALiq,        AIce,        VAir,                &
                              Temp,        CldFr,       pHCloud,             &
                              clConc_CldA, clConc_CldC, clConc_Cldg,         &
                              brConc_CldA, brConc_CldC, brConc_Cldg,         &
                              2,           H                                )
      ! not on SALC
      HET(ind_ClNO2, 2) = kIIR1Ltd( C(ind_ClNO2), C(ind_SALCCL), kITemp     )

      !----------------------------------------------------------------
      ! ClNO2 + dissolved HCl in cloud
      !----------------------------------------------------------------
      kITemp = HETClNO2_TCld( XDenA,       rLiq,        rIce,                &
                              ALiq,        AIce,        VAir,                &
                              Temp,        CldFr,       pHCloud,             &
                              clConc_CldA, clConc_CldC, clConc_Cldg,         &
                              brConc_CldA, brConc_CldC, brConc_Cldg,         &
                              3,           H                                )

      HET(ind_ClNO2, 3) = kIIR1Ltd( C(ind_ClNO2), C(ind_HCl), kITemp        )

      !----------------------------------------------------------------
      ! Reaction of ClNO2 with sea-salt Br-, XW 8/8/18)
      !----------------------------------------------------------------

      ! First consider reaction in troposphere cloud
      kITemp = HETClNO2_TCld( XDenA,       rLiq,        rIce,                &
                              ALiq,        AIce,        VAir,                &
                              Temp,        CldFr,       pHCloud,             &
                              clConc_CldA, clConc_CldC, clConc_Cldg,         &
                              brConc_CldA, brConc_CldC, brConc_Cldg,         &
                              4,           H                                )

      ! Then calculate reaction on aerosols
      kITemp = kITemp +                                                      &
               HETClNO2( XDenA,       AClRADI,     (1-CldFr)*AClAREA,        &
                         SSAlk(1),    Temp,        pHSSA(1),                 &
                         clConc_SALA, brConc_SALA, 2,                        &
                         H                                                  )

      HET(ind_ClNO2, 4) = kIIR1Ltd( C(ind_ClNO2), C(ind_BrSALA), kITemp     )

      ! First consider reaction in troposphere cloud
      kITemp = HETClNO2_TCld( XDenA,       rLiq,        rIce,                &
                              ALiq,        AIce,        VAir,                &
                              Temp,        CldFr,       pHCloud,             &
                              clConc_CldA, clConc_CldC, clConc_Cldg,         &
                              brConc_CldA, brConc_CldC, brConc_Cldg,         &
                              5,           H                                )

      ! Then calculate reaction on aerosols
      kITemp = kITemp +                                                      &
               HETClNO2( XDenA,       xRadi(12),   (1-CldFr)*xArea(12),      &
                         SSAlk(2),    Temp,        pHSSA(2),                 &
                         clConc_SALC, brConc_SALC, 2,                        &
                         H                                                  )

      HET(ind_ClNO2, 5) = kIIR1Ltd( C(ind_ClNO2), C(ind_BrSALC), kITemp     )

      !----------------------------------------------------------------
      ! ClNO2 + dissolved HBr in cloud
      !----------------------------------------------------------------
      kITemp = HETClNO2_TCld( XDenA,       rLiq,        rIce,                &
                              ALiq,        AIce,        VAir,                &
                              Temp,        CldFr,       pHCloud,             &
                              clConc_CldA, clConc_CldC, clConc_Cldg,         &
                              brConc_CldA, brConc_CldC, brConc_Cldg,         &
                              6,           H                                )

      HET(ind_ClNO2, 6) = kIIR1Ltd( C(ind_ClNO2), C(ind_HBr), kITemp        )

      !=====================================================================
      ! Iodine chemistry
      !=====================================================================

      ! Uptake reactions (forming AERI, ISALA and ISALC)
      HET(ind_HI,   1) = HETIUptake( H%HI%MW_g,   0.10_fp,  8, Input_Opt )
      HET(ind_HI,   2) = HETIUptake( H%HI%MW_g,   0.10_fp, 11, Input_Opt )
      HET(ind_HI,   3) = HETIUptake( H%HI%MW_g,   0.10_fp, 12, Input_Opt )
      HET(ind_I2O2, 1) = HETIUptake( H%I2O2%MW_g, 0.02_fp,  8, Input_Opt )
      HET(ind_I2O2, 2) = HETIUptake( H%I2O2%MW_g, 0.02_fp, 11, Input_Opt )
      HET(ind_I2O2, 3) = HETIUptake( H%I2O2%MW_g, 0.02_fp, 12, Input_Opt )
      HET(ind_I2O3, 1) = HETIUptake( H%I2O3%MW_g, 0.02_fp,  8, Input_Opt )
      HET(ind_I2O3, 2) = HETIUptake( H%I2O3%MW_g, 0.02_fp, 11, Input_Opt )
      HET(ind_I2O3, 3) = HETIUptake( H%I2O3%MW_g, 0.02_fp, 12, Input_Opt )
      HET(ind_I2O4, 1) = HETIUptake( H%I2O4%MW_g, 0.02_fp,  8, Input_Opt )
      HET(ind_I2O4, 2) = HETIUptake( H%I2O4%MW_g, 0.02_fp, 11, Input_Opt )
      HET(ind_I2O4, 3) = HETIUptake( H%I2O4%MW_g, 0.02_fp, 12, Input_Opt )

      ! These uptake reactions require non-acidic aerosol
      ! Fine sea salt first
      IF (SSAlk(1).gt.0.05) THEN
         HET(ind_HOI,  1) = HETIUptake( H%HOI%MW_g,   0.01_fp, 11, Input_Opt )
         HET(ind_IONO, 1) = HETIUptake( H%IONO%MW_g,  0.02_fp, 11, Input_Opt )
         HET(ind_IONO2,1) = HETIUptake( H%IONO2%MW_g, 0.01_fp, 11, Input_Opt )
      ENDIF

      ! Now coarse sea salt
      IF (SSAlk(2).gt.0.05) THEN
         HET(ind_HOI,  2) = HETIUptake( H%HOI%MW_g,   0.01_fp, 12, Input_Opt )
         HET(ind_IONO, 2) = HETIUptake( H%IONO%MW_g,  0.02_fp, 12, Input_Opt )
         HET(ind_IONO2,2) = HETIUptake( H%IONO2%MW_g, 0.01_fp, 12, Input_Opt )
      ENDIF

      ! Breakdown of iodine compounds on sea-salt
      kITemp           = HETIXCycleSSA( H%HOI%MW_g, 0.01_fp, SSAlk, 1       )
      HET(ind_HOI,  3) = kIIR1Ltd( C(ind_HOI), C(ind_BrSALA), kITemp        )

      kITemp           = HETIXCycleSSA( H%HOI%MW_g, 0.01_fp, SSAlk, 2       )
      HET(ind_HOI,  4) = kIIR1Ltd( C(ind_HOI), C(ind_BrSALC), kITemp        )

      kITemp           = HETIXCycleSSA( H%HOI%MW_g, 0.01_fp, SSAlk, 3       )
      HET(ind_HOI,  5) = kIIR1Ltd( C(ind_HOI), C(ind_SALACL), kITemp      )

      kITemp           = HETIXCycleSSA( H%HOI%MW_g, 0.01_fp, SSAlk, 4       )
      HET(ind_HOI,  6) = kIIR1Ltd( C(ind_HOI), C(ind_SALCCL), kITemp        )

      kITemp           = HETIXCycleSSA( H%IONO%MW_g, 0.02_fp, SSAlk, 1      )
      HET(ind_IONO, 3) = kIIR1Ltd( C(ind_IONO), C(ind_BrSALA), kITemp       )

      kITemp           = HETIXCycleSSA( H%IONO%MW_g, 0.02_fp, SSAlk, 2      )
      HET(ind_IONO, 4) = kIIR1Ltd( C(ind_IONO), C(ind_BrSALC), kITemp       )

      kITemp           = HETIXCycleSSA( H%IONO%MW_g, 0.02_fp, SSAlk, 3      )
      HET(ind_IONO, 5) = kIIR1Ltd( C(ind_IONO), C(ind_SALACL), kITemp       )

      kITemp           = HETIXCycleSSA( H%IONO%MW_g, 0.02_fp, SSAlk, 4      )
      HET(ind_IONO, 6) = kIIR1Ltd( C(ind_IONO), C(ind_SALCCL), kITemp       )

      kITemp            = HETIXCycleSSA( H%IONO2%MW_g, 0.01_fp, SSAlk, 1    )
      HET(ind_IONO2, 3) = kIIR1Ltd( C(ind_IONO2), C(ind_BrSALA), kITemp     )

      kITemp            = HETIXCycleSSA( H%IONO2%MW_g, 0.01_fp, SSAlk, 2    )
      HET(ind_IONO2, 4) = kIIR1Ltd( C(ind_IONO2), C(ind_BrSALC), kITemp     )

      kITemp            = HETIXCycleSSA( H%IONO2%MW_g, 0.01_fp, SSAlk, 3    )
      HET(ind_IONO2, 5) = kIIR1Ltd( C(ind_IONO2), C(ind_SALACL), kITemp     )

      kITemp            = HETIXCycleSSA( H%IONO2%MW_g, 0.01_fp, SSAlk, 4    )
      HET(ind_IONO2, 6) = kIIR1Ltd( C(ind_IONO2), C(ind_SALCCL), kITemp     )

      ! Nullify pointers
      H      => NULL()

    END SUBROUTINE SET_HET
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CloudHet
!
! !DESCRIPTION: Function CloudHet calculates the loss frequency (1/s) of gas
!  species due to heterogeneous chemistry on clouds in a partially cloudy grid
!  cell. The function uses the "entrainment limited uptake" equations of
!  Holmes et al. (2019). Both liquid and ice water clouds are treated.
!
!  For gasses that are that are consumed in multiple aqueous reactions with different
!  products, CloudHet can provide the loss frequency for each reaction branch using the
!  optional branch ratios (branchLiq, branchIce) as arguments.
!
!  Holmes, C.D., Bertram, T. H., Confer, K. L., Ronan, A. C., Wirks, C. K., Graham, K. A.,
!    Shah, V. (2019) The role of clouds in the tropospheric NOx cycle: a new modeling
!    approach for cloud chemistry and its global implications, Geophys. Res. Lett. 46,
!    4980-4990, https://doi.org/10.1029/2019GL081990
!\\
!\\
! !INTERFACE:
!
    FUNCTION CloudHet( fc,    Aliq,      Aice,      rLiq,      rIce,      &
                       T,     airNumDen, molmass_g, gammaLiq,  gammaIce,  &
                              branchLiq, branchIce     ) &
                       RESULT( kHet )
!
! !INPUT PARAMETERS:
!
      REAL(fp),         INTENT(IN) :: fc         ! Cloud Fraction [0-1]
      REAL(fp),         INTENT(IN) :: Aliq       ! Surface area density of
      REAL(fp),         INTENT(IN) :: AIce       !  cloud liquid & ice, cm2/cm3
                                                 !  (grid average, not in-cloud)
      REAL(fp),         INTENT(IN) :: rLiq       ! Effective radius for liquid
      REAL(fp),         INTENT(IN) :: rIce       !  and ice clouds, cm
      REAL(fp),         INTENT(IN) :: T          ! Temperature, K
      REAL(fp),         INTENT(IN) :: airNumDen  ! Air number density,
                                                 !  molec/cm3
      REAL(fp),         INTENT(IN) :: molmass_g  ! molar weight, g/mole
      REAL(fp),         INTENT(IN) :: gammaLiq   ! Gamma values on liquid water
      REAL(fp),         INTENT(IN) :: gammaIce   !   and water ice
      REAL(fp),         OPTIONAL   :: branchLiq  ! Fraction of reactant consumed in a particular reaction branch
      REAL(fp),         OPTIONAL   :: branchIce  !   in liquid and ice, fraction [0-1]
!
! !RETURN VALUE:
!
      REAL(fp)                     :: kHet ! Grid-average loss frequency, 1/s
!
! !REVISION HISTORY:
!  23 Aug 2018 - C. D. Holmes - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Residence time of air in clouds, s
      real(fp), parameter :: tauc = 3600
!
! !LOCAL VARIABLES:
!
      real(fp) :: kI, gam, rd, area
      real(fp) :: kk, ff, xx, branch, kIb, ktmp
      integer  :: K
!
!------------------------------------------------------------------------------
!
      ! If cloud fraction < 0.0001 (0.01%) or there is zero cloud surface area,
      ! then return zero uptake
      if ( (fc < 0.0001) .or. ((ALiq + AIce) <= 0) ) then
         kHet = 0.0_fp
         return
      endif

      !------------------------------------------------------------------------
      ! Loss frequency inside cloud
      !
      ! Assume both water and ice phases are inside the same cloud, so mass
      ! transport to both phases works in parallel (additive)
      !------------------------------------------------------------------------

      ! initialize loss, 1/s
      kI  = 0 ! total loss rate of a gas in cloud
      kIb = 0 ! loss rate for specific reaction branch

      ! Loop over water (K=1) and ice (K=2)
      do K=1, 2

         ! initialize reactions partitioning factor
         branch = 1.0e0_fp

         ! Liquid water cloud
         if (K==1) then

            ! gamma, unitless
            gam = gammaLiq

            ! Liquid particle radius, cm
            rd = rLiq

            ! Liquid surface area density, cm2/cm3
            area = Aliq

            if ( present(branchLiq) ) branch = branchLiq

         ! Ice water cloud
         elseif (K==2) then

            ! gamma, unitless
            gam = gammaIce

            ! Ice particle effective radius, cm
            rd = rIce

            ! Ice surface area density, cm2/cm3
            area = Aice

            if ( present(branchIce) ) branch = branchIce

         else

            print*, 'CloudHet: index value exceeded'
            call GEOS_CHEM_STOP

         endif

         ! Skip calculation if there is no surface area
         if (area <= 0) cycle

         ! Convert grid-average cloud condensate surface area density to in-cloud surface area density
         area = safe_div(area, fc, 0.e+0_fp)

         ! In-cloud loss frequency, combining ice and liquid in parallel, 1/s
         ! Pass radius in cm and mass in g.
         ktmp = arsl1k( area, rd, airnumden, gam, sqrt(T), sqrt(molmass_g))
         kI   = kI + ktmp

         ! In-cloud loss frequency for particular reaction branch, 1/s
         kIb = kIb + ktmp * branch

      end do

      ! Mean branch ratio for reaction of interest in cloud (averaged over ice and liquid)
      branch = kIb / kI

!      !------------------------------------------------------------------------
!      ! Grid-average loss frequency; Add in-cloud and entrainment rates in series
!      !
!      ! APPROXIMATE expression for entrainment-limited uptake
!      !   Approximation error in loss frequency is typically <2% and always <50%.
!      !------------------------------------------------------------------------
!
!      ! Requires declaring...
!      ! real(fp) :: kIinv, kEinv
!
!      ! Entrainment rate, inverse, s
!      kEinv = safe_div( tauc * ( 1e+0_fp - fc ), fc, 1e+30_fp )
!
!      ! In-cloud loss rate, inverse, s
!      kIinv = safe_div( 1e+0_fp, fc*kI, 1e+30_fp )
!
!      ! Overall heterogeneous loss rate, grid average, 1/s
!      kHet = safe_div( 1e+0_fp, ( kEinv + kIinv ), 0e+0_fp )
!
      !------------------------------------------------------------------------
      ! Grid-average loss frequency
      !
      ! EXACT expression for entrainment-limited uptake
      !------------------------------------------------------------------------

      ! Ratio (in cloud) of heterogeneous loss to detrainment, s/s
      kk = kI * tauc

      ! Ratio of volume inside to outside cloud
      ! ff has a range [0,+inf], so cap it at 1e30
      ff = safe_div( fc, (1e0_fp - fc), 1e30_fp )
      ff = max( ff, 1e30_fp )

      ! Ratio of mass inside to outside cloud
      ! xx has range [0,+inf], but ff is capped at 1e30, so this shouldn't overflow
      xx =     ( ff - kk - 1e0_fp ) / 2e0_fp + &
           sqrt( 1e0_fp + ff**2 + kk**2 + 2*ff**2 + 2*kk**2 - 2*ff*kk ) / 2e0_fp

      ! Overall heterogeneous loss rate, grid average, 1/s
      ! kHet = kI * xx / ( 1d0 + xx )
      !  Since the expression ( xx / (1+xx) ) may behave badly when xx>>1,
      !  use the equivalent 1 / (1 + 1/x) with an upper bound on 1/x
      kHet = kI / ( 1e0_fp + safe_div( 1e0_fp, xx, 1e30_fp ) )

      ! Overall loss rate in a particular reaction branch, 1/s
      kHet = kHet * branch
      ! Note: CloudHet currently requires calling CloudHet N times for N reaction branches.
      ! Returning both total loss frequency and mean branch ratio would allow calculation of
      ! all loss rates with N-1 calls to CloudHet (more efficient) (C.D. Holmes)

    end function CloudHet
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
    FUNCTION kIIR1Ltd( concGas, concEduct, kISource ) RESULT( kII )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(f8), INTENT(IN) :: concGas
      REAL(f8), INTENT(IN) :: concEduct
      REAL(dp), INTENT(IN) :: kISource
!
! !RETURN VALUE:
!
      REAL(dp)             :: kII
!
! !REMARKS:
!  Uses HetMinLife and HetMinRate constants from gckpp_Global.F90
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(dp) :: kIGas, kIEduct, lifeA, lifeB

      ! Copy kI as calculated assuming no limitation
      kIGas   = kISource
      kIEduct = 0.0_dp
      kII     = 0.0_dp

      IF ( concEduct < 100.0_dp ) THEN
         kIGas   = 0.0_dp
         kIEduct = 0.0_dp
         kII     = 0.0_dp
      ELSE
         ! Safe division here is probably overkill - may remove this
         IF ( Is_Safe_Div( concGas*kIGas, concEduct ) ) THEN
            kIEduct = kIGas *concGas / concEduct
            kII     = kIGas          / concEduct
         ELSE
            kIGas   = 0.0_dp
            kIEduct = 0.0_dp
            kII     = 0.0_dp
         ENDIF
      ENDIF

      ! Enforce a minimum lifetime?
      IF ( kIGas > 0.0_dp ) THEN

         ! Calculate lifetime of each reactant against removal
         lifeA = Safe_Div( 1.0_dp, kIGas,   0.0_dp )
         lifeB = Safe_Div( 1.0_dp, kIEduct, 0.0_dp )

         ! Check if either lifetime is "too short"
         IF ( ( lifeA < lifeB ) .and. ( lifeA < HetMinLife ) ) THEN
            IF ( Is_Safe_Div( concGas*kIGas, concEduct) ) THEN
               kIGas = HetMinRate
               kII   = kIGas  / concEduct
            ELSE
               kIGas = 0.0_dp
               kII   = 0.0_dp
            ENDIF
         ELSE IF ( lifeB < HetMinLife ) THEN
            IF ( Is_Safe_Div( concEduct*kIEduct, concGas ) ) THEN
               kIEduct = HetMinRate
               kII     = kIEduct / concGas
            ELSE
               kIEduct = 0.0_dp
               kII     = 0.0_dp
            ENDIF
         ENDIF
      ENDIF

    END FUNCTION kIIR1Ltd
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
! !USES:
!
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN)    :: spcVec(:)
      INTEGER,  INTENT(IN)    :: indGasA
      INTEGER,  INTENT(IN)    :: indGasB
      REAL(fp), INTENT(IN)    :: kIASource
      REAL(fp), INTENT(IN)    :: kIBSource
!
! !RETURN VALUE:
!
      REAL(fp)                :: kII

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
      REAL(fp) :: concGasA, concGasB, kIA, kIB
      REAL(fp) :: R_GasA, R_GasB
      LOGICAL  :: nonZeroRate

      ! Get the base concentrations
      concGasA = spcVec(indGasA)
      concGasB = spcVec(indGasB)

      ! Copy the first estimates of kI for each species
      kIA = kIASource
      kIB = kIBSource

      ! Assume for now that the reaction will not proceed
      nonZeroRate = .False.
      kII = 0.0e+0_fp

      ! Prevent reaction if either concentration is too low
      IF ((concGasA.gt.100.0e+0_fp).and.(concGasB.gt.100.0e+0_fp).and.&
          (kIA.gt.0.0e+0_fp).and.(kIB.gt.0.0e+0_fp)) THEN
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
         kII = 0.0e+0_fp
      ENDIF

    END FUNCTION kIIR1R2L
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetIXCycleSSA
!
! !DESCRIPTION: Set the iodine reaction rate on sea salt, assuming a fixed ratio
! of ICl and IBr (85:15) is produced.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETIXCycleSSA( A, B, SSAlk, X ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
      ! Sea salt alkalinity
      REAL(fp), INTENT(IN) :: SSAlk(2)
      ! Br or Cl
      INTEGER, INTENT(IN) :: X
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  24 Dec 2016 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: XSTKCF, ADJUSTEDRATE, r_gp
      INTEGER  :: N, NAer, N1, N2

      ! Initialize
      kISum        = 0.0_fp

      IF (X == 1) THEN
         !for BrSALA
         N1 = 1
         N2 = 1
         r_gp = 0.15e+0_fp
      ELSEIF (X == 2) THEN
         !for BrSALC
         N1 = 2
         N2 = 2
         r_gp = 0.15e+0_fp
      ELSEIF (X == 3) THEN
         !for SALACl
         N1 = 1
         N2 = 1
         r_gp = 0.85e+0_fp
      ELSEIF (X == 4) THEN
         !for SALCCl
         N1 = 2
         N2 = 2
         r_gp = 0.85e+0_fp
      ENDIF


      DO N = N1, N2
         NAer = N + 10
         ! Only allow reaction on acidic aerosol
         IF (SSAlk(N).le.0.05e+0_fp) THEN
            ! Reaction rate for surface of aerosol
            AdjustedRate = ARSL1K(XAREA(NAer),XRADI(NAer),XDENA,B,XTEMP,(A**0.5_fp))
            kISum = kISum + AdjustedRate * r_gp
         ENDIF
      ENDDO

    END FUNCTION HETIXCycleSSA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetIUptake
!
! !DESCRIPTION: Set the uptake rate for iodine species.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETIUptake( A, B, N, Input_Opt ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: A, B       ! Rate coefficients
      INTEGER,        INTENT(IN) :: N          ! Which aerosol?
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  24 Dec 2016 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      !REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      kISum        = 0.0_fp

      ! Reaction rate for surface of aerosol
      kISum = ARSL1K(XAREA(N),XRADI(N),XDENA,B,XTEMP,(A**0.5_fp))

      IF ( Input_Opt%LUCX ) THEN
         ! For UCX-based mechanisms also allow reaction on stratospheric
         ! sulfate (N=13) if tropospheric sulfate is requested (N=8)
         IF (N.eq.8) THEN
            kISum = kISum + ARSL1K(XAREA(13),XRADI(13),XDENA,B,XTEMP, &
                    (A**0.5_fp))
         ENDIF
      ENDIF

    END FUNCTION HETIUptake
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetNO3
!
! !DESCRIPTION: Set the heterogenous chemistry rate for NO3.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETNO3( A, B ) RESULT( HET_NO3 )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_NO3
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
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_NO3      = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         XSTKCF = B

         ! Set uptake coefficients
         select case (N)
            case (1:7)
               ! dust
               xstkcf = 0.01_fp
            case (8)
               ! sulfate
               !if ( relhum < 40.0_fp ) then
               !   xstkcf = 0.001_fp
               !else
               !   xstkcf = 0.002_fp
               !endif
               ! Follow NO3 hydrolysis on SSA+SNA
               xstkcf = 0.0_fp
            case (9)
               ! BC
               if ( relhum < 50.0_fp ) then
                  xstkcf = 2e-4_fp
               else
                  xstkcf = 1e-3_fp
               endif
            case (10)
               ! OC
               xstkcf = 0.005_fp
            case (11:12)
               ! sea salt
               !if ( relhum < 40.0_fp ) then
               !   xstkcf = 0.05_fp
               !elseif ( relhum > 70.0_fp ) then
               !   xstkcf = 0.002_fp
               !else
               !   xstkcf = 0.05_fp + (0.002_fp - 0.05_fp)                    &
               !          * (relhum - 40.0_fp) / (70.0_fp - 40.0_fp)
               !endif
               ! Follow NO3 hydrolysis on SSA+SNA
               xstkcf = 0.0_fp
         end select

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSEIF ((N.eq.11).or.(N.eq.12).or.(N.eq.8)) THEN
            ! Follow combined hydrolysis with Cl-, xnw
            ADJUSTEDRATE = 0.0e0_fp
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_NO3 = HET_NO3 + ADJUSTEDRATE

      ENDDO

    END FUNCTION HETNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetNO2
!
! !DESCRIPTION: Set the heterogenous chemistry rate for NO2.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETNO2( A, B ) RESULT( HET_NO2 )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_NO2
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
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_NO2      = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         XSTKCF = B

         ! Set uptake coefficients
         select case (N)
            case (1:7)
               ! dust
               xstkcf = 1e-8_fp
            case (8)
               ! sulfate
               xstkcf = 5e-6_fp
            case (9)
               ! BC
               xstkcf = 1e-4_fp
            case (10)
               ! OC
               xstkcf = 1e-6_fp
            case (11:12)
               ! sea salt
               if ( relhum < 40 ) then
                  xstkcf = 1e-8_fp
               elseif ( relhum > 70 ) then
                  xstkcf = 1e-4_fp
               else
                  xstkcf = 1e-8_fp + (1e-4_fp-1e-8_fp) * (relhum-40)/(70-40)
               endif
         end select

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_fp))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_NO2 = HET_NO2 + ADJUSTEDRATE

      END DO

    END FUNCTION HETNO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHO2
!
! !DESCRIPTION: Set the heterogenous chemistry rate for HO2.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHO2( A, B ) RESULT( HET_HO2 )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_HO2
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
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_HO2      = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         XSTKCF = 0e+0_fp
         IF (N.gt.12) THEN
            XSTKCF = TINY(1e+0_fp)
         ELSE
            XSTKCF = GAMMA_HO2
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_HO2 = HET_HO2 + ADJUSTEDRATE

      ENDDO

    END FUNCTION HETHO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHBr
!
! !DESCRIPTION: Set the heterogeneous rate for HBr.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHBr( A, B ) RESULT( HET_HBr )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_HBr
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
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_HBr      = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Only apply PSC rate adjustment if at high altitude
      DO_EDUCT     = STRATBOX

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! jpp, 3/22/11: set the sticking coefficient to
         !  ~0 for aerosol types we don't want reactions on
         !  for the HBr and HOBr surface reaction
         XSTKCF = 0e+0_fp
         ! Select proper aerosol type
         IF ( (N == 8) .OR. (N == 11) .OR. (N == 12)) THEN
            ! sulfate, 2 modes of sea-salt
            XSTKCF = B
         ELSEIF ( N == 13 ) THEN
            XSTKCF = KHETI_SLA(11)
         ELSEIF ( N == 14 ) THEN
            XSTKCF = 0.1e+0_fp
         ELSE
            XSTKCF = 0e+0_fp
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_HBr = HET_HBr + ADJUSTEDRATE

      ENDDO

    END FUNCTION HETHBr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetN2O5
!
! !DESCRIPTION: Set heterogenous chemistry rate for N2O5.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETN2O5( A, B, CldFr) RESULT( HET_N2O5 )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
      ! Cloud fraction
      REAL(fp), INTENT(IN) :: CldFr
!
! !RETURN VALUE:
!
      ! HET_N2O5(1) = rate coefficient; HET_N2O5(2) = SA-weighted gamma
      REAL(fp)             :: HET_N2O5(2)
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
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      REAL(fp) :: output(4)
      REAL(fp) :: ClNO2_yield, kClNO2, Rp, SA, SA_sum
      REAL(fp) :: TMP1,   TMP2
      REAL(fp) :: GAM_N2O5, r_gp

!------------------------------------------------------------------------------
      ! Initialize
      HET_N2O5      = 0.0_fp    !Output, holds ADJUSTEDRATE and XSTCKF results
      ADJUSTEDRATE  = 0.0_fp    !Total N2O5 loss rate coefficient for all aerosol types
      XSTKCF        = 0.0_fp    !Total N2O5 gamma for all aerosol types
      ClNO2_yield   = 0.0_fp    !ClNO2 production yield
      kClNO2        = 0.0_fp
      output        = 0.0_fp    !temporary variable
      Rp            = 0.0_fp    !Total SNA+ORG radius
      SA            = 0.0_fp    !Total SNA+ORG surface area
      SA_sum        = 0.0_fp    ! Used for weighted mean
      TMP1         = 0.0_fp
      TMP2         = 0.0_fp
      GAM_N2O5     = 0.0_fp

      ! Always apply PSC rate adjustment
      DO_EDUCT     = .TRUE.

      !NOTE: This calculations follows these general steps:
      !Loop over all aerosol types
      ! A) Calculate gamma for each aerosol type
      ! B) Use gamma to calculate N2O5 loss rate coefficient for each aerosol
      !    type
      ! C) Add results from all aerosol types to get running sum of ADJUSTEDRATE
      !    and surface-area-weighted gamma
      !
      !After looping through all aerosol types:
      ! D) Divide gamma by SA to get single SA-weighted gamma

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE
         ! A) Calculate gamma for each aerosol type

         XSTKCF = 0e+0_fp
         ! Get GAMMA for N2O5 hydrolysis, which is
         ! a function of aerosol type, temp, and RH
         IF (N.eq.14) THEN
            IF (NATSURFACE) THEN
               XSTKCF = 4.0e-4_fp ! NAT
            ELSE
               XSTKCF = 0.02e+0_fp ! Ice
            ENDIF
         ELSEIF (N.eq.13) THEN
            ! Stratospheric aerosol
            XSTKCF = KHETI_SLA(1)
         ELSEIF ((N==8) .or. (N==10) .or. (N==11) .or. (N==12)) THEN
            XSTKCF = 0.0e+0_fp
         ELSE
            ! For UCX-based mechanisms ABSHUMK is set to Spc(I,J,L,id_H2O)
            XSTKCF = N2O5( N, TEMP, RELHUM )
         ENDIF

         ! B) Use gamma to calculate N2O5 loss rate coefficient for each
         !    aerosol type
         IF (N.eq.8) THEN

            ! Properties of inorganic (sulfate-nitrate-ammonium-sea salt) coated
            ! with organics
            output = N2O5_InorgOrg( AClVOL, XVOL(10), XH2O(8), XH2O(10), &
                 AClRADI,  C(ind_NIT), C(ind_SALACL),  TEMP,    RELHUM )

            ! Gamma
            XSTKCF = output(1)
            ! ClNO2 yield, fraction [0,1]
            ! Reduce ClNO2 production yield on inorganic+organic aerosol by 75%,
            ! this is based on the analysis of WINTER observations in McDuffie et al., JGR, 2018
            ClNO2_yield = output(2)*0.25_fp
            ! Particle radius with coating, cm
            Rp     = output(3)
            ! Surface area of coated particles, cm2/cm3
            SA     = output(4)

            !For SNA aerosol...
            ! Total loss rate of N2O5 (kN2O5) on SNA+ORG+fineSSA aerosol
            ADJUSTEDRATE=ARSL1K( (1.0_fp-CldFr)*SA, Rp, XDENA, XSTKCF, XTEMP, (A**0.5_FP) )

            !Calculate ClNO2 yield on SNA(+ORG) aerosol using kN2O5 from above
            ! phi = kClNO2/kN2O5 (from ClNO2 and HNO3 production pathways)
            kClNO2 = ClNO2_yield * ADJUSTEDRATE

            ! reduce the rate of this HNO3 pathway in accordance with the yield
            ADJUSTEDRATE = ADJUSTEDRATE - kClNO2
         ELSEIF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSEIF (N.eq.12) THEN

            ! Properties of coarse mode SSA
            output = N2O5_InorgOrg( XVOL(12), 0e+0_fp, XH2O(12), 0e+0_fp, &
                 XRADI(12), C(ind_NITs), C(ind_SALCCL), TEMP, RELHUM )

            ! Gamma
            XSTKCF = output(1)
            ! ClNO2 yield, fraction [0,1]
            ClNO2_yield = output(2)
            ! Particle radius with coating, cm
            Rp     = output(3)
            ! Surface area of coated particles, cm2/cm3
            SA     = output(4)

            !For coarse mode SSA aerosol ...
            ADJUSTEDRATE=ARSL1K( (1.0_fp-CldFr)*SA, Rp, XDENA, XSTKCF, XTEMP, (A**0.5_FP))

            !Calculate ClNO2 yield using kN2O5 from above
            ! phi = kClNO2/kN2O5 (from ClNO2 and HNO3 production pathways)
            kClNO2 = ClNO2_yield * ADJUSTEDRATE

            ! reduce the rate of this HNO3 pathway in accordance with the yield
            ADJUSTEDRATE = ADJUSTEDRATE - kClNO2

         ELSE
            !For all other aerosol types...
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K((1-CldFr)*XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         ! C) Add results from all aerosol types to get running sum of
         !    ADJUSTEDRATE and surface-area-weighted gamma
         HET_N2O5(1) = HET_N2O5(1) + ADJUSTEDRATE         !loss rate coefficient

         IF (N.eq.8) THEN
              HET_N2O5(2) = HET_N2O5(2) + ( XSTKCF * SA ) !SA-weighted gamma
              SA_sum = SA_sum + SA                        !total SA in cm2/cm3
         ELSEIF ( (N.eq.10) .or. (N.eq.11) )THEN
              !don't include ORG and fine SSA SA since it has already been included above
              SA_sum = SA_sum
         ELSE
              HET_N2O5(2) = HET_N2O5(2) + ( XSTKCF * XAREA(N))!SA-weighted gamma
              SA_sum = SA_sum + XAREA(N)                      !total SA
         ENDIF

      END DO

      ! D) Divide gamma by SA to get SA-weighted gamma
      HET_N2O5(2) = safe_div( HET_N2O5(2), SA_sum, 0.0e+0_fp )

    END FUNCTION HETN2O5
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: N2O5_InorgOrg
!
! !DESCRIPTION: Function N2O5\_inorg\_org computes the GAMMA reaction probability
!     for N2O5 loss in inorganic (sulfate-nitrate-ammonium-sea salt) aerosols
!     with organic coatings, based on the recommendation of McDuffie (2018) JGR.
!     The inorganic core is based on Bertram and Thornton ACP (2009).
!\\
!\\
! !INTERFACE:
!
    FUNCTION N2O5_InorgOrg( volInorg, volOrg, H2Oinorg, H2Oorg, Rcore, &
                            NIT,      Cl,     T,        RH ) &
         RESULT ( output )
!
! !USES:
!
  USE PhysConstants,      ONLY : AVO, RGASLATM, RSTARG, PI
!
! !INPUT PARAMETERS:
!
      real(fp), intent(in) :: volInorg ! volume of wet inorganic aerosol core
                                       !  [cm3(aerosol)/cm3(air)]
      real(fp), intent(in) :: volOrg   ! volume of wet organic aerosol coating
                                       !  [cm3(aerosol)/cm3(air)]
      real(fp), intent(in) :: H2Oinorg ! volume of H2O in inorganic core
                                       !  [cm3(H2O)/cm3(air)]
      real(fp), intent(in) :: H2Oorg   ! volume of H2O in organic coating
                                       !  [cm3(H2O)/cm3(air)]
      real(fp), intent(in) :: Rcore    ! radius of inorganic core [cm]
      real(fp), intent(in) :: NIT      ! aerosol nitrate concentration
                                       !  [molecule/cm3(air)
      real(fp), intent(in) :: Cl       ! aerosol chloride concentration
                                       !  [molecule/cm3(air)
      real(fp), intent(in) :: T        ! air temperature [K]
      real(fp), intent(in) :: RH       ! relative humidity [%]
!
! !RETURN VALUE:
!
      real(fp) :: output(4) ! output(1) = gamma;
                            ! output(2) = ClNO2_yield
                            ! output(2) = particle radius, cm
                            ! output(4) = surface area, cm2/cm3
!
! !REMARKS:
!
! !REVISION HISTORY:
!  13 Dec 2018 - E. McDuffie - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Parameters from Bertram and Thornton (2009) ACP and McDuffie 2018 JGR
      real(fp), parameter :: KH    = 5.1e+1_fp   !unitless
      real(fp), parameter :: k3k2b = 4.0e-2_fp   !unitless
      real(fp), parameter :: beta  = 1.15e+6_fp  ![s-1]
      real(fp), parameter :: delta = 1.3e-1_fp   ![M-1]

      ! Organic Parameters from Antilla et al., 2006 and Riemer et al., 2009
      ! aq. Henry's Law coef [mol m-3 atm-1] (Antilla 2006)
      real(fp), parameter :: Haq  = 5e+3_fp
      ! aq. diffusion coef [m2 s-1] (Riemer, 2009)
      real(fp), parameter :: Daq  = 1e-9_fp
      ! Gas constant [m3 atm K-1 mol-1]
      real(fp), parameter :: R    = RGASLATM * 1e-3_fp
!
! !LOCAL VARIABLES
!
      real(fp)            :: k2f   ! [s-1]
      real(fp)            :: A     ! [s]
      real(fp)            :: speed ! mean molecular speed of N2O5 [m/s]

      real(fp) :: gamma, gamma_core, gamma_coat
      real(fp) :: volTotal, H2Ototal, areaTotal, volRatioDry
      real(fp) :: Rp, l, eps, OCratio
      real(fp) :: M_H2O, M_NIT, M_Cl, ClNO2_yield

      !------------------------------------------------------------------------
      ! Particle size & volume, coating thickness, molar concentrations
      !------------------------------------------------------------------------

      ! Total volume (organic + inorganic), cm3(aerosol)/cm3(air)
      volTotal = volInorg + volOrg

      ! Total H2O (organic + inorganic), cm3(H2O)/cm3(air)
      H2Ototal = H2Oinorg + H2Oorg

      ! Ratio of inorganic to total (organic+inorganic) volumes when dry, unitless
      volRatioDry = safe_div(max(volInorg - H2Oinorg, 0.0_fp),               &
                             max(volTotal - H2Ototal, 0.0_fp), 0.0_fp)

      ! Particle radius, cm
      ! Derived from spherical geometry
      ! [note: The radius and surface area of a wet particle are
      ! properly calculated from the wet volume volume ratio (including water).
      ! We use the dry volume ratio here because McDuffie et al. (2018) fitted
      ! the N2O5 gamma parameters to field data in a model using the
      ! dry ratio. cdholmes 7/22/2019]
      Rp = safe_div( Rcore, volRatioDry**(1e+0_fp/3e+0_fp), Rcore)

      ! Coating thickness, cm
      l = Rp - Rcore

      ! mean molecular speed [m s-1]
      ! sqrt( 8RT / (pi M) )
      speed = sqrt( 8e+0_fp * RSTARG * T / ( PI * 0.108+0_fp) )

      ! H2O molar concentration, mol/L
      M_H2O = H2Ototal / 18e+0_fp / volTotal * 1e+3_fp
      ! Nitrate molar concentration, mol/L
      M_NIT = NIT / volTotal / AVO * 1e+3_fp
      ! Chloride molar concentration, mol/L
      M_Cl = Cl   / volTotal / AVO * 1e+3_fp

      !------------------------------------------------------------------------
      ! Gamma for the organic shell
      ! Implements recommendations by McDuffie (2018) JGR,
      !------------------------------------------------------------------------

      !O:C ratio from Eq. 10 of Canagaratna et al., 2015 (ACP)
      ! Take average OM/OC ratio from /GeosCore/aerosol_mod.F90
      OCratio = ( ( ( OMOC_POA + OMOC_OPOA ) / 2 ) - 1.17e+0_fp ) / 1.29e+0_fp

      ! organic scaling factor (eps(Haq*Daq) = Horg*Dorg)
      ! from McDuffie (2018) JGR
      eps = 1.5e-1_fp * OCratio + 1.6e-3_fp * RH

      ! Gamma for coating
      ! [Rcore, Rp, and l converted cm -> m here]
      IF ( l <= 0.0e+0_fp ) THEN
         gamma_coat = 0.0e+0_fp
      ELSE
         gamma_coat = ( 4e+0_fp * R * T * eps * Haq * Daq * Rcore/1e+2_fp ) / &
                      ( speed * l/1e+2_fp * Rp/1e+2_fp )
      ENDIF

      ! Total particle surface area, cm2/cm3
      areaTotal = 3e+0_fp * volTotal / Rp

      !------------------------------------------------------------------------
      ! Gamma for the inorganic core
      ! Implements recommendations by McDuffie (2018) JGR,
      ! following the general approach from Bertram and Thornton ACP (2009).
      !------------------------------------------------------------------------

      ! Select dry or deliquesed aerosol based on molar concentration of H2O
      IF ( M_H2O < 0.1e+0_fp ) THEN

         ! When H2O is nearly zero, use dry aerosol value
         gamma_core = 0.005e+0_fp

      ELSE

         ! mean molecular speed [cm/s]
         speed = speed * 1e+2_fp

         ! A factor from Bertram and Thornton (2009), s
         ! Their paper suggested an approximated value of A = 3.2D-8
         A = ( ( 4 * volTotal ) / ( speed * areaTotal ) ) * KH

         ! k2f - reaction rate constant of N2O5 with H2O
         ! From McDuffie (2018): k2f = 2.14D5 * H2O
         ! This linear water dependence is not accurate at large
         ! (>20 M) aerosol water concentrations. Therefore, k2f is
         ! calculated following Bertram and Thornton ACP (2009).
         ! Eq 11 from Bertram and Thronton (2009):
         ! Modified to avoid underflow when exp(-delta*H2O) ~1
         IF ( delta * M_H2O < 1e-2_fp ) THEN
            k2f = beta * ( delta * M_H2O )
         ELSE
            k2f = beta * ( 1e+0_fp - exp( -delta * M_H2O ) )
         ENDIF

         ! Eq 12 from Bertram and Thornton (2009)
         ! Use safe_div to avoid overflow when NIT ~ 0
         gamma_core = A * k2f * ( 1e+0_fp - &
            1e+0_fp / ( 1e+0_fp + safe_div( k3k2b * M_H2O, M_NIT, 1e+30_fp ) ) )

      ENDIF

      !------------------------------------------------------------------------
      ! Gamma for overall uptake
      !------------------------------------------------------------------------

      IF (gamma_coat <= 0.0e+0_fp ) THEN
         gamma = gamma_core
      ELSEIF (gamma_core <= 0.0e+0_fp ) THEN
         gamma = 0.0e+0_fp
      ELSE
         gamma = 1e+0_fp / (( 1e+0_fp / gamma_core ) + ( 1e+0_fp / gamma_coat ))
      ENDIF

      !------------------------------------------------------------------------
      ! ClNO2 yield
      !------------------------------------------------------------------------

      ! Calculate the ClNO2 yield following Bertram and Thornton 2009 ACP
      ClNO2_yield = ClNO2_BT( M_Cl, M_H2O )

      output(1) = gamma
      output(2) = ClNO2_yield
      output(3) = Rp
      output(4) = areaTotal

    END FUNCTION N2O5_InorgOrg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ClNO2_BT
!
! !DESCRIPTION: Function ClNO2\_BT computes the PHI production yield of ClNO2
!  from N2O5 loss in sulfate-nitrate-ammonium (SNA) aerosols based on the
!  recommendation of Bertram and Thornton (2009) ACP.
!\\
!\\
! !INTERFACE:
!
    FUNCTION ClNO2_BT( Cl, H2O ) RESULT( PHI )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: Cl         !Aerosol chloride, mol/L
      REAL(fp), INTENT(IN) :: H2O        ! Aerosol H2O, mol/L

!
! !RETURN VALUE:
!
      REAL(fp)             :: PHI
!
! !REMARKS:
!
! !REVISION HISTORY:
!  13 Dec 2018 - E. McDuffie - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Parameters from Bertram and Thornton (2009) ACP
      REAL(fp), parameter :: k2k3  = 1e+0_fp / 4.5e+2_fp

      ! Initialize
      PHI     = 0.0_fp

      !When H2O is nearly zero, assign phi accordingly
      IF ( H2O < 0.1 ) THEN
           IF ( Cl > 1e-3_fp ) THEN
               PHI = 1e+0_fp
           ELSE
               PHI = 0e+0_fp
           ENDIF
      ELSE
           ! Eq from Bertram and Thronton (2009)
           ! Use safe_div to avoid overflow when Cl ~ 0
           PHI = 1e+0_fp / ( 1e+0_fp + k2k3 * ( safe_div( H2O, Cl, 1e+30_fp ) ) )

     ENDIF

    END FUNCTION ClNO2_BT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetGLYX
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for GLYX.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETGLYX( A, B ) RESULT( HET_GLYX )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_GLYX
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_GLYX     = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Default value
         XSTKCF = TINY(1e+0_fp)

         ! Only consider inorganic aqueous aerosols with RH > 35%.
         ! Uptake during the day (higher uptake than night)
         ! (Sumner et al., 2014):
         IF ( N == 8 .and. RELHUM >= CRITRH ) THEN

            ! Define gamma for GLYX:
            IF ( SUNCOS .gt. 0 ) THEN

               ! Uptake during the day (use Liggio et al., 2005):
               ! XSTKCF = 2.9e-3_fp ! Prior to 3/2/18
               XSTKCF = 4.4e-3_fp

            ELSE

               ! Uptake at night (lower uptake than day)
               ! Value is within the range 1d-5 to 1d-6
               ! (Faye McNeill personal communication, eam, 2015):
               ! XSTKCF = 5.0e-6_fp ! Prior to 3/2/18
               XSTKCF = 8.0e-6_fp

            ENDIF

         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_GLYX = HET_GLYX + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETGLYX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetMGLY
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for MGLY.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETMGLY( A, B ) RESULT( HET_MGLY )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_MGLY
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_MGLY     = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Default value
         XSTKCF = TINY(1e+0_fp)

         ! Only consider inorganic aqueous aerosols with RH > 35%.
         IF ( N == 8 .and. RELHUM >= CRITRH ) THEN

            ! Define gamma for MGLY:
            ! Obtained by scaling gamma GLYX by the
            ! ratio of effective Henry's law constants
            ! for GLYX (3d7) and MGLY (3.7d3) (eam, 02/2015):
            XSTKCF = 3.6e-7_fp

         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_MGLY = HET_MGLY + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETMGLY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetIEPOX
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for IEPOX.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETIEPOX( A, B ) RESULT( HET_IEPOX )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_IEPOX
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      REAL(fp) :: HSTAR, K_HPLUS, K_NUC, K_HSO4, K_HYDRO

      ! Initialize
      HET_IEPOX    = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Default value
         XSTKCF = TINY(1e+0_fp)

         ! Only consider inorganic aqueous aerosols with RH > 35%.
         IF ( N == 8 .and. RELHUM >= CRITRH ) THEN

            ! Define Henry's Law constant
            ! Changes H* for IEPOX again to accommodate
            ! reduction in yields of RIP, precursor
            ! of IEPOX (eam, 07/2015):
            HSTAR = HSTAR_EPOX    ! (Nguyen et al., 2014)

            ! Define first-order particle phase reaction rates
            ! specific to IEPOX (from Gaston et al., 2014):
            K_HPLUS = 3.6e-2_fp   ! Alternate: 1.2d-3 (Edding)
            K_NUC   = 2.e-4_fp    ! Alternate: 5.2d-1 (Piletic)
            K_HSO4  = 7.3e-4_fp
            K_HYDRO = 0.0e+0_fp

            ! Get GAMMA for IEPOX hydrolysis:
            XSTKCF = EPOXUPTK( XAREA(N), XRADI(N),            &
                               TEMP,     (A**0.5_fp),         &
                               HSTAR,    K_HPLUS,    H_PLUS,  &
                               K_NUC,    MSO4,       MNO3,    &
                               K_HSO4,   MHSO4,      K_HYDRO )

         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_IEPOX = HET_IEPOX + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETIEPOX
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetIMAE
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for IMAE.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETIMAE( A, B ) RESULT( HET_IMAE )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_IMAE
!
! !REMARKS:
! Here use the same values as are read in for IEPOX, but scale down gamma by
! a factor of 30 to get the value for IMAE. Gamma for the two species are
! similar under neutral conditions, so use IEPOX gamma when [H+] <= 8d-5.
! Implemented by (eam, 01/2015) using lab study findings from Riedel et al.,
! EST, 2015.
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      REAL(fp) :: HSTAR, K_HPLUS, K_NUC, K_HSO4, K_HYDRO

      ! Initialize
      HET_IMAE     = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Default value
         XSTKCF = TINY(1e+0_fp)

         ! Only consider inorganic aqueous aerosols with RH > 35%.
         IF ( N == 8 .and. RELHUM >= CRITRH ) THEN

            ! Define Henry's Law constant.
            ! Changes H* for IEPOX again to accommodate
            ! reduction in yields of RIP, precursor
            ! of IEPOX (eam, 07/2015):
            HSTAR = HSTAR_EPOX   ! (Nguyen et al., 2014)

            ! Define first-order particle phase reaction rates
            ! specific to IEPOX (from Gaston et al., 2014):
            K_HPLUS = 3.6e-2_fp   ! Alternate: 1.2d-3 (Edding)
            K_NUC   = 2.6e-4_fp   ! Alternate: 5.2d-1 (Piletic)
            K_HSO4  = 7.3e-4_fp
            K_HYDRO = 0.e+0_fp

            ! Get GAMMA for IMAE hydrolysis:
            XSTKCF = EPOXUPTK( XAREA(N), XRADI(N),            &
                               TEMP,    (A**0.5_fp),         &
                               HSTAR,    K_HPLUS,    H_PLUS,  &
                               K_NUC,    MSO4,       MNO3,    &
                               K_HSO4,   MHSO4,      K_HYDRO )

            ! Scale down gamma if H+ > 8d-5 (30x less than gamma for IEPOX)
            ! (Riedel et al., 2015)
            IF ( H_PLUS .gt. 8.e-5_fp ) THEN
               XSTKCF = XSTKCF / 30.e+0_fp
            ENDIF

         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_IMAE = HET_IMAE + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETIMAE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetLVOC
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for LVOC: condensation of
! low-volatility ISOPOOH oxidation products.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETLVOC( A, B ) RESULT( HET_LVOC )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_LVOC
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_LVOC     = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma

         IF (N.eq.8 .and. RELHUM >= CRITRH) THEN
            XSTKCF = B
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_LVOC = HET_LVOC + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETLVOC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetISN1OG
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for ISN1OG: uptake of 2nd
! generation organic nitrates formed from ISOP+NO3 reaction (eam, 02/2015).
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETISN1OG( A, B ) RESULT( HET_ISN1OG )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_ISN1OG
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_ISN1OG   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         IF (N.eq.8 .and. RELHUM >= CRITRH) THEN
            XSTKCF = B
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_ISN1OG = HET_ISN1OG + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETISN1OG
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetISOPND
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for ISOPND.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETISOPND( A, B ) RESULT( HET_ISOPND )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_ISOPND
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_ISOPND   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         IF (N.eq.8 .and. RELHUM >= CRITRH) THEN
            XSTKCF = B
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_ISOPND = HET_ISOPND + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETISOPND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetISOPNB
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for ISOPNB.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETISOPNB( A, B ) RESULT( HET_ISOPNB )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_ISOPNB
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_ISOPNB   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         IF (N.eq.8 .and. RELHUM >= CRITRH) THEN
            XSTKCF = B
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_ISOPNB = HET_ISOPNB + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETISOPNB
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetMACRN
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for MACRN.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETMACRN( A, B ) RESULT( HET_MACRN )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_MACRN
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_MACRN    = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         IF (N.eq.8 .and. RELHUM >= CRITRH) THEN
            XSTKCF = B
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_MACRN = HET_MACRN + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETMACRN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetMVKN
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for MVKN.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETMVKN( A, B ) RESULT( HET_MVKN )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_MVKN
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_MVKN     = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         IF (N.eq.8 .and. RELHUM >= CRITRH) THEN
            XSTKCF = B
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_MVKN = HET_MVKN + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETMVKN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetR4N2
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for R4N2.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETR4N2( A, B ) RESULT( HET_R4N2 )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_R4N2
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_R4N2   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         IF (N.eq.8 .and. RELHUM >= CRITRH) THEN
            XSTKCF = B
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_R4N2 = HET_R4N2 + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETR4N2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetDHDN
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for DHDN.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETDHDN( A, B ) RESULT( HET_DHDN )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_DHDN
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_DHDN     = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         IF (N.eq.8 .and. RELHUM >= CRITRH) THEN
            XSTKCF = B
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_DHDN = HET_DHDN + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETDHDN
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetMONITS
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for MONITS
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETMONITS( A, B ) RESULT( HET_MONITS )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_MONITS
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Jul 2017 - M. Sulprizio- Initial version based on SEAC4RS code and Fisher
!                              et al. 2016.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_MONITS   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         IF (N.eq.8 .and. RELHUM >= CRITRH) THEN
            XSTKCF = B
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_MONITS = HET_MONITS + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETMONITS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetMONITU
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for MONITU
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETMONITU( A, B ) RESULT( HET_MONITU )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_MONITU
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Jul 2017 - M. Sulprizio- Initial version based on SEAC4RS code and Fisher
!                              et al. 2016.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_MONITU   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         IF (N.eq.8 .and. RELHUM >= CRITRH) THEN
            XSTKCF = B
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_MONITU = HET_MONITU + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETMONITU
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHONIT
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for HONIT
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHONIT( A, B ) RESULT( HET_HONIT )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_HONIT
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Jul 2017 - M. Sulprizio- Initial version based on SEAC4RS code and Fisher
!                              et al. 2016.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_HONIT    = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma
         IF (N.eq.8 .and. RELHUM >= CRITRH) THEN
            XSTKCF = B
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         HET_HONIT = HET_HONIT + ADJUSTEDRATE
      END DO

    END FUNCTIOn HETHONIT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetIONITA
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for IONITA: Aerosol-phase
!  organic nitrate formed from monoterpene precursors.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETIONITA( A, B ) RESULT( HET_IONITA )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_IONITA
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Jul 2017 - M. Sulprizio- Initial version based on SEAC4RS code and Fisher
!                              et al. 2016.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_IONITA   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma for IONITA
         ! Imposed lifetime = 1 hour (Fisher et al., 2016)
         XSTKCF = 2.78e-4_fp

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         !HET_IONITA = HET_IONITA + ADJUSTEDRATE

      END DO
      HET_IONITA = 2.78e-4_fp

    END FUNCTIOn HETIONITA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetMONITA
!
! !DESCRIPTION: Sets the heterogenous chemistry rate for MONITA: Aerosol-phase
!  organic nitrate formed from monoterpene precursors.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETMONITA( A, B ) RESULT( HET_MONITA )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: HET_MONITA
!
! !REMARKS:
!
! !REVISION HISTORY:
!  14 Jul 2017 - M. Sulprizio- Initial version based on SEAC4RS code and Fisher
!                              et al. 2016.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      HET_MONITA   = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp

      ! Don't do PSC rate adjustment
      DO_EDUCT     = .FALSE.

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Define gamma for MONITA
         ! Imposed lifetime = 1 hour (Fisher et al., 2016)
         XSTKCF = 2.78e-4_fp

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               (A**0.5_FP))
         ENDIF

         IF ( DO_EDUCT .and. N > 12 ) THEN
            ! PSC reaction - prevent excessive reaction rate
            IF (ADJUSTEDRATE.gt.(1.e+0_fp/HetMinLife)) THEN
               ADJUSTEDRATE = 1.e+0_fp/HetMinLife
            ENDIF
         ENDIF

         ! Add to overall reaction rate
         !HET_MONITA = HET_MONITA + ADJUSTEDRATE
      END DO

      HET_MONITA = 2.78e-4_fp

    END FUNCTIOn HETMONITA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetOH
!
! !DESCRIPTION: Set heterogenous chemistry rate for N2O5 with Cl-.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETOH(denAir, rAer, AAer, TK, C_X) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp), INTENT(IN) :: rAer        ! Radius of aerosol (cm)
      REAL(fp), INTENT(IN) :: AAer        ! Area of aerosol (cm2/cm3)
      REAL(fp), INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp), INTENT(IN) :: C_X         ! Cl- concentration (mol/L)
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  12 Mar 2018 - X. Wang  - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF
      Real(fp), Parameter :: XMolWeight=17.0e+0_fp
      Real(fp), Parameter :: XSQM=SQRT(XMolWeight)
!
! !DEFINED PARAMETERS:
!
      ! Initialize
      kISum        = 0.0_fp
      XSTKCF       = 0.0_fp

      XSTKCF = 0.04_fp * C_X   !Knipping and Dabdub (2002)

      ! Reaction rate
      kISum = ARSL1K( AAer, rAer, denAir, XSTKCF, XTEMP, XSQM )

    END FUNCTION HETOH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetN2O5_SS
!
! !DESCRIPTION: Set heterogenous chemistry rate for N2O5 on Cl- containing
! aerosols. This reaction follows the N2O5 + Cl- channel.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETN2O5_SS(CldFr, X) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: CldFr       ! Cloud fraction
      INTEGER,  INTENT(IN) :: X           ! X = 1 fine; 2 coarse
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum(3) !(1) = rate coefficient, (2) = gamma, (3) = phi
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
      INTEGER  :: N
      REAL(fp) :: XSTKCF, output(4)
      REAL(fp) :: GAM_N2O5, r_gp
      REAL(fp) :: ClNO2_yield, kClNO2, Rp, SA
      Real(fp), Parameter :: XMolWeight=108.0e+0_fp
      Real(fp), Parameter :: XSQM=SQRT(XMolWeight)
!
! !DEFINED PARAMETERS:
!
      ! Initialize
      kISum        = 0.0_fp
      XSTKCF       = 0.0_fp
      ClNO2_yield   = 0.0_fp    !ClNO2 production yield (for future update)
      kClNO2        = 0.0_fp
      output        = 0.0_fp    !temporary variable
      Rp            = 0.0_fp    !Total SNA+ORG radius
      SA            = 0.0_fp    !Total SNA+ORG surface area

      IF (X .EQ. 1) THEN
         ! Properties of inorganic (sulfate-nitrate-ammonium-sea salt) coated
         ! with organics
         output = N2O5_InorgOrg( AClVOL, XVOL(10), XH2O(8), XH2O(10), &
                  AClRadi, C(ind_NIT), C(ind_SALACL), TEMP, RELHUM )
      ELSEIF (X .eq. 2) THEN
         output = N2O5_InorgOrg( XVOL(12), 0e+0_fp, XH2O(12), 0e+0_fp, &
                  XRADI(12), C(ind_NITs), C(ind_SALCCL), TEMP, RELHUM )
      ENDIF

      ! Gamma
      XSTKCF = output(1)
      ! ClNO2 yield, fraction [0,1]
      IF (X .EQ. 1) THEN
         ! Reduce ClNO2 production yield on fine inorganic+organic aerosol by 75%,
         ! this is based on the analysis of WINTER observations in McDuffie
         ! et al., JGR, 2018
         ClNO2_yield = output(2)*0.25_fp
      ELSEIF (X .eq. 2) THEN
         ClNO2_yield = output(2)
      ENDIF
      ! Particle radius with coating, cm
      Rp     = output(3)
      ! Surface area of coated particles, cm2/cm3
      SA     = output(4)

      ! Total loss rate of N2O5 (kN2O5) on SNA+ORG+SSA aerosol
      kISum(1) = ARSL1K( (1-CldFr)*SA, Rp, XDENA, XSTKCF, XTEMP, XSQM )
      kISum(1) = ClNO2_yield * kISum(1)

      kISum(2) = XSTKCF
      kISum(3) = ClNO2_yield

    END FUNCTION HETN2O5_SS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetN2O5_HCl
!
! !DESCRIPTION: Set heterogenous chemistry rate for N2O5(g) + HCl(l,s)
!  in polar stratospheric clouds and on tropospheric sulfate aerosol.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETN2O5_HCl( A, B, Input_Opt ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: A, B       ! Rate coefficients
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms.
!
! !REVISION HISTORY:
!  29 Jan 2016 - M. Sulprizio- Initial version, adapted from code previously
!                              in calcrate.F
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE, SA_total

      ! Initialize
      kISum        = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp
      SA_total     = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Assume zero
         XStkCf = 0.0e+0_fp

	 ! restore stratosphere only limitation - TMS 17/04/10
         IF ( STRATBOX ) THEN
            IF (N.eq.8) THEN
               ! Fixed gamma?
               !XSTKCF = 0.1e-4_fp ! Sulfate
               ! RH dependence
               ! Note gamma calculation on sulfate in troposphere uses
               ! the McDuffie parameterization (func: HetN2O5())
      	       XSTKCF = N2O5( N, TEMP, RELHUM )
	    ENDIF
         ENDIF

         ! For UCX-based mechanisms only consider PSC reactions in strat
         IF ( Input_Opt%LUCX .and. STRATBOX ) THEN
            IF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(2)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.003e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.03e+0_fp ! Ice
               ENDIF
            ENDIF
         ENDIF

         IF (XStkCf.gt.0.0e+0_fp) THEN
            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            ! Add to overall reaction rate
            kISum = kISum + ADJUSTEDRATE
         ENDIF

      END DO

    END FUNCTION HETN2O5_HCl
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
      REAL(fp), INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp), INTENT(IN) :: rAer        ! Radius of aerosol (cm)
      REAL(fp), INTENT(IN) :: AAer        ! Area of aerosol (cm2/cm3)
      REAL(fp), INTENT(IN) :: TK          ! Temperature (K)
      INTEGER,  INTENT(IN) :: X           ! 1: Cl-, 2: Br-
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
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
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE, XSqM
      Real(fp), Parameter :: XMolWeightHCl=36.5e+0_fp
      Real(fp), Parameter :: XSqMHCl=SQRT(XMolWeightHCl)
      Real(fp), Parameter :: XMolWeightHBr=81.0e+0_fp
      Real(fp), Parameter :: XSqMHBr=SQRT(XMolWeightHBr)

      ! Initialize
      kISum        = 0.0_fp

      ! Select between halogens
      IF (X.eq.1) THEN
         ! This would never be used since HCl uptake is
         ! handled by ISSOROPIA now, xnw 1/25/18
         XSqM = XSqMHCl
      ELSEIF (X.eq.2) THEN
         XSqM = XSqMHBr
      ENDIF

      XStkCf = Gamma_HX_Uptake( rAer, denAir, X, TK )

      ! Reaction rate for surface of aerosol
      kISum = Arsl1K(AAer,rAer,denAir,XStkCf,XTemp,XSqM)

    END FUNCTION HETHXUptake
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetO3_SS
!
! !DESCRIPTION: Sets the O3 + Br- (in sea salt) rate using Johan
!  Schmidt's updated code.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETO3_SS( denAir, rAer, AAer, alkAer, TK, halConc, O3Conc, H  ) &
                     RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: denAir   ! Density of air (#/cm3)
      REAL(fp),       INTENT(IN) :: rAer     ! Radius of aerosol (cm)
      REAL(fp),       INTENT(IN) :: AAer     ! Area of aerosol (cm2/cm3)
      REAL(fp),       INTENT(IN) :: alkAer   ! Aerosol alkalinity (?)
      REAL(fp),       INTENT(IN) :: TK       ! Temperature (K)
      REAL(fp),       INTENT(IN) :: halConc  ! Halide concentration (mol/L)
      REAL(fp),       INTENT(IN) :: O3Conc   ! Ozone concentration (#/cm3)
      TYPE(HetState), POINTER    :: H        ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum    ! Rxn rate O3 + Br- in sea salt
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: XSTKCF, XSQM

      ! Initialize
      kISum = 0.0_fp
      XSQM  = SQRT( H%O3%MW_g )

      ! Reaction can only proceed on acidic aerosol
      IF ( alkAer > 0.05_fp ) THEN
         XStkCf = 0.0_fp
      ELSE
         XStkCf = Gamma_O3_Br( rAer, denAir, TK, halConc, O3Conc, H )
      ENDIF

      ! Reaction rate for surface of aerosol
      kISum = Arsl1K( AAer, rAer, denAir, XStkCf, XTemp, XSqM )

    END FUNCTION HETO3_SS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetNO3_Cl
!
! !DESCRIPTION: Sets the NO3(g) hypsis rate on Cl-.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETNO3_Cl( denAir, rAer, AAer, TK, clConc, X, H ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: denAir   ! Density of air (#/cm3)
      REAL(fp),       INTENT(IN) :: rAer     ! Radius of aerosol (cm)
      REAL(fp),       INTENT(IN) :: AAer     ! Area of aerosol (cm2/cm3)
      REAL(fp),       INTENT(IN) :: TK       ! Temperature (K)
      REAL(fp),       INTENT(IN) :: clConc   ! Cloride concentration (mol/L)
      INTEGER,        INTENT(IN) :: X
      TYPE(HetState), POINTER    :: H        ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum       ! NO3(g) reaction rate on Cl-
!
! !REVISION HISTORY:
!  14 Mar 2018 - X. Wang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: XStkCf, XSqM

      ! Initialize
      kISum = 0.0_fp
      XSQM   = SQRT( H%NO3%MW_g )

      ! Compute reactive uptake coefficient
      XStkCf = GAMMA_NO3( rAer, TK, clConc, X, H ) * 0.01_fp

      ! Reaction rate for surface of aerosol
      kISum = Arsl1K( AAer, rAer, denAir, XStkCf, XTemp, XSqM )

    END FUNCTION HETNO3_Cl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetClNO2_TCld
!
! !DESCRIPTION: Sets the rate of the multiphase reaction ClNO2 + Cl-/Br- in
!  troposphere cloud
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETClNO2_TCld( denAir,   rLiq,     rIce,     ALiq,              &
                            AIce,     VAir,     TK,       CldFr,             &
                            pH,       clConc_A, clConc_C, clConc_g,          &
                            brConc_A, brConc_C, brConc_g, X,                 &
                            H                                              ) &
                          RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp),       INTENT(IN) :: rLiq        ! Radius of liquid cloud
                                                !  droplets (cm)
      REAL(fp),       INTENT(IN) :: rIce        ! Radius of ice cloud
                                                !  crystals (cm)
      REAL(fp),       INTENT(IN) :: ALiq        ! Area of liquid cloud
                                                !  droplets (cm2/cm3)
      REAL(fp),       INTENT(IN) :: AIce        ! Area of ice cloud
                                                !  crystals (cm2/cm3)
      REAL(fp),       INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(fp),       INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp),       INTENT(IN) :: CldFr       ! Cloud fraction
      REAL(fp),       INTENT(IN) :: clConc_A    ! Fine Chloride
                                                !  concentration (mol/L)
      REAL(fp),       INTENT(IN) :: clConc_C    ! Coarse Chloride
                                                !  concentration (mol/L)
      REAL(fp),       INTENT(IN) :: clConc_g
      REAL(fp),       INTENT(IN) :: brConc_A    ! Fine Bromide
                                                !  concentration (mol/L)
      REAL(fp),       INTENT(IN) :: brConc_C    ! Coarse Bromide
                                                !  concentration (mol/L)
      REAL(fp),       INTENT(IN) :: brConc_g
      REAL(fp),       INTENT(IN) :: pH
      INTEGER,        INTENT(IN) :: X           ! 1=fineCl;2=coarseCl;
                                                ! 3=HCl,4=fineBr;5=coarseBr;
                                                ! 6=HBr
      TYPE(HetState), POINTER    :: H           ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum       ! Rxn rate
!
 !REMARKS:
!
! !REVISION HISTORY:
!  08 Aug 2018 - X. Wang      - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: X1, X2, ADJUSTEDRATE
      REAL(fp) :: GAM_ClNO2, r_gp, clConc, r_ac, brConc, r1
      INTEGER  :: Y

      ! Initialize the return value
      kISum = 0.0_fp

      ! Return if we are in the stratosphere
      ! This skips unneccesary computations & function calls
      IF ( StratBox ) RETURN

      ! Continue initializing
      X1    = 0.0_fp

      ! Cloud halide concentration, put fine and coarse mode together
      clConc = clConc_A + clConc_C + clConc_g
      brConc = brConc_A + brConc_C + brConc_g

      If (X == 1) THEN
         r_ac = clConc_A / clConc
         Y = 1
      ELSEIF (X == 2) THEN
         r_ac = clConc_C / clConc
         Y = 1
      ElSEIF (X == 3) THEN
         r_ac = clConc_g / clConc
         Y = 1
      ELSEIF (X == 4) THEN
         r_ac = brConc_A / brConc
         Y = 2
      ELSEIF (X == 5) THEN
         r_ac = brConc_C / brConc
         Y = 2
      ELSEIF (X == 6) THEN
         r_ac = brConc_g / brConc
         Y = 2
      ENDIF

      CALL GAMMA_ClNO2( rLiq,   denAir, TK,        pH,   clConc,              &
                        brConc, Y,      GAM_ClNO2, r_gp, H                   )
      X1 = GAM_ClNO2
      r1 = r_gp * r_ac

      X2 = 0

      kISum = CloudHet( CldFr, Aliq,   Aice,         rLiq, rIce,             &
                        TK,    denAir, H%ClNO2%MW_g, X1,   X2,  r1               )

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
      REAL(fp),       INTENT(IN) :: denAir   ! Density of air (#/cm3)
      REAL(fp),       INTENT(IN) :: rAer     ! Radius of aerosol (cm)
      REAL(fp),       INTENT(IN) :: AAer     ! Area of aerosol (cm2/cm3)
      REAL(fp),       INTENT(IN) :: alkAer   ! Aerosol alkalinity (?)
      REAL(fp),       INTENT(IN) :: TK       ! Temperature (K)
      REAL(fp),       INTENT(IN) :: clConc   ! Cloride concentration (mol/L)
      REAL(fp),       INTENT(IN) :: brConc   ! Bromide concentration (mol/L)
      REAL(fp),       INTENT(IN) :: pH       ! Aerosol pH
      INTEGER,        INTENT(IN) :: X        ! 1: Cl-, 2: Br-
      TYPE(HetState), POINTER    :: H        ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum
!
! !REVISION HISTORY:
!  14 Mar 2018 - X. Wang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: XSQM, GAM_ClNO2, r_gp

      ! Initialize
      kISum     = 0.0_fp
      GAM_ClNO2 = 0.0_fp
      r_gp      = 0.0_fp
      XSQM      = SQRT( H%ClNO2%MW_g )

      CALL GAMMA_ClNO2( rAer,   denAir, TK,        pH,   clConc,             &
                        brConc, X,      GAM_ClNO2, r_gp, H                  )

      ! Reaction rate for surface of aerosol
      kISum = Arsl1K( AAer, rAer, denAir, GAM_ClNO2, XTemp, XSqM ) * r_gp

    END FUNCTION HETClNO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gamma_NO3
!
! !DESCRIPTION: Function GAMMA\_NO3 calculates reactive uptake coef.
!               for NO3 on salts and water
!\\
!\\
! !INTERFACE:
!
    FUNCTION GAMMA_NO3( Radius, T, C_X, X, H ) RESULT( GAM_NO3 )
!
! !USES:
!
      USE PhysConstants, ONLY : Pi, RStarG
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: Radius
      REAL(fp),       INTENT(IN) :: T        ! Temperature (K)
      REAL(fp),       INTENT(IN) :: C_X      ! Cl- Concentration (mol/L)
      INTEGER,        INTENT(IN) :: X        ! X = 1 fine; 2 coarse
      TYPE(HetState), POINTER    :: H        ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(fp)                   :: GAM_NO3  ! Reactive uptake coefficient [1]

!
! !REVISION HISTORY:
!  24 Sep 2015 - J. Schmidt  - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!
      ! Scalars
      REAL(fp) :: ab, M_X, k_tot, H_X, k1, k2
      REAL(fp) :: cavg, D_l, gb, l_r
      REAL(fp) :: WaterC, Vol
      REAL(fp) :: H_K0_O3

      ! Henry's law [M/atm]
      H_K0_O3 = H%O3%K0 * con_atm_bar
      H_X     = H_K0_O3 * EXP( H%O3%CR * ( 1.0e0_fp/T - 1.0e0_fp/H%O3%TK ) )

      ! O3 mol wt (kg/mol)
      M_X = H%O3%MW_g * 1e-3_fp

      WaterC = AWATER(X) / 18e12_fp !mol/cm3 air
      IF (X == 1) THEN
         Vol = AClAREA * AClRADI * 1.0e-3_fp / 3.0e0_fp !L/cm3 air
      ELSE
         Vol = XAREA(12) * XRADI(12) * 1.0e-3_fp / 3.0e0_fp !L/cm3 air
      ENDIF

      WaterC = WaterC / Vol !mol/L aerosol

      ! HNO3 mol wt (kg/mol)
      M_X = H%NO3%MW_g * 1.0e-3_fp

      ! Mass accommodation coefficient
      ab = 1.3e-2_fp

      cavg = SQRT(8.0e+0_fp*RStarG*T/(Pi*M_X)) *1.0e2_fp ! thermal velocity (cm/s)

      ! Liquid phase diffusion coefficient [cm2/s] for NO3
      ! (Ammann et al., Atmos. Chem. Phys., 2013)
      D_l  = 1.0e-5_fp

      k1 = 2.76e+6_fp !M-1s-1
      k2 = 23e+0_fp !M-1s-1

      k_tot = k1*C_X + k2*WaterC

      H_X = 0.6e+0_fp !M atm-1
      H_X = H_X * con_atm_bar !M/bar
      l_r = SQRT(D_l / k_tot)
!
!      IF (K_Tot .EQ. 0.0) THEN
!         K_tot = 1.0e-2_fp
!      ENDIF

      gb = 4.0e0_fp * H_X * con_R * T * l_r * k_tot / cavg

      gb = gb * REACTODIFF_CORR( Radius, l_r)

      GAM_NO3 = 1.0e0_fp / (1.0e0_fp/ab  +  1.0e0_fp/gb)

    END FUNCTION GAMMA_NO3
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
      REAL(fp),       INTENT(IN)  :: Radius   ! Radius (cm)
      REAL(fp),       INTENT(IN)  :: n_air    ! n_air (#/cm)
      INTEGER,        INTENT(IN)  :: X        ! 1=Cl-,2=HSO3-/SO3--
      REAL(fp),       INTENT(IN)  :: T        ! Temperature (K)
      REAL(fp),       INTENT(IN)  :: C_Y1     ! Concentration (mol/L)
      REAL(fp),       INTENT(IN)  :: C_Y3     ! Concentration (mol/L)
      REAL(fp),       INTENT(IN)  :: C_Y4     ! Concentration (mol/L)
      REAL(fp),       INTENT(IN)  :: C_Hp     ! Concentration (mol/L)
      TYPE(HetState), POINTER     :: H        ! Hetchem species metadata
!
! !OUTPUT PARAMETERS:
!
      ! Reactive uptake coefficient (unitless)
      REAL(fp),       INTENT(OUT) :: GAM_HOCl
      REAL(fp),       INTENT(OUT) :: r_gp
!
! !REVISION HISTORY:
!  01 Nov 2018 - X. Wang   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      REAL(fp) :: ab,     M_X,  k1,  k2,  H_X,   C_Y2
      REAL(fp) :: gb_tot, cavg, D_l, l_r, k_tot, H_HOCl

      C_Y2     = C_Y3 + C_Y4

      ! MW of HOCl (kg/mol)
      M_X      = H%HOCl%MW_g * 1.0e-3_fp

      ! Mass accommodation coefficient
      ab       = 0.8_fp

      ! thermal velocity (cm/s)
      cavg     = SQRT( 8.0_fp * RStarG * T / ( Pi * M_X) ) *1.0e2_fp

      ! Liquid phase diffusion coefficient [cm2/s] for HOCl
      ! (Ammann et al., Atmos. Chem. Phys., 2013)
      D_l      = 2.0e-5_fp

      k1       = 1.5e+4_fp !M-1 s-1
      k2       = 2.8e+5_fp !M-1 s-1 Liu and Abbatt, Geophys. Res. Lett., 2020
      k_tot    = k1 * C_Hp * C_Y1                                            &
               + k2 * C_Y2

      ! Henry's law
      H_HOCl   = H%HOCl%K0 * con_atm_bar
      H_X      = H_HOCl * EXP( H%HOCl%CR * ( 1.0_fp/T - 1.0_fp/H%HOCl%TK ) )

      l_r      = SQRT( D_l / k_tot )
      gb_tot   = 4.0_fp * H_X * con_R * T * l_r * k_tot / cavg
      gb_tot   = gb_tot * REACTODIFF_CORR( Radius, l_r )

      ! Reactive uptake coefficient [unitless]
      GAM_HOCl = 1.0_fp / ( 1.0_fp/ab  +  1.0_fp/gb_tot )

      ! turn off HOCl+S(IV), Xuan Wang
      !gb2 = 0.0e0_fp

      IF ( X==1 ) THEN
         r_gp = k1 * C_Hp * C_Y1  /k_tot
      ELSEIF ( X==2 ) THEN
         r_gp = k2 * C_Y2 / k_tot
      ENDIF

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
      REAL(fp),       INTENT(IN) :: Radius   ! Radius (cm)
      REAL(fp),       INTENT(IN) :: n_air    ! n_air (#/cm)
      REAL(fp),       INTENT(IN) :: T        ! Temperature (K)
      REAL(fp),       INTENT(IN) :: C_H      ! H+ concentration
      REAL(fp),       INTENT(IN) :: C_X      ! Cl- Concentration (mol/L)
      TYPE(HetState), POINTER    :: H        ! Hetchem species metadata
!
! !RETURN VALUE:
!
      ! Reactive uptake coefficient (unitless)
      REAL(fp)                   :: GAM_HOCl
!
! !REVISION HISTORY:
!  16 Mar 2018 - X. Wang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      REAL(fp) :: ab, M_X, k_ter, H_X, H_HOCl, cavg, D_l, gb, l_r

      ! MW of HOCl (kg/mol)
      M_X      = H%HOCl%MW_g * 1.0e-2_fp

      ! Mass accommodation coefficient
      ab       = 0.8_fp

      ! thermal velocity (cm/s)
      cavg     = SQRT( 8.0_fp * RStarG * T / ( Pi * M_X ) ) *1.0e2_fp

      ! Liquid phase diffusion coefficient [cm2/s] for HOCl
      ! (Ammann et al., Atmos. Chem. Phys., 2013)
      D_l      = 2.0e-5_fp

      k_ter    = 1.5e+4_fp                 ! M-1s-1
      H_HOCl   = H%HOCl%K0 * con_atm_bar   ! M/bar
      H_X      = H_HOCl*EXP( H%HOCl%CR *( 1.0_fp/T - 1.0_fp/H%HOCl%TK ) )

      l_r = SQRT(D_l / (k_ter * C_H * C_X))
      gb       = 4.0_fp * H_X * con_R * T * l_r * k_ter * C_H * C_X / cavg
      gb       = gb * REACTODIFF_CORR( Radius, l_r)

      GAM_HOCl = 1.0_fp / (1.0_fp/ab  +  1.0_fp/gb)

      !GAM_HOCl = MIN(GAM_HOCl, 2.0e-4_fp)

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
      REAL(fp),       INTENT(IN)   :: Radius   ! Radius (cm)
      REAL(fp),       INTENT(IN)   :: n_air    ! n_air (#/cm)
      REAL(fp),       INTENT(IN)   :: T        ! Temperature (K)
      REAL(fp),       INTENT(IN)   :: C_X1     ! Cl-/Br- Concentration (mol/L)
      REAL(fp),       INTENT(IN)   :: C_X2     ! Cl-/Br- Concentration (mol/L)
      INTEGER,        INTENT(IN)   :: X        ! 1=Cl-,2=Br-
      REAL(fp),       INTENT(IN)   :: pH
      TYPE(HetState), POINTER      :: H        ! Hetchem species metadata
!
! !OUTPUT PARAMETERS:
!
      ! Reactive uptake coefficient (unitless)
      REAL(fp),       INTENT(OUT)  :: GAM_ClNO2
      REAL(fp),       INTENT(OUT)  :: r_gp
!
! !REVISION HISTORY:
!  25 Jan 2018 - X. Wang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!
      ! Scalars
      REAL(fp) :: ab, M_X, H_X, fCl, k_tot
      REAL(fp) :: cavg, D_l, k_b1, k_b2, l_r, gb_tot

      ! MW of ClNO2 (kg/mol)
      M_X       = H%ClNO2%MW_g * 1.0e-3_fp

      ! Mass accommodation coefficient
      ab        = 0.01_fp

      ! thermal velocity (cm/s)
      cavg      = SQRT( 8.0_fp * RStarG * T / ( Pi * M_X ) ) *1.0e2_fp

      ! Liquid phase diffusion coefficient [cm2/s] for ClNO2
      ! (Ammann et al., Atmos. Chem. Phys., 2013)
      D_l       = 1.0e-5_fp
      H_X       = 4.5e-2_fp         !M atm-1
      H_X       = H_X * con_atm_bar !M/bar

      !ClNO2+Cl-
      k_b1      = 1.0e+7_fp !M-1s-1
      IF (pH >= 2) THEN
         k_b1   = 0.0e0_fp
      ENDIF

      !ClNO2+Br-
      k_b2      = 1.01e-1_fp / (H_X*H_X*D_l)

      k_tot     = k_b1*C_X1 + k_b2*C_X2
      l_r       = SQRT(D_l / k_tot)
      gb_tot    = 4.0_fp * H_X * con_R * T * l_r * k_tot / cavg
      gb_tot    = gb_tot * REACTODIFF_CORR( Radius, l_r )

      GAM_ClNO2 = 1.0_fp / ( 1.0_fp/ab  +  1.0_fp/gb_tot )

      IF ( X==1 ) THEN
         r_gp = k_b1*C_X1 / k_tot
      ELSEIF ( X==2 ) THEN
         r_gp = k_b2*C_X2 / k_tot
      ENDIF

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
    FUNCTION HETClNO3_SS( denAir, rAer,   AAer, alkAer, TK,                  &
                          clConc, brConc, X,    M,      H   ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: denAir   ! Density of air (#/cm3)
      REAL(fp),       INTENT(IN) :: rAer     ! Radius of aerosol (cm)
      REAL(fp),       INTENT(IN) :: AAer     ! Area of aerosol (cm2/cm3)
      REAL(fp),       INTENT(IN) :: alkAer   ! Aerosol alkalinity (?)
      REAL(fp),       INTENT(IN) :: TK       ! Temperature (K)
      REAL(fp),       INTENT(IN) :: clConc   ! Cloride concentration (mol/L)
      REAL(fp),       INTENT(IN) :: brConc   ! Bromide concentration (mol/L)
      Integer,        INTENT(IN) :: X        ! 1: Cl-, 2: Br-
      Integer,        INTENT(IN) :: M        ! 1: fine, 2: coarse
      TYPE(HetState), POINTER    :: H        ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum    ! Rxn rate ClNO3 + Br- in sea salt
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
      ! Scalars
      INTEGER  :: N
      REAL(fp) :: GAM_ClNO3, r_gp, XSQM

      ! Initialize
      kISum = 0.0_fp
      XSQM  = SQRT( H%ClNO3%MW_g )

      CALL Gamma_ClNO3_AER( rAer,  denAir,  X,         TK,   M,              &
                            clConc, brConc, GAM_ClNO3, r_gp, H              )

      ! Reaction rate for surface of aerosol
      kISum = Arsl1K( AAer, rAer, denAir, GAM_ClNO3, XTemp, XSqM ) * r_gp

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
      REAL(fp),       INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp),       INTENT(IN) :: rAer        ! Radius of aerosol (cm)
      REAL(fp),       INTENT(IN) :: AAer        ! Area of aerosol (cm2/cm3)
      REAL(fp),       INTENT(IN) :: alkAer      ! Aerosol alkalinity (?)
      REAL(fp),       INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp),       INTENT(IN) :: hConc       ! H+ concentration (mol/L)
      REAL(fp),       INTENT(IN) :: clConc      ! Cloride concentration (mol/L)
      REAL(fp),       INTENT(IN) :: brConc      ! Bromide concentration (mol/L)
      INTEGER,        INTENT(IN) :: X           ! 1: Cl-, 2: Br-
      TYPE(HetState), POINTER    :: H           ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: XSTKCF, GAM_HOBr, r_gp, XSQM

      ! Initialize
      kISum = 0.0_fp
      XSQM  = SQRT( H%HOBr%MW_g )

      ! Reaction can only proceed on acidic aerosol
      IF ( alkAer > 0.05_fp ) THEN
         XStkCf = 0.0_fp
         r_gp   = 0.0_fp
      ELSE
         CALL Gamma_HOBr_AER( rAer,   denAir,  X,       TK,   clConc,        &
                              brConc, hConc,  GAM_HOBr, r_gp, H             )
         XStkCf = GAM_HOBr
      ENDIF

      ! Reaction rate for surface of aerosol
      kISum = Arsl1K( AAer, rAer, denAir, XStkCf, XTemp, XSqM ) * r_gp

    END FUNCTION HETHOBr_SS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetO3_TCld
!
! !DESCRIPTION: Sets the O3 + Br- rate in troposphere cloud
!
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETO3_TCld( denAir,   rLiq,   rIce,  ALiq,     AIce,            &
                         VAir,     TK,     CldFr, brConc_a, brConc_c,        &
                         brConc_g, O3Conc, X,     H                        ) &
                         RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp),       INTENT(IN) :: rLiq        ! Radius of liquid cloud
                                                !  droplets (cm)
      REAL(fp),       INTENT(IN) :: rIce        ! Radius of ice cloud crystals
                                                !  (cm)
      REAL(fp),       INTENT(IN) :: ALiq        ! Area of liquid cloud droplets
                                                !  (cm2/cm3)
      REAL(fp),       INTENT(IN) :: AIce        ! Area of ice cloud crystals
                                                !  (cm2/cm3)
      REAL(fp),       INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(fp),       INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp),       INTENT(IN) :: CldFr       ! Cloud fraction
      REAL(fp),       INTENT(IN) :: brConc_a    ! Bromide concentration (mol/L)
      REAL(fp),       INTENT(IN) :: brConc_c    ! Bromide concentration (mol/L)
      REAL(fp),       INTENT(IN) :: brConc_g    ! Bromide concentration (mol/L)
      REAL(fp),       INTENT(IN) :: O3Conc      ! Ozone concentration (mol/L)
      INTEGER,        INTENT(IN) :: X           ! X=1 _a; 2 _c; 3 _g
      TYPE(HetState), POINTER    :: H           ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum
!
! !REVISION HISTORY:
!  08 Jun 2019 - X. Wang     - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER  :: N
      REAL(fp) :: brConc, r_ac, X1, X2

      ! Initialize the return value
      kISum  = 0.0_fp

      ! Return if we are in the stratosphere
      ! This skips unneeded computations & function calls
      IF ( StratBox ) RETURN

      ! Continue initializing
      r_ac   = 0.0_fp
      brConc = 0.0_fp
      X1     = 0.0_fp

      brConc = brConc_a + brConc_g + brConc_c

      IF (X.EQ.1) THEN
         r_ac = brConc_a / brConc
      ELSE IF (X .EQ. 2) THEN
         r_ac = brConc_c / brConc
      ELSE
         r_ac = brConc_g / brConc
      ENDIF

      ! Reaction on liquid clouds (tropospheric only)
      X1 = Gamma_O3_Br( rLiq, denAir, TK, brConc, O3Conc, H )
      X2 = 0.0_fp
      kISum = CloudHet( CldFr, Aliq,   Aice,      rLiq, rIce,               &
                        TK,    denAir, H%O3%MW_g, X1,   X2,   r_ac               )

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
      REAL(fp),       INTENT(IN) :: Radius   ! Radius (cm)
      REAL(fp),       INTENT(IN) :: n_air    ! Number density of air (#/cm)
      REAL(fp),       INTENT(IN) :: T        ! Temparature (K)
      REAL(fp),       INTENT(IN) :: C_Y      ! Concentration
      REAL(fp),       INTENT(IN) :: C_X_g    ! Gas-phase concentration
      TYPE(HetState), POINTER    :: H        ! Hetchem species metadata
!
! !RETURN VALUE
!
      REAL(fp)                   :: GAM      ! Reactive uptake coeff. (1)
!
! !REVISION HISTORY:
!  24 Sept 2015 - J. Schmidt - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      REAL(fp) :: ab,     gb,  gd,       gs
      REAL(fp) :: cavg,   H_X, H_O3,     M_X
      REAL(fp) :: KLangC, k_s, C_Y_surf, Nmax
      REAL(fp) :: k_b,    D_l, l_r

      ! Henry's law
      H_O3     = H%O3%K0 * con_atm_bar
      H_X      = H_O3*EXP( H%O3%CR * ( 1.0_fp/T - 1.0_fp/H%O3%TK ) )

      ! Molwt of O3 (kg/mol)
      M_X      = H%O3%MW_g * 1.0e-3_fp

      ! Thermal velocity (cm/s)
      cavg     = SQRT( 8 * RStarG * T / ( Pi * M_X ) ) *1.0e2_fp

      Nmax     = 3.0e14_fp  ! #/cm2
      KLangC   = 1.0e-13_fp !cm3
      k_s      = 1.0e-16_fp !cm2s-1, from ks*Nmax=0.03s-1

     ! [Br-(surf)] = 3.41E14 cm-2/M * [Br-(bulk)], but not gt Nmax.
      C_Y_surf = MIN( 3.41e14_fp * C_Y, Nmax )
      gs       = ( 4.0_fp * k_s * C_Y_surf * KLangC * Nmax )                 &
               / ( cavg * ( 1.0_fp + KLangC * C_X_g )      )

      k_b      = 6.3e8_fp *  EXP(-4.45e3_fp / T) !M-1 s-1
      D_l      = 8.9e-6_fp !cm2 s-1.
      l_r      = SQRT( D_l / (k_b * C_Y ) )! cm
      gb       = 4.0_fp * H_X * con_R * T * l_r * k_b * C_Y / cavg
      gb       = gb * REACTODIFF_CORR( Radius, l_r)

      GAM      = gb + gs

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
      REAL(fp),       INTENT(IN) :: A, B
      TYPE(OptInput), INTENT(IN) :: Input_Opt
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum
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
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp) :: XSQM

      ! Initialize
      kISum        = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp
      XSQM         =sqrt(A)

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Get the aerosol type
         XStkCf = 0.0e+0_fp

         IF ( N == 8 ) THEN

            ! sulfate aerosol
            ! This seems not to happen

         ELSEIF ( Input_Opt%LUCX .and. STRATBOX ) THEN
            IF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(5)
               ! Turn off het rates in stratosphere
               IF ( Input_Opt%TurnOffHetRates ) XSTKCF = 0.0e+0_fp       ! EF turn off ClNO3 + HBr on sulfates as in GMI
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.3e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.3e+0_fp ! Ice
               ENDIF
               ! Turn off het rates in stratosphere
               IF ( Input_Opt%TurnOffHetRates ) XSTKCF = 0.0e+0_fp       ! EF turn off ClNO3 + HBr on ALL PSCs as in GMI
            ENDIF

         ENDIF

         IF (XStkCf.gt.0.0e+0_fp) THEN
            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP,XSQM)
            ENDIF

            ! Add to overall reaction rate
            kISum = kISum + ADJUSTEDRATE
         ENDIF
      END DO

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
    FUNCTION HETClNO3( denAir, TK, clConc, brConc, CldFr, H ) &
             RESULT( HET_ClNO3 )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: denAir   ! Density of air (#/cm3)
      REAL(fp),       INTENT(IN) :: TK       ! Temperature (K)
      REAL(fp),       INTENT(IN) :: clconc   ! Cl- & Br- concentration in
      REAL(fp),       INTENT(IN) :: brconc   !  fine mode (M)
      REAL(fp),       INTENT(IN) :: CldFr    ! Cloud fraction
      TYPE(HetState), POINTER    :: H        ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(fp)                   :: HET_ClNO3  ! Hydrol. rate for ClNO3
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
      ! Scalars
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE, GAM_ClNO3, r_gp, XSQM

      ! Initialize
      HET_ClNO3    = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp
      GAM_ClNO3    = 0.0_fp
      r_gp         = 0.0_fp
      XSQM         = SQRT( H%ClNO3%MW_g )

      ! Only apply PSC rate adjustment if at high altitude
      DO_EDUCT = STRATBOX

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE
         XSTKCF = 0.0_fp

         ! Get the aerosol type
         IF ( N == 8 ) THEN
            ! Follow ClNO3 + Cl- channel
            XSTKCF = 0.0_fp
         ELSEIF (N == 11) THEN
            ! Follow ClNO3 + Cl- channel, xnw 1/25/18
            CALL GAMMA_ClNO3_AER( AClRADI, denAir, 3,         TK,   1,       &
                                  clconc,  brconc, GAM_ClNO3, r_gp, H       )
            XSTKCF = GAM_ClNO3
         ELSEIF (N == 12) THEN
            ! Follow ClNO3 + Cl- channel
            XSTKCF = 0.0_fp
         ELSEIF (N.eq.13) THEN
            XSTKCF = KHETI_SLA(3)
         ELSEIF (N.eq.14) THEN
            IF (NATSURFACE) THEN
               XSTKCF = 0.004_fp ! NAT
            ELSE
               XSTKCF = 0.3_fp ! Ice
            ENDIF
         ELSE
            XSTKCF = 0.0_fp
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE IF (N .eq. 11) THEN
            ADJUSTEDRATE =ARSL1K((1-CldFr)*AClAREA,AClRADI,XDENA,XSTKCF,XTEMP,XSQM)*r_gp
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP,XSQM)
         ENDIF

         ! Add to overall reaction rate
         HET_ClNO3 = HET_ClNO3 + ADJUSTEDRATE
      END DO

    END FUNCTION HETClNO3
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
      REAL(fp),       INTENT(IN) :: A, B
      TYPE(OptInput), INTENT(IN) :: Input_Opt
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  21 Dec 2016 - S. D. Eastham - Generated code based on Johan Schmidt's work
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp) :: XSQM
      REAL(fp) :: GAM_HOBr, r_gp

      ! Initialize
      kISum        = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp
      XSQM         = SQRT(A)

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         XSTKCF = 0e+0_fp
         ! Get the aerosol type

         IF ( Input_Opt%LUCX .and. STRATBOX) THEN
            ! add limitation to stratosphere, xnw 1/25/18
            IF ( N == 8 ) THEN
               XSTKCF = 0.2e+0_fp ! Sulfate, [Hanson and Ravishankara, 1995]
            ELSEIF (N.eq.13) THEN
               ! SSA/STS
               XSTKCF = KHETI_SLA(10)
            ELSEIF (N.eq.14) THEN
               ! Ice/NAT PSC
               IF (NATSURFACE) THEN
                  XSTKCF = 0.1e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.3e+0_fp ! Ice
               ENDIF
            ELSE
               XSTKCF = 0e+0_fp
            ENDIF
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSEIF (XStkCf.gt.0.0e+0_fp) THEN
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP,XSQM)
         ENDIF

         ! Add to overall reaction rate
         kISum = kISum + ADJUSTEDRATE
      END DO

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
      REAL(fp),       INTENT(IN) :: A, B        ! Rate coefficients
      TYPE(OptInput), INTENT(IN) :: Input_Opt
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  21 Dec 2016 - S. D. Eastham - Generated code based on Johan Schmidt's work
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp) :: XSQM
      Real(fp) :: SADen

      ! Initialize
      kISum        = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp
      XSQM=SQRT(A)

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE
         XSTKCF = 0e+0_fp

         ! Get the aerosol type
         IF ( Input_Opt%LUCX .and. STRATBOX ) THEN
            ! add limitation to stratosphere, xnw 1/25/18
            IF ( N == 8 ) THEN
               XSTKCF = 0.25e+0_fp ! Sulfate, [Abbatt, 1995]
            ELSEIF ( N == 13 ) THEN
               ! SSA/STS
               XSTKCF = KHETI_SLA(6)
            ELSEIF ( N == 14 ) THEN
               ! Ice/NAT PSC
               IF (NATSURFACE) THEN
                  XSTKCF = 0.001e+0_fp
               ELSE
                  XSTKCF = 0.3e+0_fp
               ENDIF
            ELSE
               XSTKCF = 0e+0_fp
            ENDIF
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSEIF (XStkCf.gt.0.0e+0_fp) THEN
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTemp,XSQM)
         ENDIF

         ! Add to overall reaction rate
         kISum = kISum + ADJUSTEDRATE
      END DO

    END FUNCTION HETHOBr_HBr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetClNO3_TCld
!
! !DESCRIPTION: Sets the rate of the multiphase reaction ClNO3 + Cl-/Br- in
!  troposphere cloud
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETClNO3_TCld( denAir,   rLiq,     rIce,     ALiq,              &
                            AIce,     VAir,     TK,       CldFr,             &
                            clConc_A, clConc_C, clConc_g, brConc_A,          &
                            brConc_C, brConc_g, hno3_th,  hcl_th,            &
                            hbr_th,   X,        H                          ) &
                          RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp),       INTENT(IN) :: rLiq        ! Radius of liquid cloud
                                                !  droplets (cm)
      REAL(fp),       INTENT(IN) :: rIce        ! Radius of ice cloud
                                                !  crystals (cm)
      REAL(fp),       INTENT(IN) :: ALiq        ! Area of liquid cloud
                                                !  droplets (cm2/cm3)
      REAL(fp),       INTENT(IN) :: AIce        ! Area of ice cloud
                                                !  crystals (cm2/cm3)
      REAL(fp),       INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(fp),       INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp),       INTENT(IN) :: CldFr       ! Cloud fraction
      REAL(fp),       INTENT(IN) :: clConc_A    ! Fine Chloride
                                                !  concentration (mol/L)
      REAL(fp),       INTENT(IN) :: clConc_C    ! Coarse Chloride
                                                !  concentration (mol/L)
      REAL(fp),       INTENT(IN) :: clConc_g
      REAL(fp),       INTENT(IN) :: brConc_A    ! Fine Bromide
                                                !  concentration (mol/L)
      REAL(fp),       INTENT(IN) :: brConc_C    ! Coarse Bromide
                                                !  concentration (mol/L)
      REAL(fp),       INTENT(IN) :: brConc_g
      REAL(fp),       INTENT(IN) :: hno3_th
      REAL(fp),       INTENT(IN) :: hcl_th
      REAL(fp),       INTENT(IN) :: hbr_th
      INTEGER,        INTENT(IN) :: X           ! 1=fineCl;2=coarseCl;3=fineBr;
                                                ! 4=coarseBr;5=HCl;6=HBr;7=H2O
      TYPE(HetState), POINTER    :: H           ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
 !REMARKS:
!
! !REVISION HISTORY:
!  01 Feb 2018 - X. Wang      - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      INTEGER  :: N, Y
      REAL(fp) :: X1, X2, GAM_ClNO3, r1, r2, clConc, r_ac, brConc, r_gp

      ! Initialize the return value
      kISum = 0.0_fp

      ! Return if we are in the stratosphere
      ! This will prevent unnecessary computations & function calls
      IF ( StratBox ) RETURN

      ! Continue initializing
      X1    = 0.0_fp
      X2    = 0.0_fp

      ! Cloud halide concentration, put fine and coarse mode together
      clConc = clConc_A + clConc_C + clConc_g
      brConc = brConc_A + brConc_C + brConc_g

      If (X == 1) THEN
         Y = 1
         r_ac = clConc_A / clConc
      ELSEIF (X == 2) THEN
         Y = 1
         r_ac = clConc_C / clConc
      ElSEIF (X == 3) THEN
         Y = 2
         r_ac = brConc_A / brConc
      ELSEIF (X == 4) THEN
         Y = 2
         r_ac = brConc_C / brConc
      ELSEIF (X == 5) THEN
         Y = 1
         r_ac = clConc_g / clConc
      ELSEIF (X == 6) THEN
         Y = 2
         r_ac = brConc_g / brConc
      ELSEIF (X == 7) THEN
         Y = 3
         r_ac = 1.0_fp
      ENDIF

      CALL GAMMA_ClNO3_AER( rLiq,   denAir, Y,         TK,   3,              &
                            clConc, brConc, GAM_ClNO3, r_gp, H              )
      X1 = GAM_ClNO3
      r1 = r_gp*r_ac

      CALL GAMMA_ClNO3_ICE( Y, TK, hno3_th, hcl_th,hbr_th, GAM_ClNO3, r_gp)
      X2 = GAM_ClNO3
      IF (X >= 5) THEN
         r2 = r_gp
      ELSE
         r2 = 0.0_fp
      ENDIF

      kISum = CloudHet( CldFr,  Aliq,   Aice,         rLiq, rIce,             &
                        TK,     denAir, H%ClNO3%MW_g, X1,   X2,   r1,   r2         )

    END FUNCTION HetClNO3_TCld
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOBr_TCld
!
! !DESCRIPTION: Sets the rate of the multiphase reaction HOBr + Cl-/Br-/S(IV) in
!  troposphere cloud
!\\
!\\
! !INTERFACE:
!
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
      REAL(fp),       INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp),       INTENT(IN) :: rLiq        ! Radius of liquid cloud
                                                !  droplets (cm)
      REAL(fp),       INTENT(IN) :: rIce        ! Radius of ice cloud crystals
                                                !  (cm)
      REAL(fp),       INTENT(IN) :: ALiq        ! Area of liquid cloud droplets
                                                !  (cm2/cm3)
      REAL(fp),       INTENT(IN) :: AIce        ! Area of ice cloud crystals
                                                !  (cm2/cm3)
      REAL(fp),       INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(fp),       INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp),       INTENT(IN) :: CldFr       ! Cloud fraction
      REAL(fp),       INTENT(IN) :: hConc_Sul   ! Sulfate H+ concentration
      REAL(fp),       INTENT(IN) :: hConc_LCl   ! Liquid cloud H+ concentration
      REAL(fp),       INTENT(IN) :: hConc_ICl   ! Ice cloud H+ concentration
      REAL(fp),       INTENT(IN) :: clConc_A    ! Fine Chloride concentration
                                                !  (mol/L)
      REAL(fp),       INTENT(IN) :: clConc_C    ! Coarse Chloride concentration
                                                !  (mol/L)
      REAL(fp),       INTENT(IN) :: clConc_g
      REAL(fp),       INTENT(IN) :: brConc_A    ! Fine Bromide concentration
                                                !  (mol/L)
      REAL(fp),       INTENT(IN) :: brConc_C    ! Coarse Bromide concentration
                                                ! (mol/L)
      REAL(fp),       INTENT(IN) :: brConc_g
      REAL(fp),       INTENT(IN) :: hso3Conc    ! HSO3-    concentration (mol/L)
      REAL(fp),       INTENT(IN) :: so3Conc     ! SO3--    concentration (mol/L)
      REAL(fp),       INTENT(IN) :: hno3_th
      REAL(fp),       INTENT(IN) :: hcl_th
      REAL(fp),       INTENT(IN) :: hbr_th
      INTEGER,        INTENT(IN) :: X           ! 1=fineCl;2=coarseCl;
                                                ! 3=fineBr;4=coarseBr;5=HCl;
                                                ! 6=HBr;7=HSO3;8=SO3
      TYPE(HetState), POINTER    :: H           ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum
!
! !REVISION HISTORY:
!  01 Feb 2018 - X. Wang      - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N,        Y
      REAL(fp) :: X1,       X2,   ADJUSTEDRATE, XSQM
      REAL(fp) :: GAM_HOBr, r_gp, clConc,       r_ac
      REAL(fp) :: brConc,   r1,   r2

      ! Initialize return value
      kISum        = 0.0_fp

      ! Return if we are in the stratosphere
      ! This skips unneccessary computations & function calls
      IF ( StratBox ) RETURN

      ! Continue initializing
      ADJUSTEDRATE = 0.0_fp
      X1           = 0.0_fp
      X2           = 0.0_fp
      XSqM         = SQRT( H%HOBr%MW_g )

      ! Cloud halide concentration, put fine and coarse mode together
      clConc = clConc_A + clConc_C + clConc_g
      brConc = brConc_A + brConc_C + brConc_g

      If (X == 1) THEN
         Y = 1
         r_ac = clConc_A / clConc
      ELSEIF (X == 2) THEN
         Y = 1
         r_ac = clConc_C / clConc
      ElSEIF (X == 3) THEN
         Y = 2
         r_ac = brConc_A / brConc
      ELSEIF (X == 4) THEN
         Y = 2
         r_ac = brConc_C / brConc
      ELSEIF (X == 5) THEN
         Y = 1
         r_ac = clConc_g / clConc
      ELSEIF (X == 6) THEN
         Y = 2
         r_ac = brConc_g / brConc
      ELSEIF (X == 7) THEN
         Y = 3
         r_ac = 1.0_fp
      ELSEIF (X == 8) THEN
         Y = 4
         r_ac = 1.0_fp
      ENDIF

      CALL Gamma_HOBr_CLD( rLiq,      denAir,   Y,        TK,                &
                           clConc,    brConc,   hso3Conc, so3Conc,           &
                           hConc_LCl, GAM_HOBr, r_gp,     H                 )
      X1 = GAM_HOBr
      r1 = r_gp * r_ac

      CALL Gamma_HOBr_ICE( Y,       TK,      hno3_th, hcl_th,                &
                           hbr_th, GAM_HOBr, r_gp                           )
      X2 = GAM_HOBr
      IF ( (X >= 5) .AND. (X<=6) ) THEN
         r2 = r_gp
      ELSE
         r2 = 0.0_fp
      ENDIF

      kISum = CloudHet( CldFr,  Aliq,   Aice,        rLiq, rIce,              &
                        TK,     denAir, H%HOBr%MW_g, X1,   X2,   r1,   r2          )

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
      REAL(fp)                         :: GAM
! !INPUT PARAMETERS:
!
      ! Radius (cm), n_air (#/cm), and X (1 for Cl and 2 for Br)
      REAL(fp), INTENT(IN)             :: Radius, n_air
      INTEGER, INTENT(IN)              :: X
      REAL(fp), INTENT(IN)             :: T
!
! !REVISION HISTORY:
!  24 Sept 2015 - J. Schmidt - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!
      REAL(fp)       :: ab

      ! 1: Cl-, 2: Br-
      IF (X==1) THEN
         ! This would never be used since HCl uptake is
         ! handled by ISORROPIA now, xnw 1/25/18
         ab = 4.4e-6_fp * EXP( 2898.0e0_fp / T ) ! ab(RT) = 0.069
      ELSE
         ab = 1.3e-8_fp * EXP( 4290.0e0_fp / T ) ! ab(RT) = 0.021
      ENDIF

      GAM = ab

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
      REAL(fp),       INTENT(IN)  :: Radius   ! Radius (cm)
      REAL(fp),       INTENT(IN)  :: n_air    ! n_air (#/cm)
      INTEGER,        INTENT(IN)  :: X        ! 1=Cl-,2=Br-,3=HSO3-,4=SO3--
      REAL(fp),       INTENT(IN)  :: T        ! Temperature (K)
      REAL(fp),       INTENT(IN)  :: C_Y1     ! Concentration (mol/L)
      REAL(fp),       INTENT(IN)  :: C_Y2     ! Concentration (mol/L)
      REAL(fp),       INTENT(IN)  :: C_Y3     ! Concentration (mol/L)
      REAL(fp),       INTENT(IN)  :: C_Y4     ! Concentration (mol/L)
      REAL(fp),       INTENT(IN)  :: C_Hp     ! Concentration (mol/L)
      TYPE(HetState), POINTER     :: H        ! Hetchem species metadata
!
! !OUTPUT PARAMETERS:
!
      ! Reactive uptake coefficient (unitless)
      REAL(fp),       INTENT(OUT) :: GAM_HOBr
      REAL(fp),       INTENT(OUT) :: r_gp
!
! !REVISION HISTORY:
!  30 Nov 2017 - Q.J. Chen   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      REAL(fp) :: ab,     gd,   M_X,  cavg,  H_X
      REAL(fp) :: gb_tot, k_b1, k_b2, k_b3,  k_b4
      REAL(fp) :: D_l,    ybr2, l_r,  k_tot, C_Hp1
      REAL(fP) :: H_HOBr, C_Hp2

      ! Reaction rate coefficient for HOBr + Cl- [M-2 s-1]
      ! (Liu and Margerum, Environ. Sci. Tech., 2001)
      k_b1   = 2.3e+10_fp ! (qjc, 12/28/16)

      ! Reaction rate coefficient for HOBr + Br- [M-2 s-1]
      k_b2   = 1.6e+10_fp

      ! Reaction rate coefficient for HOBr + HSO3- [M-2 s-1]
      ! (Liu and Abbatt, Geophys. Res. Let., 2020)
      k_b3   = 2.6e+7_fp

      ! Reaction rate coefficient for HOBr + HSO3-- [M-2 s-1]
      ! (Troy and Margerum, Inorg. Chem., 1991)
      k_b4   = 5.0e+9_fp

      ! Liquid phase diffusion coefficient [cm2/s] for HOBr
      ! (Ammann et al., Atmos. Chem. Phys., 2013)
      D_l    = 1.4e-5_fp

      ! Henry's law
      H_HOBr = H%HOBr%K0 * con_atm_bar
      H_X    = H_HOBr * EXP( H%HOBr%CR *( 1.0e0_fp/T - 1.0e0_fp/H%HOBr%TK ) )
      M_X    = H%HOBr%MW_g * 1e-3_fp

      ! Mass accommodation coefficient
      ab     = 0.6e0_fp

      ! Thermal velocity [cm/s]
      cavg   = SQRT(8*RStarG*T/(pi*M_X)) *1.0e2_fp

      ! Follow Roberts et al, (2014)
      C_Hp1  = min(C_Hp, 1.0e-6)
      C_Hp2  = min(C_Hp, 1.0e-2)
      C_Hp1  = max(C_Hp1, 1.0e-9)
      C_Hp2  = max(C_Hp2, 1.0e-6)

      k_tot  = k_b1 * C_Y1 * C_Hp1                                            &
             + k_b2 * C_Y2 * C_Hp2                                            &
             + k_b3 * C_Y3                                                    &
             + k_b4 * C_Y4

      ! l_r is diffusive length scale [cm];
      !gb is Bulk reaction coefficient [unitless]
      l_r    = SQRT( D_l / k_tot )
      gb_tot = 4.0e0_fp * H_X * con_R * T * l_r * k_tot / cavg
      gb_tot = gb_tot * REACTODIFF_CORR( Radius, l_r)

      ! Reactive uptake coefficient [unitless]
      GAM_HOBr = 1.0e0_fp / (1.0e0_fp/ab  +  1.0e0_fp/gb_tot)

      ybr2     = 0.41e0*LOG10(C_Y2/C_Y1)+2.25        ! yield of Br2
      ybr2     = MIN(ybr2, 0.9e0)
      ybr2      = MAX(ybr2, TINY(1.0e0))

      IF ( X==1 ) THEN

         r_gp = (k_b1 * C_Y1 * C_Hp1 + k_b2 * C_Y2 * C_Hp2) / k_tot

         IF (C_Y2/C_Y1>5.e-4) THEN
            r_gp = 0.1e0 * r_gp
         ELSE
            r_gp = r_gp * (1.e0 - ybr2)
         ENDIF

      ELSEIF ( X==2 ) THEN

         r_gp = (k_b1 * C_Y1 * C_Hp1 + k_b2 * C_Y2 * C_Hp2) / k_tot

         IF (C_Y2/C_Y1>5.e-4) THEN
            r_gp = 0.9e0 * r_gp
         ELSE
            r_gp = r_gp * ybr2
         ENDIF

      ELSEIF ( X==3 ) THEN

         r_gp = (k_b3 * C_Y3) / k_tot

      ELSEIF ( X==4 ) THEN

         r_gp = (k_b4 * C_Y4) / k_tot

      ENDIF

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
      REAL(fp), INTENT(IN)  :: T        ! Temperature (K)
      REAL(fp), INTENT(IN)  :: hno3_th
      REAL(fp), INTENT(IN)  :: hcl_th
      REAL(fp), INTENT(IN)  :: hbr_th
!
! !OUTPUT PARAMETERS:
!
      ! Reactive uptake coefficient (unitless)
      REAL(fp), INTENT(OUT) :: GAM_HOBr
      REAL(fp), INTENT(OUT) :: r_gp
!
! !REVISION HISTORY:
!  05 Jul 2019 - X. Wang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     REAL(fp)    :: rgs, g1, g2

     rgs  = 0.0e0_fp
     g1   = 0.0e0_fp
     g2   = 0.0e0_fp

     !HOBr + HCl
     rgs = 0.25
     g1 = rgs * hcl_th

     !HOBr + HBr
     rgs = 4.8e-4 * exp(1240_fp/T)
     g2 = rgs * hbr_th

     GAM_HOBr = g1 + g2

     IF ( GAM_HOBr == 0) THEN
        r_gp=0.0_fp
     ELSEIF ( X==1 ) THEN
         r_gp = g1 / GAM_HOBr
     ELSEIF ( X==2 ) THEN
         r_gp = g2 / GAM_HOBr
     ELSE
         r_gp = 0.0_fp
     ENDIF


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
      REAL(fp),       INTENT(IN)  :: Radius   ! Radius (cm)
      REAL(fp),       INTENT(IN)  :: n_air    ! n_air (#/cm)
      INTEGER,        INTENT(IN)  :: X        ! 1=Cl-,2=Br-,3=HSO3-,4=SO3--
      REAL(fp),       INTENT(IN)  :: T        ! Temperature (K)
      REAL(fp),       INTENT(IN)  :: C_Y1     ! Concentration (mol/L)
      REAL(fp),       INTENT(IN)  :: C_Y2     ! Concentration (mol/L)
      REAL(fp),       INTENT(IN)  :: C_Hp     ! Concentration (mol/L)
      TYPE(HetState), POINTER     :: H        ! Hetchem species metadata
!
! !OUTPUT PARAMETERS
!
      ! Reactive uptake coefficient (unitless)
      REAL(fp),       INTENT(OUT) :: GAM_HOBr
      REAL(fp),       INTENT(OUT) :: r_gp
!
! !REVISION HISTORY:
!  30 Nov 2017 - Q.J. Chen   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      REAL(fp) :: ab, gd, M_X
      REAL(fp) :: cavg, H_X, H_HOBr
      REAL(fp) :: k_tot, gb_tot
      REAL(fp) :: k_b1, k_b2
      REAL(fp) :: D_l, ybr2
      REAL(fp) :: l_r, C_Hp1, C_Hp2

      ! Reaction rate coefficient for HOBr + Cl- [M-2 s-1]
      !k_b  = 5.9e+9_fp
      ! (Liu and Margerum, Environ. Sci. Tech., 2001)
      k_b1  = 2.3e+10_fp ! (qjc, 12/28/16)
      ! Reaction rate coefficient for HOBr + Br- [M-2 s-1]
      k_b2  = 1.6e+10_fp

      ! Liquid phase diffusion coefficient [cm2/s] for HOBr
      ! (Ammann et al., Atmos. Chem. Phys., 2013)
      D_l    = 1.4e-5_fp

      ! Henry's law
      H_HOBr = H%HOBr%K0 * con_atm_bar
      H_X    = H_HOBr*EXP(H%HOBr%CR*(1.0e0_fp/T - 1.0e0_fp/H%HOBr%TK))

      ! Molwt of HOBr (kg/mol)
      M_X    = H%HOBr%MW_g * 1e-3_fp

      ! Mass accommodation coefficient
      ab     = 0.6e0_fp

      ! Thermal velocity [cm/s]
      cavg   = SQRT( 8 * RStarG * T / ( pi * M_X) ) * 1.0e2_fp

      ! Follow Roberts et al, (2014)
      C_Hp1  = min(C_Hp, 1.0e-6)
      C_Hp2  = min(C_Hp, 1.0e-2)
      C_Hp1  = max(C_Hp1, 1.0e-9)
      C_Hp2  = max(C_Hp2, 1.0e-6)

      k_tot = k_b1 * C_Y1 * C_Hp1                                            &
            + k_b2 * C_Y2 * C_Hp2

      ! l_r is diffusive length scale [cm];
      ! gb is Bulk reaction coefficient [unitless]
      l_r    = SQRT( D_l / k_tot )
      gb_tot = 4.0_fp * H_X * con_R * T * l_r * k_tot / cavg
      gb_tot = gb_tot * REACTODIFF_CORR( Radius, l_r)


      ! Reactive uptake coefficient [unitless]
      GAM_HOBr = 1.0_fp / ( 1.0_fp/ab  +  1.0_fp/gb_tot )

      ybr2 = 0.41e0*LOG10(C_Y2/C_Y1)+2.25        ! yield of Br2
      ybr2 = MIN(ybr2, 0.9e0)
      ybr2 = MAX(ybr2, TINY(1.0e0))

      IF ( X==1 ) THEN

         IF (C_Y2/C_Y1>5.e-4) THEN
            r_gp = 0.1e0
         ELSE
            r_gp = 1.e0 - ybr2
         ENDIF

      ELSEIF ( X==2 ) THEN

         IF (C_Y2/C_Y1>5.e-4) THEN
            r_gp = 0.9e0
         ELSE
            r_gp = ybr2
         ENDIF

      ENDIF

    END SUBROUTINE GAMMA_HOBr_AER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gamma_ClNO3_AER
!
! !DESCRIPTION: Function GAMMA\_ClNO3\_AER calculates reactive uptake coef.
!               for ClNO3 + Cl-/Br-
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE GAMMA_ClNO3_AER( Radius, n_air, X,         T,      M,         &
                                C_Y1,  C_Y2,   GAM_ClNO3, r_gp,   H         )
!
! !USES:
!
      USE PhysConstants, ONLY : Pi, RStarG
!
! !INPUT PARAMETERS:
!
      ! Radius (cm), n_air (#/cm), and X (1 for Cl and 2 for Br)
      REAL(fp),       INTENT(IN) :: Radius   ! Radius (cm)
      REAL(fp),       INTENT(IN) :: n_air    ! n_air (#/cm)
      INTEGER,        INTENT(IN) :: X        ! 1=Cl-,2=Br-,3=hydrosis
      REAL(fp),       INTENT(IN) :: T        ! Temperature (K)
      REAL(fp),       INTENT(IN) :: C_Y1     ! Concentration (mol/L)
      REAL(fp),       INTENT(IN) :: C_Y2     ! Concentration (mol/L)
      INTEGER,        INTENT(IN) :: M        ! 1=fine, 2= coarse, 3=cloud
      TYPE(HetState), POINTER    :: H        ! Hetchem species metadata

!
! !OUTPUT PARAMETER:
!
      ! Reactive uptake coefficient (unitless)
      REAL(fp)                   :: GAM_ClNO3
      REAl(fp)                   :: r_gp

!
! !REVISION HISTORY:
!  25 Jan 2018 - X. Wang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      REAL(fp) :: ab, M_X, fCl
      REAL(fp) :: cavg, D_l, k_0, k_2, k_tot
      REAL(fp) :: gb1, gb2, gb_tot, gb0, gbr


      ! Mol wt of ClNO3 (kg/mol)
      M_X = H%ClNO3%MW_g * 1.0e-3_fp

      ! Mass accommodation coefficient
      ab = 0.108e0_fp

      ! Calculate gb1 for ClNO3 + Cl-

      ! Following [Deiber et al., 2004], gamma is not significantly different
      ! from ClNO3 + H2O (gamma = 0.0244) independent of Cl- concentration,
      ! but Cl2 rather than HOCl formed. gb2 can be calculated reversely from
      ! gamma and ab:
      gb1 = 0.0_fp

      ! hydrolysis
      gb0 = 0.032_fp

      ! Calculate gb2 for ClNO3 + Br-

      cavg = SQRT(8.0_fp*RStarG*T/(Pi*M_X)) *1.0e2_fp ! thermal velocity (cm/s)

      D_l  = 5.0e-6_fp !cm2 s-1.
      gb2   = 4.0e0_fp * con_R * T * 1.0e6_fp * SQRT(C_Y2*D_l) / cavg ! H*sqrt(kb)=10^6 (M/s)^½ s-1

      k_2 = (1.0e6_fp ** 2.0_fp) * C_Y2 !H2k2Br

      ! Calculate gb1 for ClNO3 + Cl-

      ! Following [Deiber et al., 2004], gamma is not significantly different
      ! from ClNO3 + H2O (gamma = 0.0244) independent of Cl- concentration,
      ! but Cl2 rather than HOCl formed. gb2 can be calculated reversely from
      ! gb1 = gb0 hydrolysis
      gb0 = 4.0e0_fp * con_R * T * 1.2e5_fp * SQRT(D_l) / cavg

      k_0 = 1.2e5_fp ** 2.0_fp !H2k0

      k_tot = k_0 + k_2 !H2(k0+k2Br)

      gb_tot = 4.0e0_fp * con_R * T * SQRT(k_tot*D_l) /cavg

      gbr = k_2/k_tot

      GAM_ClNO3 = 1.0e0_fp / (1.0e0_fp/ab  +  1.0e0_fp/gb_tot)

      IF (M .EQ. 1) THEN
         fCl = C(ind_SALACL) / (C(ind_SALACL) + C(ind_NIT) + C(ind_SO4))
      ELSEIF (M .EQ. 2) THEN
         fCl = 1.0e0_fp
      ELSE
         fCl = 0.0e0_fp
      ENDIF

      IF ( X==1 ) THEN
         r_gp = (1.0e0_fp - gbr) * fCl
      ELSEIF ( X==2 ) THEN
         r_gp = gbr
      ELSEIF ( X==3 ) THEN
         r_gp = (1.0e0_fp - gbr) * (1.0e0_fp-fCl)
      ENDIF

    END SUBROUTINE GAMMA_ClNO3_AER
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Gamma_ClNO3_Ice
!
! !DESCRIPTION: Function GAMMA\_ClNO3\_ICE calculates reactive uptake coef.
!               for ClNO3 + H2O/HCl/HBr in ice clouds
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE GAMMA_ClNO3_ICE( X, T, hno3_th, hcl_th, hbr_th, GAM_ClNO3, r_gp)
!
! !USES:
!
      USE PhysConstants,      ONLY : Pi, RStarG
!
! !OUTPUT PARAMETER:
      ! Reactive uptake coefficient (unitless)
      REAL(fp)                         :: GAM_ClNO3, r_gp
! !INPUT PARAMETERS:
!
      INTEGER,  INTENT(IN)           :: X        ! 1=HCl,2=HBr,3=hydrosis
      REAL(fp), INTENT(IN)           :: T        ! Temperature (K)
      REAL(fp), INTENT(IN)           :: hno3_th, hcl_th, hbr_th
!
! !REVISION HISTORY:
!  02 Jul 2019 - X. Wang - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
     ! Universal gas consatant [bar/(mol/kg)/K]
     REAL(fp)    :: rgs, H2Os, cavg, kks, g1, g2, g3

     rgs  = 0.0e0_fp
     H2Os = 0.0e0_fp
     cavg = 0.0e0_fp
     g1   = 0.0e0_fp
     g2   = 0.0e0_fp
     g3   = 0.0e0_fp

     !ClNO3 + HCl
     rgs = 0.24
     g1 = rgs * hcl_th
     !ClNO3 + HBr
     rgs = 0.56
     g2 = rgs * hbr_th
     !ClNO3 + H2O
     cavg = SQRT(8.0e+0_fp*RStarG*T/(Pi*9.745e-2_fp)) *1.0e2_fp ! thermal velocity (cm/s)
     H2Os = 1e15_fp - 3.0*2.7e14_fp*hno3_th
     kks = 4.0_fp * 5.2e-17_fp * exp(2032_fp/T)
     g3 = 1.0_fp / (1.0_fp/0.5_fp + cavg/(kks*H2Os))

     GAM_ClNO3 = g1 + g2 + g3

     IF ( X==1 ) THEN
         r_gp = g1 / GAM_ClNO3
     ELSEIF ( X==2 ) THEN
         r_gp = g2 / GAM_ClNO3
     ELSEIF ( X==3 ) THEN
         r_gp = g3 / GAM_ClNO3
     ELSE
         r_gp = 0.0_fp
     ENDIF


    END SUBROUTINE GAMMA_ClNO3_ICE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: COTH
!
! !DESCRIPTION: COTH (Hyperbolic cotangent)
!               coth(x) = cosh(x)/sinh(x) = (1 + exp(-2x))/(1 - exp(-2x))
!\\
!\\
! !INTERFACE:
!
      REAL(fp) FUNCTION COTH( X)
!
! !INPUT PARAMETERS:
!
      REAL(fp),         INTENT(IN)  :: X           ! The argument
!
! !REVISION HISTORY:
!  24 Sept 2015 - J. Schmidt - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      REAL(fp)                 :: exp_temp

      exp_temp = EXP(-2.0e0_fp*X)
      COTH = (1.0e0_fp + exp_temp)/(1.0e0_fp - exp_temp)

      END FUNCTION COTH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: REACTODIFF_CORR
!
! !DESCRIPTION: REACTODIFF\_CORR
!     Correction =  COTH( x ) - ( 1/x )
!              x = radius / l
!     Correction approaches 1 as x becomes large, corr(x>1000)~1
!     Correction approaches x/3 as x goes towards 0
!\\
!\\
! !INTERFACE:
!
      REAL(fp) FUNCTION REACTODIFF_CORR( radius, l)
!
! !INPUT PARAMETERS:
!
      REAL(fp),         INTENT(IN)  :: radius, l           ! [cm] and [cm]
!
! !REVISION HISTORY:
!  14 Oct 2013 - J. Schmidt - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!

      REAL(fp)                 :: x

      x = radius / l

      IF (x<0.0e0_fp) THEN
         PRINT *, 'ERROR x<0, particle radius or C_Y is neg!'
      ELSEIF (x>1.0e3_fp) THEN
         REACTODIFF_CORR = 1.0e0_fp
      ELSEIF (x<1.0e-1_fp) THEN
         REACTODIFF_CORR = x/3.0e0_fp
      ELSE
         REACTODIFF_CORR = COTH(x) - (1.0e0_fp/x)
      ENDIF


      RETURN

      END FUNCTION REACTODIFF_CORR
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetBrNO3
!
! !DESCRIPTION: Sets the hydrolysis rate for BrNO3 using Johan Schmidt's
!  updated code.
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETBrNO3( denAir, TK, CldFr, H ) RESULT( HET_BrNO3 )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: denAir     ! Density of air (#/cm3)
      REAL(fp),       INTENT(IN) :: TK         ! Temperature (K)
      REAL(fp),       INTENT(IN) :: CldFr      ! Cloud fraction
      TYPE(HetState), POINTER    :: H          ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(fp)                   :: HET_BrNO3  ! Hydrolysis rate
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
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE, GAM, ab, XSQM

      ! Initialize
      HET_BrNO3    = 0.0_fp
      ADJUSTEDRATE = 0.0_fp
      XSTKCF       = 0.0_fp
      GAM          = 0.0_fp
      ab           = 0.063_fp
      XSQM         = SQRT( H%BrNO3%MW_g )

      ! Calculate temperature dependent reaction uptake rate
      ! Based on Deiber et al. (2004)
      GAM = 0.0021_fp * TK - 0.561
      GAM = MAX( GAM, 1.0e-20_fp )

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE
         XSTKCF = 0e+0_fp
         ! Get the aerosol type
         IF ( N == 8 ) THEN
            ! sulfate aerosol
            !XSTKCF = 1.0_fp / (1.0_fp/ab + 1.0_fp/GAM)
            XSTKCF = GAM
         ELSEIF ( (N == 11) .OR. ( N == 12) ) THEN
            ! 2 modes of sea-salt
            !XSTKCF = 1.0_fp / (1.0_fp/ab + 1.0_fp/GAM)
            XSTKCF = GAM
         ELSEIF ( N == 13 ) THEN
            ! SSA/STS
            XSTKCF = KHETI_SLA(6)
         ELSEIF ( N == 14 ) THEN
            ! Ice/NAT PSC
            IF (NATSURFACE) THEN
               XSTKCF = 0.001e+0_fp
            ELSE
               XSTKCF = 0.3e+0_fp
            ENDIF
         ELSE
            XSTKCF = 0e+0_fp
         ENDIF

         IF (N.eq.13) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K( (1-CldFr)*XAREA(N), XRADI(N), denAir,       &
                                 XSTKCF,             XTEMP,    XSQM         )
         ENDIF

         ! Add to overall reaction rate
         HET_BrNO3 = HET_BrNO3 + ADJUSTEDRATE
      END DO

    END FUNCTION HETBrNO3
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
      REAL(fp), INTENT(IN) :: A, B
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms.
!
! !REVISION HISTORY:
!  29 Jan 2016 - M. Sulprizio- Initial version, adapted from code previously
!                              in calcrate.F
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL  :: DO_EDUCT
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE
      Real(fp) :: XSQM

      ! Initialize
      kISum          = 0.0_fp
      ADJUSTEDRATE   = 0.0_fp
      XSTKCF         = 0.0_fp
      XSQM=SQRT(A)

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Assume zero unless proven otherwise
         XSTKCF = 0e+0_fp

!         IF (N.eq.8) THEN
!            XSTKCF = 0.1e-4_fp ! Sulfate
!         ELSEIF ( STRATBOX ) THEN
!            IF (N.eq.13) THEN
	 ! restore limitation to stratosphere - TMS 17/04/10
         IF  ( STRATBOX ) THEN
            IF (N.eq.8) THEN
               XSTKCF = 0.1e-4_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(4)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.2e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.3e+0_fp ! Ice
               ENDIF
            ENDIF
         ENDIF

         IF (STRATBOX.and.(N.eq.13)) THEN
            ! Calculate for stratospheric liquid aerosol
            ! Note that XSTKCF is actually a premultiplying
            ! factor in this case, including c-bar
            ADJUSTEDRATE = XAREA(N) * XSTKCF
         ELSE IF (XSTKCF .GT. 0.0e+0_fp) THEN
            ! Reaction rate for surface of aerosol
            ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                               XSQM)
         ENDIF

         ! Add to overall reaction rate
         kISum = kISum + ADJUSTEDRATE

      END DO

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
    FUNCTION HETBrNO3_HCl( A, B, Input_Opt ) RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      ! Rate coefficients
      REAL(fp), INTENT(IN) :: A, B
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !RETURN VALUE:
!
      REAL(fp)             :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Jan 2016 - M. Sulprizio- Initial version, adapted from code previously
!                              in calcrate.F
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      kISum         = 0.0_fp
      ADJUSTEDRATE  = 0.0_fp
      XSTKCF        = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         ! Default to zero
         XSTKCF = 0.0e+0_fp

	 ! restore limitation to stratosphere - TMS 17/04/10
	 IF ( STRATBOX ) THEN
	    IF (N.eq.8) THEN
               XSTKCF = 0.9_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(7)
               ! Turn off het rates in stratosphere
               IF ( Input_Opt%TurnOffHetRates ) XSTKCF = 0.0e+0_fp       ! EF turn off BrNO3 + HCl on sulfates as in GMI
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.3_fp ! NAT   !%%% what is the difference?
               ELSE
                  XSTKCF = 0.3_fp ! Ice
               ENDIF
            ENDIF
         ENDIF

         IF (XStkCf.gt.0.0e+0_fp) THEN
            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            ! Add to overall reaction rate
            kISum = kISum + ADJUSTEDRATE
         ENDIF

      END DO

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
      REAL(fp),       INTENT(IN) :: A, B       ! Rate coefficients
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Jan 2016 - M. Sulprizio- Initial version, adapted from code previously
!                              in calcrate.F
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      kISum         = 0.0_fp
      ADJUSTEDRATE  = 0.0_fp
      XSTKCF        = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         XSTKCF        = 0.0_fp

         ! For UCX-based mechanisms only consider PSC reactions in strat
         ! restore limitation to stratosphere - TMS 17/04/10
         IF ( Input_Opt%LUCX .and. STRATBOX) THEN
	    IF (N.eq.8) THEN
	       XSTKCF = 0.8e+0_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(8)
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.1e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.2e+0_fp ! Ice
               ENDIF
            ENDIF
         ENDIF

         IF (XStkCf.gt.0.0e+0_fp) THEN
            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            ! Add to overall reaction rate
            kISum = kISum + AdjustedRate

         ENDIF

      END DO

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
      REAL(fp),       INTENT(IN) :: A, B       ! Rate coefficients
      TYPE(OptInput), INTENT(IN) :: Input_Opt  ! Input Options object
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum
!
! !REMARKS:
!
! !REVISION HISTORY:
!  29 Jan 2016 - M. Sulprizio- Initial version, adapted from code previously
!                              in calcrate.F
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N
      REAL(fp) :: XSTKCF, ADJUSTEDRATE

      ! Initialize
      kISum         = 0.0_fp
      ADJUSTEDRATE  = 0.0_fp
      XSTKCF        = 0.0_fp

      ! Loop over aerosol types
      DO N = 1, NAEROTYPE

         XSTKCF        = 0.0_fp

         ! For UCX-based mechanisms only consider PSC reactions in strat
         ! restore limitation to stratosphere - TMS 17/04/10
         IF ( Input_Opt%LUCX .and. STRATBOX ) THEN
	    IF (N.eq.8) THEN
 	       XSTKCF = 0.8e+0_fp ! Sulfate
            ELSEIF (N.eq.13) THEN
               XSTKCF = KHETI_SLA(9)
               ! Turn off het rates in stratosphere
               IF ( Input_Opt%TurnOffHetRates ) XSTKCF = 0.0e+0_fp       ! EF turn off HOCl + HBr on sulfates as in GMI
            ELSEIF (N.eq.14) THEN
               IF (NATSURFACE) THEN
                  XSTKCF = 0.3e+0_fp ! NAT
               ELSE
                  XSTKCF = 0.3e+0_fp ! Ice
               ENDIF
               ! Turn off het rates in stratosphere
               IF ( Input_Opt%TurnOffHetRates ) XSTKCF = 0.0e+0_fp       ! EF turn off HOCl + HBr on ALL PSCs as in GMI
            ENDIF
         ENDIF

         IF (XStkCf.gt.0.0e+0_fp) THEN
            IF (N.eq.13) THEN
               ! Calculate for stratospheric liquid aerosol
               ! Note that XSTKCF is actually a premultiplying
               ! factor in this case, including c-bar
               ADJUSTEDRATE = XAREA(N) * XSTKCF
            ELSE
               ! Reaction rate for surface of aerosol
               ADJUSTEDRATE=ARSL1K(XAREA(N),XRADI(N),XDENA,XSTKCF,XTEMP, &
                                  (A**0.5_FP))
            ENDIF

            ! Add to overall reaction rate
            kISum = kISum + AdjustedRate

         ENDIF

      END DO

    END FUNCTION HETHOCl_HBr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HetHOCl_TCld
!
! !DESCRIPTION: Sets the rate of the multiphase reaction HOCl + Cl-/S(IV) in
!  troposphere cloud
!\\
!\\
! !INTERFACE:
!
    FUNCTION HETHOCl_TCld( denAir,    rLiq,      rIce,      ALiq,            &
                           AIce,      VAir,      TK,        CldFr,           &
                           hConc_Sul, hConc_LCl, hConc_ICl, clConc_A,        &
                           clConc_C,  clConc_g,  hso3Conc,  so3Conc,         &
                           hcl_th,    X,         H                         ) &
                         RESULT( kISum )
!
! !INPUT PARAMETERS:
!
      REAL(fp),       INTENT(IN) :: denAir      ! Density of air (#/cm3)
      REAL(fp),       INTENT(IN) :: rLiq        ! Radius of liquid cloud
                                                !  droplets (cm)
      REAL(fp),       INTENT(IN) :: rIce        ! Radius of ice cloud
                                                !  crystals (cm)
      REAL(fp),       INTENT(IN) :: ALiq        ! Area of liquid cloud
                                                !  droplets (cm2/cm3)
      REAL(fp),       INTENT(IN) :: AIce        ! Area of ice cloud
                                                !  crystals (cm2/cm3)
      REAL(fp),       INTENT(IN) :: VAir        ! Box volume (cm3)
      REAL(fp),       INTENT(IN) :: TK          ! Temperature (K)
      REAL(fp),       INTENT(IN) :: CldFr       ! Updated to new cloud process
      REAL(fp),       INTENT(IN) :: hConc_Sul   ! Sulfate H+ concentration
      REAL(fp),       INTENT(IN) :: hConc_LCl   ! Liquid cloud H+ concentration
      REAL(fp),       INTENT(IN) :: hConc_ICl   ! Ice cloud H+ concentration
      REAL(fp),       INTENT(IN) :: clConc_A    ! Fine Chloride
                                                !  concentration (mol/L)
      REAL(fp),       INTENT(IN) :: clConc_C    ! Coarse Chloride
                                                !  concentration (mol/L)
      REAL(fp),       INTENT(IN) :: clConc_g
      REAL(fp),       INTENT(IN) :: hso3Conc    ! HSO3-  concentration (mol/L)
      REAL(fp),       INTENT(IN) :: so3Conc     ! SO3--  concentration (mol/L)
      REAL(fp),       INTENT(IN) :: hcl_th
      INTEGER,        INTENT(IN) :: X           ! 1=fineCl;2=coarseCl;
                                                ! 3=HCl;4=HSO3;5=SO3
      TYPE(HetState), POINTER    :: H           ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum       ! Reaction rate
!
! !REVISION HISTORY:
!  08 Nov 2018 - X. Wang      - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER  :: N,      Y
      REAL(fp) :: X1,     X2,   GAM_HOCl, r_gp
      REAL(fp) :: clConc, r_ac, r1,       r2,   rgs

      ! Initialize return value
      kISum = 0.0_fp

      ! Return if we are in the stratosphere
      ! This skips unneccessary computations & function calls
      IF ( StratBox ) RETURN

      ! Continue initializing
      X1    = 0.0_fp
      X2    = 0.0_fp

      ! Cloud halide concentration, put fine and coarse mode together
      clConc = clConc_A + clConc_C + clConc_g

      If (X == 1) THEN
         Y = 1
         r_ac = clConc_A / clConc
      ELSEIF (X == 2) THEN
         Y = 1
         r_ac = clConc_C / clConc
      ElSEIF (X == 3) THEN
         Y = 1
         r_ac = clConc_g / clConc
      ELSEIF (X == 4) THEN
         Y = 2
         r_ac = hso3Conc / (hso3Conc + so3Conc)
      ELSEIF (X == 5) THEN
         Y = 2
         r_ac = so3Conc / (hso3Conc + so3Conc)
      ENDIF

      CALL Gamma_HOCl_CLD( rLiq,     denAir,   Y,       TK,                  &
                           clConc,   hso3Conc, so3Conc, hConc_LCl,           &
                           GAM_HOCl, r_gp,     H                            )
      X1 = GAM_HOCl
      r1 = r_gp*r_ac

      !For uptake on ice crystals
      rgs = 0.22_fp
      X2 = rgs * hcl_th
      IF (X == 3) THEN
         r2 = 1.0_fp
      ELSE
         r2 = 0.0_fp
      ENDIF

      kISum = CloudHet( CldFr, Aliq,   Aice,        rLiq, rIce,              &
                        TK,    denAir, H%HOCl%MW_g, X1,   X2,   r1,   r2          )

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
      REAL(fp),       INTENT(IN) :: denAir    ! Density of air (#/cm3)
      REAL(fp),       INTENT(IN) :: rAer      ! Radius of aerosol (cm)
      REAL(fp),       INTENT(IN) :: AAer      ! Area of aerosol (cm2/cm3)
      REAL(fp),       INTENT(IN) :: alkAer    ! Aerosol alkalinity (eqn)
      REAL(fp),       INTENT(IN) :: TK        ! Temperature (K)
      REAL(fp),       INTENT(IN) :: hconc     ! H+ (M)
      REAL(fp),       INTENT(IN) :: clconc    ! Cl- (M)
      TYPE(HetState), POINTER    :: H         ! Hetchem species metadata
!
! !RETURN VALUE:
!
      REAL(fp)                   :: kISum     ! Rate for HOCl on Cl- aerosol
!
! !REVISION HISTORY:
!  25 Jan 2018 - X. Wang - Initial vesion
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: XSTKCF, XSQM

      ! Initialize
      kISum = 0.0_fp
      XSQM  = SQRT( H%HOCl%MW_g )

      ! Reaction can only proceed on acidic aerosol
      IF ( alkAer > 0.05_fp ) THEN
         XStkCf = 0.0_fp
      ELSE
         XStkCf = GAMMA_HOCl_AER( rAer, denAir, TK, hconc, clconc, H )
      ENDIF

      ! Reaction rate
      kISum = ARSL1K( AAer, rAer, denAir, XSTKCF, XTEMP, XSQM )

    END FUNCTION HETHOCl_SS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: N2O5
!
! !DESCRIPTION: Internal function N2O5 computes the GAMMA sticking factor
!  for N2O5 hydrolysis. (mje, bmy, 8/7/03)
!\\
!\\
! !INTERFACE:
!
    FUNCTION N2O5( AEROTYPE, TEMP, RH ) RESULT( GAMMA )
!
! !INPUT PARAMETERS:
!
      INTEGER,   INTENT(IN) :: AEROTYPE  ! Denoting aerosol type (cf FAST_JX)
      REAL(fp),  INTENT(IN) :: TEMP      ! Temperature [K]
      REAL(fp),  INTENT(IN) :: RH        ! Relative humidity [%]
!
! !RETURN VALUE:
!
      REAL(fp)              :: GAMMA

!
! !REMARKS:
!  Taken from the old SMVGEAR function calcrate.F.
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
      REAL(fp) :: RH_P, FACT, TTEMP

      !=================================================================
      ! N2O5 begins here!
      !=================================================================

      ! RH percent max = 100%
      RH_P  = MIN( RH, 100e0_fp )

      ! Default value
      GAMMA = 0.01e+0_fp

      ! Special handling for various aerosols
      SELECT CASE ( AEROTYPE )

         !----------------
         ! Dust
         !----------------
         CASE ( 1, 2, 3, 4, 5, 6, 7 )
            GAMMA = 0.02e+0_fp

         !----------------
         ! Sulfate
         !----------------
         CASE ( 8 )
            ! This is not used when considering Cl-, XW
!            call ERROR_STOP( 'N2O5 update on sulfate should not be calculated here','gckpp_HetRates.F90')
            GAMMA = 0.00e+0_fp
         !----------------
         ! Black Carbon
         !----------------
         CASE ( 9 )

             ! From IUPAC
             GAMMA = 0.005e+0_fp

         !----------------
         ! Organic Carbon
         !----------------
         CASE ( 10 )

            ! Escoria et al (2010)
            IF ( RH_P < 30e+0_fp ) THEN
               GAMMA = 0.6e-4_fp
            ELSE
               GAMMA = 1.5e-4_fp
            ENDIF

         !----------------
         ! Sea salt
         ! accum & coarse
         !----------------
         CASE ( 11, 12 )
            ! This is not used when considering Cl-, XW
!            call ERROR_STOP( 'N2O5 update on sea salt should not be calculated here','gckpp_HetRates.F90' )
            GAMMA = 0.00e+0_fp

         !----------------
         ! Strat. aerosols
         !----------------
         CASE ( 13, 14 )

            ! These are handled outside this routine - something
            ! is wrong if AEROTYPE=13 or 14 reaches this point
            WRITE (6,*) 'Stratospheric aerosols should not '
            WRITE (6,*) 'be passed to general N2O5 het. '
            WRITE (6,*) 'chem. subroutine'
            WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
            CALL GEOS_CHEM_STOP

         !----------------
         ! Default
         !----------------
         CASE DEFAULT
            WRITE (6,*) 'Not a suitable aerosol surface '
            WRITE (6,*) 'for N2O5 hydrolysis'
            WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
            CALL GEOS_CHEM_STOP

      END SELECT

    END FUNCTION N2O5
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: HO2
!
! !DESCRIPTION: Function HO2 computes the GAMMA reaction probability
!  for HO2 loss in aerosols based on the recommendation of
!  Thornton, Jaegle, and McNeill, "Assessing Known Pathways For HO2 Loss in
!  Aqueous Atmospheric Aerosols: Regional and Global Impacts on Tropospheric
!  Oxidants" J. Geophys. Res.,  doi:10.1029/2007JD009236, 2008
!\\
!\\
! !INTERFACE:
!
      FUNCTION HO2( RADIUS,          TEMP,     DENAIR,       &
                    SQM,             HO2DENS,  AEROTYPE,     &
                    CONTINENTAL_PBL, Input_Opt            )  &
                    RESULT( GAMMA )
!
! !USES:
!
      USE Input_Opt_Mod, ONLY : OptInput
      Use PhysConstants, ONLY : AVO, RGASLATM
!
! !INPUT PARAMETERS:
!
      ! Arguments
      REAL(fp),       INTENT(IN) :: RADIUS          ! Aerosol radius [cm]
      REAL(fp),       INTENT(IN) :: TEMP            ! Temperature [K]
      REAL(fp),       INTENT(IN) :: DENAIR          ! Air density [molec/cm3]
      REAL(fp),       INTENT(IN) :: HO2DENS         ! HO2 density [molec/cm3]
      REAL(fp),       INTENT(IN) :: SQM             ! Square root of MW [g/mol]
      INTEGER,        INTENT(IN) :: AEROTYPE        ! Aerosol type (cf FAST-JX)
      INTEGER,        INTENT(IN) :: CONTINENTAL_PBL ! Flag set to 1 if the box
                                                    !  box is located in the
                                                    !  continenal boundary
                                                    !  layer, otherwise 0.
                                                    !  Also check for ICE/SNOW
                                                    !  (to disable this at
                                                    !  high latitudes).
      TYPE(OptInput), INTENT(IN) :: Input_Opt       ! Input Options object
!
! !RETURN VALUE:
!
      REAL(fp)                   :: GAMMA           ! Reaction probability

! !REMARKS:
!  Taken from the old SMVGEAR routine calcrate.F.
!  Gamma(HO2) is a function of aerosol type, radius, temperature.

!  References:
!  ---------------------------------------------------------------
!  (1) Jacob, D.J., Heterogeneous chemistry and tropospheric ozone,
!       Atmos. Environ., 34, 2131-2159, 2000. [full text (pdf)]
!  (2) J. Mao, Fan, S., Jacob, D. J., and Travis, K. R.: Radical
!       loss in the atmosphere from Cu-Fe redox coupling in aerosols,
!       Atmos. Chem. Phys., 13, 509-519, doi:10.5194/acp-13-509-2013,
!       2013.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)             :: ALPHA
      REAL(fp)             :: delG, Keq, w, H_eff
      REAL(fp)             :: A1, B1, k1, k2, A, B, C
      REAL(fp)             :: kaq, kmt, o2_ss, fluxrxn, DFKG
      REAL(fp)             :: TEST
!
! !DEFINED PARAMETERS:
!
      ! Ideal gas constant [atm cm3/mol/K]
      REAL(fp),  PARAMETER :: Raq = RGASLATM * 1e+3_fp

      !=================================================================
      ! HO2 begins here!
      !=================================================================

      ! Default value
      GAMMA = 0.0e+0_fp

      ! Error check
      IF (RADIUS.le.1e-30_fp) THEN
         RETURN
      ENDIF

      ! Special handling for various aerosols
      SELECT CASE ( AEROTYPE )

         !----------------
         ! Dust
         !----------------
         CASE ( 1, 2, 3, 4, 5, 6, 7 )

            ! Assume default gamma=0.1 on dust aerosols
            ! This is tentative as no lab measurements presently exist
            ! for gamma(HO2) on dust aerosols. We assume the rate to
            ! be fast on dust aerosols as transition metal ion induced
            ! chemistry is likely to occur in a thin aqueous surface layer.
            GAMMA = 0.1e+0_fp

         !----------------
         ! For Sulfate(8), Black Carbon (9), Organic Carbon (10),
         ! Sea-salt accum & coarse (11,12) calculate the
         ! reaction probability due to self reaction
         ! by using the algebraic expression in Thornton et al.  (2008)
         ! (equation 7) which is a function of temperature, aerosol radius,
         ! air density and HO2 concentration.
         !
         ! Transition metal ions (such as copper and iron) in sea-salt and
         ! carbonaceous aerosols are complexed to ligands and/or exist at
         ! a concentration too low to catalyze HO2 loss efficiently, so we
         ! apply the HO2 self reaction expression directly for these aerosols.
         !
         ! In the case of sulfate aerosol, the aerosols likely
         ! contain copper in the continental boundary layer and
         ! HO2 uptake proceeds rapidly. To account for the metal catalyzed
         ! uptake, we assume gamma(HO2)=0.07 (in the mid-range of the recommended
         ! 0.04-0.1 by Thornton et al, based on observed copper concentrations
         ! in the US boundary layer). Outside the continental boundary layer, we
         ! use the HO2-only algebraic expression.
         !
         ! SDE 04/18/13: Added stratospheric sulfur aerosols
         !
         !----------------
         CASE ( 8, 9, 10, 11, 12, 13 )

            ! Mean molecular speed [cm/s]
            w = 14550.5e+0_fp * sqrt(TEMP/(SQM*SQM))

            ! DFKG = Gas phase diffusion coeff [cm2/s]
            DFKG  = 9.45E+17_fp/DENAIR * SQRT(TEMP) *      &
                    SQRT(3.472E-2_fp + 1.E+0_fp/(SQM*SQM))

            !calculate T-dependent solubility and aq. reaction rate constants
            ! hydronium ion concentration
            ! A1 = 1.+(Keq/hplus)
            ! with Keq = 2.1d-5 [M], Equilibrium constant for
            ! HO2aq = H+ + O2- (Jacob, 2000)
            !      hplus=10.e+0_fp^(-pH), with pH = 5
            ! B1 = Req * TEMP
            ! with Req = 1.987d-3 [kcal/K/mol], Ideal gas constant
            ! Note that we assume a constant pH of 5.
            A1 = 1.+ (2.1e-5_fp / (10.e+0_fp**(-5) ) )
            B1 = 1.987e-3_fp * TEMP

            ! Free energy change for solvation of HO2 (Hanson 1992, Golden 1991)
            ! in [kcal/mol]:
            ! delG = -4.9-(TEMP-298e+0_fp)*delS
            ! with delS=-0.023  [kcal/mol/K],  Entropy change for solvation of HO2
            delG  = -4.9e+0_fp - (TEMP-298.e+0_fp) * (-0.023)
            H_eff = exp( -delG / B1 ) * A1

            ! Estimated temp dependent value for HO2 + O2- (k1) and
            ! HO2+HO2 (see Jacob 1989)
            k1  =   1.58e+10_fp * exp( -3. / B1 )
            k2  =   2.4e+9_fp   * exp( -4.7 / B1 )
            kaq = ( k1 * (A1 - 1.e+0_fp) + k2) / (A1**2)

            ! Calculate the mass transfer rate constant and s.s. conc. of
            ! total HO2 in the aqueous phase:
            ! kmt = (RADIUS/DFKG + 4e+0_fp/w/alpha)^(-1)
            ! with alpha = mass accomodation coefficient, assumed
            ! to be 1 (Thornton et al.)
            kmt = 1.e+0_fp/( RADIUS/DFKG + 4e+0_fp/w/1. )

            !use quadratic formula to obtain [O2-] in particle of radius RADIUS
            A = -2e+0_fp * kaq
            B = -3e+0_fp * kmt / RADIUS / (H_eff * 0.082 * TEMP)
            C =  3e+0_fp * kmt * HO2DENS * 1000e+0_fp / RADIUS / AVO

            ! Error check that B^2-(4e+0_fp*A*C) is not negative
            TEST= B**2-(4e+0_fp*A*C)
            IF ( TEST < 0e+0_fp ) THEN
                GAMMA = 0e+0_fp
            ELSE
                ! Calculate the concentration of O2- in the aerosol
                o2_ss= ( -B  -sqrt(B**2-(4e+0_fp*A*C)) )/(2e+0_fp*A)

                ! Calculate the reactive flux
                fluxrxn = kmt*HO2DENS - o2_ss*AVO*kmt/H_eff/Raq/TEMP

                IF ( fluxrxn <= 0e0_fp ) THEN
                   GAMMA = 0e+0_fp
                ELSE
                   ! Gamma for HO2 at TEMP, ho2, and RADIUS given
                   GAMMA = 1./( ( ( HO2DENS/fluxrxn ) -              &
                                  ( RADIUS/DFKG ) ) * w / 4.e+0_fp )
                ENDIF
            ENDIF
            ! For sulfate aerosols, check whether we are in
            ! the continental boundary layer, in which case
            ! copper catalyzed HO2 uptake likely dominates and
            ! speeds up the reaction: we assume gamma=0.07,
            ! which is in the middle of the 0.04-0.1 range recommended
            ! by Thornton et al. (2008)
            !
            IF ( AEROTYPE == 8 .and. CONTINENTAL_PBL == 1) THEN
                GAMMA = 0.07
            ENDIF

         !----------------
         ! NAT/ice (SDE 04/18/13)
         !----------------
         CASE ( 14 )

            GAMMA = 0.e+0_fp

         !----------------
         ! Default
         !----------------
         CASE DEFAULT
            WRITE (6,*) 'Not a suitable aerosol surface '
            WRITE (6,*) 'for HO2 uptake'
            WRITE (6,*) 'AEROSOL TYPE =',AEROTYPE
            CALL GEOS_CHEM_STOP

      END SELECT

      ! If negative value is calculated, set it to zero
      IF ( GAMMA  <= 0e+0_fp ) GAMMA = 0e+0_fp

      ! This is for the improved HO2 uptake (J. Mao)
      GAMMA = Input_Opt%GAMMA_HO2

    END FUNCTION HO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: EPOXUPTK
!
! !DESCRIPTION: Function EPOXUPTK computes the GAMMA sticking factor for
! EPOXUPTK hydrolysis to form 2-methyltetrols (AITET). (eam, 2014).
!\\
!\\
! !INTERFACE:
!
    FUNCTION EPOXUPTK( AERAREA, AERRAD,  TEMP,  SQMW,              &
                       HENRY,   KHPLUS,  HPLUS, KNUC,  SULF, NITR, &
                       KGACID,  BISULF,  KHYDRO ) &
             RESULT( GAMMA )
!
! !USES:
!
      USE Input_Opt_Mod, ONLY : OptInput
      USE PhysConstants, ONLY : RGASLATM
      USE ERROR_MOD,     ONLY : IT_IS_NAN
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN) :: AERRAD   ! Aerosol radius [cm]
      REAL(fp), INTENT(IN) :: AERAREA  ! Aerosol surf. area [cm2/cm3]
      REAL(fp), INTENT(IN) :: TEMP     ! Temperature [K]
      REAL(fp), INTENT(IN) :: SQMW     ! Square root of the molecular weight
      REAL(fp), INTENT(IN) :: HENRY    ! Henry's Law constant [M/atm]
      REAL(fp), INTENT(IN) :: KHPLUS   ! 1st order rxn rate (acid-catalyzed ring
                                       ! opening)
      REAL(fp), INTENT(IN) :: HPLUS    ! Proton activity [unitless] and [H+] [M]
      REAL(fp), INTENT(IN) :: KNUC     ! 1st order rxn rate due to specific
                                       ! nucleophiles (SO4, NO3)
      REAL(fp), INTENT(IN) :: SULF     ! Sulfate concentration [M]
      REAL(fp), INTENT(IN) :: NITR     ! Nitrate concentration [M]
      REAL(fp), INTENT(IN) :: KGACID   ! 1st order rxn rate due to general acids
                                       ! (bisulfate in this case)
      REAL(fp), INTENT(IN) :: BISULF   ! Bisulfate concentration [M]
      REAL(fp), INTENT(IN) :: KHYDRO   ! Hydrolysis rate of alkylnitrates [1/s]
!
! !RETURN VALUE:
!
      REAL(fp)             :: GAMMA    ! Reaction probability

! !REMARKS:
! Calculation is only done for inorganic aqueous phase aerosols
!                                                                             .
! This calculation uses the parameterization of Gaston et al., EST, 2014.
!                                                                             .
! Redistribution of products (e.g. AITET) to yield organosulfates and
! organonitrates is done in SOA_CHEMISTRY in carbon_mod.F.
! This is only done for IEPOX and IMAE if it's an SOA simulation
!
! !REVISION HISTORY:
!  15 Jun 2017 - M. Sulprizio- Initial version based on calcrate.F from E.Marais
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Local variables
      REAL(fp)             :: AERVOL   ! Aerosol volume [cm3/cm3]
      REAL(fp)             :: KPART    ! Particle-phase reaction rate [1/s]
      REAL(fp)             :: XMMS     ! Mean molecular speed [cm/s]
      REAL(fp)             :: VAL1, VAL2, VAL3  ! Terms for calculating GAMMA
      REAL(fp)             :: VALTMP
!
! !DEFINED PARAMETERS:
!
      ! Gas-phase diffusion constant [cm2/s]:
      REAL(fp), PARAMETER  :: DIFF_N2O5_STD = 1.0e-1_fp

      ! Mass accommodation coefficient [unitless]:
      REAL(fp), PARAMETER  :: MACOEFF = 1.0e-1_fp

      !=================================================================
      ! EPOXUPTK begins here!
      !=================================================================

      ! Initialize
      GAMMA  = 0.0_fp
      AERVOL = 0.0_fp
      KPART  = 0.0_fp
      XMMS   = 0.0_fp
      VAL1   = 0.0_fp
      VAL2   = 0.0_fp
      VAL3   = 0.0_fp
      VALTMP = 0.0_fp

      ! Calculate aerosol volume (use formula in aerosol_mod.F):
      AERVOL = (AERAREA * AERRAD)/3.0e+0_fp

      ! Calculate mean molecular speed [cm/s]:
      XMMS = SQRT( (2.117e+8_fp * TEMP) / (SQMW * SQMW) )

      ! Calculate first-order particle-phase reaction rate:
      ! (assume [H+] = proton activity)
      ! KHYDRO is only important for alkylnitrates (not currently used).
      KPART = ( KHPLUS*HPLUS )               + &
              ( KNUC*HPLUS*( NITR + SULF ) ) + &
              ( KGACID*BISULF )              + &
              ( KHYDRO )

      ! Calculate the first uptake parameterization term:
      VAL1 = ( AERRAD * XMMS )/( 4.e+0_fp * DIFF_N2O5_STD )

      ! Calculate the second uptake parameterization term:
      VAL2 = ( 1.e+0_fp/MACOEFF )

      ! Calculate the third uptake parameterization term:
      IF ( AERAREA > 0.0_fp .and. XMMS > 0.0_fp ) THEN
         VALTMP = ( 4.e+0_fp * AERVOL * RGASLATM * TEMP * HENRY * KPART ) / &
                  ( AERAREA * XMMS )
      ENDIF
      IF ( VALTMP .GT. 0 ) THEN
         VAL3 = 1.e+0_fp / VALTMP
      ELSE
         VAL3 = 0.0e+0_fp
      ENDIF

      ! Account for small reaction rates:
      IF ( KPART .LT. 1.e-8_fp ) THEN

         GAMMA = TINY(1e+0_fp)

      ELSE

         ! Calculate the uptake coefficient:
         GAMMA = 1.e+0_fp/( VAL1 + VAL2 + VAL3 )

      ENDIF

      ! Fail safes for negative, very very small, and NAN GAMMA values:
      IF ( GAMMA  .lt. 0.0e+0_fp )    GAMMA = TINY(1e+0_fp)
      IF ( IT_IS_NAN( GAMMA ) )       GAMMA = TINY(1e+0_fp)
      IF ( GAMMA .lt. TINY(1e+0_fp) ) GAMMA = TINY(1e+0_fp)

      END FUNCTION EPOXUPTK
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
    SUBROUTINE CLD_PARAMS( I,    J,    L,    State_Met, rLiq,                &
                           ALiq, VLiq, rIce, AIce,      VIce, CLDFr         )
!
! !USES:
!
      USE GcKpp_Global,  ONLY : DENAIR => XDENA
      USE GcKpp_Global,  ONLY : T      => TEMP
      USE GcKpp_Global,  ONLY : QI     => QIce
      USE GcKpp_Global,  ONLY : QL     => QLiq
      USE GcKpp_Global,  ONLY : Vair
      USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)  :: I         ! Longitude index
      INTEGER,        INTENT(IN)  :: J         ! Latitude  index
      INTEGER,        INTENT(IN)  :: L         ! Altitude  index
      TYPE(MetState), INTENT(IN)  :: State_Met ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      REAL(fp),       INTENT(OUT) :: rLiq     ! Radius of liquid cloud droplets (cm)
      REAL(fp),       INTENT(OUT) :: rIce     ! Radius of ice cloud crystals (cm)
      REAL(fp),       INTENT(OUT) :: ALiq     ! Sfc area of liq. cloud (cm2/cm3)
      REAL(fp),       INTENT(OUT) :: AIce     ! Sfc area of ice cloud (cm2/cm3)
      REAL(fp),       INTENT(OUT) :: VLiq     ! Volume of liq. cloud condensate (cm3/cm3)
      REAL(fp),       INTENT(OUT) :: VIce     ! Volume of ice cloud condensate (cm3/cm3)
      REAL(fp),       INTENT(OUT) :: CLDFr    ! cloud fraction
!
! !REMARKS:
!  References:
!   Heymsfield, A. J., Winker, D., Avery, M., et al. (2014). Relationships between ice water content and volume extinction coefficient from in situ observations for temperatures from 0° to –86°C: implications for spaceborne lidar retrievals. Journal of Applied Meteorology and Climatology, 53(2), 479–505. https://doi.org/10.1175/JAMC-D-13-087.1
!
!   Schmitt, C. G., & Heymsfield, A. J. (2005). Total Surface Area Estimates for Individual Ice Particles and Particle Populations. Journal of Applied Meteorology, 44(4), 467–474. https://doi.org/10.1175/JAM2209.1
!
! !REVISION HISTORY:
!  21 Dec 2016 - S. D. Eastham - Adapted from CLD1K_BrNO3
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! Cloud droplet radius in continental warm clouds [cm]
      REAL(fp), PARAMETER :: xCLDR_CONT =  6.e-4_fp

      ! Cloud droplet radius in marine warm clouds [cm]
      REAL(fp), PARAMETER :: xCLDR_MARI = 10.e-4_fp

      ! Ice cloud droplet radius [cm]
      !REAL(fp), PARAMETER :: XCLDrIce = 75.e-4_fp
      REAL(fp), PARAMETER :: XCLDrIce = 38.5e-4_fp

      ! Density of H2O liquid [kg/cm3]
      REAL(fp), PARAMETER :: dens_liq = 0.001e+0_fp

      ! Density of H2O ice [kg/cm3]
      REAL(fp), PARAMETER :: dens_ice = 0.91e-3_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)            :: nu         ! Mean molecular speed
      REAL(fp)            :: RADIUS     ! Radius of cloud droplet      [cm]
      REAL(fp)            :: SQM        ! Square root of molec. weight [g/mol]
      REAL(fp)            :: STK        ! Square root of temperature   [K]
      REAL(fp)            :: DFKG       ! Gas diffusion coefficient    [cm2/s]
      REAL(fp)            :: AREA_L     ! Surface area (liquid)        [cm2/cm3]
      REAL(fp)            :: AREA_I     ! Surface area (ice) )         [cm2/cm3]
      REAL(fp)            :: Vcl, Vci   ! Volume of the cloud (liq and ice) [cm3]

      LOGICAL             :: IS_LAND, IS_ICE, Is_Warm
      REAL(fp)            :: alpha,   beta

      ! Pointers
      REAL(fp), POINTER   :: AD(:,:,:)
      REAL(fp), POINTER   :: CLDF(:,:,:)
      REAL(fp), POINTER   :: FRLAND(:,:)
      REAL(fp), POINTER   :: FROCEAN(:,:)

      !=================================================================
      ! CLD_PARAMS begins here!
      !=================================================================

      ! Initialize pointers
      AD      => State_Met%AD
      CLDF    => State_Met%CLDF
      FRLAND  => State_Met%FRLAND
      FROCEAN => State_Met%FROCEAN

      CLDFr = CLDF(I,J,L)
      CLDFr = Max(CLDFr, 0.0e+0_fp)
      CLDFr = Min(CLDFr, 1.0e+0_fp)

      ! Quick test - is there any cloud?
      IF ( ( (QL+QI) <= 0.0e+0_fp) .or. (CLDF(I,J,L) <= 0.0e+0_fp) ) THEN
         rLiq = xCldR_Cont
         rIce = xCLDrIce
         ALiq = 0.0e+0_fp
         VLiq = 0.0e+0_fp
         AIce = 0.0e+0_fp
         VIce = 0.0e+0_fp
         Return
      ENDIF

      ! ----------------------------------------------
      ! In GC 12.0 and earlier, the liquid water volume was
      ! set to zero at temperatures colder than 258K and
      ! over land ice (Antarctica & Greenland). That
      ! was likely legacy code from GEOS-4, which provided
      ! no information on cloud phase. As of GC 12.0,
      ! all met data sources provide explicit liquid and
      ! ice condensate amounts, so we use those as provided.
      ! (C.D. Holmes)
      ! ----------------------------------------------

      !-----------------------------------------------
      ! Liquid water clouds
      !
      ! Droplets are spheres, so
      ! Surface area = 3 * Volume / Radius
      !
      ! Surface area density = Surface area / Grid volume
      !-----------------------------------------------

      IF ( FRLAND(I,J) > FROCEAN(I,J) ) THEN
         ! Continental cloud droplet radius [cm]
         rLiq = XCLDR_CONT
      ELSE
         ! Marine cloud droplet radius [cm]
         rLiq = XCLDR_MARI
      ENDIF

      ! get the volume of cloud condensate [cm3(condensate)/cm3(air)]
      ! QL is [g/g]
      VLiq = QL * AD(I,J,L) / dens_liq / Vair
      VIce = QI * AD(I,J,L) / dens_ice / Vair

      ALiq = 3.e+0_fp * VLiq / rLiq

      !-----------------------------------------------
      ! Ice water clouds
      !
      ! Surface area calculation requires information about
      ! ice crystal size and shape, which is a function of
      ! temperature. Use Heymsfield (2014) empirical relationships
      ! between temperature, effective radius, surface area
      ! and ice water content.
      !
      ! Schmitt and Heymsfield (2005) found that ice surface area
      ! is about 9 times its cross-sectional area.
      ! For any shape,
      ! Cross section area = pi * (Effective Radius)^2, so
      ! Cross section area = 3 * Volume / ( 4 * Effective Radius ).
      !
      ! Thus, for ice
      ! Surface area = 9 * Cross section area
      !              = 2.25 * 3 * Volume / Effective Radius
      ! (C.D. Holmes)
      !-----------------------------------------------

      ! Heymsfield (2014) ice size parameters
      if (T < 273e+0_fp - 71e+0_fp ) then
          alpha = 83.3e+0_fp
          beta  = 0.0184e+0_fp
      elseif ( T < 273e+0_fp - 56e+0_fp ) then
          alpha = 9.1744e+4_fp
          beta = 0.117e+0_fp
      else
          alpha = 308.4e+0_fp
          beta  = 0.0152e+0_fp
      endif

      ! Effective radius, cm
      rIce = 0.5e+0_fp * alpha * exp( beta * (T-273.15e+0_fp) ) / 1e+4_fp

      ! Ice surface area density, cm2/cm3
      AIce = 3.e+0_fp * VIce / rIce * 2.25e+0_fp

      ! Free Pointers
      NULLIFY( AD      )
      NULLIFY( CLDF    )
      NULLIFY( FRLAND  )
      NULLIFY( FROCEAN )

    END SUBROUTINE Cld_Params
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
      SUBROUTINE GET_HALIDE_CLDCONC( HBr, HCl, VLiq, VIce, VAir, TK, CF, pHCloud, br_conc, cl_conc )
!
! !USES:
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: HCl, HBr  ! Number density [#/cm3]
      REAL(fp),  INTENT(IN) :: VAir    ! Volume of air [cm3]
      REAL(fp),  INTENT(IN) :: VLiq, VIce ! Volume of the cloud (liq and ice) [cm3(condensate)/cm3(air)]
      REAL(fp),  INTENT(IN) :: TK      ! Air temperature [K]
      REAL(fp),  INTENT(IN) :: CF      ! Cloud fraction
      REAL(fp),  INTENT(IN) :: pHCloud ! Cloud pH
!
! !RETURN VALUE:
!
      REAL(fp), INTENT(OUT) :: cl_conc, br_conc ! Liq. phase molar concentration [mol/kg-water]
!
! !REMARKS:
!
! !REVISION HISTORY:
!  21 Dec 2016 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)            :: n_br, n_cl ! dissolved bromide and chloride [#/cm3(air)]
      REAL(fp)            :: V_tot, dr_ratio, t2l !
      REAL(fp)            :: L2G, F_L, pH

      !=================================================================
      ! GET_HALIDE_CLDCONC_G begins here!
      !=================================================================

      !---------------------------------------------------------------
      ! jas, 07/30/2014 (SETUP d/r ratio for ice cloud droplets)
      ! V_liq = 4pi/3 ( r^3 - (r - r*(d/r))^3 = (r^3 - r^3*(1 - d/r)^3) = r^3 (1
      ! - (1 - d/r)^3
      ! V_tot / V_liq = 1 / (1 - (1 - d/r)^3))
      DR_RATIO = 2e-2_fp
      T2L = 1.0e0_fp / ( 1.0e0_fp - (1.0e0_fp - DR_RATIO)**3.0e0_fp )
      !---------------------------------------------------------------


      !V_tot = VLiq + (VIce / T2L) ! (cm3(liq)/cm3(air)
      V_tot = VLiq ! (cm3(liq)/cm3(air)
      V_tot = SAFE_DIV(V_tot, CF, 0.e+0_fp) ! only consider in cloud


      IF (V_tot.lt.1.0e-20) THEN
         br_conc = 1.0e-20_fp
         cl_conc = 1.0e-20_fp
         Return
      ENDIF


      ! Chloride
      CALL COMPUTE_L2G_LOCAL( 1.0e0_fp, 9000.0e0_fp, -6.3e+0_fp, TK, V_tot, L2G, pH)
      F_L = L2G/(1.0e0_fp + L2G)
      cl_conc = F_L * HCl / (V_tot * AVO * 1.0e-3_fp) ! [Cl-] in (mol/L)
      !cl_conc = min(cl_conc,5.0e0_fp)
      !cl_conc = max(cl_conc,1.0e-20_fp)

      ! Bromide
      CALL COMPUTE_L2G_LOCAL( 7.5e-1_fp, 10200.0e0_fp, -9.0e+0_fp, TK, V_tot, L2G, pH)
      F_L = L2G/(1.0e0_fp + L2G)
      br_conc = F_L * HBr / (V_tot * AVO * 1.0e-3_fp) ! [Br-] in (mol/L)

      !br_conc = min(br_conc,5.0e0_fp)
      !br_conc = max(br_conc,1.0e-20_fp)

      END SUBROUTINE GET_HALIDE_CLDCONC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Theta_Ice
!
! !DESCRIPTION: Subroutine GET_THETA_ICE returns theta values for HNO3, HCl, and
! HBr for ice uptake calculations
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_THETA_ICE( HNO3, HCl, HBr, TK, hno3_th, hcl_th, hbr_th)


! !USES:
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: HNO3, HCl, HBr  ! Number density [#/cm3]
      REAL(fp),  INTENT(IN) :: TK      ! Air temperature [K]
!
! !RETURN VALUE:
!
      REAL(fp), INTENT(OUT) :: hno3_th, hbr_th, hcl_th

! !REMARKS:
!
! !REVISION HISTORY:
!  02 Jul 2019 - X. Wang       - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)            :: Nmax, KLangC1, KLangC2, KlinC

      !=================================================================
      ! GET_THETA_ICE begins here!
      !=================================================================

      !HNO3
      Nmax = 2.7e14_fp !molecule cm-2
      KlinC = 7.5e-5_fp * exp(4585_fp/TK)  !cm-1
      KLangC1 = KlinC / Nmax !molecule-1 cm3

      !HCl
      Nmax = 3.0e14_fp !molecule cm-2
      KlinC = 1.3e-5_fp * exp(4600_fp/TK)  !cm-1
      KLangC2 = KlinC / Nmax !molecule-1 cm3

      hno3_th = KLangC1*HNO3 / (1.0e0_fp + KLangC1*HNO3 + KLangC2*HCl)
      hcl_th  = KLangC2*HCl / (1.0e0_fp + KLangC1*HNO3 + KLangC2*HCl)

      !HBr
      hbr_th = 4.14e-10 * (HBr**0.88)
      hbr_th = min(hbr_th, 1.0e0_fp)

      END SUBROUTINE GET_THETA_ICE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get_Halide_SSAConc
!
! !DESCRIPTION: Function GET\_HALIDE\_SSACONC calculates concentration of a
!               halide in sea salt aerosol.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GET_HALIDE_SSACONC( n_x, surf_area, r_w, conc_x )
!
! !DESCRIPTION: Function GET\_HALIDE\_SSACONC calculates concentration of a
!               halide in sea salt aerosol.
!
! !OUTPUT PARAMETER:
      ! concentration of X- in SALX (mol/L)
      REAL(fp)                         :: conc_x
! !INPUT PARAMETERS:
      ! n_x = X-(ssa) number density (#/cm3), surf_area = AERO surface area
      ! conc (cm2/cm3), r_w = AERO wet radius (cm)
      REAL(fp), INTENT(IN)             :: n_x, surf_area, r_w

!
! !REVISION HISTORY:
!  25 Jul 2014 - J. Schmidt - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
!
      !REAL(fp),  PARAMETER :: con_NA      = 6.0221413e23_fp ! #/mol
      REAL(fp)             :: V_tot

      V_tot = surf_area * r_w * 0.3333333e0_fp * 1e-3_fp ! L(liq)/cm3(air)
      IF (V_tot .le. 1.0e-20) THEN
         conc_x = 1.0e-20_fp
         Return
      ELSE
         ! This calculation can be used for both SSA X- concentration and for
         ! those out of cloud only. For X- out of cloud only, V_tot =
         ! V_tot*(1-CF), n_x = n_x*(1-CF), so (1-CF) is canceled.
         ! xnw, 02/05/18
         conc_x =  (n_x / AVO) / V_tot ! mol/L
         !conc_x = MIN(conc_x,5.0e0_fp)
         conc_x = MAX(conc_x,1.0e-20_fp)
      ENDIF

      END SUBROUTINE GET_HALIDE_SSACONC
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
      SUBROUTINE COMPUTE_L2G_LOCAL( K0, CR, pKa, TK, H2OLIQ, L2G, pH )
!
! !USES:
!
      USE Henry_Mod, ONLY : Calc_KH
      USE Henry_Mod, ONLY : Calc_Heff
!
! !INPUT PARAMETERS:
!
      REAL(fp), INTENT(IN)  :: K0     ! Henry's solubility constant [M/atm]
      REAL(fp), INTENT(IN)  :: CR     ! Henry's volatility constant [K]
      REAL(fp), INTENT(IN)  :: pKa    ! Henry's pH correction factor [1]
      REAL(fp), INTENT(IN)  :: TK     ! Temperature [K]
      REAL(fp), INTENT(IN)  :: H2OLIQ ! Liquid water content [cm3 H2O/cm3 air]
      REAL(fp), INTENT(IN)  :: pH     ! Liquid water pH
!
! !OUTPUT PARAMETERS:
!
      REAL(fp), INTENT(OUT) :: L2G    ! Cliq/Cgas ratio [1]
!
! !REMARKS:
!  The ratio Cliq / Cgas is obtained via Henry's law.  The appropriate
!  values of Kstar298 and H298_R must be supplied for each species.
!  (cf Jacob et al 2000, p. 3)
!
! !REVISION HISTORY:
!  23 Feb 2000 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
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
      !pH = 4.5_fp

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

      END SUBROUTINE COMPUTE_L2G_LOCAL
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
      REAL(fp),       Intent(IN) :: DENAIR   ! Density of air [#/cm3]
      REAL(fp),       Intent(In) :: TK       ! Air temperature [K]
      REAL(fp),       Intent(In) :: rLiq     ! Radius of liquid cloud drops [cm]
      REAL(fp),       Intent(In) :: rIce     ! Radius of ice cloud crystals [cm]
      REAL(fp),       Intent(In) :: ALiq     ! Surface area (liquid) [cm2/cm3]
      REAL(fp),       Intent(In) :: AIce     ! Surface area (ice) ) [cm2/cm3]
      REAL(fp),       Intent(IN) :: MX_gmol  ! Molecular mass of XNO3 [g/mol]
      REAL(fp),       Intent(IN) :: AlphaX   ! XNO3 accomodation coef [unitless]
!
! !RETURN VALUE:
!
      REAL(fp)                   :: cld1k    ! Rate constant for heterogeneous
                                             ! cycling of BrNO3 off of cloud
                                             ! particles
!
! !REMARKS:
!
! !REVISION HISTORY:
!  27 Feb 2011 - J. Parrella - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)             :: SQM        ! Square root of molec. weight [g/mol]
      REAL(fp)             :: STK        ! Square root of temperature   [K]
      REAL(fp)             :: DFKG       ! Gas diffusion coefficient    [cm2/s]

      !=================================================================
      ! CLD1K_XNO3 begins here!
      !=================================================================

      ! Quick test - is there any cloud?
      IF ((ALiq.le.0.0e+0_fp).and.(AIce.le.0.0e+0_fp)) THEN
         cld1k = 0.0e+0_fp
         Return
      ENDIF

      ! ------------------------------------------------------------
      !   Calculate the 1st order rate constant for XNO3 hydrolysis.
      !
      !   (a) calculate the gas phase diffusion coefficient;
      !
      !   (b) calculate the hydrolysis rxn rate.
      ! ------------------------------------------------------------
      SQM = sqrt(MX_gmol)    ! square root of molar mass [g/mole]
      STK = sqrt(TK) ! square root of temperature [K]

      ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
      DFKG  = 9.45E+17_fp/DENAIR * STK * SQRT(3.472E-2_fp + 1.E+0_fp/(SQM*SQM))

      ! Compute ARSL1K according to the formula listed above
      ! Sum contribution from ice and liquid clouds
      cld1k = ALiq / ( rLiq/DFKG + 2.749064E-4 * SQM/(ALPHAX*STK) )
      cld1k = AIce / ( rIce/DFKG + 2.749064E-4 * SQM/(ALPHAX*STK) ) + cld1k

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
      REAL(fp), INTENT(IN) :: XCARBN
!
! !RETURN VALUE:
!
      REAL(fp)             :: FC_RO2HO2
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

      FC_RO2HO2 = 1E+0_fp - EXP( -0.245E+0_fp * XCARBN )

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
      REAL(fp), INTENT(IN) :: ZDNUM   ! Air density   [molec/cm3 ]
      REAL(fp), INTENT(IN) :: TT      ! Temperature   [K         ]
!
! !RETURN VALUE:
!
      REAL(fp)             :: FY_HORO
!
! !REVISION HISTORY:
!  (1 ) Branching ratio calculation (tmf, 2/6/05).
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: K1, K2, O2DNUM

      !=================================================================
      ! FYHORO begins here!
      !=================================================================
      O2DNUM  = ZDNUM * 0.21E+0_fp
      K1      = 6.0E-14_fp * EXP(-550.E+0_fp/TT) * O2DNUM
      K2      = 9.5E+13_fp * EXP(-5988.E+0_fp/TT)

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
      REAL(fp), INTENT(IN) :: XCARBN   ! Number of C atoms in RO2
      REAL(fp), INTENT(IN) :: ZDNUM    ! Air density   [molec/cm3 ]
      REAL(fp), INTENT(IN) :: TT       ! Temperature   [K         ]
!
! !RETURN VALUE:
!
      REAL(fp)             :: FYR_NO3
!
! !REVISION HISTORY:
!  (1 ) Original code from Larry Horowitz, Jinyou Liang, Gerry Gardner,
!        and Daniel Jacob circa 1989/1990.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: YYYN, XXYN,  AAA,  RARB, ZZYN
      REAL(fp) :: XF,   ALPHA, Y300, BETA, XMINF, XM0

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
! !ROUTINE: Arsl1k
!
! !DESCRIPTION: Function ARSL1K calculates the 1st-order loss rate of species
!  on wet aerosol surface.
!\\
!\\
! !INTERFACE:
!
    FUNCTION ARSL1K( AREA, RADIUS, DENAIR, STKCF, STK, SQM ) &
         RESULT( ARS_L1K )
!
! !INPUT PARAMETERS:
!
      ! Surface  area of wet aerosols/volume of air [cm2/cm3]
      REAL(fp), INTENT(IN) :: AREA

      ! Radius of wet aerosol [cm], order of 0.01-10 um;
      ! Note that radius here is Rd, not Ro
      REAL(fp), INTENT(IN) :: RADIUS

      ! Density of air [#/cm3]
      REAL(fp), INTENT(IN) :: DENAIR

      ! Sticking coefficient [unitless], order of 0.1
      REAL(fp), INTENT(IN) :: STKCF

      ! Square root of temperature [K]
      REAL(fp), INTENT(IN) :: STK

      ! Square root of molecular weight [g/mole]
      REAL(fp), INTENT(IN) :: SQM
!
! !RETURN VALUE:
!
      REAL(fp)             :: ARS_L1K
!
! !REMARKS:
!  The 1st-order loss rate on wet aerosol (Dentener's Thesis, p. 14)
!  is computed as:
!                                                                             .
!      ARSL1K [1/s] = area / [ radius/dfkg + 4./(stkcf * xmms) ]
!                                                                             .
!  where XMMS = Mean molecular speed [cm/s] = sqrt(8R*TK/pi/M) for Maxwell
!        DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)

! !REVISION HISTORY:
!  01 Jul 1994 - lwh, jyl, gmg, djj - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp) :: DFKG

      !=================================================================
      ! ARSL1K begins here!
      !=================================================================
      IF ( AREA   < 0e0_fp    .or. DENAIR < 1e-30_fp  .or.  &
           RADIUS < 1e-30_fp  .or. SQM    < 1e-30_fp  .or.  &
           STK    < 1e-30_fp  .or. STKCF  < 1e-30_fp ) THEN

         ! Use default value if any of the above values are zero
         ! This will prevent div-by-zero errors in the eqns below
         ! Value changed from 1d-3 to 1d-30 (bhh, jmao, eam, 7/18/2011)
         ARS_L1K = 1.E-30_fp

      ELSE

         ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
         DFKG  = 9.45E+17_fp/DENAIR * STK * SQRT(3.472E-2_fp + 1.E0_fp/ &
          (SQM*SQM))

         ! Compute ARSL1K according to the formula listed above
         ARS_L1K = AREA / ( RADIUS/DFKG + 2.749064E-4_fp*SQM/(STKCF*STK) )

      ENDIF

    END FUNCTION ARSL1K
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Check_Nat
!
! !DESCRIPTION: Subroutine CHECK\_NAT determines whether the solid PSC is
!  composed of ice or NAT (nitric acid trihydrate) (needed for heterogeneous
!  chemistry), or indeed if there is any direct PSC calculation at all. This
!  is important for determining whether to use the JPP or Kirner scheme for
!  ice cloud radii.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_NAT( I, J, L, IS_NAT, IS_PSC, IS_STRAT, &
                            Input_Opt, State_Met, State_Chm )
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)  :: I,J,L      ! Grid indices
      TYPE(OptInput), INTENT(IN)  :: Input_Opt  ! Input options
      TYPE(MetState), INTENT(IN)  :: State_Met  ! Meteorology State object
      TYPE(ChmState), INTENT(IN)  :: State_Chm  ! Chemistry State object
!
! !OUTPUT VARIABLES:
!
      LOGICAL,        INTENT(OUT) :: IS_NAT     ! Is surface NAT?
      LOGICAL,        INTENT(OUT) :: IS_PSC     ! Are there solid PSCs?
      LOGICAL,        INTENT(OUT) :: IS_STRAT   ! Are we in the strat?
!
! !REMARKS:
!  This routine is only activated for UCX-based mechanisms
!
! !REVISION HISTORY:
!  17 Apr 2013 - S. D. Eastham - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
!LOCAL VARIABLES:
!
      LOGICAL :: IS_TROP

      !=================================================================
      ! CHECK_NAT begins here!
      !=================================================================

      ! Check if box is in the troposphere
      IS_TROP  = ( State_Met%InTroposphere(I,J,L) )

      ! Check if box is in the stratosphere
      IS_STRAT = ( State_Met%InStratosphere(I,J,L) .and. ( .not. IS_TROP) )

      ! Check if there are solid PSCs
      IS_PSC   = ( ( Input_Opt%LPSCCHEM ) .and. &
                 ( State_Chm%STATE_PSC(I,J,L) >= 2.0 ) .and. ( IS_STRAT ) )

      ! Check if there is surface NAT
      IS_NAT   = ( ( IS_PSC ) .and. ( C(ind_NIT) .gt. TINY(1e+0_fp) ) )

    END SUBROUTINE CHECK_NAT
!EOC
  END MODULE GCKPP_HETRATES
