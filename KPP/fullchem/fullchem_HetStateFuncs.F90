!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fullchem_HetStateFuncs.F90
!
! !DESCRIPTION: Module for initializing the HetState object, which passes
!  arguments from GEOS-Chem to the heterogeneous chemistry routines.
!\\
!\\
! !INTERFACE:

MODULE fullchem_HetStateFuncs
!
! !USES:
!
  USE GcKpp_Precision

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: fullchem_SetStateHet
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
! !IROUTINE: FullChem_SetStateHet
!
! !DESCRIPTION: Initializes the State_Het object with gridbox values passed
!  from fullchem_mod.  These values are used in the heterogenous chemistry
!  reaction rate computations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE FullChem_SetStateHet( I,         J,         L,                  &
                                   Input_Opt, State_Chm, State_Met,          &
                                   H,         RC                            )
!
! !USES:
!
    USE ErrCode_Mod
    USE GcKpp_Global
    USE GcKpp_Parameters
    USE PhysConstants,    ONLY : AVO, PI
    USE Input_Opt_Mod,    ONLY : OptInput
    USE rateLawUtilFuncs
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Met_Mod,    ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I
    INTEGER,        INTENT(IN)    :: J
    INTEGER,        INTENT(IN)    :: L
    TYPE(OptInput), INTENT(IN)    :: Input_Opt
    TYPE(ChmState), INTENT(IN)    :: State_Chm
    TYPE(MetState), INTENT(IN)    :: State_Met
!
! INPUT/OUTPUT PARAMETERS:
!
    TYPE(HetState), INTENT(INOUT) :: H
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: NA

    !========================================================================
    ! Populate fields of the HetState object in gckpp_Global
    !========================================================================

    ! Initialization
    RC = GC_SUCCESS
    NA = State_Chm%nAeroType

    !========================================================================
    ! Populate fields of the HetState object in gckpp_Global
    !========================================================================

    ! Identify a box for debug printout within rate-law functions
    H%debugBox      = .FALSE.

    ! Constants (so that we can use these within KPP)
    H%AVO           = AVO
    H%PI            = PI

    ! Meteorology-related quantities
    H%CldFr         = MIN(MAX(State_Met%CLDF(I,J,L), 0.0_dp), 1.0_dp)
    H%ClearFr       = 1.0_dp - H%CldFr
    H%QICE          = State_Met%QI(I,J,L)
    H%QLIQ          = State_Met%QL(I,J,L)
    H%vAir          = State_Met%AIRVOL(I,J,L) * 1.0e6_dp

    ! Aerosol fields
    H%nAeroType     = State_Chm%nAeroType
    H%aClArea       = State_Chm%aClArea(I,J,L)
    H%aClRadi       = State_Chm%aClRadi(I,J,L)
    H%aClVol        = H%aClArea * H%aClRadi / 3.0_dp
    H%AWATER(:)     = State_Chm%IsorropAeroH2O(I,J,L,:)
    H%xArea(1:NA)   = State_Chm%AeroArea(I,J,L,1:NA)
    H%xRadi(1:NA)   = State_Chm%AeroRadi(I,J,L,1:NA)
    H%xVol(1:NA)    = H%xArea(1:NA) * H%xRadi(1:NA) / 3.0_dp
    H%wetArea(1:NA) = State_Chm%WetAeroArea(I,J,L,1:NA)
    H%xH2O(1:NA)    = State_Chm%AeroH2O(I,J,L,1:NA) * 1.0e-6_dp
    H%OMOC_POA      = State_Chm%OMOC_POA(I,J)
    H%OMOC_OPOA     = State_Chm%OMOC_OPOA(I,J)

    ! HSO3 and SO3 concentrations in cloud [mol/L]
    H%HSO3_aq       = State_Chm%HSO3_aq(I,J,L)
    H%SO3_aq        = State_Chm%SO3_aq(I,J,L)
    H%TSO3_aq       = H%HSO3_aq + H%SO3_aq
    H%frac_HSO3_aq  = SafeDiv( H%HSO3_aq, H%TSO3_aq, 0.0_dp )
    H%frac_SO3_aq   = SafeDiv( H%SO3_aq,  H%TSO3_aq, 0.0_dp )

    ! Concentrations from ISORROPIA
    H%HSO4_molal    = State_Chm%IsorropBisulfate(I,J,L)
    H%NO3_molal     = State_Chm%IsorropNitrate(I,J,L,1)
    H%SO4_molal     = State_Chm%IsorropSulfate(I,J,L)

    ! pH and alkalinity fields
    H%H_plus        = State_Chm%IsorropHplus(I,J,L,1)
    H%pHCloud       = State_Chm%pHCloud(I,J,L)
    H%pHSSA(:)      = State_Chm%IsorropAeropH(I,J,L,:)
    H%H_conc_Sul    = 10.0**( -1.0_dp * H%pHSSA(1) )
    H%H_conc_LCl    = 10.0**( -1.0_dp * H%pHCloud  )
    H%H_conc_ICl    = 10.0**( -4.5_dp              )
    H%H_conc_SSA    = H%H_conc_Sul
    H%H_conc_SSC    = 10.0**( -5.0_dp              )
    H%ssAlk         = State_Chm%SSAlk(I,J,L,:)
    H%SSA_is_Alk    = ( H%ssAlk(1) > 0.05_dp       )
    H%SSA_is_Acid   = ( .not.  H%SSA_is_Alk        )
    H%SSC_is_Alk    = ( H%ssAlk(2) > 0.05_dp       )
    H%SSC_is_Acid   = ( .not.  H%SSC_is_Alk        )

    ! Other fields
    H%gamma_HO2     = Input_Opt%gamma_HO2

    ! Correction factors for HOBr and HOCl removal by SO2 [1]
    H%fupdateHOBr  = State_Chm%fupdateHOBr(I,J,L)
    H%fupdateHOCl  = State_Chm%fupdateHOCl(I,J,L)

    ! Aqueous S(IV) in cloudwater
    !
    ! -- This is the ratio of HSO3-/SO2, both in units of molec/cm3.
    !    It allows the use of SO2 in the reactions with HOCl and HOBr,
    !    and converts SO2 to HSO3- via the reaction rate constant.
    H%HSO3m = SafeDiv( State_Chm%HSO3_aq(I,J,L) * 1.0e-3_dp *                &
                       State_Het%AVO            *                            &
                       State_Met%QL(I,J,L)      *                            &
                       State_Met%AIRDEN(I,J,L)  * 1.0e-3_dp,                 &
                       State_Met%CLDF(I,J,L),                                &
                       0.0_dp                                               )

    ! Avoid div-by-zero condition
    H%HSO3m = SafeDiv( H%HSO3m, C(ind_SO2), 0.0_dp                          )

    ! -- This is the ratio of SO3--/SO2, both in units of molec/cm3.
    !    It allows the use of SO2 in the reactions with HOCl and HOBr,
    !    and converts SO2 to SO3-- via the reaction rate constant.
    H%SO3mm = SafeDiv( State_Chm%SO3_aq(I,J,L)  * 1.0e-3_dp *                &
                       State_Het%AVO            *                            &
                       State_Met%QL(I,J,L)      *                            &
                       State_Met%AIRDEN(I,J,L)  * 1.0e-3_dp,                 &
                       State_Met%CLDF(I,J,L),                                &
                       0.0_dp                                               )

    ! Avoid div-by-zero condition
    H%SO3mm = SafeDiv( H%SO3mm, C(ind_SO2), 0.0_dp                          )

    ! Cloud fields
    CALL Cld_Params( AD      = State_Met%AD(I,J,L),                          &
                     CLDF    = State_Met%CLDF(I,J,L),                        &
                     FRLAND  = State_Met%FRLAND(I,J),                        &
                     FROCEAN = State_Met%FROCEAN(I,J),                       &
                     QI      = State_Met%QI(I,J,L),                          &
                     QL      = State_Met%QL(I,J,L),                          &
                     T       = State_Met%T(I,J,L),                           &
                     H       = H                                            )

    ! Get theta for ice cloud uptake
    CALL Get_Theta_Ice( C(ind_HNO3), C(ind_HCl), C(ind_HBr), H )

    ! Halide (Br- and Cl-) concentrations
    CALL Halide_Conc( I, J, L, H  )

    !========================================================================
    ! Copy quantities for UCX into gckpp_Global variables
    !========================================================================

    ! ... copy uptake probabilities for PSC reactions on SLA
    ! ... to the proper gckpp_Global variable
    H%KHETI_SLA(1:11) = State_Chm%KHETI_SLA(I,J,L,1:11)

    ! ... check if we are in the stratosphere
    H%stratBox = State_Met%InStratosphere(I,J,L)

    ! ... check if there are solid PSCs at this grid box
    H%pscBox  =                                                              &
         ( ( Input_Opt%LPSCCHEM                ) .and.                       &
           ( State_Chm%STATE_PSC(I,J,L) >= 2.0 ) .and. H%stratBox           )

    ! ... check if there is surface NAT at this grid box
    H%natSurface = ( H%pscBox .and. ( C(ind_NIT) > 0.0_dp )                 )

  END SUBROUTINE FullChem_SetStateHet
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
! !USES:
!
    USE Gckpp_Global, ONLY : HetState, TEMP
!
! !INPUT PARAMETERS:
!
    REAL(dp),       INTENT(IN)    :: HNO3  ! HNO3 conc [molec/cm3]
    REAL(dp),       INTENT(IN)    :: HCl   ! HCl  conc [molec/cm3]
    REAL(dp),       INTENT(IN)    :: HBr   ! HBr  conc [molec/cm3]
!
! !INPUT/OUTPUT PARAMETERS:
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
! !USES:
!
    USE Gckpp_Global
    USE GcKpp_Parameters
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
    USE gckpp_Global
    USE rateLawUtilFuncs, ONLY : SafeDiv
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
    REAL(dp) :: V_tot, dr_ratio, t2l, F_L, L2G, pH

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
    V_tot = SafeDiv( V_tot, H%CldFr, 0.0_dp )

    ! Exit if not in cloud
    IF ( V_tot < 1.0e-20_dp ) THEN
       Br_conc = 0.0_dp
       Cl_conc = 0.0_dp
       RETURN
    ENDIF

    ! Note from Viral Shah (06 Dec 2021):
    !   I believe V_tot corresponds to H2OLIQ, which is the cloud liquid
    !   water content. Whereas L2G is H_eff * H2OLIQ. Note that H_eff is
    !   dimensionless in this equation, not in the more commonly used units
    !   of M/atm.

    ! Chloride (mol/L)
    CALL Compute_L2G_Local( K0     =  1.0_dp,     CR = 9000.0_dp,            &
                            pKa    = -6.3_dp,     TK = TEMP,                 &
                            H2OLIQ =  V_tot,      pH = H%pHCloud,            &
                            L2G    =  L2G                                   )
    F_L = L2G / ( 1.0_dp + L2G )
    Cl_conc = F_L * HCl / (V_tot * H%AVO * 1.0e-3_dp)

    ! Bromide (mol/L)
    CALL Compute_L2G_Local( K0     =  7.5e-1_dp,  CR = 10200.0_dp,           &
                            pKa    = -9.0_dp,     TK = TEMP,                 &
                            H2OLIQ =  V_tot,      pH = H%pHCloud,            &
                            L2G    =  L2G                                   )
    F_L = L2G / ( 1.0_dp + L2G )
    Br_conc = F_L * HBr / ( V_tot * H%AVO * 1.0e-3_dp )

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
! !USES:
!
    USE GcKpp_Global,  ONLY : HetState
    USE PhysConstants, ONLY : AVO
!
! !INPUT PARAMETERS:
!
    REAL(dp),       INTENT(IN)  :: n_x        ! Number density     [#/cm3  ]
    REAL(dp),       INTENT(IN)  :: surf_area  ! Surface area       [cm2/cm3]
    REAL(dp),       INTENT(IN)  :: r_w        ! Aerosol wet radius [cm     ]
!
! !OUTPUT PARAMETERS:
!
    REAL(dp),       INTENT(OUT) :: conc_x     ! Halide conc in seasalt [mol/L]
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
  SUBROUTINE Compute_L2G_Local( K0, CR, pKa, TK, H2OLIQ, pH, L2G )
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
END MODULE fullchem_HetStateFuncs
