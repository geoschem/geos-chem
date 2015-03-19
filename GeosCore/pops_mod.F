!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: pops_mod
!
! !DESCRIPTION: Module POPS\_MOD contains variables and routines for the 
!  GEOS-Chem peristent organic pollutants (POPs) simulation. 
!\\
!\\
! !INTERFACE: 
!
      MODULE POPS_MOD
! 
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp, f4, f8)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC :: CHEMPOPS
      PUBLIC :: INIT_POPS
!
! !REMARKS:
!  POPs Tracers
!  ============================================================================
!  (1 ) POPG   : Gaseous POP - total tracer  
!  (2 ) POPPOC : OC-sorbed POP  - total tracer
!  (3 ) POPPBC : BC-sorbed POP  - total tracer
!
! !REVISION HISTORY:
!  20 Sep 2010 - N.E. Selin    - Initial Version
!  04 Jan 2011 - C.L. Friedman - Expansion on initial version
!  21 Aug 2014 - M. Sulprizio  - Removed emissions routines now handled by HEMCO
!  04 Mar 2015 - R. Yantosca   - Use REAL(f4) for pointer args to HCO_GetPtr 
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL(fp), PARAMETER  :: XNUMOL_OH = 6.022e+23_fp / 17e-3_fp  ! molec OH / kg OH
      REAL(fp), PARAMETER  :: CM3PERM3  = 1.e+6_fp
      REAL(fp), PARAMETER  :: SMALLNUM  = 1e-20_fp

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================

      ! Arrays
      !=================================================================
      ! TCOSZ     : Sum of COS(Solar Zenith Angle ) [unitless]
      ! TTDAY     : Total daylight time at location (I,J) [minutes]
      ! ZERO_DVEL : Array with zero dry deposition velocity [cm/s]
      ! COSZM     : Max daily value of COS(S.Z. angle) [unitless]
      !=================================================================
      REAL(fp), ALLOCATABLE :: TCOSZ(:,:)
      REAL(fp), ALLOCATABLE :: TTDAY(:,:)
      REAL(fp), ALLOCATABLE :: ZERO_DVEL(:,:)
      REAL(fp), ALLOCATABLE :: COSZM(:,:)

      ! Pointers to fields in the HEMCO data structure.
      ! These need to be declared REAL(f4), aka REAL*4.
      REAL(f4), POINTER     :: C_OC(:,:,:) => NULL()
      REAL(f4), POINTER     :: C_BC(:,:,:) => NULL()
      REAL(f4), POINTER     :: O3(:,:,:)   => NULL()
      REAL(f4), POINTER     :: OH(:,:,:)   => NULL()

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
! !IROUTINE:  chempops
!
! !DESCRIPTION: Subroutine CHEMPOPS is the driver routine for POPs chemistry 
!  (eck, 9/20/10)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEMPOPS( am_I_Root, Input_Opt,
     &                     State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version based on CHEMMERCURY
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  23 Oct 2013 - R. Yantosca - Now pass objects to GET_GLOBAL_OH routine
!  10 Sep 2014 - M. Sulprizio- Remove calls to GET_GLOBAL_OH and GET_GLOBAL_O3,
!                              these fields are now read in HEMCO
!  10 Sep 2014 - M. Sulprizio- Now get global concentrations of OH, O3, OC, and
!                              BC from HEMCO
!  03 Mar 2015 - C. Keller   - Now do dry-deposition always outside of chempops
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! For fields from Input_Opt
      INTEGER :: N_TRACERS
      LOGICAL :: LGTMM
      LOGICAL :: LNLPBL
      LOGICAL :: LPRT  

      ! For getting concentrations from HEMCO 
      LOGICAL, SAVE      :: FIRST = .TRUE.
      LOGICAL            :: aIR
      CHARACTER(LEN=255) :: LOC = 'CHEMPOPS (pops_mod.F)'

      !=================================================================
      ! CHEMPOPS begins here!
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS
      
      ! Get fields from Input_Opt
      N_TRACERS = Input_Opt%N_TRACERS
      LPRT      = Input_Opt%LPRT
      LGTMM     = Input_Opt%LGTMM
      LNLPBL    = Input_Opt%LNLPBL

      ! If it's a new 6-hr mean, then get the current average 3-D temperature

      !=================================================================
      ! Perform chemistry on POPs tracers
      !=================================================================
      
      ! am I root? 
      aIR = am_I_Root

      ! Get global concentrations of OH, O3, OC, and BC from HEMCO
      IF ( FIRST ) THEN

         CALL HCO_GetPtr( aIR, 'GLOBAL_OC', C_OC, RC )
         IF ( RC /= GIGC_SUCCESS ) 
     &      CALL ERROR_STOP ( 'Cannot get pointer to GLOBAL_OC', LOC )

         CALL HCO_GetPtr( aIR, 'GLOBAL_BC', C_BC, RC )
         IF ( RC /= GIGC_SUCCESS ) 
     &      CALL ERROR_STOP ( 'Cannot get pointer to GLOBAL_BC', LOC )

         CALL HCO_GetPtr( aIR, 'GLOBAL_OH', OH,   RC )
         IF ( RC /= GIGC_SUCCESS ) 
     &      CALL ERROR_STOP ( 'Cannot get pointer to GLOBAL_OH', LOC )

         CALL HCO_GetPtr( aIR, 'GLOBAL_O3', O3,   RC )
         IF ( RC /= GIGC_SUCCESS ) 
     &      CALL ERROR_STOP ( 'Cannot get pointer to GLOBAL_O3', LOC )

         FIRST = .FALSE.

      ENDIF

      ! Compute diurnal scaling for OH
      CALL OHNO3TIME
      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMPOPS: a OHNO3TIME' )

      !-------------------------
      ! GAS AND PARTICLE PHASE chemistry
      !-------------------------
      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMPOPS: b CHEM_GASPART' )
      
      ! Dry deposition occurs with PBL mixing,
      ! pass zero deposition frequency
      CALL CHEM_POPGP( ZERO_DVEL, ZERO_DVEL, ZERO_DVEL,
     &                 am_I_Root, Input_Opt, State_Met,
     &                 State_Chm, RC )
         
      IF ( LPRT ) CALL DEBUG_MSG( 'CHEMPOPS: a CHEM_GASPART' )
   
      END SUBROUTINE CHEMPOPS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  chem_popgp
!
! !DESCRIPTION: Subroutine CHEM\_POPGP is the chemistry subroutine for the
!  oxidation, gas-particle partitioning, and deposition of POPs.
!  (eck, clf, 1/4/2011)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHEM_POPGP ( V_DEP_G,   V_DEP_P_OC, V_DEP_P_BC,
     &                        am_I_Root, Input_Opt,  State_Met,
     &                        State_Chm, RC   )
!
! !USES:
!
      USE CMN_DIAG_MOD
      USE CMN_SIZE_MOD
#if !defined( NO_BPCH )
      USE DIAG_MOD,           ONLY : AD44
#endif
      USE DIAG53_MOD,         ONLY : AD53_PG_OC_NEG
      USE DIAG53_MOD,         ONLY : AD53_PG_BC_NEG
      USE DIAG53_MOD,         ONLY : AD53_PG_OC_POS
      USE DIAG53_MOD,         ONLY : AD53_PG_BC_POS
      USE DIAG53_MOD,         ONLY : ND53, LD53
      USE DIAG53_MOD,         ONLY : AD53_POPG_OH
      USE DIAG53_MOD,         ONLY : AD53_POPP_OC_O3
      USE DIAG53_MOD,         ONLY : AD53_POPP_BC_O3
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE ERROR_MOD,          ONLY : SAFE_DIV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GRID_MOD,           ONLY : GET_AREA_CM2
      USE PBL_MIX_MOD,        ONLY : GET_FRAC_UNDER_PBLTOP
      USE TIME_MOD,           ONLY : GET_TS_CHEM
      USE TRACERID_MOD,       ONLY : IDTPOPG,    IDTPOPPOC,  IDTPOPPBC
!
! !INPUT PARAMETERS: 
!
      ! Dry deposition frequency for gaseous POP [/s]
      REAL(fp), INTENT(IN)            :: V_DEP_G(IIPAR,JJPAR)

      !  Dry deposition frequency for OC-POP [/s]
      REAL(fp), INTENT(IN)            :: V_DEP_P_OC(IIPAR,JJPAR)

      ! Dry deposition frequency for BC-POP [/s]
      REAL(fp), INTENT(IN)            :: V_DEP_P_BC(IIPAR,JJPAR)   

      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
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
!  References:
!  ============================================================================
!  (1 ) For OH rate constant: Brubaker & Hites. 1998. OH reaction kinetics of
!  PAHs and PCDD/Fs. J. Phys. Chem. A. 102:915-921. 
!
! !REVISION HISTORY: 
!  20 Sep 2010 - N.E. Selin  - Initial Version based on CHEM_HG0_HG2
!  29 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  14 Apr 2014 - R. Yantosca - Now prevent div-by-zero in computations below
!  17 Sep 2014 - C. Keller   - Renamed EmisList_GetDataArr to HCO_GetPtr
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! R = universal gas constant for adjusting KOA for temp: 8.3145 [J/mol/K]
      REAL(fp), PARAMETER :: R          = 8.31e+0_fp  

      ! DENS_OCT = density of octanol, needed for partitioning into OC
      ! 820 [kg/m^3]
      REAL(fp), PARAMETER :: DENS_OCT   = 82e+1_fp

      ! DENS_BC = density of BC, needed for partitioning onto BC
      ! 1 [kg/L] or 1000 [kg/m^3] 
      ! From Lohmann and Lammel, Environ. Sci. Technol., 2004, 38:3793-3803.
      REAL(fp), PARAMETER :: DENS_BC    = 1e+3_fp

!
! !LOCAL VARIABLES:
!
      INTEGER           :: I, J, L
      REAL(fp)            :: DTCHEM,       SUM_F
      REAL(fp)            :: KOA_T,        KBC_T
      REAL(fp)            :: KOC_BC_T,     KBC_OC_T
      REAL(fp)            :: TK
      REAL(fp)            :: AREA_CM2
      REAL(fp)            :: F_PBL,        C_O3
      REAL(fp)            :: C_OH,         C_OC_CHEM,   C_BC_CHEM
      REAL(fp)            :: C_OC_CHEM1,   C_BC_CHEM1
      REAL(fp)            :: K_OH,         AIR_VOL,     K_O3
      REAL(fp)            :: K_OX
      REAL(fp)            :: E_KOX_T,       E_KOX_T_P
      REAL(fp)            :: K_DEPG,        K_DEPP_OC,   K_DEPP_BC
      REAL(fp)            :: OLD_POPG,      OLD_POPP_OC, OLD_POPP_BC
      REAL(fp)            :: NEW_POPG,      NEW_POPP_OC, NEW_POPP_BC
      REAL(fp)            :: POPG_BL,       POPP_OC_BL,  POPP_BC_BL
      REAL(fp)            :: POPG_FT,       POPP_OC_FT,  POPP_BC_FT
      REAL(fp)            :: TMP_POPG,      TMP_OX
      REAL(fp)            :: TMP_POPP_OC,   TMP_POPP_BC
      REAL(fp)            :: GROSS_OX,      GROSS_OX_OH, NET_OX
      REAL(fp)            :: GROSS_OX_OC,   GROSS_OX_BC
      REAL(fp)            :: DEP_POPG,      DEP_POPP_OC, DEP_POPP_BC
      REAL(fp)            :: DEP_POPG_DRY,  DEP_POPP_OC_DRY
      REAL(fp)            :: DEP_POPP_BC_DRY
      REAL(fp)            :: DEP_DRY_FLXG,  DEP_DRY_FLXP_OC
      REAL(fp)            :: DEP_DRY_FLXP_BC
      REAL(fp)            :: OLD_POP_T
      REAL(fp)            :: VR_OC_AIR,     VR_BC_AIR
      REAL(fp)            :: VR_OC_BC,      VR_BC_OC
      REAL(fp)            :: F_POP_OC,      F_POP_BC
      REAL(fp)            :: F_POP_G
      REAL(fp)            :: MPOP_OC,       MPOP_BC,     MPOP_G
      REAL(fp)            :: DIFF_G,        DIFF_OC,     DIFF_BC
      REAL(fp)            :: OC_AIR_RATIO,  OC_BC_RATIO, BC_AIR_RATIO
      REAL(fp)            :: BC_OC_RATIO,   SUM_DIFF
      REAL(fp)            :: TMP_OX_P_OC,   TMP_OX_P_BC
      REAL(fp)            :: NET_OX_OC,     NET_OX_BC
                        
      REAL(fp), SAVE      :: POP_XMW, POP_KOA, POP_KBC, POP_K_POPG_OH
      REAL(fp), SAVE      :: POP_K_POPP_O3A, POP_K_POPP_O3B
      REAL(fp), SAVE      :: POP_HSTAR, POP_DEL_H, POP_DEL_Hw

      ! Delta H for POP [kJ/mol]. Delta H is enthalpy of phase transfer
      ! from gas phase to OC. For now we use Delta H for phase transfer 
      ! from the gas phase to the pure liquid state. 
      ! For PHENANTHRENE: 
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 200, Table 6.3), or -74000 [J/mol].
      ! For PYRENE:
      ! this is taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 200, Table 6.3), or -87000 [J/mol].    
      ! For BENZO[a]PYRENE:
      ! this is also taken as the negative of the Delta H for phase transfer
      ! from the pure liquid state to the gas phase (Schwarzenbach,
      ! Gschwend, Imboden, 2003, pg 452, Prob 11.1), or -110,000 [J/mol]
      REAL(fp)            :: DEL_H

      ! KOA_298 for partitioning of gas phase POP to atmospheric OC
      ! KOA_298 = Cpop in octanol/Cpop in atmosphere at 298 K 
      ! For PHENANTHRENE:
      ! log KOA_298 = 7.64, or 4.37*10^7 [unitless]
      ! For PYRENE:
      ! log KOA_298 = 8.86, or 7.24*10^8 [unitless]
      ! For BENZO[a]PYRENE:
      ! log KOA_298 = 11.48, or 3.02*10^11 [unitless]
      ! (Ma et al., J. Chem. Eng. Data, 2010, 55:819-825).
      REAL(fp)            :: KOA_298

      ! KBC_298 for partitioning of gas phase POP to atmospheric BC
      ! KBC_298 = Cpop in black carbon/Cpop in atmosphere at 298 K
      ! For PHENANTHRENE:
      ! log KBC_298 = 10.0, or 1.0*10^10 [unitless]
      ! For PYRENE:
      ! log KBC_298 = 11.0, or 1.0*10^11 [unitless]
      ! For BENZO[a]PYRENE:
      ! log KBC_298 = 13.9, or 7.94*10^13 [unitless]
      ! (Lohmann and Lammel, EST, 2004, 38:3793-3802)
      REAL(fp)            :: KBC_298

      ! K for reaction POPG + OH  [cm3 /molecule /s]
      ! For PHENANTHRENE: 2.70d-11
      ! (Source: Brubaker & Hites, J. Phys Chem A 1998)
      ! For PYRENE: 5.00d-11
      ! Calculated with AOPWIN
      ! For BENZO[a]PYRENE: 5.00d-11
      ! Calculated with AOPWIN 
      REAL(fp)            :: K_POPG_OH !(Gas phase)

      ! k for reaction POPP + O3 [/s] depends on fitting parameters A and B. 
      ! A represents the maximum number of surface sites available to O3, and B 
      ! represents the ratio of desorption/adsorption rate coefficients for
      ! both bulk phases (Ref: Kahan et al Atm Env 2006, 40:3448)
      ! k(obs) = A x [O3(g)] / (B + [O3(g)])
      ! For PHENANTHRENE: A = 0.5 x 10^-3 s^-1, B = 2.15 x 10^15 molec/cm3
      ! For PYRENE: A = 0.7 x 10^-3 s^-1, B = 3 x 10^15 molec/cm3
      ! for BaP: A = 5.5 x 10^-3 s^-1, B = 2.8 x 10^15 molec/cm3
      REAL(fp)            :: AK  ! s^-1
      REAL(fp)            :: BK  ! molec/cm3

      ! On-particle reaction scheme 3: According to Kwamena et al. (J. Phys.
      !  Chem. A 2004 108:11626), reaction will proceed with rate
      !  k = kmax(KO3)[O3]/(1+KO3[O3])
      ! For wet axelaic acid aerosols, kmax = 0.060 s^-1 and
      !  KO3 = 0.028 x 10^-13 cm3
!      REAL(fp), PARAMETER :: KMAX = 0.060 ! s^-1
!      REAL(fp), PARAMETER :: KO3 = 0.028d-13 ! cm^3

      ! K for reaction POPP + NO3 could be added here someday

      ! For fields from Input_Opt
      INTEGER         :: N_TRACERS
      LOGICAL         :: LNLPBL
      LOGICAL         :: LGTMM
      REAL(fp)          :: XNUMOL(Input_Opt%N_TRACERS)

      REAL(fp)          :: DENOM

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER :: STT(:,:,:,:)

      ! Grid box volume
      REAL(fp)          :: BOXVL

      !=================================================================
      ! CHEM_POPGP begins here!
      !=================================================================

      ! Assume success
      RC        = GIGC_SUCCESS

      ! Copy values from Input_Opt
      DEL_H     = Input_Opt%POP_DEL_H
      KOA_298   = Input_Opt%POP_KOA
      KBC_298   = Input_Opt%POP_KBC
      K_POPG_OH = Input_Opt%POP_K_POPG_OH
      AK        = Input_Opt%POP_K_POPP_O3A
      BK        = Input_Opt%POP_K_POPP_O3B
      N_TRACERS = Input_Opt%N_TRACERS
      LNLPBL    = Input_Opt%LNLPBL
      LGTMM     = Input_Opt%LGTMM
      XNUMOL    = Input_Opt%XNUMOL(1:N_TRACERS)

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Chemistry timestep [s]
      DTCHEM = GET_TS_CHEM() * 60e+0_fp

      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Zero concentrations in loop
         MPOP_G          = 0e+0_fp
         MPOP_OC         = 0e+0_fp
         MPOP_BC         = 0e+0_fp
         OLD_POPG        = 0e+0_fp
         OLD_POPP_OC     = 0e+0_fp
         OLD_POPP_BC     = 0e+0_fp
         OLD_POP_T       = 0e+0_fp
         NEW_POPG        = 0e+0_fp
         NEW_POPP_OC     = 0e+0_fp
         NEW_POPP_BC     = 0e+0_fp
         POPG_BL         = 0e+0_fp
         POPP_OC_BL      = 0e+0_fp
         POPP_BC_BL      = 0e+0_fp
         POPG_FT         = 0e+0_fp
         POPP_OC_FT      = 0e+0_fp
         POPP_BC_FT      = 0e+0_fp
         DIFF_G          = 0e+0_fp
         DIFF_OC         = 0e+0_fp
         DIFF_BC         = 0e+0_fp
         NET_OX          = 0e+0_fp 
         TMP_POPG        = 0e+0_fp
         TMP_OX          = 0e+0_fp      
         GROSS_OX        = 0e+0_fp
         GROSS_OX_OH     = 0e+0_fp 
         DEP_POPG        = 0e+0_fp
         DEP_POPP_OC     = 0e+0_fp
         DEP_POPP_BC     = 0e+0_fp
         DEP_POPG_DRY    = 0e+0_fp
         DEP_POPP_OC_DRY = 0e+0_fp
         DEP_POPP_BC_DRY = 0e+0_fp
         DEP_DRY_FLXG    = 0e+0_fp
         DEP_DRY_FLXP_OC = 0e+0_fp
         DEP_DRY_FLXP_BC = 0e+0_fp
         E_KOX_T         = 0e+0_fp
         K_OX            = 0e+0_fp
         K_O3            = 0e+0_fp
         GROSS_OX_OC     = 0e+0_fp
         GROSS_OX_BC     = 0e+0_fp

         ! Save local temperature in TK for convenience [K]
         TK = State_Met%T(I,J,L)

         ! Grid box volume [cm3]
         BOXVL = State_Met%AIRVOL(I,J,L) * 1e+6_fp

         ! Get monthly mean OH concentrations
         C_OH        = GET_OH( I, J, L, State_Met )

         ! Get monthly mean O3 concentrations
         ! O3 is in v/v (from HEMCO), convert to molec/cm3 (mps, 9/10/14)
         C_O3 = O3(I,J,L) * ( 6.022e+23_fp / 28.97e-3_fp ) * 
     &          State_Met%AD(I,J,L)    /  BOXVL

         ! Fraction of box (I,J,L) underneath the PBL top [dimensionless]
         F_PBL = GET_FRAC_UNDER_PBLTOP( I, J, L )

         ! Define K for the oxidation reaction with POPG [/s]
         K_OH        = K_POPG_OH * C_OH

         ! Define K for the oxidation reaction with POPPOC and POPPBC [/s]
         ! (Kahan)
         K_O3        = ( AK * C_O3) / (BK + C_O3) 

         ! Define K for the oxidation reaction with POPPOC and POPPBC [/s]
         ! (Kwamena)
!         K_O3        = ( KMAX * KO3 * C_O3) / (1 + KO3 * C_O3) 
         
         ! Could add K for oxidation by NO3 here one day [/s]

         ! Total K for oxidation [/s]
         K_OX        = K_OH !+ ...

         ! Define Ks for dry deposition of gas phase POP [/s]
         K_DEPG = V_DEP_G(I,J)

         ! Define Ks for dry deposition of particle phase POP [/s]
         K_DEPP_OC = V_DEP_P_OC(I,J)

         ! Define Ks for dry deposition of particle phase POP [/s]
         K_DEPP_BC = V_DEP_P_BC(I,J)

         ! Precompute exponential factors [dimensionless]
         ! For gas phase (OH):
         E_KOX_T  = EXP( -K_OX  * DTCHEM )
         ! For OC and BC phase (O3):
         E_KOX_T_P = EXP(-K_O3  * DTCHEM )

         !==============================================================
         ! GAS-PARTICLE PARTITIONING
         !==============================================================

         OLD_POPG = MAX( STT(I,J,L,IDTPOPG), SMALLNUM )  ![kg]
         OLD_POPP_OC = MAX( STT(I,J,L,IDTPOPPOC), SMALLNUM )  ![kg]
         OLD_POPP_BC = MAX( STT(I,J,L,IDTPOPPBC), SMALLNUM )  ![kg]

         ! Total POPs in box I,J,L 
         OLD_POP_T = OLD_POPG + OLD_POPP_OC + OLD_POPP_BC

         ! Define temperature-dependant partition coefficients:
         KOA_T = KOA_298 * EXP((-DEL_H/R) * ((1e+0_fp/TK) - 
     &                  (1e+0_fp/298e+0_fp))) 

         ! Define KBC_T, the BC-air partition coeff at temp T [unitless]
         ! TURN OFF TEMPERATURE DEPENDENCY FOR SENSITIVITY ANALYSIS
         KBC_T = KBC_298 * EXP((-DEL_H/R) * ((1e+0_fp/TK) - 
     &                  (1e+0_fp/298e+0_fp)))

         ! Define KOC_BC_T, theoretical OC-BC part coeff at temp T [unitless]
         KOC_BC_T = KOA_T / KBC_T

         ! Define KBC_OC_T, theoretical BC_OC part coeff at temp T [unitless]
         KBC_OC_T = 1e+0_fp / KOC_BC_T

         ! Get monthly mean OC and BC concentrations [kg/box] 
         C_OC_CHEM = C_OC(I,J,L)
         C_BC_CHEM = C_BC(I,J,L)

         ! Make sure OC is not negative
         C_OC_CHEM = MAX( C_OC_CHEM, 0e+0_fp )

         ! Convert to units of volume per box [m^3 OC or BC/box]
         C_OC_CHEM1 = C_OC_CHEM / DENS_OCT
         C_BC_CHEM1 = C_BC_CHEM / DENS_BC

         ! Get AIRVOL
         AIR_VOL = State_Met%AIRVOL(I,J,L)

         ! Define volume ratios:
         ! VR_OC_AIR = volume ratio of OC to air [unitless]     
         VR_OC_AIR     = C_OC_CHEM1 / AIR_VOL ! could be zero

         ! VR_OC_BC = volume ratio of OC to BC [unitless]
         VR_OC_BC      = SAFE_DIV( C_OC_CHEM1, C_BC_CHEM1, 0e+0_fp )

         ! VR_BC_AIR = volume ratio of BC to air [unitless]
         VR_BC_AIR     = SAFE_DIV( VR_OC_AIR,  VR_OC_BC,   0e+0_fp )

         ! VR_BC_OC = volume ratio of BC to OC [unitless]
         VR_BC_OC      = SAFE_DIV( 1e+0_fp,        VR_OC_BC,   0e+0_fp )

         ! Redefine fractions of total POPs in box (I,J,L) that are OC-phase, 
         ! BC-phase, and gas phase with new time step (should only change if 
         ! temp changes or OC/BC concentrations change) 
         DENOM         = KOA_T * VR_OC_AIR
         OC_AIR_RATIO  = SAFE_DIV( 1e+0_fp,        DENOM,      0e+0_fp )

         DENOM         = KOC_BC_T * VR_OC_BC
         OC_BC_RATIO   = SAFE_DIV( 1e+0_fp,        DENOM,      0e+0_fp )

         DENOM         = KBC_T * VR_BC_AIR
         BC_AIR_RATIO  = SAFE_DIV( 1e+0_fp,        DENOM,      0e+0_fp )

         DENOM         = KBC_OC_T * VR_BC_OC
         BC_OC_RATIO   = SAFE_DIV( 1e+0_fp,        DENOM,      0e+0_fp )

         ! If there are zeros in OC or BC concentrations, make sure they
         ! don't cause problems with phase fractions
         IF ( C_OC_CHEM > SMALLNUM .and. C_BC_CHEM > SMALLNUM ) THEN
            F_POP_OC  = 1e+0_fp / (1e+0_fp + OC_AIR_RATIO + OC_BC_RATIO)
            F_POP_BC  = 1e+0_fp / (1e+0_fp + BC_AIR_RATIO + BC_OC_RATIO)
         
         ELSE IF (C_OC_CHEM > SMALLNUM .and.
     &            C_BC_CHEM .le. SMALLNUM ) THEN
            F_POP_OC  = 1e+0_fp / (1e+0_fp + OC_AIR_RATIO)
            F_POP_BC  = SMALLNUM           

         ELSE IF ( C_OC_CHEM .le. SMALLNUM .and.
     &             C_BC_CHEM > SMALLNUM ) THEN
            F_POP_OC  = SMALLNUM
            F_POP_BC  = 1e+0_fp / (1e+0_fp + BC_AIR_RATIO)

         ELSE IF ( C_OC_CHEM .le. SMALLNUM .and. 
     &             C_BC_CHEM .le. SMALLNUM) THEN
            F_POP_OC = SMALLNUM
            F_POP_BC = SMALLNUM
         ENDIF

         ! Gas-phase:
         F_POP_G   = 1e+0_fp - F_POP_OC - F_POP_BC

         ! Check that sum equals 1
         SUM_F = F_POP_OC + F_POP_BC + F_POP_G
         
         ! Calculate new masses of POP in each phase [kg]
         ! OC-phase:
         MPOP_OC    = F_POP_OC * OLD_POP_T

         ! BC-phase
         MPOP_BC     = F_POP_BC * OLD_POP_T

         ! Gas-phase
         MPOP_G     = F_POP_G  * OLD_POP_T

         ! Ensure new masses of POP in each phase are positive
         MPOP_OC = MAX(MPOP_OC, SMALLNUM)
         MPOP_BC = MAX(MPOP_BC, SMALLNUM)
         MPOP_G  = MAX(MPOP_G,  SMALLNUM)     

         ! Calculate differences in masses in each phase from previous time
         ! step for storage in ND53 diagnostic

         DIFF_G = MPOP_G - OLD_POPG
         DIFF_OC = MPOP_OC - OLD_POPP_OC
         DIFF_BC = MPOP_BC - OLD_POPP_BC

         ! Sum of differences should equal zero
         SUM_DIFF = DIFF_G + DIFF_OC + DIFF_BC

         !==============================================================
         ! ND53 diagnostic: Differences in distribution of gas and
         ! particle phases between time steps [kg]
         !==============================================================

         IF ( ND53 > 0 .AND. L <= LD53 ) THEN ! LD53 is max level

            IF (DIFF_OC .lt. 0) THEN
 
               AD53_PG_OC_NEG(I,J,L) = AD53_PG_OC_NEG(I,J,L)  + 
     &                                 DIFF_OC

            ELSE IF (DIFF_OC .eq. 0 .or. DIFF_OC .gt. 0) THEN

               AD53_PG_OC_POS(I,J,L) = AD53_PG_OC_POS(I,J,L)  + 
     &                                 DIFF_OC

            ENDIF

            IF (DIFF_BC .lt. 0) THEN

               AD53_PG_BC_NEG(I,J,L) = AD53_PG_BC_NEG(I,J,L)  + 
     &                                 DIFF_BC
               
            ELSE IF (DIFF_BC .eq. 0 .or. DIFF_BC .gt. 0) THEN

               AD53_PG_BC_POS(I,J,L) = AD53_PG_BC_POS(I,J,L)  + 
     &                                 DIFF_BC

            ENDIF

         ENDIF


         !==============================================================
         ! CHEMISTRY AND DEPOSITION REACTIONS
         !==============================================================
         IF ( F_PBL < 0.05e+0_fp .OR. K_DEPG < SMALLNUM ) THEN

            !==============================================================
            ! Entire box is in the free troposphere
            ! or deposition is turned off, so use RXN without deposition
            ! for gas phase POPs
            ! For particle POPs, no rxn and no deposition
            !==============================================================

            CALL RXN_OX_NODEP( MPOP_G, K_OX,
     &                         E_KOX_T, NEW_POPG, GROSS_OX )

            CALL RXN_OX_NODEP( MPOP_OC, K_O3, 
     &                         E_KOX_T_P, NEW_POPP_OC, GROSS_OX_OC)

            CALL RXN_OX_NODEP( MPOP_BC, K_O3, 
     &                         E_KOX_T_P, NEW_POPP_BC, GROSS_OX_BC)

!            NEW_POPP_OC = MPOP_OC
!            NEW_POPP_BC = MPOP_BC

            ! No deposition occurs [kg]
            DEP_POPG = 0e+0_fp
            DEP_POPP_OC = 0e+0_fp
            DEP_POPP_BC = 0e+0_fp
               

         ELSE IF ( F_PBL > 0.95e+0_fp ) THEN 

            !==============================================================
            ! Entire box is in the boundary layer
            ! so use RXN with deposition for gas phase POPs
            ! Deposition only (no rxn) for particle phase POPs
            !==============================================================

            CALL RXN_OX_WITHDEP( MPOP_G,      K_OX,
     &                           K_DEPG,      DTCHEM,      E_KOX_T,
     &                           NEW_POPG,    GROSS_OX,    DEP_POPG    )

            CALL RXN_OX_WITHDEP( MPOP_OC,     K_O3,
     &                           K_DEPP_OC,   DTCHEM,      E_KOX_T_P,
     &                           NEW_POPP_OC, GROSS_OX_OC, DEP_POPP_OC )

            CALL RXN_OX_WITHDEP( MPOP_BC,     K_O3,
     &                           K_DEPP_BC,   DTCHEM,      E_KOX_T_P,
     &                           NEW_POPP_BC, GROSS_OX_BC, DEP_POPP_BC )

!            CALL NO_RXN_WITHDEP( MPOP_OC, K_DEPP_OC, DTCHEM,
!     &                           NEW_POPP_OC, DEP_POPP_OC )

!            CALL NO_RXN_WITHDEP( MPOP_BC, K_DEPP_BC, DTCHEM,
!     &                           NEW_POPP_BC, DEP_POPP_BC )

         ELSE

            !==============================================================
            ! Box spans the top of the boundary layer
            ! Part of the mass is in the boundary layer and subject to 
            ! deposition while part is in the free troposphere and
            ! experiences no deposition.
            !
            ! We apportion the mass between the BL and FT according to the
            ! volume fraction of the box in the boundary layer.
            ! Arguably we should assume uniform mixing ratio, instead of
            ! uniform density but if the boxes are short, the air density
            ! doesn't change much.
            ! But assuming uniform mixing ratio across the inversion layer
            ! is a poor assumption anyway, so we are just using the
            ! simplest approach.
            !==============================================================

            ! Boundary layer portion of POPG [kg]
            POPG_BL = MPOP_G * F_PBL 

            ! Boundary layer portion of POPP_OC [kg]
            POPP_OC_BL = MPOP_OC * F_PBL

            ! Boundary layer portion of POPP_BC [kg]
            POPP_BC_BL = MPOP_BC * F_PBL

            ! Free troposphere portion of POPG [kg]
            POPG_FT = MPOP_G - POPG_BL

            ! Free troposphere portion of POPP_OC [kg]
            POPP_OC_FT = MPOP_OC - POPP_OC_BL

            ! Free troposphere portion of POPP_BC [kg]
            POPP_BC_FT = MPOP_BC - POPP_BC_BL
               
            ! Do chemistry with deposition on BL fraction for gas phase
            CALL RXN_OX_WITHDEP( POPG_BL,  K_OX,
     &                           K_DEPG,   DTCHEM,   E_KOX_T,
     &                           NEW_POPG, GROSS_OX, DEP_POPG )         

            ! Do chemistry without deposition on the FT fraction for gas phase
            CALL RXN_OX_NODEP( POPG_FT, K_OX,
     &                         E_KOX_T, TMP_POPG, TMP_OX ) 

            ! Now do the same with the OC and BC phase:

            ! Do chemistry with deposition on BL fraction for OC phase
            CALL RXN_OX_WITHDEP( POPP_OC_BL,  K_O3,
     &                           K_DEPP_OC,   DTCHEM,      E_KOX_T_P,
     &                           NEW_POPP_OC, GROSS_OX_OC, DEP_POPP_OC )

            ! Do chemistry without deposition on the FT fraction for OC phase
            CALL RXN_OX_NODEP( POPP_OC_FT, K_O3,
     &                         E_KOX_T_P,  TMP_POPP_OC, TMP_OX_P_OC )

            ! Do chemistry with deposition on BL fraction for BC phase
            CALL RXN_OX_WITHDEP( POPP_BC_BL,  K_O3,
     &                           K_DEPP_BC,   DTCHEM,      E_KOX_T_P,
     &                           NEW_POPP_BC, GROSS_OX_BC, DEP_POPP_BC )

            ! Do chemistry without deposition on the FT fraction for BC phase
            CALL RXN_OX_NODEP( POPP_BC_FT, K_O3,
     &                          E_KOX_T_P, TMP_POPP_BC, TMP_OX_P_BC )   

            ! Do deposition (no chemistry) on BL fraction for particulate phase
            ! No deposition (and no chem) on the FT fraction
            ! for the particulate phase
!            CALL NO_RXN_WITHDEP( POPP_OC_BL, K_DEPP_OC, DTCHEM,  
!     &                           NEW_POPP_OC, DEP_POPP_OC)

!            CALL NO_RXN_WITHDEP( POPP_BC_BL, K_DEPP_BC, DTCHEM,  
!     &                          NEW_POPP_BC, DEP_POPP_BC)
               
            ! Recombine the boundary layer and free troposphere parts [kg]
            NEW_POPG    = NEW_POPG + TMP_POPG
            NEW_POPP_OC = NEW_POPP_OC + TMP_POPP_OC
            NEW_POPP_BC = NEW_POPP_BC + TMP_POPP_BC

            ! Total gross oxidation of gas phase in the BL and FT [kg]
            GROSS_OX = GROSS_OX + TMP_OX

            ! Total gross oxidation of particle phases in the BL and FT [kg]
            GROSS_OX_OC = GROSS_OX_OC + TMP_OX_P_OC
            GROSS_OX_BC = GROSS_OX_BC + TMP_OX_P_BC

         ENDIF

         ! Ensure positive concentration [kg]
         NEW_POPG    = MAX( NEW_POPG, SMALLNUM )
         NEW_POPP_OC = MAX( NEW_POPP_OC, SMALLNUM )
         NEW_POPP_BC = MAX( NEW_POPP_BC, SMALLNUM )

         ! Archive new POPG and POPP values [kg]
         STT(I,J,L,IDTPOPG)   = NEW_POPG
         STT(I,J,L,IDTPOPPOC) = NEW_POPP_OC
         STT(I,J,L,IDTPOPPBC) = NEW_POPP_BC

         ! Net oxidation [kg] (equal to gross ox for now)
         NET_OX = MPOP_G - NEW_POPG - DEP_POPG   
         NET_OX_OC = MPOP_OC - NEW_POPP_OC - DEP_POPP_OC
         NET_OX_BC = MPOP_BC - NEW_POPP_BC - DEP_POPP_BC                

         ! Error check on gross oxidation [kg]
         IF ( GROSS_OX < 0e+0_fp .or. GROSS_OX_OC < 0e+0_fp .or.
     &        GROSS_OX_BC < 0e+0_fp ) 
     &        CALL DEBUG_MSG('CHEM_POPGP: negative gross oxidation')

         ! Apportion gross oxidation between OH and possibly
         ! NO3 someday [kg]
         IF ( (K_OX     < SMALLNUM) .OR. (GROSS_OX < SMALLNUM) ) THEN
            GROSS_OX_OH = 0e+0_fp
!            GROSS_OX_NO3 = 0D0

         ELSE
            GROSS_OX_OH = GROSS_OX * K_OH / K_OX
!            GROSS_OX_NO3 = GROSS_OX * K_NO3 / K_OX
         ENDIF

         ! Apportion deposition [kg]
         ! Right now only using dry deposition (no sea salt) (clf, 1/27/11)
         ! If ever use dep with sea salt aerosols,
         ! will need to multiply DEP_POPG by the ratio 
         ! of K_DRYG (rate of dry dep) to K_DEPG (total dep rate).
         IF ( (K_DEPG  < SMALLNUM) .OR. (DEP_POPG < SMALLNUM) ) THEN
            DEP_POPG_DRY  = 0e+0_fp
         ELSE 
            DEP_POPG_DRY  = DEP_POPG   
         ENDIF

         IF ( (K_DEPP_OC   < SMALLNUM)   .OR. 
     &        (DEP_POPP_OC < SMALLNUM) ) THEN
            DEP_POPP_OC_DRY  = 0e+0_fp
         ELSE
            DEP_POPP_OC_DRY  = DEP_POPP_OC
         ENDIF

         IF ( (K_DEPP_BC   < SMALLNUM)   .OR. 
     &        (DEP_POPP_BC < SMALLNUM) ) THEN
            DEP_POPP_BC_DRY  = 0e+0_fp
         ELSE
            DEP_POPP_BC_DRY  = DEP_POPP_BC 
         ENDIF

#if !defined( NO_BPCH )
         !=================================================================
         ! ND44 diagnostic: drydep flux of POPG and POPP [molec/cm2/s]
         !=================================================================
         IF ( ( ND44 > 0 .OR. LGTMM ) .AND. (.NOT. LNLPBL) ) THEN
         ! Not using LGTMM right now (logical switch for using GTMM soil model)
         ! Also not using non-local PBL mode yet (clf, 1/27/2011)

            ! Grid box surface area [cm2]
            AREA_CM2 = GET_AREA_CM2( I, J, L )

            ! Amt of POPG lost to drydep [molec/cm2/s]
            DEP_DRY_FLXG  = DEP_POPG_DRY * XNUMOL(IDTPOPG) / 
     &                    ( AREA_CM2 * DTCHEM )

            ! Archive POPG drydep flux in AD44 array [molec/cm2/s]
            AD44(I,J,IDTPOPG,1) = AD44(I,J,IDTPOPG,1) +
     &                            DEP_DRY_FLXG

            ! Amt of POPPOC lost to drydep [molec/cm2/s]
            DEP_DRY_FLXP_OC = DEP_POPP_OC_DRY * 
     &                        XNUMOL(IDTPOPPOC)/( AREA_CM2 * DTCHEM )        

            ! Archive POPPOC drydep flux in AD44 array [molec/cm2/s]
            AD44(I,J,IDTPOPPOC,1) = AD44(I,J,IDTPOPPOC,1) + 
     &                              DEP_DRY_FLXP_OC

            ! Amt of POPPBC lost to drydep [molec/cm2/s] 
            DEP_DRY_FLXP_BC = DEP_POPP_BC_DRY * 
     &                        XNUMOL(IDTPOPPBC)/( AREA_CM2 * DTCHEM )        

            ! Archive POPPBC drydep flux in AD44 array [molec/cm2/s]
            AD44(I,J,IDTPOPPBC,1) = AD44(I,J,IDTPOPPBC,1) +
     &                              DEP_DRY_FLXP_BC


         ENDIF
#endif

         !==============================================================
         ! ND53 diagnostic: Oxidized POPG (OH-POPG) production [kg]
         !==============================================================

         IF ( ND53 > 0 .AND. L <= LD53 ) THEN ! LD53 is max level

            AD53_POPG_OH(I,J,L)    = AD53_POPG_OH(I,J,L) + GROSS_OX
            AD53_POPP_OC_O3(I,J,L) = AD53_POPP_OC_O3(I,J,L) + 
     &                               GROSS_OX_OC
            AD53_POPP_BC_O3(I,J,L) = AD53_POPP_BC_O3(I,J,L) + 
     &                               GROSS_OX_BC

         ENDIF

      ENDDO
      ENDDO
      ENDDO

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE CHEM_POPGP 
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  rxn_ox_nodep
!
! !DESCRIPTION: Subroutine RXN\_OX\_NODEP calculates new mass of POPG for given
! oxidation rates, without any deposition. This is for the free troposphere, or
! simulations with deposition turned off. (clf, 1/27/11, based on
! RXN\_REDOX\_NODEP in mercury\_mod.f).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RXN_OX_NODEP( OLD_POPG, K_OX, E_KOX_T,
     &                         NEW_POPG, GROSS_OX )
!
! !INPUT PARAMETERS: 
!
      REAL(fp),  INTENT(IN)  :: OLD_POPG
      REAL(fp),  INTENT(IN)  :: K_OX
      REAL(fp),  INTENT(IN)  :: E_KOX_T      
!
! !OUTPUT PARAMETERS:
!
      REAL(fp),  INTENT(OUT) :: NEW_POPG,  GROSS_OX
!
! !REVISION HISTORY: 
!  27 January 2011 - CL Friedman - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! RXN_OX_NODEP begins here!
      !=================================================================

      ! Oxidation
      IF (K_OX < SMALLNUM ) THEN

         GROSS_OX = 0e+0_fp
         NEW_POPG = OLD_POPG

      ELSE 

         ! New concentration of POPG
         NEW_POPG = OLD_POPG * E_KOX_T

         ! Gross oxidation 
         GROSS_OX = OLD_POPG - NEW_POPG
         GROSS_OX = MAX( GROSS_OX, 0e+0_fp )

      ENDIF

      END SUBROUTINE RXN_OX_NODEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!-----------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  rxn_ox_withdep
!
! !DESCRIPTION: Subroutine RXN\_OX\_WITHDEP calculates new mass of POPG for
!  given rates of oxidation and deposition. This is for the boundary layer.
!  (clf, 1/27/11, based on RXN\_REDOX\_NODEP in mercury\_mod.f).
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RXN_OX_WITHDEP( OLD_POPG, K_OX, K_DEPG, DT, E_KOX_T,
     &                           NEW_POPG, GROSS_OX, DEP_POPG )
!
! USES:
!
      USE ERROR_MOD,    ONLY : ERROR_STOP
!
! !INPUT PARAMETERS: 
!
      REAL(fp),  INTENT(IN)  :: OLD_POPG,  DT
      REAL(fp),  INTENT(IN)  :: K_OX, K_DEPG
      REAL(fp),  INTENT(IN)  :: E_KOX_T
!
! !OUTPUT PARAMETERS:
!
      REAL(fp),  INTENT(OUT) :: NEW_POPG,  GROSS_OX
      REAL(fp),  INTENT(OUT) :: DEP_POPG
!
! !REVISION HISTORY: 
!  27 January 2011 - CL Friedman - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)               :: E_KDEPG_T
      REAL(fp)               :: NEWPOPG_OX
      REAL(fp)               :: NEWPOPG_DEP

      !=================================================================
      ! RXN_OX_WITHDEP begins here!
      !=================================================================

      ! Precompute exponential factor for deposition [dimensionless]
      E_KDEPG_T = EXP( -K_DEPG * DT )

      IF (K_OX < SMALLNUM) THEN     
         
         !=================================================================
         ! No Chemistry, Deposition only
         !=================================================================

         ! New mass of POPG [kg]
         NEW_POPG = OLD_POPG * E_KDEPG_T
         
         ! Oxidation of POPG [kg]
         GROSS_OX = 0e+0_fp

         ! Deposited POPG [kg]
         DEP_POPG = OLD_POPG - NEW_POPG

      ELSE

         !=================================================================
         ! Oxidation and Deposition 
         !=================================================================

         ![POPG](t) = [POPG](0) exp( -(kOx + kDPOPG) t)
         !Ox(t)     = ( [POPG](0) - [POPG](t) ) * kOx / ( kOx + kDPOPG )
         !Dep_POPG(t)   = ( [POPG](0) - [POPG](t) - Ox(t) ) 

         ! New concentration of POPG [kg]
         NEW_POPG = OLD_POPG * E_KOX_T * E_KDEPG_T

         ! Gross oxidized gas phase mass [kg]
         GROSS_OX = ( OLD_POPG - NEW_POPG ) * K_OX / ( K_OX + K_DEPG )
         GROSS_OX = MAX( GROSS_OX, 0e+0_fp )

         ! POPG deposition [kg]
         DEP_POPG = ( OLD_POPG - NEW_POPG - GROSS_OX )       
         DEP_POPG = MAX( DEP_POPG, 0e+0_fp )

      ENDIF

      END SUBROUTINE RXN_OX_WITHDEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  no_rxn_withdep
!
! !DESCRIPTION: Subroutine NO\_RXN\_WITHDEP calculates new mass of POPP for
!  given rate of deposition. No oxidation of POPP. This is for the boundary
!  layer. (clf, 2/9/11)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE NO_RXN_WITHDEP( OLD_POPP, K_DEPP, DT,
     &                           NEW_POPP, DEP_POPP )
!
! !USES:
!
      USE ERROR_MOD,    ONLY : ERROR_STOP
!
! !INPUT PARAMETERS: 
!
      REAL(fp),  INTENT(IN)  :: OLD_POPP
      REAL(fp),  INTENT(IN)  :: K_DEPP
      REAL(fp),  INTENT(IN)  :: DT
!
! !OUTPUT PARAMETERS:
!
      REAL(fp),  INTENT(OUT) :: NEW_POPP
      REAL(fp),  INTENT(OUT) :: DEP_POPP
!
! !REVISION HISTORY: 
!  9 February 2011 - CL Friedman - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      REAL(fp)               :: E_KDEPP_T

      !=================================================================
      ! NO_RXN_WITHDEP begins here!
      !=================================================================

      ! Precompute exponential factors [dimensionless]
      E_KDEPP_T = EXP( -K_DEPP * DT )     

      !=================================================================
      ! No Chemistry, Deposition only
      !=================================================================

      ! New mass of POPP [kg]
      NEW_POPP = OLD_POPP * E_KDEPP_T

      ! POPP deposition [kg]
      DEP_POPP = OLD_POPP - NEW_POPP
      DEP_POPP = MAX( DEP_POPP, 0e+0_fp )


      END SUBROUTINE NO_RXN_WITHDEP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  get_oh
!
! !DESCRIPTION: Function GET\_OH returns monthly mean OH and imposes a diurnal
! variation. 
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_OH( I, J, L, State_Met ) RESULT( OH_MOLEC_CM3 )
!
! !USES:
!
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TIME_MOD,           ONLY : GET_TS_CHEM

      USE CMN_SIZE_MOD             ! Size parameters
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)         :: I, J, L

      TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
!
! !REMARKS:
! Copied GET_OH function from mercury_mod.f - CLF
!
! !REVISION HISTORY: 
!  03 Feb 2011 - CL Friedman - Initial Version
!  29 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  23 Sep 2014 - M. Sulprizio- Now get OH for offline aerosol sim from HEMCO
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
       REAL(fp)        :: OH_MOLEC_CM3
    
       !=================================================================
       ! GET_OH begins here!
       !=================================================================

       ! Test for sunlight...
       IF ( State_Met%SUNCOS(I,J) > 0e+0_fp 
     &    .and. TCOSZ(I,J) > 0e+0_fp ) THEN

         ! Impose a diurnal variation on OH during the day
         OH_MOLEC_CM3 = OH(I,J,L)                      *           
     &                  ( State_Met%SUNCOS(I,J) / TCOSZ(I,J) ) *
     &                  ( 1440e+0_fp        / GET_TS_CHEM() )

         ! OH is in kg/m3 (from HEMCO), convert to molec/cm3 (mps, 9/23/14)
         OH_MOLEC_CM3 = OH_MOLEC_CM3 * XNUMOL_OH / CM3PERM3

         ! Make sure OH is not negative
         OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0e+0_fp )
               
       ELSE

         ! At night, OH goes to zero
         OH_MOLEC_CM3 = 0e+0_fp

       ENDIF

       END FUNCTION GET_OH
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  ohno3time
!
! !DESCRIPTION: Subroutine OHNO3TIME computes the sum of cosine of the solar
!  zenith angle over a 24 hour day, as well as the total length of daylight. 
!  This is needed to scale the offline OH and NO3 concentrations.
!  (rjp, bmy, 12/16/02, 12/8/04)
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OHNO3TIME
!
! !USES:
!
      USE GRID_MOD, ONLY : GET_XMID,    GET_YMID_R
      USE TIME_MOD, ONLY : GET_NHMSb,   GET_ELAPSED_SEC
      USE TIME_MOD, ONLY : GET_TS_CHEM, GET_DAY_OF_YEAR, GET_GMT

      USE CMN_SIZE_MOD   ! Size parameters
      USE CMN_GCTM_MOD   ! Physical constants
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version for POPS_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE       :: FIRST = .TRUE.
      INTEGER             :: I, IJLOOP, J, L, N, NT, NDYSTEP
      REAL(fp)              :: A0, A1, A2, A3, B1, B2, B3
      REAL(fp)              :: LHR0, R, AHR, DEC, TIMLOC, YMID_R
      REAL(fp)              :: SUNTMP(MAXIJ)
      
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

      DEC = A0 - A1*cos(  R) + B1*sin(  R)
     &         - A2*cos(2*R) + B2*sin(2*R)
     &         - A3*cos(3*R) + B3*sin(3*R)

      LHR0 = int(float( GET_NHMSb() )/10000.)

      ! Only do the following at the start of a new day
      IF ( FIRST .or. GET_GMT() < 1e-5 ) THEN 
      
         ! Zero arrays
         TTDAY(:,:) = 0e+0_fp
         TCOSZ(:,:) = 0e+0_fp
         COSZM(:,:) = 0e+0_fp

         ! NDYSTEP is # of chemistry time steps in this day
         NDYSTEP = ( 24 - INT( GET_GMT() ) ) * 60 / GET_TS_CHEM()      

         ! NT is the elapsed time [s] since the beginning of the run
         NT = GET_ELAPSED_SEC()

         ! Loop forward through NDYSTEP "fake" timesteps for this day 
         DO N = 1, NDYSTEP
            
            ! Zero SUNTMP array
            SUNTMP(:) = 0e+0_fp

            ! Loop over surface grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, YMID_R, IJLOOP, TIMLOC, AHR )
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Grid box latitude center [radians]
               YMID_R = GET_YMID_R( I, J, 1 )

               ! Increment IJLOOP
               IJLOOP = ( (J-1) * IIPAR ) + I
               TIMLOC = real(LHR0) + real(NT)/3600.0 +
     &                  GET_XMID( I, J, 1 )/15.0
         
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
               SUNTMP(IJLOOP) = sin(YMID_R) * sin(DEC) +
     &                          cos(YMID_R) * cos(DEC) * cos(AHR)

               ! TCOSZ is the sum of SUNTMP at location (I,J)
               ! Do not include negative values of SUNTMP
               TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP(IJLOOP), 0e+0_fp )

               ! COSZM is the peak value of SUMTMP during a day at (I,J)
               ! (rjp, bmy, 3/30/04)
               COSZM(I,J) = MAX( COSZM(I,J), SUNTMP(IJLOOP) )

               ! TTDAY is the total daylight time at location (I,J)
               IF ( SUNTMP(IJLOOP) > 0e+0_fp ) THEN
                  TTDAY(I,J) = TTDAY(I,J) + DBLE( GET_TS_CHEM() )
               ENDIF
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! Increment elapsed time [sec]
            NT = NT + ( GET_TS_CHEM() * 60 )             
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
! !IROUTINE:  init_pops
!
! !DESCRIPTION: Subroutine INIT\_POPS allocates and zeroes all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_POPS( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE ERROR_MOD,          ONLY : DEBUG_MSG
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  20 Sep 2010 - N.E. Selin  - Initial Version based on INIT_MERCURY
!  25 Mar 2013 - R. Yantosca - Now accept Input_Opt argument
!  25 Aug 2014 - M. Sulprizio- Now accept am_I_Root and RC arguments
!  10 Sep 2014 - M. Sulprizio- Remove calls to GET_GLOBAL_OC and GET_GLOBAL_BC,
!                              these fields are now read in HEMCO
!  10 Mar 2015 - R. Yantosca - Now always allocate ZERO_DVEL
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER       :: AS

      ! For fields from Input_Opt
      LOGICAL       :: LDRYD, LNLPBL, LPRT

      !=================================================================
      ! INIT_POPS begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Copy values from Input_Opt
      LDRYD     = Input_Opt%LDRYD
      LNLPBL    = Input_Opt%LNLPBL
      LPRT      = Input_Opt%LPRT

      !=================================================================
      ! Allocate and initialize arrays
      !=================================================================
      ALLOCATE( COSZM( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'COSZM' )
      COSZM = 0e+0_fp

      ALLOCATE( TCOSZ( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TCOSZ' )
      TCOSZ = 0e+0_fp

      ALLOCATE( TTDAY( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'TTDAY' )
      TTDAY = 0e+0_fp

      ! Now always allocate ZERO_DVEL.  This is needed because the calls
      ! to CHEM_POPGP have been modified in the CHEMPOPS routine above.
      ! (bmy, 3/10/15)
      ALLOCATE( ZERO_DVEL( IIPAR, JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'ZERO_DVEL' )
      ZERO_DVEL = 0e+0_fp

      END SUBROUTINE INIT_POPS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  cleanup_pops
!
! !DESCRIPTION: Subroutine CLEANUP\_POPS deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_POPS
!
! !REVISION HISTORY: 
!  20 September 2010 - N.E. Selin - Initial Version
!EOP
!------------------------------------------------------------------------------
!BOC

      IF ( ALLOCATED( COSZM      ) ) DEALLOCATE( COSZM      )     
      IF ( ALLOCATED( TCOSZ      ) ) DEALLOCATE( TCOSZ      )
      IF ( ALLOCATED( TTDAY      ) ) DEALLOCATE( TTDAY      )
      IF ( ALLOCATED( ZERO_DVEL  ) ) DEALLOCATE( ZERO_DVEL  )

      ! Free pointers
      IF ( ASSOCIATED( C_OC      ) ) C_OC => NULL()
      IF ( ASSOCIATED( C_BC      ) ) C_BC => NULL()
      IF ( ASSOCIATED( O3        ) ) O3   => NULL()
      IF ( ASSOCIATED( OH        ) ) OH   => NULL()

      END SUBROUTINE CLEANUP_POPS
!EOC
      END MODULE POPS_MOD

