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
  USE PhysConstants

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
  SUBROUTINE fullChem_SetStateHet( I,         J,         L,                  &
                                   id_SALA,   id_SALAAL, id_SALC,            &
                                   id_SALCAL, Input_Opt, State_Chm,          &
                                   State_Met, H,         RC                 )
!
! !USES:
!
    USE ErrCode_Mod
    USE GcKpp_Global
    USE GcKpp_Parameters
    USE PhysConstants,    ONLY : AVO, PI, AIRMW
    USE Input_Opt_Mod,    ONLY : OptInput
    USE rateLawUtilFuncs
    USE State_Chm_Mod,    ONLY : ChmState, Ind_
    USE State_Met_Mod,    ONLY : MetState
    USE Henry_Mod,        ONLY : Calc_KH
    USE Henry_Mod,        ONLY : Calc_Heff
    USE Species_Mod,      ONLY : SpcConc
    
  ! Species ID flags
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I          ! Lon (or X-dim) gridbox index
    INTEGER,        INTENT(IN)    :: J          ! Lat (or Y-dim) gridbox index
    INTEGER,        INTENT(IN)    :: L          ! Vertical level index
    INTEGER ,       INTENT(IN)    :: id_SALA    ! Indices of SALA, SALAAL
    INTEGER,        INTENT(IN)    :: id_SALAAL  !  SALC, and SALCAL species
    INTEGER,        INTENT(IN)    :: id_SALC    !  in the State_Chm%Species
    INTEGER,        INTENT(IN)    :: id_SALCAL  !  object
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm  ! Chemistry State object
    TYPE(MetState), INTENT(IN)    :: State_Met  ! Meterology State object
!
! INPUT/OUTPUT PARAMETERS:
!
    TYPE(HetState), INTENT(INOUT) :: H          ! Hetchem State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure?

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: NA
    REAL(dp), PARAMETER   :: R = 0.08205e+0_dp
    INTEGER,   PARAMETER  :: SUL            = 8  ! Tropospheric Sulfate!
    REAL(dp) :: K0, TK
    REAL(dp) :: FeIII, MnII
    REAL(dp) :: XSO2aq_a, XSO3_a, XHSO3_a
    REAL(dp) :: DUST,  Mn_ant,  Mn_nat
    REAL(dp) :: Mn_tot, Mn_d,    Fe_d
    REAL(dp) :: Fe_ant, Fe_nat,  Fe_tot
    REAL(dp) :: Fe_d_ant, Fe_d_nat
    REAL(dp) :: HSO3aq_a, SO3aq_a, Hplus_a, pH_a, ALWC
    REAL(dp) :: SO2aq_a, FeIII_a, MnII_a, Mn_d_a, Fe_d_a
    REAL(dp) :: Fe_d_ant_a, Fe_d_nat_a
    REAL(dp) :: MACOEFF_SO2, MACOEFF_H2O2, MACOEFF_O3, MACOEFF_NO2
    REAL(dp) :: MACOEFF_CH2O
    REAL(dp) :: val1NO2, val1O3, val1H2O2, val1TMI, val1CH2O
    REAL(dp) :: val2NO2, val2O3, val2H2O2, val2TMI, val2CH2O
    REAL(dp) :: val3NO2, val3O3, val3H2O2, val3TMI, val3CH2O
    REAL(dp) :: kchemNO2, kchemO3, kchemH2O2, kchemTMI, kchemCH2O
    REAL(dp) :: gammaO3, gammaH2O2, gammaNO2, gammaTMI, gammaCH2O
    REAL(dp) :: Kw1, KHMS, KHMS2, dOH
    REAL(dp) :: f_x, y, corr, x, xradi, l_r
    REAL(dp) :: D_aSO2, D_aO3, D_aNO2, D_aCH2O
    REAL(dp) :: HEFF, HEFF_H2O2, speed, M_SO2
    !REAL(dp) :: Eff_Fe, Eff_Mn, Eff_H2O2
    REAL(dp) :: pKa, CR, KH_NO2, KH_O3, KH_SO2, KH_H2O2
    REAL(dp) :: IONIC_a, IONIC_aMAX, IONIC_bMAX, IONIC_cMAX, IONIC_eMAX
    REAL(dp) :: kH2O2,kH_CH2O, Khc1, kNO2, kTMI6, kTMI7
    REAL(dp) :: KO0, KO1, KO2, KspFeOH3, KspMnOH2
    REAL(dp) :: hydroxide, FeIII_max, MnII_max, SIV_a
    REAL(dp) :: Ks1, Ks2, HCSO2_a
    REAL(dp) :: XSO2g_a, PATM, SO2, CNVFAC, RHO,RHO_num, id_pFe
    REAL(dp) :: ff, k9, k10,A, B, Beta, b1
    INTEGER  :: id_NO2, id_O3, id_SO2, id_DST1, id_DST2, id_DST3, id_DST4
       ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)
    pH_a                        = State_Chm%IsorropAeropH(I,J,L,1)
    Hplus_a                     = 10**(-1.0_dp * pH_a)
    hydroxide = 10.0_fp**(-14.0_fp + pH_a) ! OH- concentration, M
    KspFeOH3 = 2.6e-38_fp
    KspMnOH2 = 1.6e-13_fp
    !========================================================================
    ! Populate fields of the HetState object in gckpp_Global
    !========================================================================

    ! Initialization
    RC = GC_SUCCESS
    NA = State_Chm%nAeroType
    Spc                         => State_Chm%Species
    RHO                         = State_Met%AIRDEN(I,J,L)
    RHO_num                     = State_Met%AIRNUMDEN(I,J,L)
    CNVFAC    = 1.E3_fp * AIRMW / ( RHO * AVO ) !mcl/cm3->v/v
    !========================================================================
    ! Populate fields of the HetState object in gckpp_Global
    !========================================================================
    id_SO2    = Ind_( 'SO2'    )
    id_DST1   = Ind_( 'DST1'   )
    id_DST2   = Ind_( 'DST2'   )
    id_DST3   = Ind_( 'DST3'   )
    id_DST4   = Ind_( 'DST4'   )
    id_pFe    = Ind_( 'pFe'    )

    ! Identify a box for debug printout within rate-law functions
    debugBox        = .FALSE.

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
    H%IONIC         = State_Chm%IsorropIONIC(I,J,L)
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

    ! Concentrations from ISORROPIA/HETP
    H%HSO4_molal    = State_Chm%IsorropBisulfate(I,J,L)
    H%NO3_molal     = State_Chm%IsorropNitrate(I,J,L,1)
    H%SO4_molal     = State_Chm%IsorropSulfate(I,J,L)

    ! pH and alkalinity fields
    !H%H_plus        = State_Chm%IsorropHplus(I,J,L,1)
    H%pHCloud       = State_Chm%pHCloud(I,J,L)
    H%pHSSA(:)      = State_Chm%IsorropAeropH(I,J,L,:)
    H%H_conc_Sul    = 10.0**( -1.0_dp * H%pHSSA(1) )
    H%H_conc_LCl    = 10.0**( -1.0_dp * H%pHCloud  )
    H%H_conc_ICl    = 10.0**( -4.5_dp              )
    H%H_conc_SSA    = H%H_conc_Sul
    H%H_conc_SSC    = 10.0**( -5.0_dp              )
    H%f_Alk_SSA     = SafeDiv( State_Chm%Species(id_SALAAL)%Conc(I,J,L),     &
                               State_Chm%Species(id_SALA  )%Conc(I,J,L),     &
                               0.0_dp                                       )
    H%f_Alk_SSA     = MAX( MIN( H%f_Alk_SSA, 1.0_dp ), 0.0_dp )
    H%f_Acid_SSA    = 1.0_dp - H%f_Alk_SSA
    H%f_Alk_SSC     = SafeDiv( State_Chm%Species(id_SALCAL)%Conc(I,J,L),     &
                               State_Chm%Species(id_SALC  )%Conc(I,J,L),     &
                               0.0_dp                                       )
    H%f_Alk_SSC     = MAX( MIN( H%f_Alk_SSC, 1.0_dp ), 0.0_dp )
    H%f_Acid_SSC    = 1.0_dp - H%f_Alk_SSC
    H%SSA_is_Alk    = ( ABS( H%f_Alk_SSA ) > 0.01_dp )
    H%SSA_is_Acid   = ( .not.  H%SSA_is_Alk          )
    H%SSC_is_Alk    = ( ABS( H%f_Alk_SSC ) > 0.01_dp )
    H%SSC_is_Acid   = ( .not.  H%SSC_is_Alk          )

    ! Aerosol fields: HSO3 and SO3 concentrations in aerosol [mol/L]
    H%HSO3_aq_a       = State_Chm%HSO3_aq_a(I,J,L)
    H%SO3_aq_a        = State_Chm%SO3_aq_a(I,J,L)
    H%SO2_aq_a        = State_Chm%SO2_aq_a(I,J,L)

    H%SO2             = C(ind_SO2)*CNVFAC
    
    ! Other fields
    H%gamma_HO2     = Input_Opt%gamma_HO2

      ! Get Henry's law parameters
    K0   = HENRY_K0(ind_SO2)!SpcInfo%Henry_K0
    CR   = HENRY_CR(ind_SO2)!SpcInfo%Henry_CR
    pKa = 1.81_dp
    TK = State_Met%T(I,J,L)
    
    ! Calculate the Henry's law constant
   ! CALL CALC_KH( K0, CR, TK, K0, RC )

    ! Calculate effective Henry's law constant, corrected for pH
    ! (for those species that have a defined pKa value)
  !  CALL CALC_HEFF( pKa, H%pHSSA(1), K0, HEFF, RC )
  !  H%HEFF_a = HEFF
    Ks1    = 1.30e-2_dp * EXP( 6.75e+0_dp * ( 298.15e+0_dp / TK - 1.e+0_dp ) )
    Ks2    = 6.31e-8_dp * EXP( 5.05e+0_dp * ( 298.15e+0_dp / TK - 1.e+0_dp ) )

   ! Henry's constant [mol/l-atm] and Effective Henry's constant for SO2
    HCSO2_a  = 1.22e+0_dp * EXP( 10.55e+0_dp * ( 298.15e+0_dp / Tk - 1.e+0_dp) )
    H%HEFF_a = HCSO2_a * (1.e+0_dp + (Ks1/Hplus_a) + (Ks1*Ks2 / (Hplus_a*Hplus_a)))

     !===============================================================
    !  Do sulfur fraction for aerosols (krt, 1/11/24)
    !=================================================================
    PATM = State_Met%PMID_DRY( I, J, L ) / ( ATM * 1.e-2_dp ) ! Press, dry [atm]
    ! SIV Fraction for aerosols
    XSO2aq_a = 1.e+0_dp/(1.e+0_dp + Ks1/Hplus_a + Ks1*Ks2/(Hplus_a*Hplus_a))
    XHSO3_a  = 1.e+0_dp/(1.e+0_dp + Hplus_a/Ks1 + Ks2/Hplus_a)
    XSO3_a   = 1.e+0_dp/(1.e+0_dp + Hplus_a/Ks2 + Hplus_a*Hplus_a/(Ks1*Ks2))

     ! For aerosols (krt)
    ALWC         = State_Chm%IsorropAeroH2O(I,J,L,1) ! ug/m3 air
    ! convert ALWC from ug/m3 to m3/m3
    ALWC = ALWC * 1e-9_dp * 1e-3_dp
    
    ! QL can sometimes be negative, so force LWC to be positive
    ALWC = MAX( 0.0_dp, ALWC )
    
    XSO2g_a  = 1.e+0_dp / ( 1.e+0_dp + ( H%HEFF_a * R * Tk * ALWC ) )
    SIV_a    = H%HEFF_a * XSO2g_a * H%SO2 * PATM

    ! Effective HSO3aq 
    HSO3aq_a = SIV_a * XHSO3_a           ! unit: M (qjc, 06/10/16)

    ! Effective SO3aq for HOBr+SO3
    SO3aq_a  = SIV_a * XSO3_a            ! unit: M (qjc, 06/10/16)

    SO2aq_a = SIV_a * XSO2aq_a
    ! == Calculate gamma for uptake of SO2 to Sulfate and organic aerosol ====
    MACOEFF_SO2  = 0.11_dp
    MACOEFF_NO2  = 2e-4_dp
    MACOEFF_O3   = 2e-3_dp
    MACOEFF_H2O2 = 0.11_dp
    MACOEFF_CH2O = 0.04_dp
    ! Ionic strength from HETP (M)
    IONIC_a = State_Chm%IsorropIONIC(I,J,L) 
    ! Molecular speed
    M_SO2 = State_Chm%SpcData(id_SO2)%Info%MW_g   * 1.0e-3_dp
    speed = ( SQRT( EIGHT_RSTARG_T / ( PI * M_SO2 ) ))*100.0_dp
    ! Aqueous phase diffusion 
    D_aSO2  = 1.32e-5_dp
    D_aO3   = 1.10e-2_dp*exp(-1896_dp/Tk)
    D_aNO2  = 2.0e-5_dp
    D_aCH2O = 1.52e-6_dp ! no idea
    
    ! Get FeIII
    FeIII_a   = 0.0
    FeIII_Max = 0.0
    MnII_a    = 0.0
    MnII_Max  = 0.0
    ! Metal catalyzed oxidation of SO2 pathway
    !--------------------------------------------------------
   ! Get dust concentrations [MND -> ng/m3]

    DUST = ( Spc(id_DST1)%Conc(I,J,L)*0.7_dp + Spc(id_DST2)%Conc(I,J,L) + &
         Spc(id_DST3)%Conc(I,J,L) + Spc(id_DST4)%Conc(I,J,L) )            &
         * 1.e+15_dp * State_Chm%SpcData(id_DST1)%Info%MW_g / AVO
      
    ! Calculate Fe and Mn natural [ng m-3]
    ! Assume that Fe is 3.5% of total dust mass based on
    ! Taylor and McLennan [1985]
    Fe_nat = DUST * 35e-3_dp
    ! and Mn is 50 times less than Fe based on Desbouefs et al.[2005]
    !Mn_nat = Fe_nat / 50e+0_dp
    ! Use Cai et al
    Mn_nat = DUST * 3E-3_dp
    
    ! Anthropogenic Fe concentrations [mcl/cm3 -> ng/m3]
    Fe_ant = Spc(id_pFe)%Conc(I,J,L) * CNVFAC * &
         1.e+12_dp * State_Met%AD(I,J,L) &
         / ( AIRMW / State_Chm%SpcData(id_pFe)%Info%MW_g ) &
         / State_Met%AIRVOL(I,J,L)
    !             Fe_ant = Spc(id_pFe)%Conc(I,J,L) * 1.e+15_fp * &
    !                  State_Chm%SpcData(id_DST1)%Info%MW_g / AVO
 
    ! Calculate Mn anthropogenic [ng m-3]
    ! assume anthropogenic Mn is 1/30 times anthropogenic Fe
    Mn_ant = Fe_ant / 10e+0_dp
    
    ! Calculate total Mn and Fe [ng m-3]
    Mn_tot = Mn_ant + Mn_nat
    Fe_tot = Fe_ant + Fe_nat
    
    ! Convert Mn and Fe [ng m-3] to [mole l-1]
    ! Below I guess we are making the assumption that we wont have LWC and ALWC at the same time in one grid box?
    ! Are the aerosols scavenged before this occurs?  Are we double counting?
    !krt, for aerosols
    
    Mn_d_a             = 0e+0_dp
    IF (ALWC > 0e+0_dp) THEN
       ! Units: ng/m3 * (g/ng) / (g/mol) / (m3 H2O / m3 air) * (m3/L)
       ! max possible dissolved Mn
       MnII_max = Mn_tot * 1.0e-9_dp / 54.94e+0_dp / ALWC * 1.0e-3_dp
       MnII_max = MnII_max *0.5e+0_fp ! Maximum fractional solubility

       Mn_d_a = KspMnOH2/(hydroxide**2.0_dp)
       ! can't have more dissolved Mn than available in gas phase
       Mn_d_a = MIN(Mn_d_a, MnII_max) ! can't have more dissolved Mn than available in gas phase
       
    ENDIF
          
    ! Solubility of Fe is 10% for anthropogenic, and 1% for dust
    Fe_d_a     = 0e+0_fp
    IF ( ALWC > 0e+0_fp ) THEN
       Fe_d_ant_a = Fe_ant * 1e-9_fp / &
            State_Chm%SpcData(id_pFe)%Info%MW_g / &
            ALWC * 1.0e-3_fp
       Fe_d_nat_a = Fe_nat * 1.0e-9_fp / &
            State_Chm%SpcData(id_pFe)%Info%MW_g / &
            ALWC * 1.0e-3_fp
       ! currently not sure how to account for different
       ! solubility between anthropogenic and natural
       ! Max possible solubility
       FeIII_Max = Fe_d_ant_a*0.1e+0_fp + Fe_d_nat_a*0.1e+0_fp

       IF ( State_Met%SUNCOS(I,J) > 0e+0_fp ) THEN
          ! Assume 10% of dissolved Fe is in Fe(III)
          !oxidation state during the daytime
          FeIII_Max = FeIII_Max * 0.1e+0_fp
       ELSE
          ! Assume 90% of dissolved Fe is in Fe(III)
          ! oxidation state during the nighttime
          FeIII_Max = FeIII_Max * 0.9e+0_fp
       ENDIF

       Fe_d_a = KspFeOH3/(hydroxide**3.0_fp)
       Fe_d_a = MIN(Fe_d_a, FeIII_Max)
    ENDIF
    
    ! Assume that dissolved Mn is in Mn(II) oxidation state all of
    ! the time
    MnII = Mn_d
    MnII_a = Mn_d_a
    FeIII_a = Fe_d_a
    
    ! Add to diagnostic
    State_Chm%FeIII_A(I,J,L) = FeIII_a 
    State_Chm%MnII_A(I,J,L)  = MnII_a
    State_Chm%FeIII_AMAX(I,J,L) = FeIII_Max
    State_Chm%MnII_AMAX(I,J,L)  = MnII_Max
    H%FeIII_a       = State_Chm%FEIII_A(I,J,L)
    H%MnII_a        = State_Chm%MNII_A(I,J,L)

    ! ==== TMI/O2 ====
    ! First order in SO2
    val1TMI = ( 1.0_dp / MACOEFF_SO2 )

    !IF (IONIC_a .gt. 2.0e+0_fp ) THEN
    !   IONIC_aMAX = 2.0e+0_fp
    !ENDIF
    IONIC_aMAX = IONIC_a
    IF (IONIC_a .gt. 1.0e+0_fp ) THEN
       IONIC_aMAX = 1.0e+0_fp
    ENDIF

    ! Ionic strenght correction
    !Eff_Mn = 10.0_dp**(-4.0_dp*(SQRT(IONIC_aMAX)/(1.0_dp+SQRT(IONIC_aMAX))))
    !Eff_Fe = 10.0_dp**(-4.0_dp*(SQRT(IONIC_aMAX)/(1.0_dp+SQRT(IONIC_aMAX))))

    ff = 10.0_dp**(-4.0_dp*(SQRT(IONIC_aMAX)/(1.0_dp+SQRT(IONIC_aMAX))))

    kTMI6 = 3.72e+7_dp*exp(-8431.6_dp*(1.0e+0_dp/Tk - 1.0e+0_dp/297.0_dp)) ! M-2 s-1
    kTMI7 = 2.51e+13_dp* exp(-8431.6_dp*(1.0e+0_dp/Tk - 1.0e+0_dp/297.0_dp))
    
    !IF ( pH_a .le. 4.2_dp ) THEN
    !   kchemTMI = kTMI6*Hplus_a**(-0.74_dp)* &
    !        (MnII_a*FeIII_a*Eff_Fe*Eff_Mn)
    !ELSE
    !   kchemTMI = kTMI7*Hplus_a**(0.67_dp)* &
    !        (MnII_a*FeIII_a*Eff_Fe*Eff_Mn)
    !ENDIF
    IF ( pH_a .le. 4.2_dp ) THEN 
       kchemTMI = kTMI6*Hplus_a**(-0.74_dp)*(MnII_a*FeIII_a)
    ELSE
       kchemTMI = kTMI7*Hplus_a**(0.67_dp)*(MnII_a*FeIII_a)
    ENDIF

    ! IONIC strength impact
    kchemTMI = kchemTMI * ff
    
    M_SO2 = State_Chm%SpcData(id_SO2)%Info%MW_g   * 1.0e-3_dp    

    ! Get Henry's law parameters
    ! No idea why this is not working so hard coding...
    K0   = 1.22_dp! HENRY_K0(ind_SO2)!SpcInfo%Henry_K0
    CR   = 3100.0_dp!HENRY_CR(ind_SO2)!SpcInfo%Henry_CR
    pKa = 1.81_dp
    D_aSO2 = 1.32e-5_dp
    
    !CALL CALC_KH( K0, CR, TK, K0, RC )
    !CALL CALC_HEFF( pKa, pH_a, K0, HEFF, RC )
    
    Ks1    = 1.30e-2_fp * EXP( 6.75e+0_fp * ( 298.15e+0_fp / TK - 1.e+0_fp ) )
    Ks2    = 6.31e-8_fp * EXP( 5.05e+0_fp * ( 298.15e+0_fp / TK - 1.e+0_fp ) )
    
    ! Henry's constant [mol/l-atm] and Effective Henry's constant for SO2
    HCSO2_a  = 1.22e+0_fp * EXP( 10.55e+0_fp * ( 298.15e+0_fp / TK - 1.e+0_fp) )
    HEFF = HCSO2_a * (1.e+0_fp + (Ks1/Hplus_a) + (Ks1*Ks2 / (Hplus_a*Hplus_a)))
    
    val2TMI = speed/((FOUR_RGASLATM_T * HEFF)*SQRT(D_aSO2*kchemTMI))
    
    l_r = SQRT(D_aSO2/kchemTMI)
    xRadi = State_Chm%AeroRadi(I,J,L,SUL) 
    val3TMI = ReactoDiff_Corr(xRadi, l_r)
    
    gammaTMI = 1.0_dp/(val1TMI + val2TMI*1/val3TMI)      ! Uptake by various aerosol types
    
    ! ===== H2O2 =====
    IONIC_eMAX = IONIC_a
    IF (IONIC_eMAX .gt. 5.0e+0_fp ) THEN
       IONIC_eMAX = 5.0e+0_fp
    ENDIF
    
    ! First order in SO2
    val1H2O2 = (1.0e+0_dp / MACOEFF_SO2)

    kH2O2 = 7.45e+7_dp *exp(-4430_dp*(1.0_dp/Tk - 1.0_dp/298.0_dp))
    kchemH2O2 =(kH2O2*Hplus_a*HSO3aq_a)/(1.0_dp+13.0_dp*Hplus_a)

    ! IONIC strength impact (Cai et al., 2024)
    A = 0.509_dp
    B = 0.17_dp
    Beta = 0.18_dp
    ff = 10.0_dp**(-(2.0_dp*A*sqrt(IONIC_eMAX))/(1+B*sqrt(IONIC_eMAX)) + 2.0_dp*Beta*IONIC_eMAX)
    kchemH2O2 = kchemH2O2 * ff

    l_r = SQRT(D_aSO2/kchemH2O2)
    xRadi = State_Chm%AeroRadi(I,J,L,SUL) 

    val3H2O2 = ReactoDiff_Corr(xRadi, l_r)

    ! Calculate the Henry's law constant
    K0   = 8.3e+4_dp! HENRY_K0(ind_H2O2)!SpcInfo%Henry_K0
    CR   = 7400.0_dp!HENRY_CR(ind_H2O2)!SpcInfo%Henry_CR
    pKa  = 11.75_dp
    CALL CALC_KH( K0, CR, TK, KH_H2O2, RC )
    ! Calculate effective Henry's law constant, corrected for pH
    ! (for those species that have a defined pKa value)
    CALL CALC_HEFF( pKa, pH_a, KH_H2O2, HEFF_H2O2, RC )

    ! Calculate the second uptake parameterization term:
    val2H2O2 = speed/((FOUR_RGASLATM_T * HEFF_H2O2)*SQRT(D_aSO2*kchemH2O2))
    gammaH2O2 = 1.0_dp/(val1H2O2 + val2H2O2*1/val3H2O2)      ! Uptake by various aerosol ty
    
    ! ==== O3 ====
    IONIC_cMAX = IONIC_a
    IF (IONIC_cMAX .gt. 1.2e+0_fp ) THEN
       IONIC_cMAX = 1.2e+0_fp
    ENDIF
    
    ! First order in O3
    val1O3 = ( 1.0e+0_dp / MACOEFF_O3 )

    K0   = 0.0101325e+0_dp!SpcInfo%Henry_K0
    CR   = 2800e+0_dp !SpcInfo%Henry_CR

    CALL CALC_KH( K0, CR, TK, KH_O3, RC )

    KO0 = 2.40e+4_dp
    kO1 = 3.49e+12_dp * EXP( -4.83e+3_dp / Tk )
    kO2 = 7.32e+14_dp * EXP( -4.03e+3_dp / Tk )

    kchemO3 = KO0*SO2aq_a + KO1*HSO3aq_a + KO2*SO3aq_a
    ! IONIC strength impact (Cai et al., 2024)
    b1 = 1.94_dp
    ff = 1.0_dp + b1 * IONIC_cMAX
    kchemO3 = kchemO3 * ff
    
    val2O3 = speed/((FOUR_RGASLATM_T * KH_O3)*SQRT(D_aO3*kchemO3))

    l_r = SQRT(D_aO3/kchemO3)

    xRadi = State_Chm%AeroRadi(I,J,L,SUL) 

    val3O3 = ReactoDiff_Corr(xRadi, l_r)

    gammaO3 = 1.0_dp/(val1O3 + val2O3*1/val3O3)      ! Uptake by various aerosol ty
    
    ! ==== NO2 ====
    ! First order in NO2
    val1NO2 = (1.0e+0_dp / MACOEFF_NO2)

    K0   = 0.012159e+0_dp!SpcInfo%Henry_K0
    CR   = 2400e+0_dp !SpcInfo%Henry_CR

    CALL CALC_KH( K0, CR, TK, KH_NO2, RC )

    kNO2 = 0.0_dp
    IF ( pH_a .lt. 5.0_dp ) THEN
       kNO2 = 2.0e+6_dp
    ELSE IF (pH_a .gt. 5.8_dp ) THEN
       kNO2 = 1.67e+7_dp
    ELSE
       ! linearly interpolate
       kNO2 = 2.0E+6_dp + (pH_a - 5.0_dp) * &
            (1.67E+7_dp-2.0E+6_dp)/(5.8_dp-5.0_dp)
    ENDIF
    
    kchemNO2 =(kNO2*SIV_a)

    l_r = SQRT(D_aNO2/kchemNO2)

    xRadi = State_Chm%AeroRadi(I,J,L,SUL) 

    val3NO2 = ReactoDiff_Corr(xRadi, l_r)

    val2NO2 = speed/((FOUR_RGASLATM_T * KH_NO2)*SQRT(D_aNO2*kchemNO2))

    gammaNO2 = 1.0_dp/(val1NO2 + val2NO2*1/val3NO2)      ! Uptake by various aerosol ty

      ! ==== CH2O ====
    ! First order in CH2O
    IONIC_bMAX = IONIC_a
    IF (IONIC_bMAX .gt. 16.0e+0_fp ) THEN
       IONIC_bMAX = 16.0e+0_fp
    ENDIF
    
    val1CH2O = ( 1.0e+0_dp / MACOEFF_CH2O )

    !K0   = 3.24e+3_dp!SpcInfo%Henry_K0
    !CR   = 6800.0e+0_dp !SpcInfo%Henry_CR
    !CALL CALC_KH( K0, CR, TK, KH_CH2O, RC )
  
    ! I was using too high a value, not considering that most HCHO is the diol
    KH_CH2O  = 2.5e+0_fp * EXP( 21.6e+0_fp &
         *  (298.15e+0_fp / Tk - 1.e+0_fp) )
    
    !=================================================================
    ! Ratio of methanediol to aqeuous HCHO if I want to include this chemistry
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
         * ( 298.15e+0_fp / Tk - 1.e+0_fp ) )

    ! not sure where I got this from!
    ! Initial rates are from Cai et al. But T dependence?
    k9  = 7.9e+2_dp * EXP( -4.9e+3_dp / Tk ) ! M/s
    k10 = 2.5e+7_dp * EXP( -1.8e+3_dp / Tk ) ! M/s
    ! Now taken from HMS implementation in fullchem_SulfurChemFuncs.F90
    !k9 = 7.9e+2_fp * EXP( -16.44e+0_fp & ! L/mol/s
    !     * ( 298.15e+0_fp / Tk - 1.e+0_fp ) )
    !k10 = 2.5e+7_fp * EXP( -6.04e+0_fp  & ! L/mol/s
    !     * ( 298.15e+0_fp / Tk - 1.e+0_fp ) )
    
    kchemCH2O =  k9*HSO3aq_a + k10*SO3aq_a
    ! IONIC strength impact (Cai et al., 2024)
    A = 0.509_dp
    B = 0.17_dp
    Beta = 0.18_dp
    ff = 10.0_dp**(-(2.0_dp*A*sqrt(IONIC_bMAX))/ &
         (1.0_dp+B*sqrt(IONIC_bMAX)) + 2.0_dp*Beta*IONIC_bMAX)
    kchemCH2O = kchemCH2O * ff
    
    val2CH2O = speed/((FOUR_RGASLATM_T * KH_CH2O)*SQRT(D_aCH2O*kchemCH2O))

    l_r = SQRT(D_aCH2O/kchemCH2O)

    xRadi = State_Chm%AeroRadi(I,J,L,SUL) 

    val3CH2O = ReactoDiff_Corr(xRadi, l_r)

    gammaCH2O = 1.0_dp/(val1CH2O + val2CH2O*1/val3CH2O)      ! Uptake by various aerosol ty

    ! Now do rate for HMS decomposition back to SO2 + CH2O
    Kw1 = 1e-14_fp * EXP( -22.51e+0_fp &
         *  ( 298.15e+0_fp / Tk - 1.e+0_fp ) )
    
    ! Conversion rate from HMS to SO2 via reaction with OH-
    ! (jmm, 06/15/18; MSL 1/18/22)
    KHMS   = 3.6e+3_fp * EXP( -15.09e+0_fp & ! L/mol/s
         * ( 298.15e+0_fp / Tk - 1.e+0_fp ) )
    H%KaqHMS = KHMS * ( Kw1 / Hplus_a ) ! unit is allegedly [s-1]

    ! Now do rate of HMS + OH --> 2SO4 + CH2O - SO2
    KHMS2  = 2.65e+8_fp * EXP( -5.03e+0_fp & ! L/mol/s
         * ( 298.15e+0_fp / Tk - 1.e+0_fp ) )

    dOH = hydroxide/RHO_num ! M cm3 molec -1
    H%KaqHMS2 = KHMS2 * dOH ! [cm3/molec/s]

    IF ( RELHUM >= 50.0_dp ) THEN
       State_Chm%GammaSO2(I,J,L,1) = gammaTMI ! TMI
       State_Chm%GammaSO2(I,J,L,2) = gammaO3 ! O3
       State_Chm%GammaSO2(I,J,L,3) = gammaH2O2 ! H2O2
       State_Chm%GammaSO2(I,J,L,4) = gammaNO2 ! NO2
       State_Chm%GammaSO2(I,J,L,5) = gammaCH2O ! CH2O
    ENDIF

    H%gamma_SO2_TMI   = gammaTMI
    H%gamma_SO2_O3    = gammaO3
    H%gamma_SO2_H2O2  = gammaH2O2
    H%gamma_SO2_NO2   = gammaNO2
    H%gamma_SO2_CH2O  = gammaCH2O
    ! ===========================================================
    
    ! Correction factors for HOBr and HOCl removal by SO2 [1]
    H%fupdateHOBr  = State_Chm%fupdateHOBr(I,J,L)
    H%fupdateHOCl  = State_Chm%fupdateHOCl(I,J,L)

    ! ! mean molecular speed [cm s-1]
    ! sqrt( 8RT / (pi M) )    M_SO2  = MW(ind_SO2) * 1.0e-3_dp
    H%speed = ( SQRT( EIGHT_RSTARG_T / ( H%PI * M_SO2 ) )) *100.0_dp

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

    ! Flag to turn off heterogeneous reactions in stratosphere
    H%TurnOffHetRates = Input_Opt%TurnOffHetRates

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
    IF ((  C(ind_SALACL) + C(ind_NIT) + C(ind_SO4)) > 0.0_dp) THEN
       H%frac_SALACL = C(ind_SALACL) / ( C(ind_SALACL) + C(ind_NIT) + C(ind_SO4) )
    ELSE 
       H%frac_SALACL = 0.0_dp
    ENDIF

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
