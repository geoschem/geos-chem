!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ocean_toolbox_mod
!
! !DESCRIPTION: Module OCEAN\_TOOLBOX\_MOD contains functions and routines to
! calculate the ocean exchange velocity for any gas, according to Johnson, 2010.
!
! M.T. Johnson: "A numerical scheme to calculate temperature and salinity
! dependent air-water transfer velocities for any gas", Ocean Sci. 6, 913-932,
! 2010.
!\\
!\\
! !INTERFACE: 
!
MODULE OCEAN_TOOLBOX_MOD
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: CALC_KG
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: CALC_KA
  PRIVATE :: CALC_KL
  PRIVATE :: N_SW
  PRIVATE :: P_SW
  PRIVATE :: V_SW
  PRIVATE :: D_WC
  PRIVATE :: D_HM
  PRIVATE :: SCHMIDT_W
  PRIVATE :: SCHMIDT_SALTZMANN
  PRIVATE :: SCHMIDT_ACET
  PRIVATE :: N_AIR
  PRIVATE :: P_AIR
  PRIVATE :: V_AIR
  PRIVATE :: D_AIR
  PRIVATE :: SCHMIDT_G
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CALC_KG
!
! !DESCRIPTION: Subroutine CALC\_KG is the wrapper routine to calculate the 
! exchange velocity Kg used for calculating the ocean-air flux.
! F = Kg ( Cg - H * Cl )
! where Cg and Cl are the bulk gas and liquid concentrations and H is the 
! Henry constant (H= Cgs/Cls).
!
! 1/Kg = 1/ka + H/Kl = Ra + Rl.
!
! Note that Kg is returned in m/s and not cm h-1, as is usually reported for
! exchange velocities! 
!
! Reference: Liss and Slater: Flux of gases across the air-sea interface,
! Nature, 247, 1974.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_KG( T, P, V, SALT, H, VB, MW, SCW, KG, RA_OVER_RL, VERBOSE )
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(IN   )  :: T    ! Surface temperature     [C]
    REAL*8,  INTENT(IN   )  :: P    ! Surface pressure        [Pa]
    REAL*8,  INTENT(IN   )  :: V    ! Surface wind speed      [m/s]
    REAL*8,  INTENT(IN   )  :: SALT ! Salinity                [PSU]
    REAL*8,  INTENT(IN   )  :: H    ! Henry constant          [-]
    REAL*8,  INTENT(IN   )  :: VB   ! Liquid mol. volume      [cm3/mol]
    REAL*8,  INTENT(IN   )  :: MW   ! Molecular weight        [g/mol]
    INTEGER, INTENT(IN   )  :: SCW  ! Parameterization type 
                                    ! for Schmidt number in water 
    REAL*8,  INTENT(  OUT)  :: KG   ! Exchange velocity       [ms-1] 
    REAL*8,  INTENT(  OUT), OPTIONAL  :: RA_OVER_RL ! Ra/Rl   [-]
    LOGICAL, INTENT(IN   ), OPTIONAL  :: VERBOSE
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!  21 May 2013 - C. Keller: SCW added to argument list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8                 :: RA,RL
    LOGICAL                :: VERB

    !=================================================================
    ! CALC_KG begins here!
    !=================================================================

    ! Set verbose flag
    IF ( PRESENT ( VERBOSE ) ) THEN
       VERB = VERBOSE
    ELSE
       VERB = .FALSE.
    ENDIF

    ! Calculate air resistence RA
    RA = 1d0 / CALC_KA(T,P,V,MW,VB,VERB)

    ! Calculate water resistence RL
    RL = H / CALC_KL(T,V,SALT,VB,SCW,VERB)

    ! Calculate transfer velocity Kg
    KG = 1d0 / (RA + RL)

    ! Ratio of RA / RL
    IF ( PRESENT(RA_OVER_RL) ) RA_OVER_RL = RA / RL

  END SUBROUTINE CALC_KG
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: CALC_KA
!
! !DESCRIPTION: CALC_KA returns the air exchange velocity KA.
!\\
!\\
! !INTERFACE:
!
  FUNCTION CALC_KA(T,P,V,MW,VB,VERB) RESULT(KA)
!
! !ARGUMENTS:
!
    REAL*8, INTENT(IN)  :: T,P,V, MW, VB!T in C, P in Pa
    LOGICAL, INTENT(IN) :: VERB
    REAL*8              :: KA
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    REAL*8             :: SC, USTAR,CD
    REAL*8, PARAMETER  :: KAPPA=0.4d0

    !=================================================================
    ! CALC_KA begins here!
    !=================================================================

    ! Get Schmidt number in the air
    SC = SCHMIDT_G(T,P,MW,VB)

    !drag coefficient
    CD = 0.61d-3 + 0.063d-3 * V

    !friction velocity
    USTAR = V * SQRT(CD)

    ! Calculate KA
    KA = 1d-3 + USTAR / & 
      ( 13.3d0 * SQRT(SC) + 1d0/SQRT(CD) - 5d0 + LOG(SC)/(2d0*KAPPA) )

    IF ( VERB ) THEN
       WRITE(*,*) 'Schmidt number in air: ', SC
       WRITE(*,*) 'Drag coefficient     : ', CD
       WRITE(*,*) 'Friction velocity    : ', USTAR
       WRITE(*,*) 'Airside resistance   : ', KA
    ENDIF

  END FUNCTION CALC_KA
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: CALC_KL
!
! !DESCRIPTION: CALC_KL calculates the water exchange velocity Kl following
! Nightingale et al., Geophysical Research Letters, 2000.
!\\
!\\
! !INTERFACE:
!
  FUNCTION CALC_KL(T,V,S,VB,SCW,VERB) RESULT(K)
!
! !ARGUMENTS:
!
    REAL*8,  INTENT(IN) :: T,S,V,VB
    INTEGER, INTENT(IN) :: SCW
    LOGICAL, INTENT(IN) :: VERB
    REAL*8              :: K
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8             :: SC,ScCO2

    !=================================================================
    ! CALC_KL begins here!
    !=================================================================

    ! Get Schmidt number in water according to parameterization type
    IF ( SCW == 1 ) THEN
       SC = SCHMIDT_W(T,S,VB)
    ELSEIF ( SCW == 2 ) THEN
       SC = SCHMIDT_SALTZMANN(T)
    ELSEIF ( SCW == 3 ) THEN
       SC = SCHMIDT_ACET(T)
    ENDIF

    ! Schmidt number for CO2
    ScCO2 = 644.7d0 + T * ( -6.16d0 + T * ( 0.11d0 ) )

    ! KL in cm/h according to Nightingale, 2000
    K = V * ( 0.24d0 * V + 0.061d0) / SQRT( SC / ScCO2 )
    
    ! Convert from cm/h to m/s
    K = K / 3600d0 / 100d0

    IF ( VERB ) THEN
       WRITE(*,*) 'Schmidt number in water: ', SC
       WRITE(*,*) 'Schmidt number of CO2  : ', ScCO2
       WRITE(*,*) 'Waterside resistance   : ', K
    ENDIF

  END FUNCTION CALC_KL
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: N_SW
!
! !DESCRIPTION: N_SW returns the dynamic seawater viscosity following 
! Laliberte, 2007.
! Laliberte, M: "Model for calculating the viscosity of aqueous solutions",
! Journal of Chemical & Engineering Data, 52, 2007.
!\\
!\\
! !INTERFACE:
!
  FUNCTION N_SW(T,S) RESULT(N)
!
! !ARGUMENTS
!
    REAL*8, INTENT(IN) :: T,S !temperature (C) and salinity 
    REAL*8             :: N   ! Dynamic viscosity
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8  :: n_0, ln_n_m, w_i_ln_n_i_tot, NI,W_I_TOT,W_I
    INTEGER :: I

    !salt in the order NaCl,KCl,CaCl2,MgCl2,MgSO4
    REAL*8, PARAMETER :: MASS_FRACTION(5) = (/ 0.798D0,0.022D0, &
                                               0.033D0,0.047D0,0.1D0 /)

    REAL*8, PARAMETER :: V1(5) = (/ 16.22D0,      6.4883D0, &
                                    32.028D0,    24.032D0, 72.269D0/)

    REAL*8, PARAMETER :: V2(5) = (/  1.3229D0 ,   1.3175D0, &
                                    0.78792D0,    2.2694D0, 2.2238D0/)

    REAL*8, PARAMETER :: V3(5) = (/  1.4849D0 ,  -0.7785D0, &
                                    -1.1495D0,    3.7108D0, 6.6037D0/)

    REAL*8, PARAMETER :: V4(5) = (/  0.0074691D0, 0.09272D0, &
                                     0.0026995D0, 0.021853D0, 0.0079004D0/)

    REAL*8, PARAMETER :: V5(5) = (/ 30.78D0 ,    -1.3D0, &
                                    780860D0,    -1.1236D0, 3340.1D0/)

    REAL*8, PARAMETER :: V6(5) = (/ 2.0583D0,     2.0811D0, &
                                    5.8442D0,     0.14474D0, 6.1304D0/)

    !=================================================================
    ! N_SW begins here!
    !=================================================================

    ! Init
    W_I_TOT = 0d0
    w_i_ln_n_i_tot = 0d0

    DO I=1,5
       W_I = MASS_FRACTION(I) * S / 1000D0
       W_I_TOT = W_I + W_I_TOT
    ENDDO

    DO I=1,5
       W_I = MASS_FRACTION(I) * S / 1000d0       
       NI = exp( ((v1(I)*w_i_tot**v2(I))+v3(I)) / ((v4(I)*T)+1) ) / &
            ( (v5(I)*(w_i_tot**v6(I))) + 1d0 )
       w_i_ln_n_i_tot = w_i_ln_n_i_tot + (W_I * log(NI))
    ENDDO

    n_0 = (T+246d0)/(137.37d0+(5.2842d0*T)+(0.05594d0*(T**2d0)))

    ln_n_m = (1d0 - w_i_tot) * log(n_0) + w_i_ln_n_i_tot

    N = exp(ln_n_m)

  END FUNCTION N_SW
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: P_SW
!
! !DESCRIPTION: P_SW returns the seawater density following Millero and Poisson,
! 1981. 
! Millero and Poisson: "International one-atmosphere equation of state of
! seawater", Deep Sea Res. Pt A, 28, 1981.
!\\
!\\
! !INTERFACE:
!
  FUNCTION P_SW(T,S) RESULT(P)
!
! !ARGUMENTS
!
    REAL*8, INTENT(IN) :: T,S
    REAL*8             :: P
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8            :: A, B, C

    !=================================================================
    ! P_SW begins here!
    !=================================================================

    !A = 0.824493D0-(4.0899D-3*T)+(7.6438D-5*(T**2d0))-(8.2467D-7*(T**3d0))+(5.3875D-9*(T**4d0))
!    B = -5.72466D-3+(1.0277D-4*T)-(1.6546D-6*(T**2d0))
    A =  0.824493d0 + T * (-4.0899d-3 + T * (7.6438d-5 + T * (-8.2467d-7 + T * 5.3875d-9)))
    B = -5.72466D-3 + T * ( 1.0277D-4 + T * -1.6546D-6)
    C =  4.8314D-4

    ! Density of pure water
!    P = 999.842594D0+(6.793952D-2*T)-(9.09529D-3*(T**2d0))+(1.001685D-4*(T**3d0))-(1.120083D-6*(T**4d0))+(6.536332D-9*(T**5d0))
    P = 999.842594D0 + T * (6.793952D-2 + T * (-9.09529D-3 + T * (1.001685D-4 + T * (-1.120083D-6 + T * 6.536332D-9 ))))

    ! Salinity correction
    P = (P+(A*S)+(B*(S**(1.5D0)))+(C*S))

  END FUNCTION P_SW
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: V_SW
!
! !DESCRIPTION: V_SW returns the ???
!\\
!\\
! !INTERFACE:
!
  FUNCTION V_SW(T,S) RESULT(V)
!
! !ARGUMENTS
!
    REAL*8, INTENT(IN) :: T,S
    REAL*8             :: V
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8   :: N, P

    !=================================================================
    ! V_SW begins here!
    !=================================================================

    N=N_SW(T,S)*1D-3
    P=P_SW(T,S)

    V = 1D4*N/P

  END FUNCTION  V_SW
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: D_WC
!
! !DESCRIPTION: D_WC returns the (water) diffusion coefficient following Wilke
!  and Chang, 1955.
! Wilke and Chang: "Correlation of diffusion coefficients in dilute solutions",
! AIChE Journal, 1, 1955.
!\\
!\\
! !INTERFACE:
!
  FUNCTION D_WC(T,S,VB) RESULT(D)
!
! !ARGUMENTS
!
    REAL*8, INTENT(IN) :: T,S,VB
    REAL*8             :: D
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8, PARAMETER  :: PHI = 2.6D0

    !=================================================================
    ! D_WC begins here!
    !=================================================================

    D = ((T+273.15D0)*7.4D-8*(PHI*18.01D0)**0.5D0)/((N_SW(T,S))*(VB**0.6D0))

  END FUNCTION D_WC
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: D_HM
!
! !DESCRIPTION: D_HM returns the (water) diffusivity following Hayduk and
! Minhas, 1982.
! Hayduk and Minhas, "Correlations for prediction of molecular diffusivities in
! liquids, Can. J. Chem. Eng., 60, 1982,
!\\
!\\
! !INTERFACE:
!
  FUNCTION D_HM(T,S,VB) RESULT(D)
!
! !ARGUMENTS
!
    REAL*8, INTENT(IN) :: T,S,VB
    REAL*8             :: D
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8     :: EpsilonStar

    !=================================================================
    ! D_HM begins here!
    !=================================================================

    EpsilonStar = (9.58D0/VB)-1.12D0
    D=1.25D-8*(VB**(-0.19D0)-0.292D0)*((T+273.15D0)**(1.52D0))*((N_SW(T,S))**EpsilonStar)

  END FUNCTION D_HM
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: SCHMIDT_W
!
! !DESCRIPTION: SCHMIDT_W returns the Schmidt number of the gas in the water
! following Johnson, 2010. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION SCHMIDT_W(T,S,VB) RESULT(SC)
!
! !ARGUMENTS
!
    REAL*8, INTENT(IN) :: T,S,VB
    REAL*8             :: SC
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! SCHMIDT_W begins here!
    !=================================================================

    SC = 2D0 * V_SW(T,S) / ( D_HM(T,S,VB) + D_WC(T,S,VB) )

  END FUNCTION SCHMIDT_W
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: SCHMIDT_SALTZMANN
!
! !DESCRIPTION: SCHMIDT_SALTZMANN returns the Schmidt number of the gas in 
! the water calculated according to Saltzmann et al., 1993.
!\\
!\\
! !INTERFACE:
!
  FUNCTION SCHMIDT_SALTZMANN(T) RESULT(SC)
!
! !ARGUMENTS
!
    REAL*8, INTENT(IN) :: T   ! Temperature in C
    REAL*8             :: SC
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! SCHMIDT_SALTZMANN begins here!
    !=================================================================

    SC = 2674.0d0 + T * ( -147.12d0 + T * ( 3.726d0 - T * 0.038d0 ) )  

  END FUNCTION SCHMIDT_SALTZMANN
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: SCHMIDT_ACET
!
! !DESCRIPTION: SCHMIDT_ACET returns the Schmidt number of acetone. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION SCHMIDT_ACET(T) RESULT(SC)
!
! !ARGUMENTS
!
    REAL*8, INTENT(IN) :: T   ! Temperature in C
    REAL*8             :: SC
!
! !REVISION HISTORY:
!  11 Aug 2013 - C. Keller: Initial version 
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! SCHMIDT_ACET begins here!
    !=================================================================

    SC = 3287.687d0 + T * ( -136.2176d0 + T * ( 2.20642d0 - T*0.01410642d0 ) )  

  END FUNCTION SCHMIDT_ACET
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: N_AIR
!
! !DESCRIPTION: N_AIR returns the dynamic air viscosity. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION N_AIR(T) RESULT(N)
!
! !ARGUMENTS
!
    REAL*8, INTENT(IN) :: T
    REAL*8             :: N
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    REAL*8             :: SV_0,SV_1,SV_2,SV_3,SV_4

    !=================================================================
    ! N_AIR begins here!
    !=================================================================

    SV_0 =  1.715747771D-5
    SV_1 =  4.722402075D-8
    SV_2 = -3.663027156D-10
    SV_3 =  1.873236686D-12
    SV_4 = -8.050218737D-14

    ! in N.s/m^2 (Pa.s)
!    N = SV_0+(SV_1*T)+(SV_2*T**2d0)+(SV_3*T**3d0)+(SV_4*T**4d0)
    N = SV_0 + T * ( SV_1 + T * ( SV_2 + T * (SV_3 + T * SV_4) ) )

  END FUNCTION N_AIR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: P_AIR
!
! !DESCRIPTION: P_AIR returns the kinematic air viscosity. 
!\\
!\\
! !INTERFACE:
!
  FUNCTION P_AIR(T) RESULT(P)
!
! !ARGUMENTS
!
    REAL*8, INTENT(IN) :: T
    REAL*8             :: P
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    REAL*8             :: SD_0,SD_1,SD_2,SD_3

    !=================================================================
    ! P_AIR begins here!
    !=================================================================

    SD_0 =  1.293393662D0
    SD_1 = -5.538444326D-3
    SD_2 =  3.860201577D-5
    SD_3 = -5.2536065D-7
    P = SD_0 + T * (SD_1 + T * (SD_2 + T * SD_3 ))
!    P = SD_0+(SD_1*T)+(SD_2*T**2)+(SD_3*T**3)

  END FUNCTION P_AIR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: V_AIR
!
! !DESCRIPTION: V_AIR returns the kinematic air viscosity (m2/s).
!\\
!\\
! !INTERFACE:
!
  FUNCTION V_AIR(T) RESULT(V)
!
! !ARGUMENTS
!
    REAL*8, INTENT(IN) :: T
    REAL*8             :: V
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! V_AIR begins here!
    !=================================================================

    V = N_AIR(T)/P_AIR(T)

  END FUNCTION V_AIR
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: D_AIR
!
! !DESCRIPTION: D_AIR returns the gas phase diffusion coefficient according to
! Fuller et al., 1966.
! E. Fuller et al.: "New method for prediction of binary gas-phase diffusion
! coefficients", Industrial & Engineering Chemistry, 58, 1966.
!\\
!\\
! !INTERFACE:
!
  FUNCTION D_AIR(T,P,MW,VB) RESULT(D)
!
! !ARGUMENTS
!
    REAL*8, INTENT(IN) :: T  !T in C
    REAL*8, INTENT(IN) :: P  !P in Pa
    REAL*8, INTENT(IN) :: MW !MW in g/mol
    REAL*8, INTENT(IN) :: VB !Liq. molar volume (cm3/mol)
    REAL*8             :: D
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8            :: Pa
    REAL*8, PARAMETER :: Ma = 29.97d0 !Mw air in g/mol
    REAL*8, PARAMETER :: Va = 20.1d0  !cm3/mol (diffusion volume for air)
    
    !=================================================================
    ! D_AIR begins here!
    !=================================================================

    ! Convert P to atm
    PA = 9.8692D-6*P

    ! Calculate diffusion coefficient
    D  = 1D-3 * (T+273.15D0)**(1.75D0)*SQRT(1D0/Ma +     & 
        1D0/MW)/(PA*(VA**(1D0/3D0)+VB**(1D0/3D0))**2D0)

    !D is in cm2/s convert to m2/s
    D = D * 1D-4

  END FUNCTION D_AIR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !FUNCTION: SCHMIDT_G
!
! !DESCRIPTION: SCHMIDT_G returns the schmidt number of the gas in the air.
!\\
!\\
! !INTERFACE:
!
  FUNCTION SCHMIDT_G(T,P,MW,VB) RESULT(SC)
!
! !ARGUMENTS:
!
    REAL*8, INTENT(IN) :: T, P, MW, VB
    REAL*8             :: SC
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8             :: D,V

    !=================================================================
    ! SCHMIDT_G begins here!
    !=================================================================

    V=V_AIR(T)
    D=D_AIR(T,P,MW,VB)

    SC = V / D

  END FUNCTION SCHMIDT_G
!EOC
END MODULE OCEAN_TOOLBOX_MOD
!EOM
