!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ocean_toolbox_mod.F90
!
! !DESCRIPTION: Module Ocean\_ToolBox\_Mod contains functions and routines to
! calculate the ocean exchange velocity for any gas, according to Johnson, 2010.
!
! References:
! \begin{itemize}
! \item  M.T. Johnson: "A numerical scheme to calculate temperature and
!        salinity dependent air-water transfer velocities for any gas",
!        Ocean Sci. 6, 913-932, 2010.
! \item  Liss and Slater: Flux of gases across the air-sea interface,
!        Nature, 247, 1974.
! \item  Laliberte, M: "Model for calculating the viscosity of aqueous
!        solutions", Journal of Chemical \& Engineering Data, 52, 2007.
! \item  Millero and Poisson: "International one-atmosphere equation of
!        state of seawater", Deep Sea Res. Pt A, 28, 1981.
! \item  Wilke and Chang: "Correlation of diffusion coefficients in dilute
!        solutions", AIChE Journal, 1, 1955.
! \item  Hayduk and Minhas, "Correlations for prediction of molecular
!        diffusivities in liquids, Can. J. Chem. Eng., 60, 1982,
! \item  E. Fuller et al.: "New method for prediction of binary gas-phase
!        diffusion coefficients", Industrial \& Engineering Chemistry, 58,
!        1966.
! \item Saltzman et al.: Experimental determination of the diffusion
!    coefficient of dimethylsulfide in water, J. Geophys. Res., 98, 1993.
! \end{itemize}
!
! !INTERFACE:
!
MODULE Ocean_ToolBox_Mod
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: Calc_Kg
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Calc_Ka
  PRIVATE :: Calc_Kl
  PRIVATE :: N_SW
  PRIVATE :: P_SW
  PRIVATE :: V_SW
  PRIVATE :: D_WC
  PRIVATE :: D_HM
  PRIVATE :: Schmidt_W
  PRIVATE :: Schmidt_Saltzmann
  PRIVATE :: Schmidt_Acet
  PRIVATE :: Schmidt_Ald2
  PRIVATE :: N_Air
  PRIVATE :: P_Air
  PRIVATE :: V_Air
  PRIVATE :: D_Air
  PRIVATE :: Schmidt_G
!
! ! PARAMETER
!
  INTEGER, PARAMETER :: OC_SUCCESS = 0
  INTEGER, PARAMETER :: OC_FAIL    = -999
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!  03 Oct 2-14 - C. Keller: Added error trap for negative Schmidt numbers
!  10 Mar 2017 - M. Sulprizio- Add function Schmidt_Ald2
!EOP
!------------------------------------------------------------------------------
!BOC
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Kg
!
! !DESCRIPTION: Subroutine Calc\_Kg is the wrapper routine to calculate the
! exchange velocity Kg used for calculating the ocean-air flux (cf. Liss \&
! Slater 1974) as:
!\\
!\\
! F = Kg ( Cg - H * Cl )
!\\
!\\
! where Cg and Cl are the bulk gas and liquid concentrations and H is the
! Henry constant (H= Cgs/Cls).
!\\
!\\
! 1/Kg = 1/ka + H/Kl = Ra + Rl.
!\\
!\\
! Note that Kg is returned in m/s and not cm h-1, as is usually reported for
! exchange velocities!
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_Kg( T, P, V, SALT, H, VB, MW, SCW, KG, RC, RA_OVER_RL, VERBOSE )
!
! !INPUT PARAMETERS:
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
    LOGICAL, INTENT(IN   ), OPTIONAL  :: VERBOSE   ! turn on verbose output
!
! !OUTPUT PARAMETERS:
!
    REAL*8,  INTENT(  OUT)  :: KG   ! Exchange velocity       [ms-1]
    INTEGER, INTENT(  OUT)  :: RC   ! Error code
    REAL*8,  INTENT(  OUT), OPTIONAL  :: RA_OVER_RL ! Ra/Rl   [-]
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller: Adapted from F. Paulot
!  21 May 2013 - C. Keller: SCW added to argument list
!  15 Aug 2014 - C. Keller: Now limit temperature to -40 degC to avoid overflow
!                           error. Also added error trap for temperatures
!                           between -10.7 and -10.9 degrees that cause a div-zero
!                           error in subroutine N_SW.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL :: VERB
    REAL*8  :: RA, RL, TMP, KL, KA
    REAL*8, PARAMETER :: TMAX = -40.0d0

    !=================================================================
    ! CALC_KG begins here!
    !=================================================================

    ! Fail by default
    RC = OC_FAIL

    ! Set verbose flag
    IF ( PRESENT ( VERBOSE ) ) THEN
       VERB = VERBOSE
    ELSE
       VERB = .FALSE.
    ENDIF

    ! Set KG to zero and return if winds are 0
    IF ( V == 0.d0 ) THEN
       KG = 0d0
       RC = OC_SUCCESS
       RETURN
    ENDIF

    ! Surface temperature must be greater than -40 degC. Otherwise, an
    ! overflow error may occur!
    TMP = MAX(T,TMAX)

    ! Don't allow values between -10.7 and -10.9 to avoid div-zero error
    ! in N_SW!
    IF ( TMP > -10.7d0 .AND. TMP < -10.9d0 ) THEN
       TMP = -10.7d0
    ENDIF

    ! Calculate air resistence RA
    KA = CALC_KA(TMP,P,V,MW,VB,VERB)
    IF ( KA < 0.0d0 ) RETURN
    RA = 1d0 / KA

    ! Calculate water resistence RL
    KL = CALC_KL(TMP,V,SALT,VB,SCW,VERB)
    IF ( KL < 0.0d0 ) RETURN
    RL = H / KL

    ! Calculate transfer velocity Kg
    KG = 1d0 / (RA + RL)

    ! Ratio of RA / RL
    IF ( PRESENT(RA_OVER_RL) ) RA_OVER_RL = RA / RL

    ! Return w/ success
    RC = OC_SUCCESS

  END SUBROUTINE Calc_Kg
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Ka
!
! !DESCRIPTION: Calc\_Ka returns the air exchange velocity KA.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Calc_Ka(T,P,V,MW,VB,VERB) RESULT(KA)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN)  :: T,P,V, MW, VB!T in C, P in Pa
    LOGICAL, INTENT(IN) :: VERB
!
! !RETURN VALUE:
!
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

  END FUNCTION Calc_Ka
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Kl
!
! !DESCRIPTION: Calc\_Kl calculates the water exchange velocity Kl following
! Nightingale et al., Geophysical Research Letters, 2000.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Calc_Kl(T,V,S,VB,SCW,VERB) RESULT(K)
!
! !INPUT PARAMETERS:
!
    REAL*8,  INTENT(IN) :: T,S,V,VB
    INTEGER, INTENT(IN) :: SCW
    LOGICAL, INTENT(IN) :: VERB
!
! !RETURN VALUE:
!
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
    ELSEIF ( SCW == 4 ) THEN
       SC = SCHMIDT_ALD2(T)
    ENDIF

    ! Schmidt number for CO2
    ScCO2 = 644.7d0 + T * ( -6.16d0 + T * ( 0.11d0 ) )

    ! Error trap: Schmidt numbers MUST be positive
    IF ( SC < 0.0d0 .OR. ScCO2 < 0.0d0 ) THEN
       WRITE(*,*) 'Negative Schmidt number!'
       WRITE(*,*) 'SC, ScCO2, T, S, VB: ', SC, ScCO2, T, S, VB
       K = -999.0d0
       RETURN
    ENDIF

    ! KL in cm/h according to Nightingale, 2000
    K = V * ( 0.24d0 * V + 0.061d0) / SQRT( SC / ScCO2 )

    ! Convert from cm/h to m/s
    K = K / 3600d0 / 100d0

    IF ( VERB ) THEN
       WRITE(*,*) 'Schmidt number in water: ', SC
       WRITE(*,*) 'Schmidt number of CO2  : ', ScCO2
       WRITE(*,*) 'Waterside resistance   : ', K
    ENDIF

  END FUNCTION Calc_Kl
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: N_SW
!
! !DESCRIPTION: N\_SW returns the dynamic seawater viscosity following
! Laliberte, 2007.
!\\
!\\
! !INTERFACE:
!
  FUNCTION N_SW(T,S) RESULT(N)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: T,S !temperature (C) and salinity
!
! !RETURN VALUE:
!
    REAL*8             :: N   ! Dynamic viscosity
!
! !REVISION HISTORY:
!  11 Apr 2013 - C. Keller   - Adapted from F. Paulot
!   1 Jul 2014 - R. Yantosca - Bug fix: Don't take LOG(NI) if NI is zero
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8  :: n_0, ln_n_m, w_i_ln_n_i_tot, NI, W_I_TOT, W_I
    INTEGER :: I

    !salt in the order NaCl, KCl, CaCl2, MgCl2, MgSO4
    REAL*8, PARAMETER :: MASS_FRACTION(5) = &
       (/  0.798D0,     0.022D0,    0.033D0,     0.047D0,    0.1D0       /)

    REAL*8, PARAMETER :: V1(5) = &
       (/ 16.22D0,      6.4883D0,  32.028D0,    24.032D0,   72.269D0     /)

    REAL*8, PARAMETER :: V2(5) = &
       (/  1.3229D0,    1.3175D0,   0.78792D0,   2.2694D0,   2.2238D0    /)

    REAL*8, PARAMETER :: V3(5) = &
       (/  1.4849D0,   -0.7785D0,  -1.1495D0,    3.7108D0,   6.6037D0    /)

    REAL*8, PARAMETER :: V4(5) = &
       (/  0.0074691D0, 0.09272D0,  0.0026995D0, 0.021853D0, 0.0079004D0 /)

    REAL*8, PARAMETER :: V5(5) = &
       (/ 30.78D0,     -1.3D0,      780860D0,   -1.1236D0,   3340.1D0    /)

    REAL*8, PARAMETER :: V6(5) = &
       (/  2.0583D0,    2.0811D0,   5.8442D0,    0.14474D0,  6.1304D0    /)

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

       NI = exp(                                          &
                 ( ( v1(I) * w_i_tot**v2(I) ) + v3(I) ) / &
                 ( ( v4(I) * T ) + 1d0                )   &
               ) / (                                      &
                     ( v5(I) * (w_i_tot**v6(I)) ) + 1d0   &
                   )

       !-------------------------------------------------------------
       ! Prior to 7/1/14:
       ! NOTE: PRESERVE ORIGINAL CODE HERE:
       !w_i_ln_n_i_tot = w_i_ln_n_i_tot + (W_I * log(NI))
       !-------------------------------------------------------------

       ! BUG FIX: Don't take the log unless NI > 0!! (bmy, 7/1/14)
       IF ( NI > 0d0 ) THEN
          w_i_ln_n_i_tot = w_i_ln_n_i_tot + (W_I * log(NI))
       ENDIF
    ENDDO

    n_0 = (T+246d0)/(137.37d0+(5.2842d0*T)+(0.05594d0*(T**2d0)))

    ln_n_m = (1d0 - w_i_tot) * log(n_0) + w_i_ln_n_i_tot

    N = exp(ln_n_m)

  END FUNCTION N_SW
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: P_SW
!
! !DESCRIPTION: P\_SW returns the seawater density following Millero and
! Poisson,  1981.
!\\
!\\
! !INTERFACE:
!
  FUNCTION P_SW(T,S) RESULT(P)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: T,S
!
! !RETURN VALUE:
!
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
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: V_SW
!
! !DESCRIPTION: V\_SW returns the ???
!\\
!\\
! !INTERFACE:
!
  FUNCTION V_SW(T,S) RESULT(V)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: T,S
!
! !RETURN VALUE:
!
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
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: D_WC
!
! !DESCRIPTION: D\_WC returns the (water) diffusion coefficient following Wilke
!  and Chang, 1955.
!\\
!\\
! !INTERFACE:
!
  FUNCTION D_WC(T,S,VB) RESULT(D)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: T,S,VB
!
! !RETURN VALUE:
!
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
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: D_HM
!
! !DESCRIPTION: D\_HM returns the (water) diffusivity following Hayduk and
! Minhas, 1982.
!\\
!\\
! !INTERFACE:
!
  FUNCTION D_HM(T,S,VB) RESULT(D)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: T,S,VB
!
! !RETURN VALUE:
!
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
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Schmidt_W
!
! !DESCRIPTION: Schmidt\_W returns the Schmidt number of the gas in the water
! following Johnson, 2010.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Schmidt_W(T,S,VB) RESULT(SC)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: T,S,VB
!
! !RETURN VALUE:
!
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

  END FUNCTION Schmidt_W
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Schmidt_Saltzmann
!
! !DESCRIPTION: Schmidt\_Saltzmann returns the Schmidt number of the gas in
! the water calculated according to Saltzmann et al., 1993.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Schmidt_Saltzmann(T) RESULT(SC)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: T   ! Temperature in C
!
! !RETURN VALUE:
!
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

  END FUNCTION Schmidt_Saltzmann
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Schmidt_Acet
!
! !DESCRIPTION: Schmidt\_Acet returns the Schmidt number of acetone.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Schmidt_Acet(T) RESULT(SC)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: T   ! Temperature in C
!
! !RETURN VALUE:
!
    REAL*8             :: SC
!
! !NOTES:
!  Coefficients for fitting the Schmidt number for acetone [unitless]
!    A0 =  3287.687d0
!    A1 = -136.2176d0
!    A2 =  2.20642d0
!    A3 = -0.01410642d0
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

  END FUNCTION Schmidt_Acet
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Schmidt_Ald2
!
! !DESCRIPTION: Schmidt\_Ald2 returns the Schmidt number of acetaldehyde.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Schmidt_Ald2(T) RESULT(SC)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: T   ! Temperature in C
!
! !RETURN VALUE:
!
    REAL*8             :: SC
!
! !NOTES:
!  Coefficients for fitting the Schmidt number for acetaldehyde [unitless]
!  Derived using polynomial fit (code provided by qli, same as used
!  for acetone, methanol)
!  and partial molal volume of acetaldehyde at its normal boiling
!  temperature (51.8 cm3/g/mole) calculated using Le Bas method
!  see "The Properties of Gases and Liquids", Reid, Prausnitz, Sherwood.
!    A0 = 2581.709d0
!    A1 = -106.9671d0
!    A2 = 1.73263d0
!    A3 = -0.0110773d0
!
! !REVISION HISTORY:
!  10 Mar 2017 - M. Sulprizio- Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! SCHMIDT_ALD2 begins here!
    !=================================================================

    SC = 2581.709d0 + T * ( -106.9671d0 + T * ( 1.73263d0 - T*0.0110773d0 ) )

  END FUNCTION Schmidt_Ald2
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: N_Air
!
! !DESCRIPTION: N\_Air returns the dynamic air viscosity.
!\\
!\\
! !INTERFACE:
!
  FUNCTION N_Air(T) RESULT(N)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: T
!
! !RETURN VALUE:
!
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

  END FUNCTION N_Air
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: P_Air
!
! !DESCRIPTION: P\_Air returns the kinematic air viscosity.
!\\
!\\
! !INTERFACE:
!
  FUNCTION P_Air(T) RESULT(P)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: T
!
! !RETURN VALUE:
!
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

  END FUNCTION P_Air
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: V_Air
!
! !DESCRIPTION: V\_Air returns the kinematic air viscosity (m2/s).
!\\
!\\
! !INTERFACE:
!
  FUNCTION V_Air(T) RESULT(V)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: T
!
! !RETURN VALUE:
!
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

  END FUNCTION V_Air
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: D_Air
!
! !DESCRIPTION: D\_Air returns the gas phase diffusion coefficient according to
! Fuller et al., 1966.
!\\
!\\
! !INTERFACE:
!
  FUNCTION D_Air(T,P,MW,VB) RESULT(D)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: T  !T in C
    REAL*8, INTENT(IN) :: P  !P in Pa
    REAL*8, INTENT(IN) :: MW !MW in g/mol
    REAL*8, INTENT(IN) :: VB !Liq. molar volume (cm3/mol)
!
! !RETURN VALUE:
!
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

  END FUNCTION D_Air
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Schmidt_G
!
! !DESCRIPTION: Schmidt\_G returns the schmidt number of the gas in the air.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Schmidt_G(T,P,MW,VB) RESULT(SC)
!
! !INPUT PARAMETERS:
!
    REAL*8, INTENT(IN) :: T, P, MW, VB
!
! !RETURN VALUE:
!
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
    REAL*8 :: D,V

    !=================================================================
    ! SCHMIDT_G begins here!
    !=================================================================

    V=V_AIR(T)
    D=D_AIR(T,P,MW,VB)

    SC = V / D

  END FUNCTION Schmidt_G
!EOC
END MODULE Ocean_ToolBox_Mod
