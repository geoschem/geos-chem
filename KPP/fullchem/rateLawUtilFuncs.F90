!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: rateLawUtilFuncs
!
! !DESCRIPTION: Provides common functions for computing reaction rates.
!\\
!\\
! !INTERFACE:
!
MODULE rateLawUtilFuncs
!
! !USES:
!
  USE gckpp_Global
  USE gckpp_Precision

  IMPLICIT NONE
  PUBLIC
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: coth
!
! !DEFINED PARAMETERS:
!
  ! Minimum heterogeneous chemistry lifetime and reaction rate
  REAL(dp), PRIVATE, PARAMETER :: HET_MIN_LIFE = 1.e-3_dp
  REAL(dp), PRIVATE, PARAMETER :: HET_MIN_RATE = 1.0_dp / HET_MIN_LIFE
!EOP
!-----------------------------------------------------------------------------
!BOC
CONTAINS

  !#########################################################################
  !#####                   ARRHENIUS FUNCTIONS                         #####
  !#########################################################################

  FUNCTION GCARR_ab( a0, b0 ) RESULT( k )
    ! Arrhenius function, skipping computation of EXP( c0/T ),
    ! which evaluates to 1 when c0=0.  This avoids excess CPU
    ! cycles. (bmy, 12/18/20)
    !
    REAL(dp), INTENT(IN) :: a0, b0
    REAL(dp)             :: k
    !
    k = a0 * K300_OVER_TEMP**b0
  END FUNCTION GCARR_ab

  FUNCTION GCARR_ac( a0, c0 ) RESULT( k )
    ! Arrhenius function, skipping computation of ( 300/T )**b0,
    ! which evaluates to 1 when b0=0.  This avoids excess CPU
    ! cycles (bmy, 12/18/20)
    !
    REAL(dp), INTENT(IN) :: a0, c0
    REAL(dp)             :: k
    !
    k = a0 * EXP( c0 / TEMP )
  END FUNCTION GCARR_ac

  FUNCTION GCARR_abc( a0, b0, c0 ) RESULT( k )
    ! Arrhenius function, using all 3 terms.
    ! Use this when a0, b0, c0 are all nonzero.
    !
    REAL(dp), INTENT(IN) :: a0, b0, c0
    REAL(dp)             :: k
    !
    k = a0 * EXP( c0 / TEMP ) * K300_OVER_TEMP**b0
  END FUNCTION GCARR_abc

  !#########################################################################
  !#####         COMMON FUNCTIONS FOR COMPUTING UPTAKE RATES           #####
  !#########################################################################

  FUNCTION Ars_L1k( area, radius, gamma, srMw ) RESULT( k )
    !
    ! Calculates the 1st-order loss rate of species on wet aerosol surface.
    !
    REAL(dp), INTENT(IN) :: area, radius, gamma, srMw
    REAL(dp)             :: k,    dfkg
    !
    ! If gamma or radius is very small, set rate to zero and return
    IF ( gamma < 1.0e-30_dp .or. radius < 1.0e-30_dp ) THEN
       k = 0.0_dp
       RETURN
    ENDIF
    !
    ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
    dfkg = ( 9.45E+17_dp / NUMDEN ) * SR_TEMP *                              &
           SQRT( 3.472E-2_dp + 1.0_dp / ( srMw * srMw ) )
    !
    ! Compute ArsL1k according to the formula listed above
    k = area / ( (radius / dfkg) + 2.749064E-4_dp * srMw / (gamma * SR_TEMP) )
  END FUNCTION Ars_L1k

  FUNCTION kIIR1Ltd( concGas, concEduct, kISource ) RESULT( kII )
    !
    ! Determine removal rates for both species in an uptake reaction.
    ! - Assume that the 1st reactant (concGas) is limiting.
    ! - Assume that the 2nd reactant (concEduct) is "abundant".
    ! - Calculate the overall rate (kII) based only on the uptake
    !   rate of the first reactant.
    ! NOTE: Rewritten for computational efficiency (bmy, 5/13/21)
    !
    REAL(dp), INTENT(IN) :: concGas, concEduct, kISource
    REAL(dp)             :: kIGas,   kIEduct,   lifeA,   lifeB, kII
    !
    kIGas   = 0.0_dp
    kIEduct = 0.0_dp
    kII     = 0.0_dp
    !
    ! Prevent div by zero.   NOTE: Now use 1.0e-8 as the error trap for
    ! concEduct.  100 (the previous criterion) was too large.  We should
    ! never have have a concentration as low as 1e-8 molec/cm3, as that
    ! will definitely blow up the division. (bmy, 11/1/21)
    IF ( concEduct < 1.0e-8_dp                           ) RETURN
    IF ( .not. Is_SafeDiv( concGas*kISource, concEduct ) ) RETURN
    !
    ! Compute rates
    kIGas   = kISource
    kIEduct = kIGas * concGas / concEduct
    kII     = kIGas           / concEduct
    !
    ! Enforce a minimum lifetime?
    IF ( kIGas > 0.0_dp ) THEN
       !
       ! Calculate lifetime of each reactant against removal
       lifeA = SafeDiv( 1.0_dp, kIGas,   0.0_dp )
       lifeB = SafeDiv( 1.0_dp, kIEduct, 0.0_dp )
       !
       ! Check if either lifetime is "too short"
       IF ( ( lifeA < lifeB ) .and. ( lifeA < HET_MIN_LIFE ) ) THEN
          kII = SafeDiv( HET_MIN_RATE, concEduct, 0.0_dp )
       ELSE IF ( lifeB < HET_MIN_LIFE ) THEN
          kII = SafeDiv( HET_MIN_RATE, concGas, 0.0_dp )
       ENDIF
    ENDIF
  END FUNCTION kIIR1Ltd

  FUNCTION CloudHet( H, srMw, gamLiq, gamIce, brLiq, brIce ) RESULT( kHet )
    !
    ! Function CloudHet calculates the loss frequency (1/s) of gas species
    ! due to heterogeneous chemistry on clouds in a partially cloudy grid
    ! cell. The function uses the "entrainment limited uptake" equations of
    ! Holmes et al. (2019). Both liquid and ice water clouds are treated.
    !
    ! For gasses that are that are consumed in multiple aqueous reactions
    ! with different products, CloudHet can provide the loss frequency for
    ! each reaction branch using the branch ratios (branchLiq, branchIce).
    !
    ! Reference:
    ! Holmes, C.D., Bertram, T. H., Confer, K. L., Ronan, A. C., Wirks,
    !   C. K., Graham, K. A., Shah, V. (2019) The role of clouds in the
    !   tropospheric NOx cycle: a new modeling approach for cloud chemistry
    !   and its global implications, Geophys. Res. Lett. 46, 4980-4990,
    !   https://doi.org/10.1029/2019GL081990
    !
    TYPE(HetState), INTENT(IN) :: H              ! Hetchem State object
    REAL(dp),       INTENT(IN) :: srMw           ! SQRT( mol wt [g/mole] )
    REAL(dp),       INTENT(IN) :: gamLiq         ! Rxn prob, liquid [1]
    REAL(dp),       INTENT(IN) :: gamIce         ! Rxn prob, ice [1]
    REAL(dp),       INTENT(IN) :: brLiq          ! Frac of reactant consumed
    REAL(dp),       INTENT(IN) :: brIce          !  in liq & ice branches [0-1]
    REAL(dp)                   :: kHet           ! Grid-avg loss frequency [1/s]
    !
    REAL(dp),       PARAMETER  :: tauc = 3600.0_dp
    REAL(dp)                   :: kI, gam, rd, area
    REAL(dp)                   :: kk, ff, xx, branch, kIb, ktmp
    LOGICAL                    :: isCloud

    ! If cloud fraction < 0.0001 (0.01%) or there is zero cloud surface
    ! area, then return zero uptake
    IF ( ( H%CldFr < 0.0001_dp ) .or. ( H%aLiq + H%aIce <= 0.0_dp ) ) THEN
       kHet = 0.0_dp
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Loss frequency inside cloud
    !
    ! Assume both water and ice phases are inside the same cloud, so mass
    ! transport to both phases works in parallel (additive)
    !-----------------------------------------------------------------------

    ! initialize
    kI   = 0.0_dp
    kIb  = 0.0_dp
    ktmp = 0.0_dp
    kHet = 0.0_dp

    !-----------------------------------------------------------------------
    ! Liquid branch (skip if the liquid branching ratio is zero)
    !-----------------------------------------------------------------------
    IF ( brLiq > 0.0_dp ) THEN

       ! Convert grid-average cloud condensate surface area density
       ! to in-cloud surface area density
       area = SafeDiv( H%aLiq, H%CldFr, 0.0_dp )

       ! Skip if no area
       IF ( area > 0.0_dp ) THEN

          ! In-cloud loss frequency [1/s]
          ktmp = Ars_L1K( area, H%rLiq, gamLiq, srMw )
          kI   = kI + ktmp

          ! In-cloud loss frequency for liquid rxn branch [1/s]
          kIb  = kIb + ( ktmp * brLiq )
       ENDIF
    ENDIF

    !------------------------------------------------------------------
    ! Ice branch (skip if the ice branching ratio is zero)
    !------------------------------------------------------------------
    IF ( brIce > 0.0_dp ) THEN

       ! Convert grid-average cloud condensate surface area density
       ! to in-cloud surface area density
       area = SafeDiv( H%aIce, H%CldFr, 0.0_dp )

       ! Skip if no area
       IF ( area > 0.0_dp ) THEN

          ! In-cloud loss frequency [1/s]
          ktmp = Ars_L1K( area, H%rIce, gamIce, srMw )
          kI   = kI + ktmp

          ! In-continue loud loss frequency for ice rxn branch [1/s]
          kIb  = kIb + ( ktmp * brIce )
       ENDIF
    ENDIF

    !------------------------------------------------------------------
    ! Mean branch ratio for reaction of interest in cloud
    ! (averaged over ice and liquid)
    !
    ! If the division can't be done, set kHet = 0 and return
    !------------------------------------------------------------------
    branch = SafeDiv( kiB, kI, 0.0_dp )
    IF ( .not. branch > 0.0_dp ) THEN
       kHet = 0.0_dp
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Grid-average loss frequency
    !
    ! EXACT expression for entrainment-limited uptake
    !------------------------------------------------------------------------

    ! Ratio (in cloud) of heterogeneous loss to detrainment, s/s
    kk = kI * tauc

    ! Ratio of volume inside to outside cloud
    ! ff has a range [0,+inf], so cap it at 1e30
    ff = SafeDiv( H%CldFr, H%ClearFr, 1.0e+30_dp )
    ff = MIN( ff, 1.0e+30_dp )

    ! Ratio of mass inside to outside cloud
    ! xx has range [0,+inf], but ff is capped at 1e30, so shouldn't overflow.
    xx =     ( ff        - kk        - 1.0_dp       ) / 2.0_dp +             &
         SQRT( 1.0_dp    + ff*ff     + kk*kk                   +             &
               2.0_dp*ff + 2.0_dp*kk - 2.0_dp*ff*kk ) / 2.0_dp

    ! Overall heterogeneous loss rate, grid average, 1/s
    ! kHet = kI * xx / ( 1d0 + xx )
    !  Since the expression ( xx / (1+xx) ) may behave badly when xx>>1,
    !  use the equivalent 1 / (1 + 1/x) with an upper bound on 1/x
    kHet = kI / ( 1.0_dp + SafeDiv( 1.0_dp, xx, 1.0e+30_dp ) )

    ! Overall loss rate in a particular reaction branch, 1/s
    kHet = kHet * branch
  END FUNCTION CloudHet

  !#########################################################################
  !#####         COMMON FUNCTIONS FOR COMPUTING UPTAKE RATES           #####
  !#########################################################################

  FUNCTION coth( x ) RESULT( f_x )
    !
    ! Hyperbolic cotangent = [1 + exp(-2x)] / [1 - exp(-2x)]
    !
    REAL(dp), INTENT(IN) :: x
    REAL(dp)             :: y, f_x
    !
    y   = EXP( -2.0_dp * x )
    f_x = ( 1.0_dp + y ) / ( 1.0_dp - y )
  END FUNCTION coth

  FUNCTION ReactoDiff_Corr( radius, l ) RESULT( corr )
    !
    ! For x = radius / l, correction =  COTH( x ) - ( 1/x )
    ! Correction approaches 1 as x becomes large, corr(x>1000)~1
    ! Correction approaches x/3 as x goes towards 0
    !
    REAL(dp), INTENT(IN)  :: l, radius           ! [cm] and [cm]
    REAL(dp)              :: x, corr
    !
    x = radius / l
    IF ( x > 1000.0_dp ) THEN
       corr = 1.0_dp
       RETURN
    ENDIF
    IF ( x < 0.1_dp ) THEN
       corr = x / 3.0_dp;
       RETURN
    ENDIF
    corr = coth(x) - ( 1.0_dp / x )
  END FUNCTION ReactoDiff_Corr

  !#########################################################################
  !#####   COMMON FUNCTIONS FOR ENFORCING SAFE NUMERICAL OPERATIONS    #####
  !#########################################################################

  FUNCTION SafeDiv( num, denom, alt ) RESULT( quot )
    !
    ! Performs "safe division", that is to prevent overflow, underlow,
    ! NaN, or infinity errors.  An alternate value is returned if the
    ! division cannot be performed.
    REAL(dp), INTENT(IN) :: num, denom, alt
    REAL(dp)             :: ediff, quot
    !
    ! Exponent difference (base 2)
    ! For REAL*8, max exponent = 1024 and min = -1021
    ediff = EXPONENT( num ) - EXPONENT( denom )
    !
    IF ( ediff > 1023 .OR. denom == 0.0_dp ) THEN
       quot = alt
    ELSE IF ( ediff < -1020 ) THEN
       quot = 0.0_dp
    ELSE
       quot = num / denom
    ENDIF
  END FUNCTION SafeDiv

  FUNCTION Is_SafeDiv( num, denom ) RESULT( safe )
    !
    ! Returns TRUE if a division can be performed safely.
    REAL(dp), INTENT(IN) :: num, denom
    LOGICAL              :: safe
    REAL(dp)             :: ediff
    !
    ! Exponent difference (base 2)
    ! For REAL*8, max exponent = 1024 and min = -1021
    safe  = .TRUE.
    ediff = EXPONENT( num ) - EXPONENT( denom )
    !
    IF ( ediff < -1020 .or. ediff > 1023 .or. denom == 0.0_dp ) THEN
       safe = .FALSE.
    ENDIF
  END FUNCTION Is_SafeDiv

  FUNCTION IsSafeExp( x ) RESULT( safe )
    !
    ! Returns TRUE if an exponential can be performed safely
    !
    REAL(dp), INTENT(IN) :: x
    LOGICAL              :: safe
    !
    ! Note EXP( 708 ) = 8.2e+307 and EXP( -708 ) = 3.3e-308, which are
    ! very close to the maximum representable values at double precision.
    safe = ( ABS( x ) < 709.0_dp )
  END FUNCTION IsSafeExp

  FUNCTION SafeExp( x, alt ) RESULT( y )
    !
    ! Performs a "safe exponential", that is to prevent overflow, underflow,
    ! underlow, NaN, or infinity errors when taking the value EXP( x ).  An
    ! alternate value is returned if the exponential cannot be performed.
    !
    REAL(dp), INTENT(IN) :: x, alt
    REAL(dp)             :: y
    !
    y = alt
    IF ( ABS( X ) < 709.0_dp ) y = EXP( x )
  END FUNCTION SafeExp

END MODULE rateLawUtilFuncs
!EOC
