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
    ! Prevent div by zero.  Now use 1.0 as the error trap for concEduct.
    ! 100 and 1e-8 (the previous criteria) were too large and too small,
    ! respectively.  See https://github.com/geoschem/geos-chem/issues/1115.
    !  -- Seb Eastham, Bob Yantosca (09 Feb 2022)
    IF ( concEduct < 1.0_dp                              ) RETURN
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

    ! Do not let xx go negative, as this can cause numerical instability.
    ! See https://github.com/geoschem/geos-chem/issues/1205
    xx = MAX( xx, 0.0_dp )

    ! Overall heterogeneous loss rate, grid average, 1/s
    ! kHet = kI * xx / ( 1d0 + xx )
    !  Since the expression ( xx / (1+xx) ) may behave badly when xx>>1,
    !  use the equivalent 1 / (1 + 1/x) with an upper bound on 1/x
    kHet = kI / ( 1.0_dp + SafeDiv( 1.0_dp, xx, 1.0e+30_dp ) )

    ! Overall loss rate in a particular reaction branch, 1/s
    kHet = kHet * branch
  END FUNCTION CloudHet

  SUBROUTINE Cld_Params( AD, CLDF, FRLAND, FROCEAN, QI, QL, T, H )
    !
    ! Returns ice and liquid cloud parameters (based on State_Met)
    ! for cloud particles.
    !
    ! References:
    !  Heymsfield, A. J., Winker, D., Avery, M., et al. (2014). Relationships
    !   between ice water content and volume extinction coefficient from in
    !   situ observations for temperatures from 0° to –86°C: implications
    !   for spaceborne lidar retrievals. Journal of Applied Meteorology and
    !   Climatology, 53(2), 479–505. https://doi.org/10.1175/JAMC-D-13-087.1
    !
    !  Schmitt, C. G., & Heymsfield, A. J. (2005). Total Surface Area Estimates
    !   for Individual Ice Particles and Particle Populations. Journal of
    !   Applied Meteorology, 44(4), 467–474. https://doi.org/10.1175/JAM2209.1
    !
    REAL(dp),       INTENT(IN)    :: AD          ! Air mass [kg]
    REAL(dp),       INTENT(IN)    :: CLDF        ! Cloud fraction [1]
    REAL(dp),       INTENT(IN)    :: FRLAND      ! Land fraction [1]
    REAL(dp),       INTENT(IN)    :: FROCEAN     ! Ocean fraction [1]
    REAL(dp),       INTENT(IN)    :: QI          ! Ice mixing ratio [kg/kg]
    REAL(dp),       INTENT(IN)    :: QL          ! Liquid mixing ratio [kg/kg]
    REAL(dp),       INTENT(IN)    :: T           ! Temperature [K]
!
! !OUTPUT PARAMETERS:
!
    TYPE(HetState), INTENT(INOUT) :: H           ! Hetchem State object
!
! !REMARKS:
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
    IF ( ( QL + QI <= 0.0_dp ) .or. ( CLDF <= 0.0_dp ) ) THEN
       H%rLiq = CLDR_CONT
       H%rIce = CLDR_ICE
       H%ALiq = 0.0_dp
       H%VLiq = 0.0_dp
       H%AIce = 0.0_dp
       H%VIce = 0.0_dp
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
    IF ( FRLAND > FROCEAN ) THEN
       H%rLiq = CLDR_CONT      ! Continental cloud droplet radius [cm]
    ELSE
       H%rLiq = CLDR_MARI      ! Marine cloud droplet radius [cm]

    ENDIF

    ! get the volume of cloud condensate [cm3(condensate)/cm3(air)]
    ! QL is [g/g]
    H%VLiq = QL * AD / DENS_LIQ / H%vAir
    H%VIce = QI * AD / DENS_ICE / H%vAir
    H%ALiq = 3.0_dp * H%vLiq / H%rLiq

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
    IF ( T < 202.0_dp ) THEN          ! -71 C
       alpha = 83.3_dp
       beta  = 0.0184_dp
    ELSE IF ( T < 217.0_dp ) THEN     ! -56 C
       alpha = 9.1744e+4_dp
       beta  = 0.117_dp
    ELSE
       alpha = 308.4_dp
       beta  = 0.0152_dp
    ENDIF

    ! Effective radius, cm
    H%rIce = 0.5_dp * alpha * EXP( beta * ( T - 273.15_dp ) ) / 1e+4_dp

    ! Ice surface area density, cm2/cm3
    H%aIce = 3.0_dp * H%vIce / H%rIce * 2.25_dp

  END SUBROUTINE Cld_Params

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
!EOC
END MODULE rateLawUtilFuncs
