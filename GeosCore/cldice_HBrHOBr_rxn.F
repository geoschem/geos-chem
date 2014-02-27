!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cldice_hbrhobr_rxn
!
! !DESCRIPTION: Subroutine CLDICE\_HBrHOBr\_RXN calculates the rate constants
! for HBr and HOBr pseudo-reactions with ice.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLDICE_HBrHOBr_RXN( I,      J,      L,    DENAIR, 
     &                               QI,     hbr,    hobr,
     &                               k_hbr,  k_hobr, AREA, State_Met )
! 
! !USES:
!
      USE ERROR_MOD,          ONLY : IS_SAFE_DIV, IT_IS_NAN
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE GIGC_State_Met_Mod, ONLY : MetState

      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      INTEGER,        INTENT(IN)  :: I         ! Longitude index
      INTEGER,        INTENT(IN)  :: J         ! Latitude  index
      INTEGER,        INTENT(IN)  :: L         ! Altitude  index
      REAL*8,         INTENT(IN)  :: DENAIR    ! Density of air         [#/cm3]
      REAL*8,         INTENT(IN)  :: QI        ! Cloud ice mixing ratio [kg/kg]
      REAL*8,         INTENT(IN)  :: hbr       ! Concentration of HBr   [#/cm3]
      REAL*8,         INTENT(IN)  :: hobr      ! Concentration of HOBr  [#/cm3]
      TYPE(MetState), INTENT(IN)  :: State_Met ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      REAL*8,         INTENT(OUT) :: k_hbr     ! Rate constant for HBr  + ice
                                               !   pseudo-rxn [cm3/s]
      REAL*8,         INTENT(OUT) :: k_hobr    ! Rate constant for HOBr + ice
                                               !   pseudo-rxn [cm3/s]
      REAL*8,         INTENT(OUT) :: AREA      ! Surface area [cm2/cm3]
!
! !REMARKS:
!  The rate constant is calculated assuming:
!                                                                             .
!    1. A sticking coefficient of 0.1 [JPL 2006], Abbatt [1994],
!       Chai et al. [2000]
!    2. An effective radius is assumed as a function of (i) temperature and
!       ice water content (IWC). This relationship is taken from Wyser [1998].
!                                                                             .
!  ** Calculations of a 1st order rate constent are borrowed from the
!     subroutine arsl1k.F. Below are comments from that code:
!                                                                             .
!       The 1st-order loss rate on wet aerosol (Dentener's Thesis, p. 14)
!       is computed as:
!                                                                             .
!         ARSL1K [1/s] = area / [ radius/dfkg + 4./(stkcf * nu) ]        
!                                                                             .
!       where nu   = Mean molecular speed [cm/s] = sqrt(8R*TK/pi/M) for Maxwell
!             DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
!
! !REVISION HISTORY:
!  16 Jun 2011 - J. Parrella - Initial version
!  22 May 2012 - M. Payer    - Added ProTeX headers
!  26 Sep 2012 - R. Yantosca - For now, comment out debug print statements
!  23 Oct 2012 - R. Yantosca - Add better error checks to prevent div-by-zero
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  05 Sep 2013 - R. Yantosca - Now exit if IWC <= 0, this will cause the
!                              logarithm to choke
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      REAL*8, PARAMETER :: R        = 8.314472                  ! [J/mol/K]
      REAL*8, PARAMETER :: pi       = 3.14159265358979323846d0  ! [unitless]
      REAL*8, PARAMETER :: dens_ice = 0.9167d-3                 ! [kg/cm3]
      REAL*8, PARAMETER :: mw_hbr   = 0.081                     ! [kg/mol]
      REAL*8, PARAMETER :: mw_hobr  = 0.097                     ! [kg/mol]
!
! !LOCAL VARIABLES:
!
      REAL*8            :: MW_spec  ! Molecular weight of species
      REAL*8            :: RADIUS   ! Radius of ice particles        [cm]
      REAL*8            :: STK      ! Square root of the temperature [K]
      REAL*8            :: DFKG     ! Gas diffusion coefficient      [cm2/s]
      REAL*8            :: Vc       ! Volume of ice in cloud         [cm3]
      REAL*8            :: SQM_hbr  ! Square root of molec. weight   [g/mol]
      REAL*8            :: SQM_hobr ! Square root of molec. weight   [g/mol]
      REAL*8            :: b_param
      REAL*8            :: XAIRCM3  ! Volume of air                  [cm3]
      REAL*8            :: iwc
      REAL*8            :: gamma
      REAL*8            :: gamma_hobr, gamma_hbr
      REAL*8            :: hbr_rtemp,  hobr_rtemp
      REAL*8            :: cld1k_hbr,  cld1k_hobr
      LOGICAL           :: yn_div_safe
      LOGICAL           :: yn_nan
      LOGICAL           :: yn_stop

      ! Pointers
      REAL*8, POINTER   :: AD(:,:,:)
      REAL*8, POINTER   :: AIRVOL(:,:,:)
      REAL*8, POINTER   :: CLDF(:,:,:)
      REAL*8, POINTER   :: T(:,:,:)

      !=================================================================
      ! CLDICE_HBrHOBr_RXN begins here!
      !=================================================================

      ! Initialize pointers
      AD     => State_Met%AD
      AIRVOL => State_Met%AIRVOL
      CLDF   => State_Met%CLDF
      T      => State_Met%T

      ! ----------------------------------------------
      ! 1.
      !   Calculate the ice water content (IWC) for
      !   the box. Will be used to calculate the
      !   ice effective radius by parameterization.
      !
      !   For parameterization, we want IWC in [g/m3].
      !   Using QI which is in [kg_ice/kg_air]
      ! ----------------------------------------------

      ! IWC in [g/cm3]
#if defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
      iwc = CLDF(L,I,J) * QI * AD(I,J,L) / AIRVOL(I,J,L)
      iwc = iwc / 1.d3
#else
      ! Convert cm3/cm3 to g/cm3
      ! (cm3 ice * 0.92 g ice/cm3 H2O ) / cm3 air
      iwc = CLDF(L,I,J) * ( QI * 1d3 * dens_ice )
#endif

      !%%% ERROR CHECK!  Do not let IWC<=0 because it will cause the
      !%%% LOG10 statement in the upcoming section (bmy, 10/23/12)
!-----------------------------------------------------------------------------
! Prior to 9/5/13:
! Amend error trap, IWC < 0 will also choke because it's an argument in
! a logarithm.  Now exit unless IWC is a positive number. (bmy, 9/5/13)
!      IF ( ( .not. ( ABS( IWC ) > 0d0 ) ) ) THEN
!-----------------------------------------------------------------------------
      IF ( .not. ( IWC > 0d0 ) ) THEN
         gamma  = 0.d0
         k_hbr  = 0.d0
         k_hobr = 0.d0
         RETURN
      ENDIF

      ! ----------------------------------------------
      ! calculate the temperature dependent reactive
      ! uptake coefficient for HBr + HOBr + ice.
      !  - uses Chaix et al. [2000] for 180K < T < 205 K
      !                                 0.44 > g > 0.15
      !
      !  - uses Abbatt [1994] for T = 228 K (pure ice)
      !                           g = 0.12
      !
      !  ** using g = 0.1 for cold and mixed clouds
      ! ----------------------------------------------
      IF ( ( T(I,J,L) >= 180.d0) .and. ( T(I,J,L) <= 268.d0) ) THEN
         gamma = 0.1d0
      ELSE
         ! if temperature moves above 228K, then
         ! turn off the reaction and do nothing
         gamma  = 0.d0
         k_hbr  = 0.d0
         k_hobr = 0.d0
         RETURN
      ENDIF

      ! set the sticking coefficients for HBr and
      ! HOBr independently, following JPL [2011]
      ! ** this is only used for when we treat the
      !    ice uptake of HOBr and HBr as a loss
      !    instead of recycling.
      gamma_hobr = 0.003d0
      gamma_hbr  = 0.03d0
      ! ----------------------------------------------
      ! 2.
      !   calculate the surface area of cloud droplets
      !   in the given grid box, assuming 1 of 2
      !   conditions:
      !     a. marine warm cloud
      !       or
      !     b. continental warm cloud
      !
      !
      !   * Calculation for area is derived follows,
      !     assuming that RADIUS is constant:
      !
      !                         4/3 (pi) (RADIUS)**3
      !  1) FC = Vc / Vb = N  -------------------------
      !                                  Vb
      !
      !
      !       where N      = number of ice particles in cloud
      !             RADIUS = radius of ice particles
      !             Vc     = volumn of ice in cloud
      !             Vb     = volumn of the box = AIRVOL (in GEOS-Chem)
      !
      !
      !                     Vb
      !  2) N = FC --------------------
      !            4/3 (pi) (RADIUS)**3
      !
      !
      !  So the surface area [m2] is calculated as
      !
      !  3) total surface A = N * 4 * (pi) * (RADIUS)**2
      !
      !                  3*Vb
      !          = FC ----------
      !                 RADIUS
      !
      !  4) for this routine though we want
      !     AREA in [cm2/cm3], surface area to volume air:
      !
      !                   3
      !     AREA = FC ---------
      !                RADIUS (in cm)
      !
      !
      !    or    
      !                   3 x Vc
      !     AREA =  -----------------
      !              AIRVOL x RADIUS      (in cm)
      ! ----------------------------------------------

      ! calculate the effective radius of ice in
      ! this cloud given (a) the temperature, and
      ! (b) the ice water content of the box.
      ! This parameterization is taken from
      ! Wyser [1998]. (jpp, 6/15/2011)

      ! a. B = equation 14 in Wyser [1998]
      b_param = -2.d0 + 1.d-3 * (273.d0 - T(I,J,L))**(1.5d0) *
     &     log10(iwc/50.d0)

      ! radius is parameterized in [um]
      RADIUS  = 377.4d0 + 203.3d0 * b_param +
     &     37.91d0  * b_param**2.0 +
     &     2.3696d0 * b_param**3.0

      ! now convert radius to [cm]
      RADIUS = RADIUS * 1.d-4

!------------------------------------------------------------------------------
!      ! make sure there's enough ice in the box we're
!      ! looking at to do chemistry:
!      IF ( (iwc == 0.0) .or. (QI == 0.0) .or.
!     &     (it_is_nan( log10(iwc/50.0))) .or.
!     &     (CLDF == 0.0) .or. (RADIUS <= 0.0) ) THEN
!         gamma  = 0.d0
!         k_hbr  = 0.d0
!         k_hobr = 0.d0
!         RETURN
!      ENDIF
!------------------------------------------------------------------------------
      !%%% ERROR CHECK!  Do not let RADIUS=0 because it will cause the
      !%%% LOG10 statement in the upcoming section (bmy, 10/23/12)
      IF ( (              RADIUS   < 0     )        .or.
     &     ( .not. ( ABS( RADIUS ) > 0d0 ) )        .or.
     &     ( .not. ( ABS( CLDF(L,I,J) ) > 0d0 ) ) ) THEN
         gamma  = 0.d0
         k_hbr  = 0.d0
         k_hobr = 0.d0
         RETURN
      ENDIF

      ! store the volume of air [cm3]
      XAIRCM3 = AIRVOL(I,J,L) * (100.d0)**3

#if defined( GEOS_5 ) || defined( MERRA ) || defined( GEOS_FP )
      ! get the volume of cloud ice [cm3]
      Vc = CLDF(L,I,J) * QI * AD(I,J,L) / dens_ice
#else
      ! Doesn't matter, Vc is not used.
#endif

      ! now calculate the cloud ice surface area [cm2_ice/ cm3_air]
!jpt      AREA    = 3.d0 * (Vc/XAIRCM3) / (RADIUS) ! keep Radius in [cm]

      ! Convert cross-sectional area to surface area, making assumption
      ! given in Lawrence and Crutzen [1998]
      ! iwc must be in g/m3, not g/cm3 for Lawrence param.
      AREA    = 1.d-4 * (iwc * (100.d0)**3)**(0.9d0)
      AREA    = 2.d0 * AREA ! in [cm2/cm3]

      !### Debug
      !IF ( AREA > 1.0d-4 ) THEN
      !   print *, 'jpp: debugging cloud ice area'
      !   print *, 'AREA = ', AREA, ' cm2/cm3'
      !   print *, 'IWC  = ', iwc
      !ENDIF

      ! ----------------------------------------------------
      ! 3.
      !   Now finish calculating the 1st order rate
      !   constant, first for HBr, and then for HOBr
      !   to get rate constants for both pseudo-reactions.
      !
      !   (a) calculate the gas phase diffusion coefficient;
      !
      !   (b) calculate the hydrolysis rxn rate.
      ! ----------------------------------------------------
      SQM_hbr  = sqrt(mw_hbr  * 1.d3)      ! square root of molar mass [g/mole]
      SQM_hobr = sqrt(mw_hobr * 1.d3)      ! square root of molar mass [g/mole]
      STK      = sqrt(T(I,J,L)      )      ! square root of temperature [K]

      ! ----------------------
      ! i.) Deal with HBr
      ! ----------------------
      ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
      DFKG  = 9.45D17/DENAIR * STK * SQRT(3.472D-2 + 
     &     1.D0/(SQM_hbr*SQM_hbr))

      ! Compute ARSL1K according to the formula listed above
      cld1k_hbr = AREA / ( RADIUS/DFKG + 2.749064E-4 
     &     * SQM_hbr/(gamma*STK) )
!jp_loss     &     * SQM_hbr/(gamma_hbr*STK) )

      ! ----------------------
      ! ii.) Deal with HOBr
      ! ----------------------
      ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
      DFKG  = 9.45D17/DENAIR * STK * SQRT(3.472D-2 + 
     &     1.D0/(SQM_hobr*SQM_hobr))

      ! Compute ARSL1K according to the formula listed above
      cld1k_hobr = AREA / ( RADIUS/DFKG + 2.749064E-4 
     &     * SQM_hobr/(gamma*STK) )
!jp_loss     &     * SQM_hobr/(gamma_hobr*STK) )

      ! ----------------------------------------------------
      ! 4.
      !   Now test which loss rate (HOBr or HBr) is
      !   limiting.
      ! ----------------------------------------------------
      ! initial loss rates
      hbr_rtemp  = cld1k_hbr  * hbr
      hobr_rtemp = cld1k_hobr * hobr

      ! ---------------------------------------------
      ! kludging the rates to be equal to one another
      ! to avoid having to keep setting equality in
      ! SMVGEAR solver. (jpp, 5/10/2011)
      ! ---------------------------------------------
      IF ( hbr_rtemp > hobr_rtemp ) THEN

         ! 1. is it safe to divide?
         yn_div_safe = is_safe_div( 
     &        cld1k_hobr * hobr, hbr )
         IF (yn_div_safe) THEN
            ! 2. if it is safe, then go ahead
            cld1k_hbr = hobr_rtemp / hbr
         ELSE
            !    if not, then set rates really small...
            !    b/c the largest contributor is very small.
            cld1k_hobr = TINY(1.d0)
            cld1k_hbr  = TINY(1.d0)
         ENDIF

      ELSE ! if HOBr rate is larger than HBr rate
           ! 1. is it safe to divide?
         yn_div_safe = is_safe_div( 
     &        cld1k_hbr * hbr, hobr )
         IF (yn_div_safe) THEN
            ! 2. if it is safe, then go ahead
            cld1k_hobr = hbr_rtemp / hobr
         ELSE
            !    if not, then set rates really small...
            !    b/c the largest contributor is very small.
            cld1k_hobr = TINY(1.d0)
            cld1k_hbr  = TINY(1.d0)
         ENDIF

      ENDIF

      ! store the rate constants
      k_hbr  = cld1k_hbr
      k_hobr = cld1k_hobr

      ! test if they're are NaN's calculated
      yn_stop = .FALSE.
      yn_nan = it_is_nan( k_hbr  )
      IF ( yn_nan ) yn_stop = .TRUE.

      yn_nan = it_is_nan( k_hobr )
      IF ( yn_nan ) yn_stop = .TRUE.

      IF (yn_stop) THEN
         print*, 'stopping inside of cldice_hbrhobr_rxn().'
         print*, 'Calculated NaN rate constants for ice chem.'
         print*, 'debugging values:'
         print*, 'k_hobr =', k_hobr
         print*, 'k_hbr =', k_hbr
         print*, 'hobr =', hobr
         print*, 'hbr  =', hbr
         print*, 'radius =', radius
         print*, 'area   =', area
         print*, 'b_param =', b_param
         print*, 'iwc =', iwc
!         read(*,*) yn_stop
         CALL GEOS_CHEM_STOP
      ENDIF

      ! Free pointers
      NULLIFY( AD     )
      NULLIFY( AIRVOL )
      NULLIFY( CLDF   )
      NULLIFY( T      )

      RETURN
!EOC
      END SUBROUTINE CLDICE_HBrHOBr_RXN
