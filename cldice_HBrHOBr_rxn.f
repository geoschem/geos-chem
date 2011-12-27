      SUBROUTINE cldice_HBrHOBr_rxn( DENAIR, airvol, !I
     &     temp, QI, CLDF, AD, hbr, hobr,            !I
     &     k_hbr, k_hobr, AREA)                      !O
! -------------------------------------------------------------------------
! This subroutine calculates the rate constant
! for heterogeneous cycling of BrNO3 off of
! cloud particles, assuming:
!
!  1. A sticking coefficient of 0.1 [JPL 2006], Abbatt [1994], Chai et al. [2000]
!  2. An effective radius is assumed as a function of
!     (i) temperature and (ii) ice water content (IWC). This
!     relationship is taken from Wyser [1998].
!
!
! ** Calculation of a 1st order rate constent barrowed from 
!    the subroutine arsl1k.f. Below are comments from that
!    code.
!
! !REMARKS:
!  The 1st-order loss rate on wet aerosol (Dentener's Thesis, p. 14)
!  is computed as:
!                                                                             .
!      ARSL1K [1/s] = area / [ radius/dfkg + 4./(stkcf * nu) ]        
!                                                                             .
!  where nu   = Mean molecular speed [cm/s] = sqrt(8R*TK/pi/M) for Maxwell 
!        DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
!
! Input Variables:
! `````````````````
! ( 1) I, J, K  :: spatial indices
! ( 2) DENAIR   :: air density [#/cm3]
! ( 3) hbr      :: concentration of HBr  in [#/cm3]
! ( 4) hobr     :: concentration of HOBr in [#/cm3]
! ( 5) AIRVOL   :: volume of air in a box [m3]
! ( 6) temp     :: temperature in a given box [K]
! ( 7) QI       :: cloud ice mixing ratio [kg/kg]
! ( 8) AD       :: dry air mass [kg]
! ( 9) CLDF     :: 3D cloud fraction
!
! Output Variables:
! `````````````````
! ( 1) k_hbr  :: rate constant for HBr  + ice pseudo-rxn [cm3/s]
! ( 2) k_hobr :: rate constant for HOBr + ice pseudo-rxn [cm3/s]
!
!
! Variables:
! ```````````
! ( 1) SQM      ::  square root of the molecular weight [g/mol]
! ( 2) STK      ::  square root of the temperature [K]
! ( 3) DFKG     ::  gas diffusion coefficient [cm2/s]
! ( 4) 
!
! -------------------------------------------------------------------------
!     Justin Parrella, 6/16/2011
!     parrella@fas.harvard.edu
!
! -------------------------------------------------------------------------

      USE ERROR_MOD,       ONLY : is_safe_div, it_is_nan
      USE ERROR_MOD,       ONLY : geos_chem_stop

      IMPLICIT NONE

      ! --------------------
      ! input the variables
      ! --------------------
      REAL*8,  INTENT(IN)  :: hbr, hobr
      REAL*8,  INTENT(IN)  :: DENAIR, AIRVOL, temp, QI
      REAL*8,  INTENT(IN)  :: CLDF, AD

   
      ! --------------------
      ! the output variables
      ! --------------------
      REAL*8,  INTENT(OUT)   :: k_hbr, k_hobr, AREA
!      REAL*8,  INTENT(INOUT) :: AREA
   
      ! --------------------
      ! local variables
      ! --------------------
      real*8  :: MW_spec ! molecular weight of species
      real*8  :: RADIUS, STK, DFKG, Vc
      real*8  :: SQM_hbr, SQM_hobr, b_param
      real*8  :: XAIRCM3
      real*8  :: iwc, gamma, gamma_hobr, gamma_hbr
      real*8  :: hbr_rtemp, hobr_rtemp
      real*8  :: cld1k_hbr, cld1k_hobr
      logical :: yn_div_safe
      logical :: yn_nan, yn_stop
   
      ! --------------------
      ! parameters
      ! --------------------
      REAL*8, PARAMETER :: R = 8.314472 ! J /mol /K
      real*8, parameter :: pi = 3.14159265358979323846d0
      real*8, parameter :: dens_ice = 0.9167d-3 ! kg/cm3
      real*8, parameter :: mw_hbr  = 0.081 ! kg/mol
      real*8, parameter :: mw_hobr = 0.097 ! kg/mol


      ! ----------------------------------------------
      ! 1.
      !   Calculate the ice water content (IWC) for
      !   the box. Will be used to calculate the
      !   ice effective radius by parameterization.
      !
      !   For parameterization, we want IWC in [g/m3].
      !   Using QI which is in [kg_ice/kg_air]
      ! ----------------------------------------------

      ! IWC in [kg_ice/m3]; note, QI is [kg/kg] m.r.
      iwc = CLDF * QI * AD / AIRVOL   ! AIRVOL is in [m3], AD in [kg]
      iwc = iwc / 1.d3  ! IWC in [g_ice/cm3]

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
      if ( (temp >= 180.d0) .and. (temp <= 268.d0) ) then
         gamma = 0.1d0
      else ! if temperature moves above 228K, then
           ! turn off the reaction and do nothing
         gamma  = 0.d0
         k_hbr  = 0.d0
         k_hobr = 0.d0
         RETURN
      endif

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
      b_param = -2.d0 + 1.d-3 * (273.d0 - temp)**(1.5d0) *
     &     log10(iwc/50.d0)

      ! radius is parameterized in [um]
      RADIUS  = 377.4d0 + 203.3d0 * b_param +
     &     37.91d0  * b_param**2.0 +
     &     2.3696d0 * b_param**3.0

      ! now convert radius to [cm]
      RADIUS = RADIUS * 1.d-4

      ! make sure there's enough ice in the box we're
      ! looking at to do chemistry:
      if ( (iwc == 0.0) .or. (QI == 0.0) .or.
     &     (it_is_nan( log10(iwc/50.0))) .or.
     &     (CLDF == 0.0) .or. (RADIUS <= 0.0) ) then
         gamma  = 0.d0
         k_hbr  = 0.d0
         k_hobr = 0.d0
         RETURN
      endif

      ! store the volume of air [cm3]
      XAIRCM3 = AIRVOL * (100.d0)**3

      ! get the volume of cloud ice [cm3]
      Vc = CLDF * QI * AD / dens_ice

      ! now calculate the cloud ice surface area [cm2_ice/ cm3_air]
!jpt      AREA    = 3.d0 * (Vc/XAIRCM3) / (RADIUS) ! keep Radius in [cm]

      ! Convert cross-sectional area to surface area, making assumption
      ! given in Lawrence and Crutzen [1998]
      AREA    = 1.d-4 * (iwc * (100.d0)**3)**(0.9d0) ! iwc must be in g/m3, not g/cm3 for Lawrence param.
      AREA    = 2.d0 * AREA ! in [cm2/cm3]

      if ( AREA > 1.0d-4 ) then
         print *, 'jpp: debugging cloud ice area'
         print *, 'AREA = ', AREA, ' cm2/cm3'
         print *, 'IWC  = ', iwc
      endif

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
      SQM_hbr  = sqrt(mw_hbr  * 1.d3) ! square root of molar mass [g/mole]
      SQM_hobr = sqrt(mw_hobr * 1.d3) ! square root of molar mass [g/mole]
      STK = sqrt(temp)                ! square root of temperature [K]

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
      if ( hbr_rtemp > hobr_rtemp ) then

         ! 1. is it safe to divide?
         yn_div_safe = is_safe_div( 
     &        cld1k_hobr * hobr, hbr )
         if (yn_div_safe) then
            ! 2. if it is safe, then go ahead
            cld1k_hbr = hobr_rtemp / hbr
         else
            !    if not, then set rates really small...
            !    b/c the largest contributor is very small.
            cld1k_hobr = TINY(1.d0)
            cld1k_hbr  = TINY(1.d0)
         endif

      else ! if HOBr rate is larger than HBr rate
           ! 1. is it safe to divide?
         yn_div_safe = is_safe_div( 
     &        cld1k_hbr * hbr, hobr )
         if (yn_div_safe) then
            ! 2. if it is safe, then go ahead
            cld1k_hobr = hbr_rtemp / hobr
         else
            !    if not, then set rates really small...
            !    b/c the largest contributor is very small.
            cld1k_hobr = TINY(1.d0)
            cld1k_hbr  = TINY(1.d0)
         endif

      endif

      ! store the rate constants
      k_hbr  = cld1k_hbr
      k_hobr = cld1k_hobr

      ! test if they're are NaN's calculated
      yn_stop = .false.
      yn_nan = it_is_nan( k_hbr  )
      if ( yn_nan ) yn_stop = .true.

      yn_nan = it_is_nan( k_hobr )
      if ( yn_nan ) yn_stop = .true.

      if (yn_stop) then
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
         call geos_chem_stop
      endif



      RETURN

      end subroutine cldice_HBrHOBr_rxn
