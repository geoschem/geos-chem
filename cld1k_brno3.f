      REAL*8 function CLD1K_BrNO3(I, J, L, DENAIR) RESULT(cld1k)
! -------------------------------------------------------------------------
! This subroutine calculates the rate constant
! for heterogeneous cycling of BrNO3 off of
! cloud particles, assuming:
!
!  1. A sticking coefficient of 0.3 [Yang et al. 2005]
!  2. uniform cloud droplet size for 2 types of clouds
!     - continental warm clouds: r =  6d-4 [cm]
!     - marine warm clouds:      r = 10d-4 [cm]
!     * no distributions are assumed
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
!
! Variables:
! ```````````
! ( 1) SQM      ::  square root of the molecular weight [g/mol]
! ( 2) STK      ::  square root of the temperature [K]
! ( 3) DFKG     ::  gas diffusion coefficient [cm2/s]
! ( 4) 
!
! -------------------------------------------------------------------------
!     implemented by Justin Parrella, 2/27/2011
!     parrella@fas.harvard.edu
!
! -------------------------------------------------------------------------


      ! use the following
      USE DAO_MOD,    ONLY : T, AIRVOL, FRLAND, FROCEAN
      USE DAO_MOD,    ONLY : QL   ! cloud water mixing ratio [kg/kg]
      USE DAO_MOD,    ONLY : CLDF ! 3D cloud fraction of the box
      USE DAO_MOD,    ONLY : AD   ! dry air mass [kg]

      IMPLICIT NONE

      ! --------------------
      ! input the variables
      ! --------------------
      INTEGER, intent(in) :: I, J, L
      ! Density of air [#/cm3]
      REAL*8,  INTENT(IN) :: DENAIR  
   
      ! --------------------
      ! the output variables
      ! --------------------
!      real*8  :: cld1k
   
      ! --------------------
      ! local variables
      ! --------------------
      real*8  :: nu ! mean molecular speed
      real*8  :: RADIUS, SQM, STK, AREA, DFKG, Vc
      real*8  :: XAIRM3, FC, temp ! , calc_vcldf
      logical :: yn_continue
   
      ! --------------------
      ! parameters
      ! --------------------
      REAL*8, PARAMETER :: XCLDR_CONT =  6.d-4 ! Cloud droplet radius in continental warm clouds [cm]
      REAL*8, PARAMETER :: XCLDR_MARI = 10.d-4 ! Cloud droplet radius in marine warm clouds [cm]
      REAL*8, PARAMETER :: R = 8.314472 ! J /mol /K
      real*8, parameter :: mw_brno3 = 0.142 ! kg/mol
      real*8, parameter :: pi = 3.14159265358979323846d0
      real*8, parameter :: alpha = 0.3 ! sticking coefficient
      real*8, parameter :: dens_h2o = 0.001d0 ! kg/cm3

      ! ----------------------------------------------
      ! 1.
      !   calculate the mean molecular speed of the
      !   molecules given the temperature.
      ! ----------------------------------------------
      temp = T( I, J, L ) ! [K]
      nu   = sqrt( 8.d0 * R * temp / (mw_brno3 * pi) )


      ! ----------------------------------------------
      ! Test conditions to see if we want to continue
      ! or set the cloud rate equal to zero.
      ! ----------------------------------------------

      ! Volume cloud fraction (Sundqvist et al 1989) [unitless]
!jpt      FC      = CALC_VCLDF(I,J,L)

      ! continental or marine clouds only...
      IF ( (FRLAND(I, J) > 0) .or. (FROCEAN(I,J) > 0) ) then
         ! do we have clouds? and do we have warm temperatures?
         IF ( (FC > 0) .and. (temp > 258.0) ) then
            yn_continue = .true.
         else
            yn_continue = .false.
         endif
      else
         yn_continue = .false.
      endif

      ! test
      if ( .not. yn_continue ) then
         ! nothing to calculate...
         cld1k = 0.d0
         return
      endif


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
      !       where N      = number of cloud droplets
      !             RADIUS = radius of cloud droplet
      !             Vc     = volumn of the cloud
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
      IF ( FRLAND(I, J) > FROCEAN(I, J) ) THEN
         ! Continental cloud droplet radius [cm]
         RADIUS = XCLDR_CONT
      ELSE
         ! Marine cloud droplet radius [cm]
         RADIUS = XCLDR_MARI
      ENDIF

      ! store the volume of air [m3]
      XAIRM3 = AIRVOL(I,J,L)
      ! convert to [cm3]
      XAIRM3 = XAIRM3 * (100.d0)**3

      ! get the volume of cloud [cm3]
!jpt      Vc = CLDF(L,I,J) * QL(I,J,L) * AD(I,J,L) / dens_h2o
      Vc = QL(I,J,L) * AD(I,J,L) / dens_h2o

      ! now calculate the cloud droplet surface area
!jpt      AREA    = 3.d0 * FC * XAIRM3 / (RADIUS * 1.d-2) ! convert radius to [m]
      AREA    = 3.d0 * (Vc/XAIRM3) / (RADIUS) ! keep Radius in [cm]

      ! ----------------------------------------------------
      ! 3.
      !   Now finish calculating the 1st order rate
      !   constant for BrNO3 hydrolysis.
      !
      !   (a) calculate the gas phase diffusion coefficient;
      !
      !   (b) calculate the hydrolysis rxn rate.
      ! ----------------------------------------------------
      SQM = sqrt(mw_brno3 * 1.d3) ! square root of molar mass [g/mole]
      STK = sqrt(temp)            ! square root of temperature [K]

      ! DFKG = Gas phase diffusion coeff [cm2/s] (order of 0.1)
      DFKG  = 9.45D17/DENAIR * STK * SQRT(3.472D-2 + 1.D0/(SQM*SQM))

      ! Compute ARSL1K according to the formula listed above
      cld1k = AREA / ( RADIUS/DFKG + 2.749064E-4 
     &     * SQM/(alpha*STK) )


      ! jpp, debugging
      print*, 'jpp inside of cld1k_brno3.f:'
      print*, 'cld1k =', cld1k

      end function CLD1K_BrNO3


!jpt!------------------------------------------------------------------------------
!jpt
!jpt      FUNCTION CALC_VCLDF(I, J, L) RESULT(VCLDF)
!jpt!
!jpt!******************************************************************************
!jpt!  Subroutine GET_VCLDF computes the volume cloud fraction for SO2 chemistry.
!jpt!  (rjp, bdf, bmy, 9/23/02)
!jpt!
!jpt!  References:
!jpt!  ============================================================================
!jpt!  (1) Sundqvist et al. [1989]
!jpt!
!jpt!  NOTES:
!jpt!  (1 ) Copied from 'sulfate_mod.f' for cloud uptake of GLYX and MGLY (tmf, 2/26/07)
!jpt!
!jpt!       ** jpp copied this from 'carbon_mod.f' (2/27/2011)
!jpt!******************************************************************************
!jpt!
!jpt      ! References to F90 modules 
!jpt      USE DAO_MOD,      ONLY : RH
!jpt      USE PRESSURE_MOD, ONLY : GET_PCENTER, GET_PEDGE
!jpt
!jpt      IMPLICIT NONE
!jpt
!jpt#     include "CMN_SIZE"   ! Size parameters
!jpt
!jpt      ! input variables
!jpt      INTEGER, intent(in)  :: I,    J,    L
!jpt
!jpt      ! output variables
!jpt      REAL*8               :: VCLDF
!jpt
!jpt      ! Local variables
!jpt      REAL*8               :: PRES, PSFC, RH2, R0, B0
!jpt
!jpt      ! Parameters
!jpt      REAL*8,  PARAMETER   :: ZRT = 0.60d0, ZRS = 0.99d0
!jpt		
!jpt      !=================================================================
!jpt      ! GET_VCLDF begins here!
!jpt      !=================================================================
!jpt	
!jpt      ! Surface pressure
!jpt      PSFC = GET_PEDGE(I,J,1)
!jpt
!jpt      ! Pressure at the center of the grid box
!jpt      PRES = GET_PCENTER(I,J,L)
!jpt
!jpt      ! RH (from "dao_mod.f") is relative humidity [%]
!jpt      ! Convert to fraction and store in RH2
!jpt      RH2  = RH(I,J,L) * 1.0d-2
!jpt
!jpt      ! Terms from Sundqvist ???
!jpt      R0   = ZRT + ( ZRS - ZRT ) * EXP( 1d0 - ( PSFC / PRES )**2.5 )
!jpt      B0   = ( RH2 - R0 ) / ( 1d0 - R0 )
!jpt        
!jpt      ! Force B0 into the range 0-1
!jpt      IF ( RH2 < R0  ) B0 = 0d0
!jpt      IF ( B0  > 1d0 ) B0 = 1d0
!jpt
!jpt      ! Volume cloud fraction
!jpt      VCLDF = 1d0 - SQRT( 1d0 - B0 )
!jpt
!jpt      ! Return to calling program
!jpt      END FUNCTION CALC_VCLDF
