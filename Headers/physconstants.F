!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: physconstants.F
!
! !DESCRIPTION: PhysConstants contains GEOS-Chem specific PHYSICAL CONSTANTS
!  and DERIVED QUANTITIES.
!\\
!\\
! !INTERFACE: 
!
      MODULE PHYSCONSTANTS
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS: 
!
      ! AIRMW : Average molecular weight of dry air [g/mol]
      REAL(fp), PARAMETER :: AIRMW    = 28.97e+0_fp

      ! H2OMW : Molecular weight of water [g/mol]
      REAL(fp), PARAMETER :: H2OMW   = 18.016e+0_fp

      ! AVO   : Avogadro's number [particles/mol]
      ! Now use more precise value 6.022140857e+23 instead of 6.022e+23?
      ! Source: NIST, 2014 (ewl, 1/7/16) (NEED TO CHANGE HEMCO)
      REAL(fp), PARAMETER :: AVO      = 6.022140857e+23_fp

      ! g0    : Acceleration due to gravity at earth's surface [m/s^2]
      ! Now use more precise value of 9.80665 instead of 9.8 (ewl, 1/7/16)
      ! Source: NIST, 2014 (NEED TO CHANGE HEMCO)
      REAL(fp), PARAMETER :: g0       =   9.80665e+0_fp                 

      ! g0_100 : 100 / g0
      REAL(fp), PARAMETER :: g0_100   = 100.e+0_fp / g0

      ! PI    : Double-Precision value of PI          
      REAL(fp), PARAMETER :: PI       = 3.14159265358979323e+0_fp 

      ! PI_180 : Number of radians per degree
      REAL(fp), PARAMETER :: PI_180   = PI / 180e+0_fp

      ! Re    : Radius of Earth [m] 
      REAL(fp), PARAMETER :: Re       = 6.375e+6_fp               

      ! Rd    : Gas Constant in Dry Air [J/K/kg] 
      REAL(fp), PARAMETER :: Rd       = 287.0e+0_fp 

      ! Rv    : Gas Constant for water vapor [J/K/kg] 
      REAL(fp), PARAMETER :: Rv       = 461.0e+0_fp   

      ! Rdg0   = Rd    / g0
      REAL(fp), PARAMETER :: Rdg0     = Rd / g0

      ! SCALE_HEIGHT : Scale height of atmosphere [m]
      REAL(fp), PARAMETER :: SCALE_HEIGHT = 7600.e+0_fp

      ! VON_KARMAN : Von Karman's constant [.]
      REAL(fp), PARAMETER :: VON_KARMAN   = 0.4e+0_fp

      ! RSTARG : Molar gas constant [J/K/mol]
      ! Now use more precise value 8.3144598 instead of 8.31450? (ewl, 1/7/16)
      ! Source: NIST, 2014 (NEED TO CHANGE HEMCO)
      REAL(fp), PARAMETER :: RSTARG   = 8.3144598e+0_fp
      
      ! XNUMOLAIR : Molecules dry air per kg dry air
      REAL(fp), PARAMETER :: XNUMOLAIR = AVO / ( AIRMW * 1.e-3_fp)

      ! BOLTZ : Boltzmann's constant [J/K]  (Source: NIST, 2014)
      REAL(fp), PARAMETER :: BOLTZ   = 1.38064852e-23_fp

      ! ATM : Standard atmosphere [Pa]  (Source: NIST, 2014)
      REAL(fp), PARAMETER :: ATM     = 1.01325e+5_fp

      ! Condensation vapor pressure
      ! ** NEED SOURCE **
      !  We think 6.1078 hPa is the saturation vapor pressure at 273.16 K, the
      !   triple point of water, but this needs to be confirmed (mps, 4/21/16)
      !  Use BOLTZ [J/K] rather than BOLTG [ergs/K] from comode_loop_mod
      !   (ewl, 1/4/16)
      REAL(fp), PARAMETER :: CONSVAP  = 6.1078e+03_fp /
     &                                  ( BOLTZ * 1e+7_fp )

      ! Gas constant in: [L.atm/K.mole]
      REAL(fp), PARAMETER :: RGASLATM = 8.2057e-2_fp

      ! Molecular weight of carbon (kg/mol)
      REAL(fp), PARAMETER :: MWCARB    = 12.01e-3_fp

!
! !REFERENCES:
! (1) NIST, 2014. Website: http://physics.nist.gov/cuu/Constants/index.html
!
! !REVISION HISTORY:
!  25 Jun 2002 - R. Yantosca - Initial version
!  23 Aug 2011 - M. Long     - Converted to Module from Header file
!  23 Jul 2014 - R. Yantosca - Add Von Karman's constant here
!  23 Jul 2014 - R. Yantosca - Add Avogadro's number here
!  23 Jul 2014 - R. Yantosca - List constants in alphabetical order
!  02 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!  06 Jan 2014 - E. Lundgren - Added constants and clarified comments
!  24 Mar 2015 - E. Lundgren - Add xnumolair, previously in tracer_mod
!  08 Jan 2016 - E. Lundgren - Update values to NIST 2014 and use AIRMW in
!                              definition of XNUMOLAIR. Add BOLTZ.
!  08 Jan 2016 - E. Lundgren - Rename this file to physconstants.F. Previously
!                              it was called CMN_GCTM_mod.F
!  21 Apr 2016 - M. Sulprizio- Move CONSVAP here from comode_loop_mod.F
!EOP
!------------------------------------------------------------------------------
!BOC
      END MODULE PHYSCONSTANTS
!EOC
