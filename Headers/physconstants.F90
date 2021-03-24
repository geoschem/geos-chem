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
  ! Average molecular weight of dry air [g/mol]
  REAL(fp), PARAMETER :: AIRMW = 28.9644_fp                 ! was 28.97

  ! Molecular weight of water [g/mol]
  REAL(fp), PARAMETER :: H2OMW = 18.016_fp

  ! Avogadro's number [particles/mol] (Source: NIST, 2014)
  REAL(fp), PARAMETER :: AVO = 6.022140857e+23_fp

  ! Acceleration due to gravity at earth's surface [m/s^2]
  ! (Source: NIST, 2014)
  REAL(fp), PARAMETER :: g0     = 9.80665e+0_fp
  REAL(fp), PARAMETER :: g0_100 = 100.0_fp / g0

  ! Double-Precision value of PI (and radians per degree)
  REAL(fp), PARAMETER :: PI     = 3.14159265358979323_fp
  REAL(fp), PARAMETER :: PI_180 = PI / 180.0_fp

  ! Radius of Earth [m]
  REAL(fp), PARAMETER :: Re = 6.3710072e+6_fp               ! was 6.375e+6_fp

  ! Gas Constant in Dry Air [J/K/kg] (and divided by g)
  REAL(fp), PARAMETER :: Rd   = 287.0_fp
  REAL(fp), PARAMETER :: Rdg0 = Rd / g0

  ! Gas Constant for water vapor [J/K/kg]
  REAL(fp), PARAMETER :: Rv = 461.00_fp

  ! Scale height of atmosphere [m]
  REAL(fp), PARAMETER :: SCALE_HEIGHT = 7600.0_fp

  ! Von Karman's constant [.]
  REAL(fp), PARAMETER :: VON_KARMAN = 0.4_fp

  ! Molar gas constant [J/K/mol] (Source: NIST, 2014)
  ! NOTE: Also be sure to update con_R in gckpp_Global if you update this!
  REAL(fp), PARAMETER :: RSTARG = 8.3144598_fp

  ! XNUMOLAIR : Molecules dry air per kg dry air
  REAL(fp), PARAMETER :: XNUMOLAIR = AVO / ( AIRMW * 1.e-3_fp )

  ! BOLTZ : Boltzmann's constant [J/K]  (Source: NIST, 2014)
  REAL(fp), PARAMETER :: BOLTZ = 1.38064852e-23_fp

  ! ATM : Standard atmosphere [Pa]  (Source: NIST, 2014)
  REAL(fp), PARAMETER :: ATM = 1.01325e+5_fp

  ! Condensation vapor pressure
  ! ** NEED SOURCE **
  !  We think 6.1078 hPa is the saturation vapor pressure at 273.16 K, the
  !   triple point of water, but this needs to be confirmed (mps, 4/21/16)
  !  Use BOLTZ [J/K] rather than BOLTG [ergs/K] from comode_loop_mod
  !   (ewl, 1/4/16)
  REAL(fp), PARAMETER :: CONSVAP = 6.1078e+03_fp / ( BOLTZ * 1e+7_fp )

  ! Gas constant in: [L.atm/K.mole]
  REAL(fp), PARAMETER :: RGASLATM = 8.2057e-2_fp

  ! Molecular weight of carbon (kg/mol)
  REAL(fp), PARAMETER :: MWCARB = 12.01e-3_fp
!
! !REFERENCES:
! (1) NIST, 2014. Website: http://physics.nist.gov/cuu/Constants/index.html
!
! !REVISION HISTORY:
!  25 Jun 2002 - R. Yantosca - Initial version
!  07 Jan 2016 - E. Lundgren - Updated to NIST 2014 values
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
END MODULE PHYSCONSTANTS
!EOC
