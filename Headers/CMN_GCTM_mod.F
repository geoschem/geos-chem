!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_GCTM_mod
!
! !DESCRIPTION: CMN\_GCTM\_mod contains GEOS-Chem specific PHYSICAL CONSTANTS
!  and DERIVED QUANTITIES.
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_GCTM_MOD
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS: 
!
      ! AIRMW : Molecular weight of air [28.97 g/mole]
      REAL(fp), PARAMETER :: AIRMW    = 28.97e+0_fp

      ! Avogadro's # [#/mol]
      REAL(fp), PARAMETER :: AVO      = 6.022e+23_fp

      ! g0    : Gravity at Surface of Earth [9.8 m/s^2]
      REAL(fp), PARAMETER :: g0       =   9.8e+0_fp                 

      ! g0_100 = 100.0 / g0
      REAL(fp), PARAMETER :: g0_100   = 100e+0_fp / g0

      ! PI    : Double-Precision value of PI          
      REAL(fp), PARAMETER :: PI       = 3.14159265358979323e+0_fp 

      ! PI_180 = PI    / 180.0
      REAL(fp), PARAMETER :: PI_180   = PI / 180e+0_fp

      ! Re    : Radius of Earth [m] 
      REAL(fp), PARAMETER :: Re       = 6.375e+6_fp               

      ! Rd    : Gas Constant (R) in Dry Air [287 J/K/kg] 
      REAL(fp), PARAMETER :: Rd       = 287.0e+0_fp                 

      ! Rdg0   = Rd    / g0
      REAL(fp), PARAMETER :: Rdg0     = Rd / g0

      ! Scale height of atmosphere (7.6 km = 7600m)
      REAL(fp), PARAMETER :: SCALE_HEIGHT = 7600e+0_fp

      ! Von Karman's constant
      REAL(fp), PARAMETER :: VON_KARMAN   = 0.4e+0_fp
!
! !REVISION HISTORY:
!  25 Jun 2002 - R. Yantosca - Initial version
!  23 Aug 2011 - M. Long     - Converted to Module from Header file
!  23 Jul 2014 - R. Yantosca - Add Von Karman's constant here
!  23 Jul 2014 - R. Yantosca - Add Avogadro's number here
!  23 Jul 2014 - R. Yantosca - List constants in alphabetical order
!  02 Dec 2014 - M. Yannetti - Added PRECISION_MOD
!EOP
!------------------------------------------------------------------------------
!BOC
      END MODULE CMN_GCTM_MOD
!EOC
