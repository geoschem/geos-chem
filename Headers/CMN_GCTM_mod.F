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
      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS: 
!
      ! AIRMW : Molecular weight of air [28.97 g/mole]
      REAL*8, PARAMETER :: AIRMW        = 28.97d0

      ! Avogadro's # [#/mol]
      REAL*8, PARAMETER :: AVO          = 6.022d23

      ! g0    : Gravity at Surface of Earth [9.8 m/s^2]
      REAL*8, PARAMETER :: g0           =   9.8d0                 

      ! g0_100 = 100.0 / g0
      REAL*8, PARAMETER :: g0_100       = 100d0 / g0

      ! PI    : Double-Precision value of PI          
      REAL*8, PARAMETER :: PI           = 3.14159265358979323d0 

      ! PI_180 = PI    / 180.0
      REAL*8, PARAMETER :: PI_180       = PI / 180d0

      ! Re    : Radius of Earth [m] 
      REAL*8, PARAMETER :: Re           = 6.375d6               

      ! Rd    : Gas Constant (R) in Dry Air [287 J/K/kg] 
      REAL*8, PARAMETER :: Rd           = 287.0d0                 

      ! Rdg0   = Rd    / g0
      REAL*8, PARAMETER :: Rdg0         = Rd / g0

      ! Scale height of atmosphere (7.6 km = 7600m)
      REAL*8, PARAMETER :: SCALE_HEIGHT = 7600d0

      ! Von Karman's constant
      REAL*8, PARAMETER :: VON_KARMAN   = 0.4d0
!
! !REVISION HISTORY:
!  25 Jun 2002 - R. Yantosca - Initial version
!  23 Aug 2011 - M. Long     - Converted to Module from Header file
!  23 Jul 2014 - R. Yantosca - Add Von Karman's constant here
!  23 Jul 2014 - R. Yantosca - Add Avogadro's number here
!  23 Jul 2014 - R. Yantosca - List constants in alphabetical order
!EOP
!------------------------------------------------------------------------------
!BOC
      END MODULE CMN_GCTM_MOD
!EOC
