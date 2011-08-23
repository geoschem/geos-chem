!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
! !MODULE: CMN_GCTM_MOD
!
! !DESCRIPTION: CMN_GCTM contains GEOS-CHEM specific PHYSICAL CONSTANTS
! and DERIVED QUANTITIES (bmy, 6/25/02, 6/23/03)
!
! REMARKS:
! (1 ) Changed RCS ID tag comment character from "C" to "!$" to 
! allow freeform compilation.  Also convert PARAMETERS to F90 
! syntax. Updated comments, cosmetic changes. (bmy, 6/25/02)
! (2 ) Added atmospheric scale height (7.6 km) (bmy, 6/23/03)
!
!  !REVISION HISTORY:
!  23 Aug 2011 - M. Long   - Converted to Module from Header file
!=================================================================
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_GCTM_MOD
!
! !USES:
!
      IMPLICIT NONE
!      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS: KITCHEN SINK!
!
      PUBLIC
	SAVE
! !DEFINED PARAMETERS: 
!
      ! AIRMW : Molecular weight of air [28.97 g/mole]
      REAL*8, PARAMETER :: AIRMW  =  28.97d0

      ! g0    : Gravity at Surface of Earth [9.8 m/s^2]
      REAL*8, PARAMETER :: g0     =   9.8d0                 

      ! PI    : Double-Precision value of PI          
      REAL*8, PARAMETER :: PI     =   3.14159265358979323d0 
 
      ! Re    : Radius of Earth [m] 
      REAL*8, PARAMETER :: Re     =   6.375d6               

      ! Rd    : Gas Constant (R) in Dry Air [287 J/K/kg] 
      REAL*8, PARAMETER :: Rd     = 287.0d0                 

      ! g0_100 = 100.0 / g0
      REAL*8, PARAMETER :: g0_100 = 100d0 / g0

      ! PI_180 = PI    / 180.0
      REAL*8, PARAMETER :: PI_180 =    PI / 180d0

      ! Rdg0   = Rd    / g0
      REAL*8, PARAMETER :: Rdg0   =    Rd / g0

      ! Scale height of atmosphere (7.6 km = 7600m)
      REAL*8, PARAMETER :: SCALE_HEIGHT = 7600d0
	END MODULE CMN_GCTM_MOD
