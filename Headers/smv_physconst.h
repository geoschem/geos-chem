!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !INCLUDE: smv_physconst.h
!
! !DESCRIPTION: This include file contains physical constants for the
!  GEOS-Chem column chemistry code. 
!\\
!\\
! !DEFINED PARAMETERS: 
!
      ! Molecular weight of air [28.97e-3 kg/mol]
      REAL*8, PARAMETER :: MW_AIR       = 28.97d-3

      ! Avogadro's # [#/mol]
      REAL*8, PARAMETER :: AVO          = 6.022d23

      ! g0    : Gravity at Surface of Earth [9.8 m/s^2]
      REAL*8, PARAMETER :: g0           = 9.8d0                 

      ! PI    : Double-Precision value of PI          
      REAL*8, PARAMETER :: PI           = 3.14159265358979323d0 
 
      ! Re    : Radius of Earth [m] 
      REAL*8, PARAMETER :: Re           = 6.375d6               

      ! Rd    : Gas Constant (R) in Dry Air [287 J/K/kg] 
      REAL*8, PARAMETER :: Rd           = 287.0d0                 

      ! g0_100 = 100.0 / g0
      REAL*8, PARAMETER :: g0_100       = 100d0 / g0

      ! PI_180 = PI    / 180.0
      REAL*8, PARAMETER :: PI_180       = PI / 180d0

      ! Rdg0   = Rd    / g0
      REAL*8, PARAMETER :: Rdg0         = Rd / g0

      ! Scale height of atmosphere (7.6 km = 7600m)
      REAL*8, PARAMETER :: SCALE_HEIGHT = 7600d0

      ! Cp = 1000 J / kg / K = specific heat of air at constant P
      REAL*8, PARAMETER :: Cp           = 1000.0d0

      ! Von Karman's constant
      REAL*8, PARAMETER :: VON_KARMAN   = 0.4d0
!
! !REMARKS:
!  In older sections of code, AIRMW may be replaced by (MW_AIR*1d3).
!
! !REVISION HISTORY: 
!  14 Dec 2009 - R. Yantosca - Initial version, adapted from CMN_GCTM
!EOP
!------------------------------------------------------------------------------
