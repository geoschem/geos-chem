!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ccycle_RateLawFuncs
!
! !DESCRIPTION: Provides rate-law functions used by the "ccycle" chemical
!  mechanism.  This will be referenced from within subroutine Update_RCONST.
!\\
!\\
! !INTERFACE:
!
MODULE ccycle_RateLawFuncs
!
! !USES:
!
  USE gckpp_Global
  USE gckpp_Parameters
  USE gckpp_Precision
  USE rateLawUtilFuncs

  IMPLICIT NONE

CONTAINS

  FUNCTION GC_OHCO() RESULT( k )
    ! Reaction rate for:
    !    OH + CO = HO2 + CO2 (cf. JPL 15-10)
    !
    ! For this reaction, these Arrhenius law terms evaluate to 1:
    !    (300/T)**b0 * EXP(c0/T)
    ! because b0 = c0 = 0.  Therefore we can skip computing these
    ! terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    !
    REAL(dp)             :: klo1,   klo2,   khi1,  khi2
    REAL(dp)             :: xyrat1, xyrat2, blog1, blog2,   fexp1
    REAL(dp)             :: fexp2,  kco1,   kco2,  TEMP300, k
    !
    klo1   = 5.9E-33_dp * K300_OVER_TEMP
    khi1   = 1.1E-12_dp * K300_OVER_TEMP**(-1.3_dp)
    xyrat1 = klo1 * NUMDEN / khi1
    blog1  = LOG10( xyrat1 )
    fexp1  = 1.0_dp / ( 1.0_dp + blog1*blog1 )
    kco1   = klo1 * NUMDEN * 0.6_dp**fexp1 / ( 1.0_dp + xyrat1 )
    klo2   = 1.5E-13_dp
    khi2   = 2.1E+09_dp * K300_OVER_TEMP**(-6.1_dp)
    xyrat2 = klo2 * NUMDEN / khi2
    blog2  = LOG10( xyrat2 )
    fexp2  = 1.0_dp / ( 1.0_dp + blog2*blog2 )
    kco2   = klo2 * 0.6_dp**fexp2 / ( 1.0_dp + xyrat2 )
    k      = kco1 + kco2
  END FUNCTION GC_OHCO

END MODULE ccycle_RateLawFuncs
