MODULE Hg_RateLawFuncs

  USE gckpp_Global
  USE gckpp_Parameters
  USE gckpp_Precision
  USE rateLawUtilFuncs

  IMPLICIT NONE
  PUBLIC

CONTAINS

  FUNCTION GCJPLPR_abab( a1, b1, a2, b2, fv ) RESULT( k )
    ! Third body effect for pressure dependence of rate coefficients.
    ! a1, b1 are the Arrhenius parameters for the lower-limit rate.
    ! a2, b2 are the Arrhenius parameters for the upper-limit rate.
    ! fv     is the falloff curve paramter, (see ATKINSON ET. AL (1992)
    !        J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
    !
    ! For these reactions, these Arrhenius law terms evaluate to 1:
    !    EXP(c1/T)
    !    EXP(c2/T)
    ! because c1 = c2 = 0.  Therefore we can skip computing these
    ! terms.  Also, fct1 = fct2 = 0, so we will skip computing these
    ! terms as well.  This is more computationally efficient.
    ! (bmy, 06 Jan 2022)
    !
    REAL(dp), INTENT(IN) :: a1,   b1,    a2,    b2,   fv
    REAL(dp)             :: rlow, rhigh, xyrat, blog, fexp, k
    !
    rlow  = a1 * ( K300_OVER_TEMP**b1 ) * NUMDEN
    rhigh = a2 * ( K300_OVER_TEMP**b2 )
    xyrat = rlow / rhigh
    blog  = LOG10( xyrat )
    fexp  = 1.0_dp / ( 1.0_dp + ( blog * blog ) )
    k     = rlow * ( fv**fexp ) / ( 1.0_dp + xyrat )
  END FUNCTION GCJPLPR_abab

END MODULE Hg_RateLawFuncs
