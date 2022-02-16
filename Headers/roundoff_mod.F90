!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: roundoff_mod.F90
!
! !DESCRIPTION: Contains routines to round floating point values to a
!  given number of decimal places.
!\\
!\\
! !INTERFACE:
!
MODULE Roundoff_Mod
!
! !USES:
!
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Cast_and_RoundOff
  PUBLIC :: Roundoff
  INTERFACE RoundOff
     MODULE PROCEDURE RoundOff_F4
     MODULE PROCEDURE RoundOff_F8
  END INTERFACE
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: Roundoff_F4
  PRIVATE :: Roundoff_F8
!
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RoundOff_f4
!
! !DESCRIPTION: Rounds a number X to N decimal places of precision.
!\\
!\\
! !INTERFACE:
!
  FUNCTION RoundOff_f4( X, N ) RESULT( Y )
!
! !INPUT PARAMETERS:
!
    REAL(f4), INTENT(IN) :: X   ! Number to be rounded
    INTEGER,  INTENT(IN) :: N   ! Number of decimal places to keep
!
! !RETURN VALUE:
!
    REAL(f4)             :: Y   ! Number rounded to N decimal places
!
! !REMARKS:
!  The algorithm to round X to N decimal places is as follows:
!  (1) Multiply X by 10**(N+1)
!  (2) If X < 0, then add -5 to X; otherwise add 5 to X
!  (3) Round X to nearest integer
!  (4) Divide X by 10**(N+1)
!  (5) Truncate X to N decimal places: INT( X * 10**N ) / 10**N
!                                                                             .
!  Rounding algorithm from: Hultquist, P.F, "Numerical Methods for Engineers
!   and Computer Scientists", Benjamin/Cummings, Menlo Park CA, 1988, p. 20.
!                                                                             .
!  Truncation algorithm from: http://en.wikipedia.org/wiki/Truncation
!                                                                             .
!  The two algorithms have been merged together for efficiency.
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Round and truncate X to N decimal places
    Y = INT( NINT( X*(10.0_f4**(N+1)) + SIGN( 5.0_f4, X ) ) / 10.0_f4 ) / (10.0_f4**N)

  END FUNCTION RoundOff_f4
!EOC
!------------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1 and      !
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: RoundOff_f8
!
! !DESCRIPTION: Rounds a number X to N decimal places of precision.
!\\
!\\
! !INTERFACE:
!
  FUNCTION RoundOff_f8( X, N ) RESULT( Y )
!
! !INPUT PARAMETERS:
!
    REAL(f8), INTENT(IN) :: X   ! Number to be rounded
    INTEGER,  INTENT(IN) :: N   ! Number of decimal places to keep
!
! !RETURN VALUE:
!
    REAL(f8)             :: Y   ! Number rounded to N decimal places
!
! !REMARKS:
!  The algorithm to round X to N decimal places is as follows:
!  (1) Multiply X by 10**(N+1)
!  (2) If X < 0, then add -5 to X; otherwise add 5 to X
!  (3) Round X to nearest integer
!  (4) Divide X by 10**(N+1)
!  (5) Truncate X to N decimal places: INT( X * 10**N ) / 10**N
!                                                                             .
!  Rounding algorithm from: Hultquist, P.F, "Numerical Methods for Engineers
!   and Computer Scientists", Benjamin/Cummings, Menlo Park CA, 1988, p. 20.
!                                                                             .
!  Truncation algorithm from: http://en.wikipedia.org/wiki/Truncation
!                                                                             .
!  The two algorithms have been merged together for efficiency.
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Round and truncate X to N decimal places
    Y = INT( NINT( X*(10.0_f8**(N+1)) + SIGN( 5.0_f8, X ) ) / 10.0_f8 ) / (10.0_f8**N)

  END FUNCTION RoundOff_f8
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cast_and_roundoff
!
! !DESCRIPTION: Casts a 4-byte variable to 8-byte, and then rounds off
!  to 2 decimal places.  Used for species database fields.
!\\
!\\
! !INTERFACE:
!
  FUNCTION Cast_and_RoundOff( v_real, places ) RESULT( v_dble )
!
! !INPUT PARAMETERS:
!
    REAL(f4), INTENT(IN) :: v_real   ! Input, 4-byte real
    INTEGER,  INTENT(IN) :: places   ! Keep this many places
!
! !RETURN VALUE:
!
    REAL(f8)             :: v_dble   ! Output, 8-byte real
!EOP
!------------------------------------------------------------------------------
!BOC

    ! If v_real is a missing value, return with 8-byte missing value
    ! Assume missing value is -999.0
    IF ( v_real < -900.0_f4 ) THEN
       v_dble = -999.0_fp
       RETURN
    ENDIF

    ! Cast to real*8 and roundoff (if the number isn't too large)
    v_dble = RoundOff( DBLE( v_real ), places )

  END FUNCTION Cast_And_RoundOff
!EOC
END MODULE Roundoff_Mod
