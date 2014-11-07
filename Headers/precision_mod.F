!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: precision_mod
!
! !DESCRIPTION: Module precision\_mod is used to change the precision of
!  many variables throughout GEOS-Chem at compile-time.
!\\
!\\
! !INTERFACE:
!
      MODULE PRECISION_MOD
!
! !USES:
! 
      IMPLICIT NONE
      PRIVATE
!
! !DEFINED PARAMETERS:
! 
#if defined( USE_REAL8 )

      ! Use 8-byte floating point precision when asked.
      INTEGER, PARAMETER, PUBLIC :: fp = KIND( 0.d0 )

#else

      ! Use 4-byte floating point by default.
      INTEGER, PARAMETER, PUBLIC :: fp = KIND( 0.0 )

#endif

! !REMARKS:
!  This module is designed to help avoid hard-coding precision.
!
! !REVISION HISTORY:
!  (1 ) Created. (myannetti, 11/04/14)

      END MODULE PRECISION_MOD
