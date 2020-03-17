!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: precision_mod.F90
!
! !DESCRIPTION: Module PRECISION\_MOD is used to change the precision of
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
  !=================================================================
  ! Set parameters for floating precision
  !
  ! FP will be set to either 4-byte or 8-byte precision at compile
  ! time.  Most variables can now  declared with REAL(fp).
  !=================================================================
#ifdef USE_REAL8

  ! Use 8-byte floating point precision when asked.
  INTEGER, PARAMETER, PUBLIC :: fp = KIND( REAL( 0.0, 8 ) )

#else

  ! Use 4-byte floating point by default.
  INTEGER, PARAMETER, PUBLIC :: fp = KIND( REAL( 0.0, 4 ) )

#endif

  !=================================================================
  ! Set parameters for fixed precision
  !
  ! Not all variables can be converted into the flexible precision.
  ! Some may have to be still declared as either 4-byte or 8-byte
  ! floating point.  Use these parameters for such variables.
  !=================================================================

  ! KIND parameter for 4-byte precision
  INTEGER, PARAMETER, PUBLIC :: f4 = KIND( REAL( 0.0, 4 ) )

  ! KIND parameter for 8-byte precision
  INTEGER, PARAMETER, PUBLIC :: f8 = KIND( REAL( 0.0, 8 ) )
!
! !REMARKS:
!  This module is designed to help avoid hard-coding precision.
!
! !REVISION HISTORY:
!  04 Nov 2014 - M. Yannetti - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!-----------------------------------------------------------------------------
!BOC
END MODULE PRECISION_MOD
!EOC
