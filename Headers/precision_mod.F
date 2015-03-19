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
      !=================================================================
      ! Set parameters for floating precision
      !
      ! FP will be set to either 4-byte or 8-byte precision at compile 
      ! time.  Most variables can now  declared with REAL(fp).
      !=================================================================
#if defined( USE_REAL8 )

      ! Use 8-byte floating point precision when asked.
      INTEGER, PARAMETER, PUBLIC :: fp = KIND( 0.d0 )

#else

      ! Use 4-byte floating point by default.
      INTEGER, PARAMETER, PUBLIC :: fp = KIND( 0.0 )

#endif

      !=================================================================
      ! Set parameters for fixed precision
      !
      ! Not all variables can be converted into the flexible precision.  
      ! Some may have to be still declared as either 4-byte or 8-byte 
      ! floating point.  Use these parameters for such variables.
      !=================================================================

      ! KIND parameter for 4-byte precision
      INTEGER, PARAMETER, PUBLIC :: f4 = KIND( 0.0 )
      
      ! KIND parameter for 8-byte precision
      INTEGER, PARAMETER, PUBLIC :: f8 = KIND( 0.d0 )
!
! !REMARKS:
!  This module is designed to help avoid hard-coding precision.
!
! !REVISION HISTORY:
!  (1 ) Created. (myannetti, 11/04/14)
!EOP
!-----------------------------------------------------------------------------
!BOC
      END MODULE PRECISION_MOD
!EOC
