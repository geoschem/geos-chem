!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: registry_params_mod.F90
!
! !DESCRIPTION: Contains parameters that are used to denote the types
!  of pointers arrays (e.g. REAL(fp), REAL(f4), INTEGER) used in the
!  GEOS-Chem Registry and History routines, as well as the vertical
!  location.
!\\
!\\
! !INTERFACE:
!
MODULE Registry_Params_Mod
!
! !USES:
!
  USE Precision_Mod

  IMPLICIT NONE
!
! !DEFINED PARAMETERS:
!
  !==========================================================================
  ! Numerical type parameters used in registering & comparing variables
  !==========================================================================
  INTEGER, PUBLIC, PARAMETER :: KINDVAL_I4 = 0           ! 4-byte integer
  INTEGER, PUBLIC, PARAMETER :: KINDVAL_F4 = 1           ! 4-byte real
  INTEGER, PUBLIC, PARAMETER :: KINDVAL_F8 = 2           ! 8-byte real
#ifdef USE_REAL8
  INTEGER, PUBLIC, PARAMETER :: KINDVAL_FP = KINDVAL_F8  ! Flex = 8-byte real
#else
  INTEGER, PUBLIC, PARAMETER :: KINDVAL_FP = KINDVAL_F4  ! Flex = 4-byte real
#endif

  !==========================================================================
  ! Vertical location parameters
  !==========================================================================
  INTEGER, PUBLIC, PARAMETER :: VLocationNone   = 0
  INTEGER, PUBLIC, PARAMETER :: VLocationEdge   = 1
  INTEGER, PUBLIC, PARAMETER :: VLocationCenter = 2
!
! !REVISION HISTORY:
!  14 Jul 2017 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
END MODULE Registry_Params_Mod
!EOC
