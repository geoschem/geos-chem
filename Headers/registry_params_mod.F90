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
  ! Denotes REAL(fp), aka "flexible precision
  INTEGER, PUBLIC, PARAMETER :: KINDVAL_FP = 1

  ! Denotes REAL(f8), aka REAL*8 or 8-byte floating point
  INTEGER, PUBLIC, PARAMETER :: KINDVAL_F8 = 2

  ! Denotes REAL(f4), aka REAL*4, or 4-byte floating point
  INTEGER, PUBLIC, PARAMETER :: KINDVAL_F4 = 3

  ! Denotes INTEGER
  INTEGER, PUBLIC, PARAMETER :: KINDVAL_I4 = 4

  ! Vertical locations
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
