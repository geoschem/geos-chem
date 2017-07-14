!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: registry_params_mod.F90
!
! !DESCRIPTION: Contains parameters that are used to denote the types
!  of pointers arrays (e.g. REAL(fp), REAL(f4), INTEGER) used in the 
!  GEOS-Chem Registry and History routines.
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

  ! Denotes REAL(f4), aka REAL*4, or 4-byte floating point
  INTEGER, PUBLIC, PARAMETER :: KINDVAL_F4 = 2

  ! Denotes INTEGER
  INTEGER, PUBLIC, PARAMETER :: KINDVAL_I4 = 3
!
! !REVISION HISTORY:
!  14 Jul 2017 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
END MODULE Registry_Params_Mod
!EOC
