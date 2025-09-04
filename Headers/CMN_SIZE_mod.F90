!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_SIZE_mod.F90
!
! !DESCRIPTION: CMN\_SIZE contains size parameters for GEOS-Chem arrays.
!\\
!\\
! !INTERFACE:
!
MODULE CMN_SIZE_MOD
!
! !USES:
!
  IMPLICIT NONE
  PUBLIC
!
! !DEFINED PARAMETERS:
!
  ! Model top pressure (mb)
#if defined( MODEL_CESM )
  REAL*8             :: PTOP
#else
  REAL*8,  PARAMETER :: PTOP = 0.01d0
#endif

  ! Maximum number of surface types: 73 olson
  INTEGER, PARAMETER :: NSURFTYPE = 73

  ! Maximum number of veg types in a CTM grid box
#if defined( EXTERNAL_GRID ) || defined( EXTERNAL_TYPE )
  INTEGER, PARAMETER :: NTYPE = 50
#else
  INTEGER, PARAMETER :: NTYPE = 25
#endif

  !Number of coefficients for polynomial fits
  INTEGER, PARAMETER :: NPOLY = 20

  ! Number of FAST-J aerosol size bins (rvm, bmy, 11/15/01)
  INTEGER, PARAMETER :: NDUST = 7

  ! Number of aerosols undergoing hygroscopic growth
  INTEGER, PARAMETER :: NRHAER = 5

  ! Number of stratospheric aerosols (SDE 04/17/13)
  INTEGER, PARAMETER :: NSTRATAER = 2

  ! Number of other aerosol categories, include stratospheric aerosols
  INTEGER, PARAMETER :: NAER = NRHAER + NSTRATAER

  ! NRH -- number of relative humidity bins (rvm, bmy, 2/27/02)
  INTEGER,  PARAMETER :: NRH = 5

  ! Number of dust size bins for transport (tdf, bmy, 3/31/04)
#if defined(TOMAS)
#if defined(TOMAS40)
  INTEGER, PARAMETER :: NDSTBIN = 40
#elif defined(TOMAS15)
  INTEGER, PARAMETER :: NDSTBIN = 15
#elif defined(TOMAS12)
  INTEGER, PARAMETER :: NDSTBIN = 12
#else
  INTEGER, PARAMETER :: NDSTBIN = 30
#endif
#else
  INTEGER, PARAMETER :: NDSTBIN = 4
#endif
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!EOC
END MODULE CMN_SIZE_MOD
