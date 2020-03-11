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
  REAL*8,  PARAMETER :: PTOP = 0.01d0

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

  ! Number of other aerosol categories (rvm, bmy, 2/27/02)
  ! Now set in Init_CMN_SIZE below (mps, 1/3/18)
  INTEGER             :: NAER

  ! NRH -- number of relative humidity bins (rvm, bmy, 2/27/02)
  INTEGER,  PARAMETER :: NRH = 5

  ! Number of dust size bins for transport (tdf, bmy, 3/31/04)
#ifdef TOMAS
# if defined( TOMAS40 )
  INTEGER, PARAMETER :: NDSTBIN = 40
# elif defined( TOMAS15 )
  INTEGER, PARAMETER :: NDSTBIN = 15
# elif defined( TOMAS12 )
  INTEGER, PARAMETER :: NDSTBIN = 12
# else
  INTEGER, PARAMETER :: NDSTBIN = 30
# endif
#else
  INTEGER, PARAMETER :: NDSTBIN = 4
#endif
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Cmn_Size
!
! !DESCRIPTION: Routine INIT\_CMN\_SIZE initializes the grid dimension values
!  in module CMN\_SIZE\_mod.F90.
!\\
!\\
! !INTERFACE:

  SUBROUTINE Init_CMN_SIZE( Input_Opt, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  15 Oct 2012 - M. Long     - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! Set values for former variables declared as parameters
    !=================================================================

    ! Number of aerosol categories
    IF ( Input_Opt%LUCX ) THEN
       ! UCX-based mechanisms include stratospheric aerosols
       NAER = NRHAER + NSTRATAER
    ELSE
       NAER = NRHAER
    ENDIF

    ! Return w/ success
    RC = GC_SUCCESS

  END SUBROUTINE Init_CMN_SIZE
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Cleanup_Cmn_Size
!
! !DESCRIPTION: Subroutine CLEANUP\_CMN\_SIZE deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_CMN_SIZE( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC        ! Success or failure?
!
! !REVISION HISTORY:
!  03 Dec 2012 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Return successfully
    RC = GC_SUCCESS

  END SUBROUTINE Cleanup_CMN_SIZE
!EOC
END MODULE CMN_SIZE_MOD
