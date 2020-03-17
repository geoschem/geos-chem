#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diag_mod.F90
!
! !DESCRIPTION: Module DIAG\_MOD contains declarations for allocatable arrays
!  for use with GEOS-CHEM diagnostics.
!\\
!\\
! !INTERFACE:
!
MODULE DIAG_MOD
!
! !USES:
!
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PUBLIC
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: CLEANUP_DIAG
!
! !PUBLIC DATA MEMBERS:
!
#ifdef TOMAS
  ! For ND44 -- Dry deposition fluxes & velocities
  REAL*4,  ALLOCATABLE :: AD44     (:,:,:,:)
  REAL*4,  ALLOCATABLE :: AD59_NUMB(:,:,:,:)
  REAL*4,  ALLOCATABLE :: AD59_SULF(:,:,:,:)
  REAL*4,  ALLOCATABLE :: AD59_SALT(:,:,:,:)
  REAL*4,  ALLOCATABLE :: AD59_ECIL(:,:,:,:)
  REAL*4,  ALLOCATABLE :: AD59_ECOB(:,:,:,:)
  REAL*4,  ALLOCATABLE :: AD59_OCIL(:,:,:,:)
  REAL*4,  ALLOCATABLE :: AD59_OCOB(:,:,:,:)
  REAL*4,  ALLOCATABLE :: AD59_DUST(:,:,:,:)

  ! For ND60 -- TOMAS condensation rate diagnostic
  REAL*4,  ALLOCATABLE :: AD60_COND (:,:,:,:)
  REAL*4,  ALLOCATABLE :: AD60_COAG (:,:,:,:)
  REAL*4,  ALLOCATABLE :: AD60_NUCL (:,:,:,:)
  REAL*4,  ALLOCATABLE :: AD60_AQOX (:,:,:,:)
  REAL*4,  ALLOCATABLE :: AD60_ERROR(:,:,:,:)
  REAL*4,  ALLOCATABLE :: AD60_SOA  (:,:,:,:)

  ! For ND61 -- 3D TOMAS rate diagnostic
  REAL*4,  ALLOCATABLE :: AD61(:,:,:,:)
  REAL*4,  ALLOCATABLE :: AD61_inst(:,:,:,:)

  ! Prod/loss rates
  REAL*4,  ALLOCATABLE :: AD65(:,:,:,:)
#endif
#ifdef RRTMG
  ! For ND72 -- Radiation output diagnostic
  REAL*4,  ALLOCATABLE :: AD72(:,:,:)
#endif
!
! !REVISION HISTORY:
!  30 Nov 1999 - A. Fiore - Initial version
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
! !IROUTINE: cleanup_diag
!
! !DESCRIPTION: Subroutine CLEANUP\_DIAG deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_DIAG()
!
! !REVISION HISTORY:
!  13 Dec 2002 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
    !=================================================================
    ! CLEANUP_DIAG begins here!
    !=================================================================
#ifdef TOMAS
    IF ( ALLOCATED( AD44        ) ) DEALLOCATE( AD44        )
    IF ( ALLOCATED( AD59_NUMB   ) ) DEALLOCATE( AD59_NUMB   )
    IF ( ALLOCATED( AD59_SULF   ) ) DEALLOCATE( AD59_SULF   )
    IF ( ALLOCATED( AD59_SALT   ) ) DEALLOCATE( AD59_SALT   )
    IF ( ALLOCATED( AD59_ECIL   ) ) DEALLOCATE( AD59_ECIL   )
    IF ( ALLOCATED( AD59_ECOB   ) ) DEALLOCATE( AD59_ECOB   )
    IF ( ALLOCATED( AD59_OCIL   ) ) DEALLOCATE( AD59_OCIL   )
    IF ( ALLOCATED( AD59_OCOB   ) ) DEALLOCATE( AD59_OCOB   )
    IF ( ALLOCATED( AD59_DUST   ) ) DEALLOCATE( AD59_DUST   )
    IF ( ALLOCATED( AD60_COND   ) ) DEALLOCATE( AD60_COND   )
    IF ( ALLOCATED( AD60_COAG   ) ) DEALLOCATE( AD60_COAG   )
    IF ( ALLOCATED( AD60_NUCL   ) ) DEALLOCATE( AD60_NUCL   )
    IF ( ALLOCATED( AD60_AQOX   ) ) DEALLOCATE( AD60_AQOX   )
    IF ( ALLOCATED( AD60_ERROR  ) ) DEALLOCATE( AD60_ERROR  )
    IF ( ALLOCATED( AD60_SOA    ) ) DEALLOCATE( AD60_SOA    )
    IF ( ALLOCATED( AD61        ) ) DEALLOCATE( AD61        )
    IF ( ALLOCATED( AD61_inst   ) ) DEALLOCATE( AD61_inst   )
    IF ( ALLOCATED( AD65        ) ) DEALLOCATE( AD65        )
#endif
#ifdef RRTMG
    IF ( ALLOCATED( AD72        ) ) DEALLOCATE( AD72        )
#endif
  END SUBROUTINE CLEANUP_DIAG
!EOC
END MODULE DIAG_MOD
#endif
