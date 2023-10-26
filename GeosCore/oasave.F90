#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: oasave.F90
!
! !DESCRIPTION: Subroutine OASAVE stores the concentrations of organic aerosols
!  for the ND42 diagnostic and for the timeseries and satellite diagnostics
!\\
!\\
! !INTERFACE:
!
SUBROUTINE OASAVE( State_Grid, State_Chm )
!
! !USES:
!
  USE State_Chm_Mod,      ONLY : ChmState
  USE State_Grid_Mod,     ONLY : GrdState
  USE PRECISION_MOD            ! For GEOS-Chem Precision (fp)
  USE CMN_O3_MOD,         ONLY : SAVEOA

  IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
  TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
  TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  10 Jul 2014 - E. A. Marais- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  INTEGER  :: I, J, L
  REAL(fp) :: FACTOR

  !=================================================================
  ! OASAVE begins here!
  !
  ! Save organic aerosol concentrations
  !=================================================================

  ! Conversion factor from kg/m3 --> ug/m3
  FACTOR = 1e+9_fp

  !$OMP PARALLEL DO       &
  !$OMP DEFAULT( SHARED ) &
  !$OMP PRIVATE( I, J, L )
  DO L = 1, State_Grid%NZ
  DO J = 1, State_Grid%NY
  DO I = 1, State_Grid%NX

     ! Sum of all organic aerosol [ug/m3]
     SAVEOA(I,J,L) = SAVEOA(I,J,L) + ( State_Chm%AerMass%TSOA(I,J,L) + &
                                       State_Chm%AerMass%ASOA(I,J,L) + &
                                       State_Chm%AerMass%OCPO(I,J,L) + &
                                       State_Chm%AerMass%OCPI(I,J,L) + &
                                       State_Chm%AerMass%OPOA(I,J,L) + &
                                       State_Chm%AerMass%ISOAAQ(I,J,L) ) * FACTOR

  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
END SUBROUTINE OASAVE
!EOC
#endif

