#ifdef BPCH_DIAG
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
SUBROUTINE OASAVE( State_Grid )
!
! !USES:
!
  USE State_Grid_Mod,     ONLY : GrdState
  USE AEROSOL_MOD,        ONLY : OCPI, OCPO
  USE AEROSOL_MOD,        ONLY : TSOA, ASOA, OPOA, ISOAAQ
  USE PRECISION_MOD            ! For GEOS-Chem Precision (fp)
  USE CMN_O3_MOD,         ONLY : SAVEOA

  IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
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
     SAVEOA(I,J,L) = SAVEOA(I,J,L) + ( TSOA(I,J,L) + ASOA(I,J,L) + &
                                       OCPO(I,J,L) + OCPI(I,J,L) + &
                                       OPOA(I,J,L) + ISOAAQ(I,J,L) ) * FACTOR

  ENDDO
  ENDDO
  ENDDO
  !$OMP END PARALLEL DO
END SUBROUTINE OASAVE
!EOC
#endif

