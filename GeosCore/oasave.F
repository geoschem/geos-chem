!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: oasave.F
!
! !DESCRIPTION: Subroutine OASAVE stores the concentrations of organic aerosols
!  for the ND42 diagnostic and for the timeseries and satellite diagnostics
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE OASAVE( SAVEOA,    am_I_Root, Input_Opt,
     &                   State_Met, State_Chm, RC )
!
! !USES:
!
      USE ErrCode_Mod
      USE CMN_SIZE_MOD
      USE PRECISION_MOD                            ! For GEOS-Chem Precision (fp)

      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
      USE State_Met_Mod,      ONLY : MetState
      
      IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
      REAL(fp), ALLOCATABLE, INTENT(INOUT) :: SAVEOA(:,:,:)
      LOGICAL,        INTENT(IN)           :: am_I_Root   ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN)           :: Input_Opt   ! Input Options object
      TYPE(MetState), INTENT(IN)           :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS: 
!
      INTEGER,        INTENT(OUT)          :: RC          ! Success or failure
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState), INTENT(INOUT)        :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY: 
!  10 Jul 2014 - E. A. Marais- Initial version
!  09 Jun 2017 - M. Sulprizio- Modified for consistency with ND42 updates made
!                              in v11-01
!  08 Aug 2018 - H.P. Lin    - Accept arguments from carbon_mod.F and
!                              Moved aerosol in-module variables to state_chm_mod
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

      ! Initialize
      RC = GC_SUCCESS

      ! Conversion factor from kg/m3 --> ug/m3
      FACTOR = 1e+9_fp

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Sum of all organic aerosol [ug/m3]
         SAVEOA(I,J,L) = SAVEOA(I,J,L) +
     &                   ( State_Chm%TSOA(I,J,L) + 
     &                     State_Chm%ISOA(I,J,L) + 
     &                     State_Chm%ASOA(I,J,L) +
     &                     State_Chm%OCPO(I,J,L) + 
     &                     State_Chm%OCPI(I,J,L) + 
     &                     State_Chm%OPOA(I,J,L) +
     &                     State_Chm%ISOAAQ(I,J,L) ) * FACTOR

      ENDDO
      ENDDO
      ENDDO

!$OMP END PARALLEL DO

      END SUBROUTINE OASAVE
!EOC
