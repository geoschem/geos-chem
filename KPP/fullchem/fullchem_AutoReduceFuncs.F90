!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fullchem_AutoReduceFuncs
!
! !DESCRIPTION: Contains routines to abstract rosenbrock_autoreduce-specific
!  code out of fullchem_mod.F90.  This will avoid compilation errors when
!  building other KPP-generated mechanisms.
!\\
!\\
! !INTERFACE:
!
MODULE fullchem_AutoReduceFuncs
!
! !USES:
!
  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: fullchem_AR_KeepHalogensActive
  PUBLIC :: fullchem_AR_SetKeepActive
  PUBLIC :: fullchem_AR_UpdateKppDiags
!
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
! !IROUTINE: fullchem_AR_KeepHalogensActive
!
! !DESCRIPTION: Sets halogen species to "fast" for the rosenbrock_autoreduce 
!  integrator.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE fullchem_AR_KeepHalogensActive( doPrint )
!
! !USES:
!
    USE gckpp_Precision
    USE gckpp_Parameters
    USE gckpp_Global,    ONLY : keepSpcActive
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN) :: doPrint   ! Print informational message
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Informational printout
    IF ( doPrint ) THEN
       WRITE( 6, 100 )
 100   FORMAT('Setting halogen species to "fast" for rosenbrock_autoreduce"')
    ENDIF

    ! New halogens auto-reduce list hplin 01/27/22, 03/02/22
    ! based off Shen et al. 2020 GMD Table 1, lines 10, 11, 12
    keepSpcActive(ind_AERI  ) = .TRUE.  ! Iodine on aerosol
    keepSpcActive(ind_Br    ) = .TRUE.
    keepSpcActive(ind_Br2   ) = .TRUE.
    keepSpcActive(ind_BrCl  ) = .TRUE.
    keepSpcActive(ind_BrNO2 ) = .TRUE.
    keepSpcActive(ind_BrNO3 ) = .TRUE.
    keepSpcActive(ind_BrO   ) = .TRUE.
    keepSpcActive(ind_BrSALA) = .TRUE.
    keepSpcActive(ind_BrSALC) = .TRUE.
    keepSpcActive(ind_HBr   ) = .TRUE.
    keepSpcActive(ind_HOBr  ) = .TRUE.
    keepSpcActive(ind_Cl    ) = .TRUE.
    keepSpcActive(ind_Cl2   ) = .TRUE.
    keepSpcActive(ind_Cl2O2 ) = .TRUE.
    keepSpcActive(ind_ClNO2 ) = .TRUE.
    keepSpcActive(ind_ClNO3 ) = .TRUE.
    keepSpcActive(ind_ClO   ) = .TRUE.
    keepSpcActive(ind_ClOO  ) = .TRUE.
    keepSpcActive(ind_OClO  ) = .TRUE.
    keepSpcActive(ind_HCl   ) = .TRUE.
    keepSpcActive(ind_HOCl  ) = .TRUE.
    keepSpcActive(ind_I     ) = .TRUE.
    keepSpcActive(ind_I2    ) = .TRUE.
    keepSpcActive(ind_IO    ) = .TRUE.
    keepSpcActive(ind_I2O2  ) = .TRUE.
    keepSpcActive(ind_HI    ) = .TRUE.
    keepSpcActive(ind_ISALA ) = .TRUE.
    keepSpcActive(ind_ISALC ) = .TRUE.
    keepSpcActive(ind_I2O4  ) = .TRUE.
    keepSpcActive(ind_I2O3  ) = .TRUE.
    keepSpcActive(ind_INO   ) = .TRUE.
    keepSpcActive(ind_IONO  ) = .TRUE.
    keepSpcActive(ind_IONO2 ) = .TRUE.
    keepSpcActive(ind_ICl   ) = .TRUE.
    keepSpcActive(ind_IBr   ) = .TRUE.
    keepSpcActive(ind_HOI   ) = .TRUE.
    keepSpcActive(ind_SALACl) = .TRUE.
    keepSpcActive(ind_SALCCl) = .TRUE.
    keepSpcActive(ind_SALAAL) = .TRUE.
    keepSpcActive(ind_SALCAL) = .TRUE.

  END SUBROUTINE fullchem_AR_KeepHalogensActive
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fullchem_AR_SetKeepActive
!
! !DESCRIPTION:  Abstracts setting the rosenbrock_autoreduce keepActive flag
!  out of fullchem_mod.F90
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE fullchem_AR_SetKeepActive( option )
!
! !USES:
!
    USE gckpp_Precision
    USE gckpp_Global, ONLY : keepActive
!
! !INPUT PARAMETERS: 
!
    LOGICAL, INTENT(IN) :: option
!EOP
!------------------------------------------------------------------------------
!BOC
    keepActive = option

  END SUBROUTINE fullchem_AR_SetKeepActive
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fullchem_AR_UpdateKppDiags
!
! !DESCRIPTION: Updates KPP diagnostics for the rosenbrock_autoreduce solver
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE fullchem_AR_UpdateKppDiags( I, J, L, RSTATE, State_Diag )
!
! !USES:
!
    USE gckpp_Precision
    USE gckpp_Global,     ONLY : cNONZERO, rNVAR
    USE gckpp_Integrator, ONLY : NARthr
    USE State_Diag_Mod,   ONLY : DgnState
!
! !INPUT PARAMETERS: 
!
    INTEGER,        INTENT(IN)    :: I, J, L
    REAL(dp),       INTENT(IN)    :: RSTATE(20)
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag
!EOP
!------------------------------------------------------------------------------
!BOC
    ! # of species in auto-reduced mechanism
    IF ( State_Diag%Archive_KppAutoReducerNVAR ) THEN
       State_Diag%KppAutoReducerNVAR(I,J,L) = rNVAR
    ENDIF

    ! Computed threshold
    IF ( State_Diag%Archive_KppAutoReduceThres ) THEN
       State_Diag%KppAutoReduceThres(I,J,L) = RSTATE(NARthr)
    ENDIF

    ! # of nonzero elements in LU factorization of Jacobian, AR only
    IF ( State_Diag%Archive_KppcNONZERO ) THEN
       State_Diag%KppcNONZERO(I,J,L) = cNONZERO
    ENDIF

  END SUBROUTINE fullchem_AR_UpdateKppDiags
!EOC
END MODULE fullchem_AutoReduceFuncs
