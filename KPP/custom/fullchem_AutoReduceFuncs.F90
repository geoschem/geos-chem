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
  PUBLIC :: fullchem_AR_SetIntegratorOptions
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
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fullchem_AR_SetIntegratorOptions
!
! !DESCRIPTION: Defines the settings for ICNTRL and RCNTRL used for the
!  rosenbrock_autoreduce integrator.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE fullchem_AR_SetIntegratorOptions( Input_Opt, State_Chm,      &
                                               State_Met, FirstChem,      &
                                               I,         J,         L,   &
                                               ICNTRL,    RCNTRL          )
!
! !USES:
!
    USE gckpp_Parameters
    USE gckpp_Precision
    USE Input_Opt_Mod, ONLY : OptInput
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS: 
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input Options object 
    TYPE(ChmState), INTENT(IN)    :: State_Chm    ! Chemistry State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
    LOGICAL,        INTENT(IN)    :: FirstChem    ! Is it the 1st chem timestep
    INTEGER,        INTENT(IN)    :: I, J, L      ! Grid box indices
!
! !INPUT/OUTPUT PARAMETERS: 
!
    INTEGER,        INTENT(INOUT) :: ICNTRL(20)   ! Options for KPP (integer)
    REAL(dp),       INTENT(INOUT) :: RCNTRL(20)   ! Options for KPP (real   )
!
! !REMARKS:
!  This code was abstracted out of the parallel DO loop in DO_FULLCHEM
!  (in module GeosCore/fullchem_mod.F90)
!EOP
!------------------------------------------------------------------------------
!BOC
    !=====================================================================
    ! fullchem_AR_SetIntegratorOptions begins here!
    !=====================================================================
    
    ! Initialize
    ICNTRL = 0
    RCNTRL = 0.0_dp

    !=====================================================================
    ! Define settings in the ICNTRL vector
    !=====================================================================

    ! 0 - non-autonomous, 1 - autonomous
    ICNTRL(1) = 1

    ! NOTE: ICNTRL is already zeroed out so we don't need to
    ! set this to zero again.  Uncomment if you change this value.
    !! 0 - vector tolerances, 1 - scalars
    !ICNTRL(2) = 0  

    ! Select a particular integration method.
    ! For Rosenbrock, options are:
    ! = 0 :  default method is Rodas3
    ! = 1 :  method is  Ros2
    ! = 2 :  method is  Ros3
    ! = 3 :  method is  Ros4
    ! = 4 :  method is  Rodas3
    ! = 5:   method is  Rodas4
    ICNTRL(3) = 4

    ! 0 - adjoint, 1 - no adjoint
    ICNTRL(7) = 1

    ! Turn off calling Update_SUN, Update_RCONST, Update_PHOTO from within
    ! the integrator.  Rate updates are done before calling KPP.
    ICNTRL(15) = -1

    ! Specify that a threshold value will be used for auto-reduction.
    ! The threshold will be specified in RCNTRL(12), see below.
    IF ( Input_Opt%USE_AUTOREDUCE .and. .not. FIRSTCHEM ) ICNTRL(12) = 1

    ! Use the append functionality?
    IF ( Input_Opt%AUTOREDUCE_IS_APPEND ) ICNTRL(13) = 1

    !=====================================================================
    ! Define settings in the ICNTRL vector
    !=====================================================================

    ! Initialize Hstart (the starting value of the integration step
    ! size with the value of Hnew (the last predicted but not yet 
    ! taken timestep)  saved to the the restart file.
    RCNTRL(3) = State_Chm%KPPHvalue(I,J,L)

    !---------------------------------------------------------------------
    ! Auto-reduce threshold, Method 1: Pressure-dependent
    !                                            
    !   Actual_Threshold =
    !                                           Mid-Pressure at Level
    !     AUTOREDUCE_THRESHOLD (at surface) * --------------------------
    !                                          "Mid-Pressure" at Sfc.
    !
    !---------------------------------------------------------------------
    IF ( .not. Input_Opt%AUTOREDUCE_IS_KEY_THRESHOLD ) THEN
       IF ( Input_Opt%AUTOREDUCE_IS_PRS_THRESHOLD ) THEN
          RCNTRL(12) = Input_Opt%AUTOREDUCE_THRESHOLD                        & 
                     * State_Met%PMID(I,J,L)                                 & 
                     / State_Met%PMID(I,J,1)
       ENDIF
       
       IF ( .not. Input_Opt%AUTOREDUCE_IS_PRS_THRESHOLD ) THEN
          RCNTRL(12) = Input_Opt%AUTOREDUCE_THRESHOLD
       ENDIF
    ENDIF

    !---------------------------------------------------------------------
    ! Auto-reduce threshold, Method 2: Determine threshold 
    ! dynamically by scaling rates of key species.
    !---------------------------------------------------------------------
    IF ( Input_Opt%AUTOREDUCE_IS_KEY_THRESHOLD ) THEN

       !--------------------------------
       ! Daytime target species (OH)
       !--------------------------------
       ICNTRL(14) = ind_OH
       RCNTRL(14) = Input_Opt%AUTOREDUCE_TUNING_OH
       
       !--------------------------------
       ! Nighttime target species (NO2)
       !--------------------------------
       ! COMMENTS BY HAIPENG LIN:
       ! 1e6 daytime conc...testing shows 5e-5 as an offset here works best.
       ! Use JNO2 as night determination.
       ! RXN_NO2: NO2 + hv --> NO  + O
       ! JNO2 ranges from 0 to 0.02 and is order ~ 1e-4 at the terminator. 
       ! We set this threshold to be slightly relaxed so it captures the 
       ! terminator, but this needs some tweaking.
       !
       ! For some reason, RXN_NO2 as a proxy fails to propagate the sunset 
       ! terminator even though all diagnostics seem fine, and after a while 
       ! only the OH scheme applies.  Use SUNCOSmid as a proxy to fix this. 
       ! (hplin, 4/20/22)
       ! IF(ZPJ(L,RXN_NO2,I,J) .eq. 0.0_fp) THEN
       !
       IF( State_Met%SUNCOSmid(I,J) .le. -0.1391731e+0_dp ) THEN
          ICNTRL(14) = ind_NO2
          RCNTRL(14) = Input_Opt%AUTOREDUCE_TUNING_NO2
       ENDIF
    ENDIF

  END SUBROUTINE fullchem_AR_SetIntegratorOptions
!EOC
END MODULE fullchem_AutoReduceFuncs
