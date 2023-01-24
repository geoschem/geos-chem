!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: carbon_Funcs
!
! !DESCRIPTION: Module containing routines for passing data from GEOS-Chem
!  to the carbon chemical mechanism solver code.
!\\
!\\
! !INTERFACE:
!
MODULE carbon_Funcs
!
! !USES:
!
  USE gckpp_Precision
  USE gckpp_Parameters
  USE gckpp_Global
  USE Precision_Mod,   ONLY : fp
  USE rateLawUtilFuncs
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
! !IROUTINE: carbon_ConvertKgToMolecCm3
!
! !DESCRIPTION: Converts species from kg to molec/cm3 and stores into
!  the "C" concentration array used by the KPP-generated solver code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE carbon_ConvertKgToMolecCm3(                                &
             I,         J,          L,          id_CH4,     id_CO2,          &
             id_CO,     xnumol_CH4, xnumol_CO,  xnumol_CO2, State_Chm,       &
             State_Met                                                      )
!
! !USES:
!
    USE Species_Mod,   ONLY : SpcConc
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: I, J, L     ! Grid box indices
    INTEGER,        INTENT(IN) :: id_CH4      ! Species index for CH4
    INTEGER,        INTENT(IN) :: id_CO       ! Species index for CO
    INTEGER,        INTENT(IN) :: id_CO2      ! Species index for CO2
    REAL(fp),       INTENT(IN) :: xnumol_CH4  ! kg CH4 / molec CH4
    REAL(fp),       INTENT(IN) :: xnumol_CO   ! kg CO  / molec CO
    REAL(fp),       INTENT(IN) :: xnumol_CO2  ! kg CO2 / molec CO2
    TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meterorology State object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL(fp)               :: airvol_cm3

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    !========================================================================
    ! carbon_ConvertKgToMolecCm3 begins here!
    !========================================================================

    ! Point to species array
    Spc => State_Chm%Species

    ! Grid box volume [cm3]
    airvol_cm3 = State_Met%AIRVOL(I,J,L) * 1.0e+6_fp

    ! Convert species kg to molec/cm3 and store in the C array
    ! (skip if species are not present)
    IF ( id_CH4 > 0 ) THEN
       C(ind_CH4) = Spc(id_CH4)%Conc(I,J,L) * xnumol_CH4 / airvol_cm3
    ENDIF

    IF ( id_CO > 0 ) THEN
       C(ind_CO)  = Spc(id_CO)%Conc(I,J,L)  * xnumol_CO  / airvol_cm3
    ENDIF

    IF ( id_CO2 > 0 ) THEN
       C(ind_CO2) = Spc(id_CO2)%Conc(I,J,L) * xnumol_CO2 / airvol_cm3
    ENDIF

    ! Initialize placeholder species to 1 molec/cm3
    C(ind_DummyCH4)   = 1.0_dp
    C(ind_DummyNMVOC) = 1.0_dp

    ! Initialize fixed species to 1 molec/cm3
    ! These will later be set to values read via HEMCO
    C(ind_FixedCl) = 1.0_dp
    C(ind_FixedOH) = 1.0_dp

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE carbon_ConvertKgToMolecCm3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: carbon_ComputeRateConstants
!
! !DESCRIPTION: Computes the rate constants used in the carbon chemical
!  mechanism.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE carbon_ComputeRateConstants(                               &
             I,             J,                 L,                            &
             ConcClMnd,     ConcOHMnd,         LCH4_by_OH,                   &
             LCO_in_Strat,  OHdiurnalFac,      PCO_fr_CH4_use,               &
             PCO_fr_CH4,    PCO_fr_NMVOC_use,  PCO_fr_NMVOC,                 &
             PCO_in_Strat,  dtChem,            State_Chm,                    &
             State_Met                                                      )
!
! !USES:
!
    USE PhysConstants, ONLY : AVO, AIRMW
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp),       INTENT(IN) :: ConcClMnd       ! Cl conc [molec/cm3]
    REAL(fp),       INTENT(IN) :: ConcOHmnd       ! OH conc [molec/cm3]
    REAL(fp),       INTENT(IN) :: LCH4_by_OH      ! L(CH4) by OH [1/s]
    REAL(fp),       INTENT(IN) :: LCO_in_Strat    ! Strat L(CO) [molec/cm3/s]
    REAL(fp),       INTENT(IN) :: OHdiurnalFac    ! OH diurnal scale factor [1]
    LOGICAL,        INTENT(IN) :: PCO_fr_CH4_use  ! Use P(CO) fr CH4? [T/F]
    REAL(fp),       INTENT(IN) :: PCO_fr_CH4      ! P(CO) fr CH4 [molec/cm3/s]
    LOGICAL,        INTENT(IN) :: PCO_fr_NMVOC_use! Use P(CO) from NMVOC [T/F]
    REAL(fp),       INTENT(IN) :: PCO_fr_NMVOC    ! P(CO) fr NMVOC [molec/cm3/s]
    REAL(fp),       INTENT(IN) :: PCO_in_Strat    ! Strat P(CO) [molec/cm3/s]
    REAL(fp),       INTENT(IN) :: dtChem          ! Chemistry timestep [s]
    TYPE(ChmState), INTENT(IN) :: State_Chm       ! Chemistry State object
    TYPE(MetState), INTENT(IN) :: State_Met       ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    !========================================================================
    ! carbon_ComputeRateConstants begins here!
    !========================================================================

    ! Initialize
    k_Strat = 0.0_dp
    k_Trop  = 0.0_dp
    trop    = 0.0_dp

    IF ( State_Met%InTroposphere(I,J,L) ) THEN

       !=====================================================================
       ! Grid box (I,J.L) is in the troposphere
       !=====================================================================

       ! Set toggle that we are in the tropopshere
       ! (used for CH4 + FixedCl = Dummy reaction)
       trop = 1.0_dp

       !%%% TAGGED SPECIES HANDLING -- comment out for now
       !%%% State_Chm%SpcData(16)%Info%MW_g * 1.0e-3_fp ! kg/mol

       !---------------------------------------------------------------------
       ! k_Trop(1): Rate [1/s] for rxn: CH4 + OH_E = CO + CO_CH4
       !---------------------------------------------------------------------

       ! OH and Cl concentrations [molec/cm3]
       C(ind_FixedOH) = ConcOHmnd * OHdiurnalFac
       C(ind_FixedCl) = ConcClMnd

       ! CH4 + offline OH reaction rate [1/s]
       ! This is a pseudo-2nd order rate appropriate for CH4 + OH
       IF ( PCO_fr_CH4_use ) THEN
          k_Trop(1) = PCO_fr_CH4 * OHdiurnalFac
          k_Trop(1) = SafeDiv( k_Trop(1), C(ind_CH4)*C(ind_FixedOH), 0.0_dp )
       ELSE
          k_Trop(1) = 2.45E-12_dp * EXP( -1775.0E+0_dp /TEMP )  ! JPL 1997
       ENDIF

       !---------------------------------------------------------------------
       ! k_Trop(2): Rate [1/s        ] for CO + FixedOH = CO2 + CO2fromOH
       !---------------------------------------------------------------------
       k_Trop(2) = GC_OHCO()

       !---------------------------------------------------------------------
       ! k_Trop(3): Rate [molec/cm3/s] for FixedNMVOC   = CO  + COfromNMVOC
       !---------------------------------------------------------------------
       IF ( PCO_fr_NMVOC_use ) THEN
          k_Trop(3) = PCO_fr_NMVOC * OHdiurnalFac
       ENDIF

    ELSE

       !=====================================================================
       ! Grid box (I,J.L) is in the stratosphere
       !=====================================================================

       !---------------------------------------------------------------------
       ! k_Strat(1): Loss rate [1/s] for CH4 -> DUMMY
       ! k_Strat(3): Loss rate [1/s] for CO  -> CO2 + CO2fromOH
       !---------------------------------------------------------------------
       k_Strat(1) = LCH4_by_OH
       k_Strat(3) = LCO_in_Strat

       !---------------------------------------------------------------------
       ! k_Strat(2): Loss rate [molec/cm3/s] for DummyCH4 -> CO
       !
       ! NOTE: This reaction rate is in molec/cm3/s instead of 1/s because
       ! of the way the reaction is written.  CH4_E is a "dummy" species
       ! that is set to 1, so the result of the integration (using the
       ! forward Euler scheme) is
       !
       ! CO = k_Strat(2)  * CH4_E * DT
       !    = molec/cm3/s * 1     * s
       !    = molec/cm3
       !
       ! NOTE: PCO_in_Strat is in [v/v/s]; convert to [molec/cm3] below
       !---------------------------------------------------------------------
       k_Strat(2) = PCO_in_Strat                                             &
                  * AVO                                                      &
                  / ( AIRMW * 1.0e-3_dp )                                    &
                  * State_Met%AirDen(I,J,L)                                  &
                  * 1.0e-6_dp
    ENDIF

  END SUBROUTINE carbon_ComputeRateConstants
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: carbon_ConvertMolecCm3ToKg
!
! !DESCRIPTION: Converts species concentrations (after chemistry) from kg to
!  molec/cm3 and stores into the State_Chm%Species concentration array.
!\\
! !INTERFACE:
!
  SUBROUTINE carbon_ConvertMolecCm3ToKg(                                &
             I,         J,          L,        id_CH4,     id_CO,             &
             id_COch4,  id_COnmvoc, id_CO2,   xnumol_CH4, xnumol_CO2,        &
             xnumol_CO, State_Chm,  State_Met                               )
!
! !USES:
!
    USE Species_Mod,   ONLY : SpcConc
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L      ! Grid box indices
    INTEGER,        INTENT(IN)    :: id_CH4       ! Species index for CH4
    INTEGER,        INTENT(IN)    :: id_CO        ! Species index for CO
    INTEGER,        INTENT(IN)    :: id_COch4     ! Species index for COch4
    INTEGER,        INTENT(IN)    :: id_COnmvoc   ! Species index for COnmvoc
    INTEGER,        INTENT(IN)    :: id_CO2       ! Species index for CO2
    REAL(fp),       INTENT(IN)    :: xnumol_CH4   ! kg CH4 / molec CH4
    REAL(fp),       INTENT(IN)    :: xnumol_CO    ! kg CO  / molec CO
    REAL(fp),       INTENT(IN)    :: xnumol_CO2   ! kg CO2 / molec CO2
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meterorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm    ! Chemistry State object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL(fp)               :: airvol_cm3, convfac_CO

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    !========================================================================
    ! carbon__ConvertMolecCm3ToKg begins here!
    !========================================================================

    ! Point to species array
    Spc => State_Chm%Species

    ! Grid box volume [cm3]
    airvol_cm3  = State_Met%AIRVOL(I,J,L) * 1.0e+6_fp

    ! Send species from molec/cm3 to kg
    IF ( id_CH4 > 0 ) THEN
       Spc(id_CH4)%Conc(I,J,L) = C(ind_CH4) * airvol_cm3 / xnumol_CH4
    ENDIF

    IF ( id_CO > 0 ) THEN
       convfac_CO             = airvol_cm3 / xnumol_CO
       Spc(id_CO)%Conc(I,J,L) = C(ind_CO)  * convfac_CO

       IF ( id_COch4 > 0 ) THEN
          Spc(id_COch4)%Conc(I,J,L) = C(ind_COfromCH4) * convfac_CO
       ENDIF

       IF ( id_COnmvoc > 0 ) THEN
          Spc(id_COnmvoc)%Conc(I,J,L) = C(ind_COfromNMVOC) * convfac_CO
       ENDIF

    ENDIF

    IF ( id_CO2 > 0 ) THEN
       Spc(id_CO2)%Conc(I,J,L) = C(ind_CO2) * airvol_cm3 / xnumol_CO2
    ENDIF

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE carbon_ConvertMolecCm3ToKg
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: carbon_get_co_ch4_flux
!
! !DESCRIPTION: Returns the flux of CO_CH4 in molec/cm3/s for diagnostics.
!\\
!\\
! !INTERFACE:
!
  FUNCTION carbon_Get_COfromCH4_Flux( dtChem ) RESULT ( flux )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: dtChem       ! Chemistry timestep [s]
!
! !RETURN VALUE:
!
    REAL(dp)             :: flux         ! CO_CH4 flux [molec/cm3/s]
!EOP
!------------------------------------------------------------------------------
!BOC

    flux = C(ind_COfromCH4) / dtChem     ! molec/cm3 --> molec/cm3/s

  END FUNCTION carbon_Get_COfromCH4_Flux
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: carbon_get_cofromnmvoc_flux
!
! !DESCRIPTION: Returns the flux of CO_NMVOC in molec/cm3/s for diagnostics.
!\\
!\\
! !INTERFACE:
!
  FUNCTION carbon_Get_COfromNMVOC_Flux( dtChem ) RESULT ( flux )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: dtChem       ! Chemistry timestep [s]
!
! !RETURN VALUE:
!
    REAL(dp)             :: flux         ! CO_NMVOC flux [molec/cm3/s]

!EOP
!------------------------------------------------------------------------------
!BOC

    flux = C(ind_COfromNMVOC) / dtChem   ! molec/cm3 --> molec/cm3/s

  END FUNCTION carbon_Get_COfromNMVOC_Flux
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: carbon_get_co2fromoh_flux
!
! !DESCRIPTION: Returns the flux of CO2_OH in molec/cm3/s for diagnostics.
!\\
!\\
! !INTERFACE:
!
  FUNCTION carbon_Get_CO2fromOH_Flux( dtChem ) RESULT ( flux )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: dtChem       ! Chemistry timestep [s]
!
! !RETURN VALUE:
!
    REAL(dp)             :: flux         ! CO2_OH flux [molec/cm3/s]

!EOP
!------------------------------------------------------------------------------
!BOC

    flux = C(ind_CO2fromOH) / dtChem     ! molec/cm3 --> molec/cm3/s

  END FUNCTION carbon_Get_CO2fromOH_Flux
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: carbon_get_oh_e_flux
!
! !DESCRIPTION: Returns the flux of OH_E in molec/cm3/s for diagnostics.
!\\
!\\
! !INTERFACE:
!
  FUNCTION carbon_Get_FixedOH_Flux( dtChem ) RESULT ( flux )
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: dtChem       ! Chemistry timestep [s]
!
! !RETURN VALUE:
!
    REAL(dp)             :: flux         ! CO2_OH flux [kg/m2/s]

!EOP
!------------------------------------------------------------------------------
!BOC

    flux = C(ind_FixedOH) / dtChem       ! molec/cm3 --> molec/cm3/s

  END FUNCTION carbon_Get_FixedOH_Flux


  !==========================================================================
  ! Rate-law functions
  !==========================================================================

  FUNCTION GC_OHCO() RESULT( k )
    !
    ! Reaction rate for:
    !    OH + CO = HO2 + CO2 (cf. JPL 15-10)
    !
    ! For this reaction, these Arrhenius law terms evaluate to 1:
    !    (300/T)**b0 * EXP(c0/T)
    ! because b0 = c0 = 0.  Therefore we can skip computing these
    ! terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    !
    REAL(dp) :: klo1,   klo2,   khi1,  khi2
    REAL(dp) :: xyrat1, xyrat2, blog1, blog2,   fexp1
    REAL(dp) :: fexp2,  kco1,   kco2,  TEMP300, k
    !
    klo1   = 5.9E-33_dp * K300_OVER_TEMP
    khi1   = 1.1E-12_dp * K300_OVER_TEMP**(-1.3_dp)
    xyrat1 = klo1 * NUMDEN / khi1
    blog1  = LOG10( xyrat1 )
    fexp1  = 1.0_dp / ( 1.0_dp + blog1*blog1 )
    kco1   = klo1 * NUMDEN * 0.6_dp**fexp1 / ( 1.0_dp + xyrat1 )
    klo2   = 1.5E-13_dp
    khi2   = 2.1E+09_dp * K300_OVER_TEMP**(-6.1_dp)
    xyrat2 = klo2 * NUMDEN / khi2
    blog2  = LOG10( xyrat2 )
    fexp2  = 1.0_dp / ( 1.0_dp + blog2*blog2 )
    kco2   = klo2 * 0.6_dp**fexp2 / ( 1.0_dp + xyrat2 )
    k      = kco1 + kco2

  END FUNCTION GC_OHCO

END MODULE carbon_Funcs
