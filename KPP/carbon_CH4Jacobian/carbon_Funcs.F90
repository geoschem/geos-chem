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

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: carbon_ConvertKgToMolecCm3
  PUBLIC :: carbon_ComputeRateConstants
  PUBLIC :: carbon_ConvertMolecCm3ToKg
  PUBLIC :: carbon_InitCarbonKPPFuncs
  PUBLIC :: carbon_CleanupCarbonKPPFuncs
  PUBLIC :: carbon_Get_COfromCH4_Flux
  PUBLIC :: carbon_Get_COfromNMVOC_Flux
  PUBLIC :: GC_OHCO
!
! !PUBLIC TYPES:
!
  ! Species ID flags
  INTEGER  :: id_CH4
  INTEGER  :: id_CO
  INTEGER  :: id_CO2

  ! Jacobian CH4 tracers
  INTEGER, POINTER :: JacobianIDs(:)
  INTEGER          :: numJacobianTracers

  ! Kg species / molec species
  REAL(fp) :: xnumol_CH4
  REAL(fp) :: xnumol_CO
  REAL(fp) :: xnumol_CO2

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
  SUBROUTINE carbon_ConvertKgToMolecCm3( I, J, L, State_Chm, State_Met )
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
    TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meterorology State object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: N
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

    ! Comment out for now - CO2 chemistry is done explicitly in
    ! carbon_gases_mod.F90
    !IF ( id_CO2 > 0 ) THEN
    !   C(ind_CO2) = Spc(id_CO2)%Conc(I,J,L) * xnumol_CO2 / airvol_cm3
    !ENDIF

    ! Do the same for Jacobian CH4 tracers, if any
    IF ( numJacobianTracers > 0 ) THEN
       DO N = 1, numJacobianTracers
          C(JacobianIDs(N)) = Spc(JacobianIDs(N))%Conc(I,J,L) * xnumol_CH4 / airvol_cm3
       ENDDO
    ENDIF

    ! Initialize placeholder species to 1 molec/cm3
    C(ind_DummyCH4trop)  = 1.0_dp
    C(ind_DummyCH4strat) = 1.0_dp
    C(ind_DummyNMVOC)    = 1.0_dp

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
  SUBROUTINE carbon_ComputeRateConstants(                                    &
             I,             J,                 L,                            &
             ConcClMnd,     ConcOHMnd,         LCH4_in_Strat,                &
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
    REAL(fp),       INTENT(IN) :: LCH4_in_Strat   ! Strat L(CH4) [1/s]
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

       ! Set toggle that we are in the troposphere
       ! (used for CH4 + FixedCl = Dummy reaction)
       trop = 1.0_dp

       !%%% TAGGED SPECIES HANDLING -- comment out for now
       !%%% State_Chm%SpcData(16)%Info%MW_g * 1.0e-3_fp ! kg/mol

       ! OH and Cl concentrations [molec/cm3]
       C(ind_FixedOH) = ConcOHmnd
       C(ind_FixedCl) = ConcClMnd

       !---------------------------------------------------------------------
       ! k_Trop(1): Rate [1/s] for rxn for DummyCH4trop -> CO + COfromCH4
       !---------------------------------------------------------------------
       IF ( PCO_fr_CH4_use ) THEN
          k_Trop(1) = PCO_fr_CH4     * OHdiurnalFac
       ENDIF

       !---------------------------------------------------------------------
       ! k_Trop(2): Rate [1/s] for CO + FixedOH = LCObyOH
       !---------------------------------------------------------------------
       k_Trop(2) = GC_OHCO()

       !---------------------------------------------------------------------
       ! k_Trop(3): Rate [molec/cm3/s] for FixedNMVOC   = CO  + PCOfromNMVOC
       !---------------------------------------------------------------------
       IF ( PCO_fr_NMVOC_use ) THEN
          k_Trop(3) = PCO_fr_NMVOC * OHdiurnalFac
       ENDIF

    ELSE

       !=====================================================================
       ! Grid box (I,J.L) is in the stratosphere
       !=====================================================================

       !---------------------------------------------------------------------
       ! k_Strat(1): Loss rate [1/s] for CH4 -> LCH4inStrat
       ! k_Strat(3): Loss rate [1/s] for CO  -> LCOinStrat
       !---------------------------------------------------------------------
       k_Strat(1) = LCH4_in_Strat
       k_Strat(3) = LCO_in_Strat

       !---------------------------------------------------------------------
       ! k_Strat(2): Loss rate [molec/cm3/s] for DummyCH4strat -> CO
       !
       ! NOTE: This reaction rate is in molec/cm3/s instead of 1/s because
       ! of the way the reaction is written.  DummyCH4strat is a "dummy"
       ! species that is set to 1, so the result of the integration (using
       ! the forward Euler scheme) is
       !
       ! CO = k_Strat(2)  * dummyCH4strat * DT
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
  SUBROUTINE carbon_ConvertMolecCm3ToKg( I, J, L, State_Chm, State_Met )
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
    INTEGER                :: N
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
    ENDIF

    ! Comment out for now - CO2 chemistry is done explicitly in
    ! carbon_gases_mod.F90
    !IF ( id_CO2 > 0 ) THEN
    !   Spc(id_CO2)%Conc(I,J,L) = C(ind_CO2) * airvol_cm3 / xnumol_CO2
    !ENDIF

    ! Do the same for Jacobian CH4 tracers, if any
    IF ( numJacobianTracers > 0 ) THEN
       DO N = 1, numJacobianTracers
          Spc(JacobianIDs(N))%Conc(I,J,L) = C(JacobianIDs(N)) * airvol_cm3 / xnumol_CH4
       ENDDO
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

    flux = C(ind_PCOfromCH4) / dtChem     ! molec/cm3 --> molec/cm3/s

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

    flux = C(ind_PCOfromNMVOC) / dtChem   ! molec/cm3 --> molec/cm3/s

  END FUNCTION carbon_Get_COfromNMVOC_Flux
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GC_OHCO
!
! !DESCRIPTION: Returns the rate of tropospheric loss of CO due to chemical reaction with OH
!\\
!\\
! !INTERFACE:
!
  FUNCTION GC_OHCO() RESULT( k )

!
! !REMARKS:
! Tropospheric loss of CO due to chemical rxn w/ OH
!
!  DECAY RATE
!  The decay rate (KRATE) is calculated by:
!
!     OH + CO -> products (JPL 15-10)
!     k = (1 + 0.6Patm) * 1.5E-13
!
!  KRATE has units of [ molec^2 CO / cm6 / s ]^-1,
!  since this is a 2-body reaction.
!
! From JPL 2006: "The  reaction between HO and CO to yield
! H + CO2 akes place on a potential energy surface that
! contains the radical HOCO.  The yield of H and CO2 is
! diminished as the pressure rises.  The loss of reactants
! is thus the sum of two processes, an association to yield
! HOCO and the chemical activation process yielding H and
! CO2." So we now need two complicated reactions.
!
! GY( A0 = 5.9e-33, B0 = 1.,     A1 = 1.1e-12, B1 = -1.3e0,
!     A2 = 1.5e-13, B2 = 0.,     A3 = 2.1e09,  B3 = -6.1e0 )
!
!EOP
!------------------------------------------------------------------------------
!BOC
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
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: carbon_InitCarbonKPPFuncs
!
! !DESCRIPTION: Stores species indices and kg/molecule in module variables
!  for fast lookup.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE carbon_InitCarbonKPPFuncs( kgmolec_CH4, kgmolec_CO, kgmolec_CO2, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE State_Chm_Mod,  ONLY : Ind_
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN)  :: kgmolec_CH4   ! kg CH4 / molec CH4
    REAL(fp), INTENT(IN)  :: kgmolec_CO    ! kg CO  / molec CO
    REAL(fp), INTENT(IN)  :: kgmolec_CO2   ! kg CO2 / molec CO2
!
! !OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(OUT) :: RC            ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, spcName
    CHARACTER(LEN=4)   :: N_char

    !=================================================================
    ! carbon_InitCarbonKPPFuncs begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
       ' -> at carbon_InitCarbonKPPFuncs (in KPP/carbon_CH4Jacobian/carbon_InitCarbonKPPFuncs.F90'
    numJacobianTracers = num_Jtracers

    ! Allocate Jacobian tracer mapping array and assign ids
    IF ( numJacobianTracers > 0 ) THEN
       ALLOCATE( JacobianIDs(numJacobianTracers), STAT=RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error allocating Jacobian tracer ID mapping array!'
          CALL GC_Error( errMsg, RC, thisLoc )
       ENDIF
       DO N = 1, numJacobianTracers
          write( N_char, '(I4.4)' ) N
          spcName = 'CH4_jac' // TRIM(N_char) 
          JacobianIDs(N) = Ind_(TRIM(spcName))
       ENDDO
    ENDIF

    ! Define flags for species ID's
    id_CH4      = Ind_('CH4')
    id_CO       = Ind_('CO2')
    id_CO2      = Ind_('CO' )

    ! Set kg species / molec species
    xnumol_CH4  = kgmolec_CH4
    xnumol_CO   = kgmolec_CO
    xnumol_CO2  = kgmolec_CO2

  END SUBROUTINE carbon_InitCarbonKPPFuncs
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: carbon_CleanupCarbonKPPFuncs
!
! !DESCRIPTION: Deallocates module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE carbon_CleanupCarbonKPPFuncs( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(OUT) :: RC            ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    !=================================================================
    ! carbon_CleanupCarbonKPPFuncs begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS

    IF ( ALLOCATED( JacobianIDs ) ) THEN
       DEALLOCATE( JacobianIDs, STAT=RC )
       CALL GC_CharVar( 'carbon_Funcs.F90:JacobianIDs', 1, RC )
    ENDIF

  END SUBROUTINE carbon_CleanupCarbonKPPFuncs
!EOC
END MODULE carbon_Funcs
