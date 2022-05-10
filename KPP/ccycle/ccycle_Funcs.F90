!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: ccycle_Funcs
!
! !DESCRIPTION: Module containing routines for passing data from GEOS-Chem
!  to the ccycle chemical mechanism solver code.
!\\
!\\
! !INTERFACE:
!
MODULE ccycle_Funcs
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
! !IROUTINE: ccycle_ConvertKgToMolecCm3
!
! !DESCRIPTION: Converts species from kg to molec/cm3 and stores into
!  the "C" concentration array used by the KPP-generated solver code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ccycle_ConvertKgToMolecCm3( I,          J,          L,          &
                                         id_CH4,     id_CO2,     id_CO,      &
                                         xnumol_CH4, xnumol_CO2, xnumol_CO,  &
                                         State_Chm,  State_Met              )
!
! !USES:
!
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
    TYPE(MetState), INTENT(IN) :: State_Met   ! Meterorology State object
    TYPE(ChmState), INTENT(IN) :: State_Chm   ! Chemistry State object
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    REAL(fp)          :: airvol_cm3

    ! Pointers
    REAL(fp), POINTER :: Spc(:,:,:,:)

    !========================================================================
    ! ccycle_ConvertKgToMolecCm3 begins here!
    !========================================================================

    ! Point to species array
    Spc => State_Chm%Species

    ! Grid box volume [cm3]
    airvol_cm3 = State_Met%AIRVOL(I,J,L) * 1.0e+6_fp

    ! Convert species kg to molec/cm3 and store in the C array
    ! (skip if species are not present)
    IF ( id_CH4 > 0 ) THEN
       C(ind_CH4) = Spc(I,J,L,id_CH4) * xnumol_CH4 / airvol_cm3
    ENDIF

    IF ( id_CO > 0 ) THEN
       C(ind_CO)  = Spc(I,J,L,id_CO)  * xnumol_CO  / airvol_cm3
    ENDIF

    IF ( id_CO2 > 0 ) THEN
       C(ind_CO2) = Spc(I,J,L,id_CO2) * xnumol_CO2 / airvol_cm3
    ENDIF

    ! Initialize externally-read species (denoted by _E) to 1 molec/cm3
    C(ind_CH4_E)    = 1.0_dp
    C(ind_NMVOC_E)  = 1.0_dp

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE ccycle_ConvertKgToMolecCm3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ccycle_ComputeRateConstants
!
! !DESCRIPTION: Computes the rate constants used in the ccycle chemical
!  mechanism.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ccycle_ComputeRateConstants(                                    &
             I,        J,        L,           bAirDens,    bCl,              &
             bOH,      CH4loss,  GMI_Prod_CO, GMI_Loss_CO, PCO_nmVOC,        &
             PCO_CH4,  LPCO_CH4, dtChem,      tCosZ,       State_Chm,        &
             State_Met                                                      )
!
! !USES:
!
    USE PhysConstants, ONLY : AVO
    USE State_Chm_Mod, ONLY : ChmState
    USE State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN) :: I, J, L
    REAL(fp),       INTENT(IN) :: bAirDens
    REAL(fp),       INTENT(IN) :: bCl
    REAL(fp),       INTENT(IN) :: bOH
    REAL(fp),       INTENT(IN) :: CH4loss
    REAL(fp),       INTENT(IN) :: GMI_Prod_CO
    REAL(fp),       INTENT(IN) :: GMI_Loss_CO
    REAL(fp),       INTENT(IN) :: PCO_NMVOC
    LOGICAL,        INTENT(IN) :: LPCO_CH4
    REAL(fp),       INTENT(IN) :: PCO_CH4
    REAL(fp),       INTENT(IN) :: dtChem
    REAL(fp),       INTENT(IN) :: tCosZ
    TYPE(ChmState), INTENT(IN) :: State_Chm     ! Chemistry State object
    TYPE(MetState), INTENT(IN) :: State_Met     ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL(fp) :: kgs_to_mcm3sCO
    REAL(dp) :: mw_CO_kg_mol
    REAL(fp) :: fac_Diurnal

    !========================================================================
    ! ccycle_ComputeRateConstants begins here!
    !========================================================================

    ! Initialize
    facDiurnal     = 0.0_dp
    kgs_to_mcm3sCO = 0.0_dp
    k_Strat        = 0.0_dp
    k_Trop         = 0.0_dp
    trop           = 0.0_dp
    mw_CO_kg_mol   = MW(ind_CO) * 1.0e-3_dp

    IF ( State_Met%InTroposphere(I,J,L) ) THEN

       !---------------------------------------------------------------------
       ! Grid box (I,J.L) is in the troposphere
       !---------------------------------------------------------------------

       ! Strat rate
       K_STRAT = 0.d0
       TROP    = 1.0_dp  ! Toggle

       ! Scaling factor for diurnal cycles - zero at night
       IF ( State_Met%SUNCOSmid(I,J) > 0.0_fp .and. tCosZ > 0.0_fp ) THEN
          facDiurnal = ( State_Met%SUNCOSmid(I,J) / tCosZ  )                &
                     * ( 86400.0_fp               / dtChem )
       ENDIF

       ! Trop Rates
       ! State_Chm%SpcData(16)%Info%MW_g * 1.0e-3_fp ! kg/mol

       ! Conversion factor from kg/s CO to molec/cm3/s CO
       kgs_to_mcm3sCO = ( AVO                                                &
                      / ( mw_CO_kg_mol * State_Met%AIRVOL(I,J,L) ) )         &
                      * 1.0e-6_fp

       ! P(CO) from CH4
       k_Trop(1) = PCO_CH4     * facDiurnal
       IF ( .not. LPCO_CH4 ) THEN
          k_Trop(1) = 2.45E-12_dp * EXP( -1775.E0_dp /TEMP )  ! JPL 1997
       ENDIF
       K_TROP(2)    = GC_OHCO()
       K_TROP(3)    = PCO_NMVOC * facDiurnal

       ! OH concentration, as read from disk [molec/cm3]
       C(ind_OH_E) = bOH * State_Met%AIRNUMDEN(I,J,L) * facDiurnal

       ! Cl continue oncentration, as read from disk [molec/cm3]
       C(ind_Cl_E) = bCl * State_Met%AIRNUMDEN(I,J,L) * 1e-9_fp

       ! CH4 + offline OH reaction rate [1/s]
       IF ( LPCO_CH4 ) THEN
          k_Trop(1) = safediv( k_trop(1), ( C(ind_CH4) *C(ind_OH_E) ), 0.0_dp )
       ENDIF

    ELSE

       !---------------------------------------------------------------------
       ! Grid box (I,J.L) is in the stratosphere
       !---------------------------------------------------------------------

       ! Strat Rates [1/s]
       K_STRAT(1) = CH4LOSS                                  ! [1/s]
       K_STRAT(2) = GMI_PROD_CO * State_Met%AIRNUMDEN(I,J,L) ! [molec/cm3/s]
       K_STRAT(3) = GMI_LOSS_CO                              ! [1/s]

    ENDIF

  END SUBROUTINE ccycle_ComputeRateConstants
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ccycle_ConvertMolecCm3ToKg
!
! !DESCRIPTION: Converts species concentrations (after chemistry) from kg to
!  molec/cm3 and stores into the State_Chm%Species concentration array.
!\\
! !INTERFACE:
!
  SUBROUTINE ccycle_ConvertMolecCm3ToKg( I,          J,         L,           &
                                         id_CH4,     id_CO,     id_COch4,    &
                                         id_COnmvoc, id_CO2,    xnumol_CH4,  &
                                         xnumol_CO2, xnumol_CO, State_Chm,   &
                                         State_Met                          )
!
! !USES:
!
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
    REAL(fp)          :: airvol_cm3

    ! Pointers
    REAL(fp), POINTER :: Spc(:,:,:,:)

    ! Point to species array
    Spc => State_Chm%Species

    ! Grid box volume [cm3]
    airvol_cm3 = State_Met%AIRVOL(I,J,L) * 1.0e+6_fp

    ! Send species from molec/cm3 to kg
    IF ( id_CH4 > 0 ) THEN
       Spc(I,J,L,id_CH4)     = C(ind_CH4)      * airvol_cm3 / xnumol_CH4
    ENDIF

    IF ( id_CO > 0 ) THEN
       Spc(I,J,L,id_CO)      = C(ind_CO)       * airvol_cm3 / xnumol_CO
    ENDIF

    IF ( id_COch4 > 0 ) THEN
       Spc(I,J,L,id_COch4)   = C(ind_CO_CH4)   * airvol_cm3 / xnumol_CO
    ENDIF

    IF ( id_COnmvoc > 0 ) THEN
       Spc(I,J,L,id_COnmvoc) = C(ind_CO_NMVOC) * airvol_cm3 / xnumol_CO
    ENDIF

    IF ( id_CO2 > 0 ) THEN
       Spc(I,J,L,id_CO2)     = C(ind_CO2)      * airvol_cm3 / xnumol_CO2
    ENDIF

    Spc => NULL()

  END SUBROUTINE ccycle_ConvertMolecCm3ToKg

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

END MODULE ccycle_Funcs
