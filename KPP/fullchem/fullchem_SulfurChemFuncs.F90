!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: fullchem_SulfurChemFuncs
!
! !DESCRIPTION: FlexChem module for multiphase sulfate chemistry, via KPP.
!\\
!\\
! !INTERFACE:

MODULE fullchem_SulfurChemFuncs
!
! !USES:
!
  USE PhysConstants
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: fullchem_ConvertAlkToEquiv
  PUBLIC :: fullchem_ConvertEquivToAlk
  PUBLIC :: fullchem_HetDropChem
  PUBLIC :: fullchem_InitSulfurChem
  PUBLIC :: fullchem_SulfurAqChem
  PUBLIC :: fullchem_SulfurCldChem
!
! !PUBLIC TYPES:
!
  ! Species ID flags
  INTEGER                :: id_ACTA,   id_CH2O,   id_DMS,    id_DST1
  INTEGER                :: id_DST2,   id_DST3,   id_DST4,   id_H2O2
  INTEGER                :: id_HCL,    id_HCOOH,  id_HMS,    id_HNO3
  INTEGER                :: id_MSA,    id_NH3,    id_NH4,    id_NIT
  INTEGER                :: id_NITs,   id_O3,     id_OH,     id_pFe
  INTEGER                :: id_SALA,   id_SALAAL, id_SALACL, id_SALC
  INTEGER                :: id_SALCAL, id_SALCCL, id_SO2,    id_SO4
  INTEGER                :: id_SO4s
!
! !DEFINED_PARAMETERS
!
  REAL(fp),  PARAMETER   :: TCVV_S    = AIRMW / 32e+0_fp ! hard-coded MW
  REAL(fp),  PARAMETER   :: TCVV_N    = AIRMW / 14e+0_fp ! hard-coded MW
  REAL(fp),  PARAMETER   :: SMALLNUM  = 1e-20_fp
  REAL(fp),  PARAMETER   :: CM3PERM3  = 1.e6_fp

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fullchem_ConvertAlkToEquiv
!
! !DESCRIPTION: Converts sea salt alkalinity to equivalents.  Abstracted
!  out from fullchem_mod.F90 to prevent compilation conflicts for other
!  KPP chemical mechanisms
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE fullchem_ConvertAlkToEquiv()
!
! !USES:
!
    USE gckpp_Global,     ONLY : C,          MW
    USE gckpp_Parameters, ONLY : ind_SALAAL, ind_SALCAL
!EOP
!------------------------------------------------------------------------------
!BOC
    C(ind_SALAAL) = C(ind_SALAAL) * ( MW(ind_SALAAL) * 7.0e-5_fp )
    C(ind_SALCAL) = C(ind_SALCAL) * ( MW(ind_SALCAL) * 7.0e-5_fp )
  END SUBROUTINE fullchem_ConvertAlkToEquiv
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fullchem_ConvertEquivToAlk
!
! !DESCRIPTION: Converts sea salt alkalinity to equivalents.  Abstracted
!  out from fullchem_mod.F90 to prevent compilation conflicts for other
!  KPP chemical mechanisms
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE fullchem_ConvertEquivToAlk()
!
! !USES:
!
    USE gckpp_Global,     ONLY : C,          MW
    USE gckpp_Parameters, ONLY : ind_SALAAL, ind_SALCAL
!EOP
!------------------------------------------------------------------------------
!BOC
    C(ind_SALAAL) = C(ind_SALAAL) / ( MW(ind_SALAAL) * 7.0e-5_fp )
    C(ind_SALCAL) = C(ind_SALCAL) / ( MW(ind_SALCAL) * 7.0e-5_fp )
  END SUBROUTINE fullchem_ConvertEquivToAlk
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fullchem_SulfurAqchem
!
! !DESCRIPTION: Main aqueous/aerosol chemistry driver routine.  Sets up the
!  vector of aqueous chemistry rates for the KPP chemistry solver.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE fullchem_SulfurAqChem( I,         J,          L,                &
                                    Input_Opt, State_Chm,  State_Grid,       &
                                    State_Met, RC                           )
!
! !USES:
!

    USE ErrCode_Mod
    USE gckpp_Global
    USE gckpp_Parameters
    USE Input_Opt_Mod,    ONLY : OptInput
    USE rateLawUtilFuncs
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Met_Mod,    ONLY : MetState
    USE State_Grid_Mod,   ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L    ! Lon, lat, level indices
    TYPE(MetState), INTENT(IN)    :: State_Met  ! Meteorology State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(OptInput), INTENT(IN)    :: Input_Opt  ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(IN)    :: State_Chm  ! Chemistry State object
!
! OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC         ! Success or failure
!
! !REMARKS:
!
!  ! Reaction List (by K_MT() index)
!    1) SO2 + O3 + 2SALAAL --> SO4mm + O2 : From Sulfate_mod - 24 Mar 2021
!
! !REVISION HISTORY:
!  24 Mar 2021 - M. Long - Initial Version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

    ! Scalars
    LOGICAL            :: SALAAL_gt_0_1
    LOGICAL            :: SALCAL_gt_0_1
    REAL(fp)           :: k_ex

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !======================================================================
    ! fullchem_SulfurAqChem begins here!
    !======================================================================

    ! Initialize
    RC            = GC_SUCCESS
    k_ex          = 0.0_dp
    K_MT          = 0.0_dp
    SALAAL_gt_0_1 = ( C(ind_SALAAL) > 0.1_dp )
    SALCAL_gt_0_1 = ( C(ind_SALCAL) > 0.1_dp )

    !======================================================================
    ! Reaction rates [1/s] for fine sea salt alkalinity (aka SALAAL)
    !
    ! K_MT(1) : SALAAL + SO2 + O3 = SO4 - SALAAL
    ! K_MT(2) : SALAAL + HCl      = SALACL
    ! K_MT(3) : SALAAL + HNO3     = NIT
    !
    ! NOTE: SALAAL_gt_0_1 prevents div-by-zero errors
    !======================================================================

    !------------------------------------------------------------------------
    ! SALAAL + SO2 + O3 = SO4 - SALAAL
    !------------------------------------------------------------------------
    IF ( SALAAL_gt_0_1 ) THEN

       ! 1st order uptake
       k_ex = Ars_L1K( area   = State_Chm%WetAeroArea(I,J,L,11),             &
                       radius = State_Chm%AeroRadi(I,J,L,11),                &
                       gamma  = 0.11_dp,                                     &
                       srMw   = SR_MW(ind_SO2)                              )

       ! Assume SO2 is limiting, so recompute rxn rate accordingly
       K_MT(1) = kIIR1Ltd( C(ind_SO2), C(ind_SALAAL), k_ex ) / C(ind_O3)
    ENDIF

    !------------------------------------------------------------------------
    ! SALAAL + HCL = SALACL
    !------------------------------------------------------------------------
    IF ( SALAAL_gt_0_1 ) THEN

       ! 1st order uptake
       k_ex = Ars_L1K( area   = State_Chm%WetAeroArea(I,J,L,11),             &
                       radius = State_Chm%AeroRadi(I,J,L,11),                &
                       gamma  = 0.07_dp,                                     &
                       srMw   = SR_MW(ind_HCl)                              )

       ! Assume HCl is limiting, so recompute reaction rate accordingly
       K_MT(2) = kIIR1Ltd( C(ind_HCl), C(ind_SALAAL), k_ex )
    ENDIF

    !------------------------------------------------------------------------
    ! SALAAL + HNO3 = NIT
    !------------------------------------------------------------------------
    IF ( SALAAL_gt_0_1 ) THEN

       ! 1st order uptake
       k_ex = Ars_L1K( area   = State_Chm%WetAeroArea(I,J,L,11),             &
                       radius = State_Chm%AeroRadi(I,J,L,11),                &
                       gamma  = 0.5_dp,                                      &
                       srMw   = SR_MW(ind_HNO3)                             )

       ! Assume HNO3 is limiting, so recompute reaction rate accordingly
       K_MT(3) = kIIR1Ltd( C(ind_HNO3), C(ind_SALAAL), k_ex )
    ENDIF

    !========================================================================
    ! Reaction rates [1/s] for coarse sea salt alkalinity (aka SALAAL)
    !
    ! K_MT(4) : SALCAL + SO2 + O3 = SO4s - SALCAL
    ! K_MT(5) : SALCAL + HCl      = SALCCL
    ! K_MT(6) : SALCAL + HNO3     = NITs
    !
    ! NOTE: SALCAL_gt_0_1 prevents div-by-zero errors
    !========================================================================

    !------------------------------------------------------------------------
    ! SALCAL + SO2 + O3 = SO4s - SALCAL
    !------------------------------------------------------------------------
    IF ( SALCAL_gt_0_1 ) THEN

       ! 1st order uptake
       k_ex = Ars_L1K( area   = State_Chm%WetAeroArea(I,J,L,12),             &
                       radius = State_Chm%AeroRadi(I,J,L,12),                &
                       gamma  = 0.11_dp,                                     &
                       srMw   = SR_MW(ind_SO2)                              )

       ! Assume SO2 is limiting, so recompute rxn rate accordingly
       K_MT(4) = kIIR1Ltd( C(ind_SO2), C(ind_SALCAL), k_ex ) / C(ind_O3)
    ENDIF

    !------------------------------------------------------------------------
    ! SALCAL + HCl = SALCCL
    !------------------------------------------------------------------------
    IF ( SALCAL_gt_0_1 ) THEN

       ! 1st order uptake
       k_ex = Ars_L1K( area   = State_Chm%WetAeroArea(I,J,L,12),             &
                       radius = State_Chm%AeroRadi(I,J,L,12),                &
                       gamma  = 0.07_dp,                                     &
                       srMw   = SR_MW(ind_HCl)                              )

       ! Assume HCl is limiting, so recompute rxn rate accordingly
       K_MT(5) = kIIR1Ltd( C(ind_HCl), C(ind_SALCAL), k_ex )
    ENDIF

    !------------------------------------------------------------------------
    ! SALCAL + HNO3 = NITs
    !------------------------------------------------------------------------
    IF ( SALCAL_gt_0_1 ) THEN

       ! 1st order uptake
       k_ex = Ars_L1K( area   = State_Chm%WetAeroArea(I,J,L,12),             &
                       radius = State_Chm%AeroRadi(I,J,L,12),                &
                       gamma  = 0.5_dp,                                      &
                       srMw   = SR_MW(ind_HNO3)                             )

       ! Assume HNO3 is limiting, so recompute rxn rate accordingly
       K_MT(6) = kIIR1Ltd( C(ind_HNO3), C(ind_SALCAL), k_ex )
    ENDIF

  END SUBROUTINE fullchem_SulfurAqChem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fullchem_SulfurCldChem
!
! !DESCRIPTION: Routine that compute reaction rates for sulfur chemistry
!  in cloud, so that these can be passed to the KPP chemical solver.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE fullchem_SulfurCldChem( I,         J,         L,                &
                                    Input_Opt,  State_Chm, State_Diag,       &
                                    State_Grid, State_Met, size_res,         &
                                    RC                                      )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    INTEGER,        INTENT(IN)    :: I, J, L
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,        INTENT(OUT)   :: size_res    ! Should we call HetDropChem?
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL            :: prtDebug
    INTEGER            :: N
    CHARACTER(LEN=63)  :: OrigUnit

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !========================================================================
    ! fullchem_SulfurCldChem begins here!
    !========================================================================

    ! Initialize
    RC       = GC_SUCCESS
    size_res = .FALSE.
    ErrMsg   = ''
    ThisLoc  = &
  ' -> at fullchem_SulfurCldChem (in KPP/fullchem/fullchem_SulfurChemFuncs.F90'

    ! Should we print debug output?
    prtDebug             = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    !------------------------------------------------------------------------
    ! SO2 chemistry
    !------------------------------------------------------------------------
    CALL Set_SO2( I          = I,                                            &
                  J          = J,                                            &
                  L          = L,                                            &
                  Input_Opt  = Input_Opt,                                    &
                  State_Chm  = State_Chm,                                    &
                  State_Diag = State_Diag,                                   &
                  State_Grid = State_Grid,                                   &
                  State_Met  = State_Met,                                    &
                  size_res   = size_res,                                     &
                  RC         = RC                                           )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "SET_SO2"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE fullchem_SulfurCldChem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fullchem_HetDropChem
!
! !DESCRIPTION: Subroutine HET\_DROP\_CHEM estimates the in-cloud sulfate
!  production rate in heterogeneous cloud droplets based on the Yuen et al.,
!  1996 parameterization. (bec, 6/16/11)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE fullchem_HetDropChem( I,         J,         L,         SR,      &
                                   Input_Opt, State_Met, State_Chm          )
!
! !USES:
!
    USE fullchem_RateLawFuncs, ONLY : HOBrUptkByHSO3m, HOBrUptkBySO3mm
    USE fullchem_RateLawFuncs, ONLY : HOClUptkByHSO3m, HOClUptkBySO3mm
    USE gckpp_Global
    USE gckpp_Parameters
    USE gckpp_Precision
    USE Input_Opt_Mod,         ONLY : OptInput
    USE PhysConstants,         ONLY : AIRMW, AVO, PI, g0
    USE rateLawUtilFuncs
    USE State_Chm_Mod,         ONLY : ChmState, IND_
    USE State_Met_Mod,         ONLY : MetState
    USE Time_Mod,              ONLY : Get_Ts_Chem
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT)   :: SR          ! Sulfate production rate
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Dry sea-salt density [kg/m3]
    REAL(dp), PARAMETER :: SS_DEN        = 2200.0_dp

    ! sigma of the size distribution for sea-salt (Jaegle et al., 2011)
    REAL(fp), PARAMETER :: SIG_S         = 1.8e+0_dp

    ! geometric dry mean diameters [m] for computing lognormal size distribution
    REAL(dp), PARAMETER :: RG_S          = 0.4e-6_dp !(Jaegle et a., 2011)
    REAL(dp), PARAMETER :: RG_D2         = 1.5e-6_dp !(Ginoux et al., 2001)
    REAL(dp), PARAMETER :: RG_D3         = 2.5e-6_dp
    REAL(dp), PARAMETER :: RG_D4         = 4.e-6_dp

    ! To prevent multiple divisions
    REAL(dp), PARAMETER   :: THREE_FOURTHS = 3.0_dp / 4.0_dp
    REAL(dp), PARAMETER   :: NINE_HALVES   = 9.0_dp / 2.0_dp

!
! !LOCAL VARIABLES:
!
    REAL(dp)              :: alpha_NH3,  alpha_SO2, alpha_H2O2
    REAL(dp)              :: alpha_HNO3, alpha_B,   alpha_CN
    REAL(dp)              :: alpha_W,    alpha_SO4, sum_gas
    REAL(dp)              :: H,          NDss,      CN
    REAL(dp)              :: W,          K,         arg
    REAL(dp)              :: DTCHEM,     APV,       DSVI
    REAL(dp)              :: B,          NH3,       SO2
    REAL(dp)              :: H2O2,       HNO3,      SO4
    REAL(dp)              :: CNss,       MW_SO4,    MW_SALC
    REAL(dp)              :: CVF,        R1,        R2
    REAL(dp)              :: XX,         FC,        LST
    REAL(dp)              :: XX1,        XX2,       XX3
    REAL(dp)              :: XX4,        XX5,       GNH3

    ! Pointers
    REAL(fp), POINTER     :: AD(:,:,:)
    REAL(fp), POINTER     :: AIRDEN(:,:,:)
    REAL(fp), POINTER     :: AIRVOL(:,:,:)
    REAL(fp), POINTER     :: OMEGA(:,:,:)
    REAL(fp), POINTER     :: U(:,:,:)
    REAL(fp), POINTER     :: V(:,:,:)

    !=================================================================
    ! HET_DROP_CHEM begins here!
    !=================================================================

    ! Initialize pointers
    AD     => State_Met%AD
    AIRDEN => State_Met%AIRDEN
    AIRVOL => State_Met%AIRVOL
    OMEGA  => State_Met%OMEGA
    U      => State_Met%U
    V      => State_Met%V

    ! Zero output argument for safety's sake
    SR     =  0.0_dp

    ! Zero/initialize local variables for safety's sake
    arg    =  0.0_dp
    B      =  0.0_dp
    CN     =  0.0_dp
    CVF    =  1.0e3_fp * AIRMW / ( AIRDEN(I,J,L) * AVO ) ! molec/cm3 -> v/v
    DSVI   =  0.0_dp
    DTCHEM =  GET_TS_CHEM()                              ! seconds
    GNH3   =  0.0_dp
    K      =  0.0_dp
    LST    =  0.0_dp
    R1     =  0.0_dp
    R2     =  0.0_dp
    W      =  0.0_dp
    XX     =  0.0_dp
    XX1    =  0.0_dp
    XX2    =  0.0_dp
    XX3    =  0.0_dp
    XX4    =  0.0_dp
    XX5    =  0.0_dp

    ! FC is guaranteed to be > 1e-4, because HET_DROP_CHEM
    ! is not called otherwise (bmy, 07 Oct 2021)
    FC     =  State_Met%CLDF(I,J,L)

    !! <<>> SET THE INPUT UNITS! EITHER CONVERT IN THE ROUTINE OR
    !! <<>> CONVERT BEFOREHAND. BUT EVERYTHING IS CURRENTLY mcl/cm3
    !! <<>> AND HET_DROP_CHEM EXPECTS V/V

    ! XX* are calculated below to be consistent with
    ! Sulfate_Mod(). Values are different when
    ! computed with KPP-based variables. HET_DROP_CHEM()
    ! could use some attention to make is consistent with
    ! KPP.
    !
    ! NOTE: Use function SafeExp, which will prevent the exponential from
    ! blowing up.  Also if the entire expression will evaluate to zero
    ! then skip the exponential, which is more computationally efficient.
    !    -- Bob Yantosca, 14 Oct 2021
    !

    ! SO2 + H2O2
    R1  = C(ind_SO2)  * CVF
    R2  = C(ind_H2O2) * CVF
    K   = K_CLD(1)    / CVF/ FC
    Arg = ( R1 - R2 ) * ( K * DTCHEM )
    IF ( IsSafeExp( Arg ) .and. ABS( Arg ) > 0.0_dp ) THEN
       XX  = EXP( Arg )
       XX1 = ( R1 * R2 ) * ( XX - 1.0_dp ) / ( ( R1 * XX ) - R2 )
    ELSE
       XX1 = WhenExpCantBeDone( R1, R2, K, DTCHEM )
    ENDIF

    ! SO2 + O3
    R2  = C(ind_O3) * CVF
    K   = K_CLD(2)  / CVF / FC
    Arg = ( R1 - R2 ) * ( K * DTCHEM )
    IF ( IsSafeExp( Arg ) .and. ABS( Arg ) > 0.0_dp ) THEN
       XX = EXP( Arg )
       XX2 = ( R1 * R2 ) * ( XX - 1.0_dp ) / ( ( R1 * XX ) - R2 )
    ELSE
       XX2 = WhenExpCantBeDone( R1, R2, K, DTCHEM )
    ENDIF

    ! Metal catalyzed oxidation of SO2 pathway
    K   = -K_CLD(3) / FC
    Arg = K * DTCHEM
    XX3 = 0.0_dp
    IF ( IsSafeExp( Arg ) ) THEN
       XX  = EXP( Arg )
       XX3 = R1 * ( 1.0_dp - XX )
    ENDIF

    ! HSO3- + HOCl and SO3-- + HOCl
    R1  = C(ind_SO2) * CVF * State_Chm%HSO3_aq(I,J,L)
    R2  = C(ind_HOCl)  * CVF
    K   = HOClUptkByHSO3m(State_Het) / CVF
    Arg = ( R1 - R2 ) * ( K * DTCHEM )
    IF ( IsSafeExp( Arg ) .and. ABS( Arg ) > 0.0_dp ) THEN
       XX  = EXP( Arg )
       XX4 = ( R1 * R2 )  * ( XX - 1.0_dp ) / ( ( R1 * XX )  - R2 )
    ELSE
       XX4 = WhenExpCantBeDone( R1, R2, K, DTCHEM )
    ENDIF

    ! SO3-- + HOCl (add to HSO3- + HOCl rate)
    R1  = C(ind_SO2) * CVF * State_Chm%SO3_aq(I,J,L)
    K   = HOClUptkBySO3mm(State_Het) / CVF
    Arg = ( R1 - R2 ) * ( K * DTCHEM )
    IF ( IsSafeExp( Arg ) .and. ABS( Arg ) > 0.0_dp ) THEN
       XX  = EXP( Arg )
       XX4 = XX4 + ( ( R1 * R2 ) * ( XX - 1.0_fp ) / ( ( R1 * XX ) - R2 ) )
    ELSE
       XX4 = XX4 + WhenExpCantBeDone( R1, R2, K, DTCHEM )
    ENDIF

    ! HSO3- + HOBr
    R1  = C(ind_SO2) * CVF * State_Chm%HSO3_aq(I,J,L)
    R2  = C(ind_HOBr)  * CVF
    K   = HOBrUptkByHSO3m(State_Het) / CVF
    Arg = ( R1 - R2 ) * ( K * DTCHEM )
    IF ( IsSafeExp( Arg ) .and. ABS( Arg ) > 0.0_dp ) THEN
       XX  = EXP( Arg )
       XX5 = ( R1 * R2 ) * ( XX - 1.0_fp ) / ( ( R1 * XX ) - R2 )
    ELSE
       XX5 = WhenExpCantBeDone( R1, R2, K, DTCHEM )
    ENDIF

    ! SO3-- + HOBr (add to HSO3- + HOBr rate)
    R1  = C(ind_SO2) * CVF * State_Chm%SO3_aq(I,J,L)
    K   = HOBrUptkBySO3mm(State_Het) / CVF
    Arg = ( R1 - R2 ) * ( K * DTCHEM )
    IF ( IsSafeExp( Arg ) .and. ABS( Arg ) > 0.0_dp ) THEN
       XX  = EXP( Arg )
       XX5 = XX5 + ( ( R1 * R2 ) * ( XX - 1.0_dp ) / ( ( R1 * XX ) - R2 ) )
    ELSE
       XX5 = XX5 + WhenExpCantBeDone( R1, R2, K, DTCHEM )
    ENDIF

    ! Sum of all rates
    LST = XX1 + XX2 + XX3 + XX4 + XX5

    !### Debug print
    !IF (I .eq. 12 .and. J .eq. 7 .and. L .eq. 1) THEN
    !   write(*,*) '<<>> XX: ', XX1, XX2, XX3, XX4, XX5
    !ENDIF

    IF ( LST > R1 ) THEN
       XX1 = ( R1 * XX1 ) / LST
       XX2 = ( R1 * XX2 ) / LST
       XX3 = ( R1 * XX3 ) / LST
       XX4 = ( R1 * XX4 ) / LST
       XX5 = ( R1 * XX5 ) / LST
       LST = XX1 + XX2 + XX3 + XX4 + XX5
     ENDIF

    ! Convert gas phase concentrations from [v/v] to [pptv]
    NH3  = State_Chm%Species(id_NH3)%Conc(I,J,L) * CVF * 1.0e+12_dp
    SO2  = MAX( C(ind_SO2) * CVF - ( LST*FC ), 1.0e-20_dp ) * 1.0e+12_dp
    H2O2 = C(ind_H2O2)* CVF * 1.0e12_dp
    HNO3 = C(ind_HNO3)* CVF * 1.0e12_dp

    ! Set molecular weight local variables
    MW_SO4  = State_Chm%SpcData(id_SO4)%Info%MW_g
    MW_SALC = State_Chm%SpcData(id_SALC)%Info%MW_g

    ! Convert sulfate aerosol concentrations from [v/v] to [ug/m3]
    SO4 = ( C(ind_SO4) * CVF * AD(I,J,L) * 1.0e+9_dp ) /                     &
          ( ( AIRMW / MW_SO4 ) * AIRVOL(I,J,L) )

    ! Convert in cloud sulfate production rate from [v/v/timestep] to
    ! [ug/m3/timestep]
    B  = ( LST * AD(I,J,L) * 1.0e+9_dp ) /                                   &
         ( ( AIRMW / MW_SO4 ) * AIRVOL(I,J,L) )

    ! Convert coarse-mode aerosol concentrations from [v/v] to [#/cm3]
    ! based on equation in Hofmann, Science, 1990.
    ! First convert from [v/v] to [kg/m3 air]
    CNss = State_Chm%Species(id_SALC)%Conc(I,J,L)*CVF * AD(I,J,L)         &
         / ( ( AIRMW / MW_SALC ) * AIRVOL(I,J,L) )

    ! Now convert from [kg/m3 air] to [#/cm3 air]
    ! Sea-salt
    ARG  = NINE_HALVES * ( LOG( SIG_S ) )**2
    NDss = ( THREE_FOURTHS * CNss                           )               &
         / ( PI * SS_DEN * RG_S**3 * SafeExp( Arg, 0.0_dp ) )               &
         * 1.e-6_dp

    ! Total coarse mode number concentration [#/cm3]
    CN = NDss ! sea-salt

    ! Determine regression coefficients based on the local SO2 concentration
    IF ( SO2 <= 200.00_fp ) THEN
       alpha_B    =  0.5318_dp
       alpha_NH3  = -1.67e-7_dp
       alpha_SO2  =  2.59e-6_dp
       alpha_H2O2 = -1.77e-7_dp
       alpha_HNO3 = -1.72e-7_dp
       alpha_W    =  1.22e-6_dp
       alpha_CN   =  4.58e-6_dp
       alpha_SO4  = -1.00e-5_dp
    ELSE IF ( SO2 > 200.00_dp .and. SO2 <= 500.0_dp ) THEN
       alpha_B    =  0.5591_dp
       alpha_NH3  =  3.62e-6_dp
       alpha_SO2  =  1.66e-6_dp
       alpha_H2O2 =  1.06e-7_dp
       alpha_HNO3 = -5.45e-7_dp
       alpha_W    = -5.79e-7_dp
       alpha_CN   =  1.63e-5_dp
       alpha_SO4  = -7.40e-6_dp
    ELSE IF ( SO2 > 500.0_dp .and. SO2 < 1000.0_dp ) THEN
       alpha_B    =  1.1547_dp
       alpha_NH3  = -4.28e-8_dp
       alpha_SO2  = -1.23e-7_dp
       alpha_H2O2 = -9.05e-7_dp
       alpha_HNO3 =  1.73e-7_dp
       alpha_W    =  7.22e-6_dp
       alpha_CN   =  2.44e-5_dp
       alpha_SO4  =  3.25e-5_dp
    ELSE                          ! SO2 > 1000
       alpha_B    =  1.1795_dp
       alpha_NH3  =  2.57e-7_dp
       alpha_SO2  = -5.54e-7_dp
       alpha_H2O2 = -1.08e-6_dp
       alpha_HNO3 =  1.95e-6_dp
       alpha_W    =  6.14e-6_dp
       alpha_CN   =  1.64e-5_dp
       alpha_SO4  =  2.48e-6_dp
    ENDIF

    ! Updraft velocity over the oceans [cm/s]
    ! 500 cm/s is too high. Get W from the met field. (qjc, 04/10/16)
    !W = 500e+0_fp
    W = -OMEGA(I,J,L) / ( AIRDEN(I,J,L) * g0 ) * 100e+0_dp

    ! Compute H (integration time interval * air parcel velocity) [m]
    ! DTCHEM is the chemistry timestep in seconds

    ! Compute air parcel velocity [m/s]
    !APV = SQRT( (U(I,J,L) * U(I,J,L)) + (V(I,J,L) * V(I,J,L)) )
    !(qjc, 04/10/16)
    APV = SQRT( U(I,J,L)**2 + V(I,J,L)**2 ) + ( W**2 * 1.0e-4_dp )

    H   = DTCHEM * APV          ![m]

    sum_gas = ( alpha_NH3  * NH3  ) + ( alpha_SO2  * SO2  ) +                &
              ( alpha_H2O2 * H2O2 ) + ( alpha_HNO3 * HNO3 )

    DSVI = ( alpha_B * B ) +                                                 &
           ( ( ( alpha_CN * CN) + ( alpha_W * W ) + ( alpha_SO4 * SO4 ) +    &
                sum_gas ) * H )

    ! Only calculate SR when air parcel rises, in consistence with
    ! Yuen et al. (1996) (qjc, 04/10/16)
    IF ( W > 0.0_dp .and. C(ind_SO2) > 0.0_dp ) THEN

       ! additional sulfate production that can be
       ! attributed to ozone [ug/m3/timestep]
       ! Don't allow SR to be negative
       SR = MAX( ( DSVI - B ), 0.0_dp )

       ! Skip further computation if SR = 0
       IF ( SR > 0.0_dp ) THEN

          ! Convert SR from [ug/m3/timestep] to [v/v/timestep]
          SR = SR * ( AIRMW / MW_SO4 ) * 1.e-9_dp / AIRDEN(I,J,L)

          ! Don't produce more SO4 than SO2 available after AQCHEM_SO2
          ! -- SR is dSO4/timestep (v/v) continue onvert
          !    to 1st order rate
          SR = MIN( SR, SO2 / 1.0e12_dp ) / ( C(ind_SO2) * CVF * DT )
       ENDIF
    ENDIF

    ! Free pointers
    AD     => NULL()
    AIRDEN => NULL()
    AIRVOL => NULL()
    OMEGA  => NULL()
    U      => NULL()
    V      => NULL()

  END SUBROUTINE fullchem_HetDropChem
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: WhenExpCantBeDone
!
! !DESCRIPTION: Prevents floating point errors if exponential terms in routine
!  Het_Drop_Chem above can't be done.  In the case of a negative XX, R should be
!  approximated as R1, instead of R2.  In other words,
!  R1 * R2 * ( XX - 1.D0 ) / ( ( R1 * XX ) - R2 )
!  reaches different limits when XX reaches +Inf and -Inf.
!\\
!\\
! !INTERFACE:
!
  FUNCTION WhenExpCantBeDone( R1, R2, K, DT ) RESULT( R )
!
! !USES:
!
    USE gckpp_Precision, ONLY : dp
!
! !INPUT PARAMETERS:
!
    REAL(dp), INTENT(IN) :: R1   ! 1st term
    REAL(dp), INTENT(IN) :: R2   ! 2nd term
    REAL(dp), INTENT(IN) :: K    ! Rate [1/s]
    REAL(dp), INTENT(IN) :: DT   ! timesetep [s]
!
! !RETURN VALUE:
!
    REAL(dp)             :: R    ! new rate [1/s]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
    REAL(dp) :: DIFF

    DIFF = R1 - R2

    ! R1 <  R2
    IF ( DIFF < 0.0_dp ) THEN
       R = R1
       RETURN
    ENDIF

    ! R1 >  R2
    IF ( DIFF > 0.0_dp ) THEN
       R = R2
       RETURN
    ENDIF

    ! R1 == R2
    R = R1 - 1.0_dp / ( K * DT + ( 1.0_dp / R1 ) )

  END FUNCTION WhenExpCantBeDone
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_so2
!
! !DESCRIPTION: Subroutine SET\_SO2 is the SO2 chemistry subroutine.
!  (rjp, bmy, 11/26/02, 8/26/10) Adapted from CHEM_SO2() in SULFATE_MOD
!  (MSL - Spring 2021)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE SET_SO2( I,         J,          L,          Input_Opt,          &
                      State_Chm, State_Diag, State_Grid, State_Met,          &
                      Size_Res,  RC                                         )
!
! !USES:
!
    USE ErrCode_Mod
    USE gckpp_Global
    USE Input_Opt_Mod,   ONLY : OptInput
    USE rateLawUtilFuncs
    USE Species_Mod,     ONLY : SpcConc
    USE State_Chm_Mod,   ONLY : ChmState
    USE State_Diag_Mod,  ONLY : DgnState
    USE State_Grid_Mod,  ONLY : GrdState
    USE State_Met_Mod,   ONLY : MetState
    USE Time_Mod,        ONLY : Get_Ts_Chem
!
! !INPUT PARAMETERS:
!
    INTEGER,        INTENT(IN)    :: I, J, L     ! Grid box indices
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    LOGICAL,        INTENT(OUT)   :: Size_Res    ! Should we call HetDropChem?
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  Reaction List (by Rokjin Park)
!  ============================================================================
!  (1 ) SO2 production:
!       DMS + OH, DMS + NO3 (saved in CHEM_DMS)
!                                                                             .
!  (2 ) SO2 loss:
!       (a) SO2 + OH  -> SO4
!       (b) SO2       -> drydep
!       (c) SO2 + H2O2 or O3 (aq) -> SO4
!                                                                             .
!  (3 ) SO2 = SO2_0 * exp(-bt) +  PSO2_DMS/bt * [1-exp(-bt)]
!                                                                             .
!       where b is the sum of the reaction rate of SO2 + OH and the dry
!       deposition rate of SO2, PSO2_DMS is SO2 production from DMS in
!       MixingRatio/timestep.
!                                                                             .
!  If there is cloud in the gridbox (fraction = fc), then the aqueous
!  phase chemistry also takes place in cloud. The amount of SO2 oxidized
!  by H2O2 in cloud is limited by the available H2O2; the rest may be
!  oxidized due to additional chemistry, e.g, reaction with O3 or O2
!  (catalyzed by trace metal).
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(fp),  PARAMETER  :: HPLUS_45  = 3.16227766016837953e-5_fp  !pH = 4.5
    REAL(fp),  PARAMETER  :: HPLUS_50  = 1.0e-5_fp  !pH = 5.0
    REAL(fp),  PARAMETER  :: MINDAT    = 1.e-20_fp
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL               :: IS_OFFLINE
    LOGICAL               :: IS_FULLCHEM
    INTEGER               :: BULK
    INTEGER               :: IBIN
    REAL(fp)              :: K0,     Ki,      KK,     M,    L1
    REAL(fp)              :: L2,     L3,      Ld,     F,    Fc
    REAL(fp)              :: RK,     RKT,     DTCHEM, DT_T, TK
    REAL(fp)              :: F1,     RK1,     RK3,    SO20, AVO_over_LWC
    REAL(fp)              :: SO2_cd, H2O20,   L2S,    L3S
    REAL(fp)              :: LWC,    KaqH2O2, KaqO3,  PATM, RHO, CNVFAC
    REAL(fp)              :: ALK,    ALK1,    ALK2,   SO2_AfterSS
    REAL(fp)              :: AlkA,   AlkC
    REAL(fp)              :: Kt1,    Kt2
    REAL(fp)              :: PSO4E,  PSO4F,   Kt1N,    Kt2N
    REAL(fp)              :: XX, Kt1L, Kt2L
    REAL(fp)              :: HPLUS,  SO4nss, TNH3,   TNO3,  GNO3, ANIT
    REAL(fp)              :: LSTOT,  ALKdst, ALKss,  ALKds, NH3, CL, TNA
    REAL(fp)              :: SSCvv,  aSO4,   SO2_sr, SR,    TANIT
    REAL(fp)              :: TFA,  TAA,   TDCA    ! (jmm, 12/03/2018)
    REAL(fp)              :: SO2_gas,   PH2SO4d_tot
    REAL(fp)              :: H2SO4_cd,  H2SO4_gas

    ! (qjc, 04/10/16)
    REAL(fp)              :: L5,L5S,SRo3,SRhobr
    REAL(fp)              :: L5_1,L5S_1,L3_1,L3S_1,KaqO3_1
    REAL(fp)              :: HSO3aq, SO3aq
    REAL(fp)              :: SO2_AfterSS0, rSIV, fupdateHOBr_0
    REAL(fp)              :: HCO3, HCHOBr, KO3, KHOBr, f_srhobr, HOBr0
    REAL(fp)              :: TMP

    REAL(fp)              :: KaqO2, L4, L4S, MnII, FeIII
    REAL(fp)              :: DUST,  Mn_ant,  Mn_nat
    REAL(fp)              :: Mn_tot, Mn_d,    Fe_d
    REAL(fp)              :: Fe_ant, Fe_nat,  Fe_tot
    REAL(fp)              :: Fe_d_ant, Fe_d_nat

    REAL(fp)              :: L6,L6S,SRhocl,L6_1,L6S_1      !XW
    REAL(fp)              :: fupdateHOCl_0  !XW
    REAL(fp)              :: HCHOCl, KHOCl, f_srhocl, HOCl0 !XW

    REAL(fp)              :: KaqHCHO, KaqHMS, KaqHMS2, HMSc ! JMM, MSL

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)
    REAL(fp), POINTER      :: SSAlk(:)

    CHARACTER(LEN=255)     :: ErrMsg, ThisLoc

#ifdef LUO_WETDEP
    ! For Luo et al wetdep scheme
    LOGICAL                :: Is_QQ3D
#endif

    !========================================================================
    ! SET_SO2 begins here!
    !========================================================================
    IF ( id_H2O2 < 0 .or. id_SO2 < 0  ) RETURN

    ! Initialize
    RC       = GC_SUCCESS
    size_res = .FALSE.
    ErrMsg   = ''
    ThisLoc  = &
      ' -> at SET_SO2 (in module KPP/fullchem/fullchem_SulfurChemFuncs.F90)'

    ! Initialize variables
    Spc                         => State_Chm%Species
    SSAlk                       => State_Chm%SSAlk(I,J,L,:)
    State_Chm%isCloud(I,J,L)    = 0.0_fp
    State_Chm%pHCloud(I,J,L)    = 0.0_fp
    State_Chm%QLxpHCloud(I,J,L) = 0.0_fp
    State_Chm%HSO3_aq(I,J,L)    = 1.0e-32_fp
    State_Chm%SO3_aq(I,J,L)     = 1.0e-32_fp
    DTCHEM                      = GET_TS_CHEM()  ! Timestep [s]
    IS_FULLCHEM                 = Input_Opt%ITS_A_FULLCHEM_SIM
    IS_OFFLINE                  = ( .not. IS_FULLCHEM )
    Ld                          = 0.0_fp
    LSTOT                       = 0.0_fp
    RHO                         = State_Met%AIRDEN(I,J,L)
    CNVFAC                      = 1.E3_fp * AIRMW / ( RHO * AVO ) !mcl/cm3->v/v
    SO20                        = Spc(id_SO2)%Conc(I,J,L) * CNVFAC
    SO2_AfterSS                 = Spc(id_SO2)%Conc(I,J,L) * CNVFAC
    H2O20                       = Spc(id_H2O2)%Conc(I,J,L) * CNVFAC
    KaqH2O2                     = 0.0_fp
    KaqO3                       = 0.0_fp
    KaqO3_1                     = 0.0_fp
    KaqO2                       = 0.0_fp
    K_CLD                       = 0.0_fp
    HPLUS                       = 0.0_fp

    ! Factor to convert AIRDEN from [kg air/m3] to [molec air/cm3]
    F                           = 1000.e+0_fp / AIRMW * AVO * 1.e-6_fp

    ! Meteorological data
    PATM = State_Met%PMID_DRY( I, J, L ) / ( ATM * 1.e-2_fp ) ! Press, dry [atm]
    TK   = State_Met%T(I,J,L)                                 ! Temperature [K]
    FC   = State_Met%CLDF(I,J,L)                              ! Cloud frac [1]

    ! Get liquid water content [m3 H2O/m3 air] within cloud from met flds
    ! Units: [kg H2O/kg air] * [kg air/m3 air] * [m3 H2O/1e3 kg H2O]
#ifdef LUO_WETDEP
    ! Luo et al wetdep scheme
    IF ( Is_QQ3D ) THEN
       LWC = State_Met%QL(I,J,L) * State_Met%AIRDEN(I,J,L) * 1e-3_fp + &
            MAX( 0.0_fp, State_Chm%QQ3D(I,J,L) * DTCHEM )
    ELSE
       LWC = State_Met%QL(I,J,L) * State_Met%AIRDEN(I,J,L) * 1e-3_fp
    ENDIF
#else
    ! Default scheme
    LWC = State_Met%QL(I,J,L) * State_Met%AIRDEN(I,J,L) * 1e-3_fp
#endif

    ! QL can sometimes be negative, so force LWC to be positive
    LWC = MAX( 0.0_fp, LWC )

    ! LWC is a grid-box averaged quantity. To improve the representation
    ! of sulfate chemistry, we divide LWC by the cloud fraction and
    ! compute sulfate chemistry based on the LWC within the cloud.  We
    ! get the appropriate grid-box averaged mass of SO2 and sulfate by
    ! multiplying these quantities by FC AFTER computing the aqueous
    ! sulfur chemistry within the cloud. (lzh, jaf, bmy, 5/27/11)
    LWC = SafeDiv( LWC, FC, 0.0_fp )


    ! If (1) there is cloud, (2) there is SO2 present, (3) T > -15 C, and
    ! (4) liquid water content (LWC) is present (but not small enough to
    ! make divisions blow up), then compute sulfate production in cloud.
    IF ( ( FC          > 1.0e-4_fp  )   .and.                                &
         ( SO2_AfterSS > MINDAT     )   .and.                                &
#ifdef LUO_WETDEP
         ( TK          > 237.0_fp   )   .and.                                &
#else
         ( TK          > 258.0_fp   )  .and.                                 &
#endif
         ( LWC         > 1.0e-20_fp ) ) THEN

       !===========================================================
       ! NOTE...Sulfate production from aquatic reactions of SO2
       ! with H2O2 & O3 is computed here and followings are
       ! approximations or method used for analytical (integral)
       ! solution of these computations.
       !
       ! 1) with H2O2(aq)
       !      [HSO3-] + [H+] + [H2O2(aq)] => [SO4=]     (rxn)
       !      d[SO4=]/dt = k[H+][HSO3-][H2O2(aq)] (M/s) (rate)
       !
       ! we can rewrite k[H+][HSO3-] as K1 pSO2 hSO2,
       ! where pSO2 is equilibrium vapor pressure of SO2(g)
       ! in atm, and hSO2 is henry's law constant for SO2
       !
       ! Therefore, rate can be written as
       !
       !       k * K1 * pSO2 * hSO2 * pH2O2 * hH2O2,
       !
       ! where pH2O2 is equilibrium vapor pressure of H2O2(g),
       ! and hH2O2 is henry's law constant for H2O2. Detailed
       ! values are given in AQSET_SO2 routine.
       !
       ! Let us define a fraction of gas phase of A species
       ! in equilibrium with aqueous phase as
       !
       !        xA  = 1/(1+f),
       !
       ! where  f   = hA * R * T * LWC,
       !        hA  = Henry's constant,
       !        R   = gas constant,
       !        T   = temperature in kelvin,
       !        LWC = liquid water content [m3/m3]
       !
       ! As a result, the rate would become:
       !
       !    d[SO4=]
       !    ------- = k K1 hSO2 hH2O2 xSO2 xH2O2 P P [SO2][H2O2]
       !      dt
       !      ^       ^                            ^   ^    ^
       !      |       |____________________________|   |    |
       !
       !   mole/l/s               mole/l/s            v/v  v/v
       !
       !
       ! And we multiply rate by (LWC * R * T / P) in order to
       ! convert unit from mole/l/s to v/v/s
       !
       ! Finally we come to
       !
       !    d[SO4=]
       !    ------- = KaqH2O2 [SO2][H2O2],
       !      dt
       !
       ! where
       !
       !   KaqH2O2 = k K1 hSO2 hH2O2 xSO2 xH2O2 P LWC R T,
       !
       ! this new rate corresponds to a typical second order
       ! reaction of which analytical (integral) solution is
       !
       !   X  = A0 B0 ( exp[(A0-B0) Ka t] - 1 )
       !      / ( A0 exp[(A0-B0) Ka t] - B0 )
       !
       ! inserting variables into solution then we get
       ! [SO4=] =  [SO2][H2O2](exp[([SO2]-[H2O2]) KaqH2O2 t] - 1 )
       !        / ( [SO2] exp[([SO2]-[H2O2]) KaqH2O2 t] - [H2O2] )
       !
       ! Note...Exactly same method can be applied to O3 reaction
       ! in aqueous phase with different rate constants.
       !===========================================================

       ! Get concentrations for cloud pH calculation (bec, 12/23/11)

       ! <<>><<>><<>><<>>
       ! HAVE TO DO SOME UNIT CONVERSION HERE. THIS ROUTINE IS CALLED
       ! WITHIN FLEXCHEM_MOD WHERE SPECIES ARE IN MOLEC/CM3
       ! <<>><<>><<>><<>>

       ! Get sulfate concentration and convert from [v/v] to
       ! [moles/liter]
       ! Use a cloud scavenging ratio of 0.7

       SO4nss = 1.e+3 * ( Spc(id_SO4)%Conc(I,J,L) * 0.7e+0_fp &
              + Spc(id_SO4s)%Conc(I,J,L) ) / ( LWC * AVO ) ! mcl/cm3 -> mol/L

       ! Get HMS cloud concentration and convert from [v/v] to
       ! [moles/liter] (jmm, 06/13/2018)
       ! Use a cloud scavenging ratio of 0.7
       ! assume nonvolatile like sulfate for realistic cloud pH
       HMSc = 0.0_fp
       IF ( IS_FULLCHEM .and. id_HMS > 0 ) THEN
          HMSc = 1.e+3 * Spc(id_HMS)%Conc(I,J,L) * 0.7_fp / ( LWC * AVO ) ! mcl/cm3 -> mol/L
       ENDIF

       ! Get total ammonia (NH3 + NH4+) concentration [v/v]
       ! Use a cloud scavenging ratio of 0.7 for NH4+
       TNH3 = ( ( Spc(id_NH4)%Conc(I,J,L) * 0.7e+0_fp ) &
              + Spc(id_NH3)%Conc(I,J,L) ) * CNVFAC

       ! Get total chloride (SALACL + HCL) concentration [v/v]
       ! Use a cloud scavenging ratio of 0.7
       CL = ( Spc(id_SALACL)%Conc(I,J,L) * 0.7e+0_fp ) &
            + Spc(id_SALCCL)%Conc(I,J,L)
       CL = ( CL + Spc(id_HCL)%Conc(I,J,L) ) * CNVFAC

       ! Get total formic acid concentration [v/v]
       ! jmm (12/3/18)
       ! no cloud scavenging because gases?
       TFA = Spc(id_HCOOH)%Conc(I,J,L) * CNVFAC

       ! Get total acetic acid concentration [v/v]
       ! jmm (12/3/18)
       ! no cloud scavenging b/c gases?
       TAA = Spc(id_ACTA)%Conc(I,J,L) * CNVFAC

       ! Get total sea salt NVC concentration expressed as NA+ equivalents
       ! and convert from [MND] to [moles/liter]
       ! NVC is calculated to balance initial Cl- + alkalinity in
       ! seas salt. Note that we should not consider SO4ss here.
       ! Use a cloud scavenging ratio of 0.7 for fine aerosols
       TNA      = 1.e3_fp*( Spc(id_SALA)%Conc(I,J,L)*0.7e+0_fp &
            + Spc(id_SALC)%Conc(I,J,L) ) * &
            ( 31.6e+0_fp * 0.359e+0_fp / 23.e+0_fp ) / &
            ( LWC * AVO ) ! mcl/cm3 -> mol/L

       ! Get total dust cation concentration [mol/L]
       ! Use a cloud scavenging ratio of 1 for dust
       ! to be consistent for how it was calculated for
       ! metal catalyzed SO2 oxidation
       ! Use asumption of dust being 3% soluble Ca2+ and
       ! 0.6% soluble Mg2+ by mass (Fairlie et al., 2010)
       !
       ! Dust treated at non-volatile cation and charge applied in
       ! pH calculation
       !
       ! Move dust calculation from SO2 Metal catalzyed oxidation
       ! up here becasue needed for cloud pH
       ! jmm (12/3/18)
       !
       ! Get dust concentrations [MND -> ng/m3]

       DUST = ( Spc(id_DST1)%Conc(I,J,L)*0.7_fp + Spc(id_DST2)%Conc(I,J,L) + &
            Spc(id_DST3)%Conc(I,J,L) + Spc(id_DST4)%Conc(I,J,L) )            &
            * 1.e+15_fp * State_Chm%SpcData(id_DST1)%Info%MW_g / AVO

       ! Conversion from dust mass to Ca2+ and Mg2+ mol:
       !     0.071*(1/40.08)+0.011*(1/24.31) = 2.22e-3
       !     (Engelbrecht et al., 2016)
       !     1e-12_fp from m3->L & ng->g
       TDCA     = DUST * 2.22e-15_fp / LWC

       ! Get total nitrate (HNO3 + NIT) concentrations [v/v]
       ! Use a cloud scavenging ratio of 0.7 for NIT
       TNO3 = ( Spc(id_HNO3)%Conc(I,J,L) +             &
              ( Spc(id_NIT)%Conc(I,J,L)  * 0.7e+0_fp ) + &
              Spc(id_NITs)%Conc(I,J,L) ) * CNVFAC
       GNO3 = Spc(id_HNO3)%Conc(I,J,L) * CNVFAC ! For Fahey & Pandis decision algorithm

       ! Calculate cloud pH
       CALL GET_HPLUS( SO4nss, HMSc, TNH3, TNO3,  SO2_AfterSS,   CL, TNA, TDCA, &
                       TFA,    TAA,  TK,   PATM,  LWC, HPLUS_45, HPLUS  )

       ! Store the cloud pH quantities
       State_Chm%isCloud(I,J,L)    =  1.0_fp
       State_Chm%pHCloud(I,J,L)    = -1.0_fp * log10(HPLUS)
       State_Chm%QLxpHCloud(I,J,L) = State_Chm%pHCloud(I,J,L)             &
            * State_Met%QL(I,J,L)


       FeIII = 0.0_fp
       MnII  = 0.0_fp
       IF ( Input_Opt%LMETALCATSO2 ) THEN

          !--------------------------------------------------------
          ! Metal catalyzed oxidation of SO2 pathway
          !--------------------------------------------------------

          ! Get dust concentrations [v/v -> ng/m3]
#ifdef TOMAS
          ! TOMAS uses its own dust tracers and does not
          ! carry DST1-4.  Set DUST to zero here. (mps, 2/2/18)
          DUST = 0e+0_fp
#endif
          ! Calculate Fe and Mn natural [ng m-3]
          ! Assume that Fe is 3.5% of total dust mass based on
          ! Taylor and McLennan [1985]
          Fe_nat = DUST * 35e-3_fp
          ! and Mn is 50 times less than Fe based on Desbouefs et al.[2005]
          Mn_nat = Fe_nat / 50e+0_fp

          ! Anthropogenic Fe concentrations [mcl/cm3 -> ng/m3]
          IF ( id_pFe > 0 ) THEN
                Fe_ant = Spc(id_pFe)%Conc(I,J,L) * CNVFAC * &
                         1.e+12_fp * State_Met%AD(I,J,L) &
                         / ( AIRMW / State_Chm%SpcData(id_pFe)%Info%MW_g ) &
                         / State_Met%AIRVOL(I,J,L)
!             Fe_ant = Spc(id_pFe)%Conc(I,J,L) * 1.e+15_fp * &
!                  State_Chm%SpcData(id_DST1)%Info%MW_g / AVO
          ELSE
             Fe_ant = 0e+0_fp
          ENDIF

          ! Calculate Mn anthropogenic [ng m-3]
          ! assume anthropogenic Mn is 1/30 times anthropogenic Fe
          Mn_ant = Fe_ant / 10e+0_fp

          ! Calculate total Mn and Fe [ng m-3]
          Mn_tot = Mn_ant + Mn_nat
          Fe_tot = Fe_ant + Fe_nat

          ! Convert Mn and Fe [ng m-3] to [mole l-1]

          ! Assume that 50% of Mn is dissolved [Spokes et al., 1994]
          ! Hardcoded MW for Mn
          IF ( LWC > 0e+0_fp ) THEN
             ! Units: ng/m3 * (g/ng) / (g/mol) / (m3 H2O / m3 air) * (m3/L)
             Mn_d = Mn_tot * 1e-9_fp / 54.94e+0_fp / LWC * 1e-3_fp
             Mn_d = Mn_d * 0.5e+0_fp
          ELSE
             Mn_d = 0e+0_fp
          ENDIF

          ! Solubility of Fe is 10% for anthropogenic, and 1% for dust
          IF ( LWC > 0e+0_fp ) THEN
             Fe_d_ant = Fe_ant * 1e-9_fp / &
                  State_Chm%SpcData(id_pFe)%Info%MW_g / &
                  LWC * 1e-3_fp
             Fe_d_nat = Fe_nat * 1e-9_fp / &
                  State_Chm%SpcData(id_pFe)%Info%MW_g / &
                  LWC * 1e-3_fp
             Fe_d     = Fe_d_ant * 0.1e+0_fp + &
                  Fe_d_nat * 0.01e+0_fp
          ELSE
             Fe_d     = 0e+0_fp
          ENDIF

          ! Impose a dependence of Fe speciation on sunlight
          IF ( State_Met%SUNCOS(I,J) > 0e+0_fp ) THEN
             ! Assume 10% of dissolved Fe is in Fe(III)
             !oxidation state during the daytime
             FeIII = Fe_d * 0.1e+0_fp
          ELSE
             ! Assume 90% of dissolved Fe is in Fe(III)
             ! oxidation state during the nighttime
             FeIII = Fe_d * 0.9e+0_fp
          ENDIF

          ! Assume that dissolved Mn is in Mn(II) oxidation state all of
          ! the time
          MnII = Mn_d
       ENDIF

       ! Compute aqueous rxn rates for SO2
       CALL AQCHEM_SO2( I       = I,                                         &
                        J       = J,                                         &
                        L       = L,                                         &
                        LWC     = LWC,                                       &
                        T       = TK,                                        &
                        P       = PATM,                                      &
                        SO2     = SO2_AfterSS,                               &
                        H2O2    = Spc(id_H2O2)%Conc(I,J,L) * CNVFAC,         &
                        O3      = Spc(id_O3)%Conc(I,J,L)   * CNVFAC,         &
                        HCHO    = Spc(id_CH2O)%Conc(I,J,L) * CNVFAC,         &
                        Hplus   = Hplus,                                     &
                        MnII    = MnII,                                      &
                        FeIII   = FeIII,                                     &
                        KaqH2O2 = KaqH2O2,                                   &
                        KaqO3   = KaqO3,                                     &
                        KaqO3_1 = KaqO3_1,                                   &
                        KaqO2   = KaqO2,                                     &
                        HSO3aq  = HSO3aq,                                    &
                        SO3aq   = SO3aq,                                     &
                        KaqHCHO = KaqHCHO,                                   &
                        KaqHMS  = KaqHMS,                                    &
                        KaqHMS2 = KaqHMS2                                   )


       K_CLD(1) = KaqH2O2 * FC * CNVFAC   ! v/v/s --> cm3/mcl/s
       K_CLD(2) = KaqO3   * FC * CNVFAC   ! v/v/s --> cm3/mcl/s
       K_CLD(3) = KaqO2   * FC           ! 1/s
       ! vvvvvv Hold off using CloudHet2R until after initial S-chem benchmark
       !        -- MSL
       !K_CLD(1) = CloudHet2R( Spc(id_SO2)%Conc(I,J,L), &
       !                       Spc(id_H2O2)%Conc(I,J,L), FC, KaqH2O2 * CNVFAC )
       !K_CLD(2) = CloudHet2R( Spc(id_SO2)%Conc(I,J,L), &
       !                       Spc(id_O3)%Conc(I,J,L),   FC, KaqO3   * CNVFAC )
       !K_CLD(3) computed below

       ! HMS reaction rates (skip if HMS isn't defined)
       IF ( IS_FULLCHEM .and. id_HMS > 0 ) THEN
          K_CLD(4) = KaqHCHO * FC * CNVFAC
          K_CLD(5) = KaqHMS  * FC
          K_CLD(6) = KaqHMS2 * FC
          ! Leave comments here (bmy, 18 Jan 2022)
          !          CloudHet2R( Spc(id_HMS)%Conc(I,J,L), &
          !                      Spc(id_CH2O)%Conc(I,J,L), FC, KaqHCHO*CNVFAC )
          !          CloudHet1R( FC, KaqHMS ) ! KaqHMS is pseudo-1st order
          !          CloudHet2R( Spc(id_HMS)%Conc(I,J,L), &
          !                      Spc(id_OH)%Conc(I,J,L), FC, & ...)       ENDIF
       ENDIF

#ifdef TOMAS
       !%%%%%%%%%%%%%%%%% BUG FIX FOR TOMAS %%%%%%%%%%%%%%%%%%%%%%%
       ! NOTE: TOMAS uses its own dust tracers and does not
       ! carry ALKdst.  Set ALKdst to zero here. (bmy, 1/28/14)
       ALKdst = 0e+0_fp
#else

       ! For other simulations, Sum up the contributions from
       ! DST1 thru DST4 tracers into ALKdst. (bmy, 1/28/14)
       ! mcl/cm3 -> ug/m3
       ALKdst = ( Spc(id_DST1)%Conc(I,J,L) + Spc(id_DST2)%Conc(I,J,L) +      &
            Spc(id_DST3)%Conc(I,J,L) + Spc(id_DST4)%Conc(I,J,L) ) * CNVFAC * &
            1.e+9_fp * State_Met%AD(I,J,L)                       &
            / ( AIRMW / State_Chm%SpcData(id_DST1)%Info%MW_g ) &
            / State_Met%AIRVOL(I,J,L)
#endif

       ! mcl/cm3 -> ug/m3
       ALKss  = ( Spc(id_SALA)%Conc(I,J,L) &
            + Spc(id_SALC)%Conc(I,J,L) ) * CNVFAC * &
            1.e+9_fp * State_Met%AD(I,J,L)                       &
            / ( AIRMW / State_Chm%SpcData(id_SALA)%Info%MW_g ) &
            / State_Met%AIRVOL(I,J,L)

       ALKds = ALKdst + ALKss

       ! Get NH3 concentrations (v/v)
       NH3 = Spc(id_NH3)%Conc(I,J,L)*CNVFAC

       ! Initialize
       size_res= .FALSE.

       ! Fahey and Seinfeld decision algorithm
       ! NOTE: This is ugly, needs refactoring.  For now, just added
       ! whitespace to improve readability (bmy, 01 Oct 2021)
       IF ( H2O20 > SO2_afterss + 1e-9_fp ) THEN
          size_res = .FALSE.

       ELSE IF( LWC < 0.1e-6_fp ) THEN !10^-6 coversion from g/m3 --> m3/m3
          size_res = .TRUE.

       ELSE IF( gno3 > NH3 ) THEN

          IF ( So2_afterss >= 5.e-9_fp                                 .and. &
               H2O20       >= SO20              ) size_res = .FALSE.

          IF ( LWC         >= 0.3e-6_fp                                .and. &
               So2_afterss >= 3.e-9_fp                                 .and. &
               H2O20       >= So2_afterss       ) size_res = .FALSE.

          IF ( ALKds       >= 5.e+0_fp                                 .and. &
               LWC         >= 0.5e-6_fp                                .and. &
               H2O20       >= So2_afterss       ) size_res = .FALSE.

          IF ( LWC         >= 0.1e-6_fp                                .and. &
               gno3        <= (NH3 + 2.e-9_fp)  ) size_res = .FALSE.

       ELSE IF( LWC >= 0.5e-6_fp ) THEN

          IF ( H2O20       >=                                                &
               ( 0.9_fp * So2_afterss )         ) size_res = .FALSE.

          IF ( NH3         <= 1.e-9_fp                                 .and. &
               ALKds       >= 5.e+0_fp                                 .and. &
               So2_afterss <= 10.e-9_fp         ) size_res = .FALSE.

       ELSE IF( LWC >= 0.3e-6_fp ) THEN

          IF ( NH3         >= (gno3 + 5.e-9_fp)                        .and. &
               So2_afterss <= 10.e-9_fp         ) size_res = .FALSE.

          IF ( gno3        <= 1.e-9_fp                                 .and. &
               NH3         >= (gno3 + 2.e-9_fp) ) size_res = .FALSE.

          IF ( gno3        <= 7.e-9_fp                                 .and. &
               NH3         >= (gno3 + 3.e-9_fp) ) size_res = .FALSE.

          IF ( ALKds       >= 3.e+0_fp                                 .and. &
               NH3         <= 10e-9_fp                                 .and. &
               So2_afterss <= 5e-9_fp           ) size_res = .FALSE.

          IF ( ALKds       >= 5.e+0_fp                                 .and. &
               NH3         <= 10.e-9_fp                                .and. &
               So2_afterss <= 5.e-9_fp          ) size_res = .FALSE.

          IF ( So2_afterss >= 1.5e-9_fp                                .and. &
               H2O20       >= So2_afterss       ) size_res = .FALSE.

          IF ( NH3         <= 12.e-9_fp                                .and. &
               ALKds       >= 10.e+0_fp         ) size_res = .FALSE.

          IF ( NH3         <= 1.e-9_fp                                 .and. &
               ALKds       >= 4.e+0_fp                                 .and. &
               So2_afterss <= 10.e-9_fp         ) size_res = .FALSE.

          IF ( NH3         <= 5.e-9_fp                                 .and. &
               ALKds       >= 6.e+0_fp                                 .and. &
               So2_afterss <= 10.e-9_fp         ) size_res = .FALSE.

          IF ( NH3         <= 7.e-9_fp                                 .and. &
               ALKds       > -8.e+0_fp                                 .and. &
               So2_afterss <= 10.e-9_fp         ) size_res = .FALSE.

       ELSE IF( LWC >= 0.1e-6_fp ) THEN

          IF ( NH3         <= 1.e-9_fp                                 .and. &
               ALKds       >= 5.e+0_fp          ) size_res = .FALSE.

          IF ( NH3         <= 5.e-9_fp                                 .and. &
               ALKds       >= 10.e+0_fp         ) size_res = .FALSE.

          IF ( gno3        <= 1.e-9_fp                                 .and. &
               NH3         >= (gno3 + 2.e-9_fp)                        .and. &
               So2_afterss <= 7.e-9_fp          ) size_res = .FALSE.

          IF ( gno3        <= 1.e-9_fp                                 .and. &
               NH3         >= (gno3 + 2.e-9_fp)                        .and. &
               ALKds       >= 2.e+0_fp          ) size_res = .FALSE.

          IF ( gno3        <= 3.e-9_fp                                 .and. &
               NH3         >= (gno3 + 4.e-9_fp) ) size_res = .FALSE.

          IF ( gno3        <= 7.e-9_fp                                 .and. &
               NH3         >= (gno3 + 3.e-9_fp)                        .and. &
               So2_afterss <= 5.e-9_fp          ) size_res = .FALSE.

          IF ( gno3        <= 7.e-9_fp                                 .and. &
               NH3         >= (gno3 + 3.e-9_fp)                        .and. &
               ALKds       >= 4.e+0_fp                                 .and. &
               So2_afterss <= 9.e-9_fp          ) size_res = .FALSE.

          IF ( ALKds       >= 3.e+0_fp                                 .and. &
               NH3         <= 3.e-9_fp                                 .and. &
               So2_afterss <= 4.e-9_fp          ) size_res = .FALSE.

          IF ( ALKds       >= 5.e+0_fp                                 .and. &
               So2_afterss <= 5.e-9_fp                                 .and. &
               NH3         <= 7.e-9_fp          ) size_res = .FALSE.

          IF ( NH3         >= (gno3 + 2.e-9_fp)                        .and. &
               So2_afterss <= 5.e-9_fp          ) size_res = .FALSE.

          IF ( NH3         >= (gno3 + 4.e-9_fp)                        .and. &
               So2_afterss <= 10.e-9_fp         ) size_res = .FALSE.

          IF ( ALKds       >= 2.e+0_fp                                 .and. &
               NH3         <= 10.e-9_fp                                .and. &
               H2O20       >= So2_afterss       ) size_res = .FALSE.

          IF ( NH3         <= 1.e-9_fp                                 .and. &
               So2_afterss >= 3.e-9_fp                                 .and. &
               H2O20       >= So2_afterss       ) size_res = .FALSE.

       ELSE

          size_res = .TRUE.

       ENDIF


       State_Chm%HSO3_AQ(I,J,L) = HSO3aq
       State_Chm%SO3_AQ(I,J,L)  = SO3aq

    ENDIF
!>>    !=================================================================
!>>    ! HISTORY (aka netCDF diagnostics)
!>>    !=================================================================
!>>
!>>    ! P(SO4) from gas-phase oxidation [kg S/s]
!>>    IF ( State_Diag%Archive_ProdSO4fromGasPhase ) THEN
!>>       State_Diag%ProdSO4fromGasPhase(I,J,L) = &
!>>            ( L1  * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
!>>    ENDIF
!>>
!>>    ! P(SO4) from aqueous-phase oxidation with H2O2 [kg S/s]
!>>    IF ( State_Diag%Archive_ProdSO4fromH2O2inCloud ) THEN
!>>       State_Diag%ProdSO4fromH2O2inCloud(I,J,L) = &
!>>            ( L2S * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
!>>    ENDIF
!>>
!>>    ! P(SO4) from aqueous-phase oxidation with O3 [kg S/s]
!>>    IF ( State_Diag%Archive_ProdSO4fromO3InCloud ) THEN
!>>       State_Diag%ProdSO4fromO3InCloud(I,J,L) = &
!>>            ( ( L3S + SRo3 ) * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
!>>    ENDIF
!>>
!>>    ! P(SO4) from aqueous-phase oxidation with HOBr [kg S/s]
!>>    IF ( State_Diag%Archive_ProdSO4fromHOBrInCloud ) THEN
!>>       State_Diag%ProdSO4fromHOBrInCloud(I,J,L) = &
!>>            ( ( L5S + SRhobr ) * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
!>>    ENDIF
!>>
!>>    ! P(SO4) from aqueous-phase oxidation with O2 metal-catalyzed
!>>    ! [kg S/s]
!>>    IF ( State_Diag%Archive_ProdSO4fromO2InCloudMetal ) THEN
!>>       State_Diag%ProdSO4fromO2InCloudMetal(I,J,L) = &
!>>            ( L4S * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
!>>    ENDIF
!>>
!>>    ! P(SO4) from O3 in sea salt aerosol [kg S/s]
!>>    IF ( State_Diag%Archive_ProdSO4fromO3inSeaSalt ) THEN
!>>       State_Diag%ProdSO4fromO3inSeaSalt(I,J,L) = &
!>>            ( ( PSO4E + PSO4F ) * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
!>>    ENDIF
!>>
!>>    ! P(SO4) by SRo3 [kg S/s]
!>>    IF ( State_Diag%Archive_ProdSO4fromSRO3 ) THEN
!>>       State_Diag%ProdSO4fromSRO3(I,J,L) = &
!>>            ( SRo3 * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
!>>    ENDIF
!>>
!>>    ! P(SO4) by SRhobr [kg S/s]
!>>    IF ( State_Diag%Archive_ProdSO4fromSRHOBr ) THEN
!>>       State_Diag%ProdSO4fromSRHOBr(I,J,L) = &
!>>            ( SRhobr * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
!>>    ENDIF
!>>
!>>    ! P(SO4) by o3s [kg S/s]
!>>    IF ( State_Diag%Archive_ProdSO4fromO3s ) THEN
!>>       State_Diag%ProdSO4fromO3s(I,J,L) = &
!>>            ( L3S_1 * State_Met%AD(I,J,L) / TCVV_S ) / DTCHEM
!>>    ENDIF
!>>

    ! Free pointers
    Spc     => NULL()
    SSAlk   => NULL()

  END SUBROUTINE SET_SO2
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_hplus
!
! !DESCRIPTION: Subroutine GET\_HPLUS computes H+ concentrations in cloud
!  liquid water for pH dependent cloud chemistry. (bec, 4/11/11)
!\\
!\\
! !INTERFACE:
!
    SUBROUTINE GET_HPLUS( SO4nss, HMSc, TNH3, TNO3, SO2, CL, TNA, TDCA, TFA, &
                          TAA,  T, PRES, LWC,  iHPLUS, HPLUS )
!
! !USES:
!
    USE ERROR_MOD,       ONLY : IT_IS_NAN, GEOS_CHEM_STOP
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN)    :: SO4nss ! Total nss sulfate mixing ratio [M]
    REAL(fp),  INTENT(IN)    :: HMSc ! Total HMS mixing ratio [M]
    REAL(fp),  INTENT(IN)    :: TNO3   ! Total nitrate (gas+particulate) mixing
                                       ! ratio [v/v]
    REAL(fp),  INTENT(IN)    :: TNH3   ! NH3 mixing ratio [v/v]
    REAL(fp),  INTENT(IN)    :: SO2    ! SO2 mixing ratio [v/v]
    REAL(fp),  INTENT(IN)    :: CL     ! Total chloride (gas+particulate) mixing
    REAL(fp),  INTENT(IN)    :: TNA    ! Sodium (particulate) [v/v]
    REAL(fp),  INTENT(IN)    :: TDCA   ! Total Ca2+ and Mg2+ mixing ratio [M] ! jmm 12/3/18
    REAL(fp),  INTENT(IN)    :: TAA    ! Acetic acid mixing ratio [v/v] ! jmm 12/3/18
    REAL(fp),  INTENT(IN)    :: TFA    ! Formic acid mixing ratio [v/v] ! jmm 12/3/18
    REAL(fp),  INTENT(IN)    :: T      ! Temperature [K]
    REAL(fp),  INTENT(IN)    :: PRES   ! Dry air partial ressure [atm]
    REAL(fp),  INTENT(IN)    :: LWC    ! Cloud liquid water content [m3/m3]
    REAL(fp),  INTENT(IN)    :: iHPLUS ! Initial [H+] [M]
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),  INTENT(OUT)   :: HPLUS  ! Calculated [H+] [M]
! !REMARKS:
!  Calculation:
!  ============================================================================
!  Solve the following electroneutrality equation:
!  [H+] = 2[SO4--] + [Cl-] + [OH-] + [HCO3-] + 2[CO3--] + [HSO3-] + 2[SO3--] +!
!          [NO3-] + [HCOO-] + [CH3COO-] - [Na] - 2[Ca] - [NH4]
!  Uses Newton's method to solve the equation:
!     x_1 = x_0 -f(x_0)/f'(x_0)
!     iterate until converge
!
!  Let concentrations of [HCO3], [CO3], [HSO3], [SO3], [NO3] and [NH4] evolve
!  according to Henry's law equilibrium.
!
!  To add new species:
!    - Add species not affected by HPLUS to the "D' term
!    - Add species that disassociate once using the kHNO3 and dkHNO3
!    functions
!      as a template
!    - Add species that disassociate twice using the kSO21 and dkSO21
!    functions
!      as a template for the single charged ion and kSO22 and dkSO22
!      functions for
!      the double charged ion

!  Assume [S(VI)] = [SO4]nss (this applies for pH > 3)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Water dissociation constants
    REAL(fp),  PARAMETER   :: Kw   = 1.0e-14_fp
    REAL(fp),  PARAMETER   :: DhrKw = -6710.e+0_fp
    REAL(fp),  PARAMETER   :: MINVAL = 0.01
!
! !LOCAL VARIABLES:
!
    REAL(fp)               :: D, Kw_T, ipH, newpH, nHPLUS
    REAL(fp)               :: fHCO3, fCO3
    REAL(fp)               :: fHSO3, fSO3
    REAL(fp)               :: fHNO3, fNH4, fHCl
    REAL(fp)               :: dHCO3, dCO3
    REAL(fp)               :: dHSO3, dSO3
    REAL(fp)               :: dHNO3, dNH4, dHCl
    REAL(fp)               :: fAA, fFA, dAA, dFA
    REAL(fp)               :: f, df, nnHPLUS, fCa, dCa
    INTEGER                :: count

    !=================================================================
    ! GET_HPLUS begins here!
    !=================================================================

    ! Initial pH guess
    ipH = -log10(iHPLUS)

    ! Non-volatile aerosol concentration [M]
    ! For now sulfate is the only non-volatile species
    D = (2.e+0_fp*SO4nss) - TNA - (2.e+0_fp*TDCA) + (1.e+0_fp * HMSc)

    ! Temperature dependent water equilibrium constant
    Kw_T = Kw*exp(DhrKw*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    ! Initialize
    newpH   = 0.0
    COUNT = 0

    DO WHILE ( ABS(ipH-newpH) .gt. MINVAL )

       COUNT = COUNT+1

       IF ( COUNT .EQ. 1 ) THEN
          ipH = ipH
       ELSE
          ipH = newpH
       ENDIF

       nHPLUS = 10.e+0_fp**(-ipH)

       ! Get f(x) terms
       fHCO3  = kCO21 ( PRES, T, LWC, nHPLUS )

       fCO3 = kCO22 ( PRES, T, LWC, nHPLUS )

       fHSO3  = kSO21 ( PRES, T, LWC, nHPLUS, SO2 )

       fSO3 = kSO22 ( PRES, T, LWC, nHPLUS, SO2 )

       fHNO3 = kHNO3 ( PRES, T, LWC, nHPLUS, TNO3 )

       fNH4  = kNH3  ( PRES, T, LWC, nHPLUS, TNH3, Kw_T )

       ! include HCl in cloud pH calculations, xnw 10/17/17
       fHCl  = kHCl  ( PRES, T, LWC, nHPLUS, CL  )

       fFA   = kFA   ( PRES, T, LWC, nHPLUS, TFA ) ! jmm 12/3/18

       fAA   = kAA   ( PRES, T, LWC, nHPLUS, TAA ) ! jmm 12/3/18

       ! Get f'(x) terms
       dHCO3  = dkCO21 ( PRES, T, LWC, nHPLUS )

       dCO3 = dkCO22 ( PRES, T, LWC, nHPLUS )

       dHSO3  = dkSO21 ( PRES, T, LWC, nHPLUS, SO2 )

       dSO3 = dkSO22 ( PRES, T, LWC, nHPLUS, SO2 )

       dHNO3 = dkHNO3 ( PRES, T, LWC, nHPLUS, TNO3 )

       dNH4  = dkNH3  ( PRES, T, LWC, nHPLUS, TNH3, Kw_T )

       dHCl = dkHCl ( PRES, T, LWC, nHPLUS, CL )

       dFA   = dkFA   ( PRES, T, LWC, nHPLUS, TFA ) ! jmm 12/3/18

       dAA   = dkAA   ( PRES, T, LWC, nHPLUS, TAA ) ! jmm 12/3/18
       ! Calculate [Ca2+] in equilibrium with CaCO3(s)
       CALL CaCO3_PRECIP ( PRES, T, nHPLUS, fCa, dCa )

       ! if [Ca2+] in equilibrium with CacO3(s) is greater than total [Ca2+]
       ! then all Ca is dissolved else [Ca2+] varies with [H+]
       IF ( fCa .ge. TDCA ) THEN
          ! Non-volatile aerosol concentration [M]
          D = (2.e+0_fp*SO4nss) - (TNA+2.e+0_fp*TDCA)

          ! Define f(x)
          f = D - nHPLUS + Kw/nHPLUS + fHCO3 + 2.e+0_fp * &
               fCO3 + fHSO3 + 2.e+0_fp * fSO3 + fHNO3 - fNH4 + &
               fHCl + fFA + fAA

          ! Define f'(x)
          df = - 1.d0 - Kw/nHPLUS/nHPLUS + dHCO3 + 2.e+0_fp * &
               dCO3 + dHSO3 + 2.e+0_fp * dSO3 + dHNO3 - dNH4 + &
               dHCl + dFA + dAA

       ELSE
          ! Non-volatile aerosol concentration [M]
          D = (2.e+0_fp * SO4nss) - TNA

          ! Define f(x)
          f = D - nHPLUS + Kw/nHPLUS + fHCO3 + 2.e+0_fp * fCO3 + &
               fHSO3 + 2.e+0_fp * fSO3 + fHNO3 - fNH4 + &
               fHCl + fFA + fAA - 2.e+0_fp * fCa
          ! Define f'(x)
          df = - 1.d0 - Kw/nHPLUS/nHPLUS + dHCO3 + 2.e+0_fp * dCO3 + &
               dHSO3 + 2.e+0_fp * dSO3 + dHNO3 - dNH4 + &
               dHCl + dFA + dAA - 2.e+0_fp * dCa
       ENDIF

       ! Apply Newton's method
       nnHPLUS = nHPLUS - f/df

       ! Set minimum [H+] = 1.d-14 (pH = 14)
       nnHPLUS = MAX(nnHPLUS,1.0e-14_fp)

       ! Set maximum [H+] = 1.d-1 (pH = 1)
       nnHPLUS = MIN(nnHPLUS,1.0e-1_fp)

       ! If solution does not converge after 50 iterations
       ! average last 2 pH calculations
       IF (count > 50) THEN
          newpH = ((-log10(nnHPLUS)) + (-log10(nHPLUS))) / 2.0e+0_fp

          IF (IT_IS_NAN( newpH )) THEN
             write(6,*) 'newpH = ', newpH
             write(6,*) 'nnHPLUS = ', nnHPLUS
             write(6,*) 'nHPLUS = ', nHPLUS
             CALL GEOS_CHEM_STOP
          ENDIF

          EXIT
       ELSE
          newpH = -log10(nnHPLUS)

          IF (IT_IS_NAN( newpH )) THEN
             write(6,*) 'newpH = ', newpH
             write(6,*) 'nnHPLUS = ', nnHPLUS
             CALL GEOS_CHEM_STOP
          ENDIF

       ENDIF

    ENDDO

    HPLUS = 10.0e+0_fp**(-newpH)

  END SUBROUTINE GET_HPLUS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kCO21
!
! !DESCRIPTION: Function kCO21
!\\
!\\
! !INTERFACE:
!
  FUNCTION kCO21 ( P, T, LWC, HPLUS ) RESULT ( KCO2p )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KCO2p, KCO2p2
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! CO2 dissociation constants
    REAL(fp),  PARAMETER  :: Kc1 = 4.3e-7
    REAL(fp),  PARAMETER  :: Kc2 =4.68e-11
    REAL(fp),  PARAMETER  :: DhrKc1 = -1000.
    REAL(fp),  PARAMETER  :: DhrKc2 = -1760.
    REAL(fp),  PARAMETER  :: Hco2 = 3.4e-2
    REAL(fp),  PARAMETER  :: Dhco2 = 2.44e+3_fp
    ! CO2 concentration [v/v]
    REAL(fp),  PARAMETER  :: CO2 = 390.0e-6_fp
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: Hco2_T, Kc1_T, Kc2_T
    REAL(fp)              :: Hco2eff, xCO2, pCO2

    !=================================================================
    ! kCO21 begins here!
    !=================================================================

    !CO2 dissolution constants
    Hco2_T  = Hco2*exp(Dhco2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Kc1_T   = Kc1*exp(DhrKc1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Kc2_T   = Kc2*exp(DhrKc2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    !CO2 dissolution
    Hco2eff = Hco2_T*(1.e+0_fp+(Kc1_T/HPLUS)+((Kc1_T*Kc2_T)/(HPLUS*HPLUS)))
    xCO2    = 1.e+0_fp / ( 1.e+0_fp + ( Hco2eff * 0.08205e+0_fp * T * LWC ) )
    pCO2    = CO2 * P * xCO2

    KCO2p  = Hco2_T / HPLUS * Kc1_T * pCO2

  END FUNCTION kCO21
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkCO21
!
! !DESCRIPTION: Function dkCO21
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkCO21 ( P, T, LWC, HPLUS ) RESULT ( KCO2p )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KCO2p, KCO2p2
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  22 Mar 2017 - M. Sulprizio- Dhco2 value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (V. Shah)
!  15 Feb 2019 - J. Moch     - updated function to make output
!  derivative of [HCO3-]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! CO2 dissociation constants
      REAL(fp),  PARAMETER  :: Kc1 = 4.3e-7
      REAL(fp),  PARAMETER  :: Kc2 =4.68e-11
      REAL(fp),  PARAMETER  :: DhrKc1 = -1000.
      REAL(fp),  PARAMETER  :: DhrKc2 = -1760.
      REAL(fp),  PARAMETER  :: Hco2 = 3.4e-2
      REAL(fp),  PARAMETER  :: Dhco2 = 2.44e+3_fp
      ! CO2 concentration [v/v]
      REAL(fp),  PARAMETER  :: CO2 = 390.0e-6_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hco2_T, Kc1_T, Kc2_T

! !REMARKS:

      !=================================================================
      ! dkCO21 begins here!
      !=================================================================

      !CO2 dissolution constants
      Hco2_T = Hco2*exp(Dhco2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc1_T = Kc1*exp(DhrKc1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc2_T = Kc2*exp(DhrKc2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      !CO2 dissolution

      KCO2p  = Kc1_T * Hco2_T * CO2 * P * ( Kc1_T * Kc2_T * Hco2_T * &
          0.08205e+0_fp * T * LWC - Hco2_T * 0.08205e+0_fp * T *     &
          LWC * HPLUS * HPLUS - HPLUS * HPLUS) / (Kc1_T * Kc2_T *    &
          Hco2_T * 0.08205e+0_fp * T * LWC + Kc1_T * Hco2_T *        &
          0.08205e+0_fp * T * LWC * HPLUS + Hco2_T * 0.08205e+0_fp * &
          T * LWC * HPLUS * HPLUS + HPLUS * HPLUS)**2

      END FUNCTION dkCO21
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kCO22
!
! !DESCRIPTION: Function kCO22
!\\
!\\
! !INTERFACE:
!
  FUNCTION kCO22 ( P, T, LWC, HPLUS ) RESULT ( KCO2p2 )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KCO2p, KCO2p2
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! CO2 dissociation constants
    REAL(fp),  PARAMETER  :: Kc1 = 4.3e-7
    REAL(fp),  PARAMETER  :: Kc2 =4.68e-11
    REAL(fp),  PARAMETER  :: DhrKc1 = -1000.
    REAL(fp),  PARAMETER  :: DhrKc2 = -1760.
    REAL(fp),  PARAMETER  :: Hco2 = 3.4e-2
    REAL(fp),  PARAMETER  :: Dhco2 = 2.44e+3_fp
    ! CO2 concentration [v/v]
    REAL(fp),  PARAMETER  :: CO2 = 390.0e-6_fp
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: Hco2_T, Kc1_T, Kc2_T
    REAL(fp)              :: Hco2eff, xCO2, pCO2

    !=================================================================
    ! kCO22 begins here!
    !=================================================================

    !CO2 dissolution constants
    Hco2_T  = Hco2*exp(Dhco2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Kc1_T   = Kc1*exp(DhrKc1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Kc2_T   = Kc2*exp(DhrKc2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    !CO2 dissolution
    Hco2eff = Hco2_T*(1.e+0_fp+(Kc1_T/HPLUS)+((Kc1_T*Kc2_T)/(HPLUS*HPLUS)))
    xCO2    = 1.e+0_fp / ( 1.e+0_fp  + ( Hco2eff * 0.08205e+0_fp * T * LWC ) )
    pCO2    = CO2 * P * xCO2

    KCO2p2 = Kc1_T * Kc2_T * Hco2_T / HPLUS / HPLUS * pCO2

  END FUNCTION kCO22
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkCO22
!
! !DESCRIPTION: Function dkCO22
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkCO22 ( P, T, LWC, HPLUS ) RESULT ( KCO2p2 )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KCO2p, KCO2p2
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  22 Mar 2017 - M. Sulprizio- Dhco2 value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (V. Shah)
!  15 Feb 2019 - J. Moch     - updated function to make output deriviate
!  of [CO3--]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! CO2 dissociation constants
      REAL(fp),  PARAMETER  :: Kc1 = 4.3e-7
      REAL(fp),  PARAMETER  :: Kc2 =4.68e-11
      REAL(fp),  PARAMETER  :: DhrKc1 = -1000.
      REAL(fp),  PARAMETER  :: DhrKc2 = -1760.
      REAL(fp),  PARAMETER  :: Hco2 = 3.4e-2
      REAL(fp),  PARAMETER  :: Dhco2 = 2.44e+3_fp
      ! CO2 concentration [v/v]
      REAL(fp),  PARAMETER  :: CO2 = 390.0e-6_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hco2_T, Kc1_T, Kc2_T

      !=================================================================
      ! dkCO22 begins here!
      !=================================================================

      !CO2 dissolution constants
      Hco2_T = Hco2*exp(Dhco2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc1_T = Kc1*exp(DhrKc1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc2_T = Kc2*exp(DhrKc2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      !CO2 dissolution

      KCO2p2 = -1.e+0_fp * Kc1_T * Kc2_T * Hco2_T * CO2 * P * ( Kc1_T * &
           Hco2_T * 0.08205e+0_fp * T * LWC + 2.0e+0_fp * Hco2_T *      &
           0.08205e+0_fp * T * LWC * HPLUS + 2.0e+0_fp * HPLUS ) /      &
           ( Kc1_T * Kc2_T * Hco2_T * 0.08205e+0_fp * T * LWC +         &
           Kc1_T * Hco2_T * 0.08205e+0_fp * T * LWC * HPLUS +           &
           Hco2_T *0.08205e+0_fp * T * LWC * HPLUS * HPLUS +            &
           HPLUS * HPLUS )**2

      END FUNCTION dkCO22
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kSO21
!
! !DESCRIPTION: Function kSO21
!\\
!\\
! !INTERFACE:
!
  FUNCTION kSO21 ( P, T, LWC, HPLUS, SO2 ) RESULT ( KSO2p )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, SO2
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KSO2p, KSO2p2
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! SO2 dissociation constants
    REAL(fp),  PARAMETER  :: Ks1 = 1.3e-2
    REAL(fp),  PARAMETER  :: Ks2 = 6.6e-8
    REAL(fp),  PARAMETER  :: Hso2 = 1.23
    REAL(fp),  PARAMETER  :: Dhso2 = 3.14e+3_fp
    REAL(fp),  PARAMETER  :: DhrKso21 = 1960.
    REAL(fp),  PARAMETER  :: DhrKso22 = 1500.
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: Hso2_T, Ks1_T, Ks2_T
    REAL(fp)              :: Hso2eff, xSO2, pSO2

    !=================================================================
    ! kSO21 begins here!
    !=================================================================

    ! SO2 dissolution constants
    Hso2_T  = Hso2*exp(Dhso2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Ks1_T   = Ks1*exp(DhrKso21*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Ks2_T   = Ks2*exp(DhrKso22*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    ! SO2 dissolution
    Hso2eff = Hso2_T*(1.e+0_fp+(Ks1_T/HPLUS)+((Ks1_T*Ks2_T)/(HPLUS*HPLUS)))
    xSO2    = 1.e+0_fp / ( 1.e+0_fp  + ( Hso2eff * 0.08205e+0_fp * T * LWC ) )
    pSO2    = SO2 * P * xSO2

    KSO2p   = Hso2_T * Ks1_T * pSO2 / HPLUS

  END FUNCTION kSO21
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkSO21
!
! !DESCRIPTION: Function dkSO21
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkSO21 ( P, T, LWC, HPLUS, SO2 ) RESULT ( KSO2p )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, SO2
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KSO2p, KSO2p2
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  22 Mar 2017 - M. Sulprizio- Dhso2 value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (V. Shah)
!  15 Feb 2019 - J. Moch     - updated function to make output
!  derivative of [HSO3-]

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! SO2 dissociation constants
      REAL(fp),  PARAMETER  :: Ks1 = 1.3e-2
      REAL(fp),  PARAMETER  :: Ks2 = 6.6e-8
      REAL(fp),  PARAMETER  :: Hso2 = 1.23
      REAL(fp),  PARAMETER  :: Dhso2 = 3.14e+3_fp
      REAL(fp),  PARAMETER  :: DhrKso21 = 1960.
      REAL(fp),  PARAMETER  :: DhrKso22 = 1500.
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hso2_T, Ks1_T, Ks2_T

      !=================================================================
      ! dkSO21 begins here!
      !=================================================================



      ! SO2 dissolution constants
      Hso2_T = Hso2*exp(Dhso2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Ks1_T = Ks1*exp(DhrKso21*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Ks2_T = Ks2*exp(DhrKso22*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))


      KSO2p  = Ks1_T * Hso2_T * SO2 * P * ( Ks1_T * Ks2_T * Hso2_T *  &
           0.08205e+0_fp * T * LWC - Hso2_T * 0.08205e+0_fp * T *     &
           LWC * HPLUS * HPLUS - HPLUS * HPLUS) / (Ks1_T * Ks2_T *    &
           Hso2_T * 0.08205e+0_fp * T * LWC + Ks1_T * Hso2_T *        &
           0.08205e+0_fp * T * LWC * HPLUS + Hso2_T * 0.08205e+0_fp * &
           T * LWC * HPLUS * HPLUS + HPLUS * HPLUS)**2

      END FUNCTION dkSO21
!EOC

!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kSO22
!
! !DESCRIPTION: Function kSO22
!\\
!\\
! !INTERFACE:
!
  FUNCTION kSO22 ( P, T, LWC, HPLUS, SO2 ) RESULT ( KSO2p2 )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, SO2
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KSO2p, KSO2p2
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! SO2 dissociation constants
    REAL(fp),  PARAMETER  :: Ks1 = 1.3e-2
    REAL(fp),  PARAMETER  :: Ks2 = 6.6e-8
    REAL(fp),  PARAMETER  :: Hso2 = 1.23
    REAL(fp),  PARAMETER  :: Dhso2 = 3.14e+3_fp
    REAL(fp),  PARAMETER  :: DhrKso21 = 1960.
    REAL(fp),  PARAMETER  :: DhrKso22 = 1500.
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: Hso2_T, Ks1_T, Ks2_T
    REAL(fp)              :: Hso2eff, xSO2, pSO2

    !=================================================================
    ! kSO22 begins here!
    !=================================================================

    ! SO2 dissolution constants
    Hso2_T  = Hso2*exp(Dhso2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Ks1_T   = Ks1 *exp(DhrKso21*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Ks2_T   = Ks2 *exp(DhrKso22*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    !SO2 dissolution
    Hso2eff = Hso2_T*(1.e+0_fp+(Ks1_T/HPLUS)+((Ks1_T*Ks2_T)/(HPLUS*HPLUS)))
    xSO2    = 1.e+0_fp / ( 1.e+0_fp + ( Hso2eff * 0.08205e+0_fp * T * LWC ) )
    pSO2    = SO2 * P * xSO2

    KSO2p2 = Ks1_T * Ks2_T * Hso2_T / HPLUS / HPLUS * pSO2

  END FUNCTION kSO22
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkSO22
!
! !DESCRIPTION: Function dkSO22
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkSO22 ( P, T, LWC, HPLUS, SO2 ) RESULT ( KSO2p2 )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, SO2
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KSO2p, KSO2p2
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  22 Mar 2017 - M. Sulprizio- Dhso2 value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (V. Shah)
!  15 Feb 2019 - J. Moch     - updated function to make output
!  derivative [SO3--]

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! SO2 dissociation constants
      REAL(fp),  PARAMETER  :: Ks1 = 1.3e-2
      REAL(fp),  PARAMETER  :: Ks2 = 6.6e-8
      REAL(fp),  PARAMETER  :: Hso2 = 1.23
      REAL(fp),  PARAMETER  :: Dhso2 = 3.14e+3_fp
      REAL(fp),  PARAMETER  :: DhrKso21 = 1960.
      REAL(fp),  PARAMETER  :: DhrKso22 = 1500.
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hso2_T, Ks1_T, Ks2_T

      !=================================================================
      ! dkSO22 begins here!
      !=================================================================
      ! SO2 dissolution constants
      Hso2_T = Hso2*exp(Dhso2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Ks1_T  = Ks1 *exp(DhrKso21*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Ks2_T  = Ks2 *exp(DhrKso22*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      KSO2p2 = -1.e+0_fp * Ks1_T * Ks2_T * Hso2_T * SO2 * P * ( Ks1_T * &
           Hso2_T * 0.08205e+0_fp * T * LWC + 2.0e+0_fp * Hso2_T *      &
           0.08205e+0_fp * T * LWC * HPLUS + 2.0e+0_fp * HPLUS ) /      &
           ( Ks1_T * Ks2_T * Hso2_T * 0.08205e+0_fp * T * LWC +         &
           Ks1_T * Hso2_T * 0.08205e+0_fp * T * LWC * HPLUS +           &
           Hso2_T *0.08205e+0_fp * T * LWC * HPLUS * HPLUS +            &
           HPLUS * HPLUS )**2

      END FUNCTION dkSO22
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kHNO3
!
! !DESCRIPTION: Function kNO3
!\\
!\\
! !INTERFACE:
!
  FUNCTION kHNO3 ( P, T, LWC, HPLUS, HNO3 ) RESULT ( KHNO3p )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, HNO3
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KHNO3p
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! HNO3 dissociation constants
    REAL(fp),  PARAMETER  :: Kn1 = 15.4
    REAL(fp),  PARAMETER  :: Hhno3 = 2.1e5
    REAL(fp),  PARAMETER  :: Dhhno3 = 0.
    REAL(fp),  PARAMETER  :: DhrKn1 = 8700.
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: Hhno3_T, Kn1_T
    REAL(fp)              :: Hhno3eff, xHNO3, pHNO3

    !=================================================================
    ! kHNO3 begins here!
    !=================================================================

    ! HNO3 dissolution constants
    Hhno3_T  = Hhno3*exp(Dhhno3*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Kn1_T    = Kn1*exp(DhrKn1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    ! HNO3 dissolution
    ! The original Hhno3eff expression is valid for 298K (Seinfeld and Pandis
    ! 2006, pp 299-301), and Kn1 has a strong temperature dependence. The
    ! fix follows Eq. 7.59 of Seinfeld and Pandis (2006, pp 301).
    !Hhno3eff = 3.2e6/HPLUS
    Hhno3eff = Hhno3_T*(1.0e+0_fp+(Kn1_T/HPLUS))
    xHNO3    = 1.e+0_fp / ( 1.e+0_fp + ( Hhno3eff * 0.08205e+0_fp * T * LWC ) )
    pHNO3    = HNO3 * P * xHNO3

    kHNO3p = Hhno3_T * Kn1_T * pHNO3 / HPLUS

  END FUNCTION kHNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkHNO3
!
! !DESCRIPTION: Function dkNO3
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkHNO3 ( P, T, LWC, HPLUS, HNO3 ) RESULT ( KHNO3p )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, HNO3
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KHNO3p
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  22 Mar 2017 - M. Sulprizio- Add fix for Hhno3eff from V. Shah
!  15 Feb 2019 - J. Moch     - updated function to make output
!  derivative of [HNO3-]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! HNO3 dissociation constants
      REAL(fp),  PARAMETER  :: Kn1 = 15.4
      REAL(fp),  PARAMETER  :: Hhno3 = 2.1e5
      REAL(fp),  PARAMETER  :: Dhhno3 = 0.
      REAL(fp),  PARAMETER  :: DhrKn1 = 8700.
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hhno3_T, Kn1_T

      !=================================================================
      ! dkHNO3 begins here!
      !=================================================================

      ! HNO3 dissolution constants
      Hhno3_T = Hhno3*exp(Dhhno3*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kn1_T = Kn1*exp(DhrKn1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      ! HNO3 dissolution
      ! The original Hhno3eff expression is valid for 298K (Seinfeld and
      ! Pandis
      ! 2006, pp 299-301), and Kn1 has a strong temperature dependence.
      ! The
      ! fix follows Eq. 7.59 of Seinfeld and Pandis (2006, pp 301).

      kHNO3p = -1.0e+0_fp * Kn1_T * Hhno3_T * HNO3 * P * &
          ( 1.0e+0_fp + Hhno3_T * 0.08205e+0_fp * T * LWC ) / &
          ( Kn1_T * Hhno3_T * 0.08205e+0_fp * T * LWC + &
          Hhno3_T * 0.08205e+0_fp * T * LWC * HPLUS + &
          HPLUS )**2

      END FUNCTION dkHNO3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kHCl
!
! !DESCRIPTION: Function kHCl
!\\
!\\
! !INTERFACE:
!
  FUNCTION kHCl ( P, T, LWC, HPLUS, Cl ) RESULT ( KHClp )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, Cl
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KHClp
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! HNO3 dissociation constants
    REAL(fp),  PARAMETER  :: Kcl = 1.74e+6_fp
    REAL(fp),  PARAMETER  :: Hcl = 1.5e+3_fp
    REAL(fp),  PARAMETER  :: Dhcl = 2.3e+3_fp
    REAL(fp),  PARAMETER  :: DhrKcl = 6900.e+0_fp
!
! !LOCAL VARIABLES:
!
    REAL(fp)              :: Hcl_T, Kcl_T
    REAL(fp)              :: Hcleff, xCl, pHCl

    !=================================================================
    ! kHCl begins here!
    !=================================================================

    ! HCl dissolution constants
    HCl_T  = Hcl*exp(Dhcl*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))
    Kcl_T  = Kcl*exp(DhrKcl*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))

    !HCl dissolution
    Hcleff = Hcl_T*(1.0e+0_fp+(Kcl_T/HPLUS))
    xCl    = 1.0e+0_fp / ( 1.0e+0_fp + ( Hcleff * 0.08205e+0_fp * T * LWC ) )
    pHCl   = Cl * P * xCl

    kHClp  = Hcl_T * Kcl_T * pHCl / HPLUS

  END FUNCTION kHCl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkHCl
!
! !DESCRIPTION: Function dkHCl
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkHCl ( P, T, LWC, HPLUS, Cl ) RESULT ( KHClp )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, Cl
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KHClp
!
! !REVISION HISTORY:
!  03 Apr 2019 - X. Wang    - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! HCl dissociation constants
      REAL(fp),  PARAMETER  :: Kcl = 1.74e+6_fp
      REAL(fp),  PARAMETER  :: Hcl = 1.5e+3_fp
      REAL(fp),  PARAMETER  :: Dhcl = 2.3e+3_fp
      REAL(fp),  PARAMETER  :: DhrKcl = 6900.e+0_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hcl_T, Kcl_T

      !=================================================================
      ! dkHCl begins here!
      !=================================================================

      ! HCl dissolution constants
      Hcl_T = Hcl*exp(Dhcl*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))
      Kcl_T = Kcl*exp(DhrKcl*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))

      ! HCl dissolution
      ! The fix follows Eq. 7.59 of Seinfeld and Pandis (2006, pp 301).

      kHClp = -1.0e+0_fp * Kcl_T * Hcl_T * Cl * P *          &
           ( 1.0e+0_fp + Hcl_T * 0.08205e+0_fp * T * LWC ) / &
           ( Kcl_T * Hcl_T * 0.08205e+0_fp * T * LWC +       &
           Hcl_T * 0.08205e+0_fp * T * LWC * HPLUS +         &
           HPLUS )**2

      END FUNCTION dkHCl
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kNH3
!
! !DESCRIPTION: Function kNH3
!\\
!\\
! !INTERFACE:
!
  FUNCTION kNH3 ( P, T, LWC, HPLUS, NH3, Kw ) RESULT ( KNH3p )
!
! !INPUT PARAMETERS:
!
    REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, NH3, Kw
!
! !OUTPUT PARAMETERS:
!
    REAL(fp)              :: KNH3p
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! NH3 dissociation contants
    REAL(fp),  PARAMETER  :: Ka1 = 1.7e-5
    REAL(fp),  PARAMETER  :: Hnh3 = 60.
    REAL(fp),  PARAMETER  :: Dhnh3 = 4200e+0_fp
    REAL(fp),  PARAMETER  :: DhrKa1 = -450.

    ! Variables
    REAL(fp)              :: Hnh3_T, Ka1_T
    REAL(fp)              :: Hnh3eff, xNH3, pNH3

    !=================================================================
    ! kNH3 begins here!
    !=================================================================

    !NH3 dissolution constants
    Hnh3_T  = Hnh3*exp(Dhnh3*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
    Ka1_T   = Ka1*exp(DhrKa1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

    !NH3 dissolution
    Hnh3eff = Hnh3_T*(1.e+0_fp+((Ka1_T* HPLUS) / Kw))
    xNH3    = 1.e+0_fp / ( 1.e+0_fp + ( Hnh3eff * 0.08205e+0_fp * T * LWC ) )
    pNH3    = NH3 * P * xNH3

    KNH3p   = HPLUS * Hnh3_T * Ka1_T * pNH3 / Kw

  END FUNCTION kNH3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkNH3
!
! !DESCRIPTION: Function dkNH3
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkNH3 ( P, T, LWC, HPLUS, NH3, Kw ) RESULT ( KNH3p )
!

! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, NH3, Kw
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KNH3p
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  22 Mar 2017 - M. Sulprizio- Dhnh3 value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (V. Shah)
!  15 Feb 2019 - J. Moch     - updated function to make output
!  derivative of [NH4+]
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! NH3 dissociation contants
      REAL(fp),  PARAMETER  :: Ka1 = 1.7e-5
      REAL(fp),  PARAMETER  :: Hnh3 = 60.
      REAL(fp),  PARAMETER  :: Dhnh3 = 4200e+0_fp
      REAL(fp),  PARAMETER  :: DhrKa1 = -450.

      ! Variables
      REAL(fp)              :: Hnh3_T, Ka1_T

      !=================================================================
      ! dkNH3 begins here!
      !=================================================================

      !NH3 dissolution constants
      Hnh3_T = Hnh3*exp(Dhnh3*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Ka1_T = Ka1*exp(DhrKa1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      !NH3 dissolutionnyn

      KNH3p = Ka1_T * Hnh3_T * NH3 * Kw * P * ( 1.0e+0_fp +    &
           Hnh3_T * 0.08205e+0_fp * T * LWC ) /                &
           ( Hnh3_T * 0.08205e+0_fp * T * LWC * ( Kw + Ka1_T * &
           HPLUS ) + Kw)**2

      END FUNCTION dkNH3
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kFA
!
! !DESCRIPTION: Function kFA
!\\
!\\
! !INTERFACE:
!
      FUNCTION kFA ( P, T, LWC, HPLUS, FA ) RESULT ( kFAp )
!

! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, FA
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KFAp
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  17 Oct 2017 - M. Sulprizio- Dhck value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (Qianjie Chen)
!  03 Dec 2018 - J. Moch     - Modified for formic acid (HCOOH). Values
!                              taken from Sienfeld and Pandis. Made it
!                              to output is [FA]
!  01 May 2020 - V. Shah     - Use correct equilibrium constants
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! HCOOH dissociation constants
      REAL(fp),  PARAMETER  :: Kformate = 1.8e-4_fp ! equib const
      REAL(fp),  PARAMETER  :: Hfa = 8800e+0_fp ! henry const
      REAL(fp),  PARAMETER  :: Dhfa = 6100e+0_fp ! henry temp
      REAL(fp),  PARAMETER  :: DhrKfa = 151.e+0_fp ! equib temp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hfa_T, Kfa_T
      REAL(fp)              :: Hfaeff, xFA, pFA

      !=================================================================
      ! kFA begins here!
      !=================================================================

      ! Formic acid dissolution constants
      HFA_T = Hfa*exp(Dhfa*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))
      Kfa_T = Kformate*exp(DhrKfa*((1.0e+0_fp/T) &
           - (1.0e+0_fp/298.0e+0_fp)))

      !HCOOH  dissolution
      Hfaeff = Hfa_T*(1.0e+0_fp+(Kfa_T/HPLUS))
      xFA = 1.0e+0_fp / ( 1.0e+0_fp &
          + ( Hfaeff * 0.08205e+0_fp * T * LWC ) )
      pFA = FA * P * xFA

      kFAp = Hfa_T * Kfa_T * pFA / HPLUS

      END FUNCTION kFA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkFA
!
! !DESCRIPTION: Function dkFA
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkFA ( P, T, LWC, HPLUS, FA ) RESULT ( kFAp )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, FA
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KFAp
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  17 Oct 2017 - M. Sulprizio- Dhck value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (Qianjie Chen)
!  03 Dec 2018 - J. Moch     - Modified for formic acid (HCOOH). Values
!  taken from
!                              Sienfeld and Pandis. Made it to output is
!                              [FA]
!  01 May 2020 - V. Shah     - Use correct equilibrium constants
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! HCOOH dissociation constants
      REAL(fp),  PARAMETER  :: Kformate = 1.8e-4_fp ! equib const
      REAL(fp),  PARAMETER  :: Hfa = 8800e+0_fp ! henry const
      REAL(fp),  PARAMETER  :: Dhfa = 6100e+0_fp ! henry temp
      REAL(fp),  PARAMETER  :: DhrKfa = 151.e+0_fp ! equib temp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Hfa_T, Kfa_T

      !=================================================================
      ! dkFA begins here!
      !=================================================================

      ! Formic acid dissolution constants
      HFA_T = Hfa*exp(Dhfa*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))
      Kfa_T = Kformate*exp(DhrKfa*((1.0e+0_fp/T) &
           - (1.0e+0_fp/298.0e+0_fp)))

      !HCOOH  dissolution

      kFAp = -1.0e+0_fp * Kfa_T * HFA_T * FA * P *           &
           ( 1.0e+0_fp + HFA_T * 0.08205e+0_fp * T * LWC ) / &
           ( Kfa_T * HFA_T * 0.08205e+0_fp * T * LWC +       &
           HFA_T * 0.08205e+0_fp * T * LWC * HPLUS +         &
           HPLUS )**2

      END FUNCTION dkFA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: kAA
!
! !DESCRIPTION: Function kAA
!\\
!\\
! !INTERFACE:
!
      FUNCTION kAA ( P, T, LWC, HPLUS, AA ) RESULT ( kAAp )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, AA
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KAAp
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  17 Oct 2017 - M. Sulprizio- Dhck value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (Qianjie Chen)
!  03 Dec 2018 - J. Moch     - Modified for acetic acid (CH3COOH).
!  Values taken from
!                              Sienfeld and Pandis, value of [HCOOH]
!  01 May 2020 - V. Shah     - Use correct equilibrium constants
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! CH3HCOOH dissociation constants
      REAL(fp),  PARAMETER  :: Kacetate = 1.75e-5_fp
      REAL(fp),  PARAMETER  :: Haa = 4100e+0_fp
      REAL(fp),  PARAMETER  :: Dhaa = 6200e+0_fp
      REAL(fp),  PARAMETER  :: DhrKaa = 50.0e+0_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Haa_T, Kaa_T
      REAL(fp)              :: Haaeff, xAA, pAA
      !=================================================================
      ! kAA begins here!
      !=================================================================

      ! Formic acid dissolution constants
      HAA_T = Haa*exp(Dhaa*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))
      Kaa_T = Kacetate*exp(DhrKaa*((1.0e+0_fp/T) &
          - (1.0e+0_fp/298.0e+0_fp)))

      !HCOOH  dissolution
      Haaeff = Haa_T*(1.0e+0_fp+(Kaa_T/HPLUS))
      xAA = 1.0e+0_fp / ( 1.0e+0_fp &
         + ( Haaeff * 0.08205e+0_fp * T * LWC ) )
      pAA = AA * P * xAA

      kAAp = Haa_T * Kaa_T * pAA / HPLUS

      END FUNCTION kAA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: dkAA
!
! !DESCRIPTION: Function kdAA
!\\
!\\
! !INTERFACE:
!
      FUNCTION dkAA ( P, T, LWC, HPLUS, AA ) RESULT ( kAAp )
!
! !INPUT PARAMETERS:
!
      REAL(fp),  INTENT(IN) :: T, P, LWC, HPLUS, AA
!
! !OUTPUT PARAMETERS:
!
      REAL(fp)              :: KAAp
!
! !REVISION HISTORY:
!  25 Jan 2012 - M. Payer    - Added ProTeX headers
!  28 Apr 2015 - E. Lundgren - Input pressure is now dry air partial
!  pressure
!  17 Oct 2017 - M. Sulprizio- Dhck value is from Table 7.3 of Seinfeld
!  and
!                              Pandis (2006, pp 289) and should be
!                              positive for
!                              consistency with the way it is used here.
!                              Also,
!                              the value of R, in units of kcal mol-1
!                              K-1, is
!                              1.986x10^-3, not 0.04. (Qianjie Chen)
!  03 Dec 2018 - J. Moch     - Modified for acetic acid (CH3COOH).
!  Values taken from
!                              Sienfeld and Pandis. Output is
!                              derivative.
!  01 May 2020 - V. Shah     - Use correct equilibrium constants
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
      ! HCOOH dissociation constants
      REAL(fp),  PARAMETER  :: Kacetate = 1.75e-5_fp
      REAL(fp),  PARAMETER  :: Haa = 4100e+0_fp
      REAL(fp),  PARAMETER  :: Dhaa = 6200e+0_fp
      REAL(fp),  PARAMETER  :: DhrKaa = 50.0e+0_fp
!
! !LOCAL VARIABLES:
!
      REAL(fp)              :: Haa_T, Kaa_T
      !=================================================================
      ! kAA begins here!
      !=================================================================

      ! Formic acid dissolution constants
      HAA_T = Haa*exp(Dhaa*((1.0e+0_fp/T)-(1.0e+0_fp/298.0e+0_fp)))
      Kaa_T = Kacetate*exp(DhrKaa*((1.0e+0_fp/T) &
          - (1.0e+0_fp/298.0e+0_fp)))

      !HCOOH  dissolution
      kAAp =  -1.0e+0_fp * Kaa_T * HAA_T * AA * P * &
          ( 1.0e+0_fp + HAA_T * 0.08205e+0_fp * T * LWC ) / &
          ( Kaa_T * HAA_T * 0.08205e+0_fp * T * LWC + &
          HAA_T * 0.08205e+0_fp * T * LWC * HPLUS + &
          HPLUS )**2


      END FUNCTION dkAA
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CaCO3_PRECIP
!
! !DESCRIPTION: Subroutine CaCO3 to calculate [Ca++] in equilibrium with
! CaCO3(s) (dust particles) depending on [H+]
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CaCO3_PRECIP ( P,  T, HPLUS, fCa, dCa )
!
! !INPUT PARAMETERS:
!
      REAL(fp),        INTENT(IN) :: T, P, HPLUS
!
! !OUTPUT PARAMETERS:
!
      REAL(fp),  INTENT(OUT):: fCa, dCa ! [Ca2+] and d([Ca2+])/d[H+]
!
! !REVISION HISTORY:
!  25 Dec 2019 - V. Shah - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
! !DEFINED PARAMETERS:
!
      REAL(fp),  PARAMETER  :: Kc1 = 4.3e-7_fp
      REAL(fp),  PARAMETER  :: Kc2 = 4.68e-11_fp
      REAL(fp),  PARAMETER  :: DhrKc1 = -1000.
      REAL(fp),  PARAMETER  :: DhrKc2 = -1760.
      REAL(fp),  PARAMETER  :: Hco2 = 3.4e-2_fp
      REAL(fp),  PARAMETER  :: Dhco2 = 2.44e+3_fp
      ! CO2 concentration [v/v]
      REAL(fp),  PARAMETER  :: CO2 = 390.0e-6_fp
      REAL(fp),  PARAMETER  :: Ksp = 3.3e-9_fp
      REAL(fp),  PARAMETER  :: DHrKsp = -1200e+0_fp

! !LOCAL VARIABLES:
      REAL(fp)              :: HCO2_T, Kc1_T, Kc2_T, Ksp_T

! !REMARKS:

      !=================================================================
      ! CaCO3_PRECIP begins here!
      !=================================================================
      !Temperature adjusted eq. constants
      !CO2 dissolution constants
      Hco2_T = Hco2*exp(Dhco2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc1_T = Kc1*exp(DhrKc1*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))
      Kc2_T = Kc2*exp(DhrKc2*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      ! CaCO3 eq constants
      Ksp_T = Ksp*exp(DhrKsp*((1.e+0_fp/T)-(1.e+0_fp/298.e+0_fp)))

      !Ca concentrations [M]
      fCa = Ksp_T * HPLUS * HPLUS / (Kc1_T * Kc2_T * Hco2_T * CO2 * P)
      !derivative d[Ca2+]/dH+
      dCa  = 2e+0_fp * Ksp_T * HPLUS / (Kc1_T * Kc2_T * Hco2_T * CO2 * P)

      END SUBROUTINE CaCO3_PRECIP
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: aqchem_so2
!
! !DESCRIPTION: Subroutine AQCHEM\_SO2 computes the reaction rates for aqueous
! SO2 chemistry. (rjp, bmy, 10/31/02, 12/12/02)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE AQCHEM_SO2( I,      J,       L,      LWC,     T,      P,        &
                         SO2,    H2O2,    O3,     HCHO,    Hplus,  MnII,     &
                         FeIII,  KaqH2O2, KaqO3,  KaqO3_1, KaqO2,  HSO3aq,   &
                         SO3aq,  KaqHCHO, KaqHMS, KaqHMS2                   )
!
! !INPUT PARAMETERS:
!
    INTEGER,  INTENT(IN)  :: I, J, L ! Coordinates, for diagnostic use -- MSL

    REAL(fp), INTENT(IN)  :: LWC     ! Liq water content [m3/m3]=1.E-6*L [g/m3]
    REAL(fp), INTENT(IN)  :: T       ! Temperature [K]
    REAL(fp), INTENT(IN)  :: P       ! Dry air partial pressure [atm]
    REAL(fp), INTENT(IN)  :: SO2     ! SO2  mixing ratio [v/v]
    REAL(fp), INTENT(IN)  :: H2O2    ! H2O2 mixing ratio [v/v]
    REAL(fp), INTENT(IN)  :: O3      ! O3   mixing ratio [v/v]
    REAL(fp), INTENT(IN)  :: HPLUS   ! Concentration of H+ ion (i.e. pH) [v/v]
    REAL(fp), INTENT(IN)  :: MnII    ! Concentration of MnII [mole/l]
    REAL(fp), INTENT(IN)  :: FeIII   ! Concentration of FeIII [mole/l]
    REAL(fp), INTENT(IN)  :: HCHO    ! HCHO   mixing ratio [v/v] (jmm, 06/13/18)

!
! !OUTPUT PARAMETERS:
!
    REAL(fp), INTENT(OUT) :: KaqH2O2 ! Reaction rate for H2O2
    REAL(fp), INTENT(OUT) :: KaqO3   ! Reaction rate for O3
    REAL(fp), INTENT(OUT) :: KaqO3_1 ! only the SO3-- oxidation, (qjc, 04/10/16)
    REAL(fp), INTENT(OUT) :: KaqO2   ! Reaction rate for O2 (metal cat)
    REAL(fp), INTENT(OUT) :: KaqHCHO ! Reaction rate for SO2 and HCHO (jmm, 06/13/18)
    REAL(fp), INTENT(OUT) :: KaqHMS  ! Reaction rate for HMS and OH- (jmm, 06/13/18)
    REAL(fp), INTENT(OUT) :: KaqHMS2 ! Reaction rate for HMS and OH(aq) (jmm, 06/28/18)
    REAL(fp), INTENT(OUT) :: HSO3aq  ! Cloud bisulfite [mol/l] (qjc, 06/10/16)
    REAL(fp), INTENT(OUT) :: SO3aq   ! Cloud sulfite   [mol/l] (qjc, 06/10/16)
!
! !REMARKS:
!  Chemical Reactions:
!  ============================================================================
!  (R1) HSO3- + H2O2(aq) + H+ => SO4-- + 2H+ + H2O [Jacob, 1986]
!                                                                             .
!      d[S(VI)]/dt = k[H+][H2O2(aq)][HSO3-]/(1 + K[H+])
!      [Seinfeld and Pandis, 1998, page 366]
!                                                                             .
!  (R2) SO2(aq) + O3(aq) =>
!       HSO3-   + O3(aq) =>
!       SO3--   + O3(aq) =>
!       [Jacob, 1986; Jacobson, 1999]
!                                                                             .
!       d[S(VI)]/dt = (k0[SO2(aq)] + k1[HSO3-] + K2[SO3--])[O3(aq)]
!       [Seinfeld and Pandis, 1998, page 363]
!                                                                             .
!  (R3) HSO3-   + HCHO(aq) => HMS
!       SO3--   + HCHO(aq) => HMS + OH-
!       [Moch et al., 2018; Olson and Hoffman, 1986]
!                                                                             .
!       d[S(HMS)]/dt = (k1[HSO3-] + k2[SO3--])[HCHO(aq)]
!       [Seinfeld and Pandis, 2016, 309]
!
!  (R4) HMS + OH- => HCHO(aq) + SO3--
!       [Moch et al., 2018; Deister et al., 1986]
!        (note treated as 1st order in contrast to other reactions here)
!
!  (R5) HMS + OH(aq) =(SO2,HO2,O2)=> HCHO + 2SO4-- + O2 + 3H+ + 2H2O
!       [Jacob et al, 1986, Olson and Fessenden, 1992;
!        Seinfeld and Pandis, 2016, Table 7A.7]
!          Net reaction (R5):
!           HMS + OH(aq) =(O2)=> SO5- + HCHO + H2O
!           HO2 <=> H+ + O2-
!           SO5- + O2- =(H2O)=> HSO5- + OH- + O2
!           SO2(aq) <=> HSO3- + H+
!           H+ + OH- <=> H2O
!           HSO5- + HSO3- => 2SO4-- + 2H+
!
!  Reaction rates can be given as
!       Ra     = k [H2O2(ag)] [S(IV)]  [mole/liter*s]  OR
!       Krate  = Ra LWC R T / P        [1/s]
!                                                                             .
!  Where:
!       LWC = Liquid water content(g/m3)*10-6 [m3(water)/m3(gas)]
!       R   = 0.08205  (atm L / mol-K), Universal gas const.
!       T   = Temperature (K)
!       P   = Pressure (atm)
!                                                                             .
!  Procedure:
!  ============================================================================
!  (a ) Given [SO2] which is assumed to be total SO2 (gas+liquid) in
!        equilibrium between gas and liquid phase.
!                                                                             .
!  (b ) We can compute SO2(g) using Henry's law
!          P(so2(g)) = Xg * [SO2]
!          Xg = 1/(1 + Faq), Fraction of SO2 in gas
!       where:
!          Faq   = Kheff * R * T * LWC,
!          KHeff = Effective Henry's constant
!                                                                             .
!  (c ) Then Calculate Aquous phase, S[IV] concentrations
!        S[IV] = Kheff * P(so2(g) in atm) [M]
!                                                                             .
!  (d ) The exact same procedure is applied to calculate H2O2(aq) and HCHO(aq)
!
! !REVISION HISTORY:
!  (1 ) Updated by Rokjin Park (rjp, bmy, 12/12/02)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    REAL(fp), PARAMETER  :: R   = 0.08205e+0_fp
    REAL(fp), PARAMETER  :: dOH = 1.0e-19_fp ! [M cm^3 molec^-1]
!
! !LOCAL VARIABLES:
!
    REAL(fp)             :: KH2O2,   RA,     KS1,    KS2,    HCSO2
    REAL(fp)             :: FHCSO2,  XSO2G,  SIV,    HSO3,   XSO2AQ
    REAL(fp)             :: XHSO3,   XSO3,   KH1,    HCH2O2, FHCH2O2
    REAL(fp)             :: XH2O2G,  H2O2aq, KO0,    KO1,    KO2
    REAL(fp)             :: HCO3,    XO3g,   O3aq,   XHCHOg, HCHCHO
    REAL(fp)             :: FHCHCHO, KHCHO1, KHCHO2, KHMS,   KW1
    REAL(fp)             :: KHC1,    KHMS2

    !=================================================================
    ! AQCHEM_SO2 begins here!
    !
    ! Aqueous reaction rate
    ! HSO3- + H2O2 + H+ => SO4-- + 2H+ + H2O [Jacob, 1986]
    !=================================================================

    ! [Jacob, 1986]
    KH2O2 = 6.31e+14_fp * EXP( -4.76e+3_fp / T )

    !! [Jacobson, 1999]
    !KH2O2 = 7.45e+0_fp7 * EXP( -15.96e+0_fp * ( (298.15/T) - 1.) ) / &
    !        ( 1.e+0_fp + 13.e+0_fp * Hplus)

    !=================================================================
    ! Equilibrium reaction of SO2-H2O
    !    SO2 + H2O = SO2(aq)        (s0)
    !    SO2(ag)   = HSO3- + H+     (s1)
    !    HSO3-     = SO3-- + H+     (s2)
    !
    ! Reaction constant for Aqueous chemistry -- No big difference
    ! between Jacob and Jacobson, choose one of them.
    !
    ! Reaction rate dependent on Temperature is given
    !   H = A exp ( B (T./T - 1) )
    !
    ! For equilibrium reactions of SO2:
    !            As1      Bs1   As2      Bs2
    !  Seinfeld  1.30d-2  7.02  6.60d-8  3.76   [1998]
    !  Jacob     1.30d-2  6.75  6.31d-8  5.05   [1986]
    !  Jacobson  1.71d-2  7.04  5.99d-8  3.74   [1996]
    !=================================================================
    Ks1    = 1.30e-2_fp * EXP( 6.75e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp ) )
    Ks2    = 6.31e-8_fp * EXP( 5.05e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp ) )

    ! SIV Fraction
    XSO2aq = 1.e+0_fp/(1.e+0_fp + Ks1/Hplus + Ks1*Ks2/(Hplus*Hplus))
    XHSO3  = 1.e+0_fp/(1.e+0_fp + Hplus/Ks1 + Ks2/Hplus)
    XSO3   = 1.e+0_fp/(1.e+0_fp + Hplus/Ks2 + Hplus*Hplus/(Ks1*Ks2))

    ! Henry's constant [mol/l-atm] and Effective Henry's constant for SO2
    HCSO2  = 1.22e+0_fp * EXP( 10.55e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp) )
    FHCSO2 = HCSO2 * (1.e+0_fp + (Ks1/Hplus) + (Ks1*Ks2 / (Hplus*Hplus)))

    XSO2g  = 1.e+0_fp / ( 1.e+0_fp + ( FHCSO2 * R * T * LWC ) )
    SIV    = FHCSO2 * XSO2g * SO2 * P
    !HSO3   = Ks1 * HCSO2 * XSO2g * SO2 * P

    ! Effective HSO3aq for HOBr+HSO3
    HSO3aq = SIV * XHSO3           ! unit: M (qjc, 06/10/16)

    ! Effective SO3aq for HOBr+SO3
    SO3aq  = SIV * XSO3            ! unit: M (qjc, 06/10/16)

    !=================================================================
    ! H2O2 equilibrium reaction
    ! H2O2 + H2O = H2O2.H2O
    ! H2O2.H2O   = HO2- + H+   1)
    !
    ! Reaction rate dependent on Temperature is given
    !   H = A exp ( B (T./T - 1) )
    !
    ! For equilibrium reactions of SO2
    !            Ah1       Bh1
    !  Jacob     1.58E-12  -12.49  [1986]
    !  Jacobson  2.20E-12  -12.52  [1996]
    !=================================================================
    Kh1 = 2.20e-12_fp * EXP( -12.52e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp ) )

    ! Henry's constant [mol/l-atm] and Effective Henry's constant for H2O2
    ! [Seinfeld and Pandis, 1998]
    ! HCH2O2  = 7.45D4 * EXP( 24.48e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp) )

    ! [Jacobson,1999]
    HCH2O2  = 7.45e+4_fp * EXP( 22.21e+0_fp * (298.15e+0_fp / T - 1.e+0_fp) )
    FHCH2O2 = HCH2O2 * (1.e+0_fp + (Kh1 / Hplus))

    XH2O2g  = 1.e+0_fp / ( 1.e+0_fp + ( FHCH2O2 * R * T * LWC ) )
    !H2O2aq  = FHCH2O2 * XH2O2g * H2O2 * P

    ! Conversion rate from SO2 to SO4 via reaction with H2O2
    KaqH2O2  = kh2o2 * Ks1 * FHCH2O2 * HCSO2 * XH2O2g * XSO2g &
               * P * LWC * R * T            ! [v/v/s]

    !=================================================================
    !  Aqueous reactions of SO2 with O3
    !  SO2(aq) + O3 =>                       (0)
    !  HSO3-   + O3 => SO4-- + H+ + O2       (1)
    !  SO3--   + O3 => SO4-- + O2            (2)
    !
    ! NOTE
    ! [Jacob, 1986]
    !    KO1  = 3.49E12 * EXP( -4.83E3 / T )
    !    KO2  = 7.32E14 * EXP( -4.03E3 / T )
    !
    ! [Jacobson, 1999]
    !    KO0  = 2.40E+4
    !    KO1  = 3.70E+5 * EXP( -18.56 * ((298.15/T) - 1.))
    !    KO2  = 1.50E+9 * EXP( -17.72 * ((298.15/T) - 1.))
    !
    ! Rate constants from Jacobson is larger than those of Jacob
    ! and results in faster conversion from S(IV) to S(VI)
    ! We choose Jacob 1) 2) and Jacobson 0) here
    !=================================================================
    KO0 = 2.40e+4_fp
    KO1 = 3.49e+12_fp * EXP( -4.83e+3_fp / T )
    KO2 = 7.32e+14_fp * EXP( -4.03e+3_fp / T )

    !=================================================================
    ! H2O2 equilibrium reaction
    ! O3 + H2O = O3.H2O
    !  HCO3  = 1.13E-2 * EXP( 8.51 * (298.15/T -1.) ), S & P
    !  HCO3  = 1.13E-2 * EXP( 7.72 * (298.15/T -1.) ), Jacobson
    !=================================================================

    ! Calculate Henry's Law constant for atmospheric temperature
    HCO3  = 1.13e-2_fp * EXP( 8.51e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp ) )

    XO3g  = 1.e+0_fp / ( 1.e+0_fp + ( HCO3 * R * T * LWC ) )
    !O3aq  = HCO3 * XO3g * O3 * P

    ! Conversion rate from SO2 to SO4 via reaction with O3
    KaqO3 = (KO0*XSO2AQ + KO1*XHSO3 + KO2*XSO3) * FHCSO2 * XSO2g &
            * P * HCO3 * XO3g * LWC * R * T   ! [v/v/s]

    !(qjc, 04/10/16)
    KaqO3_1 = KO2*XSO3 * FHCSO2 * XSO2g &
              * P * HCO3 * XO3g * LWC * R * T   ! [v/v/s]

    ! ===================================================================
    ! Metal (Fe, Mn) catalyzed O2 oxidation (bec, 7/12/04)
    ! R = d[S(VI)]/dt = 750*[Mn(II)]*[S(IV)] + 2600*[Fe(III)]*[S(IV)] +
    !               1.d10*[Mn(II)]*[Fe(III)]*[S(IV)]
    ! from Seinfeld and Pandis, 1998 pg. 371
    ! S(IV) = HFCSO2 * XSO2*P*[SO2]
    ! R = KaqO2*[SO2] (v/v/s)
    ! KaqO2 = FHCSO2 * XSO2g * P *
    !        ((750*[Mn(II)])+(2600[Fe(III)])+(1.d10*[Mn(II)]*[Fe(III)]))
    ! in units of [M/s]
    ! KaqO2 = FHCSO2 * XSO2g * P *
    !        ((750*[Mn(II)])+(2600[Fe(III)])+(1.d10*[Mn(II)]*[Fe(III)])) *
    !        LWC * R * T/P
    ! in units of [v/v/s]
    ! ===================================================================

    ! Conversion rate from SO2 to SO4 via reaction with O2 (met cat)
    KaqO2 = FHCSO2 * XSO2g * ( (750e+0_fp * MnII ) + &
            ( 2600e+0_fp * FeIII ) + (1e+10_fp * MnII * FeIII ) ) * &
            LWC * R * T   ! [s-1]

    !=================================================================
    !  Aqueous reactions of SO2 with HCHO
    !     HSO3-   + HCHO(aq) => HMS + OH-           (1)
    !     SO3--   + HCHO(aq) => HMS                 (2)
    !
    !     NOTE:
    !     [Boyce and Hoffman, 1984]
    !        KHCHO1  = 7.9E2 * EXP( -16.435 * ((298.15/T) - 1.))
    !        KHCHO2  = 2.5E7 * EXP( -6.037 * ((298.15/T) - 1.))
    !
    !
    !  Aqueous reaction of HMS with OH-
    !    HMS + OH- => HCHO(aq) + SO3--             (3)
    !
    !     NOTE: unclear where B (E/R) value in Seinfeld and Pandis from,
    !     but close to Deister. Using Seinfeld and Pandis value for now
    !     [Deister et al., 1986]
    !        KHMS    = 3.6E3 * EXP( -22.027 * ((298.15/T) - 1.))
    !     [Seinfeld and Pandis, 2016; Munger et al., 1986]
    !        KHMS    = 3.6E3 * EXP( -15.09 * ((298.15/T) - 1.))
    !
    !
    !  Aqueous reaction of HMS with OH(aq)
    !    HMS + OH(aq) =(SO2,O2,HO2)=> 2SO4-- + HCHO + O2 + 3H+ + 2H2O  (4)
    !
    !    NOTE: O2, SO2, and HO2 particpate in the stoichiometry but not kinetics.
    !          Assume steady state for sulfur radicals and the following reaction chain:
    !            HMS + OH(aq) =(O2)=> SO5- + HCHO + H2O [Olsen and Fessenden, 1992]
    !            HO2 <=> H+ + O2-                       [Jacob, 1986]
    !            SO5- + O2- =(H2O)=> HSO5- + OH- + O2   [Jacob, 1986]
    !            SO2(aq) <=> HSO3- + H+
    !            H+ + OH- <=> H2O
    !            HSO5- + HSO3- => 2SO4-- + 2H+          [Jacob, 1986]
    !       Instead of assuming Henry's law for OH, use the parameter from
    !       Jacob et al, 2005 that relates gas phase OH to aqueous phase OH
    !       accounting for the HO2(aq)/O2- cylcing in cloud droplets:
    !        dOH = 1E-19 [M cm^3 molec^-1]
    !     [Olson and Fessenden, 1992]
    !        KHMS2    = 6.2E8 * EXP( -5.03 * ((298.15/T) -1.))
    !
    !
    ! (jmm, 06/28/18)
    !=================================================================
    KHCHO1 = 7.9e+2_fp * EXP( -16.44e+0_fp & ! L/mol/s
         * ( 298.15e+0_fp / T - 1.e+0_fp ) )
    KHCHO2 = 2.5e+7_fp * EXP( -6.04e+0_fp  & ! L/mol/s
         * ( 298.15e+0_fp / T - 1.e+0_fp ) )
    KHMS   = 3.6e+3_fp * EXP( -15.09e+0_fp & ! L/mol/s
         * ( 298.15e+0_fp / T - 1.e+0_fp ) )
    KHMS2  = 2.65e+8_fp * EXP( -5.03e+0_fp & ! L/mol/s
         * ( 298.15e+0_fp / T - 1.e+0_fp ) )

    !=================================================================
    ! HCHO equilibrium reaction
    ! HCHO(aq) + H2O   = HCH(OH)2
    !
    ! Reaction rate dependent on Temperature is given
    !   H = A exp ( B (T./T - 1) )
    !
    ! For equilibrium reactions of HCHO
    !                             Ah1       Bh1
    !  Sienfeld and Pandis      2.53E3    13.48  [2016]
    !
    ! (jmm, 06/15/18)
    !=================================================================
    Khc1 = 2.53e+3_fp * EXP( 13.48e+0_fp &
         * ( 298.15e+0_fp / T - 1.e+0_fp ) )

    !=================================================================
    ! H2O equilibrium reaction
    !  H2O   = H+ + OH-
    !
    ! Reaction rate dependent on Temperature is given
    !   H = A exp ( B (T./T - 1) )
    !
    ! For equilibrium reactions of HCHO
    !                             Ah1       Bh1
    !  Sienfeld and Pandis       1E-14     -22.51  [2016]
    !
    ! (jmm, 06/15/18)
    !=================================================================
    Kw1 = 1e-14_fp * EXP( -22.51e+0_fp &
         *  ( 298.15e+0_fp / T - 1.e+0_fp ) )

    ! Henry's constant [mol/l-atm] and Effective Henry's constant for HCHO
    ! [Seinfeld and Pandis, 2016]
    ! HCHCHO  = 2.5 * EXP( 21.6e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp) )
    ! (jmm, -6/15/18)
    HCHCHO  = 2.5e+0_fp * EXP( 21.6e+0_fp &
         *  (298.15e+0_fp / T - 1.e+0_fp) )
    FHCHCHO = HCHCHO * (1.e+0_fp + Khc1 )

    XHCHOg  = 1.e+0_fp / ( 1.e+0_fp + ( FHCHCHO * R * T * LWC ) )


    ! Conversion rate from SO2 to HMS via reaction with HCHO
    ! (jmm, 06/15/18)
    KaqHCHO = (KHCHO1*XHSO3 + KHCHO2*XSO3) * FHCSO2 * XSO2G &
         * P * HCHCHO * XHCHOg * LWC * R * T    ! [v/v/s]

    ! Conversion rate from HMS to SO2 via reaction with OH-
    ! (jmm, 06/15/18; MSL 1/18/22)
    KaqHMS = KHMS * ( Kw1 / Hplus ) * 0.7e0_fp  ! 70% scavenged by clouds; units [1/s]

    ! Conversion rate from HMS to SO42- & HCHO via reaction with OH(aq)
    ! (jmm, 06/28/18; MSL, 01/14/22)
    KaqHMS2 = KHMS2 * dOH * 0.7e0_fp ! 70% scavenged by clouds; units [cm3/mcl/s]

  END SUBROUTINE AQCHEM_SO2
!EOC

  SUBROUTINE SET_2R_CLD( T, LWC, FC, HPLUS, CNVFAC, P, SO2, H2O2, KaqH2O2 )

    REAL(fp), PARAMETER   :: R = 0.08205e+0_fp
    REAL(FP) :: KH2O2, KS1, KS2, T, XSO2aq, LWC, FC
    REAL(FP) :: XHSO3, XSO3, HCSO2, HPLUS, FHCSO2
    REAL(FP) :: XSO2g, KH1, HCH2O2, FHCH2O2, XH2O2g
    REAL(FP) :: KaqH2O2
    REAL(FP) :: SO2, H2O2, A, B, KAB, CNVFAC, P

    ! [Jacob, 1986]
    KH2O2  = 6.31e+14_fp * EXP( -4.76e+3_fp / T )
    Ks1    = 1.30e-2_fp * EXP( 6.75e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp ) )
    Ks2    = 6.31e-8_fp * EXP( 5.05e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp ) )

    ! SIV Fraction
    XSO2aq = 1.e+0_fp/(1.e+0_fp + Ks1/Hplus + Ks1*Ks2/(Hplus*Hplus))
    XHSO3  = 1.e+0_fp/(1.e+0_fp + Hplus/Ks1 + Ks2/Hplus)
    XSO3   = 1.e+0_fp/(1.e+0_fp + Hplus/Ks2 + Hplus*Hplus/(Ks1*Ks2))

    ! Henry's constant [mol/l-atm] and Effective Henry's constant for SO2
    HCSO2  = 1.22e+0_fp * EXP( 10.55e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp) )
    FHCSO2 = HCSO2 * (1.e+0_fp + (Ks1/Hplus) + (Ks1*Ks2 / (Hplus*Hplus)))

    XSO2g  = 1.e+0_fp / ( 1.e+0_fp + ( FHCSO2 * R * T * LWC ) )
    Kh1 = 2.20e-12_fp * EXP( -12.52e+0_fp * ( 298.15e+0_fp / T - 1.e+0_fp ) )

    ! [Jacobson,1999]
    HCH2O2  = 7.45e+4_fp * EXP( 22.21e+0_fp * (298.15e+0_fp / T - 1.e+0_fp) )
    FHCH2O2 = HCH2O2 * (1.e+0_fp + (Kh1 / Hplus))

    XH2O2g  = 1.e+0_fp / ( 1.e+0_fp + ( FHCH2O2 * R * T * LWC ) )

!    KaqH2O2  = kh2o2 * Ks1 * FHCH2O2 * HCSO2 * XH2O2g * XSO2g &

    A = SO2
    B = H2O2

    KAB = kh2o2 * Ks1 * FHCH2O2 * HCSO2 * XH2O2g * XSO2g &
               * P * LWC * R * T * CNVFAC ! cm2/mcl/s

  END SUBROUTINE SET_2R_CLD

  FUNCTION CloudHet2R( A, B, FC, KAB )  RESULT( KX )
!
! !USES
!
    USE rateLawUtilFuncs
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN) :: FC         ! Cloud Fraction [0-1]
    REAL(fp), INTENT(IN) :: A, B       ! Reactant Abundances
    REAL(fp), INTENT(IN) :: KAB        ! Bimolecular Rate Constant
!
! !RETURN VALUE:
!
    REAL(fp)             :: KX ! Grid-average loss frequency, cm3/mcl/s
!
! !REVISION HISTORY:
!  23 Aug 2018 - C. D. Holmes - Initial version
!  15 May 2021 - M. Long      - Revision for two reactants
!  See https://github.com/geoschem/geos-chem for complete history
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Residence time of air in clouds, s
    REAL(FP), PARAMETER :: TAUC = 3600.0_fp
!
! !LOCAL VARIABLES:
!
    REAL(FP) :: KAO, KBO, FF, denom, term1, term2
!EOC
!------------------------------------------------------------------------------
!

    ! Ratio of volume inside to outside cloud
    ! FF has a range [0,+inf], so cap it at 1e30
    FF = SafeDiv( FC, ( 1.0_fp - FC ), 1e30_fp )
    FF = MIN( FF, 1.0e30_fp )

    ! Avoid div by zero for the TAUC/FF term
    term1 = 0.0_fp
    IF ( ff > 0.0_fp ) term1 = tauc / ff

    ! Compute KAO and avoid div by zero
    !             term 1      term 2
    ! KAO = 1 / ( (TAUC/FF) + ( 1/(FC*KAB*B) ) ), units: 1/s
    denom = FC * KAB * B
    term2 = SafeDiv( 1.0_fp, denom, 0.0_fp )
    denom = term1 + term2
    KAO   = SafeDiv( 1.0_fp, denom, 0.0_fp )

    ! Compute KBO and avoid div by zero
    !             term 1      term 2
    ! KBO = 1 / ( (TAUC/FF) + ( 1/(FC*KAB*A) ) ), units: 1/s
    denom = FC * KAB * A
    term2 = SafeDiv( 1.0_fp, denom, 0.0_fp )
    denom = term1 + term2
    KBO   = SafeDiv( 1.0_fp, denom, 0.0_fp )

    IF ( KAO*A <= KBO*B ) THEN
       !KX = KAO/B
       KX = SafeDiv( KAO, B, 0.0_fp )
    ELSE
       !KX = KBO/A
       KX = SafeDiv( KBO, A, 0.0_fp )
    ENDIF

  END FUNCTION CloudHet2R
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: CloudHet1R
!
! !DESCRIPTION: Function CloudHet calculates the loss frequency (1/s) of gas
!  species due to heterogeneous chemistry on clouds in a partially cloudy grid
!  cell. The function uses the "entrainment limited uptake" equations of
!  Holmes et al. (2019). Both liquid and ice water clouds are treated.
!
!  For gasses that are that are consumed in multiple aqueous reactions with
!  different products, CloudHet can provide the loss frequency for each reaction
!  branch using the optional branch ratios (branchLiq, branchIce) as arguments.
!
!  Holmes, C.D., Bertram, T. H., Confer, K. L., Ronan, A. C., Wirks, C. K.,
!    Graham, K. A., Shah, V. (2019) The role of clouds in the tropospheric
!    NOx cycle: a new modeling approach for cloud chemistry and its global
!    implications, Geophys. Res. Lett. 46, 4980-4990,
!    https://doi.org/10.1029/2019GL081990
!\\
!\\
! !INTERFACE:
!
  FUNCTION CloudHet1R( fc, rate ) RESULT( kHet )
!
! !USES:
!
    USE rateLawUtilFuncs
!
! !INPUT PARAMETERS:
!
    REAL(fp), INTENT(IN) :: fc     ! Cloud Fraction [0-1]
    REAL(fp), INTENT(IN) :: rate   ! 1st order reaction rate (1/s)
!
! !RETURN VALUE:
!
    REAL(fp)             :: kHet   ! Grid-average loss frequency, 1/s
!
! !REVISION HISTORY:
!  23 Aug 2018 - C. D. Holmes - Initial version
!  25 May 2021 - M. S. Long - Modified for 1st order aqueous reaction
!                             where diffusion is not limiting
!  See https://github.com/geoschem/geos-chem for complete history
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    ! Residence time of air in clouds, s
    real(fp), parameter :: tauc = 3600.0_fp
!
! !LOCAL VARIABLES:
!
    real(fp) :: kI, gam
    real(fp) :: kk, ff, xx
!
!------------------------------------------------------------------------------
!
    ! If cloud fraction < 0.0001 (0.01%) or there is zero cloud surface area,
    ! then return zero uptake
    if ( ( fc < 0.0001_fp ) .or. ( rate <= 0.0_fp ) ) then
       kHet = 0.0_fp
       return
    endif

    !------------------------------------------------------------------------
    ! Loss frequency inside cloud
    !
    ! Assume both water and ice phases are inside the same cloud, so mass
    ! transport to both phases works in parallel (additive)
    !------------------------------------------------------------------------

    ! initialize loss, 1/s
    kI  = rate ! total loss rate of a gas in cloud

    !------------------------------------------------------------------------
    ! Grid-average loss frequency
    !
    ! EXACT expression for entrainment-limited uptake
    !------------------------------------------------------------------------

    ! Ratio (in cloud) of heterogeneous loss to detrainment, s/s
    kk = kI * tauc

    ! Ratio of volume inside to outside cloud
    ! ff has a range [0,+inf], so cap it at 1e30
    ff = SafeDiv( fc, (1e0_fp - fc), 1e30_fp )
    ff = MIN( ff, 1.0e30_fp )

    ! Ratio of mass inside to outside cloud
    ! xx has range [0,+inf], but ff is capped at 1e30, so shouldn't overflow
    xx =     ( ff        - kk        - 1.0_fp       ) / 2.0_fp +             &
         sqrt( 1.0_fp    + ff*ff     + kk*kk                   +             &
               2.0_fp*ff + 2.0_fp*kk - 2.0_fp*ff*kk ) / 2.0_fp

    ! Overall heterogeneous loss rate, grid average, 1/s
    ! kHet = kI * xx / ( 1d0 + xx )
    !  Since the expression ( xx / (1+xx) ) may behave badly when xx>>1,
    !  use the equivalent 1 / (1 + 1/x) with an upper bound on 1/x
    kHet = kI / ( 1e0_fp + SafeDiv( 1e0_fp, xx, 1e30_fp ) )

  END FUNCTION CloudHet1R
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: fullchem_InitSulfurChem
!
! !DESCRIPTION: Stores species indices in module variables for fast lookup.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE fullchem_InitSulfurChem( RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE State_Chm_Mod,  ONLY : Ind_
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  02 Jun 2000 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I, J, N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! INIT_SULFATE begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = &
 ' -> at fullchem_InitSulfurChem (in KPP/fullchem/fullchem_SulfurChemFuncs.F90'

    ! Define flags for species ID's
    id_ACTA   = Ind_( 'ACTA'   )
    id_CH2O   = Ind_( 'CH2O'   )
    id_DMS    = Ind_( 'DMS'    )
    id_DST1   = Ind_( 'DST1'   )
    id_DST2   = Ind_( 'DST2'   )
    id_DST3   = Ind_( 'DST3'   )
    id_DST4   = Ind_( 'DST4'   )
    id_H2O2   = Ind_( 'H2O2'   )
    id_HCL    = Ind_( 'HCL'    )
    id_HCOOH  = Ind_( 'HCOOH'  )
    id_HMS    = Ind_( 'HMS'    )
    id_HNO3   = Ind_( 'HNO3'   )
    id_MSA    = Ind_( 'MSA'    )
    id_NH3    = Ind_( 'NH3'    )
    id_NH4    = Ind_( 'NH4'    )
    id_NIT    = Ind_( 'NIT'    )
    id_NITs   = Ind_( 'NITs'   )
    id_O3     = Ind_( 'O3'     )
    id_OH     = Ind_( 'OH'     )
    id_pFe    = Ind_( 'pFe'    )
    id_SALA   = Ind_( 'SALA'   )  ! Sea salt aerosol     (fine mode  )
    id_SALAAL = Ind_( 'SALAAL' )  ! Sea salt alkalinity  (fine mode  )
    id_SALACL = Ind_( 'SALACL' )  ! Cl- on sea salt      (fine mode  )
    id_SALC   = Ind_( 'SALC'   )  ! Sea salt aerosol     (coarse mode)
    id_SALCAL = Ind_( 'SALCAL' )  ! Sea salt alkalinity  (coarse mode)
    id_SALCCL = Ind_( 'SALCCL' )  ! Cl- on sea salt      (coarse mode)
    id_SO2    = Ind_( 'SO2'    )
    id_SO4    = Ind_( 'SO4'    )
    id_SO4s   = Ind_( 'SO4s'   )

  END SUBROUTINE fullchem_InitSulfurChem
!EOC
END MODULE fullchem_SulfurChemFuncs
