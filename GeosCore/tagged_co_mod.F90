!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tagged_co_mod.F90
!
! !DESCRIPTION: Module TAGGED\_CO\_MOD contains variables and routines
!  used for the  geographically tagged CO simulation.
!\\
!\\
! !INTERFACE:
!
MODULE TAGGED_CO_MOD
!
! !USES:
!
  USE PhysConstants
  USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

  IMPLICIT NONE
  PRIVATE
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: CALC_DIURNAL
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: CHEM_TAGGED_CO
  PUBLIC  :: INIT_TAGGED_CO
  PUBLIC  :: CLEANUP_TAGGED_CO
!
! !REMARKS:
!  Tagged CO Species (you can modify these as needs be!)
!  ============================================================================
!  (1 ) Total CO
!  (2 ) CO from North American fossil fuel
!  (3 ) CO from European fossil fuel
!  (4 ) CO from Asian fossil fuel
!  (5 ) CO from fossil fuel from everywhere else
!  (6 ) CO from South American biomass burning
!  (7 ) CO from African biomass burning
!  (8 ) CO from Southeast Asian biomass burning
!  (9 ) CO from Oceania biomass burning
!  (10) CO from European biomass burning
!  (11) CO from North American biomass burning
!  (12) CO chemically produced from Methane
!  (13) CO from Biofuel burning (whole world)
!  (14) CO chemically produced from Isoprene
!  (15) CO chemically produced from Monoterpenes
!  (16) CO chemically produced from Methanol (CH3OH)
!  (17) CO chemically produced from Acetone
!
! !REVISION HISTORY:
!  28 Jul 2000- R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Arrays
  REAL(fp), ALLOCATABLE :: SUMACETCO  (:,:  )   ! P(CO) from Acetone
  REAL(fp), ALLOCATABLE :: SUMCH3OHCO (:,:  )   ! P(CO) from CH3OH
  REAL(fp), ALLOCATABLE :: SUMISOPCO  (:,:  )   ! P(CO) from Isoprene
  REAL(fp), ALLOCATABLE :: SUMMONOCO  (:,:  )   ! P(CO) from Monoterpenes
  REAL(fp), ALLOCATABLE :: TCOSZ      (:,:  )   ! Daily sum COS(SZA)

  ! Pointers to fields in the HEMCO data structure.
  ! These need to be declared REAL(f4), aka REAL*4.
  REAL(f4), POINTER     :: OH         (:,:,:) => NULL() ! Global OH
  REAL(f4), POINTER     :: GMI_PROD_CO(:,:,:) => NULL() ! Global P(CO)
  REAL(f4), POINTER     :: GMI_LOSS_CO(:,:,:) => NULL() ! Global L(CO)
  REAL(f4), POINTER     :: PCO_CH4    (:,:,:) => NULL() ! CH4 P(CO)
  REAL(f4), POINTER     :: PCO_NMVOC  (:,:,:) => NULL() ! NMVOC P(CO)
  REAL(f4), POINTER     :: SFC_CH4    (:,:  ) => NULL() ! Global sfc CH4

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_tagged_co
!
! !DESCRIPTION: Subroutine CHEM\_TAGGED\_CO performs CO chemistry on
!  geographically "tagged" CO species.  Loss is via reaction with OH.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEM_TAGGED_CO( Input_Opt,  State_Chm, State_Diag, State_Grid, &
                             State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD,          ONLY : CHECK_VALUE
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE HCO_Error_Mod
    USE HCO_Interface_Mod,  ONLY : HcoState, GetHcoID
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants,      ONLY : AVO
    USE State_Chm_Mod,      ONLY : ChmState, Ind_
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE TIME_MOD,           ONLY : GET_TS_CHEM,     GET_TS_EMIS
    USE TIME_MOD,           ONLY : GET_MONTH,       GET_YEAR
    USE TIME_MOD,           ONLY : ITS_A_NEW_MONTH, ITS_A_NEW_YEAR
!
! !INPUT PARAMETERS:
!
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
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  19 Oct 1999 - Q. Li, B. Duncan, B. Field - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                  :: FOUND,         LSPLIT
    LOGICAL                  :: LPCO_CH4,      LPCO_NMVOC
    INTEGER                  :: nAdvect,       HcoID,       NA
    INTEGER                  :: I,             J,           L
    INTEGER                  :: N,             MONTH,       YEAR
    REAL(fp)                 :: ALPHA_ISOP,    DTCHEM,      GCO
    REAL(fp)                 :: STTCO,         KRATE,       CH4
    REAL(fp)                 :: CO_CH4,        CO_ISOP,     CO_MONO
    REAL(fp)                 :: CO_CH3OH,      CO_OH,       CO_ACET
    REAL(fp)                 :: CO_NMVOC
    REAL(fp)                 :: CH4RATE,       DENS,        CORATE
    REAL(fp)                 :: YMID,          BOXVL,       FMOL_CO
    REAL(fp)                 :: KHI1,          KLO1,        XYRAT1
    REAL(fp)                 :: BLOG1,         FEXP1,       KHI2
    REAL(fp)                 :: KLO2,          XYRAT2,      BLOG2
    REAL(fp)                 :: FEXP2,         KCO1,        KCO2
    REAL(fp)                 :: OH_MOLEC_CM3,  FAC_DIURNAL, SUNCOS
    REAL(fp)                 :: PCO_CH4_MCM3S, PCO_NMVOC_MCM3S
    REAL(fp)                 :: kgs_to_atomsC, Emis,        DTEMIS
    REAL(fp)                 :: kgm3_to_mcm3OH,kgm3_to_mcm3sCO

    ! Strings
    CHARACTER(LEN=255)       :: ERR_VAR
    CHARACTER(LEN=255)       :: ERR_MSG
    CHARACTER(LEN=63)        :: DgnName
    CHARACTER(LEN=255)       :: ErrMsg
    CHARACTER(LEN=255)       :: ThisLoc

    ! Arrays
    INTEGER                  :: ERR_LOC(4)

    ! Pointers
    REAL(fp),        POINTER :: AD    (:,:,:  )
    REAL(fp),        POINTER :: AIRVOL(:,:,:  )
    REAL(fp),        POINTER :: Spc   (:,:,:,:)
    REAL(fp),        POINTER :: T     (:,:,:  )

    ! SAVED scalars
    LOGICAL,   SAVE          :: FIRST = .TRUE.
    REAL(fp),  SAVE          :: A3090S, A0030S, A0030N, A3090N
    INTEGER,   SAVE          :: IDch4   = -1
    INTEGER,   SAVE          :: IDnmvoc = -1
    INTEGER,   SAVE          :: IDisop  = -1
    INTEGER,   SAVE          :: IDch3oh = -1
    INTEGER,   SAVE          :: IDmono  = -1
    INTEGER,   SAVE          :: IDacet  = -1
!
! !DEFINED PARAMETERS:
!
    ! Switch to scale yield of isoprene from NOx concentration or not
    LOGICAL,   PARAMETER     :: ALPHA_ISOP_FROM_NOX = .FALSE.

    ! Yield of CO from CH4
    REAL(fp),  PARAMETER     :: ALPHA_CH4  = 1.0_fp

    ! Yield of CO from monoterpenes
    REAL(fp),  PARAMETER     :: ALPHA_MONO = 2e-1_fp

    ! Yield of CO from acetone
    REAL(fp),  PARAMETER     :: ALPHA_ACET = 2.0_fp / 3.0_fp

    !=================================================================
    ! CHEM_TAGGED_CO begins here!
    !=================================================================

    ! Assume success
    RC        = GC_SUCCESS
    ErrMsg    = ''
    ThisLoc   = ' -> at CHEM_TAGGED_CO (in GeosCore/tagged_co_mod.F90)'

    ! Number of advected species
    nAdvect   = State_Chm%nAdvect

    ! Get fields from Input_Opt
    LSPLIT     = Input_Opt%LSPLIT
    LPCO_CH4   = Input_Opt%LPCO_CH4
    LPCO_NMVOC = Input_Opt%LPCO_NMVOC

    ! DTCHEM is the number of seconds per chemistry timestep
    DTCHEM    = GET_TS_CHEM()

    ! DTEMIS is the number of seconds per emission timestep
    ! used here for conversion from HEMCO
    DTEMIS    = GET_TS_EMIS()

    ! Initialize pointers
    AD        => State_Met%AD
    AIRVOL    => State_Met%AIRVOL
    Spc       => State_Chm%Species
    T         => State_Met%T

    ! Get the molecular weight of CO [kg/mol] from the species database
    ! (All tagged species are CO, so we can use the value for species #1)
    FMOL_CO   =  State_Chm%SpcData(1)%Info%MW_g * 1.0e-3_fp

    ! Factor to convert OH from kg/m3 (from HEMCO) to molec/cm3
    kgm3_to_mcm3OH = ( AVO / 17.0e-3_fp ) * 1.0e-6_fp

    ! Factor to convert P(CO) from kg/m3 (from HEMCO) to molec/cm3/s
    ! HEMCO uses emission timestep to convert from kg/m3/s to kg/m3;
    ! apply that here.
    kgm3_to_mcm3sCO = ( AVO / (FMOL_CO * DTEMIS) ) * 1.0e-6_fp

    ! Compute diurnal cycle for OH every day (check for new day inside
    ! subroutine) - jaf 7/10/14
    CALL CALC_DIURNAL( State_Grid )

    ! Check HEMCO state object
    IF ( .NOT. ASSOCIATED(HcoState) ) THEN
       ErrMsg = 'HcoState object is not associated!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Zero diagnostic archival arrays to make sure that we don't have any
    ! leftover values from the last timestep near the top of the chemgrid
    IF ( State_Diag%Archive_Loss ) State_Diag%Loss = 0.0_f4
    IF ( State_Diag%Archive_ProdCOfromCH4 ) THEN
       State_Diag%ProdCOfromCH4 = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdCOfromNMVOC ) THEN
       State_Diag%ProdCOfromNMVOC = 0.0_f4
    ENDIF

    !=================================================================
    ! Get pointers from HEMCO for OH, P(CO), and L(CO) fields
    !
    ! NOTES:
    ! (1) We only need to get the pointers on the very first
    !      timestep.  HEMCO will update the targets automatically.
    ! (2) These calls have to be placed here instead of in routine
    !      INIT_TAGGED_CO.  This is because when INIT_TAGGED_CO is
    !      called, the HEMCO_Config file has not yet been read in.
    !=================================================================
    IF ( FIRST ) THEN

       ! Get species IDs
       IDch4    = Ind_( 'COch4'   )
       IDnmvoc  = Ind_( 'COnmvoc' )
       IDisop   = Ind_( 'COisop'  )
       IDch3oh  = Ind_( 'COch3oh' )
       IDmono   = Ind_( 'COmono'  )
       IDacet   = Ind_( 'COacet'  )

       ! Get a pointer to the OH field from the HEMCO list (bmy, 3/11/15)
       CALL HCO_GetPtr( HcoState, 'GLOBAL_OH', OH,   RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to GLOBAL_OH!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Get pointer to strat P(CO) from GMI
       CALL HCO_GetPtr( HcoState, 'GMI_PROD_CO', GMI_PROD_CO, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to GMI_PROD_CO!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Get pointer to strat L(CO) from GMI
       CALL HCO_GetPtr( HcoState, 'GMI_LOSS_CO', GMI_LOSS_CO, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to GMI_LOSS_CO!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Get pointer to trop P(CO) from CH4 if needed
       IF ( LPCO_CH4 ) THEN
          CALL HCO_GetPtr( HcoState, 'PCO_CH4', PCO_CH4, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Cannot get pointer to PCO_CH4!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       IF ( LPCO_NMVOC ) THEN
          CALL HCO_GetPtr( HcoState, 'PCO_NMVOC', PCO_NMVOC, RC )
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Cannot get pointer to PCO_NMVOC!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
       ENDIF

       ! Get pointer to surface CH4 data
       CALL HCO_GetPtr( HcoState, 'NOAA_GMD_CH4', SFC_CH4, RC )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to NOAA_GMD_CH4!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    !=================================================================
    ! Read emissions from HEMCO into SUMACETCO, SUMISOPCO, SUMMONOCO
    ! arrays.  These are needed to update the chemically-produced
    ! tagged CO species: COch4, COisop, COacet, etc. (bmy, 6/1/16)
    !
    ! NOTE: HEMCO returns 3-D emissions, so we will sum in the
    ! vertical because the code expects SUMACETCO, SUMISOPCO, and
    ! SUMMONOCO to be 2-D arrays.  We may change that later.
    ! (bmy, 6/1/16)
    !
    ! Also note, skip this section if emissions are turned off,
    ! which will keep the SUM*CO arrays zeroed out. (bmy, 10/11/16)
    !=================================================================
    IF ( Input_Opt%LEMIS ) THEN

       ! Conversion factor from [kg/s] --> [atoms C]
       ! (atoms C /mole C) / (kg C /mole C) * chemistry timestep [s]
       kgs_to_atomsC = ( AVO / 12e-3_fp ) * DTCHEM

       ! SUMACETCO (convert [kgC/m2/s] to [atoms C])
       HcoId = GetHcoId( 'ACET' )
       IF ( HcoId > 0 ) THEN
          SUMACETCO = SUM( HcoState%Spc(HcoID)%Emis%Val, 3 )    ! kgC/m2/s
          SUMACETCO = SUMACETCO * HcoState%Grid%AREA_M2%Val     ! kgC/s
          SUMACETCO = SUMACETCO * kgs_to_atomsC                 ! atoms C
       ELSE
          ErrMsg = 'ACET not turned on in the MEGAN!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! SUMISOPCO (convert [kgC/m2/s] to [atoms C])
       HcoId = GetHcoId( 'ISOP' )
       IF ( HcoId > 0 ) THEN
          SUMISOPCO = SUM( HcoState%Spc(HcoID)%Emis%Val, 3 )    ! kgC/m2/s
          SUMISOPCO = SUMISOPCO * HcoState%Grid%AREA_M2%Val     ! kgC/s
          SUMISOPCO = SUMISOPCO * kgs_to_atomsC                 ! atoms C
       ELSE
          ErrMsg = 'ISOP not turned on in MEGAN!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! SUMMONOCO (Total monoterpene = MTPA + LIMO + MTPO)
       HcoId = GetHcoId( 'MTPA' )
       IF ( HcoId > 0 ) THEN
          ! kgC/m2/s
          SUMMONOCO = SUM( HcoState%Spc(HcoID)%Emis%Val, 3 )
       ELSE
          ErrMsg = 'MTPA not turned on in Megan_Mono !'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       HcoId = GetHcoId( 'LIMO' )
       IF ( HcoId > 0 ) THEN
          ! kgC/m2/s
          SUMMONOCO = SUMMONOCO + SUM(HcoState%Spc(HcoID)%Emis%Val,3)
       ELSE
          ErrMsg = 'LIMO not turned on in Megan_Mono !'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       HcoId = GetHcoId( 'MTPO' )
       IF ( HcoId > 0 ) THEN
          ! kgC/m2/s
          SUMMONOCO = SUMMONOCO + SUM(HcoState%Spc(HcoID)%Emis%Val,3)
       ELSE
          ErrMsg = 'MTPO not turned on in Megan_Mono !'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
       ! SUMMONOCO (convert [kgC/m2/s] to [atoms C])
       SUMMONOCO = SUMMONOCO * HcoState%Grid%AREA_M2%Val ! kgC/s
       SUMMONOCO = SUMMONOCO * kgs_to_atomsC ! atoms C

    ENDIF

    !=================================================================
    ! Do tagged CO chemistry -- Put everything within a large
    ! DO-loop over all grid boxes to facilitate parallelization
    !=================================================================
    !$OMP PARALLEL DO &
    !$OMP DEFAULT( SHARED                                                  ) &
    !$OMP PRIVATE( I,       J,            L,           N,          STTCO   ) &
    !$OMP PRIVATE( GCO,     DENS,         CH4RATE,     CO_CH4,     CH4     ) &
    !$OMP PRIVATE( KRATE,   CO_ISOP,      CO_CH3OH,    CO_MONO,    CO_ACET ) &
    !$OMP PRIVATE( CORATE,  CO_OH,        YMID,        ALPHA_ISOP, KHI1    ) &
    !$OMP PRIVATE( KLO1,    XYRAT1,       BLOG1,       FEXP1,      KHI2    ) &
    !$OMP PRIVATE( KLO2,    XYRAT2,       BLOG2,       FEXP2,      KCO1    ) &
    !$OMP PRIVATE( KCO2,    OH_MOLEC_CM3, FAC_DIURNAL, SUNCOS,     ERR_LOC ) &
    !$OMP PRIVATE( PCO_CH4_MCM3S,         PCO_NMVOC_MCM3S,         CO_NMVOC) &
    !$OMP PRIVATE( ERR_VAR, ERR_MSG,      BOXVL                            )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Latitude of grid box
       YMID     = State_Grid%YMid(I,J)

       ! Grid box volume [cm3]
       BOXVL    = State_Met%AIRVOL(I,J,L) * 1e+6_fp

       !==============================================================
       ! (0) Define useful quantities
       !==============================================================

       ! STTCO [molec CO/cm3/kg CO] converts [kg CO] --> [molec CO/cm3]
       ! kg CO/box * box/cm3 * mole/0.028 kg CO * Avog.#/mole
       STTCO = 1.0_fp  / AIRVOL(I,J,L) / 1e+6_fp / FMOL_CO * AVO

       ! GCO is CO concentration in [molec CO/cm3]
       GCO   = Spc(I,J,L,1) * STTCO

       ! DENS is the number density of air [molec air/cm3]
       DENS  = AD(I,J,L) * 1000.e+0_fp / BOXVL * AVO / AIRMW

       ! Cosine of the solar zenith angle [unitless]
       SUNCOS = State_Met%SUNCOSmid(I,J)

       ! Scaling factor for diurnal cycles - zero at night
       IF ( SUNCOS > 0.0_fp .and. TCOSZ(I,J) > 0.0_fp ) THEN
          FAC_DIURNAL = ( SUNCOS    / TCOSZ(I,J)    ) * &
                        ( 86400.0_fp / GET_TS_CHEM() )
       ELSE
          FAC_DIURNAL = 0.0_fp
       ENDIF

       ! Now impose a diurnal cycle on OH.
       ! This is done in other offline simulations but was
       ! missing from tagged CO (jaf, 3/12/14)
       !
       ! NOTE: HEMCO brings in OH in kg/m3, so we need to also
       ! apply a conversion to molec/cm3 here. (bmy, 10/12/16)
       OH_MOLEC_CM3 = ( OH(I,J,L) * kgm3_to_mcm3OH ) * FAC_DIURNAL

       ! Make sure OH is not negative
       OH_MOLEC_CM3 = MAX( OH_MOLEC_CM3, 0e+0_fp )

       ! Also impose diurnal cycle on P(CO) from CH4, NMVOC
       IF ( LPCO_CH4 ) THEN

          ! HEMCO brings in PCO in kg/m3, so we need to also
          ! apply a conversion to molec/cm3/s
          PCO_CH4_MCM3S = PCO_CH4(I,J,L) * kgm3_to_mcm3sCO * FAC_DIURNAL

          ! Make sure PCO_CH4 is not negative
          PCO_CH4_MCM3S = MAX( PCO_CH4_MCM3S, 0e+0_fp )

       ENDIF

       IF ( LPCO_NMVOC ) THEN

          ! HEMCO brings in PCO in kg/m3/s, so we need to also
          ! apply a conversion to molec/cm3/s
          PCO_NMVOC_MCM3S = PCO_NMVOC(I,J,L) * kgm3_to_mcm3sCO * FAC_DIURNAL

          ! Make sure PCO_NMVOC is not negative
          PCO_NMVOC_MCM3S = MAX( PCO_NMVOC_MCM3S, 0e+0_fp )

       ENDIF

       !==============================================================
       ! (1a) Production of CO by reaction with CH4
       !==============================================================

       ! Initialize
       CO_CH4 = 0e+0_fp

       ! Test level for stratosphere or troposphere
       IF ( State_Met%InStratosphere(I,J,L) ) THEN

          !===========================================================
          ! (1a-1) Production of CO from CH4 in the stratosphere
          !===========================================================

          ! Call GET_PCO_LCO_STRAT to get the P(CO) rate from CH4
          CH4RATE = GMI_PROD_CO(I,J,L)

          ! Convert units of CH4RATE from [v/v/s] to [molec CO/cm3]
          CO_CH4  = CH4RATE * DTCHEM * DENS

       ELSE

          !===========================================================
          ! (1a-2) Production of CO from CH4 in the troposphere
          !===========================================================

          ! CH4 concentration from HEMCO [ppbv]
          CH4 = SFC_CH4(I,J)

          ! Convert CH4 from [ppbv] to [molec CH4/cm3]
          CH4 = CH4 * 1e-9_fp * DENS

          ! Use rates saved from full chemistry
          IF ( LPCO_CH4 ) THEN

             ! Diurnal cycle & unit conversion applied above
             ! Multiply by DTCHEM to convert to [molec CO/cm3]
             CO_CH4 = PCO_CH4_MCM3S * DTCHEM

          ! Original behaviour based on latitude bands
          ELSE

             ! CH4 concentration [ppbv] for the given latitude band
             ! (bmy, 1/2/01)
             CH4 = A3090S
             IF ( YMID >= -30.0 .and. YMID < 0.0  ) CH4 = A0030S
             IF ( YMID >=   0.0 .and. YMID < 30.0 ) CH4 = A0030N
             IF ( YMID >=  30.0                   ) CH4 = A3090N

             ! Convert CH4 from [ppbv] to [molec CH4/cm3]
             CH4 = CH4 * 1e-9_fp * DENS

             ! Calculate updated rate constant [s-1] (bnd, bmy, 1/2/01)
             KRATE = 2.45e-12_fp * EXP( -1775.e+0_fp / T(I,J,L) )

             ! Production of CO from CH4 = alpha * k * [CH4] * [OH] * dt
             ! Units are [molec CO/cm3]
             CO_CH4 = ALPHA_CH4 * KRATE * CH4 * OH_MOLEC_CM3 * DTCHEM

          ENDIF

       ENDIF

       ! Check CO_CH4 for NaN or Infinity
       ERR_LOC = (/ I, J, L, 0 /)
       ERR_VAR = 'CO_CH4'
       ERR_MSG = 'STOP at tagged_co_mod:1'
       CALL CHECK_VALUE( CO_CH4, ERR_LOC, ERR_VAR, ERR_MSG )

       ! Use rates saved from full chemistry
       IF ( LPCO_NMVOC) THEN
          !===========================================================
          ! (1b) Production of CO from NMVOCs (all but CH4)
          !===========================================================

          ! Initialize
          CO_NMVOC = 0e+0_fp

          ! CO production only happens in the troposphere. However,
          ! this is already taken into account in the input files,
          ! which are filled with zeros above the flexchem levels.

          ! Diurnal cycle & unit conversion applied above
          ! Multiply by DTCHEM to convert to [molec CO/cm3]
          CO_NMVOC = PCO_NMVOC_MCM3S * DTCHEM

          ! Make sure it is not negative
          CO_NMVOC = MAX( CO_NMVOC, 0e+0_fp )

          ! Check CO_NMVOC for NaN or Infinity
          ERR_LOC = (/ I, J, L, 0 /)
          ERR_VAR = 'CO_NMVOC'
          ERR_MSG = 'STOP at tagged_co_mod:1'
          CALL CHECK_VALUE( CO_NMVOC, ERR_LOC, ERR_VAR, ERR_MSG )

          ! Set individual VOC contributions to zero to avoid
          ! diagnostic problems
          CO_ISOP  = 0e+0_fp
          CO_CH3OH = 0e+0_fp
          CO_ACET  = 0e+0_fp
          CO_MONO  = 0e+0_fp

       ! Otherwise use original behaviour
       ELSE
          !===========================================================
          ! (1b) Production of CO from ISOPRENE and METHANOL (CH3OH)
          !===========================================================

          ! Initialize
          CO_ISOP  = 0e+0_fp
          CO_CH3OH = 0e+0_fp

          ! Isoprene is emitted only into the surface layer
          IF ( L == 1 ) THEN

             !========================================================
             ! Yield of CO from ISOP: 30%, from Miyoshi et al., 1994.
             ! They estimate globally 105 Tg C/yr of CO is produced
             ! from isoprene oxidation.
             !--------------------------------------------------------
             ! We need to scale the Isoprene flux to get the CH3OH
             ! (methanol) flux.  Currently, the annual isoprene flux in
             ! GEOS-CHEM is ~ 397 Tg C.
             !
             ! Daniel Jacob recommends a flux of 100 Tg/yr CO from CH3OH
             ! oxidation based on Singh et al. 2000 [JGR 105, 3795-3805]
             ! who estimate a global methanol source of 122 Tg yr-1, of
             ! which most (75 Tg yr-1) is "primary biogenic".  He also
             ! recommends for now that the CO flux from CH3OH oxidation
             ! be scaled to monthly mean isoprene flux.
             !
             ! To get CO from METHANOL oxidation, we must therefore
             ! multiply the ISOPRENE flux by the following scale factor:
             !  ( 100 Tg CO / 397 Tg C ) * ( 12 g C/mole / 28 g CO/mole )
             !-----------------------------------------------------------
             ! We now call GET_ALPHA_ISOP to get the yield factor of
             ! CO produced from isoprene, as a function of NOx, or
             ! as a constant. (bnd, bmy, 6/14/01)
             !=======================================================

             ! Get CO yield from ISOPRENE
             ! Use a 30% yield from Miyoshi et al., 1994.
             ! They estimate globally 105 Tg C/yr of CO is produced
             ! from isoprene oxidation.
             ! ALPHA_ISOP = (0.3 molec CO/atoms C) x (5 atoms C/molec ISOP)
             ALPHA_ISOP = 1.5e+0_fp

             ! P(CO) from Isoprene Flux = ALPHA_ISOP * Flux(ISOP)
             ! Convert from [molec ISOP/box] to [molec CO/cm3]
             !
             ! Units of SUMISOPCO are [atoms C/box/time step].
             ! Division by 5 is necessary to convert to
             ! [molec ISOP/box/timestep].
             !
             ! Units of ALPHA_ISOP are [molec CO/molec ISOP]
             ! Units of CO_ISOP are [molec CO/cm3]
             CO_ISOP  = SUMISOPCO(I,J)   / BOXVL &
                      / 5.0_fp           * ALPHA_ISOP

             ! P(CO) from CH3OH is scaled to Isoprene Flux (see above)
             ! Units are [molec CO/cm3]
             CO_CH3OH = ( SUMISOPCO(I,J) / BOXVL    ) &
                      * ( 100.0_fp       / 397.0_fp ) &
                      * ( 12.0_fp        / 28.0_fp  )

             ! Zero SUMISOPCO and SUMCH3OHCO for the next emission step
             SUMISOPCO(I,J)  = 0e+0_fp
             SUMCH3OHCO(I,J) = 0e+0_fp

             ! Check CO_ISOP for NaN or Infinity
             ERR_LOC = (/ I, J, L, 0 /)
             ERR_VAR = 'CO_ISOP'
             ERR_MSG = 'STOP at tagged_co_mod:2'
             CALL CHECK_VALUE( CO_ISOP,  ERR_LOC, ERR_VAR, ERR_MSG )

             ! Check CO_CH3OH for NaN or Infinity
             ERR_VAR = 'CO_CH3OH'
             ERR_MSG = 'STOP at tagged_co_mod:3'
             CALL CHECK_VALUE( CO_CH3OH, ERR_LOC, ERR_VAR, ERR_MSG )
          ENDIF

          !===========================================================
          ! (1c) Production of CO from MONOTERPENE oxidation
          !===========================================================

          ! Initialize
          CO_MONO = 0.e+0_fp

          ! Monoterpenes are emitted only into the surface layer
          IF ( L == 1 ) THEN

             !=======================================================
             ! Assume the production of CO from monoterpenes is
             ! instantaneous even though the lifetime of intermediate
             ! species may be on the order of hours or days.  This
             ! assumption will likely cause CO from monoterpene
             ! oxidation to be too high in the box in which the
             ! monoterpene is emitted.
             !-------------------------------------------------------
             ! The CO yield here is taken from:
             !   Hatakeyama et al. JGR, Vol. 96, p. 947-958 (1991)
             !   Vinckier et al. Fresenius Env. Bull., Vol. 7, p.361-368
             !     (1998)
             !
             ! Hatakeyama:  "The ultimate yield of CO from the
             !   tropospheric oxidation of terpenes (including both O3
             !   and OH reactions) was estimated to be 20% on the carbon
             !   number basis."  They studied ALPHA- & BETA-pinene.
             !
             ! Vinckier  :  "R(CO)=1.8+/-0.3" : 1.8/10 is about 20%.
             !--------------------------------------------------------
             ! Calculate source of CO per time step from monoterpene
             ! flux (assume lifetime very short) using the C number basis:
             !
             !   CO [molec CO/cm3] = Flux [atoms C from MONO/box] /
             !                       Grid Box Volume [cm^-3]       *
             !                       ALPHA_MONO
             !
             ! where ALPHA_MONO = 0.2 as explained above.
             !========================================================

             ! P(CO) from Monoterpene Flux =  alpha * Flux(Mono)
             ! Units are [molec CO/cm3]
             CO_MONO = ( SUMMONOCO(I,J) / BOXVL ) * ALPHA_MONO

             ! Zero SUMMONOCO for the next emission step
             SUMMONOCO(I,J) = 0e+0_fp

             ! Check CO_MONO for NaN or Infinity
             ERR_LOC = (/ I, J, L, 0 /)
             ERR_VAR = 'CO_MONO'
             ERR_MSG = 'STOP at tagged_co_mod:4'
             CALL CHECK_VALUE( CO_MONO, ERR_LOC, ERR_VAR, ERR_MSG )
          ENDIF

          !===========================================================
          ! (1d) Production of CO from oxidation of ACETONE
          !
          ! ALPHA_ACET = 2/3 to get a yield for CO.  This accounts
          ! for acetone loss from reaction with OH And photolysis.
          ! The acetone sources taken into account are:
          !
          ! (a) Primary emissions of acetone from biogenic sources
          ! (b) Secondary production of acetone from monoterpene
          !      oxidation
          ! (c) Secondary production of acetone from ALK4 and
          !      propane oxidation
          ! (d) Direct emissions of acetone from biomass burning and
          !      fossil fuels
          ! (e) direct emissions from ocean
          !
          ! Calculate source of CO per time step from biogenic acetone
          ! # molec CO/cc = ALPHA * ACET Emission Rate * dt
          !===========================================================

          ! Initialize
          CO_ACET = 0.e+0_fp

          ! Biogenic acetone sources are emitted only into the surface layer
          IF ( L == 1 ) THEN

             ! Units are [molec CO/cc]
             CO_ACET = SUMACETCO(I,J) / BOXVL * ALPHA_ACET

             ! Zero SUMACETCO for the next emission step
             SUMACETCO(I,J) = 0e+0_fp

             ! Check CO_ACET for NaN or Infinity
             ERR_LOC = (/ I, J, L, 0 /)
             ERR_VAR = 'CO_ACET'
             ERR_MSG = 'STOP at tagged_co_mod:5'
             CALL CHECK_VALUE( CO_ACET, ERR_LOC, ERR_VAR, ERR_MSG )
          ENDIF

          ! Add individual NMVOC contributions together to get total
          ! NMVOC contribution
          CO_NMVOC = CO_ISOP + CO_CH3OH + CO_MONO + CO_ACET

       ENDIF !Saved rates vs surface fluxes

       !==============================================================
       ! (1e) Add production of CO into the following tagged species:
       !
       ! (a) Species #12: CO produced from CH4
       ! (a) Species #13: CO produced from NMVOC
       ! (b) Species #14: CO produced from ISOPRENE - old only
       ! (c) Species #15: CO produced from MONOTERPENES - old only
       ! (d) Species #16: CO produced from METHANOL (CH3OH) - old only
       ! (e) Species #17: CO produced from ACETONE - old only
       !
       ! %%% NOTE: If you are modifying the tagged CO simulation,
       ! %%% and your simulation has less than 12 species, then
       ! %%% then comment out this section.  If you don't you can
       ! %%% get an array-out-of-bounds error (bmy, 6/11/08)
       !==============================================================
       IF ( LSPLIT ) THEN
          IF ( IDch4 >= 0 ) &
               Spc(I,J,L,IDch4) = Spc(I,J,L,IDch4) + CO_CH4   / STTCO
          IF ( IDnmvoc >= 0 ) &
               Spc(I,J,L,IDnmvoc) = Spc(I,J,L,IDnmvoc) + CO_NMVOC /STTCO
          IF (.not. LPCO_NMVOC) THEN
             IF ( IDisop >= 0 ) &
                  Spc(I,J,L,IDisop) = Spc(I,J,L,IDisop) + CO_ISOP /STTCO
             IF ( IDmono >= 0 ) &
                  Spc(I,J,L,IDmono) = Spc(I,J,L,IDmono) + CO_MONO /STTCO
             IF ( IDch3oh >= 0 ) &
                  Spc(I,J,L,IDch3oh) = Spc(I,J,L,IDch3oh)+CO_CH3OH/STTCO
             IF ( IDacet >= 0 ) &
                  Spc(I,J,L,IDacet) = Spc(I,J,L,IDacet) + CO_ACET /STTCO
          ENDIF
       ENDIF

       !==============================================================
       ! (2a) Loss of CO due to chemical reaction w/ OH
       !==============================================================

       ! Select out tropospheric or stratospheric boxes
       IF ( State_Met%InStratosphere(I,J,L) ) THEN

          !===========================================================
          ! (2a-1) Stratospheric loss of CO due to chemical rxn w/ OH
          !===========================================================

          ! Get the L(CO) rate in the stratosphere in [s-1]
          CORATE  = GMI_LOSS_CO(I,J,L)

          ! CO_OH is the fraction of CO lost to OH [unitless]
          CO_OH   = CORATE * DTCHEM

          ! Check CO_OH for NaN or Infinity
          ERR_LOC = (/ I, J, L, 0 /)
          ERR_VAR = 'CO_OH'
          ERR_MSG = 'STOP at tagged_co_mod:6'
          CALL CHECK_VALUE( CO_OH, ERR_LOC, ERR_VAR, ERR_MSG )

          ! Handle strat loss by OH for regional CO species
          IF ( LSPLIT ) THEN

             ! Loop over regional CO species
             DO NA = 2, nAdvect

                ! Advected species ID
                N = State_Chm%Map_Advect(NA)

                !------------------------------------------------------
                ! NOTE: The proper order should be:
                !   (1) Calculate CO loss rate
                !   (2) Update AD65 array
                !   (3) Update the SPC array using the loss rate
                !
                ! Therefore, we have now moved the computation of the
                ! ND65 diagnostic before we apply the loss to the
                ! tagged CO concentrations stored in the SPC array.
                !
                !    -- Jenny Fisher (27 Mar 2017)
                !
                !------------------------------------------------------

                !-----------------------------------------------------
                ! HISTORY (aka netCDF diagnostics)
                !
                ! Loss of CO by OH for "tagged" species
                !-----------------------------------------------------

                !Units: [kg/s]
                IF ( State_Diag%Archive_Loss ) THEN
                   State_Diag%Loss(I,J,L,N) = ( CORATE * Spc(I,J,L,N) )
                ENDIF

                ! Loss
                Spc(I,J,L,N) = Spc(I,J,L,N) * ( 1.0_fp - CO_OH )

                ! Species shouldn't be less than zero
                IF ( Spc(I,J,L,N) < 0.0_fp ) THEN
                   Spc(I,J,L,N) = 0.0_fp
                ENDIF

                ! Error check
                ERR_LOC = (/ I, J, L, N /)
                ERR_VAR = 'Spc (points to State_Chm%Species)'
                ERR_MSG = 'STOP at tagged_co_mod:7'
                CALL CHECK_VALUE( Spc(I,J,L,N), ERR_LOC, ERR_VAR, ERR_MSG )

             ENDDO
          ENDIF

          ! CO_OH above is just the fraction of CO lost by OH.  Here
          ! we multiply it by GCO (the initial value of Spc in molec/cm3)
          ! to convert it to an amount of CO lost by OH [molec/cm3]
          ! (bmy, 2/19/02)
          CO_OH = GCO * CO_OH

       ELSE

          !===========================================================
          ! (2a-2) Tropospheric loss of CO due to chemical rxn w/ OH
          !
          !  DECAY RATE
          !  The decay rate (KRATE) is calculated by:
          !
          !     No change from JPL '97 to JPL '03 (jaf, 2/27/09)
          !     OH + CO -> products (JPL '03)
          !     k = (1 + 0.6Patm) * 1.5E-13
          !
          !  KRATE has units of [ molec^2 CO / cm6 / s ]^-1,
          !  since this is a 2-body reaction.
          !
          ! Updated rate constant from JPL 2006; now more complicated.
          ! From JPL 2006: "The  reaction between HO and CO to yield
          ! H + CO2 akes place on a potential energy surface that
          ! contains the radical HOCO.  The yield of H and CO2 is
          ! diminished as the pressure rises.  The loss of reactants
          ! is thus the sum of two processes, an association to yield
          ! HOCO and the chemical activation process yielding H and
          ! CO2." So we now need two complicated reactions. The code
          ! is more or less copied from calcrate.f as implemented by
          ! jmao (jaf, 3/4/09)
          !
          ! Update rate constant for JPL 15-10 (mps, 4/24/17):
          ! GY( A0 = 5.9e-33, B0 = 1.4e0,  A1 = 1.1e-12, B1 = -1.3e0,
          !     A2 = 1.5e-13, B2 = -0.6e0, A3 = 2.1e09,  B3 = -6.1e0 )
          !
          ! is now
          !
          ! GY( A0 = 5.9e-33, B0 = 1.,     A1 = 1.1e-12, B1 = -1.3e0,
          !     A2 = 1.5e-13, B2 = 0.,     A3 = 2.1e09,  B3 = -6.1e0 )
          !===========================================================

          ! Decay rate
          ! NOTE: This code is nearly identical to function GC_OHCO
          !  found in KPP/Standard/gckpp_Rates.F90)
          ! new JPL 2006 version (jaf, 3/4/09)
          ! KLO1 = k_0(T) from JPL Data Eval (page 2-1)
          KLO1   = 5.9e-33_fp * ( 300 / T(I,J,L) )**(1.e+0_fp)
          ! KHI1 = k_inf(T) from JPL Data Eval (page 2-1)
          KHI1   = 1.1e-12_fp * ( 300 / T(I,J,L) )**(-1.3e+0_fp)
          XYRAT1 = KLO1 * DENS / KHI1
          BLOG1  = LOG10(XYRAT1)
          FEXP1  = 1.e+0_fp / ( 1.e+0_fp + BLOG1 * BLOG1 )
          ! KCO1 = k_f([M],T) from JPL Data Eval (page 2-1)
          KCO1   = KLO1 * DENS * 0.6**FEXP1 / ( 1.e+0_fp + XYRAT1 )
          ! KLO2 = k_0(T) from JPL Data Eval (page 2-1)
          KLO2   = 1.5e-13_fp * ( 300 / T(I,J,L) )**(0.e+0_fp)
          ! KHI2 = k_inf(T) from JPL Data Eval (page 2-1)
          KHI2   = 2.1e+09_fp * ( 300 / T(I,J,L) )**(-6.1e+0_fp)
          XYRAT2 = KLO2 * DENS / KHI2
          BLOG2  = LOG10(XYRAT2)
          FEXP2  = 1.e+0_fp / ( 1.e+0_fp + BLOG2 * BLOG2 )
          ! KCO2 = k_f^ca([M],T) from JPL Data Eval (page 2-2)
          KCO2   = KLO2 * 0.6**FEXP2 / ( 1.e+0_fp + XYRAT2 )

          ! KRATE is the sum of the two.
          KRATE  = KCO1 + KCO2

          ! CO_OH = Tropospheric loss of CO by OH [molec/cm3]
          ! Now use OH_MOLEC_CM3, which includes a diurnal cycle.
          CO_OH = KRATE * GCO * OH_MOLEC_CM3 * DTCHEM

          ! Handle trop loss by OH for regional CO species
          IF ( LSPLIT ) THEN

             ! Loop over regional CO species
             DO NA = 2, nAdvect

                ! Advected species ID
                N = State_Chm%Map_Advect(NA)

                !-----------------------------------------------------
                ! NOTE: The proper order should be:
                !   (1) Calculate CO loss rate
                !   (2) Update AD65 array
                !   (3) Update the SPC array using the loss rate
                !
                ! Therefore, we have now moved the computation of the
                ! ND65 diagnostic before we apply the loss to the
                ! tagged CO concentrations stored in the SPC array.
                !
                !    -- Jenny Fisher (27 Mar 2017)
                !-----------------------------------------------------

                !-----------------------------------------------------
                ! HISTORY (aka netCDF diagnostics)
                !
                ! Loss of CO by OH for "tagged" species
                !-----------------------------------------------------

                ! Units: [kg/s]
                IF ( State_Diag%Archive_Loss ) THEN
                   State_Diag%Loss(I,J,L,N) = ( KRATE &
                                              *   OH_MOLEC_CM3 &
                                              *   Spc(I,J,L,N) )
                ENDIF

                ! Use tropospheric rate constant
                Spc(I,J,L,N) = Spc(I,J,L,N) * &
                               ( 1e+0_fp - KRATE * OH_MOLEC_CM3 * DTCHEM )

                ! Error check
                ERR_LOC = (/ I, J, L, N /)
                ERR_VAR = 'Spc (points to State_Chm%Species)'
                ERR_MSG = 'STOP at tagged_co_mod:8'
                CALL CHECK_VALUE( Spc(I,J,L,N), ERR_LOC, ERR_VAR, ERR_MSG )
             ENDDO
          ENDIF

       ENDIF

       !==============================================================
       ! Save the total chemical production from various sources
       ! into the total CO species Spc(I,J,L,1)
       !==============================================================

       ! GCO is the total CO before chemistry was applied [molec CO/cm3]
       ! Add to GCO the sources and sinks listed above
       GCO = GCO + CO_CH4 + CO_NMVOC - CO_OH

       ! Convert net CO from [molec CO/cm3] to [kg] and store in
       ! State_Chm%Species
       Spc(I,J,L,1) = GCO / STTCO

       !==============================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Production of CO species
       !==============================================================

       ! Units: [kg/s] Production of CO from CH4
       IF ( State_Diag%Archive_ProdCOfromCH4 ) THEN
          State_Diag%ProdCOfromCH4(I,J,L) = CO_CH4 / STTCO / DTCHEM
       ENDIF

       ! Units: [kg/s] Production of CO from NMVOCs
       IF ( State_Diag%Archive_ProdCOfromNMVOC ) THEN
          State_Diag%ProdCOfromNMVOC(I,J,L) = CO_NMVOC / STTCO / DTCHEM
       ENDIF

       !==============================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Loss of total CO species
       !==============================================================

       ! Units: [kg/s]
       IF ( State_Diag%Archive_Loss ) THEN
          State_Diag%Loss(I,J,L,1) = ( CO_OH / STTCO / DTCHEM )
       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointers
    AD       => NULL()
    AIRVOL   => NULL()
    Spc      => NULL()
    T        => NULL()

  END SUBROUTINE CHEM_TAGGED_CO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: calc_diurnal
!
! !DESCRIPTION: Subroutine CALC\_DIRUNAL computes the sume of the cosine
!  of the solar zenith angle over a 24 hour day as well as the total
!  length of daylight to scale the offline OH concentrations.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CALC_DIURNAL( State_Grid )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
    USE TIME_MOD, ONLY : ITS_A_NEW_DAY
    USE TIME_MOD, ONLY : GET_MINUTE,    GET_SECOND,      GET_HOUR
    USE TIME_MOD, ONLY : GET_TS_CHEM,   GET_DAY_OF_YEAR, GET_LOCALTIME
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN) :: State_Grid  ! Grid State object
!
! !REVISION HISTORY:
!  12 Mar 2014 - J. Fisher - Initial version, Copied from OHNO3TIME in
!                            carbon_mod and COSSZA in dao_mod
!  See https://github.com/geoschem/geos-chem for complete history

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE       :: FIRST = .TRUE.
    INTEGER             :: I, J, N, NDYSTEP
    INTEGER             :: SECOND,  MINUTE, TS_SUN
    REAL*8              :: GMT_MID, TIMLOC, FACTOR
    REAL*8              :: R,       AHR,    DEC
    REAL*8              :: YMID_R,  SUNTMP_MID
!
! !DEFINED PARAMETERS:
!
    ! Coefficients for solar declination angle
    REAL*8,  PARAMETER :: A0 = 0.006918d0
    REAL*8,  PARAMETER :: A1 = 0.399912d0
    REAL*8,  PARAMETER :: A2 = 0.006758d0
    REAL*8,  PARAMETER :: A3 = 0.002697d0
    REAL*8,  PARAMETER :: B1 = 0.070257d0
    REAL*8,  PARAMETER :: B2 = 0.000907d0
    REAL*8,  PARAMETER :: B3 = 0.000148d0

    !=================================================================
    ! CALC_DIURNAL begins here!
    !=================================================================

    ! Only do at the start of a new day
    IF ( FIRST .or. ITS_A_NEW_DAY() ) THEN

       ! Zero array
       TCOSZ = 0d0

       ! Get time for central chemistry timestep
       TS_SUN = GET_TS_CHEM()                     ! Chemistry interval
       SECOND = GET_SECOND()                      ! Current seconds
       MINUTE = GET_MINUTE()                      ! Current minutes
       FACTOR = ( MINUTE * 60 + SECOND ) / TS_SUN ! Multiplying factor

       ! GMT at the midpoint of the chemistry time interval for first
       ! timestep of the day
       GMT_MID  = ( DBLE( GET_HOUR()        )        ) &
                + ( DBLE( TS_SUN * FACTOR ) / 3600d0 ) &
                + ( DBLE( TS_SUN / 2      ) / 3600d0 )

       ! Solar declination angle (low precision formula):
       ! Path length of earth's orbit traversed since Jan 1 [radians]
       R = ( 2d0 * PI / 365d0 ) * FLOAT( GET_DAY_OF_YEAR() - 1 )
       DEC = A0 - A1*COS(    R) + B1*SIN(    R) &
                - A2*COS(2d0*R) + B2*SIN(2d0*R) &
                - A3*COS(3d0*R) + B3*SIN(3d0*R)

       ! NDYSTEP is # of chemistry time steps
       NDYSTEP = INT( 24d0 * 3600d0 / GET_TS_CHEM() )

       ! Loop forward through NDYSTEP "fake" timesteps for this day
       DO N = 1, NDYSTEP

          ! Increment GMT (hours) to midpoint of next timestep
          IF ( N > 1 ) GMT_MID = GMT_MID + TS_SUN / 3600d0

          ! Loop over surface grid boxes
          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( I, J, YMID_R, TIMLOC, AHR, SUNTMP_MID )
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Zero SUNTMP_MID
             SUNTMP_MID = 0d0

             ! Grid box latitude center [radians]
             YMID_R = State_Grid%YMid_R(I,J)

             ! Local time at box (I,J) [hours]
             TIMLOC = GET_LOCALTIME( I, J, 1, State_Grid, GMT=GMT_MID)

             ! Hour angle at box (I,J) [radians]
             AHR = ABS( TIMLOC - 12d0 ) * 15d00 * PI_180

             !===========================================================
             ! The cosine of the solar zenith angle (SZA) is given by:
             !
             !  cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)
             !
             ! where LAT = the latitude angle,
             !       DEC = the solar declination angle,
             !       AHR = the hour angle, all in radians.
             !
             ! If SUNCOS < 0, then the sun is below the horizon, and
             ! therefore does not contribute to any solar heating.
             !===========================================================

             ! Compute Cos(SZA)
             SUNTMP_MID = sin(YMID_R) * sin(DEC) + &
                          cos(YMID_R) * cos(DEC) * cos(AHR)

             ! TCOSZ is the sum of SUNTMP_MID at location (I,J)
             ! Do not include negative values of SUNTMP_MID
             TCOSZ(I,J) = TCOSZ(I,J) + MAX( SUNTMP_MID, 0d0 )

          ENDDO
          ENDDO
          !$OMP END PARALLEL DO
       ENDDO

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

  END SUBROUTINE CALC_DIURNAL
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tagged_co
!
! !DESCRIPTION: Subrogutine INIT\_TAGGED\_CO allocates memory to module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_TAGGED_CO( Input_Opt, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(DgnState), INTENT(IN)    :: State_Diag  ! Diagnostic state object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid state object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  19 Jul 2000 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! INIT_TAGGED_CO begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Init_Tagged_CO (in module GeosCore/tagged_co_mod.F90)'

    ! Allocate SUMISOPCO -- array for CO from isoprene
    ALLOCATE( SUMISOPCO( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'tagged_co_mod.F90:SUMISOPCO', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SUMISOPCO = 0.0_fp

    ! Allocate SUMMONOCO -- array for CO from monoterpenes
    ALLOCATE( SUMMONOCO( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'tagged_co_mod.F90:SUMMONOCO', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SUMMONOCO = 0.0_fp

    ! Allocate SUMCH3OH -- array for CO from CH3OH
    ALLOCATE( SUMCH3OHCO( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'tagged_co_mod.F90:SUMCH3OHCO', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SUMCH3OHCO = 0.0_fp

    ! Allocate SUMACETCO -- array for CO from isoprene
    ALLOCATE( SUMACETCO( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'tagged_co_mod.F90:SUMACETCO', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    SUMACETCO = 0.0_fp

    ! Allocate TCOSZ -- array for sum of COS(SZA)
    ALLOCATE( TCOSZ( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'tagged_co_mod.F90:TCOSZ', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    TCOSZ = 0.0_fp

  END SUBROUTINE INIT_TAGGED_CO
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_tagged_co
!
! !DESCRIPTION: Subroutine CLEANUP\_TAGGED\_CO deallocates memory from
!  previously allocated module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CLEANUP_TAGGED_CO( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  19 Jul 2000 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
    ! Assume success
    RC = GC_SUCCESS

    ! Deallocate variables
    IF ( ALLOCATED( SUMISOPCO ) ) THEN
       DEALLOCATE( SUMISOPCO, STAT=RC )
       CALL GC_CheckVar( 'tagged_co_mod.F90:SUMISOPCO', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( SUMMONOCO ) ) THEN
       DEALLOCATE( SUMMONOCO , STAT=RC )
       CALL GC_CheckVar( 'tagged_co_mod.F90:SUMMONOCO', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( SUMCH3OHCO ) ) THEN
       DEALLOCATE( SUMCH3OHCO, STAT=RC )
       CALL GC_CheckVar( 'tagged_co_mod.F90:SUMCH3OHCO', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( SUMACETCO ) ) THEN
       DEALLOCATE( SUMACETCO, STAT=RC )
       CALL GC_CheckVar( 'tagged_co_mod.F90:SUMACETCO', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( TCOSZ ) ) THEN
       DEALLOCATE( TCOSZ, STAT=RC )
       CALL GC_CheckVar( 'tagged_co_mod.F90:TCOSZ', 2, RC )
       RETURN
    ENDIF

    ! Free pointers
    GMI_PROD_CO => NULL()
    GMI_LOSS_CO => NULL()
    OH          => NULL()
    SFC_CH4     => NULL()

  END SUBROUTINE CLEANUP_TAGGED_CO
!EOC
END MODULE TAGGED_CO_MOD
