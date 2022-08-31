!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: carboncycle_mod.F90
!
! !DESCRIPTION: Module CARBONCYCLE_MOD contains variables and routines
! for simulating CH4, CO, CO2, and OCS with an online calculation of the
! chemistry between them using KPP. It was adapted directly from
! the module CH4_CO_CO2_MOD.F provided by Beata Bukosa.
!\\
!\\
! !INTERFACE:
!
MODULE CarbonCycle_Mod
!
! !USES:
!
  USE Error_Mod,     ONLY : Safe_Div
  USE Hco_Error_Mod, ONLY : HCO_SUCCESS, HCO_FAIL, HCO_WARNING, hp
  USE PhysConstants
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Emiss_CarbonCycle
  PUBLIC :: Chem_CarbonCycle
  PUBLIC :: Init_CarbonCycle
  PUBLIC :: Cleanup_CarbonCycle
!
! !PUBLIC DATA MEMBERS:
!
  REAL(fp), ALLOCATABLE, PUBLIC :: CH4_EMIS_J(:,:,:)      ! [kg/m2/s]
  REAL(fp), PARAMETER,   PUBLIC :: XNUMOL_CH4 = AVO/16d-3 ! molec CH4 / kg CH4
!
! !REVISION HISTORY:
!  04 Apr 2022 - M.S. Long   - Initial version, based on work by B. Bukosa
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Scalars
  INTEGER               :: id_CH4,    id_CO,      id_COch4,  id_COnmvoc
  INTEGER               :: id_COisop, id_COch3oh, id_COmono, id_COacet
  INTEGER               :: id_CO2,    id_CO2ch,   id_OCS,    id_OH
  REAL(fp)              :: TROPOCH4

  ! Arrays
  REAL(fp), ALLOCATABLE :: SUMACETCO(:,:)             ! P(CO) from Acetone
  REAL(fp), ALLOCATABLE :: SUMCH3OHCO(:,:)            ! P(CO) from CH3OH
  REAL(fp), ALLOCATABLE :: SUMISOPCO(:,:)             ! P(CO) from Isoprene
  REAL(fp), ALLOCATABLE :: SUMMONOCO(:,:)             ! P(CO) from Monoterp's
  REAL(fp), ALLOCATABLE :: TCOSZ(:,:)                 ! Daily sum COS(SZA)
!
! !DEFINED PARAMETERS:
!
  INTEGER,  PARAMETER   :: N_CH4_DIAGS = 15           ! # of CH4 manual diags
  REAL(fp), PARAMETER   :: CM2PERM2    = 1.0e+4_fp    ! cm2 in 1 m2
  REAL(fp), PARAMETER   :: CM3PERM3    = 1.0e+6_fp    ! cm3 in 1 m3
!------------------------------------------------------------------------------
! NOTE: Hardwired MW's, we should use values from species_database.yml!!
  REAL(fp), PARAMETER   :: FMOL_CO2    = 44e-3_fp     ! kg CO2 / mole CO2
  REAL(fp), PARAMETER   :: FMOL_C      = 12e-3_fp     ! kg C   / mole C
  REAL(fp), PARAMETER   :: FMOL_CO     = 28e-3_fp     ! kg CO  / mole CO
  REAL(fp), PARAMETER   :: XNUMOL_CO2  = AVO/FMOL_CO2 ! molec CO2 / kg CO2
  REAL(fp), PARAMETER   :: XNUMOL_C    = AVO/FMOL_C   ! molec C   / kg C
  REAL(fp), PARAMETER   :: XNUMOL_CO   = AVO/FMOL_CO  ! molec CO  / kg CO
  REAL(fp), PARAMETER   :: XNUMOL_OH   = AVO/17e-3_fp ! molec OH  / kg OH
!------------------------------------------------------------------------------
  LOGICAL,  PARAMETER   :: ALPHA_ISOP_FROM_NOX = .FALSE.
  REAL(fp), PARAMETER   :: ALPHA_CH4   = 1.0_fp
  REAL(fp), PARAMETER   :: ALPHA_MONO  = 2e-1_fp
  REAL(fp), PARAMETER   :: ALPHA_ACET  = 2.0_fp / 3.0_fp
  INTEGER,  PARAMETER   :: NODAYS(12)  = (/ 31, 28, 31, 30, 31, 30,          &
                                            31, 31, 30, 31, 30, 31         /)

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_carboncycle
!
! !DESCRIPTION: Places emissions of CH4, CO, CO2, OCS [kg] into the
!  chemical species array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Emiss_CarbonCycle( Input_Opt,  State_Chm, State_Diag,            &
                                State_Grid, State_Met, RC                    )
!
! !USES:
!
    USE HCO_State_Mod,        ONLY : Hco_GetHcoId
    USE HCO_State_GC_Mod,     ONLY : HcoState
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_GetDiagn
    USE ErrCode_Mod
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Species_Mod,          ONLY : SpcConc
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  WARNING: Soil absorption has to be the 15th field in CH4_EMIS
!  Also: the ND58 diagnostics have now been removed.  We still need to
!  read the HEMCO manual diagnostics into CH4_EMIS for the analytical
!  inversion.  Therefore, we will keep EmissCh4 for the time-being
!  but only remove the bpch diagnostic.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                :: prtDebug
    INTEGER                :: I,       J,      L,      N
    REAL(fp)               :: AREA_M2, DTSRCE, E_CO2

    ! Strings
    CHARACTER(LEN=255)     :: ThisLoc
    CHARACTER(LEN=512)     :: ErrMsg

    ! Arrays
    REAL(fp)               :: Prod(State_Grid%NX,State_Grid%NY,State_Grid%NZ)
    REAL(hp)               :: CH4scale(N_CH4_DIAGS)
    CHARACTER(LEN=15)      :: CH4diag(N_CH4_DIAGS)

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)
    REAL(f4),      POINTER :: Ptr2D(:,:)

    !========================================================================
    ! Emiss_CarbonCycle begins here!
    !========================================================================

    ! Initialize
    RC       =  GC_SUCCESS
    prtDebug =  ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    Ptr2D    => NULL()
    Spc      => NULL()
    ErrMsg   =  ''
    ThisLoc  =  &
     ' -> at Emiss_CarbonCycle (in module GeosCore/carboncycle_mod.F90)'

    ! Exit with error if we can't find the HEMCO state object
    IF ( .NOT. ASSOCIATED( HcoState ) ) THEN
       ErrMsg = 'The HcoState object is undefined!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Emission timestep
    DTSRCE = HcoState%TS_EMIS

    !========================================================================
    ! CH4 emissions
    !
    ! --> All emission calculations are now done through HEMCO
    ! HEMCO stores emissions of all species internally in the HEMCO
    ! state object. Here, we pass these emissions into module array
    ! CH4_EMIS in units kg/m2/s. These values are then either added to
    ! the species array (full mixing scheme) or used later on in
    ! vdiff_mod.F90 if the non-local PBL mixing scheme is used.
    !
    ! The CH4_EMIS array is mostly used for backwards compatibility
    ! (especially the diagnostics). It is also used to ensure that
    ! in a multi-species simulation, species 1 (total CH4) is properly
    ! defined.
    !                                              (ckeller, 9/12/2013)
    !========================================================================
    IF ( id_CH4 > 0 ) THEN

       ! Initialize
       CH4_EMIS_J   = 0.0_fp
       CH4scale     = 1.0_hp
       CH4diag(1)   = 'CH4'
       CH4diag(2)   = 'CH4_OIL'
       CH4diag(3)   = 'CH4_GAS'
       CH4diag(4)   = 'CH4_COAL'
       CH4diag(5)   = 'CH4_LIVESTOCK'
       CH4diag(6)   = 'CH4_LANDFILLS'
       CH4diag(7)   = 'CH4_WASTEWATER'
       CH4diag(8)   = 'CH4_RICE'
       CH4diag(9)   = 'CH4_ANTHROTHER'
       CH4diag(10)  = 'CH4_BIOMASS'
       CH4diag(11)  = 'CH4_WETLAND'
       CH4diag(12)  = 'CH4_SEEPS'
       CH4diag(13)  = 'CH4_LAKES'
       CH4diag(14)  = 'CH4_TERMITES'
       CH4diag(15)  = 'CH4_SOILABSORB'    ! CH4 soilabsorb values are negative!
       CH4scale(15) = -1.0_hp             ! Need to convert to positive

       ! Loop over manual CH4 diagnostics
       DO N = 2, N_CH4_DIAGS

          ! Get a pointer to the emissions
          CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, CH4diag(N),           &
                                .FALSE.,   RC,         Ptr2D=Ptr2D          )

          ! Trap potential errors
          IF ( RC /= HCO_SUCCESS ) THEN
             ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(CH4diag(N))
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
          IF ( .not. ASSOCIATED( Ptr2D ) ) THEN
             ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(CH4diag(N))
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! Store emissions in CH4_EMIS_J [kg/m2/s]
          ! CH4scale is either -1 (for soil absorption) or 1 (everything else)
          CH4_EMIS_J(:,:,N) = Ptr2D * CH4scale(N)

          ! Free pointer for next iteration
          Ptr2D => NULL()
       ENDDO

       !---------------------------------------------------------------------
       ! Total emission: sum of all emissions - (2*soil absorption)
       ! We have to substract soil absorption twice because it is added
       ! to other emissions in the SUM function. (ccc, 7/23/09)
       !---------------------------------------------------------------------
       CH4_EMIS_J(:,:,1) = SUM( CH4_EMIS_J, 3 )                              &
                         - ( 2.0_fp * CH4_EMIS_J(:,:,15) )

       IF ( prtDebug ) THEN
          WRITE(*,*) 'CH4_EMIS (kg/m2/s):'
          WRITE(*,*) 'Total        : ', SUM( CH4_EMIS_J(:,:,1 ) )
          WRITE(*,*) 'Oil          : ', SUM( CH4_EMIS_J(:,:,2 ) )
          WRITE(*,*) 'Gas          : ', SUM( CH4_EMIS_J(:,:,3 ) )
          WRITE(*,*) 'Coal         : ', SUM( CH4_EMIS_J(:,:,4 ) )
          WRITE(*,*) 'Livestock    : ', SUM( CH4_EMIS_J(:,:,5 ) )
          WRITE(*,*) 'Landfills    : ', SUM( CH4_EMIS_J(:,:,6 ) )
          WRITE(*,*) 'Wastewater   : ', SUM( CH4_EMIS_J(:,:,7 ) )
          WRITE(*,*) 'Rice         : ', SUM( CH4_EMIS_J(:,:,8 ) )
          WRITE(*,*) 'Other anth   : ', SUM( CH4_EMIS_J(:,:,9 ) )
          WRITE(*,*) 'Biomass burn : ', SUM( CH4_EMIS_J(:,:,10) )
          WRITE(*,*) 'Wetlands     : ', SUM( CH4_EMIS_J(:,:,11) )
          WRITE(*,*) 'Seeps        : ', SUM( CH4_EMIS_J(:,:,12) )
          WRITE(*,*) 'Lakes        : ', SUM( CH4_EMIS_J(:,:,13) )
          WRITE(*,*) 'Termites     : ', SUM( CH4_EMIS_J(:,:,14) )
          WRITE(*,*) 'Soil absorb  : ', SUM( CH4_EMIS_J(:,:,15) )
       ENDIF

    ENDIF

    !========================================================================
    ! CO2 production from CO oxidation
    !========================================================================
    IF ( Input_Opt%LCHEMCO2 .and. id_CO2 > 0 ) THEN

       ! Point to chemical species array [kg/kg dry air]
       Spc => State_Chm%Species

       ! Evalulate the CO2 production from HEMCO
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'CO2_COPROD', Prod, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'CO_COPROD not found in HEMCO data list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Loop over all grid boxes
       !$OMP PARALLEL DO                                                    &
       !$OMP DEFAULT( SHARED                                               )&
       !$OMP PRIVATE( I, J, L, E_CO2, N                                    )&
       !$OMP COLLAPSE( 3                                                   )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Production is in [kg C/m3], convert to [molec/cm2/s]
          E_CO2 = Prod(I,J,L)                             & ! kg/m3
                  / CM3PERM3                              & ! => kg/cm3
                  * XNUMOL_C                              & ! => molec/cm3
                  / DTSRCE                                & ! => molec/cm3/s
                  * State_Met%BXHEIGHT(I,J,L) * 100.0_fp    ! => molec/cm2/s

          !------------------------------------------------------------------
          ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
          !
          ! Save production of CO2 from CO oxidation [kg/m2/s]
          !------------------------------------------------------------------
          IF ( State_Diag%Archive_ProdCO2fromCO ) THEN
             State_Diag%ProdCO2fromCO(I,J,L) = E_CO2       & ! molec/cm2/s
                                             / XNUMOL_CO2  & ! => kg/cm2/s
                                             * 1e4_fp        ! => kg/m2/s

          ENDIF

          ! Convert emissions from [molec/cm2/s] to [kg/kg dry air]
          ! (ewl, 9/11/15)
          E_CO2  =  E_CO2 * DTSRCE * CM2PERM2 / &
                    ( XNUMOL_CO2 * State_Met%DELP(I,J,L) &
                    * G0_100 * ( 1.0e+0_fp &
                    - State_Met%SPHU(I,J,L) * 1.0e-3_fp ) )

          ! Total CO2 [kg/kg dry air]
          Spc(id_CO2)%Conc(I,J,L) = Spc(id_CO2)%Conc(I,J,L) + E_CO2

          ! %%% TAGGED SPECIES HANDLING %%%
          ! Add chemical source of CO into the COchem species
          IF ( Input_Opt%LSPLIT .and. id_CO2ch > 0 ) THEN
             Spc(id_CO2ch)%Conc(I,J,L) = Spc(id_CO2ch)%Conc(I,J,L) + E_CO2
          ENDIF

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       ! Free pointer
       Spc => NULL()
    ENDIF

    ! Free pointers for safety's sake
    Spc   => NULL()
    Ptr2D => NULL()

  END SUBROUTINE Emiss_CarbonCycle
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_carboncycle
!
! !DESCRIPTION: Computes the chemical loss of carbon species (sources - sinks)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Chem_CarbonCycle( Input_Opt,  State_Met,  State_Chm,            &
                               State_Grid, State_Diag, RC                   )
!
! !USES:
!
    USE carboncycle_Funcs
    USE gckpp_Global
    USE gckpp_Integrator,     ONLY : INTEGRATE
    USE gckpp_Monitor,        ONLY : SPC_NAMES
    USE gckpp_Parameters
    USE gckpp_precision
    USE gckpp_Rates,          ONLY : UPDATE_RCONST
    USE ErrCode_Mod
    USE HCO_State_Mod,        ONLY : Hco_GetHcoId
    USE HCO_State_GC_Mod,     ONLY : HcoState
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE Input_Opt_Mod,        ONLY : OptInput
    USE PhysConstants,        ONLY : AVO
    USE rateLawUtilFuncs,     ONLY : SafeDiv
    USE Species_Mod,          ONLY : SpcConc
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Chm_Mod,        ONLY : ChmState, Ind_
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Met_Mod,        ONLY : MetState
    USE TIME_MOD,             ONLY : GET_TS_CHEM,     GET_TS_EMIS
    USE TIME_MOD,             ONLY : GET_MONTH,       GET_YEAR
    USE TIME_MOD,             ONLY : ITS_A_NEW_MONTH, ITS_A_NEW_YEAR
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
! !REMARKS:
!  CH4 SOURCES
!  ============================================================================
!  (1 ) Oxidation of methane, isoprene and monoterpenes (SRCO_fromHCs).
!  (2 ) Direct emissions of CO from fossil fuel combustion, biomass
!        burning and wood (for fuel) burning (SR SETEMIS).
!  (3 ) Emissions.
!                                                                             .
!  CH4 SINKS:
!  ============================================================================
!  (1 ) Removal of CO by OH (SR OHparam & CO_decay).
!  (2 ) CO uptake by soils (neglected).
!  (3 ) Transport of CO to stratosphere from troposphere
!        (in dynamical subroutines).
!  (4 ) Removal by OH (Clarissa's OH--climatol_OH.f and CO_decay.f)
!  (5 ) Transport of CH4 between troposphere and stratosphere, and
!        destruction in strat (CH4_strat.f).
!  (6 ) Removel by Cl
!
! !REVISION HISTORY:
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL                  :: failed,          found
    INTEGER                  :: HcoID,           I,            J
    INTEGER                  :: L,               NA,           N
    INTEGER                  :: NN,              MONTH,        YEAR
    INTEGER                  :: IERR
    REAL(fp)                 :: DTCHEM,          FAC_DIURNAL,  kgs_to_atomsC
    REAL(fp)                 :: Emis,            DTEMIS,       kgm3_to_mcm3OH
    REAL(fp)                 :: kgm3_to_mcm3sCO, TIN,          TOUT

    ! Strings
    CHARACTER(LEN=63)        :: DgnName
    CHARACTER(LEN=255)       :: ErrMsg
    CHARACTER(LEN=255)       :: ThisLoc

    ! Arrays
    INTEGER                  :: ICNTRL(20)
    INTEGER                  :: ISTATUS(20)
    REAL(dp)                 :: RCNTRL(20)
    REAL(dp)                 :: RSTATE(20)
    REAL(fp)                 :: BOH(State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Grid%NZ)
    REAL(fp)                 :: BCl(State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Grid%NZ)
    REAL(fp)                 :: CH4LOSS(State_Grid%NX,                       &
                                        State_Grid%NY,                       &
                                        State_Grid%NZ)
    REAL(fp)                 :: GMI_LOSS_CO(State_Grid%NX,                   &
                                            State_Grid%NY,                   &
                                            State_Grid%NZ)
    REAL(fp)                 :: GMI_PROD_CO(State_Grid%NX,                   &
                                            State_Grid%NY,                   &
                                            State_Grid%NZ)
    REAL(fp)                 :: OH(State_Grid%NX,                            &
                                   State_Grid%NY,                            &
                                   State_Grid%NZ)
    REAL(fp)                 :: PCO_CH4(State_Grid%NX,                       &
                                        State_Grid%NY,                       &
                                        State_Grid%NZ)
    REAL(fp)                 :: PCO_NMVOC(State_Grid%NX,                     &
                                          State_Grid%NY,                     &
                                          State_Grid%NZ)
    REAL(fp)                 :: PREVCH4(State_Grid%NX,                       &
                                        State_Grid%NY,                       &
                                        State_Grid%NZ)
    REAL(fp)                 :: SFC_CH4(State_Grid%NX,                       &
                                        State_Grid%NY)

    ! Pointers
    REAL(fp),      POINTER   :: AD(:,:,:)
    REAL(fp),      POINTER   :: AIRVOL(:,:,:)
    REAL(fp),      POINTER   :: T(:,:,:)
    TYPE(SpcConc), POINTER   :: Spc(:)

    !========================================================================
    ! Chem_Carboncycle begins here!
    !========================================================================

    ! Initialize
    RC      =  GC_SUCCESS
    DTCHEM  =  GET_TS_CHEM()
    Spc     => State_Chm%Species
    ErrMsg  = ''
    ThisLoc = &
     ' -> at Chem_CarbonCycle (in module GeosCore/carboncycle_mod.F90)'

    ! Exit with error if the HEMCO state object has not been initialized
    IF ( .NOT. ASSOCIATED(HcoState) ) THEN
       ErrMsg = 'HcoState object is not associated!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Zero diagnostic archival arrays to make sure that we don't have any
    ! leftover values from the last timestep near the top of the chemgrid.
    !========================================================================
    IF ( State_Diag%Archive_OHconcAfterChem ) THEN
       State_Diag%OHconcAfterChem = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_Loss ) THEN
       State_Diag%Loss = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdCOfromCH4 ) THEN
       State_Diag%ProdCOfromCH4 = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_ProdCOfromNMVOC ) THEN
       State_Diag%ProdCOfromNMVOC = 0.0_f4
    ENDIF

    !========================================================================
    ! Get fields that were archived from a fullchem simulation via HEMCO
    !========================================================================

    ! Loss frequencies of CH4 [1/s]
    DgnName = 'CH4_LOSS'
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName, CH4LOSS, RC        )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Cl concentration [mol/mol dry air]
    DgnName = 'GLOBAL_Cl'
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName, BCl, RC            )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! OH concentration [molec/cm3]
    DgnName = 'GLOBAL_OH'
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName, BOH, RC            )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! P(CO from GMI
    DgnName = 'GMI_PROD_CO'
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName, GMI_PROD_CO, RC    )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to ' // TRIM( DgnName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! L(CO) from GMI
    DgnName = 'GMI_LOSS_CO'
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName, GMI_LOSS_CO, RC    )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! P(CO) from CH4 [molec/cm3/s]
    IF ( Input_Opt%LPCO_CH4 ) THEN
       DgnName = 'PCO_CH4'
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName, PCO_CH4, RC    )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! L(CO) from CH4 [molec/cm3/s]
    IF ( Input_Opt%LPCO_NMVOC ) THEN
       DgnName = 'PCO_NMVOC'
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName, PCO_NMVOC, RC  )
       IF ( RC /= HCO_SUCCESS ) THEN
          ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Get pointer to surface CH4 data
    DgnName = 'NOAA_GMD_CH4'
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName, SFC_CH4, RC )
    IF ( RC /= HCO_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! Compute OH diurnal cycle and archive OH for diagnostics
    !========================================================================

    ! Compute OH diurnal cycle
    CALL CALC_DIURNAL( State_Grid )

    ! Diagnostic archival of OH [molec/cm3]
    IF ( State_Diag%Archive_OHconcAfterChem ) THEN

       ! NOTE: SUNCOS is declared THREADPRIVATE in gckpp_Global.F90,
       ! so we do not need to declare it PRIVATE again here.
       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( I, J, L, FAC_DIURNAL                                  )&
       !$OMP COLLAPSE( 3                                                    )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Initialize loop variables
          FAC_DIURNAL = 0.0_fp
          SUNCOS      = State_Met%SUNCOSmid(I,J)

          ! Scaling factor for OH diurnal cycles - zero at night
          IF ( SUNCOS > 0.0_fp .and. TCOSZ(I,J) > 0.0_fp ) THEN
             FAC_DIURNAL = ( SUNCOS     / TCOSZ(I,J)    ) * &
                           ( 86400.0_fp / GET_TS_CHEM() )
          ENDIF

          ! Archive OH if we are in the chemistry grid
          IF ( State_Met%InChemGrid(I,J,L) ) THEN
             IF ( State_Diag%Archive_OHconcAfterChem ) THEN
                State_Diag%OHconcAfterChem(I,J,L) = &
                 ( BOH(I,J,L) * XNUMOL_OH / CM3PERM3 * FAC_DIURNAL )
             ENDIF
          ENDIF

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    ! Save OH concentrations [molec/cm3] for metrics
    CALL CH4_OhSave_CarbonCycle( State_Met, State_Chm, State_Grid, BOH )

    !========================================================================
    ! TAGGED SPECIES HANDLING
    ! Store total CH4 and update chemically-produced CO arrays
    !========================================================================
    IF ( Input_Opt%LSPLIT ) THEN

       ! If there are multiple CH4 species, store the total CH4 concentration
       ! so that we can distribute the sink after the chemistry.
       IF ( id_CH4 > 0 ) THEN
          PrevCH4 = Spc(id_CH4)%Conc
       ENDIF

       !---------------------------------------------------------------------
       ! Update the chemically-produced tagged CO species, which are
       ! COch4, COisop, COacet, COch4, COnmvoc.
       !
       ! Because HEMCO returns 3-D emissions, we need to sum the
       ! SUMACETCO, SUMISOPCO, and SUMMONOCO arrays in the vertical.
       !
       ! Also note, skip this section if emissions are turned off,
       ! which will keep the SUM*CO arrays zeroed out.
       !---------------------------------------------------------------------
       IF ( Input_Opt%DoEmissions .and. id_CO > 0 ) THEN

          ! Conversion factor from [kg/s] --> [atoms C]
          ! (atoms C /mole C) / (kg C /mole C) * chemistry timestep [s]
          kgs_to_atomsC = ( AVO / 12e-3_fp ) * DTCHEM

          ! SUMACETCO (convert [kgC/m2/s] to [atoms C])
          HcoId = HCO_GetHcoId( 'ACET', HcoState )
          IF ( HcoId > 0 ) THEN
             SUMACETCO = SUM( HcoState%Spc(HcoID)%Emis%Val, 3 )    !kgC/m2/s
             SUMACETCO = SUMACETCO * HcoState%Grid%AREA_M2%Val     !kgC/s
             SUMACETCO = SUMACETCO * kgs_to_atomsC                 !atoms

          ELSE
             ErrMsg = 'ACET not turned on in the MEGAN!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! SUMISOPCO (convert [kgC/m2/s] to [atoms C])
          HcoId = HCO_GetHcoID( 'ISOP', HcoState )
          IF ( HcoId > 0 ) THEN
             SUMISOPCO = SUM( HcoState%Spc(HcoID)%Emis%Val, 3 )    !kgC/m2/s
             SUMISOPCO = SUMISOPCO * HcoState%Grid%AREA_M2%Val     !kgC/s
             SUMISOPCO = SUMISOPCO * kgs_to_atomsC                 !atoms

          ELSE
             ErrMsg = 'ISOP not turned on in MEGAN!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

          ! SUMMONOCO (Total monoterpene = MTPA + LIMO + MTPO)
          HcoId = HCO_GetHcoId( 'MTPA', HcoState )
          IF ( HcoId > 0 ) THEN
             ! kgC/m2/s
             SUMMONOCO = SUM( HcoState%Spc(HcoID)%Emis%Val, 3 )
          ELSE
             ErrMsg = 'MTPA not turned on in Megan_Mono !'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
          HcoId = HCO_GetHcoId( 'LIMO', HcoState )
          IF ( HcoId > 0 ) THEN
             ! kgC/m2/s
             SUMMONOCO = SUMMONOCO + SUM(HcoState%Spc(HcoID)%Emis%Val,3)
          ELSE
             ErrMsg = 'LIMO not turned on in Megan_Mono !'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF
          HcoId = HCO_GetHcoId( 'MTPO', HcoState )
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
    ENDIF

    !%%%%% TIMESTEPS %%%%%
    TIN        = 0
    TOUT       = DTCHEM

    !%%%%% Fine-tune the integrator %%%%%
    ICNTRL     =  0
    ICNTRL(1)  =  1   ! Verbose error output
    ICNTRL(2)  =  0   ! Stop model on negative values
    ICNTRL(15) = -1   ! Do not call Update_SUN, Update_RCONST w/in integrator

    ! Set a flag to denote if the chemistry failed
    failed     = .FALSE.

    !========================================================================
    ! Do chemistry -- Put everything within a large DO-loop
    ! over all grid boxes to facilitate parallelization
    !
    ! NOTE: SUNCOS is held THREADPRIVATE in gckpp_Global.F90
    !========================================================================
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J, L, N, FAC_DIURNAL                                  )&
    !$OMP COLLAPSE( 3                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 24                                             )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize PRIVATE and THREADPRIVATE variables
       C              = 0.0_dp                   ! Species conc. [molec/cm3]
       CFACTOR        = 1.0_dp                   ! Not used, set = 1
       fac_Diurnal    = 0.0_dp                   ! Diurnal scaling factor [1]
       K_STRAT        = 0.0_dp                   ! Rate in stratosphere [1/s]
       K_TROP         = 0.0_dp                   ! Rate in troposphere  [1/s]
       TROP           = 0.0_dp                   ! Toggle
       TEMP           = State_Met%T(I,J,L)       ! Temperature [K]
       INV_TEMP       = 1.0_dp / TEMP            ! 1/T  term for equations
       TEMP_OVER_K300 = TEMP / 300.0_dp          ! T/300 term for equations
       K300_OVER_TEMP = 300.0_dp / TEMP          ! 300/T term for equations
       SUNCOS         = State_Met%SUNCOSmid(I,J) ! Cos(SZA) ) [1]

       ! Convert CO, CO2, CH4 to molec/cm3 for the KPP solver
       CALL carboncycle_ConvertKgtoMolecCm3(                               &
            I          = I,                                                &
            J          = J,                                                &
            L          = L,                                                &
            id_CH4     = id_CH4,                                           &
            id_CO      = id_CO,                                            &
            id_CO2     = id_CO2,                                           &
            xnumol_CO  = xnumol_CO,                                        &
            xnumol_CH4 = xnumol_CH4,                                       &
            xnumol_CO2 = xnumol_CO2,                                       &
            State_Met  = State_Met,                                        &
            State_Chm  = State_Chm                                        )

       ! Compute the rate constants that will be used
       ! return the diurnal factor
       CALL carboncycle_ComputeRateConstants(                              &
            I            = I,                                              &
            J            = J,                                              &
            L            = L,                                              &
            LPCO_CH4     = Input_Opt%LPCO_CH4,                             &
            dtChem       = dtChem,                                         &
            bCl          = bCL(I,J,L),                                     &
            bOH          = bOH(I,J,L),                                     &
            CH4loss      = CH4loss(I,J,L),                                 &
            GMI_Prod_CO  = GMI_Prod_CO(I,J,L),                             &
            GMI_Loss_CO  = GMI_Loss_CO(I,J,L),                             &
            PCO_NMVOC    = PCO_NMVOC(I,J,L),                               &
            PCO_CH4      = PCO_CH4(I,J,L),                                 &
            TCOSZ        = TCOSZ(I,J),                                     &
            State_Met    = State_Met,                                      &
            State_Chm    = State_Chm                                      )

       !===================================================================
       ! Update reaction rates
       !===================================================================

       ! Update the array of rate constants for the KPP solver
       CALL Update_RCONST( )

       ! Call the KPP integrator
       CALL Integrate( TIN      = TIN,                                     &
                       TOUT     = TOUT,                                    &
                       ICNTRL_U = ICNTRL,                                  &
                       IERR_U   = IERR                                    )

       ! Trap potential errors
       IF ( IERR /= 1 ) failed = .TRUE.

       ! Convert CO, CO2, CH4 to molec/cm3 for the KPP solver
       CALL carboncycle_ConvertMolecCm3ToKg(                               &
            I            = I,                                              &
            J            = J,                                              &
            L            = L,                                              &
            id_CH4       = id_CH4,                                         &
            id_CO        = id_CO,                                          &
            id_COch4     = id_COch4,                                       &
            id_COnmvoc   = id_COnmvoc,                                     &
            id_CO2       = id_CO2,                                         &
            xnumol_CO    = xnumol_CO,                                      &
            xnumol_CH4   = xnumol_CH4,                                     &
            xnumol_CO2   = xnumol_CO2,                                     &
            State_Met    = State_Met,                                      &
            State_Chm    = State_Chm                                      )

       !===================================================================
       ! TAGGED SPECIES HANDLING
       ! Handle trop loss by OH for regional CO species
       !===================================================================
       IF ( Input_Opt%LSPLIT .and. id_CO > 0 ) THEN

          !%%% NOTE: Un-hardwire the species IDs!
          ! Loop over regional CO species
          DO NA = 16, State_Chm%nAdvect-11

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

             ! Update regional species
           !<<Not   SUBROUTINE re if this is correct - MSL>>
             IF (NA .ne. 16)                                     &
                  Spc(N)%Conc(I,J,L) = Spc(N)%Conc(I,J,L) *      &
                  ( 1e+0_fp - K_TROP(2) * C(ind_OH_E) * DTCHEM )

             !-----------------------------------------------------
             ! HISTORY (aka netCDF diagnostics)
             !
             ! Loss of CO by OH for "tagged" species
             !-----------------------------------------------------

             ! Units: [kg/s]
             IF ( State_Diag%Archive_Loss ) THEN
                State_Diag%Loss(I,J,L,N) = Spc(N)%Conc(I,J,L) * K_TROP(2) &
                     * C(ind_OH_E) * DTCHEM
                !   C(ind_CO2_OH) / DTCHEM  &
                !   * State_Met%AIRVOL(I,J,L) * 1e+6_fp / XNUMOL_CO
             ENDIF
          ENDDO
       ENDIF

       !=====================================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Production and loss of CO
       !
       ! NOTE: Call functions in KPP/carboncycle/carbonccycle_Funcs.F90 so
       ! that we avoid bringing in KPP species indices into this module.
       ! This avoids compile-time dependency errors.
       !=====================================================================

       ! Production of CO2 from CO oxidation [molec/cm3/s]
       IF ( Input_Opt%LCHEMCO2 ) THEN
          IF ( State_Diag%Archive_ProdCO2fromCO ) THEN
             State_Diag%ProdCO2fromCO(I,J,L) =                               &
                carboncycle_Get_CO2_OH_Flux( dtChem )
          ENDIF
       ENDIF

       ! Production of CO from CH4
       IF ( State_Diag%Archive_ProdCOfromCH4 ) THEN
          State_Diag%ProdCOfromCH4(I,J,L) =                                  &
             carboncycle_Get_CO_CH4_Flux( dtChem )
       ENDIF

       ! Units: [kg/s] Production of CO from NMVOCs
       IF ( State_Diag%Archive_ProdCOfromNMVOC ) THEN
          State_Diag%ProdCOfromNMVOC(I,J,L) =                                &
             carboncycle_Get_CO_NMVOC_Flux( dtChem )
       ENDIF

! For now, comment out tagged species handling (Bob Yantosca, 30 Aug 2022)
!       ! Loss of CO by OH -- tagged species
!       ! Units: [kg/s]
!       IF ( State_Diag%Archive_Loss ) THEN
!          State_Diag%Loss(I,J,L,16) = ( CO_OH / STTCO / DTCHEM )
!       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    IF ( failed ) THEN
       errMsg = 'KPP integration failed!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !========================================================================
    ! TAGGED SPECIES HANDLING
    ! Allocate the CH4 chemistry sink to different tagged CH4 species
    !========================================================================
    IF ( Input_Opt%LSPLIT .and. id_CH4 > 0 ) THEN
       CALL CH4_Distrib_CarbonCycle( PREVCH4,   Input_Opt,                   &
                                     State_Chm, State_Grid )
    ENDIF

  END SUBROUTINE Chem_CarbonCycle
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_ohsave_carboncycle
!
! !DESCRIPTION: Subroutine CH4\_OHSAVE archives the CH3CCl3 lifetime from the
!  OH used in the CH4 simulation. (bnd, jsw, bmy, 1/16/01, 7/20/04)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CH4_OhSave_CarbonCycle( State_Met, State_Chm, State_Grid, BOH )
!
! !USES:
!
    USE Species_Mod,        ONLY : SpcConc
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Met_Mod,      ONLY : MetState
    USE State_Grid_Mod,     ONLY : GrdState
    USE TIME_MOD,           ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    TYPE(MetState), INTENT(IN)  :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(IN)  :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN)  :: State_Grid  ! Grid State object
    REAL(fp),       INTENT(IN)  :: BOH(:,:,:)
!
! !REMARKS:
!  We now use function ITS_IN_THE_CHEMGRID from chemgrid_mod.F to diagnose
!  if box (I,J,L) is in the troposphere or stratosphere.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                :: I,       J,           L
    REAL(fp)               :: MASST,   AREA_M2
    REAL(fp)               :: KCLO,    LOSS,        OHMASS
    REAL(fp)               :: KCH4,    CH4LOSE,     CH4MASS
    REAL(fp)               :: CH4EMIS, CH4TROPMASS, BOXVL
    REAL(fp)               :: C_OH, DT, FAC_DIURNAL, SUNCOS

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    !=================================================================
    ! CH4_OHSAVE begins here!
    !=================================================================

    ! Chemistry timestep in seconds
    DT = GET_TS_CHEM()

    ! Point to chemical species array [kg/kg dry]
    Spc => State_Chm%Species

    ! Calculate OH mass and total air mass
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I,           J,       L,       BOXVL,   C_OH             )&
    !$OMP PRIVATE( OHMASS,      MASST,   KCLO,    LOSS,    KCH4             )&
    !$OMP PRIVATE( CH4TROPMASS, CH4MASS, CH4LOSE, CH4EMIS, AREA_M2          )&
    !$OMP PRIVATE( FAC_DIURNAL, SUNCOS                                      )&
    !$OMP COLLAPSE( 3                                                       )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !---------------------------------------------------------------------
       ! Initialize loop variables
       !---------------------------------------------------------------------
       BOXVL       = 0.0_fp
       SUNCOS      = 0.0_fp
       FAC_DIURNAL = 0.0_fp
       C_OH        = 0.0_fp
       OHMASS      = 0.0_fp
       MASST       = 0.0_fp
       KCLO        = 0.0_fp
       LOSS        = 0.0_fp
       CH4TROPMASS = 0.0_fp
       CH4LOSE     = 0.0_fp
       CH4EMIS     = 0.0_fp
       CH4MASS     = Spc(id_CH4)%Conc(I,J,L) * XNUMOL_CH4

       !---------------------------------------------------------------------
       ! Only process tropospheric boxes (bmy, 4/17/00)
       !---------------------------------------------------------------------
       IF ( State_Met%InChemGrid(I,J,L) ) THEN

          ! Grid box volume [cm3]
          BOXVL  = State_Met%AIRVOL(I,J,L) * 1e+6_fp

          ! Cosine of the solar zenith angle [unitless]
          SUNCOS = State_Met%SUNCOSmid(I,J)

          ! Scaling factor for diurnal cycles - zero at night
          IF ( SUNCOS > 0.0_fp .and. TCOSZ(I,J) > 0.0_fp ) THEN
             FAC_DIURNAL = ( SUNCOS     / TCOSZ(I,J) )                       &
                         * ( 86400.0_fp / DT         )
          ENDIF

          ! OH concentration in molec/cm3. BOH is imported
          ! from HEMCO in kg/m3 (ckeller, 9/16/2014)
          C_OH = BOH(I,J,L) * XNUMOL_OH / CM3PERM3 * FAC_DIURNAL

          ! Calculate OH mass [molec / box]
          OHMASS = C_OH * State_Met%AIRNUMDEN(I,J,L) * BOXVL

          ! Calculate total air mass [molec / box]
          MASST  = State_Met%AIRNUMDEN(I,J,L) * BOXVL

          ! Calculate CH3CCl3 + OH rate constant from JPL '06
          ! [cm3 / molec / s]
          KCLO = 1.64e-12_fp * EXP( -1520.0_fp / State_Met%T(I,J,L))

          ! Calculate Loss term [molec / box / s]
          LOSS   = KCLO                       * C_OH                         &
                 * State_Met%AIRNUMDEN(I,J,L) * BOXVL


          ! Calculate CH4 + OH rate constant from JPL '06
          ! [cm3 / molec / s]
          KCH4 = 2.45e-12_fp * EXP( -1775e+0_fp / State_Met%T(I,J,L) )

          ! Calculate CH4 mass [molec / box] from [kg / box]
          !CH4TROPMASS = Spc(id_CH4)%Conc(I,J,L) * XNUMOL_CH4
          !CH4MASS     = Spc(id_CH4)%Conc(I,J,L) * XNUMOL_CH4

          ! CH4MASS is now computed above (Bob Yantosca, 30 Aug 2022)
          CH4TROPMASS  = CH4MASS

          ! Calculate loss term  [molec /box / s]
          CH4LOSE = KCH4            * C_OH  * &
                    State_Met%AIRNUMDEN(I,J,L) * BOXVL

          ! Calculate CH4 emissions [molec / box / s]
          !   Only for surface level
          ! Grid box surface area [cm2]
          ! HEMCO update: CH4_EMIS now in kg/m2/s (ckeller, 9/12/2014)
          IF ( L == 1 ) THEN
             CH4EMIS = CH4_EMIS_J(I,J,1)                                     &
                     * State_Grid%Area_M2(I,J)                               &
                     * XNUMOL_CH4
          ENDIF
       ENDIF

       ! Pass OH mass, total mass, and loss to "diag_oh_mod.f",
       ! which calculates mass-weighted mean [OH] and CH3CCl3
       ! lifetime.
!         CALL DO_DIAG_OH_CH4( I, J, L, OHMASS, MASST, LOSS, &
!              CH4LOSE, CH4TROPMASS, CH4EMIS, CH4MASS )

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE CH4_OhSave_CarbonCycle
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_distrib_ccycle
!
! !DESCRIPTION: Allocates the chemistry sink to different emission species.
!  (Only called if there are tagged CH4 species.)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CH4_Distrib_CarbonCycle( PREVCH4,   Input_Opt,                  &
                                      State_Chm, State_Grid                 )
!
! !USES:
!
    USE ERROR_MOD,          ONLY : SAFE_DIV
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Species_Mod,        ONLY : SpcConc
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState

    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    REAL(fp)                      :: PREVCH4(State_Grid%NX,State_Grid%NY,State_Grid%NZ)! CH4 bef chem
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  07 Mar 2012 - M. Payer    - Added ProTeX headers
!  25 Mar 2013 - R. Yantosca - Now accept Input_Opt, State_Chm args
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER           :: I, J, L, N, NA

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    !========================================================================
    ! CH4_DISTRIB begins here
    !========================================================================

    ! Point to chemical species array [kg]
    Spc => State_Chm%Species

    !%%% NOTE: Need to unhardwire the species ID's
    ! Loop over the number of advected species
    DO NA = 2, State_Chm%nAdvect-24

       ! Advected species ID
       N = State_Chm%Map_Advect(NA)

       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( I, J, L                                               )&
       !$OMP COLLAPSE( 3                                                    )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
         Spc(N)%Conc(I,J,L) = &
              SAFE_DIV( Spc(N)%Conc(I,J,L), PREVCH4(I,J,L), 0.0_fp) &
              * Spc(id_CH4)%Conc(I,J,L)
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDDO

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE CH4_DISTRIB_CarbonCycle
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_carboncycle
!<
! !DESCRIPTION: Allocates and zeroes module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_CarbonCycle( Input_Opt,  State_Chm, State_Diag,             &
                               State_Grid, RC                                )
!
! !USES:
!
    USE gckpp_Global,    ONLY : MW, SR_MW, HENRY_CR, HENRY_K0
    USE ErrCode_Mod
    USE Input_Opt_Mod,   ONLY : OptInput
    USE State_Chm_Mod,   ONLY : ChmState, Ind_
    USE State_Diag_Mod,  ONLY : DgnState
    USE State_Grid_Mod,  ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag  ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  This routine is called from GC_INIT_EXTRA (in GeosCore/input_mod.f)
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: KppId, N

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !========================================================================
    ! Initialize
    !========================================================================
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
     ' -> at Init_CarbonCycle (in module GeosCore/ccyclechem_mod.F90)'

    !========================================================================
    ! Define GEOS-Chem species indices
    ! NOTE: Some of these are for tagged species, which are deactivated
    ! by default.  Interested users can add tagged species if they wish.
    !========================================================================
    id_CH4     = Ind_( 'CH4'     )
    id_CO      = Ind_( 'CO'      )
    id_COacet  = Ind_( 'COacet'  )
    id_COch3oh = Ind_( 'COch3oh' )
    id_COch4   = Ind_( 'COch4'   )
    id_COisop  = Ind_( 'COisop'  )
    id_COmono  = Ind_( 'COmono'  )
    id_COnmvoc = Ind_( 'COnmvoc' )
    id_CO2     = Ind_( 'CO2'     )
    id_CO2ch   = Ind_( 'CO2ch'   )
    id_OCS     = Ind_( 'OCS'     )

    !========================================================================
    ! Save physical parameters from the species_database.yml file into KPP
    ! arrays located in module gckpp_Global.F90.  These do not vary with
    ! (I,J,L) location, and so can be defined here in the init phase.
    !========================================================================
    DO KppId = 1, State_Chm%nKppSpc + State_Chm%nOmitted
       N                  = State_Chm%Map_KppSpc(KppId)
       IF ( N > 0 ) THEN
          MW(KppId)       = State_Chm%SpcData(N)%Info%MW_g
          SR_MW(KppId)    = SQRT( MW(KppId ) )
          HENRY_K0(KppId) = State_Chm%SpcData(N)%Info%Henry_K0
          HENRY_CR(KppId) = State_Chm%SpcData(N)%Info%Henry_CR
       ENDIF
    ENDDO

    !========================================================================
    ! Initialize variables for CH4 chemistry
    !========================================================================
    IF ( id_CH4 > 0 ) THEN

       ALLOCATE( CH4_EMIS_J( State_Grid%NX, State_Grid%NY, N_CH4_DIAGS ),    &
            STAT=RC )
       CALL GC_CheckVar( 'carboncycle_mod.F90:CH4_EMIS', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       CH4_EMIS_J = 0.0_fp

       ! Initialize tropoch4 (counts total decay of CH4 due to OH)
       TROPOCH4 = 0.0_fp

    ENDIF

    !========================================================================
    ! Initialize variables for CO chemistry
    !========================================================================
    IF ( id_CO > 0 ) THEN

       ! Allocate SUMISOPCO -- array for CO from isoprene
       ALLOCATE( SUMISOPCO( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'carboncycle_mod.F90:SUMISOPCO', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       SUMISOPCO = 0.0_fp

       ! Allocate SUMMONOCO -- array for CO from monoterpenes
       ALLOCATE( SUMMONOCO( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'carboncycle_mod.F90:SUMMONOCO', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       SUMMONOCO = 0.0_fp

       ! Allocate SUMCH3OH -- array for CO from CH3OH
       ALLOCATE( SUMCH3OHCO( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'carboncycle_mod.F90:SUMCH3OHCO', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       SUMCH3OHCO = 0.0_fp

       ! Allocate SUMACETCO -- array for CO from isoprene
       ALLOCATE( SUMACETCO( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'carboncycle_mod.F90:SUMACETCO', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       SUMACETCO = 0.0_fp

       ! Allocate TCOSZ -- array for sum of COS(SZA)
       ALLOCATE( TCOSZ( State_Grid%NX, State_Grid%NY ), STAT=RC )
       CALL GC_CheckVar( 'carboncycle_mod.F90:TCOSZ', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       TCOSZ = 0.0_fp

    ENDIF

    !========================================================================
    ! Initialize variables for CO2 chemistry
    !========================================================================
    ! none yet

    !========================================================================
    ! Initialize variables for OCS chemistry
    !========================================================================
    ! none yet

  END SUBROUTINE Init_CarbonCycle
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_carboncycle
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_CH4 deallocates module arrays.
!  (bmy, 1/16/01)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_CarbonCycle( RC )
!
! !USES:
!
    USE ErrCode_Mod
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC          ! Success or failure?
!EOP
!------------------------------------------------------------------------------
!BOC
!
    !=================================================================
    ! Cleanup_CarbonCycle begins here!
    !=================================================================

    ! Initialize
    RC = GC_SUCCESS

    ! Deallocate
    IF ( ALLOCATED( CH4_EMIS_J ) ) THEN
       DEALLOCATE( CH4_EMIS_J, STAT=RC )
       CALL GC_CheckVar( 'carboncycle_mod.F90:CH4_EMIS', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( SUMISOPCO ) ) THEN
       DEALLOCATE( SUMISOPCO, STAT=RC )
       CALL GC_CheckVar( 'carboncycle_mod.F90:SUMISOPCO', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( SUMMONOCO ) ) THEN
       DEALLOCATE( SUMMONOCO , STAT=RC )
       CALL GC_CheckVar( 'carboncycle_mod.F90:SUMMONOCO', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( SUMCH3OHCO ) ) THEN
       DEALLOCATE( SUMCH3OHCO, STAT=RC )
       CALL GC_CheckVar( 'carboncycle_mod.F90:SUMCH3OHCO', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( SUMACETCO ) ) THEN
       DEALLOCATE( SUMACETCO, STAT=RC )
       CALL GC_CheckVar( 'carboncycle_mod.F90:SUMACETCO', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( TCOSZ ) ) THEN
       DEALLOCATE( TCOSZ, STAT=RC )
       CALL GC_CheckVar( 'carboncycle_mod.F90:TCOSZ', 2, RC )
       RETURN
    ENDIF

  END SUBROUTINE Cleanup_CarbonCycle
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
          !$OMP PARALLEL DO                                                  &
          !$OMP DEFAULT( SHARED                                             )&
          !$OMP PRIVATE( I, J, YMID_R, TIMLOC, AHR, SUNTMP_MID              )&
          !$OMP COLLAPSE( 2                                                 )
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
END MODULE CarbonCycle_Mod
