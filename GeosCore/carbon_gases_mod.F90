!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: carbon_gases_mod.F90
!
! !DESCRIPTION: Module CARBON_GASES_MOD contains variables and routines
! for simulating CH4, CO, CO2, and OCS with an online calculation of the
! chemistry between them using KPP. It was adapted directly from
! the module CH4_CO_CO2_MOD.F provided by Beata Bukosa.
!\\
!\\
! !INTERFACE:
!
MODULE Carbon_Gases_Mod
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
  PUBLIC :: Emiss_Carbon_Gases
  PUBLIC :: Chem_Carbon_Gases
  PUBLIC :: Init_Carbon_Gases
  PUBLIC :: Cleanup_Carbon_Gases
!
! !PUBLIC DATA MEMBERS:
!
  REAL(fp), ALLOCATABLE, PUBLIC :: CH4_EMIS_J(:,:,:)      ! [kg/m2/s]
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
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ TAGGED SPECIES HANDLING
!@@@ Set a switch to activate tagged species
!@@@ NOTE: By default, this is not supported and needs further work!
!#define ACTIVATE_TAGGED_SPECIES
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ! Scalars
  LOGICAL               :: useGlobOH,  useGlobOHv5
  INTEGER               :: id_CH4,     id_CO,      id_COch4,   id_COnmvoc
  INTEGER               :: id_COisop,  id_COch3oh, id_COmono,  id_COacet
  INTEGER               :: id_CO2,     id_CO2ch,   id_OCS,     id_OH
  REAL(fp)              :: xnumol_CH4, xnumol_CO,  xnumol_CO2, xnumol_OH

  ! Arrays
  REAL(fp), ALLOCATABLE :: sumOfCosSza(:,:)
!
! !DEFINED PARAMETERS:
!
  INTEGER,  PARAMETER   :: N_CH4_DIAGS = 15
  REAL(fp), PARAMETER   :: CM2perM2    = 1.0e+4_fp
  REAL(fp), PARAMETER   :: CM3perM3    = 1.0e+6_fp
  REAL(fp), PARAMETER   :: toMolecCm3  = ( AVO / AIRMW ) * 1.0e-3_fp

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emiss_carbon_gases
!
! !DESCRIPTION: Places emissions of CH4, CO, CO2, OCS [kg] into the
!  chemical species array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Emiss_Carbon_Gases( Input_Opt,  State_Chm, State_Diag,            &
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
    INTEGER                :: I,        J
    INTEGER                :: L,        N
    REAL(fp)               :: dtSrce,   E_CO2

    ! Strings
    CHARACTER(LEN=255)     :: thisLoc
    CHARACTER(LEN=512)     :: errMsg

    ! Arrays
    REAL(hp)               :: CH4scale(N_CH4_DIAGS)
    REAL(fp)               :: PCO2_fr_CO(                                    &
                               State_Grid%NX,                                &
                               State_Grid%NY,                                &
                               State_Grid%NZ)

    ! String arrays
    CHARACTER(LEN=15)      :: CH4diag(N_CH4_DIAGS)

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)
    REAL(f4),      POINTER :: Ptr2D(:,:)
!
! !DEFINED PARAMETERS:
!
    ! Molecules C per kg C
    REAL(fp),    PARAMETER :: xnumol_C = AVO / 12.0e-3_fp

    !========================================================================
    ! Emiss_Carbon_Gases begins here!
    !========================================================================

    ! Initialize
    RC       =  GC_SUCCESS
    prtDebug =  ( Input_Opt%LPRT .and. Input_Opt%amIRoot )
    Ptr2D    => NULL()
    Spc      => NULL()
    errMsg   =  ''
    thisLoc  =  &
     ' -> at Emiss_Carbon_Gases (in module GeosCore/carbon_gases_mod.F90)'

    ! Exit with error if we can't find the HEMCO state object
    IF ( .NOT. ASSOCIATED( HcoState ) ) THEN
       errMsg = 'The HcoState object is undefined!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Emission timestep
    dtSrce = HcoState%TS_EMIS

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
             errMsg = 'Cannot get pointer to HEMCO field ' // TRIM(CH4diag(N))
             CALL GC_Error( errMsg, RC, thisLoc )
             RETURN
          ENDIF
          IF ( .not. ASSOCIATED( Ptr2D ) ) THEN
             errMsg = 'Cannot get pointer to HEMCO field ' // TRIM(CH4diag(N))
             CALL GC_Error( errMsg, RC, thisLoc )
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
       CALL HCO_GC_EvalFld( Input_Opt,  State_Grid, 'CO2_COPROD',            &
                            PCO2_fr_CO, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'CO_COPROD not found in HEMCO data list!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Loop over all grid boxes
       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( I, J, L, E_CO2, N                                     )&
       !$OMP COLLAPSE( 3                                                    )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Production is in [kg C/m3], convert to [molec/cm2/s]
          E_CO2 = PCO2_fr_CO(I,J,L)                       & ! kg/m3
                  / CM3perM3                              & ! => kg/cm3
                  * xnumol_C                              & ! => molec/cm3
                  / dtSrce                                & ! => molec/cm3/s
                  * State_Met%BXHEIGHT(I,J,L) * 100.0_fp    ! => molec/cm2/s

          !------------------------------------------------------------------
          ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
          !
          ! Save production of CO2 from CO oxidation [kg/m2/s]
          !------------------------------------------------------------------
          IF ( State_Diag%Archive_ProdCO2fromCO ) THEN
             State_Diag%ProdCO2fromCO(I,J,L) = E_CO2       & ! molec/cm2/s
                                             / xnumol_CO2  & ! => kg/cm2/s
                                             * CM2perM2      ! => kg/m2/s

          ENDIF

          ! Convert emissions from [molec/cm2/s] to [kg/kg dry air]
          ! (ewl, 9/11/15)
          E_CO2  =  E_CO2 * DTSRCE * CM2perM2 /                              &
                    ( XNUMOL_CO2 * State_Met%DELP(I,J,L)                     &
                    * G0_100 * ( 1.0e+0_fp                                   &
                    - State_Met%SPHU(I,J,L) * 1.0e-3_fp )                   )

          ! Total CO2 [kg/kg dry air]
          Spc(id_CO2)%Conc(I,J,L) = Spc(id_CO2)%Conc(I,J,L) + E_CO2

#ifdef ACTIVATE_TAGGED_SPECIES
          !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          !@@@ TAGGED SPECIES HANDLING
          !@@@ Add chemical source of CO into the COchem species
          !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
          IF ( Input_Opt%LSPLIT .and. id_CO2ch > 0 ) THEN
             Spc(id_CO2ch)%Conc(I,J,L) = Spc(id_CO2ch)%Conc(I,J,L) + E_CO2
          ENDIF
#endif

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

  END SUBROUTINE Emiss_Carbon_Gases
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chem_carbon_gases
!
! !DESCRIPTION: Computes the chemical loss of carbon species (sources - sinks)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Chem_Carbon_Gases( Input_Opt,  State_Met,  State_Chm,            &
                               State_Grid, State_Diag, RC                   )
!
! !USES:
!
    USE carbon_Funcs
    USE gckpp_Global
    USE gckpp_Integrator,     ONLY : Integrate
   !USE gckpp_Monitor,        ONLY : Spc_Names
    USE gckpp_Parameters
    USE gckpp_Precision
    USE gckpp_Rates,          ONLY : Update_Rconst
    USE ErrCode_Mod
    USE HCO_State_Mod,        ONLY : Hco_GetHcoId
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_HcoStateOK
    USE Input_Opt_Mod,        ONLY : OptInput
    USE rateLawUtilFuncs,     ONLY : SafeDiv
    USE Species_Mod,          ONLY : SpcConc
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Chm_Mod,        ONLY : ChmState, Ind_
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Met_Mod,        ONLY : MetState
    USE Time_Mod,             ONLY : Get_Ts_Chem
    USE Timers_Mod
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
    ! SAVEd scalars
    LOGICAL                :: first = .TRUE.

    ! Scalars
    LOGICAL                :: failed
    INTEGER                :: HcoID,    I
    INTEGER                :: J,        L
    INTEGER                :: NA,       N
    INTEGER                :: IERR,     threadNum
    REAL(fp)               :: dtChem,   facDiurnal
    REAL(fp)               :: tsPerDay

    ! Strings
    CHARACTER(LEN=255)     :: errMsg
    CHARACTER(LEN=255)     :: thisLoc

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    ! Arrays
    INTEGER                :: ICNTRL(20)
    INTEGER                :: ISTATUS(20)
    REAL(dp)               :: RCNTRL(20)
    REAL(dp)               :: RSTATE(20)
    REAL(fp)               :: OHdiurnalFac(State_Grid%NX, State_Grid%NY)

    ! Arrays for data read in via HEMCO
    REAL(fp) :: Global_OH(   State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    REAL(fp) :: Global_Cl(   State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    REAL(fp) :: LCH4_by_OH(  State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    REAL(fp) :: LCO_in_Strat(State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    REAL(fp) :: PCO_in_Strat(State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    REAL(fp) :: PCO_fr_CH4  (State_Grid%NX, State_Grid%NY, State_Grid%NZ)
    REAL(fp) :: PCO_fr_NMVOC(State_Grid%NX, State_Grid%NY, State_Grid%NZ)
#ifdef ACTIVATE_TAGGED_SPECIES
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !@@@ TAGGED SPECIES HANDLING
    !@@@ PrevCH4 stores CH4 before chemistry, for later distribution
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    REAL(fp) :: PREVCH4(     State_Grid%NX, State_Grid%NY, State_Grid%NZ)
#endif


#ifdef MODEL_CLASSIC
#ifndef NO_OMP
    ! For GEOS-Chem Classic timers
    INTEGER,     EXTERNAL  :: OMP_GET_THREAD_NUM
#endif
#endif

    !========================================================================
    ! Chem_Carbon_Gases begins here!
    !========================================================================

    ! Initialize
    RC       =  GC_SUCCESS
    dtChem   =  Get_Ts_Chem()
    tsPerDay =  86400.0_fp / dtChem
    Spc      => State_Chm%Species
    errMsg   = ''
    thisLoc  = &
     ' -> at Chem_Carbon_Gases (in module GeosCore/carbon_gases_mod.F90)'

    ! First-time only safety checks
    IF ( first ) THEN

       ! Exit if the HEMCO state object has not yet been initialized
       IF ( .not. HCO_GC_HcoStateOK() ) THEN
          errMsg = 'HcoState object is not associated!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Determine which OH oxidant field we are using
       CALL InquireGlobalOHversion( Input_Opt, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          errMsg = 'Error encountered in "InquireGlobalOHversion"!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Reset first-time flag
       first = .FALSE.
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
    ! Read chemical inputs (oxidant fields, concentrations) via HEMCO
    !========================================================================
    CALL ReadChemInputFields(                                                &
         Input_Opt    = Input_Opt,                                           &
         State_Grid   = State_Grid,                                          &
         State_Met    = State_Met,                                           &
         Global_OH    = Global_OH,                                           &
         Global_Cl    = Global_Cl,                                           &
         LCH4_by_OH   = LCH4_by_OH,                                          &
         LCO_in_Strat = LCO_in_Strat,                                        &
         PCO_in_Strat = PCO_in_Strat,                                        &
         PCO_fr_CH4   = PCO_fr_CH4,                                          &
         PCO_fr_NMVOC = PCO_fr_NMVOC,                                        &
         RC           = RC                                                  )

    !========================================================================
    ! Compute OH diurnal cycle scaling factor
    ! (this scales OH by the position of the sun, and zeroes it at night)
    !========================================================================
    CALL Calc_Diurnal(                                                       &
         State_Grid   = State_Grid,                                          &
         State_Met    = State_Met,                                           &
         OHdiurnalFac = OHdiurnalFac                                        )

    !========================================================================
    ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
    !
    ! Diagnostic archival of OH [molec/cm3]
    !========================================================================
    IF ( State_Diag%Archive_OHconcAfterChem ) THEN

       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( I, J, L                                               )&
       !$OMP COLLAPSE( 3                                                    )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          ! Archive OH if we are in the chemistry grid [molec/cm3]
          IF ( State_Met%InChemGrid(I,J,L) ) THEN
             IF ( State_Diag%Archive_OHconcAfterChem ) THEN
                State_Diag%OHconcAfterChem(I,J,L) = Global_OH(I,J,L)         &
                                                  * OHdiurnalFac(I,J)
             ENDIF
          ENDIF

       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

#ifdef ACTIVATE_TAGGED_SPECIES
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%% TAGGED SPECIES HANDLING
    !%%% If there are multiple CH4 species, store the total CH4 conc
    !%%% so that we can distribute the sink after the chemistry.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IF ( Input_Opt%LSPLIT .and. id_CH4 > 0 ) THEN
       PrevCH4 = Spc(id_CH4)%Conc
    ENDIF
#endif

    !========================================================================
    ! Main chemistry loop -- call KPP to integrate the mechanism forward
    !========================================================================

    ! KPP forward-Euler integrator settings
    ICNTRL     =  0
    ICNTRL(1)  =  1   ! Verbose error output
    ICNTRL(2)  =  0   ! Stop model on negative values
    ICNTRL(15) = -1   ! Do not call Update_SUN, Update_RCONST w/in integrator

    ! Set a flag to denote if the chemistry failed
    failed     = .FALSE.

    ! Loop over grid boxes
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J, L, N                                               )&
#ifdef MODEL_CLASSIC
#ifndef NO_OMP
    !$OMP PRIVATE( threadNum                                                )&
#endif
#endif
    !$OMP COLLAPSE( 3                                                       )&
    !$OMP SCHEDULE( DYNAMIC, 24                                             )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize PRIVATE and THREADPRIVATE loop variables
       C              = 0.0_dp                    ! Species conc. [molec/cm3]
       CFACTOR        = 1.0_dp                    ! Not used, set = 1
       k_Strat        = 0.0_dp                    ! Rate in stratosphere [1/s]
       k_Trop         = 0.0_dp                    ! Rate in troposphere  [1/s]
       NUMDEN         = State_Met%AIRNUMDEN(I,J,L)! Air density [molec/cm3]
       TROP           = 0.0_dp                    ! Toggle for reaction
       TEMP           = State_Met%T(I,J,L)        ! Temperature [K]
       INV_TEMP       = 1.0_dp / TEMP             ! 1/T  term for equations
       TEMP_OVER_K300 = TEMP / 300.0_dp           ! T/300 term for equations
       K300_OVER_TEMP = 300.0_dp / TEMP           ! 300/T term for equations
       SUNCOS         = State_Met%SUNCOSmid(I,J)  ! Cos(SZA) ) [1]
#ifdef MODEL_CLASSIC
#ifndef NO_OMP
       threadNum      = OMP_GET_THREAD_NUM() + 1   ! OpenMP thread number
#endif
#endif

       !=====================================================================
       ! Start KPP main timer
       !=====================================================================
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start(                                                  &
               timerName = "  -> KPP",                                       &
               inLoop    = .TRUE.,                                           &
               threadNum = ThreadNum,                                        &
               RC        = RC                                               )
       ENDIF

       !=====================================================================
       ! Convert CO, CO2, CH4 to molec/cm3 for the KPP solver
       !=====================================================================

       ! Convert units
       CALL carbon_ConvertKgtoMolecCm3(                                      &
            I          = I,                                                  &
            J          = J,                                                  &
            L          = L,                                                  &
            id_CH4     = id_CH4,                                             &
            id_CO      = id_CO,                                              &
            id_CO2     = id_CO2,                                             &
            xnumol_CH4 = xnumol_CH4,                                         &
            xnumol_CO  = xnumol_CO,                                          &
            xnumol_CO2 = xnumol_CO2,                                         &
            State_Met  = State_Met,                                          &
            State_Chm  = State_Chm                                          )

       !===================================================================
       ! Update reaction rates
       !===================================================================

       ! Start timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start(                                                  &
               timerName = "     RCONST",                                    &
               inLoop    = .TRUE.,                                           &
               threadNum = ThreadNum,                                        &
               RC        =  RC                                              )
       ENDIF

       ! Compute the rate constants that will be used
       CALL carbon_ComputeRateConstants(                                     &
            I                = I,                                            &
            J                = J,                                            &
            L                = L,                                            &
            dtChem           = dtChem,                                       &
            ConcClMnd        = Global_Cl(I,J,L),                             &
            ConcOHmnd        = Global_OH(I,J,L),                             &
            LCH4_by_OH       = LCH4_by_OH(I,J,L),                            &
            LCO_in_Strat     = LCO_in_Strat(I,J,L),                          &
            OHdiurnalFac     = OHdiurnalFac(I,J),                            &
            PCO_in_Strat     = PCO_in_Strat(I,J,L),                          &
            PCO_fr_CH4_use   = Input_Opt%LPCO_CH4,                           &
            PCO_fr_CH4       = PCO_fr_CH4(I,J,L),                            &
            PCO_fr_NMVOC_use = Input_Opt%LPCO_NMVOC,                         &
            PCO_fr_NMVOC     = PCO_fr_NMVOC(I,J,L),                          &
            State_Met        = State_Met,                                    &
            State_Chm        = State_Chm                                    )

       ! Update the array of rate constants for the KPP solver
       CALL Update_RCONST()

       ! Stop timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End(                                                    &
               timerName = "     RCONST",                                    &
               inLoop    = .TRUE.,                                           &
               threadNum = threadNum,                                        &
               RC        =  RC                                              )
       ENDIF

       !=====================================================================
       ! Call the KPP integrator
       !=====================================================================

       ! Start timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_Start(                                                  &
               timerName = "     Integrate 1",                               &
               inLoop    =  .TRUE.,                                          &
               threadNum = threadNum,                                        &
               RC        = RC                                               )
       ENDIF

       ! Integrate the mechanism forward in time
       CALL Integrate(                                                       &
            TIN      = 0.0_dp,                                               &
            TOUT     = dtChem,                                               &
            ICNTRL_U = ICNTRL,                                               &
            IERR_U   = IERR                                                 )

       ! Stop timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End(                                                    &
               timerName = "     Integrate 1",                               &
               inLoop    = .TRUE.,                                           &
               threadNum = ThreadNum,                                        &
               RC        = RC                                               )
       ENDIF

       ! Trap potential errors
       IF ( IERR /= 1 ) failed = .TRUE.

       !=====================================================================
       ! HISTORY: Archive KPP solver diagnostics
       !=====================================================================
       IF ( State_Diag%Archive_KppDiags ) THEN

          ! # of integrator calls
          IF ( State_Diag%Archive_KppIntCounts ) THEN
             State_Diag%KppIntCounts(I,J,L) = ISTATUS(1)
          ENDIF

          ! # of times Jacobian was constructed
          IF ( State_Diag%Archive_KppJacCounts ) THEN
             State_Diag%KppJacCounts(I,J,L) = ISTATUS(2)
          ENDIF

          ! # of internal timesteps
          IF ( State_Diag%Archive_KppTotSteps ) THEN
             State_Diag%KppTotSteps(I,J,L) = ISTATUS(3)
          ENDIF

          ! # of accepted internal timesteps
          IF ( State_Diag%Archive_KppTotSteps ) THEN
             State_Diag%KppAccSteps(I,J,L) = ISTATUS(4)
          ENDIF

          ! # of rejected internal timesteps
          IF ( State_Diag%Archive_KppTotSteps ) THEN
             State_Diag%KppRejSteps(I,J,L) = ISTATUS(5)
          ENDIF

          ! # of LU-decompositions
          IF ( State_Diag%Archive_KppLuDecomps ) THEN
             State_Diag%KppLuDecomps(I,J,L) = ISTATUS(6)
          ENDIF

          ! # of forward and backwards substitutions
          IF ( State_Diag%Archive_KppSubsts ) THEN
             State_Diag%KppSubsts(I,J,L) = ISTATUS(7)
          ENDIF

          ! # of singular-matrix decompositions
          IF ( State_Diag%Archive_KppSmDecomps ) THEN
             State_Diag%KppSmDecomps(I,J,L) = ISTATUS(8)
          ENDIF
       ENDIF

       ! Convert CO, CO2, CH4 to molec/cm3 for the KPP solver
       CALL carbon_ConvertMolecCm3ToKg(                                      &
            I            = I,                                                &
            J            = J,                                                &
            L            = L,                                                &
            id_CH4       = id_CH4,                                           &
            id_CO        = id_CO,                                            &
            id_COch4     = id_COch4,                                         &
            id_COnmvoc   = id_COnmvoc,                                       &
            id_CO2       = id_CO2,                                           &
            xnumol_CO    = xnumol_CO,                                        &
            xnumol_CH4   = xnumol_CH4,                                       &
            xnumol_CO2   = xnumol_CO2,                                       &
            State_Chm    = State_Chm,                                        &
            State_Met    = State_Met                                        )

       IF ( Input_Opt%useTimers ) THEN

          ! Stop main KPP timer
          CALL Timer_End(                                                    &
               timerName  = "  -> KPP",                                      &
               inLoop     = .TRUE.,                                          &
               threadNum  = threadNum,                                        &
               RC         = RC                                              )

          ! Start Prod/Loss timer
          CALL Timer_Start(                                                  &
               timerName = "  -> Prod/loss diags",                           &
               inLoop    = .TRUE.,                                           &
               threadNum = threadNum,                                        &
               RC        = RC                                               )

       ENDIF

#ifdef ACTIVATE_TAGGED_SPECIES
       !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       !@@@ TAGGED SPECIES HANDLING
       !@@@ Handle trop loss by OH for regional CO species
       !@@@ This is turned off by default -- needs further work
       !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
            IF (NA .ne. 16)                                                   &
                 Spc(N)%Conc(I,J,L) = Spc(N)%Conc(I,J,L) *                    &
                 ( 1e+0_fp - K_TROP(2) * C(FixedOH) * DTCHEM )

            !-----------------------------------------------------
            ! HISTORY (aka netCDF diagnostics)
            !
            ! Loss of CO by OH for "tagged" species
            !-----------------------------------------------------

            ! Units: [kg/s]
            IF ( State_Diag%Archive_Loss ) THEN
               State_Diag%Loss(I,J,L,N) = Spc(N)%Conc(I,J,L) * k_Trop(2)     &
                    * C(FixedOH) * DTCHEM
               !   C(ind_CO2_OH) / DTCHEM  &
               !   * State_Met%AIRVOL(I,J,L) * 1e+6_fp / XNUMOL_CO
            ENDIF
         ENDDO
      ENDIF
#endif

       !=====================================================================
       ! HISTORY (aka netCDF diagnostics)
       !
       ! Production and loss of CO
       !
       ! NOTE: Call functions in KPP/carbon/carbon_Funcs.F90 so
       ! that we avoid bringing in KPP species indices into this module.
       ! This avoids compile-time dependency errors.
       !=====================================================================

       ! Production of CO2 from CO oxidation [molec/cm3/s]
       IF ( Input_Opt%LCHEMCO2 ) THEN
          IF ( State_Diag%Archive_ProdCO2fromCO ) THEN
             State_Diag%ProdCO2fromCO(I,J,L) =                               &
                carbon_Get_CO2fromOH_Flux( dtChem )
          ENDIF
       ENDIF

       ! Production of CO from CH4
       IF ( State_Diag%Archive_ProdCOfromCH4 ) THEN
          State_Diag%ProdCOfromCH4(I,J,L) =                                  &
             carbon_Get_COfromCH4_Flux( dtChem )
       ENDIF

       ! Units: [kg/s] Production of CO from NMVOCs
       IF ( State_Diag%Archive_ProdCOfromNMVOC ) THEN
          State_Diag%ProdCOfromNMVOC(I,J,L) =                                &
             carbon_Get_COfromNMVOC_Flux( dtChem )
       ENDIF

#ifdef ACTIVATE_TAGGED_SPECIES
       !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       !@@@ TAGGED SPECIES HANDLING
       !@@@ Loss of CO by OH -- tagged species; Units: [kg/s]
       !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       IF ( State_Diag%Archive_Loss ) THEN
          State_Diag%Loss(I,J,L,16) = ( COfrom_OH / STTCO / DTCHEM )
       ENDIF
#endif

       ! Stop Prod/Loss timer
       IF ( Input_Opt%useTimers ) THEN
          CALL Timer_End(                                                    &
               timerName =  "  -> Prod/loss diags",                          &
               inLoop    = .TRUE.,                                           &
               threadNum = threadNum,                                        &
               RC        = RC                                               )
       ENDIF

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    IF ( failed ) THEN
       errMsg = 'KPP integration failed!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

#ifdef ACTIVATE_TAGGED_SPECIES
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !@@@ TAGGED SPECIES HANDLING
    !@@@ Allocate the CH4 chemistry sink to different tagged CH4 species
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    IF ( Input_Opt%LSPLIT .and. id_CH4 > 0 ) THEN
       CALL CH4_Distrib_Carbon(                                              &
            PrevCh4    = PrevCh4,                                            &
            Input_Opt  = Input_Opt,                                          &
            State_Chm  = State_Chm,                                          &
            State_Grid = State_Grid                                         )
    ENDIF
#endif

    ! Free pointers for safety's sake
    Spc => NULL()

  END SUBROUTINE Chem_Carbon_Gases
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ReadChemInputFields
!
! !DESCRIPTION: Calls HCO_EvalFld to read the various oxidant and other
!  non-emissions fields needed for the carbon simulation
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE ReadChemInputFields( Input_Opt,    State_Grid,   State_Met,     &
                                  Global_OH,    Global_Cl,    LCH4_by_OH,    &
                                  LCO_in_Strat, PCO_in_Strat, PCO_fr_CH4,    &
                                  PCO_fr_NMVOC, RC                          )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt        ! Input Options object
    TYPE(GrdState), INTENT(IN)  :: State_Grid       ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met        ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT) :: Global_OH(                                &
                                    State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Grid%NZ)  ! OH conc
    REAL(fp),       INTENT(OUT) :: Global_Cl(                                &
                                    State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Grid%NZ)  ! Cl conc
    REAL(fp),       INTENT(OUT) :: LCH4_by_OH(                               &
                                    State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Grid%NZ)  ! L(CH4) by OH
    REAL(fp),       INTENT(OUT) :: LCO_in_Strat(                             &
                                    State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Grid%NZ)  ! L(CO) in strat
    REAL(fp),       INTENT(OUT) :: PCO_in_Strat(                             &
                                    State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Grid%NZ)  ! P(CO) in strat
    REAL(fp),       INTENT(OUT) :: PCO_fr_CH4(                               &
                                    State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Grid%NZ)  ! P(CO) from CH4
    REAL(fp),       INTENT(OUT) :: PCO_fr_NMVOC(                             &
                                    State_Grid%NX,                           &
                                    State_Grid%NY,                           &
                                    State_Grid%NZ)  ! P(CO) from NMVOC
    INTEGER,        INTENT(OUT) :: RC               ! Success/failure
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: found

    ! Strings
    CHARACTER(LEN=63)  :: dgnName
    CHARACTER(LEN=255) :: thisLoc
    CHARACTER(LEN=512) :: errMsg

    !========================================================================
    !  ReadChemInputFields begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    found   = .FALSE.
    dgnName = ''
    errMsg  = ''
    thisLoc = &
     ' -> at ReadInputChemFields (in module GeosCore/carbon_gases_mod.F90)'

    !------------------------------------------------------------------------
    ! Loss frequencies of CH4
    ! Input via HEMCO ("CH4_LOSS" container) as [1/s]
    !------------------------------------------------------------------------
    DgnName = 'CH4_LOSS'
    CALL HCO_GC_EvalFld( Input_Opt,  State_Grid, DgnName,                    &
                         LCH4_by_OH, RC,         found=found                )
    IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
       errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! Cl concentration:
    ! Input via HEMCO ("SpeciesConc" collection) as [mol/mol dry]
    ! Convert to [molec/cm3] below
    !------------------------------------------------------------------------
    DgnName = 'GLOBAL_Cl'
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName,                     &
                         Global_Cl, RC,         found=found                 )
    IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
       errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Convert orignal units [mol/mol dry air] to [molec/cm3]
    Global_Cl = ( Global_Cl * State_Met%AirDen ) * toMolecCm3

    !------------------------------------------------------------------------
    ! OH concentration: from GEOS-Chem v5 or GEOS-Chem 10yr benchmark
    !------------------------------------------------------------------------
    IF ( useGlobOH ) THEN

       ! NOTE: Container name is GLOBAL_OH for both data sets!
       DgnName = 'GLOBAL_OH'
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, DgnName,               &
                            Global_OH, RC,         found=found           )
       IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
          errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       IF ( useGlobOHv5 ) THEN

          ! If we are using Global_OH from GEOS-Chem v5 (e.g. for the IMI or
          ! methane simulations) then convert OH from [kg/m3] to [molec/cm3]
          Global_OH = Global_OH * xnumol_OH / CM3perM3

       ELSE

          ! If we are using OH from a 10-year benchmark ("SpeciesConc")
          ! then convert OH [mol/mol dry air] to [molec/cm3]
          Global_OH = ( Global_OH * State_Met%AirDen ) * toMolecCm3

       ENDIF

    ENDIF

    !------------------------------------------------------------------------
    ! P(CO) from GMI:
    ! Input via HEMCO ("GMI_PROD_CO" field) as [v/v/s]
    ! Units will be converted in carbon_ComputeRateConstants
    !------------------------------------------------------------------------
    DgnName = 'GMI_PROD_CO'
    CALL HCO_GC_EvalFld( Input_Opt,    State_Grid, DgnName,                  &
                         PCO_in_Strat, RC,         found=found              )
    IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
       errMsg = 'Cannot get pointer to ' // TRIM( DgnName )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! L(CO) from GMI
    ! Input via HEMCO ("GMI_LOSS_CO" field) as [1/s]
    !------------------------------------------------------------------------
    DgnName = 'GMI_LOSS_CO'
    CALL HCO_GC_EvalFld( Input_Opt,    State_Grid, DgnName,                  &
                         LCO_in_Strat, RC,         found=found              )
    IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
       errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !------------------------------------------------------------------------
    ! P(CO) from CH4
    ! Input via HEMCO ("ProdCOfromCH4 field") as [molec/cm3/s]
    !------------------------------------------------------------------------
    IF ( Input_Opt%LPCO_CH4 ) THEN
       DgnName = 'PCO_CH4'
       CALL HCO_GC_EvalFld( Input_Opt,  State_Grid, DgnName,                  &
                            PCO_fr_CH4, RC,         found=found              )
       IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
          errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    !------------------------------------------------------------------------
    ! P(CO) from NMVOC
    ! Input via HEMCO ("ProdCOfromNMVOC" field) as [molec/cm3/s]
    !------------------------------------------------------------------------
    IF ( Input_Opt%LPCO_NMVOC ) THEN
       DgnName = 'PCO_NMVOC'
       CALL HCO_GC_EvalFld( Input_Opt,    State_Grid, DgnName,               &
                            PCO_fr_NMVOC, RC,         found=found           )
       IF ( RC /= GC_SUCCESS .or. .not. found ) THEN
          errMsg = 'Cannot get pointer to HEMCO field ' // TRIM( DgnName )
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

  END SUBROUTINE ReadChemInputFields
!EOC
#ifdef ACTIVATE_TAGGED_SPECIES
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_distrib
!
! !DESCRIPTION: Allocates the chemistry sink to different emission species.
!  (Only called if there are tagged CH4 species.)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CH4_Distrib( PrevCh4, Input_Opt, State_Chm, State_Grid )
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
    TYPE(OptInput), INTENT(IN)    :: Input_Opt       ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid      ! Grid State object
    REAL(fp)                      :: PREVCH4(                                &
                                      State_Grid%NX,                         &
                                      State_Grid%NY,                         &
                                      State_Grid%NZ) ! CH4 befire chemistry
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm       ! Chemistry State object
!
!
! !REMARKS:
!  This routine is only used with tagged CH4 species.
!
! !REVISION HISTORY:
!  See the Git version history
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

  END SUBROUTINE CH4_DISTRIB
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_carbon_gases
!
! !DESCRIPTION: Allocates and zeroes module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Carbon_Gases( Input_Opt,  State_Chm, State_Diag,             &
                                State_Grid, RC                                )
!
! !USES:
!
    USE gckpp_Global,   ONLY : MW, SR_MW, HENRY_CR, HENRY_K0
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState, Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
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
    INTEGER            :: KppId,  N

    ! Strings
    CHARACTER(LEN=255) :: errMsg
    CHARACTER(LEN=255) :: thisLoc

    !========================================================================
    ! Initialize
    !========================================================================
    RC         = GC_SUCCESS
    xnumol_CH4 = 1.0_fp
    xnumol_CO  = 1.0_fp
    xnumol_CO2 = 1.0_fp
    xnumol_OH  = 1.0_fp
    errMsg     = ''
    thisLoc    = &
     ' -> at Init_Carbon_Gases (in module GeosCore/carbon_gases_mod.F90)'

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
    id_OH      = Ind_( 'FixedOH' )

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

          ! Also define xnumol factors = molec/kg ratios
          ! These are useful in unit conversions
          IF ( N == id_CH4 ) xnumol_CH4 = AVO / ( MW(KppId) * 1.0e-3_fp )
          IF ( N == id_CO  ) xnumol_CO  = AVO / ( MW(KppId) * 1.0e-3_fp )
          IF ( N == id_CO2 ) xnumol_CO2 = AVO / ( MW(KppId) * 1.0e-3_fp )
          IF ( N == id_OH  ) xnumol_OH  = AVO / ( MW(KppId) * 1.0e-3_fp )
       ENDIF
    ENDDO

    !========================================================================
    ! Initialize variables for COchemistry
    !========================================================================
    ALLOCATE( sumOfCosSza( State_Grid%NX, State_Grid%NY ), STAT=RC )
    CALL GC_CheckVar( 'carbon_gases_mod.F90:sumOfCosSza', 0, RC )
    IF ( RC /= GC_SUCCESS ) RETURN
    sumOfCosSza = 0.0_fp

    !========================================================================
    ! Initialize variables for CH4 chemistry
    !========================================================================
    IF ( id_CH4 > 0 ) THEN
       ALLOCATE( CH4_EMIS_J( State_Grid%NX, State_Grid%NY, N_CH4_DIAGS ),    &
            STAT=RC )
       CALL GC_CheckVar( 'carbon_gases_mod.F90:CH4_EMIS', 0, RC )
       IF ( RC /= GC_SUCCESS ) RETURN
       CH4_EMIS_J = 0.0_fp
    ENDIF

    !========================================================================
    ! Initialize variables for CO2 chemistry
    !========================================================================
    ! none yet

    !========================================================================
    ! Initialize variables for OCS chemistry
    !========================================================================
    ! none yet

  END SUBROUTINE Init_Carbon_Gases
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_carbon_gases
!
! !DESCRIPTION: Subroutine CLEANUP\_GLOBAL\_CH4 deallocates module arrays.
!  (bmy, 1/16/01)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Carbon_Gases( RC )
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
    ! Cleanup_Carbon_Gases begins here!
    !=================================================================

    ! Initialize
    RC = GC_SUCCESS

    ! Deallocate
    IF ( ALLOCATED( CH4_EMIS_J ) ) THEN
       DEALLOCATE( CH4_EMIS_J, STAT=RC )
       CALL GC_CheckVar( 'carbon_gases_mod.F90:CH4_EMIS', 2, RC )
       RETURN
    ENDIF

    IF ( ALLOCATED( sumOfCosSza ) ) THEN
       DEALLOCATE( sumOfCosSza, STAT=RC )
       CALL GC_CheckVar( 'carbon_gases_mod.F90:sumOfCosSza', 2, RC )
       RETURN
    ENDIF

  END SUBROUTINE Cleanup_Carbon_Gases
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
  SUBROUTINE Calc_Diurnal( State_Grid, State_Met, OHdiurnalFac )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE TIME_MOD,       ONLY : ITS_A_NEW_DAY
    USE TIME_MOD,       ONLY : GET_MINUTE,    GET_SECOND,      GET_HOUR
    USE TIME_MOD,       ONLY : GET_TS_CHEM,   GET_DAY_OF_YEAR, GET_LOCALTIME
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid        ! Grid State object
    TYPE(MetState), INTENT(IN)  :: State_Met         ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
    REAL(fp),       INTENT(OUT) :: OHdiurnalFac(                             &
                                    State_Grid%NX,                           &
                                    State_Grid%NY)   ! OH diurnal scaling [1]
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
    ! SAVEd scalars
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Scalars
    INTEGER            :: I, J, N, NDYSTEP
    INTEGER            :: SECOND,  MINUTE, TS_SUN
    REAL*8             :: GMT_MID, TIMLOC, FACTOR
    REAL*8             :: R,       AHR,    DEC
    REAL*8             :: YMID_R,  SUNTMP_MID
    REAL*8             :: dtChem,  timestepsPerDay
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
       sumOfCosSza = 0d0

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

             ! sumOfCosSza is the sum of SUNTMP_MID at location (I,J)
             ! Do not include negative values of SUNTMP_MID
             sumOfCosSza(I,J) = sumOfCosSza(I,J) + MAX( SUNTMP_MID, 0d0 )

          ENDDO
          ENDDO
          !$OMP END PARALLEL DO
       ENDDO

       ! Reset first-time flag
       FIRST = .FALSE.
    ENDIF

    !========================================================================
    ! Calculate the OH diurnal scaling factor
    !========================================================================

    ! Chemistry timestep [s] and timesteps per day
    dtChem          = GET_TS_CHEM()
    timestepsPerDay = 86400.0_fp / dtChem

    ! Loop over surface grid boxes
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, J                                                     )&
    !$OMP COLLAPSE( 2                                                       )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Initialize loop variables
       OHdiurnalFac(I,J) = 0.0_fp

       ! Scaling factor for OH diurnal cycles - zero at night
       IF ( State_Met%SUNCOSmid(I,J) > 0.0_fp  .and.                         &
            sumOfCosSza(I,J)         > 0.0_fp ) THEN
          OHdiurnalFac(I,J) = State_Met%SUNCOSmid(I,J)                       &
                            / sumOfCosSza(I,J)                               &
                            * timestepsPerDay
       ENDIF
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Calc_Diurnal
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: InquireGlobalOHversion
!
! !DESCRIPTION: Determines if we are using global OH from a GEOS-Chem 10-year
!  benchmark, or from GEOS-Chem v5 (needed for IMI & methane simulations).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE InquireGlobalOHversion( Input_Opt, RC )
!
! !USES:
!
    USE CharPak_Mod,          ONLY : To_UpperCase
    USE ErrCode_Mod
    USE Hco_Utilities_Gc_Mod, ONLY : HCO_GC_GetOption
    USE Input_Opt_Mod,        ONLY : OptInput
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC   ! Success or failure
!
! !REMARKS:
!  This has to be called on the first chemistry timestep, because the HEMCO
!  Configuration file (HEMCO_Config.rc) will not have been read when
!  Init_Carbon_Gases is called.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: isOHon, isOHv5on

    ! Strings
    CHARACTER(LEN=255) :: optVal
    CHARACTER(LEN=255) :: thisLoc
    CHARACTER(LEN=512) :: errMsg

    !========================================================================
    ! InquireGlobalOHversion begins here!
    !========================================================================

    ! Initialize
    RC       = GC_SUCCESS
    errMsg   = ''
    thisLoc  = &
     ' -> at InquireGlobalOHversion (in module GeosCore/carbon_gases_mod.F90)'

    ! Test if any Global OH option is on
    optVal   = HCO_GC_GetOption( "GLOBAL_OH", extNr=0 )
    isOHon   = ( To_UpperCase( TRIM( optVal ) ) == 'TRUE' )

    ! Test for GEOS-Chem v5 OH
    ! =T
    optVal   = HCO_GC_GetOption( "GLOBAL_OH_GCv5", extNr=0 )
    isOHv5On = ( To_UpperCase( TRIM( optVal ) ) == 'TRUE' )

    ! Set a global variable to determine which OH to use
    useGlobOH   = isOHon
    useGlobOHv5 = isOHv5on

    IF ( Input_Opt%amIRoot ) THEN
       IF ( useGlobOH ) THEN
          IF ( useGlobOHv5 ) THEN
             WRITE( 6, 100 ) 'GLOBAL_OH_GCv5'
          ELSE
             WRITE( 6, 100 ) 'GLOBAL_OH'
          ENDIF
 100   FORMAT( 'Carbon_Gases: Using global OH oxidant field option: ', a )

       ELSE
          WRITE( 6, 110 )
 110      FORMAT( 'Carbon_Gases: Global OH is set to zero!' )
       ENDIF

    ENDIF

  END SUBROUTINE InquireGlobalOHversion
!EOC
END MODULE Carbon_Gases_Mod
