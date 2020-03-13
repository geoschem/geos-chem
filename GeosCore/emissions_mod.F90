!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: emissions_mod.F90
!
! !DESCRIPTION: Module emissions\_mod.F90 is a wrapper module to interface
! GEOS-Chem and HEMCO. It basically just calls the GEOS-Chem - HEMCO interface
! routines. For some specialty sims, a few additional steps are required that
! are also executed here.
!\\
!\\
! !INTERFACE:
!
MODULE Emissions_Mod
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Emissions_Init
  PUBLIC :: Emissions_Run
  PUBLIC :: Emissions_Final
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRIVATE :: MMR_Compute_Flux
!
! !REVISION HISTORY:
!  27 Aug 2014 - C. Keller   - Initial version.
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  ! Species ID flags
  INTEGER :: id_CH4
  LOGICAL :: doMaintainMixRatio

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissions_init
!
! !DESCRIPTION: Subroutine EMISSIONS\_INIT calls the HEMCO - GEOS-Chem
! interface initialization routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Emissions_Init( Input_Opt, State_Chm, State_Grid, State_Met,  &
                             RC,        HcoConfig )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCOI_GC_Main_Mod,   ONLY : HCoi_GC_Init
    USE HCO_Types_Mod,      ONLY : ConfigObj
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Chm_Mod,      ONLY : Ind_
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(ChmState),  INTENT(IN   )          :: State_Chm  ! Chemistry state
    TYPE(GrdState),  INTENT(IN   )          :: State_Grid ! Grid state
    TYPE(MetState),  INTENT(IN   )          :: State_Met  ! Met state
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(OptInput),  INTENT(INOUT)          :: Input_Opt  ! Input opts
    INTEGER,         INTENT(INOUT)          :: RC         ! Failure or success
    TYPE(ConfigObj), POINTER,      OPTIONAL :: HcoConfig  ! HEMCO config object
!
! !REVISION HISTORY:
!  27 Aug 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! EMISSIONS_INIT begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Emissions_Init (in module GeosCore/emissions_mod.F90)'

    ! Define species ID flags for use in routines below
    id_CH4   = Ind_('CH4'  )

    ! Are we including a species for which the global mixing ratio should
    ! remain constant?
    doMaintainMixRatio = ( Ind_('GlobEmis90dayTracer') > 0 .OR. &
                           Ind_('GlobNH90dayTracer'  ) > 0 .OR. &
                           Ind_('GlobSH90dayTracer'  ) > 0 )

    ! Initialize the HEMCO environment for this GEOS-Chem run.
    CALL HCOI_GC_Init( Input_Opt, State_Chm, State_Grid, State_Met, &
                       RC,        HcoConfig=HcoConfig )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HCOI_GC_Init"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Emissions_Init
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissions_run
!
! !DESCRIPTION: Subroutine EMISSIONS\_RUN calls the HEMCO - GEOS-Chem
! interface run routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Emissions_Run( Input_Opt, State_Chm, State_Diag, &
                            State_Grid, State_Met, EmisTime,  Phase, RC   )
!
! !USES:
!
    USE CARBON_MOD,         ONLY : EMISSCARBON
    USE CO2_MOD,            ONLY : EMISSCO2
    USE ErrCode_Mod
    USE GLOBAL_CH4_MOD,     ONLY : EMISSCH4
    USE HCOI_GC_MAIN_MOD,   ONLY : HCOI_GC_RUN
    USE Input_Opt_Mod,      ONLY : OptInput
    USE Precision_Mod
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE Time_Mod,           ONLY : Get_Ts_Emis
    Use SfcVmr_Mod,         Only : FixSfcVmr_Run
    USE MERCURY_MOD,        ONLY : EMISSMERCURY
    USE Pops_Mod,           ONLY : GetPopsDiagsFromHemco
#ifdef TOMAS
    USE CARBON_MOD,         ONLY : EMISSCARBONTOMAS !jkodros
    USE SULFATE_MOD,        ONLY : EMISSSULFATETOMAS !jkodros
#endif
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN   )  :: EmisTime   ! Emissions in this time step
    INTEGER,        INTENT(IN   )  :: Phase      ! Run phase
    TYPE(GrdState), INTENT(IN   )  :: State_Grid ! Grid State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT)  :: State_Chm  ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT)  :: State_Diag ! Diagnostics State object
    TYPE(MetState), INTENT(INOUT)  :: State_Met  ! Meteorology State object
    TYPE(OptInput), INTENT(INOUT)  :: Input_Opt  ! Input Options object
    INTEGER,        INTENT(INOUT)  :: RC         ! Failure or success
!
! !REVISION HISTORY:
!  27 Aug 2014 - C. Keller   - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! EMISSIONS_RUN begins here!
    !=================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Emissions_Run (in module GeosCore/emissions_mod.F90)'

    ! Run HEMCO. Phase 1 will only update the HEMCO clock and the
    ! HEMCO data list, phase 2 will perform the emission calculations.
    CALL HCOI_GC_Run( Input_Opt, State_Chm, State_Grid, &
                      State_Met, EmisTime,  Phase,     RC          )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HCOI_GC_Run"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Exit if Phase 0 or 1, or if it is a GEOS-Chem dry-run
    IF ( ( Phase == 1 ) .or. ( Phase == 0 ) .or. Input_Opt%DryRun ) RETURN

    ! Call carbon emissions module to make sure that sesquiterpene
    ! emissions calculated in HEMCO (SESQ) are passed to the internal
    ! species array in carbon, as well as to ensure that POA emissions
    ! are correctly treated.
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM   .or. &
         Input_Opt%ITS_AN_AEROSOL_SIM ) THEN
       CALL EmissCarbon( Input_Opt, State_Grid, State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "EmissCarbon"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

#ifdef TOMAS
    ! Call TOMAS emission routines (JKodros 6/2/15)
    CALL EmissCarbonTomas( Input_Opt, State_Chm, State_Grid, State_Met, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "EmissCarbonTomas"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    CALL EmissSulfateTomas( Input_Opt, State_Chm, State_Grid, State_Met, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "EmissSulfateTomas"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
#endif

    ! For the CO2 simulation, we manually add the chemical production of CO2
    ! from CO oxidation (which is listed as a non-chemical source in HEMCO)
    ! to State_Chm%Species.  This is done in EmissCO2.  All other CO2 emissions
    ! (as of GEOS-Chem 12.0.2) are now added via HEMCO, and diagnostics for
    ! these quantities are saved out via HEMCO diagnostics. (bmy, 10/18/18)
    IF ( Input_Opt%ITS_A_CO2_SIM ) THEN
       CALL EmissCO2( Input_Opt, State_Chm, State_Diag, State_Grid, &
                      State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "EmissCO2"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! For CH4 simulation or if CH4 is defined, call EMISSCH4.
    ! This will get the individual CH4 emission terms (gas, coal, wetlands,
    ! ...) and write them into the individual emissions arrays defined in
    ! global_ch4_mod (CH4_EMIS). Emissions are all done in mixing_mod, the
    ! call to EMISSCH4 is for backwards consistency.  This is especially
    ! needed to do the analytical inversions.  NOTE: The CH4 manual
    ! diagnostics are no longer used to force-feed the ND58 bpch diagnostics
    ! becasue we now archive the exact same quantities to the HEMCO
    ! diagnostics output. (bmy, mps, 10/19/18)
    IF ( Input_Opt%ITS_A_CH4_SIM .OR.            &
       ( id_CH4 > 0 .and. Input_Opt%LCH4EMIS ) ) THEN
       CALL EmissCh4( Input_Opt, State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "EmissCH4"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! For mercury, use old emissions code for now
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN
       CALL EmissMercury( Input_Opt, State_Chm, State_Diag, State_Grid, &
                          State_Met, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "EmissMercury"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! For the POPS simulation, copy values from several HEMCO-based manual
    ! diagnostics (defined in hcoi_gc_diagn_mod.F90) from the HEMCO state
    ! object into the State_Diag object.  This will allow us to save these
    ! fields to netCDF output via HISTORY. (bmy, 10/15/18)
    IF ( Input_Opt%ITS_A_POPS_SIM ) THEN
       CALL GetPopsDiagsFromHemco( Input_Opt, State_Diag, RC )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "GetPopsDiagsFromHemco"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Prescribe some concentrations if needed
    IF ( Input_Opt%ITS_A_FULLCHEM_SIM ) THEN

       ! Set other (non-UCX) fixed VMRs
       If ( Input_Opt%LEMIS ) Then
          CALL FixSfcVMR_Run( Input_Opt, State_Chm, State_Grid, State_Met, RC )

          ! Trap potential errors
          IF ( RC /= GC_SUCCESS ) THEN
             ErrMsg = 'Error encountered in "FixSfcVmr_Run"!'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
          ENDIF

       endif

    ENDIF

    IF ( doMaintainMixRatio ) THEN

       ! Compute the surface flux needed to restore the total burden
       CALL MMR_Compute_Flux( Input_Opt, State_Chm, State_Grid, State_Met, RC )

    ENDIF

   END SUBROUTINE EMISSIONS_RUN
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissions_final
!
! !DESCRIPTION: Subroutine EMISSIONS\_FINAL calls the HEMCO - GEOS-Chem
! interface finalization routines.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Emissions_Final( Error, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCOI_GC_Main_Mod, ONLY : HCOI_GC_Final
!
! !INPUT PARAMETERS:
!
    LOGICAL, INTENT(IN)  :: Error      ! Cleanup arrays after crash?
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: RC         ! Success or failure?
!
! !REVISION HISTORY:
!  27 Aug 2014 - C. Keller    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc

    !=================================================================
    ! EMISSIONS_FINAL begins here!
    !=================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at HCOI_GC_Final (in module GeosCore/hcoi_gc_final_mod.F90)'

    ! Finalize HEMCO
    CALL HCOI_GC_Final( Error, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "HCOI_GC_Final"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF
  END SUBROUTINE Emissions_Final
!EOC
!------------------------------------------------------------------------------
!                  Harvard-NASA Emissions Component (HEMCO)                   !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: MMR_Compute_Flux
!
! !DESCRIPTION: Subroutine MMR\_Compute\_Flux computes the surface flux
!  needed to maintain a given mixing ratio value.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE MMR_Compute_Flux( Input_Opt, State_Chm, &
                               State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_INTERFACE_MOD,  ONLY : HcoState
    USE HCO_EmisList_Mod,   ONLY : HCO_GetPtr
    USE HCO_STATE_MOD,      ONLY : HCO_GetHcoID
    USE Input_Opt_Mod,      ONLY : OptInput
    USE PhysConstants
    USE Species_Mod,        ONLY : Species
    USE State_Chm_Mod,      ONLY : ChmState
    USE State_Grid_Mod,     ONLY : GrdState
    USE State_Met_Mod,      ONLY : MetState
    USE UnitConv_Mod,       ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object!
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY:
!  31 Jan 2019 - M. Sulprizio- Initial version, modified from MMR code in
!                              TR_GridCompMod.F90 from GEOS model
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                    :: I, J, L, N, nAdvect
    REAL(fp)                   :: Total_Spc
    REAL(fp)                   :: Total_Area
    REAL(fp)                   :: YMID

    ! Arrays
    REAL(fp)                   :: Flux(State_Grid%NX,State_Grid%NY)
    REAL(fp)                   :: Mask(State_Grid%NX,State_Grid%NY)

    ! Strings
    CHARACTER(LEN=63)          :: OrigUnit
    CHARACTER(LEN=255)         :: ThisLoc
    CHARACTER(LEN=512)         :: ErrMsg

    ! Pointers
    REAL(fp),        POINTER   :: Spc(:,:,:,:)
    TYPE(Species),   POINTER   :: SpcInfo
!
! !DEFINED PARAMETERS:
!
    ! Hardcode global burden to 100 ppbv for now
    REAL(fp),        PARAMETER :: GlobalBurden = 1.0e-7_fp ! [v/v]

    !=================================================================
    ! MMR_Compute_Flux begins here!
    !=================================================================

    ! Initialize
    RC          = GC_SUCCESS
    ErrMsg      = ''
    ThisLoc     = ' -> at MMR_Compute_Flux (in module GeosCore/emissions_mod.F90)'

    ! Point to chemical species array [kg/kg dry air]
    Spc        => State_Chm%Species

    ! Number of advected species
    nAdvect     = State_Chm%nAdvect

    !=======================================================================
    ! Convert species units to v/v dry
    !=======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'v/v dry', RC, OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error (kg/kg dry -> v/v dry)'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Compute the surface flux needed to restore the total burden
    !=======================================================================

    ! Loop over species
    DO N = 1, nAdvect

       ! Point to the Species Database entry for species N
       SpcInfo => State_Chm%SpcData(N)%Info

       ! Only do calculation if this is an MMR tracer
       IF ( TRIM(SpcInfo%Name) == 'GlobEmis90dayTracer' .OR. &
            TRIM(SpcInfo%Name) == 'NHEmis90dayTracer'   .OR. &
            TRIM(SpcInfo%Name) == 'SHEmis90dayTracer'   ) THEN

          ! Initialize
          Total_Spc   = 0.0e+0_fp
          Total_Area  = 0.0e+0_fp
          Mask        = 1.0e+0_fp

          ! Loop over grid boxes
          DO L = 1, State_Grid%NZ
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Compute mol of Tracer needed to achieve the desired value
             Total_Spc = Total_Spc + &
                ( GlobalBurden - State_Chm%Species(I,J,L,N)) * &
                (State_Met%AIRNUMDEN(I,J,L)/ AVO) * State_Met%AIRVOL(I,J,1)

             ! To distribute it uniformly on the surface, compute the total
             ! area [m2]
             IF ( L == 1 ) THEN

                ! Latitude of grid box
                YMID = State_Grid%YMid(I,J)

                ! Define mask if needed
                IF (TRIM(SpcInfo%Name)=='NHEmis90dayTracer') THEN
                   IF ( YMID <  0.0 ) MASK(I,J) = 0.0e+0_fp
                ELSE IF (TRIM(SpcInfo%Name)=='SHEmis90dayTracer') THEN
                   IF ( YMID >= 0.0 ) MASK(I,J) = 0.0e+0_fp
                ENDIF
                Total_Area = Total_Area + State_Grid%Area_M2(I,J) * MASK(I,J)
             ENDIF

          ENDDO
          ENDDO
          ENDDO

          ! Compute flux [mol/m2]
          Flux(:,:) = ( Total_Spc / Total_Area ) * MASK(:,:)

          ! Update species concentrations [mol/mol]
          Spc(:,:,1,N) = Spc(:,:,1,N) + Flux(:,:) * &
             AVO / ( State_Met%BXHEIGHT(:,:,1) * State_Met%AIRNUMDEN(:,:,1) )

       ENDIF ! MMR tracer

       ! Free pointers
       SpcInfo => NULL()

    ENDDO

    !=======================================================================
    ! Convert species units back to original unit
    !=======================================================================
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit,  RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE MMR_Compute_Flux
!EOC
END MODULE Emissions_Mod
