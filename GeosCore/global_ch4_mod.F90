!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: global_ch4_mod.F90
!
! !DESCRIPTION: Module GLOBAL\_CH4\_MOD contains variables and routines for
!  simulating CH4 chemistry in the troposphere.
!\\
!\\
! !INTERFACE:
!
MODULE GLOBAL_CH4_MOD
!
! !USES:
!
  USE PhysConstants, ONLY : AVO, AIRMW
  USE PRECISION_MOD       ! For GEOS-Chem Precision (fp, f4, f8)

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: EMISSCH4
  PUBLIC :: CHEMCH4
  PUBLIC :: INIT_GLOBAL_CH4
!
! !PUBLIC DATA MEMBERS:
!
  REAL(fp), PARAMETER,   PUBLIC :: XNUMOL_CH4 = AVO / 16d-3 ! hard-coded MW
!
! !REVISION HISTORY:
!  17 Jan 2001- J. Wang, B. Duncan, R. Yantosca -- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
  !========================================================================
  ! Module Variables:
  ! XNUMOL_CH4 : Molecules CH4 / kg CH4                     [molec/kg]
  !========================================================================

  REAL(fp), PARAMETER   :: XNUMOL_OH = AVO / 17e-3_fp  ! molec OH / kg OH
                                                       ! hard-coded MW
  REAL(fp), PARAMETER   :: CM3PERM3  = 1.e+6_fp
!
! !LOCAL VARIABLES:
!
  ! Scalars
  INTEGER               :: id_CH4
  REAL(fp)              :: TROPOCH4

CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: emissch4
!
! !DESCRIPTION: Subroutine EMISSCH4 places emissions of CH4 [kg] into the
!  chemical species array.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE EMISSCH4( Input_Opt, State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_GetDiagn
    USE ErrCode_Mod
    USE Input_Opt_Mod,        ONLY : OptInput
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Met_Mod,        ONLY : MetState
    USE State_Grid_Mod,       ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid  ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
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
!
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f"
!        by Bob Yantosca. (bmy, 1/16/01)
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
    CHARACTER(LEN= 63) :: DgnName
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    ! For fields from Input_Opt
    LOGICAL            :: ITS_A_CH4_SIM
    LOGICAL            :: prtDebug
    LOGICAL, SAVE      :: FIRST = .TRUE.

    ! Arrays of state vector elements for applying emissions perturbations
    REAL(fp)           :: STATE_VECTOR(State_Grid%NX,State_Grid%NY)

    ! Array of scale factors for emissions (from HEMCO)
    REAL(fp)           :: EMIS_SF(State_Grid%NX,State_Grid%NY)

    ! Pointers
    REAL(f4), POINTER  :: Ptr2D(:,:)

    !=================================================================
    ! EMISSCH4 begins here!
    !=================================================================

    ! Nullify pointers
    Ptr2D => NULL()

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at EMISSCH4 (in GeosCore/global_ch4_mod.F90)'

    ! Copy values from Input_Opt
    ITS_A_CH4_SIM  = Input_Opt%ITS_A_CH4_SIM

    ! Do we have to print debug output?
    prtDebug = ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    IF ( ITS_A_CH4_SIM .and. prtDebug ) THEN
       print*,'BEGIN SUBROUTINE: EMISSCH4'
    ENDIF

    ! =================================================================
    ! Get fields for CH4 analytical inversions if needed
    ! =================================================================
    IF ( Input_Opt%AnalyticalInv ) THEN

       ! Evaluate the state vector field from HEMCO
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'CH4_STATE_VECTOR', &
                         STATE_VECTOR, RC)
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'CH4_STATE_VECTOR not found in HEMCO data list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    IF ( Input_Opt%UseEmisSF ) THEN

       ! Evaluate CH4 emissions scale factors from HEMCO
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'EMIS_SF', EMIS_SF, RC)
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'EMIS_SF not found in HEMCO data list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    ! =================================================================
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
    !
    !                                              (ckeller, 9/12/2013)
    ! =================================================================
    State_Chm%CH4_EMIS(:,:,:) = 0e+0_fp

    !-------------------
    ! Oil
    !-------------------
    DgnName = 'CH4_OIL'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

    ! Trap potential errors and assign HEMCO pointer to array
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSEIF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
       ErrMsg = 'Unassociated pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
    ELSE
       State_Chm%CH4_EMIS(:,:,2) =  Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()

    !-------------------
    ! Gas
    !-------------------
    DgnName = 'CH4_GAS'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

    ! Trap potential errors and assign HEMCO pointer to array
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSEIF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
       ErrMsg = 'Unassociated pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
    ELSE
       State_Chm%CH4_EMIS(:,:,3) =  Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()

    !-------------------
    ! Coal
    !-------------------
    DgnName = 'CH4_COAL'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

    ! Trap potential errors and assign HEMCO pointer to array
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSEIF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
       ErrMsg = 'Unassociated pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
    ELSE
       State_Chm%CH4_EMIS(:,:,4) =  Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()

    !-------------------
    ! Livestock
    !-------------------
    DgnName = 'CH4_LIVESTOCK'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

    ! Trap potential errors and assign HEMCO pointer to array
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSEIF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
       ErrMsg = 'Unassociated pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
    ELSE
       State_Chm%CH4_EMIS(:,:,5) =  Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()

    !-------------------
    ! Landfills
    !-------------------
    DgnName = 'CH4_LANDFILLS'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

    ! Trap potential errors and assign HEMCO pointer to array
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSEIF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
       ErrMsg = 'Unassociated pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
    ELSE
       State_Chm%CH4_EMIS(:,:,6) =  Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()

    !-------------------
    ! Wastewater
    !-------------------
    DgnName = 'CH4_WASTEWATER'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

    ! Trap potential errors and assign HEMCO pointer to array
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSEIF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
       ErrMsg = 'Unassociated pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
    ELSE
       State_Chm%CH4_EMIS(:,:,7) =  Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()

    !-------------------
    ! Rice
    !-------------------
    DgnName = 'CH4_RICE'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

    ! Trap potential errors and assign HEMCO pointer to array
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSEIF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
       ErrMsg = 'Unassociated pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
    ELSE
       State_Chm%CH4_EMIS(:,:,8) =  Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()

    !-------------------
    ! Other anthropogenic
    !-------------------
    DgnName = 'CH4_ANTHROTHER'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

    ! Trap potential errors and assign HEMCO pointer to array
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSEIF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
       ErrMsg = 'Unassociated pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Warning( ErrMsg, RC, ThisLoc = ThisLoc )
    ELSE
       State_Chm%CH4_EMIS(:,:,9) =  Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()

    !-------------------
    ! Biomass burning
    !-------------------
    DgnName = 'CH4_BIOMASS'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

    ! Trap potential errors and assign HEMCO pointer to array
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSEIF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
       ErrMsg = 'Unassociated pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
    ELSE
       State_Chm%CH4_EMIS(:,:,10) =  Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()

    !-------------------
    ! Wetland
    !-------------------
    DgnName = 'CH4_WETLAND'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

    ! Trap potential errors and assign HEMCO pointer to array
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSEIF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
       ErrMsg = 'Unassociated pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
    ELSE
       State_Chm%CH4_EMIS(:,:,11) =  Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()

    !-------------------
    ! Global seeps
    !-------------------
    DgnName = 'CH4_SEEPS'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

    ! Trap potential errors and assign HEMCO pointer to array
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSEIF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
       ErrMsg = 'Unassociated pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
    ELSE
       State_Chm%CH4_EMIS(:,:,12) =  Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()

    !-------------------
    ! Lakes
    !-------------------
    DgnName = 'CH4_LAKES'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

    ! Trap potential errors and assign HEMCO pointer to array
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSEIF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
       ErrMsg = 'Unassociated pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
    ELSE
       State_Chm%CH4_EMIS(:,:,13) =  Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()

    !-------------------
    ! Termites
    !-------------------
    DgnName = 'CH4_TERMITES'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

    ! Trap potential errors and assign HEMCO pointer to array
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSEIF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
       ErrMsg = 'Unassociated pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
    ELSE
       State_Chm%CH4_EMIS(:,:,14) =  Ptr2D(:,:)
    ENDIF
    Ptr2D => NULL()

    !-------------------
    ! Soil absorption (those are negative!)
    !-------------------
    DgnName = 'CH4_SOILABSORB'
    CALL HCO_GC_GetDiagn( Input_Opt, State_Grid, DgnName, .FALSE., RC, Ptr2D=Ptr2D )

    ! Trap potential errors and assign HEMCO pointer to array
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Cannot get pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ELSEIF ( .NOT. ASSOCIATED(Ptr2D) ) THEN
       ErrMsg = 'Unassociated pointer to HEMCO field ' // TRIM(DgnName)
       CALL GC_Warning( ErrMsg, RC, ThisLoc=ThisLoc )
    ELSE
       State_Chm%CH4_EMIS(:,:,15) =  Ptr2D(:,:) * -1.0_fp
    ENDIF
    Ptr2D => NULL()

    ! =================================================================
    ! Total emission: sum of all emissions - (2*soil absorption)
    ! We have to substract soil absorption twice because it is added
    ! to other emissions in the SUM function. (ccc, 7/23/09)
    ! =================================================================
    State_Chm%CH4_EMIS(:,:,1) = SUM(State_Chm%CH4_EMIS, 3) - (2 * State_Chm%CH4_EMIS(:,:,15))

    IF ( prtDebug ) THEN
       WRITE(*,*) 'CH4_EMIS (kg/m2/s):'
       WRITE(*,*) 'Total        : ', SUM(State_Chm%CH4_EMIS(:,:,1))
       WRITE(*,*) 'Oil          : ', SUM(State_Chm%CH4_EMIS(:,:,2))
       WRITE(*,*) 'Gas          : ', SUM(State_Chm%CH4_EMIS(:,:,3))
       WRITE(*,*) 'Coal         : ', SUM(State_Chm%CH4_EMIS(:,:,4))
       WRITE(*,*) 'Livestock    : ', SUM(State_Chm%CH4_EMIS(:,:,5))
       WRITE(*,*) 'Landfills    : ', SUM(State_Chm%CH4_EMIS(:,:,6))
       WRITE(*,*) 'Wastewater   : ', SUM(State_Chm%CH4_EMIS(:,:,7))
       WRITE(*,*) 'Rice         : ', SUM(State_Chm%CH4_EMIS(:,:,8))
       WRITE(*,*) 'Other anth   : ', SUM(State_Chm%CH4_EMIS(:,:,9))
       WRITE(*,*) 'Biomass burn : ', SUM(State_Chm%CH4_EMIS(:,:,10))
       WRITE(*,*) 'Wetlands     : ', SUM(State_Chm%CH4_EMIS(:,:,11))
       WRITE(*,*) 'Seeps        : ', SUM(State_Chm%CH4_EMIS(:,:,12))
       WRITE(*,*) 'Lakes        : ', SUM(State_Chm%CH4_EMIS(:,:,13))
       WRITE(*,*) 'Termites     : ', SUM(State_Chm%CH4_EMIS(:,:,14))
       WRITE(*,*) 'Soil absorb  : ', SUM(State_Chm%CH4_EMIS(:,:,15))
    ENDIF

    ! =================================================================
    ! Do scaling for analytical inversion
    ! =================================================================
    IF ( Input_Opt%AnalyticalInv  .or. &
         Input_Opt%UseEmisSF      .or. &
         Input_Opt%UseOHSF        ) THEN

       ! Don't optimize for soil absorption so remove from the total
       ! emissions array
       State_Chm%CH4_EMIS(:,:,1) = State_Chm%CH4_EMIS(:,:,1) + State_Chm%CH4_EMIS(:,:,15)

       ! Rescale emissions
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J)	 
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX

          !------------------------------------------------------------
          ! For applying scale factors from analytical inversion
          !------------------------------------------------------------
          IF ( Input_Opt%UseEmisSF ) THEN
             ! Scale total emissions
             State_Chm%CH4_EMIS(I,J,1) = State_Chm%CH4_EMIS(I,J,1) * EMIS_SF(I,J)
          ENDIF

          !------------------------------------------------------------
          ! Perturb emissions for analytical inversion
          !------------------------------------------------------------
          IF ( Input_Opt%AnalyticalInv ) THEN

             ! Only apply emission perturbation to current state vector
             ! element number
             IF ( Input_Opt%StateVectorElement .GT. 0 ) THEN
                IF ( STATE_VECTOR(I,J) == Input_Opt%StateVectorElement ) THEN
                   State_Chm%CH4_EMIS(I,J,1) = State_Chm%CH4_EMIS(I,J,1) * Input_Opt%PerturbEmis
                   !Print*, 'Analytical Inversion: Scaled state vector element ', &
                   !        Input_Opt%StateVectorElement, ' by ', &
                   !        Input_Opt%PerturbEmis
                ENDIF
             ENDIF
          ENDIF

       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       ! Now that we've done the emission factor scaling, add soil absorption
       ! back to the total emissions array
       State_Chm%CH4_EMIS(:,:,1) = State_Chm%CH4_EMIS(:,:,1) - State_Chm%CH4_EMIS(:,:, 15)

    ENDIF

    IF ( ITS_A_CH4_SIM .and. prtDebug ) THEN
       print*,'END SUBROUTINE: EMISSCH4'
    ENDIF

  END SUBROUTINE EMISSCH4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: chemch4
!
! !DESCRIPTION: Subroutine CHEMCH4 computes the chemical loss of CH4
!  (sources - sinks). (jsw, bnd, bmy, 6/8/00, 10/3/05)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CHEMCH4( Input_Opt,  State_Chm, State_Diag, &
                      State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
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
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (6/8/00).  Inserted into module "global_ch4_mod.f"
!        by Bob Yantosca. (bmy, 1/16/01)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    LOGICAL, SAVE      :: FIRSTCHEM = .TRUE.
    INTEGER            :: I, J, L
    REAL(fp)           :: PREVCH4(State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    ! Number of days per month
    INTEGER            :: NODAYS(12) = (/ 31, 28, 31, 30, 31, 30, &
                                          31, 31, 30, 31, 30, 31 /)

    ! For fields from Input_Opt
    LOGICAL            :: LSPLIT
    LOGICAL            :: prtDebug

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !=================================================================
    ! CHEMCH4 begins here!
    !=================================================================

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at CHEMCH4 (in module GeosCore/global_ch4_mod.F90)'

    ! Copy values from Input_Opt
    LSPLIT  = Input_Opt%LSPLIT
    prtDebug= ( Input_Opt%LPRT .and. Input_Opt%amIRoot )

    ! Point to the chemical species
    Spc     => State_Chm%Species

    IF ( prtDebug ) THEN
       WRITE( 6, '(a)' ) '% --- ENTERING CHEMCH4! ---'
    ENDIF

    !================================================================
    ! Evaluate OH and Cl fields from HEMCO. Doing this every call
    ! allows usage of HEMCO scaling and masking features.
    !================================================================

    ! Evalulate the global OH from HEMCO
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'GLOBAL_OH', State_Chm%BOH, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'GLOBAL_OH not found in HEMCO data list!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Evalulate the global Cl from HEMCO
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'GLOBAL_Cl', State_Chm%BCl, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'GLOBAL_Cl not found in HEMCO data list!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! HISTORY (aka netCDF diagnostics)
    ! OH concentration in [molec/cm3] after chemistry
    !
    ! BOH from HEMCO is in kg/m3, convert to molec/cm3
    !=================================================================
    IF ( State_Diag%Archive_OHconcAfterChem ) THEN
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          IF ( State_Met%InChemGrid(I,J,L) ) THEN
             State_Diag%OHconcAfterChem(I,J,L) = &
                  ( State_Chm%BOH(I,J,L) * XNUMOL_OH / CM3PERM3 )
          ELSE
             State_Diag%OHconcAfterChem(I,J,L) = 0.0_f4
          ENDIF
       ENDDO
       ENDDO
       ENDDO
    ENDIF

    !=================================================================
    ! HISTORY (aka netCDF diagnostics)
    ! Archive quantities for computing CH4 metrics such as global
    ! mean OH, MCF lifetime, and CH4 lifetimes.
    !=================================================================
    CALL CH4_Metrics( Input_Opt,  State_Chm, State_Diag,                     &
                      State_Grid, State_Met, RC                             )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in routine "CH4_Metrics!"'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=================================================================
    ! If multi-CH4 species, we store the CH4 total conc. to
    ! distribute the sink after the chemistry. (ccc, 2/10/09)
    !=================================================================
    IF ( LSPLIT ) THEN

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          PREVCH4(I,J,L) = Spc(1)%Conc(I,J,L)
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

    !=================================================================
    ! Calculate rate of decay of CH4 by OH oxidation.
    !=================================================================
    CALL CH4_DECAY( Input_Opt,  State_Chm, State_Diag, &
                    State_Grid, State_Met, RC )

    !=================================================================
    ! Calculate CH4 chemistry in layers above tropopause
    !=================================================================
    CALL CH4_STRAT( Input_Opt,  State_Chm, State_Diag, &
                    State_Grid, State_Met, RC )

    !=================================================================
    ! Distribute the chemistry sink from total CH4 to other CH4
    !     species. (ccc, 2/10/09)
    !=================================================================
    IF ( LSPLIT ) THEN
       CALL CH4_DISTRIB( Input_Opt, State_Chm, State_Grid, PREVCH4 )
    ENDIF

    ! Free pointer
    Spc => NULL()


    ! Set FIRSTCHEM to FALSE
    FIRSTCHEM = .FALSE.

  END SUBROUTINE CHEMCH4
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_decay
!
! !DESCRIPTION: Subroutine CH4\_DECAY calculates the decay rate of CH4 by OH.
!  OH is the only sink for CH4 considered here. (jsw, bnd, bmy, 1/16/01,
!  7/20/04)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CH4_DECAY( Input_Opt,  State_Chm, State_Diag,  &
                        State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Species_Mod,          ONLY : SpcConc
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE TIME_MOD,             ONLY : GET_TS_CHEM
    USE TIME_MOD,             ONLY : GET_MONTH
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
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
!  Monthly loss of CH4 is summed in TCH4(3)
!     TCH4(3)  = CH4 sink by OH
!
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f"
!        by Bob Yantosca. (bmy, 1/16/01)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER            :: I,  J,    L
    REAL(fp)           :: DT, GCH4, Spc2GCH4
    REAL(fp)           :: KRATE, C_OH
    REAL(fp)           :: KRATE_Cl, C_Cl
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc
    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    ! Array of scale factors for OH (from HEMCO)
    REAL(fp)           :: OH_SF(State_Grid%NX,State_Grid%NY)

    !=================================================================
    ! CH4_DECAY begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Chemistry timestep in seconds
    DT = GET_TS_CHEM()

    ! Point to the chemical species array
    Spc => State_Chm%Species

    ! =================================================================
    ! Get fields for CH4 analytical inversions if needed
    ! =================================================================
    IF ( Input_Opt%UseOHSF ) THEN

       ! Evaluate OH scale factors from HEMCO
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'OH_SF', OH_SF, RC)
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'OH_SF not found in HEMCO data list!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

    ENDIF

    !=================================================================
    ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
    !
    ! Zero the relevant diagnostic fields of State_Diag because the
    ! position of the tropopause changes from one timestep to the next
    !=================================================================
    IF ( State_Diag%Archive_LossCH4byClinTrop ) THEN
       State_Diag%LossCH4byClinTrop = 0.0_f4
    ENDIF

    IF ( State_Diag%Archive_LossCH4byOHinTrop ) THEN
       State_Diag%LossCH4byOHinTrop = 0.0_f4
    ENDIF

    !=================================================================
    ! Compute decay of CH4 by OH and Cl in the troposphere
    !
    ! The decay for CH4 is calculated by:
    !    OH + CH4 -> CH3 + H2O
    !    k = 2.45E-12 exp(-1775/T)
    !
    !    This is from JPL '97.
    !    JPL '00, '06, & '11 do not revise '97 value. (jsw, kjw, ajt)
    !
    ! The decay for CH4 by Cl is calculated by:
    !    Cl + CH4 -> HCl + CH3
    !    k = 9.6E-12 exp(-1360/T)
    !
    !    This is from Kirschke et al., Nat. Geosci., 2013.
    !=================================================================

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( L, J, I, KRATE, Spc2GCH4, GCH4, C_OH ) &
    !$OMP PRIVATE( C_Cl, KRATE_Cl ) &
    !$OMP REDUCTION( +:TROPOCH4 )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Only consider tropospheric boxes
       IF ( State_Met%InChemGrid(I,J,L) ) THEN

          ! Calculate rate coefficients
          KRATE    = 2.45e-12_fp * EXP( -1775e+0_fp / State_Met%T(I,J,L))
          KRATE_Cl = 9.60e-12_fp * EXP( -1360e+0_fp / State_Met%T(I,J,L))

          ! Conversion from [kg/box] --> [molec/cm3]
          ! [kg CH4/box] * [box/cm3] * XNUMOL_CH4 [molec CH4/kg CH4]
          Spc2GCH4 = 1e+0_fp / State_Met%AIRVOL(I,J,L) / 1e+6_fp * XNUMOL_CH4

          ! CH4 in [molec/cm3]
          GCH4 = Spc(1)%Conc(I,J,L) * Spc2GCH4

          ! OH in [molec/cm3]
          ! BOH from HEMCO in units of kg/m3, convert to molec/cm3
          C_OH = State_Chm%BOH(I,J,L) * XNUMOL_OH / CM3PERM3

          ! Apply scale factors from analytical inversion
          IF ( Input_Opt%UseOHSF ) THEN
             C_OH = C_OH * OH_SF(I,J)
             !Print*, 'Applying scale factor to OH: ', OH_SF(I,J)
          ENDIF

          ! Cl in [molec/cm3]
          ! BCl from HEMCO in units of mol/mol, convert to molec/cm3
          C_Cl = State_Chm%BCl(I,J,L) * State_Met%AIRNUMDEN(I,J,L)

          TROPOCH4 = TROPOCH4 + GCH4 * KRATE    * C_OH * DT / Spc2GCH4 &
                              + GCH4 * KRATE_Cl * C_Cl * DT / Spc2GCH4

          !-----------------------------------------------------------
          ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
          !
          ! Archive Loss of CH4 (kg/s) reactions with OH and Cl
          !-----------------------------------------------------------

          ! Loss CH4 by reaction with Cl [kg/s]
          IF ( State_Diag%Archive_LossCH4byClinTrop ) THEN
             State_Diag%LossCH4byClinTrop(I,J,L) = &
                  ( GCH4 * KRATE_Cl * C_Cl ) / Spc2GCH4
          ENDIF

          IF ( State_Diag%Archive_LossCH4byOHinTrop ) THEN
             State_Diag%LossCH4byOHinTrop(I,J,L) = &
                  ( GCH4 * KRATE * C_OH ) / Spc2GCH4
          ENDIF

          ! Calculate new CH4 value: [CH4]=[CH4](1-k[OH]*delt)
          GCH4 = GCH4 * ( 1e+0_fp - KRATE    * C_OH * DT )
          GCH4 = GCH4 * ( 1e+0_fp - KRATE_Cl * C_Cl * DT )

          ! Convert back from [molec/cm3] --> [kg/box]
          Spc(1)%Conc(I,J,L) = GCH4 / Spc2GCH4

       ENDIF
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !print*,'% --- CHEMCH4: CH4_DECAY: TROP DECAY (Tg): ',TROPOCH4/1e9
    !print*,'Trop decay should be over 1Tg per day globally'
    !print*,'    ~ 500Tg/365d ~ 1.37/d'

    ! Free pointers
    Spc => NULL()

  END SUBROUTINE CH4_DECAY
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_metrics
!
! !DESCRIPTION: Computes mass-weighted mean OH columns (full-atmosphere and
!  trop-only) that are needed to compute the overall mean OH concentration.
!  This is used as a metric as to how reactive, or "hot" the chemistry
!  mechanism is.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CH4_Metrics( Input_Opt,  State_Chm, State_Diag, &
                          State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Chm_Mod,  ONLY : ChmState
    USE PhysConstants,  ONLY : AVO
    USE PhysConstants,  ONLY : XNUMOLAIR
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt    ! Input Options object
    TYPE(ChmState), INTENT(IN)    :: State_Chm    ! Chemistry State object
    TYPE(GrdState), INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met    ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!
! !REMARKS:
!  References:
!  (1) Prather, M. and C. Spivakovsky, "Tropospheric OH and
!       the lifetimes of hydrochlorofluorocarbons", JGR,
!       Vol 95, No. D11, 18723-18729, 1990.
!  (2) Lawrence, M.G, Joeckel, P, and von Kuhlmann, R., "What
!       does the global mean OH concentraton tell us?",
!       Atm. Chem. Phys, 1, 37-49, 2001.
!  (3) WMO/UNEP Scientific Assessment of Ozone Depletion: 2010
!
! !REVISION HISTORY:
!  18 Aug 2020 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS
!
    REAL(f8), PARAMETER :: M3toCM3        = 1.0e+6_f8
    REAL(f8), PARAMETER :: MCM3toKGM3_OH  = M3toCM3 * 17.01e-3_f8 / AVO
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL,  SAVE      :: first          = .TRUE.
    INTEGER,  SAVE      :: id_CH4         = -1
    REAL(f8), SAVE      :: MCM3toKGM3_CH4 = -1.0_f8

    ! Scalars
    INTEGER             :: I,           J,            L
    REAL(f8)            :: airMass_m,   airmass_kg,   airMassFull
    REAL(f8)            :: airMassTrop, CH4conc_kgm3, CH4conc_mcm3
    REAL(f8)            :: CH4mass_kg,  CH4mass_m,    CH4massFull
    REAL(f8)            :: CH4massTrop, OHconc_kgm3,  OHconc_mcm3
    REAL(f8)            :: OHmassWgt,   OHmassFull,   OHmassTrop
    REAL(f8)            :: Ktrop,       LossOHbyCH4,  LossOHbyMCF
    REAL(f8)            :: volume

    ! Strings
    CHARACTER(LEN=255)  :: errMsg,      thisLoc

    !========================================================================
    ! Compute_Mean_OH_and_CH4 begins here!
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    errMsg  = ''
    thisLoc = ' -> at Compute_Mean_OH (in module GeosCore/diagnostics_mod.F90)'

    ! Exit if we have not turned on the Metrics collection
    IF ( .not. State_Diag%Archive_Metrics ) RETURN

    !========================================================================
    ! First-time setup
    !========================================================================
    IF ( first ) THEN

       ! Get the species ID for CH4
       id_CH4 = Ind_('CH4')
       IF ( id_CH4 < 0 ) THEN
          errMsg = 'CH4 is not a defined species in this simulation!!!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF

       ! Convert [molec CH4 cm-3] --> [kg CH4 m-3]
       MCM3toKGM3_CH4 = M3toCM3                                              &
                      * ( State_Chm%SpcData(id_CH4)%Info%MW_g * 1.0e-3_f8 )  &
                      / AVO

       ! Reset first-time flag
       first  = .FALSE.
    ENDIF

    !========================================================================
    ! Loop over surface boxes and compute mean OH in columns
    !========================================================================
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I,            J,           L,           airMass_kg       )&
    !$OMP PRIVATE( airMass_m,    airMassFull, airMassTrop, CH4conc_kgm3     )&
    !$OMP PRIVATE( CH4conc_mcm3, CH4mass_kg,  CH4massFull, CH4massTrop      )&
    !$OMP PRIVATE( Ktrop,        LossOHbyCH4, LossOHbyMCF, OHconc_kgm3      )&
    !$OMP PRIVATE( OHconc_mcm3,  OHmassWgt,   OHmassFull,  OHmassTrop       )&
    !$OMP PRIVATE( volume                                                   )&
    !$OMP SCHEDULE( DYNAMIC, 4                                              )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       !--------------------------------------------------------------------
       ! Zero column-specific quantities
       !--------------------------------------------------------------------
       airMass_kg   = 0.0_f8
       airMass_m    = 0.0_f8
       airMassFull  = 0.0_f8
       airMassTrop  = 0.0_f8
       CH4conc_kgm3 = 0.0_f8
       CH4conc_mcm3 = 0.0_f8
       CH4mass_kg   = 0.0_f8
       CH4massFull  = 0.0_f8
       CH4massTrop  = 0.0_f8
       Ktrop        = 0.0_f8
       LossOHbyCH4  = 0.0_f8
       LossOHbyMCF  = 0.0_f8
       OHconc_kgm3  = 0.0_f8
       OHconc_mcm3  = 0.0_f8
       OHmassWgt    = 0.0_f8
       OHmassFull   = 0.0_f8
       OHmassTrop   = 0.0_f8
       volume       = 0.0_f8

       !--------------------------------------------------------------------
       ! Loop over the number of levels in the chemistry grid
       ! (which for CH4 simulations is within the troposphere)
       !--------------------------------------------------------------------
       DO L = 1, State_Grid%NZ

          ! Compute box volume [cm3], and air mass ([molec] and [kg])
          ! Note: air mass in [molec] is also the atmospheric burden of
          ! methyl chloroform (aka MCF, formula=CH3CCl3), since we assume
          ! a uniform mixing ratio (=1) of MCF in air.
          volume       = State_Met%AIRVOL(I,J,L)    * M3toCM3
          airMass_m    = State_Met%AIRNUMDEN(I,J,L) * volume
          airMass_kg   = airMass_m / XNUMOLAIR

          ! CH4 mass [kg]
          CH4mass_kg   = State_Chm%Species(id_CH4)%Conc(I,J,L)

          ! CH4 concentration [kg m-3] and [molec cm-3]
          CH4conc_kgm3 = CH4mass_kg   / volume
          CH4conc_mcm3 = CH4conc_kgm3 / MCM3toKGM3_CH4

          ! OH concentration [kg m-3] and [molec cm-3]
          OHconc_kgm3  = State_Chm%BOH(I,J,L)
          OHconc_mcm3  = OHconc_kgm3 /  MCM3toKGM3_OH

          ! Airmass-weighted OH [kg air * (kg OH  m-3)]
          OHmassWgt    = airmass_kg * OHconc_kgm3

          ! Sum the air mass, mass-weighted CH4,
          ! and mass-weighted OH in the full-atm column
          airMassFull  = airMassFull + airMass_kg
          CH4massFull  = CH4MassFull + CH4mass_kg
          OHmassFull   = OHmassFull  + OHmassWgt

          !------------------------------------------------------------------
          ! Only do the following for tropospheric boxes
          !------------------------------------------------------------------
          IF ( State_Met%InTroposphere(I,J,L) ) THEN

             ! Sum the air mass, mass-weighted CH4,
             ! and mass-weighted OH in the trop-only column
             airMassTrop = airMassTrop + airMass_kg
             CH4massTrop = CH4MassTrop + CH4mass_kg
             OHmassTrop  = OHmassTrop  + OHmassWgt

             ! Compute CH4 + OH loss rate in troposphere
             ! Ktrop (Arrhenius parameter) has units [cm3/molec/s]
             ! OHconc has units [molec/cm3]
             ! AirMass has units [molec]
             ! Resultant units of CH4 loss rate = [molec/s]
             Ktrop = 2.45e-12_f8 * EXP( -1775.0_f8 / State_Met%T(I,J,L) )
             LossOHbyCH4 = LossOHbyCH4 + ( Ktrop * OHconc_MCM3 * airMass_m )

             ! Compute MCF + OH loss rate in the troposphere
             ! Ktrop (Arrhenius parameter) has units [cm3/molec/s]
             ! OHconc has units [molec/cm3]
             ! AirMass has units [molec]
             ! Resultant units of MCF loss rate = [molec/s]
             Ktrop = 1.64e-12_f8 * EXP( -1520.0_f8 / State_Met%T(I,J,L) )
             LossOHbyMCF = LossOHbyMCF + ( Ktrop * OHconc_MCM3 * airMass_m )

             !---------------------------------------------------------------
             ! HISTORY (aka netCDF diagnostics)
             !
             ! Keep track of CH4 emisisons [kg/s] for computing
             ! the various lifetime metrics in post-processing
             !---------------------------------------------------------------
             IF ( L == 1 .and. State_Diag%Archive_CH4emission ) THEN
                State_Diag%CH4emission(I,J) = State_Chm%CH4_EMIS(I,J,id_CH4)           &
                                            * State_Grid%Area_M2(I,J)
             ENDIF
          ENDIF
       ENDDO

       !---------------------------------------------------------------------
       ! HISTORY (aka netCDF diagnostics)
       ! Air mass [kg]
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_AirMassColumnFull ) THEN
          State_Diag%AirMassColumnFull(I,J) = airMassFull
       ENDIF

       IF ( State_Diag%Archive_AirMassColumnTrop ) THEN
          State_Diag%AirMassColumnTrop(I,J) = airMassTrop
       ENDIF

       !---------------------------------------------------------------------
       ! HISTORY (aka netCDF diagnostics)
       ! CH4 mass [kg], full-atmosphere and trop-only column sums
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_CH4massColumnFull ) THEN
          State_Diag%CH4massColumnFull(I,J) = CH4massFull
       ENDIF

       IF ( State_Diag%Archive_CH4massColumnTrop ) THEN
          State_Diag%CH4massColumnTrop(I,J) = CH4massTrop
       ENDIF

       !---------------------------------------------------------------------
       ! HISTORY (aka netCDF diagnostics)
       ! Mass-weighted mean OH [kg air * (kg OH m-3)]
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_OHwgtByAirMassColumnFull ) THEN
          State_Diag%OHwgtByAirMassColumnFull(I,J) = OHmassFull
       ENDIF

       IF ( State_Diag%Archive_OHwgtByAirMassColumnTrop ) THEN
          State_Diag%OHwgtByAirMassColumnTrop(I,J) = OHmassTrop
       ENDIF

       !-----------------------------------------------------------------
       ! HISTORY (aka netCDF diagnostics)
       !
       ! OH loss by CH4 + OH loss in troposphere [molec/s] and
       ! OH loss by MCF + OH loss in troposphere [molec/s]
       !----------------------------------------------------------------
       IF ( State_Diag%Archive_LossOHbyCH4columnTrop ) THEN
          State_Diag%LossOHbyCH4columnTrop(I,J) = LossOHbyCH4
       ENDIF

       IF ( State_Diag%Archive_LossOHbyMCFcolumnTrop ) THEN
          State_Diag%LossOHByMCFcolumnTrop(I,J) = LossOHbyMCF
       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE CH4_Metrics
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_strat
!
! !DESCRIPTION: Subroutine CH4\_STRAT calculates uses production rates for CH4
!  to  calculate loss of CH4 in above the tropopause. (jsw, bnd, bmy, 1/16/01,
!  7/20/04). This is only done if unified chemistry is not active.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CH4_STRAT( Input_Opt,  State_Chm, State_Diag, &
                        State_Grid, State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE Input_Opt_Mod,        ONLY : OptInput
    USE Species_Mod,          ONLY : SpcConc
    USE State_Chm_Mod,        ONLY : ChmState
    USE State_Diag_Mod,       ONLY : DgnState
    USE State_Grid_Mod,       ONLY : GrdState
    USE State_Met_Mod,        ONLY : MetState
    USE TIME_MOD,             ONLY : GET_TS_CHEM
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input options
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
!  Production (mixing ratio/sec) rate provided by Dylan Jones.
!  Only production by CH4 + OH is considered.
!
! !REVISION HISTORY:
!  (1 ) Created by Bryan Duncan (1/99).  Adapted for CH4 chemistry by
!        James Wang (7/00).  Inserted into module "global_ch4_mod.f"
!        by Bob Yantosca. (bmy, 1/16/01)
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER           :: I,  J,    L
    REAL(fp)          :: DT, GCH4, Spc2GCH4, LRATE

    ! Strings
    CHARACTER(LEN=255)    :: ThisLoc
    CHARACTER(LEN=255)    :: ErrMsg

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    ! Array for monthly average CH4 loss freq [1/s] (from HEMCO)
    REAL(fp)          :: CH4LOSS(State_Grid%NX,State_Grid%NY,State_Grid%NZ)

    !=================================================================
    ! CH4_STRAT begins here!
    !=================================================================

    ! Initialize
    RC    =  GC_SUCCESS
    ErrMsg      = ''
    ThisLoc     = ' -> at CH4_STRAT (in module GeosCore/global_ch4_mod.F90)'

    ! Point to chemical species
    Spc   => State_Chm%Species

    ! Evalulate CH4 loss frequency from HEMCO. This must be done
    ! every timestep to allow masking or scaling in HEMCO config.
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'CH4_LOSS', CH4LOSS, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'CH4_LOSS not found in HEMCO data list!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !==============================================================
    ! %%%%% HISTORY (aka netCDF diagnostics) %%%%%
    !
    ! Zero the relevant diagnostic fields of State_Diag because
    ! the position of the tropopause changes from one timestep
    ! to the next
    !==============================================================
    IF ( State_Diag%Archive_LossCH4inStrat ) THEN
       State_Diag%LossCH4inStrat = 0.0_f4
    ENDIF

    ! Chemistry timestep [s]
    DT  = GET_TS_CHEM()

    !=================================================================
    ! Loop over stratospheric boxes only
    !=================================================================
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( I, J, L, Spc2GCH4, GCH4, LRATE )
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Only proceed if we are outside of the chemistry grid
       IF ( .not. State_Met%InChemGrid(I,J,L) ) THEN

          ! Conversion factor [kg/box] --> [molec/cm3]
          ! [kg/box] / [AIRVOL * 1e6 cm3] * [XNUMOL_CH4 molec/mole]
          Spc2GCH4 = 1e+0_fp / State_Met%AIRVOL(I,J,L) / 1e+6_fp * XNUMOL_CH4

          ! CH4 in [molec/cm3]
          GCH4 = Spc(1)%Conc(I,J,L) * Spc2GCH4

          ! Loss rate [molec/cm3/s]
          LRATE = GCH4 * CH4LOSS( I,J,L )

          ! Update Methane concentration in this grid box [molec/cm3]
          GCH4 = GCH4 - ( LRATE * DT )

          ! Convert back from [molec CH4/cm3] --> [kg/box]
          Spc(1)%Conc(I,J,L) = GCH4 / Spc2GCH4

          !------------------------------------------------------------
          ! %%%%%% HISTORY (aka netCDF diagnostics) %%%%%
          !
          ! Loss of CH4 by OH above tropopause [kg/s]
          !------------------------------------------------------------
          IF ( State_Diag%Archive_LossCH4inStrat ) THEN
             State_Diag%LossCH4inStrat(I,J,L) = LRATE / Spc2GCH4
          ENDIF

       ENDIF
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE CH4_STRAT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ch4_distrib
!
! !DESCRIPTION: Subroutine CH4\_DISTRIB allocates the chemistry sink to
!  different emission species. (ccc, 10/2/09)
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE CH4_DISTRIB( Input_Opt, State_Chm, State_Grid, PREVCH4 )
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
    REAL(fp),       INTENT(IN)    :: PREVCH4(State_Grid%NX, & ! CH4 before chem
                                             State_Grid%NY, &
                                             State_Grid%NZ)
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
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
    INTEGER           :: I, J, L, N, NA, nAdvect

    ! Pointers
    TYPE(SpcConc), POINTER :: Spc(:)

    !========================================================================
    ! CH4_DISTRIB begins here
    !========================================================================

    ! Point to chemical species array [kg]
    Spc => State_Chm%Species

    ! fix nAdvect (Xueying Yu, 12/10/2017)
    nAdvect = State_Chm%nAdvect

    ! Loop over the number of advected species
    DO NA = 2, nAdvect

       ! Advected species ID
       N = State_Chm%Map_Advect(NA)

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( I, J, L )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          Spc(N)%Conc(I,J,L) = &
                         SAFE_DIV(Spc(N)%Conc(I,J,L),PREVCH4(I,J,L),0.e+0_fp) &
                         * Spc(1)%Conc(I,J,L)
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDDO

    ! Free pointer
    Spc => NULL()

  END SUBROUTINE CH4_DISTRIB
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_global_ch4
!
! !DESCRIPTION: Subroutine INIT\_GLOBAL\_CH4 allocates and zeroes module
!  arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE INIT_GLOBAL_CH4( Input_Opt, State_Chm, State_Diag, State_Grid, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,      ONLY : OptInput
    USE State_Chm_Mod,      ONLY : Ind_, ChmState
    USE State_Diag_Mod,     ONLY : DgnState
    USE State_Grid_Mod,     ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
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
! !REVISION HISTORY:
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Strings
    CHARACTER(LEN=255) :: ErrMsg
    CHARACTER(LEN=255) :: ThisLoc

    !=================================================================
    ! INIT_GLOBAL_CH4 begins here!
    !=================================================================

    ! Assume Success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> INIT_CH4 (in module GeosCore/global_ch4_mod.F90)'

    ! Exit immediately if it's a dry-run simulation
    IF ( Input_Opt%DryRun ) RETURN

    ! Define species ID flag
    id_CH4  = Ind_('CH4')

    ! Make sure CH4 is a defined species (bmy, 6/20/16)
    IF ( id_CH4 <= 0 ) THEN
       ErrMsg = 'CH4 is an undefined species!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Initialize tropoch4 (counts total decay of CH4 due to OH)
    TROPOCH4 = 0e+0_fp

  END SUBROUTINE INIT_GLOBAL_CH4
!EOC
END MODULE GLOBAL_CH4_MOD
