!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: diagnostics_mod.F90
!
! !DESCRIPTION: Module diagnostics\_mod.F90 contains subroutines for
!  setting State_Diag diagnostics arrays for the purposes of outputting
!  in netcdf format. Source code for setting diagnostics arrays for output
!  in binary format are not included in this module.
!
! !INTERFACE:
!
MODULE Diagnostics_mod
!
! !USES:
!
  USE ErrCode_Mod
  USE Precision_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Set_Diagnostics_EndofTimestep
  PUBLIC :: Zero_Diagnostics_StartofTimestep
  PUBLIC :: Compute_Budget_Diagnostics
#ifdef ADJOINT
  PUBLIC :: Set_SpcAdj_Diagnostic
#endif
!
! !PRIVATE MEMBER FUNCTIONS
!
  PRIVATE :: Set_SpcConc_Diags_VVDry
  PRIVATE :: Set_SpcConc_Diags_MND
!
! !REVISION HISTORY:
!  01 Feb 2018 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
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
! !IROUTINE: Set_Diagnostics_EndofTimestep
!
! !DESCRIPTION: Updates various diagnostics for History output at the end
!  of the GEOS-Chem "heartbeat" timestep.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Diagnostics_EndofTimestep( Input_Opt,  State_Chm,           &
                                            State_Diag, State_Grid,          &
                                            State_Met,  RC )
!
! !USES:
!
    USE Input_Opt_Mod,    ONLY : OptInput
    USE State_Met_Mod,    ONLY : MetState
    USE State_Chm_Mod,    ONLY : ChmState, Ind_
    USE State_Diag_Mod,   ONLY : DgnState, DgnMap
    USE State_Grid_Mod,   ONLY : GrdState
    USE PhysConstants,    ONLY : AIRMW,  AVO
    USE TIME_MOD,         ONLY : GET_LOCALTIME
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt      ! Input Options object
    TYPE(GrdState),   INTENT(IN)    :: State_Grid     ! Grid state object
    TYPE(MetState),   INTENT(IN)    :: State_Met      ! Meteorology state object
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm      ! Chemistry state obj
    TYPE(DgnState),   INTENT(INOUT) :: State_Diag     ! Diagnostics state obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,              INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  01 Feb 2018 - E. Lundgren - initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER                 :: I, J, L, N, S
    REAL(fp)                :: ToPptv, LT

    ! SAVEd scalars
    INTEGER, SAVE           :: id_Hg2 = -1
    INTEGER, SAVE           :: id_HgP = -1
    LOGICAL, SAVE           :: FIRST_Hg  = .TRUE.

    ! Strings
    CHARACTER(LEN=255)      :: ErrMsg, ThisLoc

    ! Objects
    TYPE(DgnMap), POINTER   :: mapData

    !========================================================================
    ! Set_Diagnostics_EndofTimestep begins here
    !========================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
      ' -> at Set_Diagnostics_EndofTimestep (in GeosCore/diagnostics_mod.F90)'

    !------------------------------------------------------------------------
    ! Set species concentration for diagnostics in units of
    ! v/v dry air = mol/mol dry air
    !------------------------------------------------------------------------
    CALL Set_SpcConc_Diags_VVDry( Input_Opt,  State_Chm, State_Diag,         &
                                  State_Grid, State_Met, RC                 )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered setting SpeciesConcVV diagnostic'
       CALL GC_ERROR( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Set species concentration for diagnostics in units of
    ! molec/cm3 (hplin, 11/21/21)
    !-----------------------------------------------------------------------
    CALL Set_SpcConc_Diags_MND  ( Input_Opt,  State_Chm, State_Diag,         &
                                  State_Grid, State_Met, RC                 )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered setting SpeciesConcMND diagnostic'
       CALL GC_ERROR( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

#ifdef ADJOINT
    !-----------------------------------------------------------------------
    ! Set species concentration diagnostic in units specified in state_diag_mod
    !-----------------------------------------------------------------------
    IF ( State_Diag%Archive_SpeciesAdj ) THEN
       CALL Set_SpcAdj_Diagnostic( Input_Opt,  State_Chm, State_Diag,         &
                                   State_Grid, State_Met, RC                 )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered setting SpeciesAdj diagnostic'
          CALL GC_ERROR( ErrMsg, RC, ThisLoc )
       ENDIF
    ENDIF
#endif

    !------------------------------------------------------------------------
    ! Set total dry deposition flux
    !------------------------------------------------------------------------
    IF ( State_Diag%Archive_DryDep ) THEN
       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( I, J, S                                               )&
       !$OMP COLLAPSE( 3                                                    )
       DO S = 1, State_Diag%Map_DryDep%nSlots
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          State_Diag%DryDep(I,J,S) = State_Diag%DryDepChm(I,J,S)             &
                                   + State_Diag%DryDepMix(I,J,S)
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !------------------------------------------------------------------------
    ! Set total dry deposition flux
    !------------------------------------------------------------------------
    IF ( State_Diag%Archive_SatDiagnDryDep ) THEN
       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( I, J, S                                               )&
       !$OMP COLLAPSE( 3                                                    )
       DO S = 1, State_Diag%Map_SatDiagnDryDep%nSlots
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          State_Diag%SatDiagnDryDep(I,J,S) = State_Diag%DryDepChm(I,J,S)  &
                                           + State_Diag%DryDepMix(I,J,S)
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !------------------------------------------------------------------------
    ! Compute fraction of time each grid box spent in the troposphere
    !------------------------------------------------------------------------
    IF ( State_Diag%Archive_FracOfTimeInTrop ) THEN
       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP SCHEDULE( DYNAMIC, 8                                           )&
       !$OMP PRIVATE( I, J, L                                               )&
       !$OMP COLLAPSE( 3                                                    )
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          IF ( State_Met%InTroposphere(I,J,L) ) THEN
             State_Diag%FracOfTimeInTrop(I,J,L) = 1.0_f4
          ELSE
             State_Diag%FracOfTimeInTrop(I,J,L) = 0.0_f4
          ENDIF
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !------------------------------------------------------------------------
    ! Diagnostics for the mercury and tagged mercury simulations
    !------------------------------------------------------------------------
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

       ! Get species indices for Hg2 and HgP
       IF ( FIRST_Hg ) THEN
          id_Hg2 = Ind_('Hg2')
          id_HgP = Ind_('HgP')
          FIRST_Hg  = .FALSE.
       ENDIF

       !--------------------------------------------
       ! Ractive gaseous mercury (RGM) [pptv]
       !--------------------------------------------
       IF ( id_Hg2 > 0 .and. State_Diag%Archive_ReactiveGaseousHg ) THEN

          ! Conversion factor to pptv
          ToPptv = ( AIRMW                                  /                &
                     State_Chm%SpcData(id_Hg2)%Info%MW_g  *                  &
                     1.0e+12_fp                               )

          ! Save into State_diag
          State_Diag%ReactiveGaseousHg = &
                     State_Chm%Species(id_Hg2)%Conc(:,:,:) * ToPptv
       ENDIF

       !--------------------------------------------
       ! Ractive particulate mercury (RGM) [pptv]
       !--------------------------------------------
       IF ( id_HgP > 0 .and. State_Diag%Archive_ParticulateBoundHg ) THEN

          ! Conversion factor to pptv
          ToPptv = ( AIRMW                                  /                &
                     State_Chm%SpcData(id_HgP)%Info%MW_g  *                  &
                     1.0e+12_fp                               )

          ! Save into State_Diag
          State_Diag%ParticulateBoundHg = &
                     State_Chm%Species(id_HgP)%Conc(:,:,:) * ToPptv
       ENDIF
    ENDIF

    !========================================================================
    ! Archive quantities for satellite diagnostics (if requested)
    !========================================================================
    IF ( State_Diag%Archive_SatDiagn ) THEN
       CALL Do_Archive_SatDiagn( Input_Opt,  State_Chm,  State_Diag,         &
                                 State_Grid, State_Met,  RC                 )
    ENDIF


    ! Error handling
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error converting species units for archiving diagnostics #2'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF


  END SUBROUTINE Set_Diagnostics_EndofTimestep
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Zero_Diagnostics_StartofTimestep
!
! !DESCRIPTION: This routine sets certain diagnostic arrays to zero. This
!  is intended for diagnostics that must be reset to zero each timestep but
!  that do not have a clear place in the source code execution for doing this,
!  generally because they are set in multiple places.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Zero_Diagnostics_StartofTimestep( Input_Opt, State_Diag, RC )
!
! !USES:
!
    USE Input_Opt_Mod,    ONLY : OptInput
    USE State_Diag_Mod,   ONLY : DgnState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt    ! Input Options object
!
! !INPUT AND OUTPUT PARAMETERS:
!
    TYPE(DgnState),   INTENT(INOUT) :: State_Diag     ! Diagnostics state obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC
!
! !REVISION HISTORY:
!  01 Feb 2018 - E. Lundgren - initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255)      :: ErrMsg, thisLoc

    !=======================================================================
    ! Zero_Diagnostics_StartofTimestep begins here
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = &
    ' -> at Zero_Diagnostics_StartofTimestep (in GeosCore/diagnostics_mod.F90)'

    ! Mercury simulation
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

       IF ( State_Diag%Archive_DryDepChm .or. State_Diag%Archive_DryDep ) THEN
          State_Diag%DryDepChm = 0.0_f4
       ENDIF

       IF ( State_Diag%Archive_EmisHg2rivers ) THEN
          State_Diag%EmisHg2rivers = 0.0_f4
       ENDIF

       IF ( State_Diag%Archive_EmisHg2snowToOcean ) THEN
          State_Diag%EmisHg2snowToOcean = 0.0_f4
       ENDIF

      IF ( State_Diag%Archive_FluxHg0fromAirToOcean ) THEN
          State_Diag%FluxHg0fromAirToOcean = 0.0_f4
       ENDIF

       IF ( State_Diag%Archive_FluxHg0fromOceantoAir ) THEN
          State_Diag%FluxHg0fromOceanToAir = 0.0_f4
       ENDIF

       IF ( State_Diag%Archive_FluxHg2HgPfromAirToOcean ) THEN
          State_Diag%FluxHg2HgPfromAirToOcean = 0.0_f4
       ENDIF

       IF ( State_Diag%Archive_FluxOCtoDeepOcean ) THEN
          State_Diag%FluxOCtoDeepOcean = 0.0_f4
       ENDIF

    ENDIF

    ! Dry deposition
    IF ( Input_Opt%LDRYD ) THEN
       ! Initialize the DryDepMix diagnostic array for the History Component.
       ! This will prevent leftover values from being carried over to this
       ! timestep. (For example, if on the last iteration, the PBL height
       ! was higher than it is now, then we will have stored drydep fluxes
       ! up to that height, so we need to zero these out.)
       IF ( State_Diag%Archive_DryDepMix .or. State_Diag%Archive_DryDep ) THEN
          State_Diag%DryDepMix = 0.0_f4
       ENDIF
    ENDIF

  END SUBROUTINE Zero_Diagnostics_StartofTimestep
!EOC
#ifdef ADJOINT
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_SpcAdj_Diagnostic
!
! !DESCRIPTION: Subroutine Set_SpcAdj\_Diagnostic sets the passed species
!  adjoint diagnostic array stored in State_Diag to the instantaneous
!  State_Chm%SpeciesAdj values converted to the diagnostic unit stored in
!  the State_Diag metadata.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_SpcAdj_Diagnostic( Input_Opt, State_Chm, State_Diag,    &
                                    State_Grid, State_Met, RC           )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnMap
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt      ! Input Options object
    TYPE(GrdState),   INTENT(IN)  :: State_Grid     ! Grid state object
    TYPE(MetState),   INTENT(IN)  :: State_Met      ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm    ! Chemistry State object
    TYPE(DgnState),   INTENT(INOUT) :: State_Diag   ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC      ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY:
!  15 Dec 2019 - C. Lee - Initial version
!  17 Dec 2020 - C. Lee - Updated to account for changes to Set_SpcConcs
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL               :: Found
    INTEGER               :: N, S
    REAL(fp)              :: LT,  GOOD

    ! Strings
    CHARACTER(LEN=255)    :: ErrMsg, ThisLoc

    ! Objects
    TYPE(DgnMap), POINTER :: mapData

    !====================================================================
    ! Set_SpcAdj_Diagnostic begins here!
    !====================================================================

    ! Assume success
    RC      =  GC_SUCCESS
    Found   = .FALSE.
    ThisLoc = ' -> Set_SpcAdj_Diagnostic (in GeosCore/diagnostics_mod.F90)'

    ! Verify that incoming State_Chm%Species units are kg/kg dry air.
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
       ErrMsg = 'Incorrect species units in Set_SpcAdj_Diags_VVDry!' // trim( State_Chm%Spc_Units)
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Copy species to SpeciesAdj (concentrations diagnostic) [v/v dry]
    !=======================================================================
    IF ( Input_Opt%Is_Adjoint ) THEN

       ! Point to mapping obj specific to SpeciesAdj diagnostic collection
       mapData => State_Diag%Map_SpeciesAdj

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( N, S   )
       DO S = 1, mapData%nSlots
          N = mapData%slot2id(S)
          State_Diag%SpeciesAdj(:,:,:,S) = State_Chm%SpeciesAdj(:,:,:,N)
       ENDDO
       !$OMP END PARALLEL DO

       ! Free pointer
       mapData => NULL()

    ENDIF

    !=======================================================================
    ! Copy species to SatDiagn (satellite diagnostic output) [v/v dry]
    !=======================================================================
    IF ( State_Diag%Archive_SatDiagnConc ) THEN

       ! Point to mapping obj specific to species boundary conditions
       mapData => State_Diag%Map_SatDiagnConc

       ! Loop over the number of advected species that we wish
       ! to save at a user-specified local time range
       ! Loop over longitudes:
       DO I = 1, State_Grid%NX

          ! Get local time in hours:
          LT = GET_LOCALTIME(I, 1, 1, State_Grid)
          IF ( LT < 0  ) LT = LT + 24e+0_fp
          ! Check if local time is during satellite overpass time:
          IF ( LT >= State_Diag%SatDiagn_StartHr .and. &
               LT <= State_Diag%SatDiagn_EndHr ) THEN

             ! GOOD = 1 if during local time range, 0 otherwise:
             GOOD = 1e+0_fp

          ELSE

             GOOD = 0e+0_fp

          ENDIF

          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( N, S   )
          DO S = 1, mapData%nSlots
             N = mapData%slot2id(S)
             State_Diag%SatDiagnConc(I,:,:,S) = TmpSpcArr(I,:,:,N) * GOOD
          ENDDO
          !$OMP END PARALLEL DO

       ENDDO

       ! Free pointer
       mapData => NULL()

    ENDIF

  END SUBROUTINE Set_SpcAdj_Diagnostic
!EOC
#endif
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_SpcConc_Diags_VVDry
!
! !DESCRIPTION: Subroutine Set_SpcConc\_DiagVVDry sets several species
!  concentration diagnostic arrays stored in State_Diag to the instantaneous
!  State_Chm%Species values (in units of "v/v, dry air").
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_SpcConc_Diags_VVDry( Input_Opt,  State_Chm, State_Diag,     &
                                      State_Grid, State_Met, RC            )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE PhysConstants,  ONLY : AIRMW
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnMap
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE Time_Mod,       ONLY : Get_LocalTime
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt    ! Input Options object
    TYPE(GrdState),   INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState),   INTENT(IN)    :: State_Met    ! Meteorology State obj
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm    ! Chemistry State object
    TYPE(DgnState),   INTENT(INOUT) :: State_Diag   ! Diagnsotics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  08 Jul 2019 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL               :: Found
    INTEGER               :: D, I, J, L, N, S
    REAL(fp)              :: TmpVal, Conv, LT, GOOD

    ! Strings
    CHARACTER(LEN=255)    :: ErrMsg, ThisLoc

    ! Objects
    TYPE(DgnMap), POINTER :: mapData

    ! Arrays
    REAL(fp)              :: TmpSpcArr(State_Grid%NX,State_Grid%NY, &
                                       State_Grid%NZ,State_Chm%nSpecies)

    !====================================================================
    ! Set_SpcConc_Diags_VVDry begins here!
    !====================================================================

    ! Assume success
    RC      =  GC_SUCCESS
    Found   = .FALSE.
    ThisLoc = &
         ' -> at Set_SpcConc_Diags_VVDry (in GeosCore/diagnostics_mod.F90)'

    ! Verify that incoming State_Chm%Species units are kg/kg dry air.
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
       ErrMsg = 'Incorrect species units in Set_SpcConc_Diags_VVDry!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Store species in v/v dry as temporary variable if diagnostics on
    IF ( State_Diag%Archive_SpeciesConcVV                               .or. &
         State_Diag%Archive_SpeciesBC                                   .or. &
         State_Diag%Archive_SpeciesRst                                  .or. &
         State_Diag%Archive_ConcAboveSfc                                .or. & 
         State_Diag%Archive_SatDiagnConc                              ) THEN

       !$OMP PARALLEL DO                                                     &
       !$OMP DEFAULT( SHARED                                                )&
       !$OMP PRIVATE( I, J, L, N                                            )&
       !$OMP COLLAPSE( 4                                                    )
       DO N = 1, State_Chm%nSpecies
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          TmpSpcArr(I,J,L,N) = State_Chm%Species(N)%Conc(I,J,L) *            &
                               ( AIRMW / State_Chm%SpcData(N)%Info%MW_g )
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

    !=======================================================================
    ! Copy species to SpeciesConc (concentrations diagnostic) [v/v dry]
    !=======================================================================
    IF ( State_Diag%Archive_SpeciesConcVV ) THEN

       ! Point to mapping obj specific to SpeciesConcVV diagnostic collection
       mapData => State_Diag%Map_SpeciesConcVV

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( N, S   )
       DO S = 1, mapData%nSlots
          N = mapData%slot2id(S)
          State_Diag%SpeciesConcVV(:,:,:,S) = TmpSpcArr(:,:,:,N)
       ENDDO
       !$OMP END PARALLEL DO

       ! Free pointer
       mapData => NULL()

    ENDIF

    !=======================================================================
    ! Copy species to SatDiagn (satellite diagnostic output) [v/v dry]
    !=======================================================================
    IF ( State_Diag%Archive_SatDiagnConc ) THEN

       ! Loop over longitudes
       !$OMP PARALLEL DO                                                    &
       !$OMP DEFAULT( SHARED                                               )&
       !$OMP PRIVATE( I, LT, GOOD, S, N                                    )
       DO I = 1, State_Grid%NX

          ! Get local time in hours
          LT = Get_LocalTime( I, 1, 1, State_Grid )
          IF ( LT < 0  ) LT = LT + 24.0_fp

          ! Check if local time is during satellite overpass time:
          ! GOOD = 1 if during local time range, 0 otherwise
          GOOD = 0e+0_fp
          IF ( LT >= State_Diag%SatDiagn_StartHr .and.                      &
               LT <= State_Diag%SatDiagn_EndHr ) GOOD = 1e+0_fp

          ! Archie into SatDiagnConc diagnostic array
          DO S = 1, State_Diag%Map_SatDiagnConc%nSlots
             N = State_Diag%Map_SatDiagnConc%slot2id(S)
             State_Diag%SatDiagnConc(I,:,:,S) = TmpSpcArr(I,:,:,N) * GOOD
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !=======================================================================
    ! Copy species to SpeciesBC (transport boundary conditions) [v/v dry]
    !=======================================================================
    IF ( State_Diag%Archive_SpeciesBC ) THEN

       ! Point to mapping obj specific to species boundary conditions
       mapData => State_Diag%Map_SpeciesBC

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( N, S   )
       DO S = 1, mapData%nSlots
          N = mapData%slot2id(S)
          State_Diag%SpeciesBC(:,:,:,S) = TmpSpcArr(:,:,:,N)
       ENDDO
       !$OMP END PARALLEL DO

       ! Free pointer
       mapData => NULL()

    ENDIF

    !=======================================================================
    ! Copy species to SpeciesRst (restart file output) [v/v dry]
    !=======================================================================
    IF ( State_Diag%Archive_SpeciesRst ) THEN
       State_Diag%SpeciesRst(:,:,:,:) = TmpSpcArr(:,:,:,:)
    ENDIF

    !=======================================================================
    ! Diagnostic for correcting species concentrations from the height
    ! of the lowest model level to the surface.
    !
    ! Use this diagnostic to correct species concentration values from
    ! (typically for O3 or HNO3) from the lowest model layer, ~60m,
    ! to the surface.
    !
    !    C(Zc) = [ 1 - Ra(Z1,Zc) * Vd(Z1) ] * C(Z1)
    !
    ! where
    !    Ra(Z1,ZC) is the aerodynamic resistance between Z1 and ZC,
    !
    !    Vd(Z1) is the ozone deposition velocity at Z1, and
    !
    !    C(Z1) is the ozone concentration at Z1.
    !
    ! Ra(Z1,Zc) is calculated to the lowest model level in drydep_mod.F90.
    ! We recalculate Ra using Z1 using a value specified in geoschem_config.yml;
    ! usually 10m, which is the height of the CASTNET measurement for O3.
    ! This new Ra is stored in State_Diag%DryDepRaALT1.
    !
    ! References:
    ! (1) Travis, K.R., et al, "Resolving vertical ozone vertical gradients
    !      in air quality models, Atmos. Chem. Phys. Disc., 2017.
    ! (2) Zhang, L.,et al, "Nitrogen deposition to the United States:
    !      distribution, sources, and processes" Atmos. Chem. Phys.,
    !      12, 4,539-4,4554, 2012.
    !=======================================================================
    IF ( State_Diag%Archive_ConcAboveSfc ) THEN

       ! Loop over the number of drydep species that we wish
       ! to save at a user-specified altitude above the surface
       !$OMP PARALLEL DO                         &
       !$OMP DEFAULT( SHARED                   ) &
       !$OMP PRIVATE( D, N, I, J, TmpVal, Conv )
       DO D = 1, State_Chm%nDryAlt

          ! Get the corresponding species index and drydep index
          N = State_Chm%Map_DryAlt(D)

          ! Loop over surface locations
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX

             ! Species concentration [v/v dry]
             TmpVal = TmpSpcArr(I,J,1,N)

             ! Conversion factor used to translate from
             ! lowest model layer (~60m) to the surface
             Conv = ( 1.0_fp                                              &
                  -   ( State_Diag%DryDepRaALT1(I,J) / 100.0_fp )         &
                  *   State_Diag%DryDepVelForALT1(I,J,D)                 )

             ! Do not let CONV go negative
             IF ( Conv < 0.0_fp ) Conv = 1.0_fp

             ! Save concentration at the user-defined altitude
             ! as defined in geoschem_config.yml (usually 10m).
             State_Diag%SpeciesConcALT1(I,J,D) = TmpVal * Conv

          ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

  END SUBROUTINE Set_SpcConc_Diags_VVDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_SpcConc_Diags_MND
!
! !DESCRIPTION: Subroutine Set_SpcConc\_Diags\_MND sets several species
!  concentration diagnostic arrays stored in State_Diag to the instantaneous
!  State_Chm%Species values (in units of "v/v, dry air").
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_SpcConc_Diags_MND  ( Input_Opt,  State_Chm, State_Diag,     &
                                      State_Grid, State_Met, RC            )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE PhysConstants,  ONLY : AVO
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnMap
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt    ! Input Options object
    TYPE(GrdState),   INTENT(IN)    :: State_Grid   ! Grid State object
    TYPE(MetState),   INTENT(IN)    :: State_Met    ! Meteorology State obj
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm    ! Chemistry State object
    TYPE(DgnState),   INTENT(INOUT) :: State_Diag   ! Diagnsotics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC           ! Success or failure?
!
! !REVISION HISTORY:
!  21 Nov 2021 - H.P. Lin    - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER               :: N, S
    REAL(fp)              :: MW_kg

    ! Strings
    CHARACTER(LEN=255)    :: ErrMsg, ThisLoc, OrigUnit

    ! Objects
    TYPE(DgnMap), POINTER :: mapData

    !====================================================================
    ! Set_SpcConc_Diags_MND begins here!
    !====================================================================

    ! Assume success
    RC      =  GC_SUCCESS
    ThisLoc = &
         ' -> at Set_SpcConc_Diags_MND (in GeosCore/diagnostics_mod.F90)'

    ! Verify that incoming State_Chm%Species units are kg/kg dry air.
    IF ( TRIM( State_Chm%Spc_Units ) /= 'kg/kg dry' ) THEN
       ErrMsg = 'Incorrect species units in Set_SpcConc_Diags_MND!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Copy species to SpeciesConcMND (concentrations diagnostic) [molec/cm3]
    !=======================================================================
    IF ( State_Diag%Archive_SpeciesConcMND ) THEN

       ! Point to mapping obj specific to SpeciesConcMND diagnostic collection
       mapData => State_Diag%Map_SpeciesConcMND

       !$OMP PARALLEL DO            &
       !$OMP DEFAULT( SHARED      ) &
       !$OMP PRIVATE( N, S, MW_kg )
       DO S = 1, mapData%nSlots
          N = mapData%slot2id(S)

          ! Molecular weight for the species [kg]
          MW_kg = State_Chm%SpcData(N)%Info%MW_g * 1.e-3_fp

          State_Diag%SpeciesConcMND(:,:,:,S) =                      &
               State_Chm%Species(N)%Conc(:,:,:) *                   &
               State_Met%AIRDEN(:,:,:) * ( AVO / MW_kg ) / 1e+6_fp
       ENDDO
       !$OMP END PARALLEL DO

       ! Free pointer
       mapData => NULL()

    ENDIF

  END SUBROUTINE Set_SpcConc_Diags_MND
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute_Column_Mass
!
! !DESCRIPTION: Subroutine Compute\_Budget\_Diagnostics calculates the
!  budget diagnostics for a given component by taking the difference of the
!  final and initial kg per grid cell and dividing by the timestep in seconds
!  to get kg/s.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Budget_Diagnostics( Input_Opt,   State_Chm, State_Grid, &
                                         State_Met,   isFull,    diagFull,   &
                                         mapDataFull, isTrop,    diagTrop,   &
                                         mapDataTrop, isPBL,     diagPBL,    &
                                         mapDataPBL,  colMass,   RC,         &
                                         timeStep,    isWetDep,  before_op  )
!
! !USES:
!
    USE Input_Opt_Mod,  Only : OptInput
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnMap
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt        ! Input options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid       ! Grid state object
    TYPE(MetState), INTENT(IN)    :: State_Met        ! Meteorology state obj
    LOGICAL,        INTENT(IN)    :: isFull           ! T if full col diag on
    TYPE(DgnMap),   POINTER       :: mapDataFull      ! Map to species indexes
    LOGICAL,        INTENT(IN)    :: isTrop           ! T if trop col diag on
    TYPE(DgnMap),   POINTER       :: mapDataTrop      ! Map to species indexes
    LOGICAL,        INTENT(IN)    :: isPBL            ! T if PBL col diag on
    TYPE(DgnMap),   POINTER       :: mapDataPBL       ! Map to species indexes
    LOGICAL,        OPTIONAL      :: isWetDep         ! T = wetdep budgets
    LOGICAL,        OPTIONAL      :: before_op        ! T = before operation
    REAL(f8),       OPTIONAL      :: timestep         ! F = after operation
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm        ! Chemistry state obj
    REAL(f8),       POINTER       :: diagFull(:,:,:)  ! ptr to full col diag
    REAL(f8),       POINTER       :: diagTrop(:,:,:)  ! ptr to trop col diag
    REAL(f8),       POINTER       :: diagPBL(:,:,:)   ! ptr to pbl col diag
    REAL(f8),       POINTER       :: colMass(:,:,:,:) ! Initial column mass
                                                      ! (I,J,spc,col region)
                                                      ! 1:full, 2:trop, 3:pbl
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC               ! Success or failure?
!
! !REVISION HISTORY:
!  28 Aug 2018 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    LOGICAL            :: after,  before, wetDep
    INTEGER            :: I,      J,      L,       N
    INTEGER            :: numSpc, region, topLev,  S
    REAL(f8)           :: colSum, dt

    ! Arrays
    REAL(f8)           :: spcMass(State_Grid%NZ)

    ! Strings
    CHARACTER(LEN=255) :: errMsg, thisLoc

    !====================================================================
    ! Compute_Budget_Diagnostics begins here!
    !====================================================================

    ! Initialize
    RC      =  GC_SUCCESS
    errMsg  = ''
    ThisLoc = ' -> at Compute_Column_Mass (in GeosCore/diagnostics_mod.F90)'
    colSum  = 0.0_f8
    spcMass = 0.0_f8

    ! Exit if concentrations are not in [kg/kg dry]
    IF ( State_Chm%Spc_Units /= 'kg/kg dry' ) THEN
       errMsg = 'State_Chm%Species units must be kg/kg dry. ' // &
                'Incorrect units: '// TRIM( State_Chm%Spc_Units )
       CALL GC_Error( errMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Set logicals to denote if we are calling this routine
    ! before the operation or after the operation
    IF ( PRESENT( before_op ) ) THEN
       before = before_op
    ELSE
       before = .FALSE.
    ENDIF
    after = ( .not. before )

    ! Test if the budgets are for wetdep species
    IF ( PRESENT( isWetDep ) ) THEN
       wetDep = isWetDep
    ELSE
       wetDep = .FALSE.
    ENDIF

    ! Make sure the timeStep argument is passed (if after operation)
    IF ( after .and. ( .not. PRESENT( timeStep ) ) ) THEN
       errMsg = 'The timeStep argument was not passed!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    ! Make sure mapDataFull and diagFull are not undefined
    IF ( isFull ) THEN
       IF ( .not. ASSOCIATED( mapDataFull ) ) THEN
          errMsg = 'The mapDataFull object is undefined!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       IF ( after .and. ( .not. ASSOCIATED( diagFull ) ) ) THEN
          errMsg = 'The diagFull array is undefined!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Make sure mapDataTrop and diagTrop are not undefined
    IF ( isTrop ) THEN
       IF ( .not. ASSOCIATED( mapDataTrop ) ) THEN
          errMsg = 'The mapDataTrop object is undefined!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       IF ( after .and. ( .not. ASSOCIATED( diagTrop ) ) ) THEN
          errMsg = 'The diagTrop array is undefined!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Make sure mapDataPBL and diagPBL are not undefined
    IF ( isPBL ) THEN
       IF ( .not. ASSOCIATED( mapDataPBL ) ) THEN
          errMsg = 'The mapDataPBL object is undefined!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       IF ( after .and. ( .not. ASSOCIATED( diagPBL ) ) ) THEN
          errMsg = 'The diagPBL array is undefined!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
    ENDIF

    ! Make sure the colMass array is not undefined
    IF ( .not. ASSOCIATED( colMass ) ) THEN
       errMsg = 'The colMass array is undefined!'
       CALL GC_Error( errMsg, RC, thisLoc )
       RETURN
    ENDIF

    !====================================================================
    ! Before operation: Compute column masses (full, trop, PBL)
    !
    ! After operation:  Compute column differences (final-initial)
    !                   and them update diagnostic arrays
    !====================================================================

    ! Zero out the column mass array if we are calling this routine
    ! before the desired operation.  This will let us compute initial mass.
    IF ( before ) THEN
       colMass = 0.0_f8
    ENDIF

    ! Loop over NX and NY dimensions
    !$OMP PARALLEL DO                                        &
    !$OMP DEFAULT( SHARED                                  ) &
    !$OMP PRIVATE( I, J, colSum, spcMass, topLev, S, N, L  )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Zero column-specific variables
       colSum  = 0.0_f8
       spcMass = 0.0_f8
       topLev  = 0

       !--------------------------------------------------------------------
       ! Full-column budget for requested species
       !--------------------------------------------------------------------
       IF ( isFull ) THEN

          ! Loop over # of diagnostic slots
          DO S = 1, mapDataFull%nSlots

             ! Initialize column-specfic variables
             colSum  = 0.0_f8
             spcMass = 0.0_f8

             ! For wetdep budgets, translate wetdep ID to modelId
             ! Otherwise, get the modelId from the slotId
             IF ( wetDep ) THEN
                N = State_Chm%Map_WetDep(mapDataFull%slot2Id(S))
             ELSE
                N = mapDataFull%slot2Id(S)
             ENDIF

             ! Compute mass at each grid box in the column [kg]
             DO L = 1, State_Grid%NZ
                spcMass(L) = State_Chm%Species(N)%Conc(I,J,L) * &
                             State_Met%AD(I,J,L)
             ENDDO

             ! Compute the full-atmosphere column mass [kg]
             colSum = SUM( spcMass(1:State_Grid%NZ)  )

             ! Before operation: Compute initial full-atm column mass
             ! After operation: Compute change in column mass (final-initial),
             ! convert to [kg/s], and store in the diagFull array.
             IF ( before ) THEN
                colMass(I,J,N,1) = colSum
             ELSE
#ifdef MODEL_GEOS
                diagFull(I,J,S) = ( colSum - colMass(I,J,N,1) ) / timeStep &
                                / State_Grid%AREA_M2(I,J)
#else
                diagFull(I,J,S) = ( colSum - colMass(I,J,N,1) ) / timeStep
#endif
             ENDIF
          ENDDO
       ENDIF

       !---------------------------------------------------------------------
       ! Troposphere-only budget for each requested species
       !---------------------------------------------------------------------
       IF ( isTrop ) THEN

          ! Top level in the column
          topLev = State_Met%TropLev(I,J)

          ! Loop over # of diagnostic slots
          DO S = 1, mapDataTrop%nSlots

             ! Initialize column-specfic variables
             colSum  = 0.0_f8
             spcMass = 0.0_f8

             ! For wetdep budgets, translate wetdep ID to modelId
             ! Otherwise, get the modelId from the slotId
             IF ( wetDep ) THEN
                N = State_Chm%Map_WetDep(mapDataTrop%slot2Id(S))
             ELSE
                N = mapDataTrop%slot2Id(S)
             ENDIF

             ! Compute mass at each grid box in the troposphere [kg]
             DO L = 1, topLev
                spcMass(L) = State_Chm%Species(N)%Conc(I,J,L) * &
                             State_Met%AD(I,J,L)
             ENDDO

             ! Compute the trop-column mass [kg]
             colSum = SUM( spcMass(1:topLev) )

             ! Before operation: Compute initial trop-column mass
             ! After operation: Compute change in column mass (final-initial),
             ! convert to [kg/s], and store in the diagTrop array.
             IF ( before ) THEN
                colMass(I,J,N,2) = colSum
             ELSE
#ifdef MODEL_GEOS
                diagTrop(I,J,S) = ( colSum - colMass(I,J,N,2) ) / timeStep &
                                / State_Grid%AREA_M2(I,J)
#else
                diagTrop(I,J,S) = ( colSum - colMass(I,J,N,2) ) / timeStep
#endif
             ENDIF
          ENDDO
       ENDIF

       !---------------------------------------------------------------------
       ! PBL-only budget for each requested species
       !---------------------------------------------------------------------
       IF ( isPBL ) THEN

          ! Top level of column is where PBL top occurs
          topLev = MAX( 1, FLOOR( State_Met%PBL_TOP_L(I,J) ) )

          ! Loop over # of diagnostic slots
          DO S = 1, mapDataPBL%nSlots

             ! Initialize column-specfic variables
             colSum  = 0.0_f8
             spcMass = 0.0_f8

             ! For wetdep budgets, translate wetdep ID to modelId
             ! Otherwise, get the modelId from the slotId
             IF ( wetDep ) THEN
                N = State_Chm%Map_WetDep(mapDataPBL%slot2Id(S))
             ELSE
                N = mapDataPBL%slot2Id(S)
             ENDIF

             ! Compute mass at each grid box in the column [kg]
             DO L = 1, topLev
                spcMass(L) = State_Chm%Species(N)%Conc(I,J,L) * &
                             State_Met%AD(I,J,L)
             ENDDO

             ! Compute column mass in PBL region [kg]
             colSum = SUM( spcMass(1:topLev) )

             ! Before operation: Compute initial PBL-column mass
             ! After operation: Compute change in column mass (final-initial),
             ! convert to [kg/s], and store in the diagPBL array.
             IF ( before ) THEN
                colMass(I,J,N,3) = colSum
             ELSE
#ifdef MODEL_GEOS
                diagPBL(I,J,S) = ( colSum - colMass(I,J,N,3) ) / timeStep &
                               / State_Grid%AREA_M2(I,J)
#else
                diagPBL(I,J,S) = ( colSum - colMass(I,J,N,3) ) / timeStep
#endif
             ENDIF
          ENDDO
       ENDIF

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

  END SUBROUTINE Compute_Budget_Diagnostics
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Archive_SatDiagn
!
! !DESCRIPTION: Masks satellite diagnostic fields by the requested local
!  time window.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Archive_SatDiagn( Input_Opt,  State_Chm,  State_Diag,        &
                                  State_Grid, State_Met,  RC                )
!
! !USES:
!
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE PhysConstants,  ONLY : AVO
    USE Species_Mod,    ONLY : SpcConc
    USE State_Chm_Mod,  ONLY : ChmState, Ind_
    USE State_Diag_Mod, ONLY : DgnState, DgnMap
    USE State_Grid_Mod, ONLY : GrdState
    USE State_Met_Mod,  ONLY : MetState
    USE Time_Mod,       ONLY : Get_LocalTime
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
    TYPE(DgnState), INTENT(INOUT) :: State_Diag   ! Diagnostic State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC           ! Success or failure?
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! SAVEd scalars
    LOGICAL, SAVE          :: first = .TRUE.
    INTEGER, SAVE          :: id_OH = -1

    ! Scalars
    INTEGER                :: I,    N,      S
    REAL(fp)               :: good, locTime

    ! Strings
    CHARACTER(LEN=255)     :: thisLoc
    CHARACTER(LEN=512)     :: errMsg

    ! Pointers & Objects
    TYPE(SpcConc), POINTER :: Spc(:)

    !=======================================================================
    ! Do_Archive_SatDiagn begins here!
    !=======================================================================

    ! Initialize
    RC      =  GC_SUCCESS
    good    =  0.0_fp
    locTime =  0.0_fp
    Spc     => State_Chm%Species
    errMsg  = ''
    thisLoc = &
     ' -> at Do_Archive_SatDiagn (in module GeosCore/diagnostics_mod.F90)'

    ! Get the species ID for OH if this is the first call
    IF ( first ) THEN
       id_OH  = Ind_('OH')
       IF ( id_OH < 0 ) THEN
          errMsg = 'OH is not a defined species in this simulation!!!'
          CALL GC_Error( errMsg, RC, thisLoc )
          RETURN
       ENDIF
       first= .FALSE.
    ENDIF

    !========================================================================
    ! Archive satellite diagnostics
    !========================================================================

    ! Loop over longitudes
    !$OMP PARALLEL DO                                                        &
    !$OMP DEFAULT( SHARED                                                   )&
    !$OMP PRIVATE( I, locTime, good, S, N                                   )
    DO I = 1, State_Grid%NX

       !---------------------------------------------------------------------
       ! Local time
       !---------------------------------------------------------------------

       ! Get local time (and make sure it isn't negative)
       locTime = Get_LocalTime( I, 1, 1, State_Grid )
       IF ( locTime < 0 ) locTime = locTime + 24.0_fp

       ! Determine whether during satellite overpass time window:
       ! good = 1 if during local time range, 0 otherwise:
       !
       !%%% TODO This should be a property of the HISTORY container
       !%%% rather than SatDiagn.  This will prevent multiple SatDiagn
       !%%% collections from being run at once.  Fix this later.
       !%%%   -- Bob Yantosca (01 Nov 2022)
       good = 0.0_fp
       IF ( locTime >= State_Diag%SatDiagn_StartHr  .and.                    &
            locTime <= State_Diag%SatDiagn_EndHr   ) good = 1.0_fp

       !---------------------------------------------------------------------
       ! SatDiagnCount: Count number of local times
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnCount ) THEN
          State_Diag%SatDiagnCount(I,:,:) = &
          State_Diag%SatDiagnCount(I,:,:) + good
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagnOH: OH concentration [molec/cm3]:
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnOH ) THEN
          State_Diag%SatDiagnOH(I,:,1:State_Grid%MaxChemLev) =               &
               ( Spc(id_OH)%Conc(I,:,1:State_Grid%MaxChemLev) * good  *      &
               State_Met%AIRDEN(I,:,1:State_Grid%MaxChemLev)  * AVO ) /      &
               ( State_Chm%SpcData(id_OH)%Info%MW_g ) / 1.0e+3_fp
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagnRH: Relative humidity [%]
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnRH ) THEN
          State_Diag%SatDiagnRH(I,:,:) = State_Met%RH(I,:,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagnAirDen: Air density [molec/cm3]
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnAirDen ) THEN
          State_Diag%SatDiagnAirDen(I,:,:) =                                 &
               State_Met%AirNumDen(I,:,:)  * good
       ENDIF

       !---------------------------------------------------------------------
       ! Grid box height [m]:
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnBoxHeight ) THEN
          State_Diag%SatDiagnBoxHeight(I,:,:) =                              &
               State_Met%BXHEIGHT(I,:,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! Pressure edges [hPa]:
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnPEdge ) THEN
          State_Diag%SatDiagnPEdge(I,:,:) = &
               State_Met%PEDGE(I,:,1:State_Grid%NZ) * good
       ENDIF

       !---------------------------------------------------------------------
       ! Tropopause pressure [hPa]:
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnTROPP ) THEN
          State_Diag%SatDiagnTROPP(I,:) = State_Met%TROPP(I,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! PBL Height [m]:
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnPBLHeight ) THEN
          State_Diag%SatDiagnPBLHeight(I,:) = State_Met%PBLH(I,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! PBL Top [m]:
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnPBLTop ) THEN
          State_Diag%SatDiagnPBLTop(I,:) = State_Met%PBL_TOP_m(I,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! Air temperature [K]: Temperature Interpolated to Current Time
       ! This temperture is interpolated from 3 h Met Field (TMPU1 and TMPU2)
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnTAir ) THEN
          State_Diag%SatDiagnTAir(I,:,:) = State_Met%T(I,:,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! Root Zone Soil Moisture (or Wetness) [fraction]:
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnGWETROOT ) THEN
          State_Diag%SatDiagnGWETROOT(I,:) = State_Met%GWETROOT(I,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! Topsoil Moisture (or Wetness) [fraction]:
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnGWETTOP ) THEN
          State_Diag%SatDiagnGWETTOP(I,:) = State_Met%GWETTOP(I,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! Direct Photosynthetically Active Radiation [W/m2]:
       ! Aka Surface downward PAR beam flux
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnPARDR ) THEN
          State_Diag%SatDiagnPARDR(I,:) = State_Met%PARDR(I,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! Diffuse Photosynthetically Active Radiation [W/m2]:
       ! Aka Surface downward PAR diffuse flux
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnPARDF ) THEN
          State_Diag%SatDiagnPARDF(I,:) = State_Met%PARDF(I,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! Total Precipitation (at surface) [mm/day]:
       ! Documentation says this variable is converted from original
       ! units of kg/m2/s
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnPRECTOT ) THEN
          State_Diag%SatDiagnPRECTOT(I,:) = State_Met%PRECTOT(I,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! Sea Level Pressure [hPa]:
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnSLP ) THEN
          State_Diag%SatDiagnSLP(I,:) = State_Met%SLP(I,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! Specific Humidity Interpolated to Current Time [g H2O/kg air]:
       ! Linearly interpolated from 3 h met field (SPHU1 and SPHU2)
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnSPHU ) THEN
          State_Diag%SatDiagnSPHU(I,:,:) = State_Met%SPHU(I,:,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! Surface Temperature at 2m [K]
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnTS ) THEN
          State_Diag%SatDiagnTS(I,:) = State_Met%TS(I,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! PBL Top Height [Levels]:
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnPBLTOPL ) THEN
          State_Diag%SatDiagnPBLTOPL(I,:) = State_Met%PBL_TOP_L(I,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! MODIS Daily LAI [m2/m2]:
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnMODISLAI ) THEN
          State_Diag%SatDiagnMODISLAI(I,:) = State_Met%MODISLAI(I,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagn Diagnostic for WetLossLS [units of kg/s as per WetLossLS]
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnWetLossLS ) THEN
          DO S = 1, State_Diag%Map_SatDiagnWetLossLS%nSlots
             State_Diag%SatDiagnWetLossLS(I,:,:,S) =                         &
             State_Diag%SatDiagnWetLossLS(I,:,:,S) * good
          ENDDO
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagn Diagnostic for WetLossConv [units of kg/s as per WetLossConv]
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnWetLossConv ) THEN
          DO S = 1, State_Diag%Map_SatDiagnWetLossConv%nSlots
             State_Diag%SatDiagnWetLossConv(I,:,:,S) =                       &
             State_Diag%SatDiagnWetLossConv(I,:,:,S) * good
          ENDDO
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagn Diagnostic for Jval [units of s-1 as per Jval]
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnJval ) THEN
          DO S = 1, State_Diag%Map_SatDiagnJval%nSlots
             State_Diag%SatDiagnJval(I,:,:,S) =                              &
             State_Diag%SatDiagnJval(I,:,:,S) * good
          ENDDO
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagn Diagnostic for JvalO3O1D [units of s-1 as per JvalO3O1D]:
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnJvalO3O1D ) THEN
          State_Diag%SatDiagnJvalO3O1D(I,:,:) =                              &
          State_Diag%SatDiagnJvalO3O1D(I,:,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagn Diagnostic for JvalO3O3P [units of s-1 as per JvalO3O3P]:
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnJvalO3O3P ) THEN
          State_Diag%SatDiagnJvalO3O3P(I,:,:) =                              &
          State_Diag%SatDiagnJvalO3O3P(I,:,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagn Diagnostic for DryDep [units of molec cm-2 s-1 as per DryDep]
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnDryDep ) THEN
          DO S = 1, State_Diag%Map_SatDiagnDryDep%nSlots
             State_Diag%SatDiagnDryDep(I,:,S) =                              &
             State_Diag%SatDiagnDryDep(I,:,S) * good
          ENDDO
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagn Diagnostic for DryDepVel [units of cm s-1 as per DryDepVel]
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnDryDepVel ) THEN
          DO S = 1, State_Diag%Map_SatDiagnDryDepVel%nSlots
             State_Diag%SatDiagnDryDepVel(I,:,S) =                           &
             State_Diag%SatDiagnDryDepVel(I,:,S) * good
          ENDDO
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagn Diagnostic for OH Reactivity [units of s-1]
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnOHreactivity ) THEN
          State_Diag%SatDiagnOHreactivity(I,:,:) =                           &
          State_Diag%SatDiagnOHreactivity(I,:,:) * good
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagn Diagnostic for Column Emissions (ColEmis) [units of kg/m2/s]:
       ! From surface to maximum vertical level for advected species
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnColEmis ) THEN
          DO S = 1, State_Diag%Map_SatDiagnColEmis%nSlots
             State_Diag%SatDiagnColEmis(I,:,S) =                             &
             State_Diag%SatDiagnColEmis(I,:,S) * good
          ENDDO
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagn Diagnostic for Total Surface Fluxes [units of kg/m2/s]:
       ! From surface to top of the PBL for Advected Species
       ! (eflx (emis) - dflx(drydep)))
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnSurfFlux ) THEN
          DO S = 1, State_Diag%Map_SatDiagnSurfFlux%nSlots
             State_Diag%SatDiagnSurfFlux(I,:,:) =                            &
             State_Diag%SatDiagnSurfFlux(I,:,:) * good
          ENDDO
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagn Diagnostic for Chemical Loss
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnLoss ) THEN
          DO S = 1, State_Diag%Map_SatDiagnLoss%nSlots
             State_Diag%SatDiagnLoss(I,:,:,S) =                              &
             State_Diag%SatDiagnLoss(I,:,:,S) * good
          ENDDO
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagn Diagnostic for Chemical Production
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnProd ) THEN
          DO S = 1, State_Diag%Map_SatDiagnProd%nSlots
             State_Diag%SatDiagnProd(I,:,:,S) =                              &
             State_Diag%SatDiagnProd(I,:,:,S) * good
          ENDDO
       ENDIF

       !---------------------------------------------------------------------
       ! SatDiagn Diagnostic for Reaction Rates
       ! SatDiagnRxnRate was previously defined in fullchem_mod.F90
       !---------------------------------------------------------------------
       IF ( State_Diag%Archive_SatDiagnRxnRate ) THEN
          DO S = 1, State_Diag%Map_SatDiagnRxnRate%nSlots
             State_Diag%SatDiagnRxnRate(I,:,:,S) =                           &
             State_Diag%SatDiagnRxnRate(I,:,:,S) * good
          ENDDO
       ENDIF

    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointers for safety's sake
    Spc => NULL()

  END SUBROUTINE Do_Archive_SatDiagn
!EOC
END MODULE Diagnostics_mod
