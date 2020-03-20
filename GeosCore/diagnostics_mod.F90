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
  PUBLIC :: Compute_Column_Mass
  PUBLIC :: Compute_Budget_Diagnostics
!
! !PRIVATE MEMBER FUNCTIONS
!
  PRIVATE :: Set_SpcConc_Diagnostic
  PRIVATE :: Set_SpcConc_Diags_VVDry
!
! !PRIVATE DATA MEMBERS:
!
    CHARACTER(LEN=255), PARAMETER :: &
         ModLoc = '(in module GeosCore/diagnostics_mod.F90)'
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
  SUBROUTINE Set_Diagnostics_EndofTimestep( Input_Opt,  State_Chm, State_Diag, &
                                            State_Grid, State_Met, RC )
!
! !USES:
!
    USE Input_Opt_Mod,    ONLY : OptInput
    USE State_Met_Mod,    ONLY : MetState
    USE State_Chm_Mod,    ONLY : ChmState, Ind_
    USE State_Diag_Mod,   ONLY : DgnState
    USE State_Grid_Mod,   ONLY : GrdState
    USE PhysConstants,    ONLY : AIRMW
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
    INTEGER                 :: I, J, L, N
    REAL(fp)                :: ToPptv

    ! SAVEd scalars
    INTEGER, SAVE           :: id_Hg2 = -1
    INTEGER, SAVE           :: id_HgP = -1
    LOGICAL, SAVE           :: FIRST  = .TRUE.

    ! Strings
    CHARACTER(LEN=255)      :: ErrMsg, ThisLoc

    !=======================================================================
    ! Set_Diagnostics_EndofTimestep begins here
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at Set_Diagnostics_EndofTimestep ' // ModLoc

! NOTE: If you need to change SpeciesConc from "v/v dry" to other units,
! then uncomment this subroutine call.  Also comment out where SpeciesConc
! gets updated in routine Set_SpcConc_Diags_VVDry below.
!    !-----------------------------------------------------------------------
!    ! Set species concentration diagnostic in units specified in state_diag_mod
!    !-----------------------------------------------------------------------
!    IF ( State_Diag%Archive_SpeciesConc ) THEN
!       CALL Set_SpcConc_Diagnostic( 'SpeciesConc', State_Diag%SpeciesConc,   &
!                                    Input_Opt,  State_Chm,                   &
!                                    State_Grid, State_Met,  RC              )
!
!       ! Trap potential errors
!       IF ( RC /= GC_SUCCESS ) THEN
!          ErrMsg = 'Error encountered setting species concentration diagnostic'
!          CALL GC_ERROR( ErrMsg, RC, ThisLoc )
!       ENDIF
!    ENDIF

    !-----------------------------------------------------------------------
    ! Set species concentration for diagnostics in units of
    ! v/v dry air = mol/mol dry air
    !-----------------------------------------------------------------------
    CALL Set_SpcConc_Diags_VVDry( Input_Opt,  State_Chm, State_Diag, &
                                  State_Grid, State_Met, RC      )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered setting species concentration diagnostic'
       CALL GC_ERROR( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !-----------------------------------------------------------------------
    ! Set total dry deposition flux
    !-----------------------------------------------------------------------
    IF ( State_Diag%Archive_DryDep ) THEN
       !$OMP PARALLEL DO          &
       !$OMP DEFAULT( SHARED  )   &
       !$OMP PRIVATE( I, J, N )
       DO N = 1, State_Chm%nDryDep
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          State_Diag%DryDep(I,J,N) = State_Diag%DryDepChm(I,J,N)             &
                                   + State_Diag%DryDepMix(I,J,N)
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !-----------------------------------------------------------------------
    ! Compute fraction of time each grid box spent in the troposphere
    !-----------------------------------------------------------------------
    IF ( State_Diag%Archive_FracOfTimeInTrop ) THEN
       !$OMP PARALLEL DO            &
       !$OMP DEFAULT( SHARED      ) &
       !$OMP SCHEDULE( DYNAMIC, 8 ) &
       !$OMP PRIVATE( I, J, L )
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

    !-----------------------------------------------------------------------
    ! Diagnostics for the mercury and tagged mercury simulations
    !-----------------------------------------------------------------------
    IF ( Input_Opt%ITS_A_MERCURY_SIM ) THEN

       ! Get species indices for Hg2 and HgP
       IF ( FIRST ) THEN
          id_Hg2 = Ind_('Hg2')
          id_HgP = Ind_('HgP')
          FIRST  = .FALSE.
       ENDIF

       !--------------------------------------------
       ! Ractive gaseous mercury (RGM) [pptv]
       !--------------------------------------------
       IF ( id_Hg2 > 0 .and. State_Diag%Archive_ReactiveGaseousHg ) THEN

          ! Conversion factor to pptv
          ToPptv = ( AIRMW                                  /                &
                     State_Chm%SpcData(id_Hg2)%Info%EmMW_g  *                &
                     1.0e+12_fp                               )

          ! Save into State_diag
          State_Diag%ReactiveGaseousHg = State_Chm%Species(:,:,:,id_Hg2)     &
                                       * ToPptv
       ENDIF

       !--------------------------------------------
       ! Ractive particulate mercury (RGM) [pptv]
       !--------------------------------------------
       IF ( id_HgP > 0 .and. State_Diag%Archive_ParticulateBoundHg ) THEN

          ! Conversion factor to pptv
          ToPptv = ( AIRMW                                  /                &
                     State_Chm%SpcData(id_HgP)%Info%EmMW_g  *                &
                     1.0e+12_fp                               )

          ! Save into State_Diag
          State_Diag%ParticulateBoundHg = State_Chm%Species(:,:,:,id_HgP)    &
                                        * ToPptv
       ENDIF
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
! !DESCRIPTION: This routine is currently not used but is available if needed
!   in the future for diagnostics that are summed over a timestep and must
!   therefore be zeroed at the start of each time in the dynamic loop.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Zero_Diagnostics_StartofTimestep( State_Chm, State_Diag, RC )
!
! !USES:
!
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Diag_Mod,   ONLY : DgnState
!
! !INPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm      ! Chemistry state obj
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
    INTEGER                 :: I, J, L, N
    CHARACTER(LEN=255)      :: ErrMsg, thisLoc

    !=======================================================================
    ! Zero_Diagnostics_StartofTimestep begins here
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> Zero_Diagnostics_StartofTimestep ' // ModLoc

    ! Zero diagnostics here

  END SUBROUTINE Zero_Diagnostics_StartofTimestep
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_SpcConc_Diagnostic
!
! !DESCRIPTION: Subroutine Set_SpcConc\_Diagnostic sets the passed species
!  concentration diagnostic array stored in State_Diag to the instantaneous
!  State_Chm%Species values converted to the diagnostic unit stored in
!  the State_Diag metadata.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_SpcConc_Diagnostic( DiagMetadataID, Ptr2Data, Input_Opt, &
                                     State_Chm, State_Grid, State_Met, RC )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState, Get_Metadata_State_Diag
    USE State_Grid_Mod, ONLY : GrdState
    USE UnitConv_Mod,   ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN)  :: DiagMetadataID ! Diagnostic id
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt      ! Input Options object
    TYPE(GrdState),   INTENT(IN)  :: State_Grid     ! Grid state object
    TYPE(MetState),   INTENT(IN)  :: State_Met      ! Meteorology state object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm         ! Chemistry state obj
    REAL(f8),         POINTER       :: Ptr2Data(:,:,:,:) ! Diagnostics array
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC      ! Success or failure?
!
! !REMARKS:
!  The name argument is used to retrieve metadata (units) for the diagnostic
!  of interest. The Ptr2Data should be of form State_Diag%xxx where xxx is
!  the name of the diagnostic array to be set.
!
!  This routine allows the freedom to easily create multiple species
!  concentration diagnostics other than the default end-of-timestep
!  diagnostic State_Diag%SpeciesConc, although this routine is used to set
!  that as well.
!
!  For example, you may create diagnostics for concentrations at different
!  phases of the GEOS-Chem run by adding new diagnostic arrays, e.g.
!  State_Diag%SpeciesConc_preChem and State_Diag%SpeciesConc_postChem, to
!  state_diag_mod.F90. Metadata could be used for 'SpeciesConc' or a new
!  metadata entry could be created.
!
!  Changing the unit string of the metadata in state_diag_mod.F90 will
!  result in different units in the species concencentration diagnostic
!  output. The units in the netcdf file metadata will reflect the new units.
!
! !REVISION HISTORY:
!  27 Sep 2017 - E. Lundgren - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, Units, OrigUnit
    LOGICAL            :: Found
    INTEGER            :: I, J, L, N

    !====================================================================
    ! Set_SpcConc_Diagnostic begins here!
    !====================================================================

    ! Assume success
    RC      =  GC_SUCCESS
    Found   = .FALSE.
    ThisLoc = ' -> Set_SpcConc_Diagnostic ' // ModLoc

    ! Exit if species concentration is not a diagnostics in HISTORY.rc
    IF ( ASSOCIATED( Ptr2Data ) ) THEN

       IF ( TRIM( DiagMetadataID ) == 'SpeciesRst' .or. &
            TRIM( DiagMetadataID ) == 'SpeciesBC' ) THEN

          ! For GEOS-Chem restart and BC collections force units to v/v dry
          Units = 'v/v dry'

       ELSE

          ! Retrieve the units of the diagnostic from the metadata
          CALL Get_Metadata_State_Diag( Input_Opt%amIRoot,     &
                                        TRIM(DiagMetadataID),  &
                                        Found, RC, Units=Units )

          ! Allow for alternate format of units
          IF ( TRIM(Units) == 'mol mol-1 dry' ) Units = 'v/v dry'
          IF ( TRIM(Units) == 'kg kg-1 dry'   ) Units = 'kg/kg dry'
          IF ( TRIM(Units) == 'kg m-2'        ) Units = 'kg/m2'
          IF ( TRIM(Units) == 'molec cm-3'    ) Units = 'molec/cm3'

       ENDIF

       ! Convert State_Chm%Species unit to diagnostic units
       CALL Convert_Spc_Units( Input_Opt, State_Chm, State_grid, State_Met, &
                               Units, RC, OrigUnit=OrigUnit )

       ! Copy species concentrations to diagnostic array
       !$OMP PARALLEL DO           &
       !$OMP DEFAULT( SHARED     ) &
       !$OMP PRIVATE( I, J, L, N )
       DO N = 1, State_Chm%nSpecies
       DO L = 1, State_Grid%NZ
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          Ptr2Data(I,J,L,N) = State_Chm%Species(I,J,L,N)
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       ! Convert State_Chm%Species back to original unit
       CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                               OrigUnit, RC )

       ! Error handling
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error converting species units for archiving diagnostics'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
       ENDIF
    ENDIF

  END SUBROUTINE Set_SpcConc_Diagnostic
!EOC
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
  SUBROUTINE Set_SpcConc_Diags_VVDry( Input_Opt,  State_Chm, State_Diag, &
                                      State_Grid, State_Met, RC )
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
    USE UnitConv_Mod,   ONLY : Convert_Spc_Units
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
! !REMARKS:
!  This subroutine was written so as to minimize the number of unit
!  conversions that occur per call (which happens once per timestep).
!  Units  are now converted to and from "v/v dry air" only once per call.
!  The prior algorithm, which used routine Set_SpcConc_Diagnostic,
!  was converting units 2 or 3 times per call, which can make run times
!  substantially longer.
!
!  The State_Diag%SpeciesConc diagnostic has units of "mol mol-1 dry",
!  which is equivalent to "v/v dry".  Therefore, we can include the
!  State_Diag%SpeciesConc diagnostic in this routine.  But if you change
!  the units of State_Diag%SpeciesConc, you should instead comment it out
!  below and call routine Set_SpcConc_Diagnostic instead.  This will
!  ensure that State_Diag%SpeciesConc will get set to the proper units.
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
    LOGICAL            :: Found
    INTEGER            :: D, I, J, L, N
    REAL(fp)           :: TmpVal, Conv

    ! Strings
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc, Units, OrigUnit

    !====================================================================
    ! Set_SpcConc_Diags_VVDry begins here!
    !====================================================================

    ! Assume success
    RC      =  GC_SUCCESS
    Found   = .FALSE.
    ThisLoc = ' -> Set_SpcConc_Diagnostics ' // ModLoc

    ! We a ssume all diagnostics are already in [v/v dry]
    ! This will allow us to minimize unit conversions
    Units   = 'v/v dry'

    ! Convert State_Chm%Species unit to [v/v dry]
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            Units, RC, OrigUnit=OrigUnit )

    ! Error handling
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error converting species units for archiving diagnostics #1'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Copy species concentrations to diagnostic arrays [v/v dry]
    !=======================================================================
    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( I, J, L, N )
    DO N = 1, State_Chm%nSpecies
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Boundary conditions for nested grid [v/v dry]
       IF ( State_Diag%Archive_SpeciesBC ) THEN
          State_Diag%SpeciesBC(I,J,L,N) = State_Chm%Species(I,J,L,N)
       ENDIF

       ! Species concentrations diagnostic [v/v dry]
       ! NOTE: If you change the units of SpeciesConc in state_diag_mod.F90,
       ! then comment this IF block out and then also uncomment the IF block
       ! in the main routine above where Set_SpcConc_Diagnostic is called.
       IF ( State_Diag%Archive_SpeciesConc ) THEN
          State_Diag%SpeciesConc(I,J,L,N) = State_Chm%Species(I,J,L,N)
       ENDIF

       ! Species concentrations for restart file [v/v dry]
       IF ( State_Diag%Archive_SpeciesRst ) THEN
          State_Diag%SpeciesRst(I,J,L,N) = State_Chm%Species(I,J,L,N)
       ENDIF

    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

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
    ! We recalculate Ra using Z1 using a value specified in input.geos;
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
             TmpVal = State_Chm%Species(I,J,1,N)

             ! Conversion factor used to translate from
             ! lowest model layer (~60m) to the surface
             Conv = ( 1.0_fp                                              &
                  -   ( State_Diag%DryDepRaALT1(I,J) / 100.0_fp )         &
                  *   State_Diag%DryDepVelForALT1(I,J,D)                 )

             ! Do not let CONV go negative
             IF ( Conv < 0.0_fp ) Conv = 1.0_fp

             ! Save concentration at the user-defined altitude
             ! as defined in input.geos (usually 10m).
             State_Diag%SpeciesConcALT1(I,J,D) = TmpVal * Conv

          ENDDO
          ENDDO
       ENDDO
       !$OMP END PARALLEL DO

    ENDIF

    ! Convert State_Chm%Species back to original unit
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit,  RC )

    ! Error handling
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error converting species units for archiving diagnostics #2'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Set_SpcConc_Diags_VVDry
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute_Column_Mass
!
! !DESCRIPTION: Subroutine Compute\_Column\_Mass calculates the
!  initial or final mass for a given region of the column for use in the
!  calculation of the budget diagnostics.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Column_Mass( Input_Opt, State_Chm, State_Grid, State_Met, &
                                  SpcMap,    isFull,    isTrop,     isPBL,     &
                                  ColMass,    RC      )
!
! !USES:
!
    USE Input_Opt_Mod,  Only : OptInput
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Grid_Mod, ONLY : GrdState
!ewl    USE UnitConv_Mod,   ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt        ! Input options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid       ! Grid state object
    TYPE(MetState), INTENT(IN)    :: State_Met        ! Meteorology state object
    INTEGER,        POINTER       :: SpcMap(:)        ! Map to species indexes
    LOGICAL,        INTENT(IN)    :: isFull           ! True if full col diag on
    LOGICAL,        INTENT(IN)    :: isTrop           ! True if trop col diag on
    LOGICAL,        INTENT(IN)    :: isPBL            ! True if PBL col diag on
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm        ! Chemistry state obj
    REAL(f8),       POINTER       :: colMass(:,:,:,:) ! column masses
                                                      ! (I,J,spc,col region)
                                                      ! 1:full, 2:trop, 3:pbl
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC                 ! Success or failure?
!
! !REMARKS:
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
    CHARACTER(LEN=255)  :: ErrMsg, ThisLoc
    INTEGER             :: I, J, L, M, N, numSpc, region, PBL_TOP
    REAL*8, ALLOCATABLE :: SpcMass(:,:,:,:)

    !====================================================================
    ! Compute_Column_Mass begins here!
    !====================================================================

    ! Initialize
    RC      =  GC_SUCCESS
    ThisLoc = ' -> Compute_Column_Mass ' // ModLoc
    numSpc = SIZE(SpcMap)
    colMass = 0.0_f8

    ! Get concentrations in units of kg. Incoming units should be kg/kg dry.
    IF (State_Chm%Spc_Units == 'kg/kg dry' ) THEN
       ALLOCATE(SpcMass(State_Grid%NX, State_Grid%NY, State_Grid%NZ, numSpc))
    ELSE
       CALL GC_Error( 'State_Chm%Species units must be kg/kg dry. '// &
                      'Incorrect units: '//TRIM(State_Chm%Spc_Units),  &
                      RC, ThisLoc )
      RETURN
    ENDIF
    !$OMP PARALLEL DO        &
    !$OMP DEFAULT( SHARED )  &
    !$OMP PRIVATE( I, J, L, M, N )
    DO M = 1, numSpc
    DO L = 1, State_Grid%NZ
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX
       N = SpcMap(M)
       SpcMass(I,J,L,M) = State_Chm%Species(I,J,L,N) * State_Met%AD(I,J,L)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Full column
    IF ( isFull ) THEN
       region = 1
       !$OMP PARALLEL DO        &
       !$OMP DEFAULT( SHARED )  &
       !$OMP PRIVATE( I, J, M, N )
       DO M = 1, numSpc
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          N = SpcMap(M)
          colMass(I,J,N,region) = SUM(SpcMass(I,J,:,M))
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    ! Troposphere
    IF ( isTrop ) THEN
       region = 2
       !$OMP PARALLEL DO        &
       !$OMP DEFAULT( SHARED )  &
       !$OMP PRIVATE( I, J, M, N )
       DO M = 1, numSpc
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          N = SpcMap(M)
          colMass(I,J,N,region) = SUM(SpcMass(I,J,1:State_Met%TropLev(I,J),M))
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    ! PBL
    IF ( isPBL ) THEN
       region = 3
       !$OMP PARALLEL DO        &
       !$OMP DEFAULT( SHARED )  &
       !$OMP PRIVATE( I, J, M, N )
       DO M = 1, numSpc
       DO J = 1, State_Grid%NY
       DO I = 1, State_Grid%NX
          N = SpcMap(M)
          PBL_TOP = MAX( 1, FLOOR( State_Met%PBL_TOP_L(I,J) ) )
          colMass(I,J,N,region) = SUM(SpcMass(I,J,1:PBL_TOP,M))
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    ! Clean up
    DEALLOCATE(SpcMass)

  END SUBROUTINE Compute_Column_Mass
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Compute_Budget_Diagnostics
!
! !DESCRIPTION: Subroutine Compute\_Budget\_Diagnostics calculates the
!  budget diagnostics for a given component by taking the difference of the
!  final and initial kg per grid cell and dividing by the timestep in seconds
!  to get kg/s.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Budget_Diagnostics( State_Grid,               &
                                         SpcMap,       TS,         &
                                         isFull,       isTrop,     &
                                         isPBL,        diagFull,   &
                                         diagTrop,     diagPBL,    &
                                         mass_initial, mass_final, &
                                         RC )
!
! !USES:
!
    USE State_Grid_Mod, ONLY : GrdState
!
! !INPUT PARAMETERS:
!
    TYPE(GrdState), INTENT(IN)  :: State_Grid      ! Grid State object
    INTEGER,        POINTER     :: SpcMap(:)       ! Map to species indexes
    REAL(fp),       INTENT(IN)  :: TS              ! timestep [s]
    LOGICAL,        INTENT(IN)  :: isFull          ! True if full col diag on
    LOGICAL,        INTENT(IN)  :: isTrop          ! True if trop col diag on
    LOGICAL,        INTENT(IN)  :: isPBL           ! True if PBL col diag on
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f8),       TARGET      :: diagFull(:,:,:)       ! ptr to full col diag
    REAL(f8),       TARGET      :: diagTrop(:,:,:)       ! ptr to trop col diag
    REAL(f8),       TARGET      :: diagPBL(:,:,:)        ! ptr to pbl col diag
    REAL(f8),       POINTER     :: mass_initial(:,:,:,:) ! ptr to initial mass
    REAL(f8),       POINTER     :: mass_final(:,:,:,:)   ! ptr to final mass
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT) :: RC              ! Success or failure?
!
! !REMARKS:
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
    INTEGER            :: I, J, M, N, R, numSpc, numRegions
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    LOGICAL            :: setDiag
    REAL(f8), POINTER  :: ptr3d(:,:,:)

    !====================================================================
    ! Compute_Budget_Diagnostics begins here!
    !====================================================================

    ! Initialize
    RC         =  GC_SUCCESS
    ThisLoc    = ' -> Compute_Budget_Diagnostics ' // ModLoc
    numSpc     = SIZE(SpcMap)
    numRegions = 3
    ptr3d      => NULL()

    ! Loop over regions and only set diagnostic for those that are on
    DO R = 1, numRegions
       setDiag = .FALSE.
       SELECT CASE ( R )
          CASE ( 1 )
             ptr3d => diagFull
             setDiag = isFull
          CASE ( 2 )
             ptr3d => diagTrop
             setDiag = isTrop
          CASE ( 3 )
             ptr3d => diagPBL
             setDiag = isPBL
          CASE DEFAULT
             ErrMsg = 'Region not defined for budget diagnostics'
             CALL GC_Error( ErrMsg, RC, ThisLoc )
             RETURN
       END SELECT

       ! Compute diagnostics as [kg/s] by taking dividing the mass
       ! difference by component dt in seconds
       !
       ! NOTE: if changing the definition of budget diagnostics below be sure
       ! to also update the budget diagnostic metadata in state_diag_mod.F90
       ! within subroutine Get_Metadata_State_Diag. If you wish to output
       ! different units, e.g. kg/m2/s instead of kg/s, update the metadata
       ! and also the unit string in the unit conversion call above in subroutine
       ! Compute_Column_Mass. (ewl, 9/26/18)
       IF ( setDiag ) THEN
          !$OMP PARALLEL DO        &
          !$OMP DEFAULT( SHARED )  &
          !$OMP PRIVATE( I, J, M, N )
          DO M = 1, numSpc
          DO J = 1, State_Grid%NY
          DO I = 1, State_Grid%NX
             N  = SpcMap(M)
             ptr3d(I,J,M) =   &
                   ( mass_final(I,J,N,R) - mass_initial(I,J,N,R) ) / TS
          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO
       ENDIF

       ! Free pointer
       ptr3d => NULL()

    ENDDO

    ! Zero the mass arrays now that diagnostics are set
    mass_initial = 0.0_f8
    mass_final   = 0.0_f8

  END SUBROUTINE Compute_Budget_Diagnostics
!EOC
END MODULE Diagnostics_mod
