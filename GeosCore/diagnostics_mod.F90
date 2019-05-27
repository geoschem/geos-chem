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
  USE CMN_Size_Mod
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
!
!
! !PRIVATE DATA MEMBERS:
!
    CHARACTER(LEN=255), PARAMETER :: &
         ModLoc = '(in module GeosCore/diagnostics_mod.F90)'
!
! !REVISION HISTORY:
!  01 Feb 2018 - E. Lundgren - Initial version
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
  SUBROUTINE Set_Diagnostics_EndofTimestep ( am_I_Root,  Input_Opt,  &
                                             State_Met,  State_Chm,  &
                                             State_Diag, RC )
!
! !USES:
!
    USE Input_Opt_Mod,    ONLY : OptInput
    USE State_Met_Mod,    ONLY : MetState
    USE State_Chm_Mod,    ONLY : ChmState, Ind_
    USE State_Diag_Mod,   ONLY : DgnState
    USE PhysConstants,    ONLY : AIRMW
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt      ! Input Options object
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

    !-----------------------------------------------------------------------
    ! Set species concentration for restart in units of mol/mol dry air
    !-----------------------------------------------------------------------
    IF ( State_Diag%Archive_SpeciesRst ) THEN
       CALL Set_SpcConc_Diagnostic( am_I_Root, 'SpeciesRst',                &
                                    State_Diag%SpeciesRst,                  &
                                    Input_Opt,  State_Met,                   &
                                    State_Chm,  RC                          )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered setting species concentration diagnostic'
          CALL GC_ERROR( ErrMsg, RC, ThisLoc )
       ENDIF
    ENDIF
    
    !-----------------------------------------------------------------------
    ! Set species concentration diagnostic in units specified in state_diag_mod
    !-----------------------------------------------------------------------
    IF ( State_Diag%Archive_SpeciesConc ) THEN
       CALL Set_SpcConc_Diagnostic( am_I_Root, 'SpeciesConc',                &
                                    State_Diag%SpeciesConc,                  &
                                    Input_Opt,  State_Met,                   &
                                    State_Chm,  RC                          )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered setting species concentration diagnostic'
          CALL GC_ERROR( ErrMsg, RC, ThisLoc )
       ENDIF
    ENDIF

    !-----------------------------------------------------------------------
    ! Set total dry deposition flux
    !-----------------------------------------------------------------------
    IF ( State_Diag%Archive_DryDep ) THEN
       !$OMP PARALLEL DO          &
       !$OMP DEFAULT( SHARED  )   &
       !$OMP PRIVATE( I, J, N )
       DO N = 1, State_Chm%nDryDep
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          State_Diag%DryDep(I,J,N) = State_Diag%DryDepChm(I,J,N)             &
                                   + State_Diag%DryDepMix(I,J,N)
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
  SUBROUTINE Zero_Diagnostics_StartofTimestep( am_I_Root,  State_Chm,        &
                                               State_Diag, RC               )
!
! !USES:
!
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Diag_Mod,   ONLY : DgnState
!
! !INPUT PARAMETERS:
!
    LOGICAL,          INTENT(IN)    :: am_I_Root
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
  SUBROUTINE Set_SpcConc_Diagnostic( am_I_Root, DiagMetadataID, Ptr2Data,    &
                                     Input_Opt, State_Met,      State_Chm,   &
                                     RC                                     ) 
!
! !USES:
!
    USE Input_Opt_Mod,  ONLY : OptInput
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState, Get_Metadata_State_Diag
    USE UnitConv_Mod,   ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS: 
!
    LOGICAL,          INTENT(IN)  :: am_I_Root      ! Are we on the root CPU?
    CHARACTER(LEN=*), INTENT(IN)  :: DiagMetadataID ! Diagnostic id
    TYPE(OptInput),   INTENT(IN)  :: Input_Opt      ! Input Options object
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
!  06 Nov 2018 - M. Sulprizio- Only allow for different units in GCHP and GEOS-5
!                              and force units to v/v for GEOS-Chem Classic to
!                              avoid issues with the GEOS-Chem restart file
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

       IF ( TRIM( DiagMetadataID ) == 'SpeciesRst' ) THEN

          ! For GEOS-Chem restart collection force units to v/v dry
          Units = 'v/v dry'

       ELSE

          ! Retrieve the units of the diagnostic from the metadata
          CALL Get_Metadata_State_Diag( am_I_Root, TRIM(DiagMetadataID), &
                                        Found, RC, Units=Units )

          ! Allow for alternate format of units
          IF ( TRIM(Units) == 'mol mol-1 dry' ) Units = 'v/v dry'
          IF ( TRIM(Units) == 'kg kg-1 dry'   ) Units = 'kg/kg dry'
          IF ( TRIM(Units) == 'kg m-2'        ) Units = 'kg/m2'
          IF ( TRIM(Units) == 'molec cm-3'    ) Units = 'molec/cm3'

       ENDIF

       ! Convert State_Chm%Species unit to diagnostic units
       CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, State_Chm, &
                               Units, RC, OrigUnit=OrigUnit )
       
       ! Copy species concentrations to diagnostic array
       !$OMP PARALLEL DO           &
       !$OMP DEFAULT( SHARED     ) &
       !$OMP PRIVATE( I, J, L, N )
       DO N = 1, State_Chm%nSpecies
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          Ptr2Data(I,J,L,N) = State_Chm%Species(I,J,L,N)
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO

       ! Convert State_Chm%Species back to original unit
       CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                               State_Chm, OrigUnit, RC )
       
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
! !IROUTINE: Compute_Column_Mass
!
! !DESCRIPTION: Subroutine Compute\_Column\_Mass calculates the
!  initial or final mass for a given region of the column for use in the 
!  calculation of the budget diagnostics.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Compute_Column_Mass( am_I_Root,  Input_Opt,  State_Met,       &
                                  State_Chm,  SpcMap,     isFull,          &
                                  isTrop,     isPBL,      colMass,         &
                                  RC      ) 
!
! !USES:
!
    USE Input_Opt_Mod,  Only : OptInput
    USE State_Met_Mod,  ONLY : MetState
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Diag_Mod, ONLY : DgnState
    USE UnitConv_Mod,   ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root        ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt        ! Input options object
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    CHARACTER(LEN=63)  :: OrigUnit
    CHARACTER(LEN=255) :: ErrMsg, ThisLoc
    INTEGER            :: I, J, M, N, numSpc, region

    !====================================================================
    ! Compute_Column_Mass begins here!
    !====================================================================

    ! Initialize
    RC      =  GC_SUCCESS
    ThisLoc = ' -> Compute_Column_Mass ' // ModLoc
    numSpc = SIZE(SpcMap)

    ! Convert species units to kg
    ! NOTE: If you wish to output different units, e.g. kg/m2/s instead of 
    ! kg/s, update the unit string passed in the conversion subroutine below.
    ! Removing the per seconds would not be done here, but in the subroutine
    ! that follows to compute the budget diagnostic. Also update the budget 
    ! diagnostic units metadata in state_diag_mod.F90 within subroutine 
    ! Get_Metadata_State_Diag (ewl, 9/26/18)
    CALL Convert_Spc_Units( am_I_Root,  Input_Opt, State_Met,   &
                            State_Chm,  'kg',      RC,          &
                            OrigUnit=OrigUnit                  )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( 'Unit conversion error', RC, ThisLoc )
       RETURN
    ENDIF

    ! Start with all zeros
    colMass = 0.0_f8

    ! Full column
    IF ( isFull ) THEN
       region = 1
       !$OMP PARALLEL DO        &
       !$OMP DEFAULT( SHARED )  &
       !$OMP PRIVATE( I, J, M, N )
       DO M = 1, numSpc
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          N = SpcMap(M)
          colMass(I,J,N,region) =    &
                  SUM( State_Chm%Species(I,J,:,N) )
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
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          N = SpcMap(M)
          colMass(I,J,N,region) =    &
             SUM( State_Chm%Species(I,J,1:State_Met%TropLev(I,J),N) )
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
       DO J = 1, JJPAR
       DO I = 1, IIPAR
          N = SpcMap(M)
          colMass(I,J,N,region) =    &
             SUM( State_Chm%Species(I,J,1:State_Met%PBL_TOP_L(I,J),N) )
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    ! Convert species conc back to original units
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met,  &
                            State_Chm, OrigUnit,  RC         )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error( 'Unit conversion error', RC, ThisLoc )
       RETURN
    ENDIF

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
  SUBROUTINE Compute_Budget_Diagnostics( am_I_Root,    SpcMap,     TS,       &
                                         isFull,       isTrop,     isPBL,    &
                                         diagFull,     diagTrop,   diagPBL,  & 
                                         mass_initial, mass_final, RC        ) 
!
! !USES:
!
!
! !INPUT PARAMETERS: 
!
    LOGICAL,  INTENT(IN)  :: am_I_Root       ! Are we on the root CPU?
    INTEGER,  POINTER     :: SpcMap(:)       ! Map to species indexes
    REAL(fp), INTENT(IN)  :: TS              ! timestep [s]
    LOGICAL,  INTENT(IN)  :: isFull          ! True if full col diag on
    LOGICAL,  INTENT(IN)  :: isTrop          ! True if trop col diag on
    LOGICAL,  INTENT(IN)  :: isPBL           ! True if PBL col diag on
!
! !INPUT/OUTPUT PARAMETERS:
!
    REAL(f8), TARGET      :: diagFull(:,:,:)       ! ptr to full col diag
    REAL(f8), TARGET      :: diagTrop(:,:,:)       ! ptr to trop col diag 
    REAL(f8), TARGET      :: diagPBL(:,:,:)        ! ptr to pbl col diag
    REAL(f8), POINTER     :: mass_initial(:,:,:,:) ! ptr to initial mass
    REAL(f8), POINTER     :: mass_final(:,:,:,:)   ! ptr to final mass
!
! !OUTPUT PARAMETERS:
!
    INTEGER,  INTENT(OUT) :: RC              ! Success or failure?
!
! !REMARKS:
!
! !REVISION HISTORY: 
!  28 Aug 2018 - E. Lundgren - Initial version
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
          DO J = 1, JJPAR
          DO I = 1, IIPAR
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
