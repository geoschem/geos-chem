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
  PUBLIC :: Init_Diagnostics_Mod
  PUBLIC :: Set_Diagnostics_EndofTimestep
  PUBLIC :: Zero_Diagnostics_StartofTimestep
!
! !PRIVATE MEMBER FUNCTIONS
!
  PRIVATE :: Set_SpcConc_Diagnostic
!
!
! !PRIVATE DATA MEMBERS:
!
#if defined( NC_DIAG )
      LOGICAL            :: Archive_SpeciesConc
      LOGICAL            :: Archive_DryDep
      CHARACTER(LEN=255) :: ModLoc
#endif
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
! !IROUTINE: Init_Diagnostics_Mod
!
! !DESCRIPTION:
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Diagnostics_Mod ( am_I_Root, State_Diag, RC )
!
! !USES:
!
    USE State_Diag_Mod, ONLY : DgnState
!
! !INPUT PARAMETERS:
!
    LOGICAL,         INTENT(IN)  :: am_I_Root
    TYPE(DgnState),  INTENT(IN)  :: State_Diag     ! Diagnostics state obj
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT) :: RC
!
! !REVISION HISTORY:
!  01 Feb 2018 - E. Lundgren - initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( NC_DIAG )
    CHARACTER(LEN=255)      :: ErrMsg, thisLoc

    !=======================================================================
    ! Init_Diagnostics_Mod begins here
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ModLoc = '(GeosCore/diagnostics_mod.F90)'
    ThisLoc = ' -> Init_Diagnostics_mod ' // ModLoc
    
    ! Determine whether diagnostics will be archived
    Archive_SpeciesConc = .FALSE.
    Archive_DryDep = .FALSE.
    IF ( ASSOCIATED( State_Diag%SpeciesConc ) ) Archive_SpeciesConc = .TRUE. 
    IF ( ASSOCIATED( State_Diag%DryDep ) ) THEN 
       IF ( ASSOCIATED( State_Diag%DryDepChm ) .AND. &
            ASSOCIATED( State_Diag%DryDepMix ) ) THEN
          Archive_DryDep = .TRUE.
       ELSE        
          ErrMsg = 'State_Diag%DryDepChm or State_Diag%DryDepMix must be ' &
                   // 'set to output the total dry deposition flux diagnostics'
          CALL GC_ERROR( ErrMsg, RC, ThisLoc )
       ENDIF
    ENDIF

#endif
  END SUBROUTINE Init_Diagnostics_Mod
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Diagnostics_EndofTimestep
!
! !DESCRIPTION:
!\\
!\\
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
    USE State_Chm_Mod,    ONLY : ChmState
    USE State_Diag_Mod,   ONLY : DgnState
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
#if defined( NC_DIAG )
    INTEGER                 :: I, J, L, N
    CHARACTER(LEN=255)      :: ErrMsg, thisLoc

    !=======================================================================
    ! Set_Diagnostics_EndofTimestep begins here
    !=======================================================================

    ! Initialize
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> Set_Diagnostics_EndofTimestep ' // ModLoc

    ! Set species concentration diagnostic in units of mol/mol dry air
    IF ( Archive_SpeciesConc ) THEN
       CALL Set_SpcConc_Diagnostic( am_I_Root, 'SpeciesConc',           &
                                    State_Diag%SpeciesConc,             &
                                    Input_Opt, State_Met, State_Chm, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered setting species concentration diagnostic'
          CALL GC_ERROR( ErrMsg, RC, ThisLoc )
       ENDIF
    ENDIF

    ! Set total dry deposition flux
    IF ( Archive_DryDep ) THEN
       !$OMP PARALLEL DO          &
       !$OMP DEFAULT( SHARED  )   &
       !$OMP PRIVATE( I, J, L )
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
       DO N = 1, State_Chm%nDryDep
          State_Diag%DryDep(I,J,L,N) = State_Diag%DryDepChm(I,J,L,N) &
                                     + State_Diag%DryDepMix(I,J,L,N)
       ENDDO
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

#endif
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
  SUBROUTINE Zero_Diagnostics_StartofTimestep ( am_I_Root,  State_Chm, &
                                                State_Diag, RC )
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
#if defined( NC_DIAG )
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

#endif
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
  SUBROUTINE Set_SpcConc_Diagnostic( am_I_Root, DiagMetadataID,       &
                                     Ptr2Data,  Input_Opt, State_Met, &
                                     State_Chm, RC ) 
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( NC_DIAG )
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

       ! Retrieve the units of the diagnostic from the metadata
       CALL Get_Metadata_State_Diag( am_I_Root, TRIM(DiagMetadataID), &
                                     Found, RC, Units=Units )

       ! Allow for alternate format of units
       IF ( TRIM(Units) == 'mol mol-1 dry' ) Units = 'v/v dry'
       IF ( TRIM(Units) == 'kg kg-1 dry'   ) Units = 'kg/kg dry'
       IF ( TRIM(Units) == 'kg m-2'        ) Units = 'kg/m2'
       IF ( TRIM(Units) == 'molec cm-3'    ) Units = 'molec/cm3'
       
       ! Convert State_Chm%Species unit to diagnostic units
       CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, State_Chm, &
                               Units, RC, OrigUnit=OrigUnit )
       
       ! Copy species concentrations to diagnostic array
       !$OMP PARALLEL DO          &
       !$OMP DEFAULT( SHARED  )   &
       !$OMP PRIVATE( I, J, L, N )
       DO L = 1, LLPAR
       DO J = 1, JJPAR
       DO I = 1, IIPAR
       DO N = 1, State_Chm%nSpecies
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

#endif
  END SUBROUTINE Set_SpcConc_Diagnostic
!EOC
END MODULE Diagnostics_mod
