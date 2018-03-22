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
  ! Diagnostic switches
  LOGICAL            :: Archive_SpeciesConc
  LOGICAL            :: Archive_DryDep
  LOGICAL            :: Archive_AerMass
  LOGICAL            :: Archive_AerMassASOA
  LOGICAL            :: Archive_AerMassBC
  LOGICAL            :: Archive_AerMassINDIOL
  LOGICAL            :: Archive_AerMassISN1OA
  LOGICAL            :: Archive_AerMassISOA
  LOGICAL            :: Archive_AerMassLVOCOA
  LOGICAL            :: Archive_AerMassNH4
  LOGICAL            :: Archive_AerMassNIT
  LOGICAL            :: Archive_AerMassOPOA
  LOGICAL            :: Archive_AerMassPOA
  LOGICAL            :: Archive_AerMassSAL
  LOGICAL            :: Archive_AerMassSO4
  LOGICAL            :: Archive_AerMassSOAGX
  LOGICAL            :: Archive_AerMassSOAIE
  LOGICAL            :: Archive_AerMassSOAME
  LOGICAL            :: Archive_AerMassSOAMG
  LOGICAL            :: Archive_AerMassTSOA
  LOGICAL            :: Archive_BetaNO
  LOGICAL            :: Archive_PM25
  LOGICAL            :: Archive_TotalOA
  LOGICAL            :: Archive_TotalOC
  LOGICAL            :: Archive_TotalBiogenicOA
  LOGICAL            :: Is_BetaNO
  LOGICAL            :: Is_POA

  ! Species ID flags
  INTEGER            :: id_POA1
  INTEGER            :: id_POA2
  INTEGER            :: id_SOAIE
  INTEGER            :: id_SOAME
  INTEGER            :: id_INDIOL
  INTEGER            :: id_SOAGX
  INTEGER            :: id_SOAMG
  INTEGER            :: id_LVOCOA
  INTEGER            :: id_ISN1OA

  ! Conversionf factors to ugC/m3 for Total Organic Carbon diagnostic
  REAL(fp)           :: Fac_INDIOL
  REAL(fp)           :: Fac_ISN1OA
  REAL(fp)           :: Fac_LVOCOA
  REAL(fp)           :: Fac_SOAGX
  REAL(fp)           :: Fac_SOAIE
  REAL(fp)           :: Fac_SOAME
  REAL(fp)           :: Fac_SOAMG

  ! Strings
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
    ModLoc  = '(in module GeosCore/diagnostics_mod.F90)'
    ThisLoc = ' -> at Init_Diagnostics_mod ' // ModLoc
    
    !=======================================================================
    ! Determine whether SpeciesConc diagnostics will be archived
    !=======================================================================
    Archive_SpeciesConc = ASSOCIATED( State_Diag%SpeciesConc )

    !=======================================================================
    ! Determine whether the total DryDep diagnostic will be archived
    !=======================================================================
    Archive_DryDep = .FALSE.
    IF ( ASSOCIATED( State_Diag%DryDep ) ) THEN 
       IF ( ASSOCIATED( State_Diag%DryDepChm ) .AND. &
            ASSOCIATED( State_Diag%DryDepMix ) ) THEN
          Archive_DryDep = .TRUE.
       ELSE        
          ErrMsg = 'State_Diag%DryDepChm or State_Diag%DryDepMix must be '// &
                   'set to output the total dry deposition flux diagnostics'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
       ENDIF
    ENDIF

    !=======================================================================
    ! Determine whether the aerosol mass diagnostics will be archived
    !=======================================================================
    Archive_AerMassASOA     = ASSOCIATED( State_Diag%AerMassASOA            )
    Archive_AerMassBC       = ASSOCIATED( State_Diag%AerMassBC              )
    Archive_AerMassINDIOL   = ASSOCIATED( State_Diag%AerMassINDIOL          )
    Archive_AerMassISN1OA   = ASSOCIATED( State_Diag%AerMassISN1OA          )
    Archive_AerMassISOA     = ASSOCIATED( State_Diag%AerMassISOA            )
    Archive_AerMassLVOCOA   = ASSOCIATED( State_Diag%AerMassLVOCOA          )
    Archive_AerMassNH4      = ASSOCIATED( State_Diag%AerMassNH4             )
    Archive_AerMassNIT      = ASSOCIATED( State_Diag%AerMassNIT             )
    Archive_AerMassOPOA     = ASSOCIATED( State_Diag%AerMassOPOA            )
    Archive_AerMassPOA      = ASSOCIATED( State_Diag%AerMassPOA             )
    Archive_AerMassSAL      = ASSOCIATED( State_Diag%AerMassSAL             )
    Archive_AerMassSO4      = ASSOCIATED( State_Diag%AerMassSO4             )
    Archive_AerMassSOAGX    = ASSOCIATED( State_Diag%AerMassSOAGX           )
    Archive_AerMassSOAIE    = ASSOCIATED( State_Diag%AerMassSOAIE           )
    Archive_AerMassSOAME    = ASSOCIATED( State_Diag%AerMassSOAME           )
    Archive_AerMassSOAMG    = ASSOCIATED( State_Diag%AerMassSOAMG           )
    Archive_AerMassTSOA     = ASSOCIATED( State_Diag%AerMassTSOA            )
    Archive_BetaNO          = ASSOCIATED( State_Diag%BetaNO                 )
    Archive_PM25            = ASSOCIATED( State_Diag%PM25                   )
    Archive_TotalOA         = ASSOCIATED( State_Diag%TotalOA                )
    Archive_TotalOC         = ASSOCIATED( State_Diag%TotalOC                )
    Archive_TotalBiogenicOA = ASSOCIATED( State_Diag%TotalBiogenicOA        )

    ! Set a flag to denote if we need to call the routine to archive
    ! the aerosol mass diagnostics.  True if any of the above are true.
    Archive_AerMass         = ( Archive_AerMassASOA     .or.                 &
                                Archive_AerMassBC       .or.                 &
                                Archive_AerMassINDIOL   .or.                 &
                                Archive_AerMassISN1OA   .or.                 &
                                Archive_AerMassISOA     .or.                 &
                                Archive_AerMassLVOCOA   .or.                 &
                                Archive_AerMassNH4      .or.                 &
                                Archive_AerMassNIT      .or.                 &
                                Archive_AerMassOPOA     .or.                 &
                                Archive_AerMassPOA      .or.                 &
                                Archive_AerMassSAL      .or.                 &
                                Archive_AerMassSO4      .or.                 &
                                Archive_AerMassSOAGX    .or.                 &
                                Archive_AerMassSOAIE    .or.                 &
                                Archive_AerMassSOAME    .or.                 &
                                Archive_AerMassSOAMG    .or.                 &
                                Archive_AerMassTSOA     .or.                 &
                                Archive_BetaNO          .or.                 &
                                Archive_PM25            .or.                 &
                                Archive_TotalOA         .or.                 &
                                Archive_TotalOC         .or.                 &
                                Archive_TotalBiogenicOA                     )

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
    ThisLoc = ' -> at Set_Diagnostics_EndofTimestep ' // ModLoc

    !-----------------------------------------------------------------------
    ! Set species concentration diagnostic in units of mol/mol dry air
    !-----------------------------------------------------------------------
    IF ( Archive_SpeciesConc ) THEN
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
    IF ( Archive_DryDep ) THEN
       !$OMP PARALLEL DO          &
       !$OMP DEFAULT( SHARED  )   &
       !$OMP PRIVATE( I, J )
       DO J = 1, JJPAR
       DO I = 1, IIPAR
       DO N = 1, State_Chm%nDryDep
          State_Diag%DryDep(I,J,N) = State_Diag%DryDepChm(I,J,N)         &
                                     + State_Diag%DryDepMix(I,J,N)
       ENDDO
       ENDDO
       ENDDO
       !$OMP END PARALLEL DO
    ENDIF

    !-----------------------------------------------------------------------
    ! Archive aerosol mass and PM2.5 diagnostics
    !-----------------------------------------------------------------------
    IF ( Archive_AerMass ) THEN
       CALL Set_AerMass_Diagnostic( am_I_Root, Input_Opt,  State_Met,        &
                                    State_Chm, State_Diag, RC               )

       ! Trap potential errors
       IF ( RC /= GC_SUCCESS ) THEN
          ErrMsg = 'Error encountered in "Set_AerMass_Diagnostic"!'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF
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

#endif
  END SUBROUTINE Set_SpcConc_Diagnostic
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_aermass_diagnostic
!
! !DESCRIPTION: Computes the aerosol mass diagnostic (formerly ND42 bpch
!  diagnostic).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_AerMass_Diagnostic( am_I_Root, Input_Opt,  State_Met,       &
                                     State_Chm, State_Diag, RC              )
!
! !USES:
!
    USE Aerosol_Mod
    USE Carbon_Mod,     ONLY : BetaNOSave
    USE CMN_Size_Mod
    USE ErrCode_Mod
    USE Input_Opt_Mod,  ONLY : OptInput
    USE Species_Mod,    ONLY : Species
    USE State_Chm_Mod,  ONLY : ChmState
    USE State_Chm_Mod,  ONLY : Ind_
    USE State_Diag_Mod, ONLY : DgnState
    USE State_Met_Mod,  ONLY : MetState
    USE PhysConstants,  ONLY : MwCarb
    USE UnitConv_Mod,   ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS: 
!
    LOGICAL,          INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
    TYPE(OptInput),   INTENT(IN)    :: Input_Opt   ! Input Options object
    TYPE(MetState),   INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object
    TYPE(DgnState),   INTENT(INOUT) :: State_Diag  ! Diagnostic State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(OUT)   :: RC          ! Success or failure?
!
! !REMARKS:
!  NOTE: This diagnostic mimics the bpch diagnostic routine "DIAG42".
!
! !REVISION HISTORY:
!  05 Feb 2018 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
#if defined( NC_DIAG )
    ! SAVEd scalars
    LOGICAL                  :: First = .TRUE.

    ! Scalars
    INTEGER                  :: I, J, L

    ! Strings
    CHARACTER(LEN=63)        :: OrigUnit
    CHARACTER(LEN=255)       :: ThisLoc
    CHARACTER(LEN=512)       :: ErrMsg

    ! Pointers      
    REAL(fp),      POINTER   :: AirVol(:,:,:  )
    REAL(fp),      POINTER   :: Spc   (:,:,:,:)
    TYPE(Species), POINTER   :: SpcInfo
!
! !DEFINED PARAMETERS:
!
    ! Convert [kg/m3] to [ug/m3]
    REAL(fp),      PARAMETER :: kgm3_to_ugm3 = 1.0e+9_fp

    ! Define number of carbon atoms in each irreversible isoprene
    ! SOA tracer species. Named according to the parent HC (same 
    ! number of carbons):
    REAL(fp),      PARAMETER :: NCIMAE   = 4e+0_fp
    REAL(fp),      PARAMETER :: NCIEPOX  = 5e+0_fp
    REAL(fp),      PARAMETER :: NCINDIOL = NCIEPOX
    REAL(fp),      PARAMETER :: NCGLYX   = 2e+0_fp
    REAL(fp),      PARAMETER :: NCGLYC   = NCGLYX
    REAL(fp),      PARAMETER :: NCMGLY   = 3e+0_fp
    REAL(fp),      PARAMETER :: NCLVOC   = NCIEPOX
    REAL(fp),      PARAMETER :: NCISN1   = NCIEPOX

    !=======================================================================
    ! Set_AerMass_Diagnostic begins here!
    !=======================================================================

    ! Initialize
    RC       = GC_SUCCESS
    ErrMsg   = ''
    ThisLoc  = ' -> at Set_AerMass_Diagnostic ' // TRIM( ModLoc )

    ! Define species ID flags for the aerosol mass diagnostics
    IF ( First ) THEN

       !--------------------------------------------------------------------
       ! Look up species indices in State_Chm%SPECIES
       !--------------------------------------------------------------------
       id_INDIOL  = Ind_( 'INDIOL' )
       id_ISN1OA  = Ind_( 'ISN1OA' )
       id_LVOCOA  = Ind_( 'LVOCOA' )
       id_POA1    = Ind_( 'POA1'   )
       id_POA2    = Ind_( 'POA2'   )
       id_SOAGX   = Ind_( 'SOAGX'  )
       id_SOAIE   = Ind_( 'SOAIE'  )
       id_SOAME   = Ind_( 'SOAME'  )
       id_SOAMG   = Ind_( 'SOAMG'  )
       Is_BetaNO  = ALLOCATED( BETANOSAVE )
       Is_POA     = ( id_POA1 > 0 .and. id_POA2 > 0 )

       ! Initialize conversion factors for total OC diagnostic
       Fac_INDIOL = 0.0_fp
       Fac_ISN1OA = 0.0_fp
       Fac_LVOCOA = 0.0_fp
       Fac_SOAGX  = 0.0_fp
       Fac_SOAIE  = 0.0_fp
       Fac_SOAME  = 0.0_fp
       Fac_SOAMG  = 0.0_fp

       !--------------------------------------------------------------------
       ! Set conversion factors for certain isoprene SOA species,
       ! or, if they aren't present, disable their diagnostics
       !--------------------------------------------------------------------
       IF ( id_INDIOL > 0 ) THEN
          SpcInfo    => State_Chm%SpcData(id_INDIOL)%Info
          Fac_INDIOL =  ( NCINDIOL  * MwCarb / ( SpcInfo%Mw_g * 1e-3_fp ) )
          SpcInfo    => NULL()
       ELSE
          IF ( Archive_AerMassINDIOL ) THEN
             Archive_AerMassINDIOL = .FALSE.
             ErrMsg = 'Disabling AerMassINDIOL diagnostic.  INDIOL is '   // &
                      'not a defined species for this simulation.'
             CALL GC_Warning( ErrMsg, RC, ThisLoc )
          ENDIF
       ENDIF

       IF ( id_ISN1OA > 0 ) THEN
          SpcInfo    => State_Chm%SpcData(id_ISN1OA)%Info
          Fac_ISN1OA =  ( NCISN1 * MwCarb / ( SpcInfo%Mw_g * 1e-3_fp ) )
          SpcInfo    => NULL()
       ELSE
          IF ( Archive_AerMassISN1OA ) THEN
             Archive_AerMassISN1OA = .FALSE.
             ErrMsg = 'Disabling AerMassISN1OA diagnostic.  ISN1OA is '   // &
                      'not a defined species for this simulation.'
             CALL GC_Warning( ErrMsg, RC, ThisLoc )
          ENDIF
       ENDIF

       IF ( id_LVOCOA > 0  ) THEN
          SpcInfo    => State_Chm%SpcData(id_LVOCOA)%Info
          Fac_LVOCOA = ( NCLVOC * MwCarb / ( SpcInfo%Mw_G * 1e-3_fp ) )
          SpcInfo    => NULL()
       ELSE
          IF ( Archive_AerMassLVOCOA ) THEN
             Archive_AerMassLVOCOA = .FALSE.
             ErrMsg = 'Disabling AerMassLVOCOA diagnostic.  LVOCOA is '   // &
                      'not a defined species for this simulation.'
             CALL GC_Warning( ErrMsg, RC, ThisLoc )
          ENDIF
       ENDIF

       IF ( id_SOAGX > 0 ) THEN
          SpcInfo    => State_Chm%SpcData(id_SOAGX)%Info
          Fac_SOAGX  = ( NCGLYX * MwCarb / ( SpcInfo%Mw_g * 1e-3_fp ) )
          SpcInfo    => NULL()
       ELSE
          IF ( Archive_AerMassSOAGX ) THEN
             Archive_AerMassSOAGX = .FALSE.
             ErrMsg = 'Disabling AerMassSOAGX diagnostic.  SOAGX is  '    // &
                      'not a defined species for this simulation.'
             CALL GC_Warning( ErrMsg, RC, ThisLoc )
          ENDIF
       ENDIF

       IF ( id_SOAIE > 0 ) THEN
          SpcInfo    => State_Chm%SpcData(id_SOAIE)%Info
          Fac_SOAIE  =  ( NCIEPOX * MwCarb / ( SpcInfo%Mw_g * 1e-3_fp ) )
          SpcInfo    => NULL()
       ELSE
          IF ( Archive_AerMassSOAIE ) THEN
             Archive_AerMassSOAIE = .FALSE.
             ErrMsg = 'Disabling AerMassSOAIE diagnostic.  SOAIE is  '    // &
                      'not a defined species for this simulation.'
             CALL GC_Warning( ErrMsg, RC, ThisLoc )
          ENDIF
       ENDIF

       IF ( id_SOAME > 0 ) THEN
          SpcInfo    => State_Chm%SpcData(id_SOAME)%Info
          Fac_SOAME  =  ( NCIMAE  * MwCarb / ( SpcInfo%Mw_g * 1e-3_fp ) )
          SpcInfo    => NULL()
       ELSE
          IF ( Archive_AerMassSOAME ) THEN
             Archive_AerMassSOAME = .FALSE.
             ErrMsg = 'Disabling AerMassLVOCOA diagnostic.  SOAME is  '   // &
                      'not a defined species for this simulation.'
             CALL GC_Warning( ErrMsg, RC, ThisLoc )
          ENDIF
       ENDIF

       IF ( id_SOAMG > 0 ) THEN
          SpcInfo    => State_Chm%SpcData(id_SOAMG)%Info
          Fac_SOAMG  =  ( NCMGLY  * MwCarb / ( SpcInfo%Mw_g * 1e-3_fp ) )
          SpcInfo    => NULL()
       ELSE
          IF ( Archive_AerMassSOAMG ) THEN
             Archive_AerMassSOAMG = .FALSE.
             ErrMsg = 'Disabling AerMassSOAMG diagnostic.  SOAMG is not ' // &
                      'a defined GEOS-Chem species for this simulation.'
             CALL GC_Warning( ErrMsg, RC, ThisLoc )
          ENDIF
       ENDIF

       ! Reset first-time flag
       First = .FALSE.
    ENDIF

    !=======================================================================
    ! Convert species units to [kg] for this routine
    !=======================================================================
    CALL Convert_Spc_Units( am_I_Root,        Input_Opt, State_Met,          &
                            State_Chm,        'kg',      RC,                 & 
                            OrigUnit=OrigUnit                               )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Convert_Spc_Units"!  Could not '   // &
                'convert units from [' // TRIM( OrigUnit ) // ' to [kg] ' // &
                'before computing aerosol mass and PM2.5 diagnostics!' 
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !=======================================================================
    ! Compute Aerosol mass and PM2.5 diagnostics using concentrations 
    ! from the end of the chemistry timestep, which should be more 
    ! consistent with the legacy ND42 bpch diagnostics
    !=======================================================================

    ! Point to fielss of State_Chm and State_Met
    Spc    => State_Chm%Species
    AirVol => State_Met%AIRVOL

    !$OMP PARALLEL DO         &
    !$OMP DEFAULT( SHARED   ) &
    !$OMP PRIVATE( I, J, L  )
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
       
       !--------------------------------------
       ! AerMassASOA [ug/m3]
       !--------------------------------------
       IF ( Archive_AerMassASOA ) THEN
          State_Diag%AerMassASOA(I,J,L)     = ASOA(I,J,L)                    &
                                            * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassBC [ug C/m3]
       !--------------------------------------
       IF ( Archive_AerMassBC ) THEN
          State_Diag%AerMassBC(I,J,L)       = ( BCPI(I,J,L)                  &
                                              + BCPO(I,J,L) )                &
                                            * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassINDIOL [ug/m3]
       !--------------------------------------
       IF ( Archive_AerMassINDIOL ) THEN
          State_Diag%AerMassINDIOL(I,J,L)   = Spc(I,J,L,id_INDIOL)           &
                                            * kgm3_to_ugm3                   &
                                            / AirVol(I,J,L)
       ENDIF

       !--------------------------------------
       ! AerMassISN1OA [ug/m3]
       !--------------------------------------
       IF ( Archive_AerMassISN1OA ) THEN
          State_Diag%AerMassISN1OA(I,J,L)   = Spc(I,J,L,id_ISN1OA)           &
                                            * kgm3_to_ugm3                   &
                                            / AirVol(I,J,L)
       ENDIF

       !--------------------------------------
       ! AerMassISOA [ug/m3]
       !--------------------------------------
       IF ( Archive_AerMassISOA ) THEN
          State_Diag%AerMassISOA(I,J,L)     = ( ISOA(I,J,L)                  &
                                              + ISOAAQ(I,J,L) )              &
                                            * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassLVOCOA [ug/m3]
       !--------------------------------------
       IF ( Archive_AerMassLVOCOA ) THEN
          State_Diag%AerMassLVOCOA(I,J,L)   = Spc(I,J,L,id_LVOCOA)           &
                                            * kgm3_to_ugm3                   &
                                            / AirVol(I,J,L)
       ENDIF

       !--------------------------------------
       ! AerMassNH4 [ug/m3]
       !--------------------------------------
       IF ( Archive_AerMassNH4 ) THEN
          State_Diag%AerMassNH4(I,J,L)      = NH4(I,J,L)                     &
                                            * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassNIT [ug/m3]
       !--------------------------------------
       IF ( Archive_AerMassNIT ) THEN
          State_Diag%AerMassNIT(I,J,L)      = NIT(I,J,L)                     &
                                            * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassOPOA [ug/m3], OA:OC=2.1
       !--------------------------------------
       IF ( Archive_AerMassOPOA ) THEN
          State_Diag%AerMassOPOA(I,J,L)     = OPOA(I,J,L)                    &
                                            * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassPOA [ug/m3], OA:OC=2.1
       !--------------------------------------
       IF ( Archive_AerMassPOA ) THEN
          IF ( Is_POA ) THEN
             State_Diag%AerMassPOA(I,J,L)   = OCPO(I,J,L)                    &
                                            * kgm3_to_ugm3
          ELSE
             State_Diag%AerMassPOA(I,J,L)   = ( OCPI(I,J,L)                  &
                                              + OCPO(I,J,L) )                &
                                            * kgm3_to_ugm3
          ENDIF
       ENDIF

       !--------------------------------------
       ! AerMassSAL [ug/m3]
       !--------------------------------------
       IF ( Archive_AerMassSAL ) THEN
          State_Diag%AerMassSAL(I,J,L)      = ( SALA(I,J,L)                  &
                                              + SALC(I,J,L) )                &
                                            * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassSO4 [ug/m3]
       !--------------------------------------
       IF ( Archive_AerMassSO4 ) THEN
          State_Diag%AerMassSO4(I,J,L)      = SO4(I,J,L)                     &
                                            * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassSOAGX [ug/m3]
       !--------------------------------------
       IF ( Archive_AerMassSOAGX ) THEN
          State_Diag%AerMassSOAGX(I,J,L)    = SOAGX(I,J,L)                   &
                                            * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassSOAIE [ug/m3]
       !--------------------------------------
       IF ( Archive_AerMassSOAIE ) THEN
          State_Diag%AerMassSOAIE(I,J,L)    = Spc(I,J,L,id_SOAIE)            &
                                            * kgm3_to_ugm3                   &
                                            / AirVol(I,J,L)
       ENDIF

       !--------------------------------------
       ! AerMassSOAME [ug/m3]
       !--------------------------------------
       IF ( Archive_AerMassSOAME ) THEN
          State_Diag%AerMassSOAME(I,J,L)    = Spc(I,J,L,id_SOAME)            &
                                            * kgm3_to_ugm3                   &
                                            / AirVol(I,J,L)
       ENDIF

       !--------------------------------------
       ! AerMassSOAMG [ug/m3]
       !--------------------------------------
       IF ( Archive_AerMassSOAMG ) THEN
          State_Diag%AerMassSOAMG(I,J,L)    = SOAMG(I,J,L)                   &
                                            * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! AerMassTSOA [ug/m3]
       !--------------------------------------
       IF ( Archive_AerMassTSOA ) THEN
          State_Diag%AerMassTSOA(I,J,L)     = TSOA(I,J,L)                    &
                                            * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! BetaNO [ug C/m3]
       !--------------------------------------
       IF ( Archive_BetaNO ) THEN
          State_Diag%BetaNO(I,J,L)          = BetaNOsave(I,J,L)
       ENDIF

       !--------------------------------------
       ! PM25 [ug/m3]
       !--------------------------------------
       IF ( Archive_BetaNO ) THEN
          State_Diag%PM25(I,J,L)            = PM25(I,J,L)                    &
                                            * kgm3_to_ugm3
       ENDIF

       !--------------------------------------
       ! Sum of all biogenic organic aerosol
       !--------------------------------------
       IF ( Archive_TotalBiogenicOA ) THEN
          State_Diag%TotalBiogenicOA(I,J,L) = ( TSOA(I,J,L)                   &
                                              + ISOA(I,J,L)                   &
                                              + ISOAAQ(I,J,L) )               &
                                            * kgm3_to_ugm3  
       ENDIF

       !--------------------------------------
       ! Sum of all organic aerosol [ug/m3]
       !--------------------------------------
       IF ( Archive_TotalOA ) THEN
          State_Diag%TotalOA(I,J,L)         = ( TSOA(I,J,L)                  &
                                              + ISOA(I,J,L)                  & 
                                              + ASOA(I,J,L)                  &
                                              + OCPO(I,J,L)                  &
                                              + OCPI(I,J,L)                  &
                                              + OPOA(I,J,L)                  &
                                              + ISOAAQ(I,J,L) )              & 
                                            * kgm3_to_ugm3

       ENDIF

       !--------------------------------------
       ! Sum of all organic carbon [ug/m3]
       !--------------------------------------
       IF ( Archive_TotalOC ) THEN
          
          IF ( Is_POA ) THEN
             State_Diag%TotalOC(I,J,L) =                              &
                  ( ( TSOA(I,J,L) + ISOA(I,J,L) + ASOA(I,J,L)         &
                    + OCPI(I,J,L) + OPOA(I,J,L) ) / OCFOPOA(I,J)    &
                    + OCPO(I,J,L) / OCFPOA(I,J) ) * kgm3_to_ugm3 
            
          ELSE
             State_Diag%TotalOC(I,J,L) =                              &
                  ( ( TSOA(I,J,L) + ISOA(I,J,L) + ASOA(I,J,L)         &
                    + OCPO(I,J,L) + OCPI(I,J,L) + OPOA(I,J,L) )       &
                    / OCFOPOA(I,J) ) * kgm3_to_ugm3
          ENDIF
         
          IF ( Input_Opt%LSOA ) THEN
             State_Diag%TotalOC(I,J,L) =  State_Diag%TotalOC(I,J,L)   &
                        + ( ( Spc(I,J,L,id_SOAIE)  * Fac_SOAIE  )     &
                          + ( Spc(I,J,L,id_SOAME)  * Fac_SOAME  )     &
                          + ( Spc(I,J,L,id_INDIOL) * Fac_INDIOL )     &
                          + ( Spc(I,J,L,id_SOAGX)  * Fac_SOAGX  )     &
                          + ( Spc(I,J,L,id_SOAMG)  * Fac_SOAMG  )     &
                          + ( Spc(I,J,L,id_LVOCOA) * Fac_LVOCOA )     &
                          + ( Spc(I,J,L,id_ISN1OA) * Fac_ISN1OA ) )   &
                        / AirVol(I,J,L) * kgm3_to_ugm3
         ENDIF

      ENDIF
!------------------------------------------------------------------------------

    ENDDO
    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Free pointers
    Spc      => NULL()
    AirVol   => NULL()

    !=======================================================================
    ! Convert species units to [kg] for this routine
    !=======================================================================
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met,                 &
                            State_Chm, OrigUnit,  RC                        )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Error encountered in "Convert_Spc_Units"!  Could not '   // &
                'convert units from [kg] back to the original unit of ['  // &
                TRIM( OrigUnit ) // '] after computing aerosol mass '     // &
                'and PM2.5 diagnostics'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

#endif
  END SUBROUTINE Set_AerMass_Diagnostic
!EOC
END MODULE Diagnostics_mod
