!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: set_global_ch4_mod.F90
!
! !DESCRIPTION: Module SET\_GLOBAL\_CH4 contains variables and routines for
!  reading and applying surface CH4 concentrations from NOAA GMD data
!  (provided by Lee Murray)
!\\
!\\
! !INTERFACE:
!
MODULE Set_Global_CH4_Mod
!
! !USES:
!
  USE PRECISION_MOD

  IMPLICIT NONE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: Set_CH4
!
! !REVISION HISTORY:
!  18 Jan 2018 - M. Sulprizio- Initial version
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
! !IROUTINE: set_ch4
!
! !DESCRIPTION: Subroutine SET\_CH4 copies monthly mean surface CH4 from
!  HEMCO and applies it to CH4 concentrations in State\_Chm%Species.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_CH4( Input_Opt, State_Chm, State_Diag, State_Grid, &
                      State_Met, RC )
!
! !USES:
!
    USE ErrCode_Mod
    USE ERROR_MOD
    USE HCO_Error_Mod
    USE HCO_Utilities_GC_Mod, ONLY : HCO_GC_EvalFld
    USE Input_Opt_Mod,     ONLY : OptInput
    USE PhysConstants,     ONLY : AIRMW
    USE State_Chm_Mod,     ONLY : ChmState, Ind_
    USE State_Diag_Mod,    ONLY : DgnState
    USE State_Grid_Mod,    ONLY : GrdState
    USE State_Met_Mod,     ONLY : MetState
    USE TIME_MOD,          ONLY : GET_TS_DYN
    USE UnitConv_Mod
!
! !INPUT PARAMETERS:
!
    TYPE(OptInput), INTENT(IN)    :: Input_Opt ! Input Options object
    TYPE(GrdState), INTENT(IN)    :: State_Grid! Grid State object
    TYPE(MetState), INTENT(IN)    :: State_Met ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm  ! Chemistry State object
    TYPE(DgnState), INTENT(INOUT) :: State_Diag ! Diagnostics State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC        ! Success or failure?
!
! !REMARKS:
!  Lee Murray wrote:
!   The monthly mean surface methane mixing ratios used here were spatially
!   interpolated from NOAA GLOBALVIEW discrete monthly mean flask data using
!   ordinary kriging for 1983-2016 (ch4_flask_surface_2017-07-28.tar.gz).
!   Surface mixing ratios are extended back to 1979 and forward to 2020 via
!   linear extrapolation of the local 1983-1990 and 2011-2016 trends,
!   respectively.
!
! !REVISION HISTORY:
!  18 Jan 2018 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: I, J, L, PBL_TOP, id_CH4, DT
    INTEGER             :: previous_units
    REAL(fp)            :: CH4, dCH4
    LOGICAL             :: FOUND

    ! Arrays
    INTEGER, TARGET     :: mapping(1)

    ! Pointers
    INTEGER, POINTER    :: theMapping(:)

    ! Strings
    CHARACTER(LEN=255)  :: ErrMsg
    CHARACTER(LEN=255)  :: ThisLoc

#if defined( MODEL_GEOS )
    REAL(hp), ALLOCATABLE :: GEOS_CH4(:,:,:)
    REAL(hp), ALLOCATABLE :: CH4_OFFSET(:,:)
    LOGICAL               :: USE_GEOS_CH4
#endif
    LOGICAL, SAVE         :: FIRST = .TRUE.
    CHARACTER(LEN=255)    :: SRCNAME

    !=================================================================
    ! SET_CH4 begins here!
    !=================================================================

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at SET_CH4 (in module GeosCore/set_global_ch4_mod.F90)'

    ! Skip unless we are doing a fullchem simulation
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
       RETURN
    ENDIF

    ! Get species ID
    id_CH4     = Ind_( 'CH4' )
    mapping(1) = id_CH4
    theMapping => mapping

    ! Get dynamic timestep
    DT = GET_TS_DYN()

    FOUND   = .FALSE.
    SRCNAME = ''
#if defined( MODEL_GEOS )
    ! Check for CH4 offset first
    ALLOCATE(CH4_OFFSET(State_Grid%NX,State_Grid%NY))
    CH4_OFFSET(:,:) = 0.0
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'CH4_OFFSET', &
                         CH4_OFFSET, RC, FOUND=FOUND )
    IF ( .NOT. FOUND ) CH4_OFFSET = 0.0
    ! Now get CH4 concentrations
    ALLOCATE(GEOS_CH4(State_Grid%NX,State_Grid%NY,State_Grid%NZ))
    GEOS_CH4(:,:,:) = 0.0
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'GEOS_CH4', &
                         GEOS_CH4, RC, FOUND=FOUND )
    USE_GEOS_CH4 = FOUND
    IF ( FOUND ) SRCNAME = 'GEOS_CH4'
#endif

    ! Use the NOAA spatially resolved data where available
    IF (.NOT. FOUND ) THEN
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'NOAA_GMD_CH4', &
                            State_Chm%SFC_CH4, RC, FOUND=FOUND )
       IF ( FOUND ) SRCNAME = 'NOAA_GMD_CH4'
    ENDIF
    IF (.NOT. FOUND ) THEN
       FOUND = .TRUE.
       ! Use the CMIP6 data from Meinshausen et al. 2017, GMD
       ! https://doi.org/10.5194/gmd-10-2057-2017a
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'CMIP6_Sfc_CH4', &
                            State_Chm%SFC_CH4, RC, FOUND=FOUND )
       IF ( FOUND ) SRCNAME = 'CMIPS_Sfc_CH4'
    ENDIF
    IF (.NOT. FOUND ) THEN
       FOUND = .TRUE.
       ! Use the CMIP6 data boundary conditions processed for GCAP 2.0
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'SfcVMR_CH4', &
                            State_Chm%SFC_CH4, RC, FOUND=FOUND )
       IF ( FOUND ) SRCNAME = 'SfcVMR_CH4'
    ENDIF
    IF (.NOT. FOUND ) THEN
       ErrMsg = 'Cannot retrieve data for NOAA_GMD_CH4, CMIP6_Sfc_CH4, or ' // &
                'SfcVMR_CH4 from HEMCO! Make sure the data source ' // &
                'corresponds to your emissions year in HEMCO_Config.rc ' // &
                '(NOAA GMD for 1978 and later; else CMIP6). To use the last year ' // &
                'available you can change the time cycle flag in HEMCO_Config.rc for ' // &
                'the inventory from RY to CY.'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Convert species to [v/v dry] aka [mol/mol dry] for this routine
    CALL Convert_Spc_Units(                                                  &
         Input_Opt      = Input_Opt,                                         &
         State_Chm      = State_Chm,                                         &
         State_Grid     = State_Grid,                                        &
         State_Met      = State_Met,                                         &
         mapping        = theMapping,                                        &
         new_units      = MOLES_SPECIES_PER_MOLES_DRY_AIR,                   &
         previous_units = previous_units,                                    &
         RC             = RC                                                )

    ! Add info to logfile
    IF ( FOUND .AND. Input_Opt%amIRoot .AND. FIRST ) THEN
       WRITE(*,*) 'Getting CH4 boundary conditions in GEOS-Chem from :'//TRIM(SRCNAME)
       FIRST = .FALSE.
    ENDIF 

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error at start of "SET_CH4"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    !$OMP PARALLEL DO                            &
    !$OMP DEFAULT( SHARED )                      &
    !$OMP PRIVATE( I, J, L, PBL_TOP, CH4, dCH4 ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Top level of boundary layer at (I,J)
       PBL_TOP = CEILING( State_Met%PBL_TOP_L(I,J) )

       ! Surface CH4 from HEMCO is in units [ppbv], convert to [v/v dry]
       CH4 = State_Chm%SFC_CH4(I,J) * 1e-9_fp

       ! Zero diagnostics
       IF ( State_Diag%Archive_CH4pseudoFlux ) THEN
          State_Diag%CH4pseudoFlux(I,J) = 0.0_fp
       ENDIF

       ! Prescribe methane concentrations throughout PBL
       DO L=1,PBL_TOP

       ! In GEOS, we may be getting CH4 from a 3D field 
#if defined( MODEL_GEOS )
          IF ( USE_GEOS_CH4 ) CH4 = GEOS_CH4(I,J,L) + CH4_OFFSET(I,J)
#endif

          ! Compute implied CH4 flux if diagnostic is on
          IF ( State_Diag%Archive_CH4pseudoFlux ) THEN
             ! v/v dry
             dCH4 = CH4 - State_Chm%Species(id_CH4)%Conc(I,J,L)
             ! Convert to kg/kg dry
             dCH4 = dCH4 * State_Chm%SpcData(id_CH4)%Info%MW_g / AIRMW
             ! Convert to kg/m2/s
             dCH4 = dCH4 * State_Met%AD(I,J,L) / State_Met%AREA_M2(I,J) / DT
             ! Accumulate statistics
             State_Diag%CH4pseudoFlux(I,J) = &
                State_Diag%CH4pseudoFlux(I,J) + dCH4
          ENDIF

          State_Chm%Species(id_CH4)%Conc(I,J,L) = CH4
       ENDDO

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Convert species back to original unit
    CALL Convert_Spc_Units(                                                  &
         Input_Opt  = Input_Opt,                                             &
         State_Chm  = State_Chm,                                             &
         State_Grid = State_Grid,                                            &
         State_Met  = State_Met,                                             &
         mapping    = theMapping,                                            &
         new_units  = previous_units,                                        &
         RC         = RC                                                    )

    ! Free pointer
    theMapping => NULL()

#if defined( MODEL_GEOS )
    ! Cleanup
    IF(ALLOCATED(GEOS_CH4))   DEALLOCATE(GEOS_CH4)
    IF(ALLOCATED(CH4_OFFSET)) DEALLOCATE(CH4_OFFSET)
#endif

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error at end of "SET_CH4"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Set_CH4
!EOC
END MODULE Set_Global_CH4_Mod
