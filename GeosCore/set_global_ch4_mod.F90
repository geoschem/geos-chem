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
    USE UnitConv_Mod,      ONLY : Convert_Spc_Units
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
    CHARACTER(LEN=63)   :: OrigUnit
    REAL(fp)            :: CH4, dCH4
    LOGICAL             :: FOUND

    ! Strings
    CHARACTER(LEN=255)  :: ErrMsg
    CHARACTER(LEN=255)  :: ThisLoc

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
    id_CH4 = Ind_( 'CH4' )

    ! Get dynamic timestep
    DT = GET_TS_DYN()

    ! Use the NOAA spatially resolved data where available
    CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'NOAA_GMD_CH4', &
                         State_Chm%SFC_CH4, RC, FOUND=FOUND )
    IF (.NOT. FOUND ) THEN
       FOUND = .TRUE.
       ! Use the CMIP6 data from Meinshausen et al. 2017, GMD
       ! https://doi.org/10.5194/gmd-10-2057-2017a
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'CMIP6_Sfc_CH4', &
                            State_Chm%SFC_CH4, RC, FOUND=FOUND )
    ENDIF
    IF (.NOT. FOUND ) THEN
       FOUND = .TRUE.
       ! Use the CMIP6 data boundary conditions processed for GCAP 2.0
       CALL HCO_GC_EvalFld( Input_Opt, State_Grid, 'SfcVMR_CH4', &
                            State_Chm%SFC_CH4, RC, FOUND=FOUND )
    ENDIF
    IF (.NOT. FOUND ) THEN
       ErrMsg = 'Cannot retrieve data for NOAA_GMD_CH4, CMIP6_Sfc_CH4, or ' // &
                'SfcVMR_CH4 from HEMCO! Make sure the data source ' // &
                'corresponds to your emissions year in HEMCO_Config.rc ' // &
                '(NOAA GMD for 1978 and later; else CMIP6).'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

    ! Convert species to [v/v dry] for this routine
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            'v/v dry', RC, OrigUnit=OrigUnit )

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
    CALL Convert_Spc_Units( Input_Opt, State_Chm, State_Grid, State_Met, &
                            OrigUnit, RC )

    ! Trap potential errors
    IF ( RC /= GC_SUCCESS ) THEN
       ErrMsg = 'Unit conversion error at end of "SET_CH4"!'
       CALL GC_Error( ErrMsg, RC, ThisLoc )
       RETURN
    ENDIF

  END SUBROUTINE Set_CH4
!EOC
END MODULE Set_Global_CH4_Mod
