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
  PUBLIC :: Cleanup_Set_Global_CH4
!
! !REVISION HISTORY:
!  18 Jan 2018 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE TYPES:
!
  REAL(f4), POINTER :: SFC_CH4(:,:) => NULL() ! Global surface CH4

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
    USE HCO_EMISLIST_MOD,  ONLY : HCO_GetPtr
    USE HCO_Error_Mod
    USE HCO_INTERFACE_MOD, ONLY : HcoState
    USE Input_Opt_Mod,     ONLY : OptInput
    USE State_Chm_Mod,     ONLY : ChmState, Ind_
    USE State_Diag_Mod,    ONLY : DgnState
    USE State_Grid_Mod,    ONLY : GrdState
    USE State_Met_Mod,     ONLY : MetState
    USE UnitConv_Mod,      ONLY : Convert_Spc_Units
#if defined( MODEL_GEOS )
    USE PhysConstants,     ONLY : AIRMW
    USE TIME_MOD,          ONLY : GET_TS_DYN
#endif
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
    INTEGER             :: I, J, L, PBL_TOP
    CHARACTER(LEN=63)   :: OrigUnit
    REAL(fp)            :: CH4
    LOGICAL             :: FOUND
#if defined( MODEL_GEOS )
    INTEGER             :: DT
    REAL(fp)            :: dCH4, MWCH4
    LOGICAL             :: PseudoFlux
#endif

    ! Strings
    CHARACTER(LEN=255)  :: ErrMsg
    CHARACTER(LEN=255)  :: ThisLoc

    ! SAVEd variables
    LOGICAL, SAVE       :: FIRST = .TRUE.
    INTEGER, SAVE       :: id_CH4

    !=================================================================
    ! SET_GLOBAL_CH4 begins here!
    !=================================================================

    ! Assume success
    RC      = GC_SUCCESS
    ErrMsg  = ''
    ThisLoc = ' -> at SET_CH4 (in module GeosCore/set_global_ch4_mod.F90)'

    ! Skip unless we are doing a fullchem simulation
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
       RETURN
    ENDIF

    IF ( FIRST ) THEN

       ! Get species ID
       id_CH4 = Ind_( 'CH4' )

       ! Use the NOAA spatially resolved data where available
       CALL HCO_GetPtr( HcoState, 'NOAA_GMD_CH4', SFC_CH4, RC, FOUND=FOUND )
       IF (.NOT. FOUND ) THEN
          FOUND = .TRUE.
          ! Use the CMIP6 data from Meinshausen et al. 2017, GMD
          ! https://doi.org/10.5194/gmd-10-2057-2017a
          CALL HCO_GetPtr( HcoState, 'CMIP6_Sfc_CH4', SFC_CH4, RC, FOUND=FOUND )
       ENDIF
       IF (.NOT. FOUND ) THEN
          ErrMsg = 'Cannot get pointer to NOAA_GMD_CH4 or CMIP6_Sfc_CH4 ' // &
                   'in SET_CH4! Make sure the data source corresponds '   // &
                   'to your emissions year in HEMCO_Config.rc (NOAA GMD ' // &
                   'for 1978 and later; else CMIP6).'
          CALL GC_Error( ErrMsg, RC, ThisLoc )
          RETURN
       ENDIF

       ! Reset first-time flag
       FIRST = .FALSE.

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

#if defined( MODEL_GEOS )
    ! Write out pseudo (implied) CH4 flux?
    PseudoFlux = ASSOCIATED(State_Diag%CH4pseudoFlux)
    MWCH4      = State_Chm%SpcData(id_CH4)%Info%emMW_g
    IF ( MWCH4 <= 0.0_fp ) MWCH4 = 16.0_fp
    DT         = GET_TS_DYN()
#endif

    !$OMP PARALLEL DO                 &
    !$OMP DEFAULT( SHARED )           &
    !$OMP PRIVATE( I, J, L, PBL_TOP, CH4 ) &
#if defined( MODEL_GEOS )
    !$OMP PRIVATE( dCH4 ) &
#endif
    !$OMP SCHEDULE( DYNAMIC )
    DO J = 1, State_Grid%NY
    DO I = 1, State_Grid%NX

       ! Top level of boundary layer at (I,J)
       PBL_TOP = CEILING( State_Met%PBL_TOP_L(I,J) )

       ! Surface CH4 from HEMCO is in units [ppbv], convert to [v/v dry]
       CH4 = SFC_CH4(I,J) * 1e-9_fp

#if defined( MODEL_GEOS )
       ! Zero diagnostics
       IF ( PseudoFlux ) State_Diag%CH4pseudoFlux(I,J) = 0.0_fp
#endif

       ! Prescribe methane concentrations throughout PBL
       DO L=1,PBL_TOP

#if defined( MODEL_GEOS )
          ! Eventually compute implied CH4 flux
          IF ( PseudoFlux ) THEN
             ! v/v dry
             dCH4 = CH4 - State_Chm%Species(I,J,L,id_CH4)
             ! Convert to kg/kg dry
             dCH4 = dCH4 * MWCH4 / AIRMW
!             ! Convert to kg/m2/s
!             dCH4 = dCH4 * State_Met%AIRDEN(I,J,L) &
!                  * State_Met%BXHEIGHT(I,J,L) / DT
              dCH4 = dCH4 * State_Met%AD(I,J,L) / State_Met%AREA_M2(I,J) / DT
             ! Accumulate statistics
             State_Diag%CH4pseudoFlux(I,J) = &
                State_Diag%CH4pseudoFlux(I,J) + dCH4
          ENDIF
#endif

          State_Chm%Species(I,J,L,id_CH4) = CH4
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
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_set_global_ch4
!
! !DESCRIPTION: Subroutine CLEANUP\_SET\_GLOBAL\_CH4 deallocates memory from
!  previously allocated module arrays.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Cleanup_Set_Global_CH4
!
! !REVISION HISTORY:
!  18 Jan 2018 - M. Sulprizio- Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Free pointers
    SFC_CH4     => NULL()

  END SUBROUTINE Cleanup_Set_Global_CH4
!EOC
END MODULE Set_Global_CH4_Mod
