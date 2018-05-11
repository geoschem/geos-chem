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
  SUBROUTINE Set_CH4( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
!
! !USES:
!
    USE CMN_SIZE_MOD
    USE ErrCode_Mod
    USE ERROR_MOD
    USE HCO_EMISLIST_MOD,  ONLY : HCO_GetPtr 
    USE HCO_Error_Mod
    USE HCO_INTERFACE_MOD, ONLY : HcoState
    USE Input_Opt_Mod,     ONLY : OptInput
    USE PBL_MIX_MOD,       ONLY : GET_PBL_TOP_L
    USE State_Chm_Mod,     ONLY : ChmState, Ind_
    USE State_Met_Mod,     ONLY : MetState
    USE UnitConv_Mod,      ONLY : Convert_Spc_Units
!
! !INPUT PARAMETERS:
!
    LOGICAL,        INTENT(IN)    :: am_I_Root ! Are we on the root CPU?
    TYPE(OptInput), INTENT(IN)    :: Input_Opt ! Input Options object
    TYPE(MetState), INTENT(IN)    :: State_Met ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
    TYPE(ChmState), INTENT(INOUT) :: State_Chm ! Chemistry State object
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    INTEGER             :: I, J, L, PBL_TOP
    CHARACTER(LEN=63)   :: OrigUnit
    CHARACTER(LEN=255)  :: LOC='SET_CH4 (set_global_ch4_mod.F90)'
    REAL(fp)            :: CH4

    ! SAVEd variables
    LOGICAL, SAVE       :: FIRST = .TRUE.
    INTEGER, SAVE       :: id_CH4

    !=================================================================
    ! SET_GLOBAL_CH4 begins here!
    !=================================================================

    ! Assume success
    RC = GC_SUCCESS

    ! Skip unless we are doing a fullchem simulation
    IF ( .not. Input_Opt%ITS_A_FULLCHEM_SIM ) THEN
       RETURN
    ENDIF

    IF ( FIRST ) THEN

       ! Get species ID
       id_CH4 = Ind_( 'CH4' )

       ! Get pointer to surface CH4 data
       CALL HCO_GetPtr( am_I_Root, HcoState, 'NOAA_GMD_CH4', SFC_CH4, RC )
       IF ( RC /= GC_SUCCESS ) THEN
          CALL ERROR_STOP ( 'Cannot get pointer to NOAA_GMD_CH4', LOC)
       ENDIF

       ! Reset first-time flag
       FIRST = .FALSE.

    ENDIF

    ! Convert species to [v/v dry] for this routine
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, 'v/v dry', RC, &
                            OrigUnit=OrigUnit )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, 'Start of SET_GLOBAL_CH4')
       RETURN
    ENDIF

    !$OMP PARALLEL DO                 &
    !$OMP DEFAULT( SHARED )           &
    !$OMP PRIVATE( I, J, L, PBL_TOP, CH4 ) &
    !$OMP SCHEDULE( DYNAMIC )
    DO J=1,JJPAR
    DO I=1,IIPAR

       ! Top level of boundary layer at (I,J)
       PBL_TOP = CEILING( GET_PBL_TOP_L(I,J) )

       ! Surface CH4 from HEMCO is in units [ppbv], convert to [v/v dry]
       CH4 = SFC_CH4(I,J) * 1e-9_fp

       ! Prescribe methane concentrations throughout PBL
       DO L=1,PBL_TOP
          State_Chm%Species(I,J,L,id_CH4) = CH4
       ENDDO

    ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    ! Convert species back to original unit
    CALL Convert_Spc_Units( am_I_Root, Input_Opt, State_Met, &
                            State_Chm, OrigUnit, RC )
    IF ( RC /= GC_SUCCESS ) THEN
       CALL GC_Error('Unit conversion error', RC, 'End of SET_GLOBAL_CH4')
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
!EOP
!------------------------------------------------------------------------------
!BOC

    ! Free pointers
    SFC_CH4     => NULL()

  END SUBROUTINE Cleanup_Set_Global_CH4
!EOC
END MODULE Set_Global_CH4_Mod
