!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: cleanup.F90
!
! !DESCRIPTION: Subroutine CLEANUP deallocates the memory assigned to
!  dynamically allocatable arrays just before exiting a GEOS-Chem simulation.
!\\
!\\
! !INTERFACE:
!
SUBROUTINE CLEANUP( Input_Opt, State_Grid, ERROR, RC )
!
! !USES:
!
  USE AEROSOL_MOD,             ONLY : CLEANUP_AEROSOL
  USE CARBON_MOD,              ONLY : CLEANUP_CARBON
  USE CO2_MOD,                 ONLY : CLEANUP_CO2
  USE CMN_FJX_Mod,             ONLY : Cleanup_CMN_FJX
  USE CMN_SIZE_Mod,            ONLY : Cleanup_CMN_SIZE
  USE DEPO_MERCURY_MOD,        ONLY : CLEANUP_DEPO_MERCURY
  USE DIAG_OH_MOD,             ONLY : CLEANUP_DIAG_OH
  USE DRYDEP_MOD,              ONLY : CLEANUP_DRYDEP
  USE DUST_MOD,                ONLY : CLEANUP_DUST
  USE ErrCode_Mod
  USE ERROR_MOD,               ONLY : DEBUG_MSG
  USE FLEXCHEM_MOD,            ONLY : CLEANUP_FLEXCHEM
  USE GLOBAL_Br_MOD,           ONLY : CLEANUP_GLOBAL_Br
  USE GLOBAL_CH4_MOD,          ONLY : CLEANUP_GLOBAL_CH4
  USE Grid_Registry_Mod,       ONLY : Cleanup_Grid_Registry
  USE History_Mod,             ONLY : History_Cleanup
  USE Input_Opt_Mod,           ONLY : OptInput
  USE ISORROPIAII_MOD,         ONLY : CLEANUP_ISORROPIAII
  USE LAND_MERCURY_MOD,        ONLY : CLEANUP_LAND_MERCURY
  USE MERCURY_MOD,             ONLY : CLEANUP_MERCURY
  USE ObsPack_Mod,             ONLY : ObsPack_SpeciesMap_Cleanup
  USE OCEAN_MERCURY_MOD,       ONLY : CLEANUP_OCEAN_MERCURY
  USE PBL_MIX_MOD,             ONLY : CLEANUP_PBL_MIX
  USE PJC_PFIX_MOD,            ONLY : CLEANUP_PJC_PFIX
  USE PLANEFLIGHT_MOD,         ONLY : CLEANUP_PLANEFLIGHT
  USE POPs_Mod,                ONLY : Cleanup_POPs
  USE PRESSURE_MOD,            ONLY : CLEANUP_PRESSURE
  USE Regrid_A2A_Mod,          ONLY : Cleanup_Map_A2a
  USE SEASALT_MOD,             ONLY : CLEANUP_SEASALT
  USE SULFATE_MOD,             ONLY : CLEANUP_SULFATE
  USE State_Grid_Mod,          ONLY : GrdState
  USE STRAT_CHEM_MOD,          ONLY : CLEANUP_STRAT_CHEM
  USE TAGGED_CO_MOD,           ONLY : CLEANUP_TAGGED_CO
  USE UCX_MOD,                 ONLY : CLEANUP_UCX
  USE EMISSIONS_MOD,           ONLY : EMISSIONS_FINAL
  USE SFCVMR_MOD,              ONLY : FixSfcVmr_Final
#ifdef BPCH_DIAG
  USE CMN_O3_Mod,              ONLY : Cleanup_CMN_O3
  USE DIAG_MOD,                ONLY : CLEANUP_DIAG
  USE DIAG03_MOD,              ONLY : CLEANUP_DIAG03
  USE DIAG51_MOD,              ONLY : CLEANUP_DIAG51
  USE DIAG53_MOD,              ONLY : CLEANUP_DIAG53
#endif
#ifdef TOMAS
  USE TOMAS_MOD,               ONLY : CLEANUP_TOMAS  !sfarina, 1/16/13
#endif
  USE TPCORE_FVDAS_MOD,        ONLY : EXIT_TPCORE
  USE TPCORE_WINDOW_MOD,       ONLY : EXIT_TPCORE_WINDOW
#if !defined( ESMF_ ) && !defined( MODEL_ )
  USE TRANSPORT_MOD,           ONLY : CLEANUP_TRANSPORT
#endif
#ifdef RRTMG
  USE RRTMG_RAD_TRANSFER_MOD,  ONLY : Cleanup_RRTMG_Rad_Transfer
#endif

  IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
  TYPE(OptInput), INTENT(IN)  :: Input_Opt    ! Input options
  TYPE(GrdState), INTENT(IN)  :: State_Grid   ! Grid state object
  LOGICAL,        INTENT(IN)  :: ERROR        ! Cleanup after error?

!
! !OUTPUT PARAMETERS:
!
  INTEGER,        INTENT(OUT) :: RC           ! Success or failure
!
! !REVISION HISTORY:
!  29 Nov 1999 - R. Yantosca - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
  ! Strings
  CHARACTER(LEN=255) :: ErrMsg, ThisLoc

  !=================================================================
  ! CLEANUP begins here!
  !=================================================================

  ! Initialize
  RC      = GC_SUCCESS
  ErrMsg  = ''
  ThisLoc = ' -> at CLEANUP (in module GeosCore/cleanup.F)'

  ! Echo info
  IF ( Input_Opt%amIRoot ) THEN
     WRITE( 6, 100 )
  ENDIF
100 FORMAT( '     - CLEANUP: deallocating arrays now...' )

  !=================================================================
  !         ***** H I S T O R Y   C L E A N U P *****
  !
  ! Finalize the History Component.
  ! Also closes all netCDF files that may still be open.
  !=================================================================

  ! Finalize the history component
  CALL History_Cleanup( RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "History_Cleanup"!'
     CALL GC_Error( ErrMsg, RC, ThisLoc )
     RETURN
  ENDIF

  ! UCX needs to be cleaned up before emissions, because the UCX
  ! restart variables needs to be passed to the HEMCO diagnostics
  ! first.
  CALL CLEANUP_UCX()

  !=================================================================
  ! Cleanup HEMCO
  !=================================================================
  CALL EMISSIONS_FINAL( ERROR, RC )

  ! Cleanup surface mixing concentration module
  CALL FixSfcVmr_Final( RC )

  !=================================================================
  ! Call cleanup routines from individual F90 modules
  !=================================================================
  CALL CLEANUP_AEROSOL()
  CALL CLEANUP_CARBON()
  CALL CLEANUP_CO2()
  CALL CLEANUP_DIAG_OH()
  CALL CLEANUP_DRYDEP()
  CALL CLEANUP_DUST()
  CALL CLEANUP_ISORROPIAII()
  CALL CLEANUP_MAP_A2A()
  CALL CLEANUP_PBL_MIX()
  CALL CLEANUP_PJC_PFIX()
  CALL CLEANUP_PRESSURE()
  CALL CLEANUP_SEASALT()
  CALL CLEANUP_SULFATE()
  CALL CLEANUP_STRAT_CHEM()

  CALL Cleanup_FlexChem( RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_FlexChem"!'
     CALL GC_Error( ErrMsg, RC, ThisLoc )
     RETURN
  ENDIF

  CALL Cleanup_Global_CH4( RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_Global_CH4"!'
     CALL GC_Error( ErrMsg, RC, ThisLoc )
     RETURN
  ENDIF

  CALL Cleanup_Grid_Registry( RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_Grid_Registry"!'
     CALL GC_Error( ErrMsg, RC, ThisLoc )
     RETURN
  ENDIF

  IF ( State_Grid%NestedGrid ) THEN
     CALL EXIT_TPCORE_WINDOW()
  ELSE
     CALL EXIT_TPCORE()
  ENDIF

  ! Cleanup Tagged CO code
  CALL CLEANUP_TAGGED_CO( RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_Tagged_CO"!'
     CALL GC_Error( ErrMsg, RC, ThisLoc )
     RETURN
  ENDIF

  CALL Cleanup_CMN_SIZE( RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_CMN_SIZE"!'
     CALL GC_Error( ErrMsg, RC, ThisLoc )
     RETURN
  ENDIF

  CALL Cleanup_CMN_FJX( RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_CMN_FJX"!'
     CALL GC_Error( ErrMsg, RC, ThisLoc )
     RETURN
  ENDIF

  CALL CLEANUP_MERCURY()
  CALL CLEANUP_OCEAN_MERCURY()
  CALL CLEANUP_DEPO_MERCURY()
  CALL CLEANUP_LAND_MERCURY()
  CALL CLEANUP_PLANEFLIGHT()

  CALL Cleanup_Global_Br( RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_Global_Br"!'
     CALL GC_Error( ErrMsg, RC, ThisLoc )
     RETURN
  ENDIF

  CALL Cleanup_POPs( RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_POPs"!'
     CALL GC_Error( ErrMsg, RC, ThisLoc )
     RETURN
  ENDIF

#ifdef BPCH_DIAG
  !=====================================================================
  ! These routines are only needed when GEOS-Chem
  ! is compiled with BPCH_DIAG=y
  !=====================================================================
  CALL CLEANUP_DIAG()
  CALL CLEANUP_DIAG03()
  CALL CLEANUP_DIAG51()
  CALL CLEANUP_DIAG53()

  CALL Cleanup_CMN_O3( RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_CMN_O3"!'
     CALL GC_Error( ErrMsg, RC, ThisLoc )
     RETURN
  ENDIF
#endif

#ifdef RRTMG
  !=====================================================================
  ! These routines are only needed when GEOS-Chem
  ! is compiled with RRTMG=y
  !=====================================================================
  CALL Cleanup_RRTMG_Rad_Transfer( RC )
  IF ( RC /= GC_SUCCESS ) THEN
     ErrMsg = 'Error encountered in "Cleanup_RRTMG_Rad_Transfer"!'
     CALL GC_Error( ErrMsg, RC, ThisLoc )
     RETURN
  ENDIF
#endif

#ifdef TOMAS
  !=====================================================================
  ! These routines are only needed when GEOS-Chem
  ! is compiled with TOMAS{12|15|30|40}=y
  !=====================================================================
  CALL CLEANUP_TOMAS()
#endif

#if !defined( ESMF_ ) && !defined( MODEL_ )
  !=====================================================================
  ! These routines are only needed when GEOS-Chem
  ! is compiled for GCHP, or for external ESMs.
  !=====================================================================
  CALL CLEANUP_TRANSPORT()
#endif

END SUBROUTINE CLEANUP
!EOC
