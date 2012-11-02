#if defined (ESMF_)
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_finalization_mod
!
! !DESCRIPTION: Module GIGC\_FINALIZATION\_MOD is the module that contains 
!  the GEOS-Chem finalize methods for the ESMF framework.
!\\
!\\
! !INTERFACE: 
!      
MODULE GIGC_Finalization_Mod
!
! !USES:
!      
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!      
  PUBLIC :: GIGC_Finalize
  PUBLIC :: GIGC_Cleanup_GeosChem
!
! !REVISION HISTORY:
!  19 Oct 2012 - R. Yantosca - Added ProTeX headers
!  22 Oct 2012 - R. Yantosca - Renamed to gigc_finalization_mod.F90
!  02 Nov 2012 - R. Yantosca - Now cleanup the Input Options object
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_finalize
!
! !DESCRIPTION: GIGC\_Finalize calls the cleanup routines which deallocate
!  memory from the State objects of the Grid-Independent GEOS-Chem (aka
!  "GIGC") code.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_Finalize( am_I_Root, Input_Opt, State_Chm, State_Met, RC )
!
! !USES:
!
    USE GIGC_ErrCode_Mod
    USE GIGC_Input_Opt_Mod
    USE GIGC_State_Chm_Mod
    USE GIGC_State_Met_Mod
!
! !INPUT PARAMETERS: 
!
    LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
!
! !INPUT/OUTPUT PARAMETERS: 
!
    TYPE(OptInput), INTENT(INOUT) :: Input_Opt   ! Input Options object
    TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
    TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !OUTPUT PARAMETERS:
!
    INTEGER,        INTENT(OUT)   :: RC          ! Success or failure
!
! !REMARKS:
! 
! !REVISION HISTORY: 
!  16 Oct 2012 - R. Yantosca - Initial version
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_chm_mod.F90
!  19 Oct 2012 - R. Yantosca - Now reference gigc_state_met_mod.F90
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_Finalize
!  02 Nov 2012 - R. Yantosca - Now reference gigc_input_opt_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Assume success
    RC = GIGC_SUCCESS

    ! Deallocate fields of the Input Options object
    CALL Cleanup_GIGC_Input_Opt( am_I_Root, Input_Opt, RC )

    ! Deallocate fields of the Chemistry State object
    CALL Cleanup_GIGC_State_Chm( am_I_Root, State_Chm, RC )

    ! Deallocate fields of the Meteorology State object
    CALL Cleanup_GIGC_State_Met( am_I_Root, State_Met, RC )

    ! Deallocate all other GEOS-Chem allocatable arrays
    CALL GIGC_CleanUp_GeosChem( am_I_Root, RC )

  END SUBROUTINE GIGC_Finalize
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gigc_cleanup_geoschem
!
! !DESCRIPTION: Routine GIGC\_CLEANUP\_GEOSCHEM deallocates arrays in the 
!  various GEOS-Chem modules.  This is an analog to subroutine CLEANUP in 
!  GeosCore.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GIGC_CleanUp_GeosChem( am_I_Root, RC )
!
! !USES:
!
    USE AEROSOL_MOD,             ONLY : CLEANUP_AEROSOL
    USE AIRCRAFT_NOX_MOD,        ONLY : CLEANUP_AIRCRAFT_NOX
    USE ARCTAS_SHIP_EMISS_MOD,   ONLY : CLEANUP_ARCTAS_SHIP
    USE BIOMASS_MOD,             ONLY : CLEANUP_BIOMASS
    USE BIOFUEL_MOD,             ONLY : CLEANUP_BIOFUEL
    USE BRAVO_MOD,               ONLY : CLEANUP_BRAVO
    USE C2H6_MOD,                ONLY : CLEANUP_C2H6
    USE CAC_ANTHRO_MOD,          ONLY : CLEANUP_CAC_ANTHRO
    USE CARBON_MOD,              ONLY : CLEANUP_CARBON
    USE CO2_MOD,                 ONLY : CLEANUP_CO2
    USE COMODE_MOD,              ONLY : CLEANUP_COMODE
    USE GCKPP_COMODE_MOD,        ONLY : CLEANUP_GCKPP_COMODE
    USE DAO_MOD,                 ONLY : CLEANUP_DAO
    USE DIAG_MOD,                ONLY : CLEANUP_DIAG
    USE DIAG03_MOD,              ONLY : CLEANUP_DIAG03
    USE DIAG04_MOD,              ONLY : CLEANUP_DIAG04
    USE DIAG41_MOD,              ONLY : CLEANUP_DIAG41
    USE DIAG50_MOD,              ONLY : CLEANUP_DIAG50
    USE DIAG51_MOD,              ONLY : CLEANUP_DIAG51
    USE DIAG_OH_MOD,             ONLY : CLEANUP_DIAG_OH
    USE DIAG_PL_MOD,             ONLY : CLEANUP_DIAG_PL
    USE DRYDEP_MOD,              ONLY : CLEANUP_DRYDEP
    USE DUST_MOD,                ONLY : CLEANUP_DUST
    USE DUST_DEAD_MOD,           ONLY : CLEANUP_DUST_DEAD
    USE EDGAR_MOD,               ONLY : CLEANUP_EDGAR
    USE EMEP_MOD,                ONLY : CLEANUP_EMEP
    USE EPA_NEI_MOD,             ONLY : CLEANUP_EPA_NEI
    USE ERROR_MOD,               ONLY : DEBUG_MSG
    USE GC_BIOMASS_MOD,          ONLY : CLEANUP_GC_BIOMASS
    USE GFED2_BIOMASS_MOD,       ONLY : CLEANUP_GFED2_BIOMASS
    USE GFED3_BIOMASS_MOD,       ONLY : CLEANUP_GFED3_BIOMASS
    USE GIGC_ErrCode_Mod
    USE GLOBAL_CH4_MOD,          ONLY : CLEANUP_GLOBAL_CH4
    USE GLOBAL_HNO3_MOD,         ONLY : CLEANUP_GLOBAL_HNO3
    USE GLOBAL_NO3_MOD,          ONLY : CLEANUP_GLOBAL_NO3
    USE GLOBAL_NOX_MOD,          ONLY : CLEANUP_GLOBAL_NOX
    USE GLOBAL_O1D_MOD,          ONLY : CLEANUP_GLOBAL_O1D
    USE GLOBAL_OH_MOD,           ONLY : CLEANUP_GLOBAL_OH
    USE GLOBAL_GRID_MOD,         ONLY : CLEANUP_GLOBAL_GRID
    USE GRID_MOD,                ONLY : CLEANUP_GRID
    USE H2_HD_MOD,               ONLY : CLEANUP_H2_HD
    USE HCN_CH3CN_MOD,           ONLY : CLEANUP_HCN_CH3CN
    USE HDF_MOD,                 ONLY : CLEANUP_HDF
    USE ISOROPIAII_MOD,          ONLY : CLEANUP_ISOROPIAII
    USE LIGHTNING_NOX_MOD,       ONLY : CLEANUP_LIGHTNING_NOX
    USE LINOZ_MOD,               ONLY : CLEANUP_LINOZ
    USE MEGAN_MOD,               ONLY : CLEANUP_MEGAN
    USE MERCURY_MOD,             ONLY : CLEANUP_MERCURY
    USE MODIS_LAI_MOD,           ONLY : CLEANUP_MODIS_LAI
    USE OCEAN_MERCURY_MOD,       ONLY : CLEANUP_OCEAN_MERCURY
    USE DEPO_MERCURY_MOD,        ONLY : CLEANUP_DEPO_MERCURY
    USE LAND_MERCURY_MOD,        ONLY : CLEANUP_LAND_MERCURY
    USE PBL_MIX_MOD,             ONLY : CLEANUP_PBL_MIX
    USE PJC_PFIX_MOD,            ONLY : CLEANUP_PJC_PFIX
    USE PLANEFLIGHT_MOD,         ONLY : CLEANUP_PLANEFLIGHT
    USE PRESSURE_MOD,            ONLY : CLEANUP_PRESSURE
    USE REGRID_1x1_MOD,          ONLY : CLEANUP_REGRID_1x1
    USE SEASALT_MOD,             ONLY : CLEANUP_SEASALT
    USE SULFATE_MOD,             ONLY : CLEANUP_SULFATE
    USE STRAT_CHEM_MOD,          ONLY : CLEANUP_STRAT_CHEM
    USE TAGGED_CO_MOD,           ONLY : CLEANUP_TAGGED_CO
    USE TOMS_MOD,                ONLY : CLEANUP_TOMS
    USE TRACER_MOD,              ONLY : CLEANUP_TRACER
    USE TROPOPAUSE_MOD,          ONLY : CLEANUP_TROPOPAUSE
    USE UVALBEDO_MOD,            ONLY : CLEANUP_UVALBEDO
    USE VISTAS_ANTHRO_MOD,       ONLY : CLEANUP_VISTAS_ANTHRO
    USE WETSCAV_MOD,             ONLY : CLEANUP_WETSCAV
    USE ICOADS_SHIP_MOD,         ONLY : CLEANUP_ICOADS_SHIP  !(cklee,7/09/09)
    USE RETRO_MOD,               ONLY : CLEANUP_RETRO
    USE BROMOCARB_MOD,           ONLY : CLEANUP_BROMOCARB    ! jpp, 6/17/09

    IMPLICIT NONE
#   include "define.h"
!
! !INPUT PARAMETERS: 
!
    LOGICAL,         INTENT(IN)    :: am_I_Root
!
! !OUTPUT PARAMETERS:
!
    INTEGER,         INTENT(OUT)   :: RC
!
! !REMARKS:
!  Add calls to cleanup routines as needed.
! 
! !REVISION HISTORY: 
!  15 Oct 2012 - R. Yantosca - Initial version
!  22 Oct 2012 - R. Yantosca - Renamed to GIGC_CleanUp_GeosChem
!  22 Oct 2012 - R. Yantosca - Now reference gigc_errcode_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    
    ! Assume success at first
    RC = GIGC_SUCCESS

    ! Call cleanup routines from individual F90 modules
    CALL CLEANUP_AEROSOL
    CALL CLEANUP_AIRCRAFT_NOX
    CALL CLEANUP_ARCTAS_SHIP
    CALL CLEANUP_BIOMASS
    CALL CLEANUP_BIOFUEL
    CALL CLEANUP_BRAVO
    CALL CLEANUP_BROMOCARB  ! jpp, 6/17/09
    CALL CLEANUP_C2H6
    CALL CLEANUP_CAC_ANTHRO
    CALL CLEANUP_CARBON
    CALL CLEANUP_CO2
    CALL CLEANUP_COMODE
    CALL CLEANUP_GCKPP_COMODE
    CALL CLEANUP_DAO
    CALL CLEANUP_DIAG
    CALL CLEANUP_DIAG03
    CALL CLEANUP_DIAG04
    CALL CLEANUP_DIAG41
    CALL CLEANUP_DIAG50
    CALL CLEANUP_DIAG51
    CALL CLEANUP_DIAG_OH
    CALL CLEANUP_DIAG_PL
    CALL CLEANUP_DRYDEP
    CALL CLEANUP_DUST_DEAD
    CALL CLEANUP_DUST
    CALL CLEANUP_EDGAR
    CALL CLEANUP_EMEP
    CALL CLEANUP_EPA_NEI
    CALL CLEANUP_GC_BIOMASS
    CALL CLEANUP_GFED2_BIOMASS
    CALL CLEANUP_GFED3_BIOMASS
    CALL CLEANUP_GLOBAL_CH4
    CALL CLEANUP_GLOBAL_HNO3
    CALL CLEANUP_GLOBAL_NO3
    CALL CLEANUP_GLOBAL_NOX
    CALL CLEANUP_GLOBAL_NO3
    CALL CLEANUP_GLOBAL_O1D
    CALL CLEANUP_GLOBAL_OH
    CALL CLEANUP_GLOBAL_GRID
    CALL CLEANUP_GRID
    CALL CLEANUP_H2_HD
    CALL CLEANUP_HCN_CH3CN
    CALL CLEANUP_HDF
    CALL CLEANUP_ISOROPIAII
    CALL CLEANUP_LIGHTNING_NOX
    CALL CLEANUP_LINOZ
    CALL CLEANUP_MEGAN
    CALL CLEANUP_MERCURY
    CALL CLEANUP_MODIS_LAI
    CALL CLEANUP_OCEAN_MERCURY
    CALL CLEANUP_DEPO_MERCURY
    CALL CLEANUP_LAND_MERCURY
    CALL CLEANUP_PBL_MIX
    CALL CLEANUP_PJC_PFIX
    CALL CLEANUP_PLANEFLIGHT
    CALL CLEANUP_PRESSURE
    CALL CLEANUP_REGRID_1x1
    CALL CLEANUP_SEASALT
    CALL CLEANUP_SULFATE
    CALL CLEANUP_STRAT_CHEM
    CALL CLEANUP_TAGGED_CO
    CALL CLEANUP_TOMS
    CALL CLEANUP_TRACER
    CALL CLEANUP_TROPOPAUSE
    CALL CLEANUP_UVALBEDO
    CALL CLEANUP_VISTAS_ANTHRO
    CALL CLEANUP_WETSCAV
    CALL CLEANUP_ICOADS_SHIP !(cklee,7/09/09)
    CALL CLEANUP_RETRO
    
  END SUBROUTINE GIGC_CleanUp_GeosChem
!EOC
END MODULE GIGC_Finalization_Mod
#endif
