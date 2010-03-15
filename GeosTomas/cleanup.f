! $Id: cleanup.f,v 1.2 2010/03/15 19:33:20 ccarouge Exp $
      SUBROUTINE CLEANUP
!
!******************************************************************************
!  Subroutine CLEANUP deallocates the memory assigned to dynamic allocatable 
!  arrays just before exiting the GEOS-CHEM model. (bmy, 11/29/99, 1/25/10)
!
!  NOTES:
!  (1 ) CLEANUP is written in Fixed-Format F90.
!  (2 ) Now calls CLEANUP_WETSCAV, which deallocates arrays from 
!        "wetscav_mod.f". (bmy, 3/9/00)
!  (3 ) Add call to CLEANUP_SULFATE, which deallocates arrays from
!        "sulfate_mod.f".  Also now deallocate ND32 arrays. (bmy, 6/6/00)
!  (4 ) Add call to CLEANUP_DAO, which deallocates arrays from "dao_mod.f".  
!        (bmy, 6/26/00)
!  (5 ) Add call to CLEANUP_TAGGED_CO and CLEANUP_COMODE, which deallocates 
!        arrays from and "comode_mod.f". (bmy, 7/19/00)
!  (6 ) Add call to CLEANUP_GLOBAL_OH and CLEANUP_COMODE, which deallocates 
!        arrays from "global_oh_mod.f". (bmy, 7/28/00)
!  (7 ) Add calls to CLEANUP_BIOMASS and CLEANUP_BIOFUEL, which deallocates 
!        arrays from "biomass_mod.f" and "biofuel_mod.f".  Also deallocate
!        the AD32_bf array for the biofuel NOx diagnostic. (bmy, 9/12/00)
!  (8 ) Add call to CLEANUP_DIAG51, to deallocate module arrays from
!        "diag51_mod.f" (bmy, 11/29/00)
!  (9 ) Removed obsolete code from 11/29/00 (bmy, 12/21/00)
!  (10) Add call to CLEANUP_CH4, to deallocate module arrays from
!        "global_ch4_mod.f" (bmy, 1/16/01)
!  (11) Now deallocate the AD34 array.  Also updated comments and
!        made some cosmetic changes. (bmy, 3/15/01)
!  (12) Now deallocate the AD12 array (bdf, bmy, 6/15/01)
!  (13) Add call to CLEANUP_ACETONE, to deallocate module arrays from 
!        "acetone_mod.f"  Also deallocate AD11 array.  Also deallocate 
!        variables from dao_mod.f last, to try to avoid bus error on 
!        SGI (bmy, 8/3/01) 
!  (14) Added call to CLEANUP_UVALBEDO from "uvalbedo_mod.f".  Also removed
!        obsolete code from 9/01.  Also only include references to CLEANUP_* 
!        subroutines in other modules for clarity. (bmy, 1/15/02)
!  (15) Added call to CLEANUP_C2H6 from "c2h6_mod.f" (bmy, 1/25/02)
!  (16) Added call to CLEANUP_AIRCRAFT_NOX from "aircraft_nox_mod.f" 
!        (bmy, 2/14/02)
!  (17) Now deallocate CTNO2, CTHO2, LTNO2, LTHO2 arrays (rvm, bmy, 2/27/02)
!  (18) Now reference CLEANUP_PLANEFLIGHT from "planeflight_mod.f".
!        Now also deallocate AD01 and AD02 arrays. (mje, bmy, 8/7/02)
!  (19) Now reference cleanup routines from "global_nox_mod.f", 
!        "global_hno3_mod.f", "global_no3_mod.f", "drydep_mod.f", and
!        "rpmares_mod.f". (bmy, 12/16/02)
!  (20) Now reference cleanup routine from "transport_mod.f" (bmy, 2/10/03)
!  (21) Now reference cleanup routine from "pjc_pfix_mod.f" and 
!        "tpcore_fvdas_mod.f90". (bmy, 5/9/03)
!  (22) Now reference cleanup routine from "toms_mod.f" (bmy, 7/14/03)
!  (23) Now reference cleanup routine from "carbon_mod.f", "dust_mod.f", and
!        "dust_dead_mod.f". (bmy, 7/14/03)
!  (23) Now references cleanup routine from "lightning__nox_mod.f" 
!        (bmy, 4/14/04)
!  (24) Now references cleanup routine from "seasalt_mod.f" (bmy, 4/26/04)
!  (25) Now references cleanup routines from new modules (bmy, 7/20/04)
!  (26) Now calls cleanup routine from "epa_nei_mod.f" (bmy, 11/5/04)
!  (27) Now call CLEANUP_MERCURY from "mercury_mod.f" (eck, bmy, 12/7/04)
!  (28) Now call CLEANUP_OCEAN_MERCURY from "ocean_mercury_mod.f".  Also
!        reordered the calling sequence. (sas, bmy, 1/21/05)
!  (29) Now call CLEANUP_PBL_MIX from "pbl_mix_mod.f".  Now call CLEANUP_DIAG41
!        from "diag41_mod.f". (bmy, 2/17/05)
!  (30) Now calls CLEANUP_HCN_CH3CN from "hcn_ch3cn_mod.f (bmy, 6/23/05)
!  (31) Now calls CLEANUP_DIAG04, CLEANUP_CO2, and CLEANUP_TROPOPAUSE 
!         (bmy, 8/15/05) 
!  (32) Now calls CLEANUP_LAI from "lai_mod.f", CLEANUP_MEGAN from
!        "megan_mod.f" and CLEANUP_REGRID_1x1 from "regrid_1x1_mod.f"
!        (tmf, bdf, bmy, 10/24/05)
!  (33) Now calls CLEANUP_EMEP from "emep_mod.f" (bdf, bmy, 11/1/05)
!  (34) Now calls CLEANUP_GC_BIOMASS and CLEANUP_GFED2_BIOMASS (bmy, 4/5/06)
!  (35) Now calls CLEANUP_DIAG56 from "diag56_mod.f" and
!        CLEANUP_LIGHTNING_NOX_NL from "lightning_nox_nl_mod.f"
!        (ltm, bmy, 5/5/06)
!  (36) Now references CLEANUP_BRAVO from "bravo_mod.f" and CLEANUP_EDGAR
!        from "edgar_mod.f" (bmy, 7/6/06)
!  (37) Now calls CLEANUP_H2_HD from "h2_hd_mod.f" and CLEANUP_GLOBAL_O1D
!        from "global_o1d_mod.f".  Remove call to CLEANUP_LIGHTNING_NOx_NL 
!        from "lightning_nox_nl_mod.f (hup, phs, bmy, 10/2/07)
!  (38) Now calls GEOS5_EXIT_TPCORE_WINDOW to finalize the TPCORE for
!        GEOS-5 nested window simulations (yxw, dan, bmy, 11/6/08)
!  (39) Now references CLEANUP_CAC_ANTHRO (amv, phs, 3/10/08)
!  (40) Now references CLEANUP_ARCTAS_SHIP (phs, 3/10/08)
!  (41) Now references CLEANUP_VISTAS_ANTHRO (phs, 3/10/08)
!  (41) Now references CLEANUP_LINOZ (phs, 10/16/09)
!  (42) Now references CLEANUP_HDF (amv, bmy, 12/21/09)
!  (43) Now references CLEANUP_TOMAS (win, bmy, 1/25/10)
!******************************************************************************
!
      ! References to F90 modules 
      USE ACETONE_MOD,             ONLY : CLEANUP_ACETONE
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
      USE GLOBAL_CH4_MOD,          ONLY : CLEANUP_GLOBAL_CH4
      USE GLOBAL_HNO3_MOD,         ONLY : CLEANUP_GLOBAL_HNO3
      USE GLOBAL_NO3_MOD,          ONLY : CLEANUP_GLOBAL_NO3
      USE GLOBAL_NOX_MOD,          ONLY : CLEANUP_GLOBAL_NOX
      USE GLOBAL_O1D_MOD,          ONLY : CLEANUP_GLOBAL_O1D
      USE GLOBAL_OH_MOD,           ONLY : CLEANUP_GLOBAL_OH
      USE H2_HD_MOD,               ONLY : CLEANUP_H2_HD
      USE HCN_CH3CN_MOD,           ONLY : CLEANUP_HCN_CH3CN
      USE HDF_MOD,                 ONLY : CLEANUP_HDF
      USE ISOROPIAII_MOD,          ONLY : CLEANUP_ISOROPIAII
      USE LAI_MOD,                 ONLY : CLEANUP_LAI
      USE LIGHTNING_NOX_MOD,       ONLY : CLEANUP_LIGHTNING_NOX
      USE LINOZ_MOD,               ONLY : CLEANUP_LINOZ
      USE MEGAN_MOD,               ONLY : CLEANUP_MEGAN
      USE MERCURY_MOD,             ONLY : CLEANUP_MERCURY
      USE OCEAN_MERCURY_MOD,       ONLY : CLEANUP_OCEAN_MERCURY
      USE PBL_MIX_MOD,             ONLY : CLEANUP_PBL_MIX
      USE PJC_PFIX_MOD,            ONLY : CLEANUP_PJC_PFIX
      USE PLANEFLIGHT_MOD,         ONLY : CLEANUP_PLANEFLIGHT
      USE PRESSURE_MOD,            ONLY : CLEANUP_PRESSURE
      USE REGRID_1x1_MOD,          ONLY : CLEANUP_REGRID_1x1
      USE SEASALT_MOD,             ONLY : CLEANUP_SEASALT
      USE SULFATE_MOD,             ONLY : CLEANUP_SULFATE
      USE TAGGED_CO_MOD,           ONLY : CLEANUP_TAGGED_CO
      USE TOMAS_MOD,               ONLY : CLEANUP_TOMAS       !(win, 1/25/10)
      USE TOMS_MOD,                ONLY : CLEANUP_TOMS
      USE TPCORE_FVDAS_MOD,        ONLY : EXIT_TPCORE
      USE TPCORE_GEOS5_WINDOW_MOD, ONLY : EXIT_GEOS5_TPCORE_WINDOW
      USE TRACER_MOD,              ONLY : CLEANUP_TRACER
      USE TRANSPORT_MOD,           ONLY : CLEANUP_TRANSPORT
      USE TROPOPAUSE_MOD,          ONLY : CLEANUP_TROPOPAUSE
      USE UVALBEDO_MOD,            ONLY : CLEANUP_UVALBEDO
      USE VISTAS_ANTHRO_MOD,       ONLY : CLEANUP_VISTAS_ANTHRO
      USE WETSCAV_MOD,             ONLY : CLEANUP_WETSCAV
      USE ICOADS_SHIP_MOD,         ONLY : CLEANUP_ICOADS_SHIP  !(cklee,7/09/09)

      IMPLICIT NONE

#     include "define.h"

      !=================================================================
      ! CLEANUP begins here!
      !=================================================================

      ! Echo info
      WRITE( 6, 100 ) 
 100  FORMAT( '     - CLEANUP: deallocating arrays now...' )

      ! Call cleanup routines from individual F90 modules
      CALL CLEANUP_ACETONE
      CALL CLEANUP_AEROSOL
      CALL CLEANUP_AIRCRAFT_NOX
      CALL CLEANUP_ARCTAS_SHIP
      CALL CLEANUP_BIOMASS
      CALL CLEANUP_BIOFUEL
      CALL CLEANUP_BRAVO
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
      CALL CLEANUP_GLOBAL_CH4
      CALL CLEANUP_GLOBAL_HNO3
      CALL CLEANUP_GLOBAL_NO3
      CALL CLEANUP_GLOBAL_NOX
      CALL CLEANUP_GLOBAL_NO3
      CALL CLEANUP_GLOBAL_O1D
      CALL CLEANUP_GLOBAL_OH
      CALL CLEANUP_H2_HD
      CALL CLEANUP_HCN_CH3CN
      CALL CLEANUP_HDF
      CALL CLEANUP_ISOROPIAII
      CALL CLEANUP_LAI
      CALL CLEANUP_LIGHTNING_NOX
      CALL CLEANUP_LINOZ
      CALL CLEANUP_MEGAN
      CALL CLEANUP_MERCURY
      CALL CLEANUP_OCEAN_MERCURY
      CALL CLEANUP_PBL_MIX
      CALL CLEANUP_PJC_PFIX
      CALL CLEANUP_PLANEFLIGHT
      CALL CLEANUP_PRESSURE
      CALL CLEANUP_REGRID_1x1
      CALL CLEANUP_SEASALT
      CALL CLEANUP_SULFATE
      CALL CLEANUP_TAGGED_CO
      CALL CLEANUP_TRANSPORT
      CALL CLEANUP_TOMS
      CALL CLEANUP_TRACER
      CALL CLEANUP_TROPOPAUSE
      CALL CLEANUP_UVALBEDO
      CALL CLEANUP_VISTAS_ANTHRO
      CALL CLEANUP_WETSCAV
      CALL CLEANUP_ICOADS_SHIP !(cklee,7/09/09)
      !if ( IDTNK1 > 0 ) CALL CLEANUP_TOMAS !(win, 7/28/09)
      CALL CLEANUP_TOMAS

#if   defined( GEOS_5 ) && defined( GRID05x0666 )
      CALL EXIT_GEOS5_TPCORE_WINDOW 
#else
      CALL EXIT_TPCORE
#endif

      ! Return to calling program
      END SUBROUTINE CLEANUP
