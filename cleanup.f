! $Id: cleanup.f,v 1.8 2004/12/16 16:52:44 bmy Exp $
      SUBROUTINE CLEANUP
!
!******************************************************************************
!  Subroutine CLEANUP deallocates the memory assigned to dynamic allocatable 
!  arrays just before exiting the GEOS-CHEM model. (bmy, 11/29/99, 12/7/04)
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
!******************************************************************************
!
      ! References to F90 modules 
      USE ACETONE_MOD,       ONLY : CLEANUP_ACETONE
      USE AEROSOL_MOD,       ONLY : CLEANUP_AEROSOL
      USE AIRCRAFT_NOX_MOD,  ONLY : CLEANUP_AIRCRAFT_NOX
      USE BIOMASS_MOD,       ONLY : CLEANUP_BIOMASS
      USE BIOFUEL_MOD,       ONLY : CLEANUP_BIOFUEL
      USE C2H6_MOD,          ONLY : CLEANUP_C2H6
      USE CARBON_MOD,        ONLY : CLEANUP_CARBON
      USE COMODE_MOD,        ONLY : CLEANUP_COMODE
      USE DAO_MOD,           ONLY : CLEANUP_DAO
      USE DIAG_MOD,          ONLY : CLEANUP_DIAG
      USE DIAG50_MOD,        ONLY : CLEANUP_DIAG50
      USE DIAG51_MOD,        ONLY : CLEANUP_DIAG51
      USE DIAG_OH_MOD,       ONLY : CLEANUP_DIAG_OH
      USE DIAG_PL_MOD,       ONLY : CLEANUP_DIAG_PL
      USE DRYDEP_MOD,        ONLY : CLEANUP_DRYDEP
      USE DUST_MOD,          ONLY : CLEANUP_DUST
      USE DUST_DEAD_MOD,     ONLY : CLEANUP_DUST_DEAD
      USE EPA_NEI_MOD,       ONLY : CLEANUP_EPA_NEI
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE GLOBAL_CH4_MOD,    ONLY : CLEANUP_GLOBAL_CH4
      USE GLOBAL_HNO3_MOD,   ONLY : CLEANUP_GLOBAL_HNO3
      USE GLOBAL_NO3_MOD,    ONLY : CLEANUP_GLOBAL_NO3
      USE GLOBAL_NOX_MOD,    ONLY : CLEANUP_GLOBAL_NOX
      USE GLOBAL_OH_MOD,     ONLY : CLEANUP_GLOBAL_OH
      USE LIGHTNING_NOX_MOD, ONLY : CLEANUP_LIGHTNING_NOX
      USE MERCURY_MOD,       ONLY : CLEANUP_MERCURY
      USE PJC_PFIX_MOD,      ONLY : CLEANUP_PJC_PFIX
      USE PLANEFLIGHT_MOD,   ONLY : CLEANUP_PLANEFLIGHT
      USE PRESSURE_MOD,      ONLY : CLEANUP_PRESSURE
      USE SEASALT_MOD,       ONLY : CLEANUP_SEASALT
      USE SULFATE_MOD,       ONLY : CLEANUP_SULFATE
      USE TAGGED_CO_MOD,     ONLY : CLEANUP_TAGGED_CO
      USE TOMS_MOD,          ONLY : CLEANUP_TOMS
      USE TPCORE_FVDAS_MOD,  ONLY : EXIT_TPCORE
      USE TRANSPORT_MOD,     ONLY : CLEANUP_TRANSPORT
      USE UVALBEDO_MOD,      ONLY : CLEANUP_UVALBEDO
      USE WETSCAV_MOD,       ONLY : CLEANUP_WETSCAV

      IMPLICIT NONE

      !=================================================================
      ! CLEANUP begins here!
      !=================================================================
      WRITE( 6, '(a)' ) '     - CLEANUP: deallocating arrays now...'

      !=================================================================
      ! Deallocate arrays in other modules (the order may be important)
      !=================================================================
      CALL EXIT_TPCORE
      CALL CLEANUP_PLANEFLIGHT
      CALL CLEANUP_PJC_PFIX
      CALL CLEANUP_TAGGED_CO
      CALL CLEANUP_TRANSPORT
      CALL CLEANUP_SULFATE
      CALL CLEANUP_DIAG50
      CALL CLEANUP_DIAG51
      CALL CLEANUP_DIAG_OH
      CALL CLEANUP_DIAG_PL
      CALL CLEANUP_TOMS
      CALL CLEANUP_GLOBAL_CH4
      CALL CLEANUP_GLOBAL_HNO3
      CALL CLEANUP_GLOBAL_NO3
      CALL CLEANUP_GLOBAL_NOX
      CALL CLEANUP_GLOBAL_NO3
      CALL CLEANUP_GLOBAL_OH
      CALL CLEANUP_DRYDEP
      CALL CLEANUP_EPA_NEI
      CALL CLEANUP_BIOMASS
      CALL CLEANUP_BIOFUEL
      CALL CLEANUP_AEROSOL
      CALL CLEANUP_ACETONE
      CALL CLEANUP_AIRCRAFT_NOX
      CALL CLEANUP_C2H6
      CALL CLEANUP_CARBON
      CALL CLEANUP_COMODE
      CALL CLEANUP_DUST_DEAD
      CALL CLEANUP_DUST
      CALL CLEANUP_LIGHTNING_NOX
      CALL CLEANUP_MERCURY
      CALL CLEANUP_SEASALT
      CALL CLEANUP_UVALBEDO
      CALL CLEANUP_WETSCAV
      CALL CLEANUP_PRESSURE
      CALL CLEANUP_DIAG
      CALL CLEANUP_DAO

      ! Return to calling program
      END SUBROUTINE CLEANUP
