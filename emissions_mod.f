! $Id: emissions_mod.f,v 1.8 2004/12/16 16:52:45 bmy Exp $
      MODULE EMISSIONS_MOD
!
!******************************************************************************
!  Module EMISSIONS_MOD is used to call the proper emissions subroutine
!  for the various GEOS-CHEM simulations. (bmy, 2/11/03, 12/7/04)
! 
!  Module Routines:
!  ============================================================================
!  (1 ) DO_EMISSIONS     : Driver which calls various emissions routines
!
!  GEOS-CHEM modules referenced by emissions_mod.f
!  ============================================================================
!  (1 ) c2h6_mod.f       : Module containing routines for C2H6 chemistry
!  (2 ) carbon_mod.f     : Module containing routines for carbon arsl emissions
!  (3 ) ch3i_mod.f       : Module containing routines for CH3I chemistry
!  (4 ) dust_mod.f       : Module containing routines for dust arsl emissions
!  (5 ) epa_nei_mod.f    : Module containing routines to read EPA/NEI99 data
!  (6 ) error_mod.f      : Module containing NaN and other error checks
!  (7 ) global_ch4_mod.f : Module containing routines for CH4 emissions
!  (8 ) Kr85_mod.f       : Module containing routines for Kr85 emissions
!  (9 ) logical_mod.f    : Module containing GEOS-CHEM logical switches
!  (10) mercury_mod.f    : Module containing routines for mercury chemistry
!  (11) RnPbBe_mod.f     : Module containing routines for Rn-Pb-Be emissions
!  (12) tagged_co_mod.f  : Module containing routines for Tagged CO emissions
!  (13) time_mod.f       : Module containing routines to compute date & time
!  (14) tracer_mod.f     : Module containing GEOS-CHEM tracer array STT etc.
!  (15) seasalt_mod.f    : Module containing routines for seasalt emissions
!  (16) sulfate_mod.f    : Module containing routines for sulfate emissions
!
!  NOTES:
!  (1 ) Now references DEBUG_MSG from "error_mod.f"
!  (2 ) Now references "Kr85_mod.f" (jsw, bmy, 8/20/03)
!  (3 ) Now references "carbon_mod.f" and "dust_mod.f" (rjp, tdf, bmy, 4/2/04)
!  (4 ) Now references "seasalt_mod.f" (rjp, bmy, bec, 4/20/04)
!  (5 ) Now references "logical_mod" & "tracer_mod.f" (bmy, 7/20/04)
!  (6 ) Now references "epa_nei_mod.f" and "time_mod.f" (bmy, 11/5/04)
!  (7 ) Now references "emissions_mod.f" (bmy, 12/7/04)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
      
      SUBROUTINE DO_EMISSIONS
!
!******************************************************************************
!  Subroutine DO_EMISSIONS is the driver routine which calls the appropriate
!  emissions subroutine for the various GEOS-CHEM simulations. 
!  (bmy, 2/11/03, 12/7/04)
!
!  NOTES:
!  (1 ) Now references DEBUG_MSG from "error_mod.f" (bmy, 8/7/03)
!  (2 ) Now calls Kr85 emissions if NSRCX == 12 (jsw, bmy, 8/20/03)
!  (3 ) Now calls EMISSCARBON and EMISSDUST for carbon aerosol and dust
!        aerosol chemistry (rjp, tdf, bmy, 4/2/04)
!  (4 ) Now calls EMISSSEASALT for seasalt emissions (rjp, bec, bmy, 4/20/04)
!  (5 ) Now use inquiry functions from "tracer_mod.f".  Now references
!        "logical_mod.f" (bmy, 7/20/04)
!  (6 ) Now references ITS_A_NEW_MONTH from "time_mod.f".  Now references
!        EMISS_EPA_NEI from "epa_nei_mod.f" (bmy, 11/5/04)
!  (7 ) Now calls EMISSMERCURY from "mercury_mod.f" (eck, bmy, 12/7/04)
!******************************************************************************
!
      ! References to F90 modules
      USE C2H6_MOD,       ONLY : EMISSC2H6
      USE CARBON_MOD,     ONLY : EMISSCARBON
      USE CH3I_MOD,       ONLY : EMISSCH3I
      USE DUST_MOD,       ONLY : EMISSDUST
      USE EPA_NEI_MOD,    ONLY : EMISS_EPA_NEI
      USE ERROR_MOD,      ONLY : DEBUG_MSG
      USE GLOBAL_CH4_MOD, ONLY : EMISSCH4
      USE Kr85_MOD,       ONLY : EMISSKr85
      USE LOGICAL_MOD
      USE MERCURY_MOD,    ONLY : EMISSMERCURY
      USE RnPbBe_MOD,     ONLY : EMISSRnPbBe
      USE SEASALT_MOD,    ONLY : EMISSSEASALT
      USE SULFATE_MOD,    ONLY : EMISSSULFATE 
      USE TIME_MOD,       ONLY : ITS_A_NEW_MONTH
      USE TRACER_MOD
      USE TAGGED_CO_MOD,  ONLY : EMISS_TAGGED_CO

#     include "CMN_SIZE"  ! Size parameters

      !=================================================================
      ! DO_EMISSIONS begins here!
      !=================================================================
      IF ( ITS_A_FULLCHEM_SIM() ) THEN

         ! Read EPA/NEI99 emissions once per month
         IF ( LNEI99 .and. ITS_A_NEW_MONTH() ) CALL EMISS_EPA_NEI

         ! NOx-Ox-HC (w/ or w/o aerosols)
         CALL EMISSDR

         ! Emissions for various aerosol types
         IF ( LSULF  ) CALL EMISSSULFATE
         IF ( LCARB  ) CALL EMISSCARBON
         IF ( LDUST  ) CALL EMISSDUST
         IF ( LSSALT ) CALL EMISSSEASALT

      ELSE IF ( ITS_A_RnPbBe_SIM() ) THEN
         
         ! Rn-Pb-Be
         CALL EMISSRnPbBe

      ELSE IF ( ITS_A_CH3I_SIM() ) THEN

         ! CH3I
         CALL EMISSCH3I

      ELSE IF ( ITS_A_HCN_SIM() ) THEN

         ! HCN - CH3CN
         CALL EMISSHCN

      ELSE IF ( ITS_A_TAGCO_SIM() ) THEN

         ! Tagged CO
         CALL EMISS_TAGGED_CO

      ELSE IF ( ITS_A_C2H6_SIM() ) THEN

         ! C2H6
         CALL EMISSC2H6

      ELSE IF ( ITS_A_CH4_SIM() ) THEN

         ! CH4
         CALL EMISSCH4

      ELSE IF ( ITS_A_MERCURY_SIM() ) THEN

         ! Mercury
         CALL EMISSMERCURY

      ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

         ! Offline Sulfate, carbon, & dust
         IF ( LSULF  ) CALL EMISSSULFATE
         IF ( LCARB  ) CALL EMISSCARBON
         IF ( LDUST  ) CALL EMISSDUST
         IF ( LSSALT ) CALL EMISSSEASALT

!-----------------------------------------------------------------------------
! Prior to 2/11/03:
! Reinstate this later in a more integrated way (bmy, 2/11/03)
!#if   defined( LGEOSCO )
!      ELSE IF ( ITS_A_COPARAM_SIM() ) THEN
!
!         ! CO w/ parameterized OH
!         CALL EMISSCO( FIRSTEMISS, NSEASON, LMN, SUNCOS )
!
!#endif
!-----------------------------------------------------------------------------

      ENDIF

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG ( '### DO_EMISSIONS: a EMISSIONS' )

      ! Return to calling program
      END SUBROUTINE DO_EMISSIONS

!------------------------------------------------------------------------------

      END MODULE EMISSIONS_MOD
