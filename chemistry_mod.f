! $Id: chemistry_mod.f,v 1.8 2004/05/03 14:46:15 bmy Exp $
      MODULE CHEMISTRY_MOD
!
!******************************************************************************
!  Module CHEMISTRY_MOD is used to call the proper chemistry subroutine
!  for the various GEOS-CHEM simulations. (bmy, 4/14/03, 4/20/04)
! 
!  Module Routines:
!  ============================================================================
!  (1 ) DO_CHEMISTRY       : Driver which calls various chemistry routines
!
!  GEOS-CHEM modules referenced by chemistry_mod.f
!  ============================================================================
!  (1 ) acetone_mod.f      : Module containing routines for ACET chemistry
!  (2 ) c2h6_mod.f         : Module containing routines for C2H6 chemistry
!  (3 ) carbon_mod.f       : Module containing routines for carbon arsl chem.
!  (4 ) ch3i_mod.f         : Module containing routines for CH3I chemistry
!  (5 ) dao_mod.f          : Module containing arrays for DAO met fields
!  (6 ) drydep_mod.f       : Module containing GEOS-CHEM drydep routines
!  (7 ) dust_mod.f         : Module containing routines for dust arsl chem.
!  (8 ) error_mod.f        : Module containing NaN and error checks
!  (9 ) global_ch4_mod.f   : Module containing routines for CH4 chemistry
!  (10) Kr85_mod.f         : Module containing routines for Kr85 chemistry
!  (11) RnPbBe_mod.f       : Module containing routines for Rn-Pb-Be chemistry
!  (12) rpmares_mod.f      : Module containing routines for arsl phase equilib.
!  (13) seasalt_mod.f      : Module containing routines for seasalt chemistry
!  (13) sulfate_mod.f      : Module containing routines for sulfate chemistry
!  (14) tagged_co_mod.f    : Module containing routines for Tagged CO chemistry
!  (15) tagged_ox_mod.f    : Module containing routines for Tagged Ox chemistry
!  (16) time_mod.f         : Module containing routines to compute time & date
!  (17) tracerid_mod.f     : Module containing pointers to tracers & emissions
!
!  NOTES:
!  (1 ) Bug fix in DO_CHEMISTRY (bnd, bmy, 4/14/03)
!  (2 ) Now references DEBUG_MSG from "error_mod.f" (bmy, 8/7/03)
!  (3 ) Now references "tagged_ox_mod.f"(bmy, 8/18/03)
!  (4 ) Now references "Kr85_mod.f" (jsw, bmy, 8/20/03)
!  (5 ) Bug fix: Now also call OPTDEPTH for GEOS-4 (bmy, 1/27/04)
!  (6 ) Now references "carbon_mod.f" and "dust_mod.f" (rjp, tdf, bmy, 4/5/04)
!  (7 ) Now references "seasalt_mod.f" (rjp, bec, bmy, 4/20/04)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
      
      SUBROUTINE DO_CHEMISTRY
!
!******************************************************************************
!  Subroutine DO_CHEMISTRY is the driver routine which calls the appropriate
!  chemistry subroutine for the various GEOS-CHEM simulations. 
!  (bmy, 2/11/03, 4/20/04)
!
!  NOTES:
!  (1 ) Now reference DELP, T from "dao_mod.f" since we need to pass this
!        to OPTDEPTH for GEOS-1 or GEOS-STRAT met fields (bnd, bmy, 4/14/03)
!  (2 ) Now references DEBUG_MSG from "error_mod.f" (bmy, 8/7/03)
!  (3 ) Removed call to CHEMO3, it's obsolete.  Now calls CHEM_TAGGED_OX !
!        from "tagged_ox_mod.f" when NSRCX==6.  Now calls Kr85 chemistry if 
!        NSRCX == 12 (jsw, bmy, 8/20/03)
!  (4 ) Bug fix: added GEOS-4 to the #if block in the call to OPTDEPTH.
!        (bmy, 1/27/04)
!  (5 ) Now calls CHEMCARBON and CHEMDUST to do carbon aerosol & dust 
!        aerosol chemistry (rjp, tdf, bmy, 4/2/04)
!  (5 ) Now calls CHEMSEASALT to do seasalt aerosol chemistry 
!        (rjp, bec, bmy, 4/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE ACETONE_MOD,     ONLY : OCEAN_SINK_ACET
      USE C2H6_MOD,        ONLY : CHEMC2H6
      USE CARBON_MOD,      ONLY : CHEMCARBON
      USE CH3I_MOD,        ONLY : CHEMCH3I
      USE DAO_MOD,         ONLY : CLMOSW,  CLROSW, DELP, 
     &                            MAKE_RH, OPTDEP, OPTD, T
      USE DRYDEP_MOD,      ONLY : DRYFLX, DRYFLXRnPbBe
      USE DUST_MOD,        ONLY : CHEMDUST
      USE ERROR_MOD,       ONLY : DEBUG_MSG
      USE GLOBAL_CH4_MOD,  ONLY : CHEMCH4
      USE Kr85_MOD,        ONLY : CHEMKr85
      USE OPTDEPTH_MOD,    ONLY : OPTDEPTH
      !USE PLANEFLIGHT_MOD, ONLY : PLANEFLIGHT
      USE RnPbBe_MOD,      ONLY : CHEMRnPbBe
      USE RPMARES_MOD,     ONLY : DO_RPMARES
      USE SEASALT_MOD,     ONLY : CHEMSEASALT
      USE SULFATE_MOD,     ONLY : CHEMSULFATE
      USE TAGGED_CO_MOD,   ONLY : CHEM_TAGGED_CO
      USE TAGGED_OX_MOD,   ONLY : CHEM_TAGGED_OX
      USE TIME_MOD,        ONLY : GET_ELAPSED_MIN, GET_TS_CHEM
      USE TRACERID_MOD,    ONLY : IDTACET

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! NSRCX
#     include "CMN_SETUP" ! LSULF, LCHEM
#     include "CMN_DIAG"  ! NDxx flags

      !=================================================================
      ! DO_CHEMISTRY begins here!
      !=================================================================

#if   defined( GEOS_1 ) || defined( GEOS_STRAT )

      ! Compute optical depths (except for CO-OH run)
      ! GEOS-1/GEOS-STRAT: compute from T, CLMO, CLRO, DELP
      IF ( NSRCX /= 5 .and. NSRCX /= 9 ) THEN
         CALL OPTDEPTH( LM, CLMOSW, CLROSW, DELP, T, OPTD )
      ENDIF

#elif defined( GEOS_3 ) || defined( GEOS_4 ) 

      ! Compute optical depths (except for CO-OH run)
      ! GEOS-2/GEOS-3: Copy OPTDEP to OPTD, also archive diagnostics
      IF ( NSRCX /= 5 .and. NSRCX /= 9 ) THEN
         CALL OPTDEPTH( LM, OPTDEP, OPTD )
      ENDIF
#endif 

      !=================================================================
      ! If LCHEM=T then call the chemistry subroutines
      !=================================================================
      IF ( LCHEM ) THEN 
                  
         ! Chemistry breakpoint
         SELECT CASE ( NSRCX )

            !---------------------------------
            ! Rn-Pb-Be
            !---------------------------------
            CASE ( 1  )
               CALL CHEMRnPbBe 
               CALL DRYFLXRnPbBe
                  
            !---------------------------------
            ! CH3I
            !---------------------------------
            CASE ( 2  ) 
               CALL CHEMCH3I

            !---------------------------------
            ! NOx-Ox-HC (w/ or w/o aerosols) 
            !---------------------------------
            CASE ( 3 )

               ! Call SMVGEAR routines
               CALL CHEMDR

               ! Also do sulfate chemistry
               IF ( LSULF ) THEN

                  ! Get relative humidity
                  CALL MAKE_RH

                  ! Do sulfate chemistry
                  CALL CHEMSULFATE

                  ! Do aerosol phase equilibrium
                  CALL DO_RPMARES
               ENDIF

               ! Do carbonaceous aerosol chemistry
               IF ( LCARB ) CALL CHEMCARBON

               ! Do dust aerosol chemistry
               IF ( LDUST ) CALL CHEMDUST

               ! Do seasalt aerosol chemistry
               IF ( LSSALT ) CALL CHEMSEASALT

               ! ND44 drydep fluxes
               CALL DRYFLX     

               ! ND43 chemical production
               CALL DIAGOH

               ! Remove acetone ocean sink
               IF ( IDTACET /= 0 ) THEN
                  CALL OCEAN_SINK_ACET( STT(:,:,1,IDTACET) ) 
               ENDIF

               ! ND65 P-L diagnostics
               IF ( ND65 > 0 ) CALL DIAGPL

               ! Save P(Ox) and L(Ox) for offline run
               IF ( ND20 > 0 ) CALL DIAG20

            !---------------------------------            
            ! HCN
            !---------------------------------
            CASE ( 4  ) 
               CALL CHEMHCN

            !---------------------------------
            ! Tagged O3
            !---------------------------------
            CASE ( 6  ) 
               CALL CHEM_TAGGED_OX

            !---------------------------------
            ! Tagged CO
            !---------------------------------
            CASE ( 7  )
               CALL CHEM_TAGGED_CO

            !---------------------------------
            ! C2H6
            !---------------------------------
            CASE ( 8  )
               CALL CHEMC2H6

            !---------------------------------
            ! CH4
            !---------------------------------
            CASE ( 9  )
               ! Only call after the first 24 hours
               IF ( GET_ELAPSED_MIN() >= GET_TS_CHEM() ) THEN
                  CALL CHEMCH4
               ENDIF

            !---------------------------------
            ! Offline sulfate simulation
            !---------------------------------
            CASE ( 10 )

               IF ( LSULF ) THEN

                  ! Get relative humidity
                  CALL MAKE_RH

                  ! Do sulfate chemistry
                  CALL CHEMSULFATE

                  ! Do aerosol phase equilibrium
                  CALL DO_RPMARES
               ENDIF

               ! Do carbonaceous aerosol chemistry
               IF ( LCARB ) CALL CHEMCARBON

               ! Do dust aerosol chemistry
               IF ( LDUST ) CALL CHEMDUST

               ! Do seasalt aerosol chemistry
               IF ( LSSALT ) CALL CHEMSEASALT

!-----------------------------------------------------------------------------
! Prior 
!#if   defined( LGEOSCO )
!                  ! Parameterized OH
!                  ! Call after 24 hours have passed
!                  CASE ( 5 )
!                     IF ( NMIN >= NCHEM ) THEN 
!                        CALL CHEMCO( FIRSTCHEM, LMN, 
!     &                       NSEASON,   NYMDb, NYMDe )
!                     ENDIF 
!#endif
!-----------------------------------------------------------------------------

            !---------------------------------
            ! Kr85   
            !---------------------------------
            CASE ( 12 )
               CALL CHEMKr85

         END SELECT
                  
         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CHEMISTRY' )
      ENDIF

      !-----------------------------------------------------------------------
      ! Prior to 4/22/04:
      ! Now call this from "main.f", after wetdep, since we need to save
      ! out soluble tracers for ICARTT (bmy, 4/22/04)
      !!=================================================================
      !! Call the planeflight diagnostic
      !!=================================================================
      !IF ( ND40 > 0 ) THEN
      !   CALL PLANEFLIGHT
      !
      !   !### Debug
      !   IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a PLANEFLIGHT' )
      !ENDIF
      !-----------------------------------------------------------------------
         
      ! Return to calling program
      END SUBROUTINE DO_CHEMISTRY

!------------------------------------------------------------------------------

      END MODULE CHEMISTRY_MOD
