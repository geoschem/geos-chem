! $Id: chemistry_mod.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      MODULE CHEMISTRY_MOD
!
!******************************************************************************
!  Module CHEMISTRY_MOD is used to call the proper chemistry subroutine
!  for the various GEOS-CHEM simulations. (bmy, 4/14/03)
! 
!  Module Routines:
!  ============================================================================
!  (1 ) DO_CHEMISTRY       : Driver which calls various chemistry routines
!
!  GEOS-CHEM modules referenced by tpcore_call_mod.f
!  ============================================================================
!  (1 ) c2h6_mod.f         : Module containing routines for C2H6 chemistry
!  (2 ) ch3i_mod.f         : Module containing routines for CH3I chemistry
!  (3 ) global_ch4_mod.f   : Module containing routines for CH4 chemistry
!  (4 ) RnPbBe_mod.f       : Module containing routines for Rn-Pb-Be chemistry
!  (5 ) tagged_co_mod.f    : Module containing routines for Tagged CO chemistry
!
!  NOTES:
!  (1 ) Bug fix in DO_CHEMISTRY (bnd, bmy, 4/14/03)
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
!  chemistry subroutine for the various GEOS-CHEM simulations. (bmy, 2/11/03)
!
!  NOTES:
!  (1 ) Now reference DELP, T from "dao_mod.f" since we need to pass this
!        to OPTDEPTH for GEOS-1 or GEOS-STRAT met fields (bnd, bmy, 4/14/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ACETONE_MOD,     ONLY : OCEAN_SINK_ACET
      USE C2H6_MOD,        ONLY : CHEMC2H6
      USE CH3I_MOD,        ONLY : CHEMCH3I
      !-----------------------------------------------------------------------
      ! Prior to 4/14/03 -- supplemental fix
      ! Need to pass DELP, T via "dao_mod.f" (bnd, bmy, 4/14/03)
      !USE DAO_MOD,         ONLY : CLMOSW, CLROSW, MAKE_RH, OPTDEP, OPTD
      !-----------------------------------------------------------------------
      USE DAO_MOD,         ONLY : CLMOSW,  CLROSW, DELP, 
     &                            MAKE_RH, OPTDEP, OPTD, T
      USE GLOBAL_CH4_MOD,  ONLY : CHEMCH4
      USE DRYDEP_MOD,      ONLY : DRYFLX, DRYFLXRnPbBe
      USE OPTDEPTH_MOD,    ONLY : OPTDEPTH
      USE PLANEFLIGHT_MOD, ONLY : PLANEFLIGHT
      USE RnPbBe_MOD,      ONLY : CHEMRnPbBe
      USE RPMARES_MOD,     ONLY : DO_RPMARES
      USE TAGGED_CO_MOD,   ONLY : CHEM_TAGGED_CO
      USE TIME_MOD,        ONLY : GET_ELAPSED_MIN, GET_TS_CHEM
      USE TRACERID_MOD,    ONLY : IDTACET
      USE SULFATE_MOD,     ONLY : CHEMSULFATE

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

#elif defined( GEOS_2 ) || defined( GEOS_3 ) 

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
            CASE ( 3  )

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
               CALL CHEMO3

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

               ! Get relative humidity
               CALL MAKE_RH

               ! Do sulfate chemistry
               CALL CHEMSULFATE

               ! Do aerosol phase equilibrium
               CALL DO_RPMARES

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

         END SELECT
                  
         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CHEMISTRY' )
      ENDIF

      !=================================================================
      ! Call the planeflight diagnostic
      !=================================================================
      IF ( ND40 > 0 ) THEN
         CALL PLANEFLIGHT

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a PLANEFLIGHT' )
      ENDIF
         
      ! Return to calling program
      END SUBROUTINE DO_CHEMISTRY

!------------------------------------------------------------------------------

      END MODULE CHEMISTRY_MOD
