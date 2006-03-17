! $Id: chemistry_mod.f,v 1.20 2006/03/17 15:30:36 bmy Exp $
      MODULE CHEMISTRY_MOD
!
!******************************************************************************
!  Module CHEMISTRY_MOD is used to call the proper chemistry subroutine
!  for the various GEOS-CHEM simulations. (bmy, 4/14/03, 3/16/06)
! 
!  Module Routines:
!  ============================================================================
!  (1 ) DO_CHEMISTRY       : Driver which calls various chemistry routines
!
!  GEOS-CHEM modules referenced by chemistry_mod.f
!  ============================================================================
!  (1 ) acetone_mod.f      : Module w/ routines for ACET chemistry
!  (2 ) c2h6_mod.f         : Module w/ routines for C2H6 chemistry
!  (3 ) carbon_mod.f       : Module w/ routines for carbon arsl chem.
!  (4 ) ch3i_mod.f         : Module w/ routines for CH3I chemistry
!  (5 ) dao_mod.f          : Module w/ arrays for DAO met fields
!  (6 ) diag_pl_mod.f      : Module w/ routines to save P(Ox), L(Ox)
!  (7 ) drydep_mod.f       : Module w/ GEOS-CHEM drydep routines
!  (8 ) dust_mod.f         : Module w/ routines for dust arsl chem.
!  (9 ) error_mod.f        : Module w/ NaN and error checks
!  (10) global_ch4_mod.f   : Module w/ routines for CH4 chemistry
!  (11) hcn_ch3cn_mod.f    : Module w/ routines for HCN and CH3CN chemistry
!  (12) Kr85_mod.f         : Module w/ routines for Kr85 chemistry
!  (13) logical_mod.f      : Module w/ GEOS-CHEM logical switches
!  (14) RnPbBe_mod.f       : Module w/ routines for Rn-Pb-Be chemistry
!  (15) rpmares_mod.f      : Module w/ routines for arsl phase equilib.
!  (16) seasalt_mod.f      : Module w/ routines for seasalt chemistry
!  (17) sulfate_mod.f      : Module w/ routines for sulfate chemistry
!  (18) tagged_co_mod.f    : Module w/ routines for Tagged CO chemistry
!  (19) tagged_ox_mod.f    : Module w/ routines for Tagged Ox chemistry
!  (20) time_mod.f         : Module w/ routines to compute time & date
!  (21) tracer_mod.f       : Module w/ GEOS-CHEM tracer array STT etc. 
!  (22) tracerid_mod.f     : Module w/ pointers to tracers & emissions
!
!  NOTES:
!  (1 ) Bug fix in DO_CHEMISTRY (bnd, bmy, 4/14/03)
!  (2 ) Now references DEBUG_MSG from "error_mod.f" (bmy, 8/7/03)
!  (3 ) Now references "tagged_ox_mod.f"(bmy, 8/18/03)
!  (4 ) Now references "Kr85_mod.f" (jsw, bmy, 8/20/03)
!  (5 ) Bug fix: Now also call OPTDEPTH for GEOS-4 (bmy, 1/27/04)
!  (6 ) Now references "carbon_mod.f" and "dust_mod.f" (rjp, tdf, bmy, 4/5/04)
!  (7 ) Now references "seasalt_mod.f" (rjp, bec, bmy, 4/20/04)
!  (8 ) Now references "logical_mod.f", "tracer_mod.f", "diag20_mod.f", and
!        "diag65_mod.f", and "aerosol_mod." (bmy, 7/20/04)
!  (9 ) Now references "mercury_mod.f" (bmy, 12/7/04)
!  (10) Updated for SO4s, NITs chemistry (bec, bmy, 4/13/05)
!  (11) Now call CHEM_HCN_CH3CN from "hcn_ch3cn_mod.f".  Also remove all
!        references to the obsolete CO-OH param simulation. (xyp, bmy, 6/24/05)
!  (12) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (13) Now we call MAKE_RH from "main.f" (bmy, 3/16/06)
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
!  (bmy, 2/11/03, 10/3/05)
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
!  (6 ) Now calls CHEMSEASALT to do seasalt aerosol chemistry 
!        (rjp, bec, bmy, 4/20/04)
!  (7 ) Now references "logical_mod.f" & "tracer_mod.f".  Now references
!        AEROSOL_CONC, AEROSOL_RURALBOX, and RDAER from "aerosol_mod.f".  
!        Now includes "CMN_DIAG" and "comode.h".  Also call READER, READCHEM, 
!        and INPHOT to initialize the FAST-J arrays so that we can save out !
!        AOD's to the ND21 diagnostic for offline runs. (bmy, 7/20/04)
!  (8 ) Now call routine CHEMMERCURY from "mercury_mod.f" for an offline
!        Hg0/Hg2/HgP simulation. (eck, bmy, 12/7/04)
!  (9 ) Now do not call DO_RPMARES if we are doing an offline aerosol run
!        with crystalline sulfur & aqueous tracers (cas, bmy, 1/7/05)
!  (10) Now use ISOROPIA for aer thermodyn equilibrium if we have seasalt 
!        tracers defined, or RPMARES if not.  Now call CHEMSEASALT before
!        CHEMSULFATE.  Now do aerosol thermodynamic equilibrium before
!        aerosol chemistry for offline aerosol runs.  Now also reference 
!        CLDF from "dao_mod.f" (bec, bmy, 4/20/05)
!  (11) Now modified for GCAP met fields.  Now call CHEM_HCN_CH3CN from 
!        "hcn_ch3cn_mod.f".  Also remove allreferences to the obsolete 
!         CO-OH param simulation. (xyp, bmy, 6/23/05)
!  (12) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE ACETONE_MOD,     ONLY : OCEAN_SINK_ACET
      USE AEROSOL_MOD,     ONLY : AEROSOL_CONC, AEROSOL_RURALBOX
      USE AEROSOL_MOD,     ONLY : RDAER,        SOILDUST
      USE C2H6_MOD,        ONLY : CHEMC2H6
      USE CARBON_MOD,      ONLY : CHEMCARBON
      USE CH3I_MOD,        ONLY : CHEMCH3I
      USE DAO_MOD,         ONLY : CLDF,    CLMOSW, CLROSW, DELP
      !-------------------------------------------------------------
      ! Prior to 3/16/06:
      !USE DAO_MOD,         ONLY : MAKE_RH, OPTDEP, OPTD,   T
      !-------------------------------------------------------------
      USE DAO_MOD,         ONLY : OPTDEP, OPTD,   T
      USE DRYDEP_MOD,      ONLY : DRYFLX, DRYFLXRnPbBe
      USE DUST_MOD,        ONLY : CHEMDUST, RDUST_ONLINE
      USE ERROR_MOD,       ONLY : DEBUG_MSG
      USE GLOBAL_CH4_MOD,  ONLY : CHEMCH4
      USE HCN_CH3CN_MOD,   ONLY : CHEM_HCN_CH3CN
      USE ISOROPIA_MOD,    ONLY : DO_ISOROPIA
      USE Kr85_MOD,        ONLY : CHEMKr85
      USE LOGICAL_MOD,     ONLY : LCARB, LCHEM,  LCRYST, LDUST
      USE LOGICAL_MOD,     ONLY : LPRT,  LSSALT, LSULF  
      USE MERCURY_MOD,     ONLY : CHEMMERCURY
      USE OPTDEPTH_MOD,    ONLY : OPTDEPTH
      USE RnPbBe_MOD,      ONLY : CHEMRnPbBe
      USE RPMARES_MOD,     ONLY : DO_RPMARES
      USE SEASALT_MOD,     ONLY : CHEMSEASALT
      USE SULFATE_MOD,     ONLY : CHEMSULFATE
      USE TAGGED_CO_MOD,   ONLY : CHEM_TAGGED_CO
      USE TAGGED_OX_MOD,   ONLY : CHEM_TAGGED_OX
      USE TIME_MOD,        ONLY : GET_ELAPSED_MIN, GET_TS_CHEM
      USE TRACER_MOD,      ONLY : N_TRACERS,       STT  
      USE TRACER_MOD,      ONLY : ITS_A_C2H6_SIM
      USE TRACER_MOD,      ONLY : ITS_A_CH3I_SIM
      USE TRACER_MOD,      ONLY : ITS_A_CH4_SIM
      USE TRACER_MOD,      ONLY : ITS_A_FULLCHEM_SIM
      USE TRACER_MOD,      ONLY : ITS_A_HCN_SIM
      USE TRACER_MOD,      ONLY : ITS_A_MERCURY_SIM
      USE TRACER_MOD,      ONLY : ITS_A_RnPbBe_SIM
      USE TRACER_MOD,      ONLY : ITS_A_TAGCO_SIM
      USE TRACER_MOD,      ONLY : ITS_A_TAGOX_SIM
      USE TRACER_MOD,      ONLY : ITS_AN_AEROSOL_SIM
      USE TRACER_MOD,      ONLY : ITS_NOT_COPARAM_OR_CH4
      USE TRACERID_MOD,    ONLY : IDTACET

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_DIAG"  ! NDxx flags
#     include "comode.h"  ! NPHOT

      ! Local variables
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER       :: N_TROP

      !=================================================================
      ! DO_CHEMISTRY begins here!
      !=================================================================

#if   defined( GEOS_1 ) || defined( GEOS_STRAT )

      ! Compute optical depths (except for CO-OH run)
      ! GEOS-1/GEOS-STRAT: compute from T, CLMO, CLRO, DELP
      IF ( ITS_NOT_COPARAM_OR_CH4() ) THEN
         CALL OPTDEPTH( LLPAR, CLMOSW, CLROSW, DELP, T, OPTD )
      ENDIF

#elif defined( GEOS_3 ) || defined( GEOS_4 ) || defined( GCAP )

      ! Compute optical depths (except for CO-OH run)
      ! GEOS-2/GEOS-3: Copy OPTDEP to OPTD, also archive diagnostics
      IF ( ITS_NOT_COPARAM_OR_CH4() ) THEN
         CALL OPTDEPTH( LLPAR, CLDF, OPTDEP, OPTD )
      ENDIF
#endif 


      !=================================================================
      ! If LCHEM=T then call the chemistry subroutines
      !=================================================================
      IF ( LCHEM ) THEN 

         !---------------------------------
         ! NOx-Ox-HC (w/ or w/o aerosols) 
         !---------------------------------
         IF ( ITS_A_FULLCHEM_SIM() ) THEN 

            ! Call SMVGEAR routines
            CALL CHEMDR

            ! Do seasalt aerosol chemistry
            IF ( LSSALT ) CALL CHEMSEASALT

            ! Also do sulfate chemistry
            IF ( LSULF ) THEN
               
               !--------------------------------------------
               ! Prior to 3/16/06:
               ! Now call this from main.f (bmy, 3/16/06)
               !! Get relative humidity
               !CALL MAKE_RH
               !--------------------------------------------

               ! Do sulfate chemistry
               CALL CHEMSULFATE

               ! Do aerosol thermodynamic equilibrium
               IF ( LSSALT ) THEN

                  ! ISOROPIA takes Na+, Cl- into account
                  CALL DO_ISOROPIA

               ELSE

                  ! RPMARES does not take Na+, Cl- into account
                  CALL DO_RPMARES

               ENDIF
               
            ENDIF

             ! Do carbonaceous aerosol chemistry
            IF ( LCARB ) CALL CHEMCARBON

            ! Do dust aerosol chemistry
            IF ( LDUST ) CALL CHEMDUST

            ! ND44 drydep fluxes
            CALL DRYFLX     

            ! ND43 chemical production
            CALL DIAGOH

            ! Remove acetone ocean sink
            IF ( IDTACET /= 0 ) THEN
               CALL OCEAN_SINK_ACET( STT(:,:,1,IDTACET) ) 
            ENDIF

         !---------------------------------
         ! Offline aerosol simulation
         !---------------------------------
         ELSE IF ( ITS_AN_AEROSOL_SIM() ) THEN

            !--------------------------------------------
            ! Prior to 3/16/06:
            ! Now call this from main.f (bmy, 3/16/06)
            !! Get relative humidity
            !CALL MAKE_RH
            !--------------------------------------------

            ! Define loop index and other SMVGEAR arrays
            ! N_TROP, the # of trop boxes, is returned
            CALL AEROSOL_RURALBOX( N_TROP )

            ! Initialize FAST-J quantities for computing AOD's
            IF ( FIRST ) THEN
               CALL READER( FIRST )
               CALL READCHEM
               CALL INPHOT( LLTROP, NPHOT )

               ! Reset NCS with NCSURBAN
               NCS     = NCSURBAN

               ! Reset NTLOOP and NTTLOOP after call to READER
               ! with the actual # of boxes w/in the ann mean trop
               NTLOOP  = N_TROP
               NTTLOOP = N_TROP

               ! Reset first-time flag
               FIRST = .FALSE.
            ENDIF

            ! Compute aerosol & dust concentrations [kg/m3]
            ! (NOTE: SOILDUST in "aerosol_mod.f" is computed here)
            CALL AEROSOL_CONC

            ! Compute AOD's and surface areas
            CALL RDAER

            !*** AEROSOL THERMODYNAMIC EQUILIBRIUM ***
            IF ( LSSALT ) THEN

               ! ISOROPIA takes Na+, Cl- into account
               CALL DO_ISOROPIA

            ELSE

               ! RPMARES does not take Na+, Cl- into account
               ! (skip for crystalline & aqueous offline run)
               IF ( .not. LCRYST ) CALL DO_RPMARES

            ENDIF

            !*** SEASALT AEROSOLS ***
            IF ( LSSALT ) CALL CHEMSEASALT

            !*** SULFATE AEROSOLS ***
            IF ( LSULF .or. LCRYST ) THEN

               ! Do sulfate chemistry
               CALL CHEMSULFATE

            ENDIF
               
            !*** CARBON AND 2NDARY ORGANIC AEROSOLS ***
            IF ( LCARB ) CALL CHEMCARBON

            !*** MINERAL DUST AEROSOLS ***
            IF ( LDUST ) THEN 

               ! Do dust aerosol chemsitry
               CALL CHEMDUST

               ! Compute dust OD's & surface areas
               CALL RDUST_ONLINE( SOILDUST )
            ENDIF

         !---------------------------------
         ! Rn-Pb-Be
         !---------------------------------                 
         ELSE IF ( ITS_A_RnPbBe_SIM() ) THEN

            CALL CHEMRnPbBe 
            CALL DRYFLXRnPbBe
                  
         !---------------------------------
         ! CH3I
         !---------------------------------
         ELSE IF ( ITS_A_CH3I_SIM() ) THEN
            CALL CHEMCH3I

         !---------------------------------            
         ! HCN
         !---------------------------------
         ELSE IF ( ITS_A_HCN_SIM() ) THEN
            CALL CHEM_HCN_CH3CN( N_TRACERS, STT )

         !---------------------------------
         ! Tagged O3
         !---------------------------------
         ELSE IF ( ITS_A_TAGOX_SIM() ) THEN 
            CALL CHEM_TAGGED_OX

         !---------------------------------
         ! Tagged CO
         !---------------------------------
         ELSE IF ( ITS_A_TAGCO_SIM() ) THEN
            CALL CHEM_TAGGED_CO

         !---------------------------------
         ! C2H6
         !---------------------------------
         ELSE IF ( ITS_A_C2H6_SIM() ) THEN
            CALL CHEMC2H6

         !---------------------------------
         ! CH4
         !---------------------------------
         ELSE IF ( ITS_A_CH4_SIM() ) THEN

            ! Only call after the first 24 hours
            IF ( GET_ELAPSED_MIN() >= GET_TS_CHEM() ) THEN
               CALL CHEMCH4
            ENDIF

         !---------------------------------
         ! Mercury
         !---------------------------------
         ELSE IF ( ITS_A_MERCURY_SIM() ) THEN

            ! Get relative humidity
            CALL MAKE_RH    

            ! Do Hg chemistry
            CALL CHEMMERCURY
               
!-----------------------------------------------------------------------------
! Prior to 7/19/04:
! Fully install Kr85 run later (bmy, 7/19/04)
!            !---------------------------------
!            ! Kr85   
!            !---------------------------------
!            CASE ( 12 )
!               CALL CHEMKr85
!-----------------------------------------------------------------------------
         ENDIF

         !### Debug
         IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a CHEMISTRY' )
      ENDIF
         
      ! Return to calling program
      END SUBROUTINE DO_CHEMISTRY

!------------------------------------------------------------------------------

      END MODULE CHEMISTRY_MOD
