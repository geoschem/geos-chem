! $Id: inphot.f,v 1.1 2009/09/16 14:06:24 bmy Exp $
      SUBROUTINE INPHOT( NLAYER, NREACS )
!
!******************************************************************************
!  Subroutine INPHOT initializes quantities for FAST-J photolysis, including
!  JPL spectral data (e.g. cross sections, quantum yields), standard O3 and T 
!  profiles, and the translation indices between GEOS-Chem and FAST-J species
!  names. (Oliver Wild, 4/99, ppm, bmy, 9/7/99, 5/10/10)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) NLAYER (INTEGER) : Number of layers for FAST-J photolysis
!  (2 ) NREACS (INTEGER) : Total # of photolysis reactions for FAST-J
!
!  NOTES: 
!  (1 ) Remove PTOP from the arg list, since it is now a 
!        parameter in "CMN_SIZE" (bmy, 2/10/00).
!  (2 ) Remove SIGE from the argument list, since we are now using
!        a hybrid pressure specification.  Now define ETAA and ETAB
!        for use in "set_prof.f". (bmy, 8/23/02)
!  (3 ) Now reference ERROR_STOP from "error_mod.f".  Updated comments and
!        made cosmetic changes (bmy, 10/15/02)
!  (4 ) Remove IPH -- now use IU_FASTJ directly (bmy, 4/8/03)
!  (5 ) Removed ETAA and ETAB arrays.  We now compute PJ directly from the 
!        GET_PEDGE routine.  Also remove reference to "pressure_mod.f".  
!        Updated comments. (bmy, 10/30/07)
!  (6 ) Read jv_spec_aod.dat file for AOD diagnostics. (clh, bmy, 5/10/10)
!******************************************************************************
!
      ! References to F90 modules (bmy, 6/27/02)
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE FILE_MOD,     ONLY : IU_FASTJ

#     include "define.h"

      USE CMN_FJ_MOD,   ONLY : JPMAX, JPPJ, JPNL
      USE CMN_SIZE_MOD, ONLY : LLPAR
      USE JV_CMN_MOD

      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN) :: NLAYER, NREACS
      
      !=================================================================
      ! INPHOT begins here!
      !=================================================================

      ! # of layers to do chemistry
      JPNL  = NLAYER             

      ! # of reactions in chemistry
      JPPJ  = NREACS + 4         

      ! Error check # of layers
      IF ( JPNL > LLPAR ) THEN 
         CALL ERROR_STOP( 'JPNL > LLPAR!', 'inphot.f' )
      ENDIF

      ! Error check # of rxns
      IF ( JPPJ > JPMAX ) THEN
         CALL ERROR_STOP( 'JPPJ > JPMAX!', 'inphot.f' )
      ENDIF

      ! Read in labels of photolysis rates required
      CALL RD_JS( IU_FASTJ, 'ratj.d' )

      ! Call JV_INDEX to translate between GEOS-Chem species 
      ! nomenclature and Fast-J species nomenclature (bmy, 9/13/99)
      CALL JV_INDEX 

      ! Read in JPL spectral data set (e.g. X-sections, quantum yields)
      CALL RD_TJPL( IU_FASTJ, 'jv_spec.dat' )

      ! Read in aerosol optics to be used to calculate AODs (clh)
      CALL RD_AOD( IU_FASTJ, 'jv_spec_aod.dat' )

      ! Read in T & O3 climatology (cf. Nagatani/92 and McPeters/91)
      CALL RD_PROF( IU_FASTJ, 'jv_atms.dat' )
 
      ! Select Aerosol/Cloud types to be used
      CALL SET_AER

      ! Return to calling program
      END SUBROUTINE INPHOT
