!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: inphot
!
! !DESCRIPTION: Subroutine INPHOT initializes quantities for FAST-J photolysis,
!  including JPL spectral data (e.g. cross sections, quantum yields), standard 
!  O3 and T profiles, and the translation indices between GEOS-Chem and FAST-J
!  species names.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INPHOT( NLAYER, NREACS, am_I_Root )
!
! !USES:
!
#     include "define.h"

      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE FILE_MOD,     ONLY : IU_FASTJ
      USE CMN_FJ_MOD,   ONLY : JPMAX, JPPJ, JPNL
      USE CMN_SIZE_MOD, ONLY : LLPAR
      USE JV_CMN_MOD

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN) :: NLAYER      ! # of layers for FAST-J photolysis
      INTEGER, INTENT(IN) :: NREACS      ! # of FAST-J photolysis reactions
      LOGICAL, INTENT(IN) :: am_I_Root   ! Is this the root CPU?
! 
! !REVISION HISTORY: 
!  01 Apr 1999 - O. Wild     - Initial version
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
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  30 Jul 2012 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now add LUN as a local variable
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      INTEGER :: LUN

      !=================================================================
      ! INPHOT begins here!
      !=================================================================

#if   !defined( DEVEL )
      ! For the traditional driver "main.F", we use IU_FASTJ as the
      ! LUN for ASCII file I/O.  (bmy, 7/30/12)
      LUN   = IU_FASTJ
#endif

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
      CALL RD_JS( LUN, 'ratj.d', am_I_Root )

      ! Call JV_INDEX to translate between GEOS-Chem species 
      ! nomenclature and Fast-J species nomenclature (bmy, 9/13/99)
      CALL JV_INDEX( am_I_Root )

      ! Read in JPL spectral data set (e.g. X-sections, quantum yields)
      CALL RD_TJPL( LUN, 'jv_spec.dat___.rc', am_I_Root )

      ! Read in aerosol optics to be used to calculate AODs (clh)
      CALL RD_AOD( LUN, 'jv_spec_aod.dat___.rc', am_I_Root )

      ! Read in T & O3 climatology (cf. Nagatani/92 and McPeters/91)
      CALL READ_JV_ATMS_DAT( 'fastj.jv_atms_dat.nc', am_I_Root )

      ! Select Aerosol/Cloud types to be used
      CALL SET_AER( am_I_Root )

      ! Return to calling program
      END SUBROUTINE INPHOT
