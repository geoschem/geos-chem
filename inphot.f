! $Id: inphot.f,v 1.1 2003/06/30 20:26:05 bmy Exp $
      SUBROUTINE INPHOT( NLAYER, NREACS )
!
!******************************************************************************
!  Routine to initialise photolysis rate data, called directly from the
!  cinit routine in ASAD. Currently use it to read the JPL spectral data
!  and standard O3 and T profiles and to set the appropriate reaction index.
!-----------------------------------------------------------------------
!  Add the following input variables (bmy, 9/7/99)
!
!  Variable  Type    Dimensn Units   Description
!  --------  ----    ------- -----   -----------
!  nlayer    int        -      -     Top layer where J-values required
!  nreacs    int        -      -     Total # of photolysis reactions
!                                    starting with bottom of lowest layer
!
!  NOTES: 
!  (1) Remove PTOP from the arg list, since it is now a 
!       parameter in "CMN_SIZE" (bmy, 2/10/00).
!  (2) Remove SIGE from the argument list, since we are now using
!       a hybrid pressure specification.  Now define ETAA and ETAB
!       for use in "set_prof.f". (bmy, 8/23/02)
!  (3) Now reference ERROR_STOP from "error_mod.f".  Updated comments and
!       made cosmetic changes (bmy, 10/15/02)
!  (4) Remove IPH -- now use IU_FASTJ directly (bmy, 4/8/03)
!------------------------------------------------------------------------------
!  Other relevant variables:
!
!     rad       Radius of Earth (cm)
!     zzht      Effective scale height above top of atmosphere (cm)
!     dtaumax   Maximum opt.depth above which a new level should be inserted
!     dtausub   No. of opt.depths at top of cloud requiring subdivision
!     dsubdiv   Number of additional levels to add at top of cloud
!     szamax    Solar zenith angle cut-off, above which to skip calculation
!
!******************************************************************************
!
      ! References to F90 modules (bmy, 6/27/02)
      USE ERROR_MOD,    ONLY : ERROR_STOP
      USE FILE_MOD,     ONLY : IU_FASTJ
      USE PRESSURE_MOD, ONLY : GET_AP, GET_BP

      IMPLICIT NONE

#     include "cmn_fj.h"
#     include "jv_cmn.h"

      ! Arguments
      INTEGER, INTENT(IN) :: NLAYER, NREACS
      
      ! Local variables
      !-----------------------------------------------
      ! Prior to 4/8/03:
      ! Remove IPH (bmy, 4/8/03)
      !INTEGER             :: IPH, L
      !-----------------------------------------------
      INTEGER             :: L

      !=================================================================
      ! INPHOT begins here!
      !=================================================================

      ! Define A and B coordinates for hybrid grid (bmy, 8/22/02)
      DO L = 1, LPAR + 1
         ETAA(L) = GET_AP(L)
         ETAB(L) = GET_BP(L)
      ENDDO

      JPNL  = NLAYER             ! # of layers to do chemistry
      JPPJ  = NREACS + 4         ! # of reactions in chemistry

      IF ( JPNL > LPAR ) THEN 
         CALL ERROR_STOP( 'JPNL > LPAR!', 'inphot.f' )
      ENDIF

      IF ( JPPJ > JPMAX ) THEN
         CALL ERROR_STOP( 'JPPJ > JPMAX!', 'inphot.f' )
      ENDIF

      !---------------------------------------------------------------------
      ! Prior to 4/8/03:
      ! Now use IU_FASTJ directly (bmy, 4/8/03)
      !! Use channel 8 to read files at the moment
      !! Use file unit IU_FASTJ (which =8) to read files (bmy, 6/27/02)
      !IPH = IU_FASTJ
      !---------------------------------------------------------------------

      ! Read in labels of photolysis rates required
      CALL RD_JS( IU_FASTJ, 'ratj.d' )

      ! Call JV_INDEX to translate between GEOS-CHEM species 
      ! nomenclature and Fast-J species nomenclature (bmy, 9/13/99)
      CALL JV_INDEX 

      ! Read in JPL spectral data set
      CALL RD_TJPL( IU_FASTJ, 'jv_spec.dat' )

      ! Read in T & O3 climatology
      CALL RD_PROF( IU_FASTJ, 'jv_atms.dat' )
 
      ! Select Aerosol/Cloud types to be used
      CALL SET_AER

      ! Return to calling program
      END SUBROUTINE INPHOT
