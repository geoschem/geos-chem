!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
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
      SUBROUTINE INPHOT( NLAYER, NREACS, Input_Opt, am_I_Root )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : ERROR_STOP
      USE CMN_FJ_MOD,         ONLY : JPMAX, JPPJ, JPNL
      USE CMN_SIZE_MOD,       ONLY : LLPAR
      USE JV_CMN_MOD
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE inquireMod,         ONLY : findFreeLUN

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER,        INTENT(IN) :: NLAYER    ! # of layers for FAST-J photolysis
      INTEGER,        INTENT(IN) :: NREACS    ! # of FAST-J photolysis reactions
      LOGICAL,        INTENT(IN) :: am_I_Root ! Is this the root CPU?
      TYPE(OptInput), INTENT(IN) :: Input_Opt ! Input Options
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
!  01 Aug 2012 - R. Yantosca - Add reference to findFreeLUN from inqure_mod.F90
!  03 Aug 2012 - R. Yantosca - Move calls to findFreeLUN out of DEVEL block
!  07 Aug 2012 - R. Yantosca - Now print LUN used to open file
!  26 Feb 2013 - M. Long     - Now accept Input_Opt as an argument
!  19 Mar 2013 - R. Yantosca - When using ESMF interface to GEOS-5, append
!                              ".rc" to filenames (instead of ___.rc)
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      INTEGER            :: LUN
      CHARACTER(LEN=255) :: FILENAME

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

      !-----------------------------------------------------------------
      ! Read in the names of FAST-J photolysis species from "ratj.d"
      !-----------------------------------------------------------------

      ! Find a free file unit
      LUN      = findFreeLUN()

      ! Name of file to read
      FILENAME = 'ratj.d'

#if defined( ESMF_ )
      ! We need to add the ".rc" extension to the "ratj.d" file
      ! so that GEOS-5 will copy it to the proper run directory.
      FILENAME = TRIM( FILENAME ) // '.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), LUN
 100     FORMAT( 'GEOSCHEM::INPHOT: Reading ', a ', on unit ', i4 )
      ENDIF

      ! Read J_value names 
      CALL RD_JS( LUN, TRIM( FILENAME ), am_I_Root )

      !-----------------------------------------------------------------
      ! Call JV_INDEX to translate between GEOS-Chem species 
      ! nomenclature and Fast-J species nomenclature (bmy, 9/13/99)
      !-----------------------------------------------------------------
      CALL JV_INDEX( am_I_Root )

      !-----------------------------------------------------------------
      ! Read in JPL spectral data set (e.g. X-sections, quantum yields)
      !-----------------------------------------------------------------

      ! Find a free file unit
      LUN      = findFreeLUN()

      ! Name of file to read
      FILENAME = 'jv_spec.dat'

#if defined( ESMF_ )
      ! We need to add the ".rc" extension to the "jv_spec.dat" file
      ! so that GEOS-5 will copy it to the proper run directory.
      FILENAME = TRIM( FILENAME ) // '.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), LUN
      ENDIF
      
      ! Read JPL data from disk
      CALL RD_TJPL( LUN, TRIM( FILENAME ), am_I_Root )

      !-----------------------------------------------------------------
      ! Read in aerosol optics to be used to calculate AODs (clh)
      !-----------------------------------------------------------------

      ! Find a free file unit
      LUN      = findFreeLUN()

      ! Name of file to read
      FILENAME = 'jv_spec_aod.dat'

#if defined( ESMF_ )
      ! We need to add the ".rc" extension to the "jv_spec_aod.dat"
      ! file so that GEOS-5 will copy it to the proper run directory.
      FILENAME = TRIM( FILENAME ) // '.rc'
#endif

      ! Echo info
      IF ( am_I_Root ) THEN
         WRITE( 6, 100 ) TRIM( FILENAME ), LUN
      ENDIF

      ! Read aerosol optics from disk
      CALL RD_AOD( LUN, TRIM( FILENAME ), am_I_Root )

      !-----------------------------------------------------------------
      ! Read in T & O3 climatology (cf. Nagatani/92 and McPeters/91)
      !-----------------------------------------------------------------
      CALL READ_JV_ATMS_DAT( am_I_Root, Input_Opt%Data_Dir_1x1 )

      !-----------------------------------------------------------------
      ! Select Aerosol/Cloud types to be used
      !-----------------------------------------------------------------
      CALL SET_AER( am_I_Root )

      END SUBROUTINE INPHOT
!EOC
