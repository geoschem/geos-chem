!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: cmn_fj_mod
!
! !DESCRIPTION: Module cmn\_fj\_mod contains parameters and global variables
!  used to interface between Harvard chemistry and UC-Irvine Fast-J photolysis 
!  programs.
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_FJ_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR, LLPAR

      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS:
!
      ! max # of photolysis rxns = 4 + IPHOT (see comode.h)
      ! FP increased JPMAX since IPHOT was increased (hotp 7/31/09)
      INTEGER, PARAMETER :: JPMAX = 89
!
! !PUBLIC DATA MEMBERS:
!
      ! Variables for number of layers and number of photolysis rxns
      INTEGER            :: JPNL, JPPJ       

      ! Branches for photolysis species
      INTEGER            :: BRANCH(JPMAX)

      ! Names of photolysis species
      ! FP increased length of RNAMES for species indistinguishable
      ! with only 4 characters (hotp 7/31/09)
      ! used in jv_index and rd_js.f
      !CHARACTER (LEN=4)  :: RNAMES
      CHARACTER (LEN=7)  :: RNAMES(JPMAX)

      ! Mapping array from Harvard species names to UCI species names
      INTEGER            :: RINDEX(JPMAX)

      ! Output J-values
      REAL*8, ALLOCATABLE :: ZPJ(:,:,:,:)
! 
! !REMARKS:
!  Based on code from Oliver Wild (9 Jul 1999)
!
! !REVISION HISTORY:
!  (1 ) Uses Fortran 90 declarations for parameters and variables
!  (2 ) Pass CTM size parameters and preprocessor switches via CMN_SIZE.
!  (3 ) Update JPMAX for new chemistry mechanism (amf, bmy, 4/20/00)
!  (4 ) Return JPMAX to original setting (bmy, 9/25/00)
!  (5 ) Return JPMAX to 55 for peroxy recycling (again) (bmy, 12/20/00)
!  (6 ) Now need to use the window parameters IIPAR,JJPAR,LLPAR (bmy, 9/25/01)
!  (7 ) Changed RCS ID tag comment character from "C" to "!" to allow freeform
!        compilation. (bmy, 6/25/02)
!  (8 ) Replaced ESIG array with ETAA and ETAB arrays for the hybrid
!        pressure formulation.  Also deleted PREST, since we don't need that
!        anymore. (bmy, 8/23/02)
!  (9 ) Removed ETAA and ETAB arrays.  We now compute PJ directly from the 
!        GET_PEDGE routine. (bmy, 10/30/07)
!  (10) Increase photolysis rxns JPMAX = 79 (tmf, 1/7/09)
!  (11) Increase photolysis rxns JPMAX = 89 for Isoprene (fp, 2/2/10)
!  (12) Increase species name length. (fp, 2/2/10)
!  23 Aug 2011 - M. Long   - Converted to Module from Header file
!  10 Aug 2012 - R. Yantosca - Replace IPAR, JPAR, LPAR w/ IIPAR, JJPAR, LLPAR
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_cmn_fj
!
! !DESCRIPTION: Subroutine INIT\_CMN\_FJ allocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Init_CMN_FJ( am_I_Root, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GIGC_SUCCESS

      ! Allocate arrays
      ALLOCATE( ZPJ( LLPAR, JPMAX, IIPAR, JJPAR ), STAT=RC )
      IF ( RC /= GIGC_SUCCESS ) THEN
         WRITE( 6, '(a)' ) 'Could not allocate ZPJ (cmn_fj_mod.F)!'
         RETURN
      ENDIF

      ! Zero arrays
      ZPJ = 0d0

      END SUBROUTINE Init_CMN_FJ
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_cmn_fj
!
! !DESCRIPTION: Subroutine CLEANUP\_CMN\_FJ deallocates all module arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE Cleanup_CMN_FJ( am_I_Root, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
!
! !INPUT PARAMETERS:
!
      LOGICAL, INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
!
! !OUTPUT PARAMETERS:
!
      INTEGER, INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  19 Nov 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      ! Assume success
      RC = GIGC_SUCCESS

      ! Deallocate arrays
      IF ( ALLOCATED( ZPJ ) ) DEALLOCATE( ZPJ )

      END SUBROUTINE Cleanup_CMN_FJ
!EOC
      END MODULE CMN_FJ_MOD

