!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      ! Global array sizes in longitude, latitude, altitude
!      INTEGER :: IPAR
!      INTEGER :: JPAR
!      INTEGER :: LPAR

      ! max # of photolysis rxns = 4 + IPHOT (see comode.h)
      ! FP increased JPMAX since IPHOT was increased (hotp 7/31/09)
      !INTEGER, PARAMETER :: JPMAX = 79
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
!EOP
!------------------------------------------------------------------------------
!BOC
      CONTAINS

      SUBROUTINE SET_CMN_FJ_MOD

      IMPLICIT NONE

      INTEGER AS

!      IPAR = IIPAR
!      JPAR = JJPAR
!      LPAR = LLPAR

      ALLOCATE(ZPJ(LLPAR,JPMAX,IIPAR,JJPAR), STAT=AS)

      END SUBROUTINE SET_CMN_FJ_MOD

      END MODULE CMN_FJ_MOD
!EOC
