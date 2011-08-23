!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
! !MODULE: CMN_MONOT_MOD
!
! !DESCRIPTION: 
!
!  !REVISION HISTORY:
!  23 Aug 2011 - M. Long   - Converted to Module from Header file
!------------------------------------------------------------------------------
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_MONOT_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD

      IMPLICIT NONE
!      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS: KITCHEN SINK!
!
      PUBLIC
	SAVE
! !DEFINED PARAMETERS: 
!
!-----COMMON BLOCKS FOR MONOTERPENE
!
! (1 ) Changed RCS ID tag comment character from "C" to "!" to allow freeform
!       compilation.  Updated comments, cosmetic changes. (bmy, 6/25/02)

      ! Baseline emission
      REAL*8 BASEMONOT(MAXIJ,NTYPE)

	END MODULE CMN_MONOT_MOD
