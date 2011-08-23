!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
! !MODULE: CMN_VEL_mod.F
!
! !DESCRIPTION: COMMON BLOCKS FOR DEPVEL
!
!  !REVISION HISTORY:
!  23 Aug 2011 - M. Long   - Converted to Module from Header file
!------------------------------------------------------------------------------
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_VEL_MOD
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
!
! !REMARKS:
! (1 ) Changed RCS ID tag comment character from "C" to "!" to allow freeform
!       compilation.  Also added & continuation characters in column 73
!       to allow header files to be included in F90 freeform files.
!       (bmy, 6/25/02)
! (2 ) Commented out internal variables that are local to DEPVEL since
!       these are now contained within "drydep_mod.f" (bmy, 11/21/02)

      ! NOTE: This should be split off into a "landtype_mod.f"
      INTEGER :: IJREG(MAXIJ),IJLAND(MAXIJ,NTYPE),                      &
     &           IJUSE(MAXIJ,NTYPE)

      ! NOTE: This should be split off into a "lai_mod.f"
      REAL*8  :: XYLAI(MAXIJ,NTYPE)
	END MODULE CMN_VEL_MOD
 
