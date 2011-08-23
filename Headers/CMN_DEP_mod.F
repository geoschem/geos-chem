!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
! !MODULE: CMN_DEP_mod.F
!
! !DESCRIPTION: COMMON BLOCKS FOR DRYDEP
!
! REMARKS:
! (1 ) Changed RCS ID tag comment character from "C" to "!" to allow freeform
!       compilation.  Changed continuation character from numbers
!       to &.  Changed comment character from C to !.  Added & continuation 
!       characters in column 73 to allow header files to be included in 
!       F90 freeform files. (bmy, 6/25/02)*
!
!  !REVISION HISTORY:
!  23 Aug 2011 - M. Long   - Converted to Module from Header file
!-----------------------------------------------------------------------------
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_DEP_MOD
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
      !  IREG    - Number of landtypes in grid square (I,J)
      !  ILAND   - Land type ID in grid square (I,J) for IREG landtypes
      !  IUSE    - Fraction ((per mil) of gridbox area occupied by
      !                          land type
      INTEGER ::      IREG(IGLOB,JGLOB),                                &
     &                ILAND(IGLOB,JGLOB,NTYPE),                         &
     &                IUSE(IGLOB,JGLOB,NTYPE)                          
      REAL*8  ::      FRCLND(IGLOB,JGLOB)

      ! XLAI    - Leaf Area Index of land type for current MONTH
      REAL*8 :: XLAI(IIPAR,JJPAR,NTYPE)

      ! XLAI2   - Leaf Area Index of land type for following MONTH
      REAL*8 :: XLAI2(IIPAR,JJPAR,NTYPE)
	END MODULE CMN_DEP_MOD
