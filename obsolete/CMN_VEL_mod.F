!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_VEL_mod.F
!
! !DESCRIPTION: COMMON BLOCKS FOR DEPVEL
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_VEL_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : MAXIJ, NTYPE

      IMPLICIT NONE
      PUBLIC

      ! NOTE: This should be split off into a "landtype_mod.f"
      INTEGER :: IJREG(MAXIJ)
      INTEGER :: IJLAND(MAXIJ,NTYPE)
      INTEGER :: IJUSE(MAXIJ,NTYPE)

      ! NOTE: This should be split off into a "lai_mod.f"
      REAL*8  :: XYLAI(MAXIJ,NTYPE)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE IS DEPRECATED.  THE ARRAYS CONTAIN THE SAME  %%%
!  %%%  INFORMATION AS IN "CMN_DEP_mod.F", BUT ARE ARCHIVED AS 1-D     %%%
!  %%%  ARRAYS (OF SIZE MAXIJ=IIPAR*JJPAR).  THIS IS BECAUSE THE VERY  %%%
!  %%%  OLD DRYDEP ROUTINE "DEPVEL" WAS ORIGINALLY WRITTEN TO EXPECT   %%%
!  %%%  1-D ARRAYS AND HAS NOT BEEN UPDATED SINCE.  THIS MODULE CAN    %%%
!  %%%  BE REMOVED WHEN THE DRYDEP MODULE IS UPDATED (IN THE FUTURE).  %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! !REVISION HISTORY:
!  23 Aug 2011 - M. Long   - Converted to Module from Header file
!EOP
!-------------------------X-----------------------------------------------------
!BOC
 	END MODULE CMN_VEL_MOD
!EOC
