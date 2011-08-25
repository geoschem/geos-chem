!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_DEP_mod
!
! !DESCRIPTION: Common blocks for dry deposition.
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_DEP_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR, NTYPE

      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC DATA MEMBERS:
!
      ! IREG   - Number of landtypes in grid square (I,J)
      INTEGER :: IREG(IIPAR,JJPAR)                                

      ! ILAND  - Land type ID in grid square (I,J) for IREG landtypes
      INTEGER :: ILAND(IIPAR,JJPAR,NTYPE)

      ! IUSE   - Fraction ((per mil) of gridbox area occupied by land type
      INTEGER :: IUSE(IIPAR,JJPAR,NTYPE)                          

      ! Fraction of land in grid box
      REAL*8  :: FRCLND(IIPAR,JJPAR)

      ! XLAI    - Leaf Area Index of land type for current MONTH
      REAL*8  :: XLAI(IIPAR,JJPAR,NTYPE)

      ! XLAI2   - Leaf Area Index of land type for following MONTH
      REAL*8  :: XLAI2(IIPAR,JJPAR,NTYPE)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE IS DEPRECATED.  IT WILL BE REMOVED WHEN THE  %%%
!  %%%  GEIA BIOGENIC EMISSIONS ARE REMOVED FROM GEOS-CHEM (v9-01-03)  %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  23 Aug 2011 - M. Long     - Converted to Module from Header file
!  25 Aug 2011 - R. Yantosca - Change IGLOB,JGLOB to IIPAR,JJPAR
!EOP
!-----------------------------------------------------------------------------
!BOC
	END MODULE CMN_DEP_MOD
!EOC
