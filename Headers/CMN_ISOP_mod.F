!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_ISOP_mod
!
! !DESCRIPTION: This include file contains variables for the old GEIA
!  biogenic emissions (ISOPRENE).
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_ISOP_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : MAXIJ, NPOLY, NTYPE

      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC DATA MEMBERS:
!
      ! polynomial fitting's coeffcients
      REAL*8 :: SOPCOEFF(NPOLY)
      REAL*8 :: BASEISOP(MAXIJ,NTYPE)

      ! Added for EMISOP_GRASS and EMISOP_MB (bdf, bmy, 8/1/01)
      REAL*8 :: BASEGRASS(MAXIJ)
      REAL*8 :: BASEMB(MAXIJ)
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE IS DEPRECATED.  IT WILL BE REMOVED WHEN THE  %%%
!  %%%  GEIA BIOGENIC EMISSIONS ARE REMOVED FROM GEOS-CHEM (v9-01-03)  %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  23 Aug 2011 - M. Long   - Converted to Module from Header file
!EOP
!------------------------------------------------------------------------------
!BOC
	END MODULE CMN_ISOP_MOD
!EOC
