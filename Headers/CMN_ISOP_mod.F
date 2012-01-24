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
#if defined( DEVEL )
      REAL*8, ALLOCATABLE :: BASEISOP(:,:)
#else
      REAL*8 :: BASEISOP(MAXIJ,NTYPE)
#endif

      ! Added for EMISOP_GRASS and EMISOP_MB (bdf, bmy, 8/1/01)
#if defined( DEVEL )
      REAL*8, ALLOCATABLE :: BASEGRASS(:)
      REAL*8, ALLOCATABLE :: BASEMB(:)
#else
      REAL*8 :: BASEGRASS(MAXIJ)
      REAL*8 :: BASEMB(MAXIJ)
#endif
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
#if defined( DEVEL )
      CONTAINS

      SUBROUTINE SET_CMN_ISOP_MOD

      IMPLICIT NONE
      
      INTEGER AS

      ALLOCATE(
     &     BASEISOP(MAXIJ,NTYPE),
     &     BASEGRASS(MAXIJ),
     &     BASEMB(MAXIJ),
     &     STAT=AS)

      END SUBROUTINE SET_CMN_ISOP_MOD
#endif      
      END MODULE CMN_ISOP_MOD
!EOC
