!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_MONOT_mod
!
! !DESCRIPTION: This include file contains variables for the old GEIA
!  biogenic emissions (MONOTERPENES).
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_MONOT_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : MAXIJ, NTYPE

      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC DATA MEMBERS:
!
      ! Baseline emission
      REAL*8, ALLOCATABLE :: BASEMONOT(:,:)
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
      CONTAINS

      SUBROUTINE SET_CMN_MONOT_MOD

      IMPLICIT NONE
      
      INTEGER AS

      ALLOCATE(
     &     BASEMONOT(MAXIJ,NTYPE),
     &     STAT=AS)

      END SUBROUTINE SET_CMN_MONOT_MOD
      
      END MODULE CMN_MONOT_MOD
!EOC
