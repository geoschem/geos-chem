!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: smv_dimension_mod.F
!
! !DESCRIPTION: This include file contains the various placeholder parameters
!  that are required to replace references to GEOS-Chem grid parameters.  
!  This is necessary because several quantities in the FAST-J and SMVGEAR
!  codes are contained in common blocks, and we need to have these parameters
!  for sizing those arrays properly.
!\\
!\\
! !INTERFACE: 
!
      MODULE SMV_DIMENSION_MOD
!
! !USES:
!
      IMPLICIT NONE
      PUBLIC
!
! !DEFINED PARAMETERS: 
!
      ! Locally defined replacement for GEOS-Chem parameter "LLPAR"
      !INTEGER, PARAMETER :: MAX_COLUMN  = 72    ! Full GEOS-5 vertical grid
      INTEGER, PARAMETER :: MAX_COLUMN  = 47    ! Reduced GEOS-5 vertical grid

      ! Locally defined replacement for GEOS-Chem parameter "NNPAR"
      INTEGER, PARAMETER :: MAX_TRACERS = 100

      ! Locally defined replacement for "comode.h" parameter "IGAS"
      INTEGER, PARAMETER :: MAX_SPECIES = 125
!
! !REMARKS:
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%  NOTE: THIS MODULE WAS ORIGINALLY DEVELOPED FOR THE COLUMN CODE.  %%%
!  %%%  THE PARAMETERS HERE CAN BE COMBINED INTO A SINGLE MODULE WITH    %%%
!  %%%  DECLARED SIZE PARAMETERS (IN THE NEAR FUTURE).                   %%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY: 
!  24 Mar 2009 - R. Yantosca - Initial version
!  16 Apr 2010 - R. Yantosca - Added MAX_SPECIES = 125
!  03 Aug 2011 - M. Long     - Converted from Header file to Module
!EOP
	END MODULE SMV_DIMENSION_MOD
!------------------------------------------------------------------------------
