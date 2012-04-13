!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: CMN_NOX_mod
!
! !DESCRIPTION: CMN\_NOX\_mod is the module file for containing NOx from soils.
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_NOX_MOD
!
! !USES:
!
      USE CMN_SIZE_MOD, ONLY : IIPAR, JJPAR

      IMPLICIT NONE
      PUBLIC
!
! !PUBLIC DATA MEMBERS:
! 
      ! NOTE: Keep this for backwards compatibility for now (bmy, 10/2/07)
      ! GEMISNOX2 = Soil Nox [molec NOx/cm3/s]
      REAL*8 :: GEMISNOX2(IIPAR,JJPAR)

! !REMARKS:
!  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  !%%%  NOTE: THIS MODULE IS DEPRECATED.  IT CAN BE REMOVED WHEN THE  %%%
!  !%%%  UPDATED SOIL NOX EMISSIONS ARE FINALLY ADDED TO GEOS-CHEM.    %%%
!  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! !REVISION HISTORY:
!  05 Mar 1998 - M. Schultz - Initial version
!  (1 ) Changed RCS ID tags from "C" to "!" to allow 
!        freeform compilation.  (bmy, 6/25/02)
!  (2 ) Moved BXHEIGHT to "dao_mod.f".  The fact that BXHEIGHT
!        was in "CMN_NOX" is historical baggage. (bmy, 9/18/02)
!  (3 ) Now everything except GEMISNOX, GEMISNOX2 is in
!        "lightning_mod.f" (bmy, 4/14/04)
!  (4 ) Remove GEMISNOX from common block (ltm, bmy, 10/2/07)
!  23 Aug 2011 - M. Long    - Converted to Module from Header file
!-----------------------------------------------------------------------------
!BOC
	END MODULE CMN_NOX_MOD
!EOC
