!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
! !MODULE: CMN_NOX_MOD
!
! !DESCRIPTION: CMN_NOX is the header file for containing NOx from aircraft,
!       lightning, and soils.  CMN_NOX also contains variables for the
!       lightning subroutines. (mgs, bey, bdf, bmy, 3/5/98, 10/2/07)
!
! REMARKS:
! (1 ) Changed RCS ID tags from "C" to "!" to allow 
! freeform compilation.  (bmy, 6/25/02)
! (2 ) Moved BXHEIGHT to "dao_mod.f".  The fact that BXHEIGHT
! was in "CMN_NOX" is historical baggage. (bmy, 9/18/02)
! (3 ) Now everything except GEMISNOX, GEMISNOX2 is in
! "lightning_mod.f" (bmy, 4/14/04)
! (4 ) Remove GEMISNOX from common block (ltm, bmy, 10/2/07)
!=================================================================
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
!      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS: KITCHEN SINK!
!
      PUBLIC
	SAVE
! !DEFINED PARAMETERS: 
!
      ! NOTE: Keep this for backwards compatibility for now (bmy, 10/2/07)
      ! GEMISNOX2 = Soil Nox [molec NOx/cm3/s]
      REAL*8            :: GEMISNOX2(IIPAR,JJPAR)
	END MODULE CMN_NOX_MOD
