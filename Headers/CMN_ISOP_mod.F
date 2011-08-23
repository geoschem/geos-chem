!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
! !MODULE: CMN_ISOP_MOD
!
! !DESCRIPTION: This include file contains physical constants for the
!  GEOS-Chem column chemistry code. 
!
!  !REVISION HISTORY:
!  23 Aug 2011 - M. Long   - Converted to Module from Header file
!------------------------------------------------------------------------------
!\\
!\\
! !INTERFACE: 
!
      MODULE CMN_ISOP_MOD
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
!-----COMMON BLOCKS FOR ISOPRENE
!
! NOTES:
! (1 ) Changed RCS ID tag comment character from "C" to "!" to allow freeform
!       compilation.  Also added & continuation characters in column 73
!       to allow header files to be included in F90 freeform files.
!       (bmy, 6/25/02)
!
      REAL*8 
             ! polynomial fitting's coeffcients
     &       SOPCOEFF(NPOLY),BASEISOP(MAXIJ,NTYPE),                     &
             ! Added for EMISOP_GRASS and EMISOP_MB (bdf, bmy, 8/1/01)
     &       BASEGRASS(MAXIJ),BASEMB(MAXIJ)
	END MODULE CMN_ISOP_MOD
