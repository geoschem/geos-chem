#if defined (ESMF_)
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gc_setvals_mod
!
! !DESCRIPTION: Module GC\_IDXMAP\_MOD contains the routine
!  to initialize the mapping between ESMF tracer indexes and
!  GEOS-Chem internal indexes.
!\\
!\\
! !INTERFACE: 
!      
MODULE GC_SETVALS_MOD

  IMPLICIT NONE
  PUBLIC
!
! !REMARKS:
!  This module appears to be not used in the current implementation of
!  the GEOS-Chem ESMF interface.
! 
! !REVISION HISTORY: 
!  16 Oct 2012 - M. Long     - Initial version
!  16 Oct 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !DESCRIPTION: 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_SETVALS()
!
! !USES:
!
    IMPLICIT NONE
!
! !REMARKS:
!  This routine appears to be not used in the current implementation of
!  the GEOS-Chem ESMF interface.
! 
! !REVISION HISTORY: 
!  16 Oct 2012 - M. Long     - Initial version
!  16 Oct 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC        
  END SUBROUTINE GC_SETVALS
!eoc
END MODULE GC_SETVALS_MOD
#endif
