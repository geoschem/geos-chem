!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: gigc_errcode_mod.F90
!
! !DESCRIPTION: Module GIGC\_ERRCODE\_MOD contains the error codes (i.e. that
!  report success or failure) returned by routines of the Grid-Independent
!  GEOS-Chem (aka "GIGC").
!\\
!\\
! !INTERFACE: 
!
MODULE GIGC_ErrCode_Mod
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !DEFINED PARAMETERS: 
!
  INTEGER, PUBLIC, PARAMETER :: GIGC_SUCCESS =  0   ! Routine returns success
  INTEGER, PUBLIC, PARAMETER :: GIGC_FAILURE = -1   ! Routine returns failure
!
! !REMARKS:
!  The error codes are returned by routines at various levels of the 
!  Grid-Independent GEOS-Chem implementation.
!
! !REVISION HISTORY: 
!  19 Oct 2012 - R. Yantosca - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
END MODULE GIGC_ErrCode_Mod
!EOC
