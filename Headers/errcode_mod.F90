!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: errcode_mod.F90
!
! !DESCRIPTION: Module ERRCODE\_MOD contains the error codes (i.e. that
!  report success or failure) returned by GEOS-Chem routines.
!\\
!\\
! !INTERFACE: 
!
MODULE ErrCode_Mod
!
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !DEFINED PARAMETERS: 
!
  INTEGER, PUBLIC, PARAMETER :: GC_SUCCESS =  0   ! Routine returns success
  INTEGER, PUBLIC, PARAMETER :: GC_FAILURE = -1   ! Routine returns failure
!
! !REMARKS:
!  The error codes are returned by routines at various levels of the 
!  Grid-Independent GEOS-Chem implementation.
!
! !REVISION HISTORY: 
!  19 Oct 2012 - R. Yantosca - Initial version
!  16 Aug 2016 - M. Sulprizio- Rename from gigc_errcode_mod.F90 to
!                              errcode_mod.F90. The "gigc" nomenclature is
!                              no longer used.
!EOP
!------------------------------------------------------------------------------
!BOC
END MODULE ErrCode_Mod
!EOC
