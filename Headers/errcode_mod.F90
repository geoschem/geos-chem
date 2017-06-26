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
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: GC_Error
!
! !DEFINED PARAMETERS: 
!
  INTEGER, PUBLIC, PARAMETER :: GC_SUCCESS =  0   ! Routine returns success
  INTEGER, PUBLIC, PARAMETER :: GC_FAILURE = -1   ! Routine returns failure
!
! !REMARKS:
!  The error codes are returned by routines at various levels of GEOS-Chem.
!
! !REVISION HISTORY: 
!  19 Oct 2012 - R. Yantosca - Initial version
!  16 Aug 2016 - M. Sulprizio- Rename from gigc_errcode_mod.F90 to
!                              errcode_mod.F90. The "gigc" nomenclature is
!                              no longer used.
!  23 Jun 2017 - R. Yantosca - Moved subroutine GC_Error here
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GC_Error
!
! !DESCRIPTION: Subroutine GC\_Error prints an error message and sets RC to 
!  GC\_FAILURE. Note that this routine does not stop a run, but it will cause 
!  a stop at a higher level if you add a catch for RC /= GC\_SUCCESS. 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GC_Error( ErrMsg, RC, ThisLoc )
!
! !INPUT PARAMETERS:
!
    CHARACTER(LEN=*), INTENT(IN   )            :: ErrMsg 
    CHARACTER(LEN=*), INTENT(IN   ), OPTIONAL  :: ThisLoc 
!
! !INPUT/OUTPUT PARAMETERS:
!
    INTEGER,          INTENT(INOUT)            :: RC 
!
! !REVISION HISTORY:
!  13 Aug 2015 - E. Lundgren - Initial version, based on C. Keller's HCO_ERROR
!  16 Aug 2016 - M. Sulprizio- Rename from GIGC_ERROR to GC_ERROR
!  23 Jun 2017 - R. Yantosca - Now moved from error_mod.F to errcode_mod.F90
!EOP
!------------------------------------------------------------------------------
!BOC

    CHARACTER(LEN=255) :: Message

    !=======================================================================
    ! GC_ERROR begins here 
    !=======================================================================

    ! Print error message to log
    Message =  'GEOS-Chem ERROR: ' // TRIM(ErrMsg)
    WRITE( 6, '(a)' ) TRIM( Message )
      
    ! Print error location to log
    IF ( PRESENT( ThisLoc ) ) THEN
       Message = 'ERROR LOCATION: ' // TRIM( ThisLoc )
       WRITE( 6, '(a)' ) TRIM( ThisLoc )
    ENDIF

    ! Return with failure, but preserve existing error code
    !IF ( RC == GC_SUCCESS ) THEN
       RC = GC_FAILURE
    !ENDIF

  END SUBROUTINE GC_Error

END MODULE ErrCode_Mod
!EOC
