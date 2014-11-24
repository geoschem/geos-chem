!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: hco_error_mod 
!
! !DESCRIPTION: Module HCO\_ERROR\_MOD contains all routines for error
! handling in HEMCO. 
! \\
! \\
! !INTERFACE: 
!
      MODULE HCO_ERROR_MOD
!
! !USES:
!
#if defined( ESMF_ )
      USE ESMF
      USE MAPL_MOD
#endif
      IMPLICIT NONE
      PRIVATE

#if defined( ESMF_ )
#     include "MAPL_Generic.h"
#endif
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC           :: HCO_STOP 
      PUBLIC           :: HCO_ERROR
      PUBLIC           :: HCO_WARNING 
!
! !MODULE VARIABLES:
!
      INTEGER, PARAMETER, PUBLIC  :: HCO_SUCCESS = 0
      INTEGER, PARAMETER, PUBLIC  :: HCO_FAIL    = -999
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
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
! !ROUTINE: HCO_STOP
!
! !DESCRIPTION: Subroutine HCO_STOP stops HEMCO at the given location.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_STOP
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC
      INTEGER :: STAT

      !======================================================================
      ! HCO_STOP begins here 
      !======================================================================

#if defined( ESMF_ )
      !IF( MAPL_VRFYt(RC,ErrMsg,ErrLoc,1,STAT) ) CALL MAPL_Abort
      CALL MAPL_Abort
#else
      ! Cleanup HEMCO
      CALL HCO_CLEANUP
      CALL EXIT( 99999 )
#endif

      END SUBROUTINE HCO_STOP
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_ERROR
!
! !DESCRIPTION: Subroutine HCO_ERROR promts an error message and sets RC to 
! -999 (/= HCO_SUCCESS). Note that this routine does not stop a run, but it 
! will cause a stop at higher level (when RC gets evaluated). 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_ERROR ( ErrMsg, ErrLoc, RC )
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*), INTENT(IN   )  :: ErrMsg 
      CHARACTER(LEN=*), INTENT(IN   )  :: ErrLoc 
      INTEGER,          INTENT(INOUT)  :: RC 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! HCO_ERROR begins here 
      !======================================================================

      ! Print error message
      WRITE( *, '(a)' ) REPEAT( '=', 79 )
      WRITE( *, '(a)' ) 'HEMCO ERROR: ' // TRIM( ErrMsg )
      WRITE( *, '(a)' ) 'LOCATION   : ' // TRIM( ErrLoc )
      WRITE( *, '(a)' ) REPEAT( '=', 79 )

      ! Return w/ error
      RC = HCO_FAIL 

      END SUBROUTINE HCO_ERROR
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: HCO_WARNING
!
! !DESCRIPTION: Subroutine HCO_WARNING promts a warning message without 
! forcing HEMCO to stop. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE HCO_WARNING ( ErrMsg, ErrLoc, RC )
!
! !ARGUMENTS:
!
      CHARACTER(LEN=*), INTENT(IN   )  :: ErrMsg 
      CHARACTER(LEN=*), INTENT(IN   )  :: ErrLoc 
      INTEGER,          INTENT(INOUT)  :: RC 
!
! !REVISION HISTORY:
!  23 Sep 2013 - C. Keller - Initialization
!
!EOP
!------------------------------------------------------------------------------
!BOC

      !======================================================================
      ! HCO_WARNING begins here 
      !======================================================================

      ! Print warning
      WRITE( *, '(a)' ) REPEAT( '=', 79 )
      WRITE( *, '(a)' ) 'HEMCO WARNING: ' // TRIM( ErrMsg )
      WRITE( *, '(a)' ) 'LOCATION   : ' // TRIM( ErrLoc )
      WRITE( *, '(a)' ) REPEAT( '=', 79 )

      ! Return w/ success
      RC = HCO_SUCCESS

      END SUBROUTINE HCO_WARNING
!EOC
      END MODULE HCO_ERROR_MOD
!EOM
