#ifdef MAPL_ESMF
#ifdef MAPL3
#include "MAPL.h"
#else
#include "MAPL_Generic.h"
#endif
#endif
!------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1    !
!------------------------------------------------------------------------
!BOP
!
! !MODULE: inquireMod
!
! !DESCRIPTION: Module inquireMod contains functions to find free and
!  unopened logical file units (LUNs) for Fortran I/O.
!
! !INTERFACE:
!
MODULE inquireMod
!
! !USES:
!
#ifdef MAPL_ESMF
  USE ESMF
#ifdef MAPL3
  USE mapl_ErrorHandlingMod, only: MAPL_Verify
#else
  USE MAPL_Mod
#endif
#endif

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: findFreeLUN
!
! !REVISION HISTORY:
!  14 Jun 2012 - E. Nielsen  - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
  CONTAINS
!EOC
!------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1    !
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: findFreeLUN
!
! !DESCRIPTION: Inquire for an existing, but unopened, logical unit number
!\\
!\\
! !INTERFACE:
!
  FUNCTION findFreeLUN( b ) RESULT( lun )
!
! !USES:
!
#if defined( MODEL_CESM )
    USE UNITS,      ONLY : GETUNIT
#endif
    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    INTEGER, INTENT(IN), OPTIONAL :: b   ! Not really used here
!
! !RETURN VALUE:
!
    INTEGER :: lun
!
! !REVISION HISTORY:
!  14 Jun 2012 - E. Nielsen  - Initial version
!  See https://github.com/geoschem/geos-chem for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER                    :: i, rc, status
    LOGICAL                    :: exists        ! File existence
    LOGICAL                    :: found         ! Detect unused logical unit
    LOGICAL                    :: open          ! Is open?
#ifdef MAPL_ESMF
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam
#else
    CHARACTER(LEN=255)         :: Iam
#endif

!
! !DEFINED PARAMETERS
!
    INTEGER, PARAMETER         :: iTop = 199     ! Maximum LUN limit

#if defined( MODEL_CESM )
    lun = GETUNIT()
#else
    !======================================================================
    ! Initialization
    !======================================================================
    Iam = "GEOS-Chem::findFreeLUN"
    status = 0
    rc     = 0

    !======================================================================
    ! Find an available logical unit
    !======================================================================
    found = .FALSE.
    i     = 11

    DO WHILE ( .NOT. found .AND. i <= iTop )
       INQUIRE( UNIT=i, EXIST=exists, OPENED=open )
       IF ( exists .AND. .NOT. open ) THEN
          found = .TRUE.
          lun = i
       ENDIF
       i = i + 1
    ENDDO

    IF ( .NOT. found ) THEN
       status = 1
       PRINT *,TRIM( Iam ) // ": No available logical units"
    ENDIF
#endif

#ifdef MAPL_ESMF
#ifdef MAPL3
    ! Comment out for now (did this in HEMCO too)
    !_VERIFY(status)
#else
    VERIFY_(status)    
#endif
#endif

  END FUNCTION findFreeLUN
!EOC
END MODULE inquireMod
