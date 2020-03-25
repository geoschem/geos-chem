#ifdef ESMF_
! We only need to refer to this include file if we are connecting
! to the GEOS-5 GCM via the ESMF/MAPL framework (bmy, 8/3/12)
#include "MAPL_Generic.h"
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
#ifdef ESMF_
  ! We only need to refer to these modules if we are connecting
  ! to the GEOS-5 GCM via the ESMF/MAPL framework (bmy, 8/3/12)
  USE ESMF
  USE MAPL_Mod
#endif

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: findFreeLUN
  PUBLIC  :: I_Am_UnOPENed
!
! !REVI<SION HISTORY:
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

#ifdef ESMF_
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam
#else
    CHARACTER(LEN=255)         :: Iam
#endif
!
! !DEFINED PARAMETERS
!
    INTEGER, PARAMETER         :: iTop = 199     ! Maximum LUN limit

    !======================================================================
    ! Initialization
    !======================================================================
    Iam = "GEOSCHEMCHEM::findFreeLUN"
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

#ifdef ESMF_
    VERIFY_(status)
#endif

  END FUNCTION findFreeLUN
!EOC
!------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1    !
!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: I_Am_UnOPENed
!
! !DESCRIPTION: Inquire as to the availability of a given logical unit
!\\
!\\
! !INTERFACE:
!
  FUNCTION I_Am_UnOPENed( n ) RESULT( TorF )
!
! !USES:
    !
    IMPLICIT NONE
!
! !INPUT PARAMETERS:
!
    INTEGER :: n     ! Logical unit # to test
!
! !RETURN VALUE:
!
    LOGICAL :: TorF  ! .TRUE. means the file is unopened
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
    INTEGER                    :: rc, status
    LOGICAL                    :: exists       ! File existence
    LOGICAL                    :: open         ! Is open?

#ifdef ESMF_
    CHARACTER(LEN=ESMF_MAXSTR) :: Iam
#else
    CHARACTER(LEN=255)         :: Iam
#endif

    !======================================================================
    ! Initialization
    !======================================================================
    Iam    = "GEOSCHEMCHEM::I_Am_UnOPENed"
    status = 0
    rc     = 0

    !======================================================================
    ! Inquire if the LUN is available
    !======================================================================
    INQUIRE( UNIT=n, EXIST=exists, OPENED=open )

    IF ( exists .AND. .NOT. open ) THEN
       TorF = .TRUE.
    ELSE
       TorF = .FALSE.
    ENDIF

  END FUNCTION I_Am_UnOPENed
!EOC
END MODULE inquireMod
