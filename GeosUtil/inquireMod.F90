#if defined( DEVEL )
#include "MAPL_Generic.h"
!------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1    !
!------------------------------------------------------------------------
!BOP
!
! !MODULE: inquireMod
!
! !INTERFACE: 
!
MODULE inquireMod
! 
! !USES:
!
  USE ESMF_Mod
  USE MAPL_Mod

  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  :: findFreeLUN
  PUBLIC  :: I_Am_UnOPENed

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
!
! !INTERFACE:
!
 FUNCTION findFreeLUN(b) RESULT(lun)
!
! !USES:
!
  IMPLICIT NONE
!
! !INPUT ARGUMENTS:
!
  INTEGER, INTENT(IN), OPTIONAL :: b

!
! !RETURN ARGUMENTS:
!
  INTEGER :: lun

! !REVISION HISTORY:
!  14 Jun 2012: Nielsen     First crack
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!
  INTEGER :: i, rc, status
  INTEGER, PARAMETER :: iTop = 99

  LOGICAL :: exists ! File existence
  LOGICAL :: found  ! Detect unused logical unit
  LOGICAL :: open   ! Is open?

  CHARACTER(LEN=ESMF_MAXSTR) :: Iam
 
  Iam = "GEOSCHEMCHEM::findFreeLUN"
  status = 0
  rc = 0

! Find an available logical unit 
! ------------------------------
  found = .FALSE.
  i = 11

  DO WHILE (.NOT. found .AND. i <= iTop)
   INQUIRE(UNIT=i, EXIST=exists, OPENED=open)
   IF(exists .AND. .NOT. open) THEN
    found = .TRUE.
    lun = i
   END IF
   i = i + 1
  END DO

  IF(.NOT. found) THEN
   status = 1
   PRINT *,TRIM(Iam)//": No available logical units"
  END IF
  VERIFY_(status)

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
!
! !INTERFACE:
!
 FUNCTION I_Am_UnOPENed(n) RESULT(TorF)
!
! !USES:
!
  IMPLICIT NONE
!
! !INPUT ARGUMENTS:
!
  INTEGER :: n
!
! !RETURN ARGUMENTS:
!
  LOGICAL :: TorF
!
! !REVISION HISTORY:
!  14 Jun 2012: Nielsen     First crack
!EOP
!------------------------------------------------------------------------------
!BOC
! 
! !LOCAL VARIABLES:
!

  INTEGER :: rc, status

  LOGICAL :: exists ! File existence
  LOGICAL :: open   ! Is open?

  CHARACTER(LEN=ESMF_MAXSTR) :: Iam
 
  Iam = "GEOSCHEMCHEM::I_Am_UnOPENed"
  status = 0
  rc = 0

  INQUIRE(UNIT=n, EXIST=exists, OPENED=open)

  IF(exists .AND. .NOT. open) THEN
   TorF = .TRUE.
  ELSE
   TorF = .FALSE.
  END IF

 END FUNCTION I_Am_UnOPENed

!EOC
!------------------------------------------------------------------------------
END MODULE inquireMod
#endif
