! ******************************************************************
! ------------------------------------------------------------------
MODULE MESSY_NCREGRID_MPI
! ------------------------------------------------------------------
! Author: Patrick Joeckel, MPICH, Mainz, April 2004
! ******************************************************************
#if defined(MPI)

  IMPLICIT NONE
  PRIVATE

  INCLUDE 'mpif.h'

  PUBLIC   :: ncregrid_abort

CONTAINS

! --------------------------------------------------------------
  SUBROUTINE ncregrid_abort (name, text, exit_no)

    CHARACTER(*)           :: name
    CHARACTER(*), OPTIONAL :: text
    INTEGER, OPTIONAL      :: exit_no
    INTEGER                :: iexit
    INTEGER                :: iou
    LOGICAL                :: opened

    IF (PRESENT(exit_no)) THEN
       iexit = exit_no
    ELSE
       iexit = 1
    END IF

    IF (iexit /=0) WRITE (*,'(/,80("*"),/)')

    IF (PRESENT(text)) THEN
       WRITE (*,'(1x,a,a,a)') TRIM(name), ': ', TRIM(text)
    ELSE
       WRITE (*,'(1x,a,a)') TRIM(name), ': '
    ENDIF

    IF (iexit /=0) WRITE (*,'(/,80("*"),/)')

#if defined(MESSY)
    ! WRITE FILE 'END' FOR BREAKING RERUN-CHAIN IN RUNSCRIPT xmessy
    DO iou=100,300
       INQUIRE(unit=iou,opened=opened)
       IF (.NOT.opened) EXIT
    END DO
    OPEN(iou, FILE='END', STATUS='unknown')
    IF (PRESENT(text)) THEN
       WRITE (iou,'(1x,a,a,a)') TRIM(name), ': ', TRIM(text)
    ELSE
       WRITE (iou,'(1x,a,a)') TRIM(name), ': '
    ENDIF
    CLOSE(iou)
#endif

    CALL p_abort

  END SUBROUTINE ncregrid_abort
! --------------------------------------------------------------

! --------------------------------------------------------------
  SUBROUTINE p_abort

    ! this routine should be used instead of abort, util_abort() or STOP
    ! in all routines for proper clean up of all PEs

    IMPLICIT NONE

    EXTERNAL :: MPI_ABORT

    INTEGER :: p_error

    CALL MPI_ABORT (MPI_COMM_WORLD, 1, p_error)

    IF (p_error /= MPI_SUCCESS) THEN
       WRITE (*,'(a)') ' MPI_ABORT failed.'
       WRITE (*,'(a,i4)') ' Error =  ', p_error
       STOP
    END IF

  END SUBROUTINE p_abort
! --------------------------------------------------------------

#endif
! ******************************************************************
END MODULE MESSY_NCREGRID_MPI
! ******************************************************************
