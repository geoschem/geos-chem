      PROGRAM MAIN

!-------------------------------------------------------------------------------
! PKUCPL.f95 provides utility to interrupt the normal DO-Loop of main.f
! and wait for a dispatcher to coordinate runs. (yanyy,6/18/14)
!-------------------------------------------------------------------------------

!Define variables
      IMPLICIT NONE
      CHARACTER*8      :: strarg
      CHARACTER        :: arg
      INTEGER          :: narg, iarg
      REAL*8           :: interval=0.5
      INTEGER,EXTERNAL :: ifsleep, ifnotend
      INTEGER,EXTERNAL :: ifsleeps, ifnotends
      LOGICAL          :: ifdone=.FALSE., ifnest(3)=.FALSE.
      CHARACTER*400    :: command,sed,combash
      CHARACTER*400    :: rundir,lockdir,codedir,pkucpl

      CALL getarg(2,codedir)
      CALL getarg(3,rundir)
      pkucpl = TRIM(codedir) // "PKUCPL/PKUCPL.sh"
      lockdir = TRIM(rundir) // "lock/"
      sed="sed '7c rundir='"
      CALL getarg(1,strarg)
      narg=LEN_TRIM(ADJUSTL(strarg))

      DO iarg=1,narg
      arg=strarg(iarg:iarg)

      SELECT CASE(arg)
      CASE("c")
              ifnest(1)=.true.
      CASE("n")
              ifnest(2)=.true.
      CASE("e")
              ifnest(3)=.true.
      CASE DEFAULT
      END SELECT

      ENDDO

      WRITE(combash,'(A,A)') TRIM(sed) // TRIM(rundir) // TRIM('' ) // TRIM(pkucpl)
      CALL SYSTEM(combash)
      WRITE(command,'(A,A)') TRIM(pkucpl) // " -",TRIM(ADJUSTL(strarg))
      CALL SYSTEM(command)

      DO WHILE( ifnotends(ifnest,lockdir) == 1 )
      DO WHILE( ifsleeps(ifnest,lockdir) == 1 &
               .AND. ifnotends(ifnest,lockdir) == 1 )
            CALL SLEEP(interval)
      ENDDO
         PRINT *,"Get done signal-----------"
         PRINT *,"------------------Let GEOS-Chem go--------------"
         CALL TOUCH_RM(ifnest)
      ENDDO

      PRINT *,"Congratulation! The simulations is complete----------"
      CALL SYSTEM("date")

!------------------------------------------------------------------------
!      Begin Subroutine Here
!------------------------------------------------------------------------

      CONTAINS

      SUBROUTINE TOUCH_RM(ifnest)
      IMPLICIT NONE
      LOGICAL, INTENT(IN)   :: ifnest(3)

      CALL SYSTEM("rm " // TRIM(lockdir) // "done*done")
      CALL SYSTEM("touch " // TRIM(lockdir) // "key.GLOBAL.key")

      IF (ifnest(1)) THEN
           CALL SYSTEM("touch " // TRIM(lockdir) // "key.NESTED_CH.key")
      END IF
      IF (ifnest(2)) THEN
           CALL SYSTEM("touch " // TRIM(lockdir) // "key.NESTED_NA.key")
      END IF
      IF (ifnest(3)) THEN
           CALL SYSTEM("touch " // TRIM(lockdir) // "key.NESTED_EU.key")
      ENDIF

      END SUBROUTINE TOUCH_RM
!                      End Subtouine Touch_rm
!-----------------------------------------------------------------------

      END PROGRAM MAIN
!                     End Program Main
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
!              Begin Funtion Here
!-----------------------------------------------------------------------

      FUNCTION ifsleeps(ifnest,lockdir) RESULT(ifsleep)
      IMPLICIT NONE
      LOGICAL, INTENT(IN)    :: ifnest(3)
      INTEGER                :: ifsleep
      LOGICAL                :: done_global=.FALSE.
      LOGICAL                :: done_nest(3)=.FALSE.,ifdone=.FALSE.
      CHARACTER*200          :: lockdir, done(4)

      done(1) = trim(lockdir) // "done.GLOBAL.done"
      done(2) = trim(lockdir) // "done.NESTED_CH.done"
      done(3) = trim(lockdir) // "done.NESTED_NA.done"
      done(4) = trim(lockdir) // "done.NESTED_EU.done"

      INQUIRE( FILE=TRIM(done(1)),EXIST=done_global )
      INQUIRE( FILE=TRIM(done(2)),EXIST=done_nest(1) )
      INQUIRE( FILE=TRIM(done(3)),EXIST=done_nest(2) )
      INQUIRE( FILE=TRIM(done(4)),EXIST=done_nest(3) )
      ifdone=done_global .AND. (ifnest(1) .EQV. done_nest(1)) .AND. &
                               (ifnest(2) .EQV. done_nest(2)) .AND. &
                               (ifnest(3) .EQV. done_nest(3))
      IF (ifdone) THEN
              ifsleep= 0
      ELSE
              ifsleep= 1
      END IF
      RETURN
      END FUNCTION ifsleeps

!----------------------------------------------------------------------

      FUNCTION ifnotends(ifnest,lockdir) RESULT(ifnotend)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: ifnest(3)
      LOGICAL             :: ifend=.false.
      INTEGER             :: ifnotend
      LOGICAL             :: check(3)=.false.,end_global=.false.
      CHARACTER*200       :: lockdir, endrun(4)

      endrun(1) = trim(lockdir) // "end.GLOBAL.end"
      endrun(2) = trim(lockdir) // "end.NESTED_CH.end"
      endrun(3) = trim(lockdir) // "end.NESTED_NA.end"
      endrun(4) = trim(lockdir) // "end.NESTED_EU.end"

      INQUIRE( FILE=TRIM(endrun(1)),EXIST=end_global )
      INQUIRE( FILE=TRIM(endrun(2)),EXIST=check(1) )
      INQUIRE( FILE=TRIM(endrun(3)),EXIST=check(2) )
      INQUIRE( FILE=TRIM(endrun(4)),EXIST=check(3) )
      ifend=end_global .AND. (ifnest(1) .EQV. check(1)) .AND. &
                             (ifnest(2) .EQV. check(2)) .AND. &
                             (ifnest(3) .EQV. check(3))

      IF (ifend) THEN
              ifnotend = 0
      ELSE
              ifnotend = 1
      ENDIF
      RETURN
      END FUNCTION ifnotends

!-------------------------------------------------------------------------
!      End Function
!-------------------------------------------------------------------------
