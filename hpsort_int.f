C $Id: hpsort_int.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      SUBROUTINE HPSORT_INT ( N, XX )
!
!-----------------------------------------------------------------------------
!  HPSORT_INT does a heap sort on integer array XX(N) in place (cf. Press)
!
!  Recreated for GEOS-CTM (bmy, 7/7/98) from original HPSORT.F  
!-----------------------------------------------------------------------------
!
      IMPLICIT NONE

! Arguments
      INTEGER, INTENT(IN   ) :: N
      INTEGER, INTENT(INOUT) :: XX(N)

! Local variables
      INTEGER :: I, J, L, IR
      INTEGER :: XXA
!
! HPSORT_INT begins here!!
!
      IF ( N .lt. 2 ) RETURN
      L  = N/2 + 1
      IR = N

   10 CONTINUE
      IF( L .gt. 1 ) THEN
         L   = L-1
         XXA = XX(L)
      ELSE
         XXA    = XX(IR)
         XX(IR) = XX(1)
         IR     = IR-1

         IF ( IR .eq. 1 ) THEN
            XX(1) = XXA
            RETURN
         ENDIF

      ENDIF

      I = L
      J = L + L
 20   IF( J .le. IR ) THEN
         IF ( J .lt. IR ) THEN
            IF ( XX(J) .lt. XX(J+1) ) J = J + 1
         ENDIF

         IF( XXA .lt. XX(J) ) THEN
            XX(I) = XX(J)
            I     = J
            J     = J + J
         ELSE
            J = IR + 1
         ENDIF

         GOTO 20
      ENDIF

      XX(I) = XXA

      GOTO 10

      END SUBROUTINE HPSORT_INT
