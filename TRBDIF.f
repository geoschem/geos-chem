C $Id: TRBDIF.f,v 1.1 2003/06/30 20:26:03 bmy Exp $
      SUBROUTINE TRBDIF ( XX1,RHOKDZ,FLXFAC,NLEV,LTYPE,EPSL,irun )
C
C**********************************************************************
C
C   ARGUMENTS ::
C
C     INPUT:
C     ------
C    XX1           -         PROPERTY TO BE DIFFUSED
C    RHOKDZ        -         RHO * K * WEIGHT / DZ AT INTERFACES
C    FLXFAC        -         G * DT / (DP*WEIGHT) AT EDGES
C    NLEV          -         NUMBER OF ATMOSPHERIC LEVELS
C    LTYPE         -         LOGICAL FLAG FOR UNDERFLOW CUTOFF
C    EPSL          -         UNDERFLOW CUTOFF VALUE
C
C     OUTPUT:
C     ------
C    XX1           -         NEW VALUE RETURNED
C
C**********************************************************************
C
C   VALUES PASSED FROM diffuse
C
C    XX1     <--  XX
C    RHOKDZ  <--  RHOKDZ
C    FLXFAC  <--  FLXFAC
C    NLEV    <--  LM
C    LTYPE   <--  LTYPE
C    EPSL    <--  EPS
C    IRUN    <--  JJPAR
C
C**********************************************************************
C  Modification notes
C  ==================
C  (1) 1-d loops over 2-d arrays changed to 2-d loops (bdf, 2/99)
C  (2) implicit none added (bdf 2/24/99)
C  (3) Now use double-precision exponents (bmy, 3/23/99)
C  (4) Now use F90 syntax for declarations (bmy, 4/6/99)
C  (5) Eliminate ISTNLV, ISTNM1, NLEVP1...these are the old
C      CRAY loop boundary limits (bmy, 4/6/99)
C  (6) Use F90 WHERE statement to handle the underflow test (bmy, 4/8/99)
C**********************************************************************

      IMPLICIT NONE

!**** Arguments
      LOGICAL, INTENT(IN)    :: LTYPE

      INTEGER, INTENT(IN)    :: IRUN, NLEV

      REAL*8,  INTENT(INOUT) :: XX1(IRUN,NLEV+1)
      REAL*8,  INTENT(IN)    :: RHOKDZ(IRUN,NLEV)
      REAL*8,  INTENT(IN)    :: FLXFAC(IRUN,NLEV+1)
      REAL*8,  INTENT(IN)    :: EPSL

!**** Local variables
      INTEGER                :: J, L

      REAL*8                 :: AA(IRUN,NLEV)
      REAL*8                 :: BB(IRUN,NLEV)
      REAL*8                 :: CC(IRUN,NLEV+1)

!**** TRBDIF.f begins here!
C
C     DEFINE MATRIX
C
      DO J = 1, IRUN
         CC(J,1) = 0d0
      ENDDO

      DO L = 1, NLEV
         DO J = 1, IRUN
            CC(J,L+1) = RHOKDZ(J,L) * FLXFAC(J,L+1)
         ENDDO
      ENDDO

      DO L = 1, NLEV
         DO J = 1, IRUN
            BB(J,L) = RHOKDZ(J,L) * FLXFAC(J,L)
            AA(J,L) = 1d0 + CC(J,L) + BB(J,L)
         ENDDO
      ENDDO

C     SOLVE MATRIX EQUATION FOR XX1
      CALL VTRI0(AA,BB,CC,XX1,XX1,NLEV,IRUN)

C
C     ELIMINATE UNDERFLOW
!**** Use F90 WHERE statement to perform test (bmy, 4/6/99)
!****      IF(LTYPE) THEN
!****         DO L = 1, NLEV
!****            DO J = 1, IRUN
!****               XX1(J,L) = amax1( EPSL, XX1(J,L) )
!****            ENDDO
!****         ENDDO
!****      ENDIF
      IF ( LTYPE ) THEN
         WHERE( XX1 <= EPSL ) XX1 = EPSL
      ENDIF

!**** Return to calling program
      RETURN
      END SUBROUTINE TRBDIF
