! $Id: diagpl.f,v 1.1 2003/06/30 20:26:05 bmy Exp $
      SUBROUTINE DIAGPL
!
!*****************************************************************************
!  Subroutine DIAGPL records production and loss for specified families 
!  for the ND65 diagnostic (bey, bmy, 3/16/00, 6/25/02)
!
!  NOTES:
!  (1 ) Original code is from Loretta Mickley and Isabelle Bey
!  (2 ) Now use allocatable array for ND65 diagnostic.
!       Also made cosmetic changes, updated comments (bmy, 3/16/00)
!  (3 ) Removed obsolete code from 3/16/00 (bmy, 8/31/00)
!  (4 ) Replaced all instances of IM with IIPAR and JM with JJPAR, in order
!        to prevent namespace confusion for the new TPCORE.  Also removed
!        obsolete, commented-out code. (bmy, 6/25/02)
!*****************************************************************************
!
      ! References to F90 modules
      USE DIAG_MOD, ONLY : AD65, FAMPL

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! NFAM
#     include "CMN_DIAG"  ! Diagnostic arrays & switches

      ! Local variables
      INTEGER :: I, J, L, N  

      !=================================================================
      ! DIAGPL begins here!
      !
      ! If ND65 is turned on, then archive P-L for specified families
      ! and store in the AD65 array. 
      !=================================================================
      IF ( ND65 > 0 ) THEN
         DO N = 1, NFAM
         DO L = 1, LD65
         DO J = 1, JJPAR
         DO I = 1, IIPAR
            AD65(I,J,L,N) = AD65(I,J,L,N) + FAMPL(I,J,L,N)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      ! Return to calling program
      END SUBROUTINE DIAGPL
