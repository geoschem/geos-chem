! $Id: plsave.f,v 1.2 2003/08/06 15:30:52 bmy Exp $
      SUBROUTINE PLSAVE  
!
!*****************************************************************************
!  Subroutine PLSAVE saves info on production and loss of families 
!  specified in prodloss.dat (bey, bmy, 3/16/00, 8/1/03)
!
!  NOTES:
!  (1 ) Original code is from Loretta Mickley and Isabelle Bey
!  (2 ) Now use allocatable array for ND65 diagnostic.
!        Also made cosmetic changes, updated comments (bmy, 3/16/00)
!  (3 ) Now reference the CSPEC array from "comode_mod.f" instead of from
!        common block header "comode.h".  Also make sure that FAMPL and
!        CSPEC arrays have been allocated. (bmy, 7/11/00)
!  (4 ) Also reference JLOP from "comode_mod.f" (bmy, 10/19/00)
!  (5 ) Removed obsolete code from 10/19/00 (bmy, 12/21/00)
!  (6 ) Updated comments, Cosmetic changes (bmy, 5/13/03)
!  (7 ) Added OpenMP parallelization commands (bmy, 8/1/03)
!*****************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD, ONLY : CSPEC, JLOP
      USE DIAG_MOD,   ONLY : FAMPL

      IMPLICIT NONE

#     include "CMN_SIZE"  ! Size parameters
#     include "comode.h"  ! SMVGEAR II arrays

      ! Local variables 
      INTEGER :: I, J, L, JLOOP, N

      !=================================================================
      ! PLSAVE begins here!
      !
      ! If ND65 is turned on, then archive P-L for specified families 
      ! and store in the AD65 array.  
      !
      ! Make sure that memory has already been allocated to arrays
      ! FAMPL, JLOP, and CSPEC.
      !=================================================================
      IF ( ALLOCATED( FAMPL ) .and. 
     &     ALLOCATED( CSPEC ) .and. 
     &     ALLOCATED( JLOP  ) ) THEN 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, JLOOP )
!$OMP+SCHEDULE( DYNAMIC ) 
         DO N = 1, NFAMILIES
         DO L = 1, NVERT
         DO J = 1, NLAT
         DO I = 1, NLONG
         
            ! JLOOP is the 1-D grid box index for SMVGEAR arrays
            JLOOP = JLOP(I,J,L)
            IF ( JLOOP == 0 ) CYCLE

            ! FAMPL stores the "fake" prodloss family, which have been
            ! appended to the SMVGEAR species list. [molecules/cm3/s]
            FAMPL(I,J,L,N)       = CSPEC(JLOOP,IFAM(N)) / CHEMINTV

            ! Zero each "fake" ND65 prod/loss family for next iteration
            CSPEC(JLOOP,IFAM(N)) = 0.0d0
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDIF

      ! Return to calling program
      END SUBROUTINE PLSAVE
