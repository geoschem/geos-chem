C $Id: sumup.f,v 1.1 2003/06/30 20:26:08 bmy Exp $
C*********************************************************
      FUNCTION SUMUP(I1,I2,J1,J2,L1,L2,K1,K2,UPDOWN)
C*********************************************************
C
      IMPLICIT NONE

#     include "CMN_SIZE"       ! Size parameters
#     include "CMN"            ! LPAUSE
#     include "CMN_CO_BUDGET"  ! TCO
C*****************************************************************************
C Store the sources/sinks of CO in TCO in total molecules
C           ( 1) = initial burden
C           ( 2) = final burden
C  SINKS
C           ( 3) = CO sink by OH
C  SOURCES
C           ( 4) = CH4 oxidation
C           ( 5) = isoprene oxidation
C           ( 6) = biomass burning
C           ( 7) = wood burning
C           ( 8) = fossil fuels
C           ( 9) = monoterpene oxidation
C
C           (10) = interhemispheric exchange (+ = northward)
C           (11) = vacant
C           (12) = vacant
C*****************************************************************************
C Called from SR CO_budget.
C
C Levels            1 <= L <= LPAUSE(I,J) - 1 are tropospheric
C         LPAUSE(I,J) <= L <= LLPAR           are stratospheric (bmy, 4/17/00)
C
      REAL*8 SUMUP
      INTEGER I1,I2,J1,J2,L1,L2,K1,K2,I,J,K,L,UPDOWN
      ! Variable to hold the minimum value of LPAUSE (bmy, 4/18/00)
      INTEGER  LPAUSE_MIN,LPAUSE_MAX

      ! Compute the minimum value of LPAUSE once for use in
      ! the DO-loops below (bmy, 4/18/00)
      LPAUSE_MIN = MINVAL( LPAUSE )
      LPAUSE_MAX = MAXVAL( LPAUSE )
      print*,'LPAUSE MIN/MAX=',LPAUSE_MIN,LPAUSE_MAX  
      print*,'L1,L2=',L1,L2

      ! Start L-loop from the lowest stratospheric level (bmy, 4/18/00)

      SUMUP=0.

      IF(UPDOWN.EQ.1) THEN
         
         DO K=K1,K2
         DO L=L1,LPAUSE_MAX
         DO J=J1,J2
         DO I=I1,I2

         ! Only process tropospheric boxes (bmy, 4/17/00)
            IF ( L < LPAUSE(I,J) ) THEN 
               SUMUP=SUMUP+TCO(I,J,L,K)
            ENDIF

         ENDDO
         ENDDO
         ENDDO
         ENDDO

      ELSE

         DO K=K1,K2
         DO L=LPAUSE_MIN,L2
         DO J=J1,J2
         DO I=I1,I2

         ! Only process stratospheric boxes (bmy, 4/17/00)
            IF ( L >= LPAUSE(I,J) ) THEN 
               SUMUP=SUMUP+TCO(I,J,L,K)
            ENDIF
            
         ENDDO
         ENDDO
         ENDDO
         ENDDO
      ENDIF

      RETURN
      END
