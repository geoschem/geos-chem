! $Id: lump.f,v 1.1 2003/06/30 20:26:04 bmy Exp $
      SUBROUTINE LUMP( NTRACER, XNUMOL, STT )
!
!******************************************************************************
!  Subroutine LUMP takes individual chemistry species and "lumps" them back 
!  into tracers after each SMVGEAR chemistry timestep. (bmy, 4/1/03)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) NTRACER (INTEGER) : Number of tracers
!  (2 ) XNUMOL  (REAL*8 ) : Array of molecules tracer / kg tracer 
!  (3 ) STT     (REAL*8 ) : Tracer concentrations [molec/cm3/box]
!
!  Arguments as Output:
!  ============================================================================
!  (3 ) STT     (REAL*8 ) : Tracer concentrations [kg/box]
!
!  NOTES:
!  (1 ) Updated comments, cosmetic changes (bmy, 4/1/03)
!******************************************************************************
!
      ! References to F90 modules
      USE COMODE_MOD,   ONLY : CSPEC,  JLOP,    VOLUME
      USE TRACERID_MOD, ONLY : IDTRMB, NMEMBER, CTRMB

      IMPLICIT NONE

#     include "CMN_SIZE"     ! Size parameters
#     include "comode.h"     ! SMVGEAR II arrays

      ! Arguments
      INTEGER, INTENT(IN)    :: NTRACER
      REAL*8,  INTENT(IN)    :: XNUMOL(NNPAR)
      REAL*8,  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,NNPAR)        

      ! Local variables
      INTEGER                :: I, J, L, N, JLOOP, KK, JJ
      REAL*8                 :: CONCTMP  

      !=================================================================
      ! LUMP begins here!
      !=================================================================
      DO N = 1, NTRACER
         
         ! Skip if not a valid tracer
         IF ( IDTRMB(N,1) == 0 ) CYCLE
       
         ! Loop over grid boxes
         DO L = 1, NPVERT
         DO J = 1, NLAT
         DO I = 1, NLONG

            ! 1-D SMVGEAR grid box index 
            JLOOP = JLOP(I,J,L)
            IF ( JLOOP == 0 ) CYCLE

            ! Compute tracer concentration [molec/cm3/box] by
            ! looping over all species belonging to this tracer
            CONCTMP = 0.d0
            DO KK = 1, NMEMBER(N)
               JJ =IDTRMB(N, KK)
               CONCTMP = CONCTMP + ( 1d0+CTRMB(N,KK) ) * CSPEC(JLOOP,JJ)
            ENDDO

            ! Save tracer concentrations back to STT
            STT(I,J,L,N) = CONCTMP

            ! Change STT from [molec/cm3/box] back to [kg/box]
            STT(I,J,L,N) = STT(I,J,L,N) * VOLUME(JLOOP) / XNUMOL(N)
         ENDDO
         ENDDO
         ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE LUMP

