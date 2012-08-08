!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: lump
!
! !DESCRIPTION: Subroutine LUMP takes individual chemistry species and 
!  "lumps" them back into tracers after each SMVGEAR chemistry timestep.

!\\
!\\
! !INTERFACE:
!
      SUBROUTINE LUMP( NTRACER, XNUMOL, STT )
!
! !USES:
!
      USE COMODE_MOD,   ONLY : CSPEC,  JLOP,    VOLUME
      USE TRACERID_MOD, ONLY : IDTRMB, NMEMBER, CTRMB
      USE CMN_SIZE_MOD                            ! Size parameters
      USE COMODE_LOOP_MOD                         ! SMVGEAR II arrays

      IMPLICIT NONE
#     include "define.h"
!
! !INPUT PARAMETERS: 
!
      INTEGER, INTENT(IN)    :: NTRACER           ! # of tracers
      REAL*8,  INTENT(IN)    :: XNUMOL(NTRACER)   ! Molec tracer / kg tracer
!
! !INPUT/OUTPUT PARAMETERS: 
!
      ! Tracer concentrations [kg/box]
      REAL*8,  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,NTRACER)
! 
! !REVISION HISTORY: 
!  01 Apr 2003 - R. Yantosca - Initial version
!  (1 ) Updated comments, cosmetic changes (bmy, 4/1/03)
!  (2 ) Added OpenMP parallelization commands (bmy, 8/1/03)
!  (3 ) Now dimension args XNUMOL, STT w/ NTRACER and not NNPAR (bmy, 7/20/04)
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  30 Jul 2012 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES
!
      INTEGER :: I, J, L, N, JLOOP, KK, JJ
      REAL*8  :: CONCTMP  

      !=================================================================
      ! LUMP begins here!
      !=================================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, JLOOP, CONCTMP, KK, JJ )
!$OMP+SCHEDULE( DYNAMIC )
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
               JJ = IDTRMB(N, KK)
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
!$OMP END PARALLEL DO

      END SUBROUTINE LUMP
!EOC

