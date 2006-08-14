! $Id: partition.f,v 1.8 2006/08/14 17:58:12 bmy Exp $
      SUBROUTINE PARTITION( NTRACER, STT, XNUMOL ) 
!
!******************************************************************************
!  Subroutine PARTITION separates GEOS-CHEM tracers into its individual
!  constituent chemistry species before each SMVGEAR chemistry timestep.
!  (bdf, bmy, 4/1/03, 7/14/06)
! 
!  Arguments as Input:
!  ============================================================================
!  (1 ) NTRACER (INTEGER) : Number of tracers
!  (2 ) STT     (REAL*8 ) : Tracer concentrations [kg/box]
!  (3 ) XNUMOL  (REAL*8 ) : Array of molecules tracer / kg tracer 
!
!  Arguments as Output:
!  ============================================================================
!  (1 ) STT     (REAL*8 ) : Updated tracer concentrations [molec/cm3/box]
!
!  NOTES:
!  (1 ) Now make CSAVE a local dynamic array.  Updated comments, cosmetic 
!        changes (bmy, 4/24/03)
!  (2 ) Add OpenMP parallelization commands (bmy, 8/1/03)
!  (3 ) Now dimension args XNUMOL, STT w/ NTRACER and not NNPAR (bmy, 7/20/04)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Resize CSAVE to save local memory, for SUN compiler. (bmy, 7/14/06)
!******************************************************************************
!
      ! References to F90 modules 
      USE COMODE_MOD,   ONLY : CSPEC,     JLOP,      VOLUME
      USE ERROR_MOD,    ONLY : ALLOC_ERR, ERROR_STOP
      USE TRACERID_MOD, ONLY : IDTOX,     IDTNOX,    IDTRMB
      USE TRACERID_MOD, ONLY : IDO3,      IDNO,      IDHNO2
      USE TRACERID_MOD, ONLY : CTRMB,     NMEMBER

      IMPLICIT NONE

#     include "CMN_SIZE"
#     include "comode.h"

      ! Arguments
      INTEGER, INTENT(IN)    :: NTRACER
      REAL*8,  INTENT(INOUT) :: STT(IIPAR,JJPAR,LLPAR,NTRACER)
      REAL*8,  INTENT(IN)    :: XNUMOL(NTRACER)

      ! Local variables
      INTEGER                :: I, J, L, N, JLOOP, IPL, JJ, KK
      INTEGER                :: CSAVEID(IGAS)
      INTEGER                :: CSAVEID_JJ(IGAS)
      INTEGER                :: CS, IDNUM, AS  
      REAL*8                 :: CONCTMP, CONCNOX, SUM, SUM1
      REAL*8                 :: CSAVE( ITLOOP, NTRACER )
      
      !=================================================================
      ! PARTITION begins here!
      !
      ! Copy values of CSPEC that need to be saved  (bdf, 3/30/99)
      !=================================================================

      ! Initialize
      IDNUM         = 0
      CSAVEID(:)    = 0
      CSAVEID_JJ(:) = 0

      ! Loop over tracers
      DO N = 1, NTRACER

         ! Skip if this is not a valid tracer
         IF ( IDTRMB(N,1) == 0 ) CYCLE

         ! Handle all other tracers except Ox 
         IF ( N /= IDTOX ) THEN
            DO KK = 1, NMEMBER(N)
               IDNUM             = IDNUM + 1
               JJ                = IDTRMB(N,KK)
               CSAVEID(JJ)       = IDNUM
               CSAVEID_JJ(IDNUM) = JJ
            ENDDO

         ! Handle Ox
         ELSE IF ( IDTOX /= 0 ) THEN
            JJ                = IDTRMB(N,1)
            IDNUM             = IDNUM + 1
            CSAVEID(JJ)       = IDNUM
            CSAVEID_JJ(IDNUM) = JJ
         ENDIF
      ENDDO

      ! Loop over tracer members and boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, JLOOP )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1, IDNUM
      DO L = 1, NPVERT
      DO J = 1, NLAT
      DO I = 1, NLONG

         ! 1-D SMVGEAR grid box index
         JLOOP = JLOP(I,J,L)
         IF ( JLOOP == 0 ) CYCLE

         ! Store into CSAVE
         CSAVE(JLOOP,N) = CSPEC(JLOOP,CSAVEID_JJ(N))
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Split each tracer up into its components (if any)
      ! Family tracers are partitioned among members according to 
      ! initial ratios. In tracer sequence, OX must be after NOX, 
      ! otherwise, adjust the code
      !=================================================================
      DO N = 1, NTRACER

         ! Skip if it's not a valid tracer
         IF ( IDTRMB(N,1) == 0 ) CYCLE

         !### Debug
         !WRITE(6,*) 'IN PARTITION N= ', N

         ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, JLOOP, CONCTMP, SUM, KK, JJ, SUM1, CONCNOX )
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, NPVERT
         DO J = 1, NLAT
         DO I = 1, NLONG

            ! 1-D SMVGEAR grid box index
            JLOOP = JLOP(I,J,L)
            IF ( JLOOP == 0 ) CYCLE

            ! Convert tracer concentration from [kg/box] to [molec/cm3/box]
            STT(I,J,L,N) = STT(I,J,L,N) / VOLUME(JLOOP) * XNUMOL(N)

            ! Store concentration in CONCTMP 
            CONCTMP = STT(I,J,L,N)

            !===========================================================
            ! First, find sum of starting concentrations
            !===========================================================

            !------------------------
            ! All tracers except Ox
            !------------------------
            IF ( N /= IDTOX ) THEN
               SUM = 0.d0

               DO KK = 1, NMEMBER(N)
                  JJ = IDTRMB(N, KK)

                  ! Error check
                  IF ( JJ == 0 ) THEN
                     PRINT *,JJ,JLOOP,N,KK,IDTRMB(N, KK)
                  ENDIF

                  SUM = SUM + CSAVE(JLOOP,CSAVEID(JJ)) * (CTRMB(N,KK)+1)
               ENDDO

            !------------------------
            ! Ox
            !------------------------
            ELSE IF ( IDTOX /= 0 ) THEN
               JJ   = IDTRMB(N,1)
               SUM  = CSAVE(JLOOP,CSAVEID(JJ)) * (CTRMB(N,1)+1)
               SUM1 = 0.d0


               ! SUM  = sum of starting values for all Ox species (incl. O3)
               ! SUM1 = sum of new values for all Ox species except O3,
               ! based on NOx partitioning
               DO KK = 2, NMEMBER(N)
                  JJ   = IDTRMB(N,KK)
                  SUM  = SUM + CSAVE(JLOOP,CSAVEID(JJ))*(CTRMB(N,KK)+1)
                  SUM1 = SUM1+ CSPEC(JLOOP,JJ) * (CTRMB(N,KK)+1)
               ENDDO

            ENDIF

            !===========================================================
            ! Now perform the partitioning
            !===========================================================

            !----------------------------------------
            ! All tracers except Ox
            !----------------------------------------
            IF ( N /= IDTOX ) THEN
               DO KK = 1, NMEMBER(N)
                  JJ = IDTRMB(N, KK)
                  CSPEC(JLOOP,JJ) = ( CSAVE(JLOOP,CSAVEID(JJ)) / SUM ) *
     &                              CONCTMP
               ENDDO

            !----------------------------------------
            ! For Ox, take O3 = Ox - SUM(NO2+NO3*2)
            !----------------------------------------
            ELSE IF ( IDTOX /= 0 .AND. IDTNOX /= 0 ) THEN

               ! Find Ox in CSPEC
               JJ              = IDO3
               CSPEC(JLOOP,JJ) = CONCTMP - SUM1

               ! If Ox partitioning is OK, then skip to next tracer
               IF ( CSPEC(JLOOP,JJ) > 0.0d0 ) GOTO 220

               ! Ox partitioning failed, we are getting a negative ozone 
               ! concentration.  Instead, try partitioning Ox before NOx
               DO KK = 1, NMEMBER(N)
                  JJ = IDTRMB(N, KK)
                  CSPEC(JLOOP,JJ) = 
     &                 ( CSAVE(JLOOP,CSAVEID(JJ)) / SUM ) * CONCTMP
               ENDDO

               ! then partition NO+HNO2 
               ! (the only NOx species not contained in Ox)
               ! SUM  = sum of starting values for NO and HNO2
               ! SUM1 = sum of new values for all NOx species except 
               ! NO and HNO2, based on Ox partitioning
               SUM  = 0.d0
               SUM1 = 0.d0

               DO KK = 1, NMEMBER(IDTNOX)
                  JJ = IDTRMB(IDTNOX, KK)

                  IF ( JJ == IDNO .OR. JJ == IDHNO2 ) THEN
                     SUM = SUM + CSAVE(JLOOP,CSAVEID(JJ)) *
     &                           (CTRMB(IDTNOX,KK)+1)
                  ELSE
                     SUM1 = SUM1 + CSPEC(JLOOP,JJ)*(CTRMB(IDTNOX,KK)+1)
                  ENDIF
               ENDDO

               ! Get NOx concentration from STT
               CONCNOX = STT(I,J,L,IDTNOX)

               ! Stop w/ error
               IF ( CONCNOX - SUM1 < 0.d0 ) THEN
                  CALL ERROR_STOP( 'STOP 30000', 'partition.f' )
               ENDIF

               DO KK = 1, NMEMBER(IDTNOX)
                  JJ = IDTRMB(IDTNOX,KK)

                  IF ( JJ == IDNO .OR. JJ == IDHNO2 ) THEN
                     CSPEC(JLOOP,JJ) = 
     &               ( CSAVE(JLOOP,CSAVEID(JJ)) / SUM ) * (CONCNOX-SUM1)
                  ENDIF
               ENDDO

               !========================================================
               ! Ox partitioning is OK
               !========================================================
 220           CONTINUE
            ENDIF
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDDO

      ! Return to calling program
      END SUBROUTINE PARTITION
