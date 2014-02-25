!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: partition
!
! !DESCRIPTION: Subroutine PARTITION separates each GEOS-Chem tracer into 
!  its individual constituent chemistry species before each SMVGEAR or KPP 
!  chemistry timestep.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE PARTITION( NTRACER, XNUMOL, am_I_Root, State_Chm ) 
!
! !USES:
!
      USE COMODE_MOD,         ONLY : CSPEC,     JLOP,       VOLUME
      USE COMODE_MOD,         ONLY : JLOP_PREVIOUS
      USE ERROR_MOD,          ONLY : ALLOC_ERR, ERROR_STOP, SAFE_DIV
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE TRACERID_MOD,       ONLY : IDTRMB
      USE TRACERID_MOD,       ONLY : IDO3,      IDNO,       IDHNO2
      USE TRACERID_MOD,       ONLY : CTRMB,     NMEMBER
      USE CMN_SIZE_MOD
      USE COMODE_LOOP_MOD

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      ! # of tracers
      INTEGER, INTENT(IN)    :: NTRACER

      ! Conversion factor: molecules tracer / kg tracer 
      REAL*8,  INTENT(IN)    :: XNUMOL(NTRACER)

      ! Is this the root CPU?
      LOGICAL, INTENT(IN)    :: am_I_Root
!
! !INPUT/OUTPUT PARAMETERS: 
!
      ! Chemistry State object
      TYPE(ChmState), INTENT(INOUT) :: State_Chm
!
! !REMARKS:
!  Warning: 
!  ========
!  Partition was written assuming NOx tracer is before Ox tracer
!  in the tracer declaration in input.geos. If you want to change this
!  order, you need to adjust the code.
! 
! 
! !REVISION HISTORY: 
!  01 Apr 2003 - B. Field, R. Yantosca - Initial version, based on older code
!  (1 ) Now make CSAVE a local dynamic array.  Updated comments, cosmetic 
!        changes (bmy, 4/24/03)
!  (2 ) Add OpenMP parallelization commands (bmy, 8/1/03)
!  (3 ) Now dimension args XNUMOL, STT w/ NTRACER and not NNPAR (bmy, 7/20/04)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Resize CSAVE to save local memory, for SUN compiler. (bmy, 7/14/06)
!  (6 ) Now do safe division to eliminate FP errors (phs, bmy, 2/26/08)
!  (7 ) Now change error stop 30000 into a warning (phs, ccc, bmy, 1/7/09)
!  27 Jun 2011 - D. Henze, J. Koo - Fix to variable tropopause by Daven Henze.
!                                   When initializing CSAVE, search downward
!                                   in the column until we find a grid box
!                                   that was in the troposphere on the previous
!                                   timestep.
!  27 Jun 2011 - R. Yantosca - Added ProTeX headers
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  14 Mar 2013 - M. Payer    - Remove partitioning of NOx and Ox. Family tracers
!                              NOx and Ox have now been replaced with tracers
!                              NO, NO2, NO3, HNO3, and O3.
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J, L, N, JLOOP, IPL, JJ, KK
      INTEGER :: CSAVEID(IGAS)
      INTEGER :: CSAVEID_JJ(IGAS)
      INTEGER :: CS, IDNUM, AS 
      INTEGER :: LTEST
      INTEGER :: JLOOP_P 
      REAL*8  :: CONCTMP, CONCNOX, SUM, SUM1
      REAL*8  :: CSAVE( ITLOOP, NTRACER )
      REAL*8  :: QTEMP

      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL*8, POINTER :: STT(:,:,:,:)

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

            DO KK = 1, NMEMBER(N)
               IDNUM             = IDNUM + 1
               JJ                = IDTRMB(N,KK)
               CSAVEID(JJ)       = IDNUM
               CSAVEID_JJ(IDNUM) = JJ
            ENDDO

      ENDDO

      ! Loop over tracer members and boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, JLOOP, JLOOP_P, LTEST )
!$OMP+SCHEDULE( DYNAMIC )
      DO N = 1, IDNUM
      DO L = 1, NPVERT
      DO J = 1, NLAT
      DO I = 1, NLONG

         ! 1-D index for SMVGEAR/KPP arrays corresponding to 
         ! grid box (I,J,L) @ current chemistry timestep
         JLOOP = JLOP(I,J,L)

         ! Exit if grid box (I,J,L) is in the stratosphere 
         ! at the current chemistry timestep
         IF ( JLOOP == 0 ) CYCLE

         ! 1-D index for SMVGEAR/KPP arrays corresponding to 
         ! grid box (I,J,L) @ previous chemistry timestep
	 JLOOP_P = JLOP_PREVIOUS(I,J,L)

         ! If the grid box (I,J,L) was in the stratosphere 
         ! at the  previous timestep, then  ...
         IF ( JLOOP_P == 0 ) THEN

            ! Start with the current vertical level
            LTEST = L

            ! Search downward in the column (I,J,:) until we find a grid box 
            ! that was in the troposphere @ the previous chemistry timestep
            DO WHILE ( JLOOP_P == 0 )

               ! Try the next lowest grid box in the column
               LTEST = LTEST - 1

               ! Error: exit if we reach the surface without finding
               ! any boxes in the troposphere.  This should never happen!
               IF ( LTEST == 0 ) THEN
!$OMP CRITICAL
                  CALL ERROR_STOP( 'could no find trop box',
     &				   'partition.f' )
!$OMP END CRITICAL
               ENDIF

               ! Test if this next lowest grid box (I,J,LTEST) was in the
               ! troposphere at the previous chemistry timestep.  If so,
               ! this will terminate the search.  Save the index JLOOP_P,
               ! which corresponds to (I,J,LTEST), for use below.
               JLOOP_P = JLOP_PREVIOUS(I,J,LTEST)

            ENDDO
            
         ENDIF
		    
         ! Archive into the CSAVE array the concentration of each species
         ! in the CSPEC array at the previous chemistry timestep at index
         ! JLOOP_P.  JLOOP_P is guaranteed to be a tropopsheric box.  This
         ! will ensure proper partitioning.
         CSAVE(JLOOP,N) = CSPEC(JLOOP_P,CSAVEID_JJ(N))

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
      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      DO N = 1, NTRACER

         ! Skip if it's not a valid tracer
         IF ( IDTRMB(N,1) == 0 ) CYCLE

         !### Debug
         !IF( am_I_Root ) WRITE(6,*) 'IN PARTITION N= ', N

         ! Loop over grid boxes
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, JLOOP, CONCTMP, SUM, KK, JJ, SUM1, CONCNOX )
!$OMP+PRIVATE( QTEMP )
!$OMP+SCHEDULE( DYNAMIC )
         DO L = 1, NPVERT
         DO J = 1, NLAT
         DO I = 1, NLONG

            ! 1-D SMVGEAR grid box index
            JLOOP = JLOP(I,J,L)
            IF ( JLOOP == 0 ) CYCLE

            ! Convert tracer concentration from [kg/box] to [molec/cm3/box]
            STT(I,J,L,N) = STT(I,J,L,N) / VOLUME(JLOOP) * XNUMOL(N)

            ! Store concentration of tracer N at grid box (I,J,L) in CONCTMP 
            CONCTMP = STT(I,J,L,N)

            !===========================================================
            ! First, find sum of starting concentrations
            !===========================================================

            SUM = 0.d0

            DO KK = 1, NMEMBER(N)
               JJ = IDTRMB(N, KK)

               ! Error check
               IF ( JJ == 0 ) THEN
!$OMP CRITICAL
                  IF ( am_I_Root ) THEN
                     PRINT *,JJ,JLOOP,N,KK,IDTRMB(N, KK)
                  ENDIF
!$OMP END CRITICAL
               ENDIF

               SUM = SUM + CSAVE(JLOOP,CSAVEID(JJ)) * (CTRMB(N,KK)+1)
            ENDDO

            !===========================================================
            ! Now perform the partitioning
            !===========================================================

            ! Loop over # of member species in this tracer
            DO KK = 1, NMEMBER(N)

               ! Index of member species for CSPEC
               JJ = IDTRMB(N, KK)

               ! QTEMP is the fraction of the given member species KK
               ! in the tracer N.  The value QTEMP*CONCTMP is the 
               ! concentration of the member species itself, and that
               ! needs to be saved into CSPEC.
               !
               ! In the partitioning, now be sure to perform a safe
               ! floating point division of CSAVE/SUM.  Return the value
               ! 1/NMEMBER(N) if the division can't be done, i.e. do a
               ! uniform paritioning among all member species of the
               ! given tracer. (phs, bmy, 2/26/08)
               QTEMP = SAFE_DIV( CSAVE(JLOOP,CSAVEID(JJ)),
     &                           SUM, 1d0/NMEMBER(N) )
                     
               ! Store the concentration of member species KK
               ! into the CSPEC array.  Do not allow underflow!
               CSPEC(JLOOP,JJ) = MAX( QTEMP*CONCTMP, SMAL2 )
            ENDDO

         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO
      ENDDO

      ! Free pointer
      NULLIFY( STT )

      END SUBROUTINE PARTITION
!EOC
