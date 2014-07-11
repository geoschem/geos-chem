!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: setemis
!
! !DESCRIPTION: Subroutine SETEMIS places emissions from the 3-D 
! emission fields in State\_Chm%Trac\_Tend into the REMIS array (used
! for SMVGEAR II chemistry). Emissions are converted from kg/m2/s to
! molec/cm3/s.
!\\
!\\
! If the full-mixing PBL scheme is used, the entire emission arrays
! of State\_Chm%Trac\_Tend are passed to REMIS. For the non-local PBL
! scheme, only the emissions above the PBL are passed since all
! emissions w/in the PBL are already considered in the non-local PBL
! mixing routine (vdiff\_mod.F90). 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SETEMIS( am_I_Root, Input_Opt, 
     &                    State_Met, State_Chm, RC )
!
! !USES:
!
      USE GIGC_ErrCode_Mod
      USE COMODE_MOD,         ONLY : JLOP, REMIS
      USE CMN_SIZE_MOD,       ONLY : IIPAR, JJPAR, LLPAR
      USE COMODE_LOOP_MOD,    ONLY : NTTLOOP, NCS, IDEMS, NEMIS
      USE TRACERID_MOD,       ONLY : CTRMB, IDEMIS
      USE CHEMGRID_MOD,       ONLY : ITS_IN_THE_NOCHEMGRID
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE PBL_MIX_MOD,        ONLY : GET_PBL_TOP_L
!
! !ARGUMENTS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root
      TYPE(OptInput), INTENT(IN)  :: Input_Opt  ! Input options
      TYPE(MetState), INTENT(IN)  :: State_Met
      TYPE(ChmState), INTENT(IN)  :: State_Chm
      INTEGER,        INTENT(OUT) :: RC

! !REVISION HISTORY:
!  13 Apr 2012 - C. Keller - Initial Version
!  31 Jan 2013 - C. Keller - Now pass chemistry and met state as argument
!  22 Aug 2013 - C. Keller - Now pass Trac_Tend2 to REMIS
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: I, J,  JLOOP, JLOOP1, LTROP
      INTEGER :: L, L1, LL, N, NN,  NBB, NBF
      INTEGER :: N_TRACERS
      REAL*8  :: COEF1, COEF2
      LOGICAL :: LNLPBL

      !=================================================================
      ! SETEMIS begins here!
      !=================================================================

      ! Assume success
      RC = GIGC_SUCCESS

      ! Get variables from Input object
      N_TRACERS = Input_Opt%N_TRACERS
      LNLPBL    = Input_Opt%LNLPBL

      ! Reset 
      REMIS = 0d0

      ! NOTE: If NCS=0 and we are NOT doing a full-chemistry run, then
      ! setting REMIS=0 can lead to a seg-fault error.  Therefore, you should
      ! add an error check in the calling routine to prevent SETEMIS from
      ! being called if we are not doing full-chemistry. (bmy, 12/21/10)
      IF ( NCS == 0 ) THEN
         RETURN
      ENDIF

      ! Loop over all emission species
      DO N = 1, NEMIS(NCS)

         ! Get CTM tracer number corresponding to emission species N
         NN = IDEMS(N)
         IF ( NN == 0 ) CYCLE

         ! Skip to next tracer if IDEMIS(NN) is not defined in
         ! order to avoid array-out-of-bounds errors (bmy, 8/6/04)
         IF ( IDEMIS(NN) <= 0 ) THEN
            IF ( am_I_Root ) THEN
               PRINT*, 'Tracer ', NN, ' is not an emitted species!'
            ENDIF 
            CYCLE
         ENDIF

         ! COEF1 = molecules of emission species / molecules of tracer
         COEF1 = 1.0 + CTRMB(NN, IDEMIS(NN))

         ! COEF2 = conversion factor from kg/m2/s to molec/cm2/s
         COEF2 = Input_Opt%XNUMOL(NN) / 1.0d4

         ! Loop over all tropospheric grid boxes
         ! --> Stratospheric emissions are ignored!
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! For non-local PBL, only consider emissions above
            ! PBL. Emissions within PBL were used within the
            ! non-local PBL scheme in vdiff_mod.F!
            IF ( LNLPBL ) THEN
               L1 = MAX( 1, FLOOR( GET_PBL_TOP_L(I,J) ) ) + 1
            ELSE
               L1 = 1
            ENDIF

            ! Loop over all vertical levels
            DO L = L1, LLPAR

               ! Proceed to next grid box if we surpassed the chemistry grid 
               IF ( ITS_IN_THE_NOCHEMGRID(I,J,L,State_Met) ) EXIT

               ! Write emissions to REMIS array. Convert from 
               ! kg(tracer)/m2/s to molec(species)/cm3/s.
               JLOOP = JLOP(I,J,L)
               IF ( JLOOP /= 0 ) THEN
                  REMIS(JLOOP,N) = State_Chm%Trac_Tend(I,J,L,NN)
     &                           / COEF1 * COEF2 
     &                           / ( State_Met%BXHEIGHT(I,J,L) * 1d2 )
               ENDIF
            ENDDO !L
         ENDDO
         ENDDO

         ! testing only
!         write(*,*) 'total REMIS of tracer ', NN , ': ', SUM(REMIS(:,N))

      ENDDO ! N

      END SUBROUTINE SETEMIS
!EOC
