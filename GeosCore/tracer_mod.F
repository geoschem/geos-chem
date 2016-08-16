!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: tracer_mod
!
! !DESCRIPTION: Module TRACER\_MOD contains GEOS-CHEM tracer array STT plus
!  various other related quantities.  TRACER\_MOD also contains inquiry
!  functions that can be used to determine the type of GEOS-CHEM simulation.
!\\
!\\
! !INTERFACE:
!
      MODULE TRACER_MOD
!
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS: 
!
      PUBLIC :: CHECK_STT
      PUBLIC :: CHECK_STT_NESTED
!
! !REMARKS:
!  Most of the functionality of this module has been superseded by
!  the Input_Opt derived type object.
!
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca -  Initial version
!  (1 ) Added function GET_SIM_NAME (bmy, 5/3/05)
!  (2 ) Removed ITS_A_COPARAM_SIM; the CO-OH param is obsolete (bmy, 6/24/05)
!  (3 ) Added ITS_A_CO2_SIM (pns, bmy, 7/25/05)
!  (4 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (5 ) Now added XNUMOL, XNUMOLAIR as module variables (bmy, 10/25/05)
!  (6 ) Added public routine ITS_A_H2HD_SIM (phs, 9/18/07)
!  (7 ) Added public routine ITS_A_POPS_SIM (eck, 9/20/10)
!  16 Feb 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!  04 Apr 2013 - R. Yantosca - Removed STT (now in State_Chm)
!  21 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  24 Mar 2015 - E. Lundgren - Move param XNUMOLAIR to CMN_GCTM_mod.F
!  24 Mar 2015 - E. Lundgren - Remove shadow copy of Input_Opt%N_TRACERS
!  24 Mar 2015 - E. Lundgren - Put init_tracer within APM ifdef block
!  24 Mar 2015 - E. Lundgren - Change tracer units from kg to kg/kg
!  23 Jun 2016 - R. Yantosca - Remove references to APM code; it is no longer
!                              compatible with the FlexChem implementation
!EOP
!------------------------------------------------------------------------------
!BOC
!
      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_stt
!
! !DESCRIPTION: Subroutine CHECK\_STT checks the STT tracer array for negative
!  values, NaN values, or Infinity values.  If any of these are found, the code
!  will stop with an error message.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_STT( State_Chm, Input_Opt, LOCATION )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE ERROR_MOD,          ONLY : IT_IS_FINITE
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN) :: LOCATION
      TYPE(ChmState),   INTENT(IN) :: State_Chm   ! Chemistry State object
      TYPE(OptInput),   INTENT(IN) :: Input_Opt   ! Input Options object
! 
! !REVISION HISTORY:
!  (1 ) CHECK_STT uses the interfaces defined above -- these will do the
!        proper error checking for either SGI or DEC/Compaq platforms.
!        (bmy, 3/8/01)
!  (2 ) Now call GEOS_CHEM_STOP to shutdown safely.  Now use logicals LNAN,
!        LNEG, LINF to flag if we have error conditions, and then stop the
!        run outside of the parallel DO-loop. (bmy, 11/27/02)
!  (3 ) Bug fix in FORMAT statement: replace missing commas (bmy, 3/23/03)
!  (4 ) Moved from "error_mod.f" to "tracer_mod.f" (bmy, 7/15/04)
!  (5 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!  25 Mar 2013 - M. Payer    - Now pass State_Chm object via the arg list
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL           :: LNEG, LNAN, LINF
      INTEGER           :: I,    J,    L,   N,  NA,  nAdvect
      
      ! Pointers
      REAL(fp), POINTER :: Spc(:,:,:,:)

      !=================================================================
      ! CHECK_STT begins here!
      !=================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ REMOVE TRACERS MODIFICATION (bmy, 6/30/16)
!@@@ Need to force State_Chm%Species = State_Chm%Tracers during development
!@@@ This can be removed later once State_Chm%Tracers is removed everywhere
!@@@
      ! Number of advected species
      nAdvect = State_Chm%nAdvect

      ! Force State_Chm%SPECIES = State_Chm%TRACERS for testing  
      DO NA = 1, nAdvect
         N                          = State_Chm%Map_Advect(NA)
         State_Chm%Species(:,:,:,N) = State_Chm%Tracers(:,:,:,N)
      ENDDO
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      ! Point to GEOS-Chem species [kg/kg]
      Spc => State_Chm%Species

      ! Initialize
      LNEG = .FALSE.
      LNAN = .FALSE.
      LINF = .FALSE.

      ! Loop over the advected species
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, NA )
      DO NA = 1, nAdvect
         
         ! Advected species ID
         N = State_Chm%Map_Advect(NA)

         ! Loop over grid boxes
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            !---------------------------
            ! Check for Negatives
            !---------------------------
            Spc(I,J,L,N) = ABS(Spc(I,J,L,N))
            IF ( Spc(I,J,L,N) < 0e+0_fp ) THEN 
!$OMP CRITICAL
               LNEG = .TRUE.
               WRITE( 6, 100 ) I, J, L, N, Spc(I,J,L,N)
               PRINT*, Spc(I,J,L,N)
!$OMP END CRITICAL

            !---------------------------
            ! Check for NaN's
            !---------------------------
            ELSE IF ( IT_IS_NAN( Spc(I,J,L,N) ) ) THEN
!$OMP CRITICAL
               LNAN = .TRUE.
               WRITE( 6, 100 ) I, J, L, N, Spc(I,J,L,N)
!$OMP END CRITICAL

            !----------------------------
            ! Check Spc's for Infinities
            !----------------------------
            ELSE IF ( .not. IT_IS_FINITE( Spc(I,J,L,N) ) ) THEN
!$OMP CRITICAL
               LINF = .TRUE.
               WRITE( 6, 100 ) I, J, L, N, Spc(I,J,L,N)
!$OMP END CRITICAL            

            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      Spc => NULL()

      !=================================================================
      ! Stop the run if any of LNEG, LNAN, LINF is true
      !=================================================================
      IF ( LNEG .or. LNAN .or. LINF ) THEN
         WRITE( 6, 120 ) TRIM( LOCATION ), LNEG, LNAN, LINF
         CALL GEOS_CHEM_STOP
      ENDIF

      !=================================================================
      ! FORMAT statements
      !=================================================================
 100  FORMAT( 'CHECK_STT: STT(',i3,',',i3,',',i3,',',i3,') = ', f13.6 )
 120  FORMAT( 'CHECK_STT: STOP at ', a , 3L2 )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ REMOVE TRACERS MODIFICATION (bmy, 6/30/16)
!@@@ Need to restore State_Chm%TRACERS = State_Chm%SPECIES for testing 
!@@@
      DO NA = 1, nAdvect
         N                          = State_Chm%Map_Advect(NA)
         State_Chm%Tracers(:,:,:,N) = State_Chm%Species(:,:,:,N)
      ENDDO
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      END SUBROUTINE CHECK_STT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_stt_nested
!
! !DESCRIPTION: Subroutine CHECK\_STT\_NESTED checks the STT tracer array for
!  negative values, NaN values, or Infinity values.  If any of these are found,
!  the STT array will be set to a specified value.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_STT_NESTED( State_Chm, Input_Opt, LOCATION )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE ERROR_MOD,          ONLY : IT_IS_FINITE
      USE Input_Opt_Mod,      ONLY : OptInput
      USE State_Chm_Mod,      ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN)    :: LOCATION
      TYPE(OptInput),   INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  05 Mar 2012 - M. Payer    - Initial version based on CHECK_STT and updates
!                              for nested grid by Yuxuan Wang.
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!  26 Oct 2015 - M. Sulprizio- Rename from CHECK_STT_05x0666 to CHECK_STT_NESTED
!                              for use with all nested grids
!  30 Jun 2016 - R. Yantosca - Remove instances of STT.  Now get the advected
!                              species ID from State_Chm%Map_Advect.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER           :: I, J, L, N, NA, nAdvect
      
      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER :: Spc(:,:,:,:)

      !=================================================================
      ! CHECK_STT_NESTED begins here!
      !=================================================================
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ REMOVE TRACERS MODIFICATION (bmy, 6/30/16)
!@@@ Need to force State_Chm%Species = State_Chm%Tracers during development
!@@@ This can be removed later once State_Chm%Tracers is removed everywhere
!@@@
      ! Number of advected species
      nAdvect = State_Chm%nAdvect

      ! Force State_Chm%SPECIES = State_Chm%TRACERS for testing  
      DO NA = 1, nAdvect
         N                          = State_Chm%Map_Advect(NA)
         State_Chm%Species(:,:,:,N) = State_Chm%Tracers(:,:,:,N)
      ENDDO
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      ! Point to chemical species array [kg/kg]
      Spc => State_Chm%Species

      ! Loop over advected species
      ! Loop over grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, NA )
      DO NA = 1, nAdvect

         ! Advected species ID
         N = State_Chm%Map_Advect(NA)

         ! Loop over grid boxes
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            !---------------------------
            ! Check for Negatives
            !---------------------------
            IF ( Spc(I,J,L,N) < 0e+0_fp ) THEN 
!$OMP CRITICAL
               WRITE( 6, 100 ) I, J, L, N, Spc(I,J,L,N)
               PRINT*, 'Neg Spc ' // TRIM( LOCATION ) //
     &                 '. SET Spc TO BE ZERO.'
               Spc(I,J,L,N) = 0e+0_fp
!$OMP END CRITICAL

            !---------------------------
            ! Check for NaN's
            !---------------------------
            ELSE IF ( IT_IS_NAN( Spc(I,J,L,N) ) ) THEN 
!$OMP CRITICAL
               WRITE( 6, 100 ) I, J, L, N, Spc(I,J,L,N)
               PRINT*, 'NaN Spc ' // TRIM( LOCATION ) //
     &              '. SET Spc TO BE LOWER LEVEL.'
               Spc(I,J,L,N) = Spc(I,J,L-1,N)
!$OMP END CRITICAL

            !----------------------------
            ! Check Spc's for Infinities
            !----------------------------
            ELSE IF ( .not. IT_IS_FINITE( Spc(I,J,L,N) ) ) THEN
!$OMP CRITICAL
               WRITE( 6, 100 ) I, J, L, N, Spc(I,J,L,N)
               PRINT*, 'Inf Spc ' // TRIM( LOCATION ) //
     &                 '. SET Spc TO BE LOWER LEVEL.'
               Spc(I,J,L,N) =  Spc(I,J,L-1,N)
!$OMP END CRITICAL

            ENDIF
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      Spc => NULL()

 100  FORMAT( ' Spc(',i3,',',i3,',',i3,',',i3,') = ', f13.6 )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@ REMOVE TRACERS MODIFICATION (bmy, 6/30/16)
!@@@ Need to restore State_Chm%TRACERS = State_Chm%SPECIES for testing 
!@@@
      DO NA = 1, nAdvect
         N                          = State_Chm%Map_Advect(NA)
         State_Chm%Tracers(:,:,:,N) = State_Chm%Species(:,:,:,N)
      ENDDO
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      END SUBROUTINE CHECK_STT_NESTED
!EOC
      END MODULE TRACER_MOD
