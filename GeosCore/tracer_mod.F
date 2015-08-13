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
      PUBLIC :: CHECK_STT_05x0666
      PUBLIC :: CHECK_STT_025x03125
      PUBLIC :: INIT_TRACER
!
! !DEFINED PARAMETERS:
!
      ! Conversion factor: molec air / kg air
      REAL(fp), PUBLIC, PARAMETER :: XNUMOLAIR = 
     &                        6.022e+23_fp / 28.9644e-3_fp
!
! !PUBLIC DATA MEMBERS:
!
#if   defined( APM )
      ! # of APM tracers
      INTEGER, PUBLIC           :: N_APMTRA   ! # of APM tracers
#endif
!
! !PRIVATE TYPES:

      ! Shadow copy: stores Input_Opt%N_TRACERS locally
      INTEGER, PRIVATE          :: N_TRACERS
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
      SUBROUTINE CHECK_STT( State_Chm, LOCATION )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : GEOS_CHEM_STOP
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE ERROR_MOD,          ONLY : IT_IS_FINITE
      USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN) :: LOCATION
!
! !INPUT PARAMETERS: 
!
      TYPE(ChmState),   INTENT(IN) :: State_Chm   ! Chemistry State object
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL                      :: LNEG, LNAN, LINF
      INTEGER                      :: I,    J,    L,   N, N_TOT_TRC
      
      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHECK_STT begins here!
      !=================================================================

      ! Pick the proper # of tracers
#if   defined( APM )
      N_TOT_TRC = N_TRACERS + N_APMTRA   ! G-C advected tracers + APM tracers
#else
      N_TOT_TRC = N_TRACERS              ! G-C advected tracers only
#endif

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Initialize
      LNEG = .FALSE.
      LNAN = .FALSE.
      LINF = .FALSE.

      ! Loop over grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TOT_TRC
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !---------------------------
         ! Check for Negatives
         !---------------------------
!clb fill nega
         STT(I,J,L,N) = ABS(STT(I,J,L,N))
         IF ( STT(I,J,L,N) < 0e+0_fp ) THEN 
!$OMP CRITICAL
            LNEG = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
            PRINT*, STT(I,J,L,N)
!$OMP END CRITICAL

         !---------------------------
         ! Check for NaN's
         !---------------------------
         ELSE IF ( IT_IS_NAN( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            LNAN = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
!$OMP END CRITICAL

         !----------------------------
         ! Check STT's for Infinities
         !----------------------------
         ELSE IF ( .not. IT_IS_FINITE( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            LINF = .TRUE.
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
!$OMP END CRITICAL            

         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

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

      END SUBROUTINE CHECK_STT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_stt_05x0666
!
! !DESCRIPTION: Subroutine CHECK\_STT\_05x0666 checks the STT tracer array for
!  negative values, NaN values, or Infinity values.  If any of these are found,
!  the STT array will be set to a specified value.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_STT_05x0666( State_Chm, LOCATION )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE ERROR_MOD,          ONLY : IT_IS_FINITE
      USE GIGC_State_Chm_Mod, ONLY : ChmState
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object
!
! !REVISION HISTORY:
!  05 Mar 2012 - M. Payer    - Initial version based on CHECK_STT and updates
!                              for nested grid by Yuxuan Wang.
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER                      :: I,    J,    L,   N
      
      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHECK_STT_05x0666 begins here!
      !=================================================================

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Loop over grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !---------------------------
         ! Check for Negatives
         !---------------------------
         IF ( STT(I,J,L,N) < 0e+0_fp ) THEN 
!$OMP CRITICAL
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
            PRINT*, 'Neg STT ' // TRIM( LOCATION ) //
     &              '. SET STT TO BE ZERO.'
            STT(I,J,L,N) = 0e+0_fp
!$OMP END CRITICAL

         !---------------------------
         ! Check for NaN's
         !---------------------------
         ELSE IF ( IT_IS_NAN( STT(I,J,L,N) ) ) THEN 
!$OMP CRITICAL
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
            PRINT*, 'NaN STT ' // TRIM( LOCATION ) //
     &              '. SET STT TO BE LOWER LEVEL.'
            STT(I,J,L,N) = STT(I,J,L-1,N)
!$OMP END CRITICAL

         !----------------------------
         ! Check STT's for Infinities
         !----------------------------
         ELSE IF ( .not. IT_IS_FINITE( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
            PRINT*, 'Inf STT ' // TRIM( LOCATION ) //
     &              '. SET STT TO BE LOWER LEVEL.'
            STT(I,J,L,N) =  STT(I,J,L-1,N)
!$OMP END CRITICAL

         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

 100  FORMAT( ' STT(',i3,',',i3,',',i3,',',i3,') = ', f13.6 )

      END SUBROUTINE CHECK_STT_05x0666
!EOC
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: check_stt_025x03125
!
! !DESCRIPTION: Subroutine CHECK\_STT\_025x03125 checks the STT tracer 
!  array for negative values, NaN values, or Infinity values.  If any of 
!  these are found, the STT array will be set to a specified value.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CHECK_STT_025x03125( State_Chm, LOCATION )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE ERROR_MOD,          ONLY : IT_IS_NAN
      USE ERROR_MOD,          ONLY : IT_IS_FINITE
!
! !INPUT PARAMETERS: 
!
      CHARACTER(LEN=*), INTENT(IN)    :: LOCATION
!
! !INPUT/OUTPUT PARAMETERS: 
!
      TYPE(ChmState),   INTENT(INOUT) :: State_Chm   ! Chemistry State object
! 
! !REVISION HISTORY:
!  05 Mar 2012 - M. Payer    - Initial version based on CHECK_STT and updates
!                              for nested grid by Yuxuan Wang.
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!  07 Jun 2013 - R. Yantosca - Now pass State_Chm object via the arg list
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER         :: I, J, L, N
      
      ! Pointers
      ! We need to define local arrays to hold corresponding values 
      ! from the Chemistry State (State_Chm) object. (mpayer, 12/6/12)
      REAL(fp), POINTER :: STT(:,:,:,:)

      !=================================================================
      ! CHECK_STT_05x0666 begins here!
      !=================================================================

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

      ! Loop over grid boxes
!$OMP PARALLEL DO 
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TRACERS
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         !---------------------------
         ! Check for Negatives
         !---------------------------
!clb         STT (I,J,L,N)=ABS(STT(I,J,L,N))
         IF ( STT(I,J,L,N) < 0e+0_fp ) THEN
!$OMP CRITICAL
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
            PRINT*, 'Neg STT ' // TRIM( LOCATION ) //
     &              '. SET STT TO BE ZERO.'
            STT(I,J,L,N) = 0e+0_fp
!$OMP END CRITICAL

         !---------------------------
         ! Check for NaN's
         !---------------------------
         ELSE IF ( IT_IS_NAN( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
            PRINT*, 'NaN STT ' // TRIM( LOCATION ) //
     &              '. SET STT TO BE LOWER LEVEL.'
            STT(I,J,L,N) = STT(I,J,L-1,N)
!$OMP END CRITICAL

         !----------------------------
         ! Check STT's for Infinities
         !----------------------------
         ELSE IF ( .not. IT_IS_FINITE( STT(I,J,L,N) ) ) THEN
!$OMP CRITICAL
            WRITE( 6, 100 ) I, J, L, N, STT(I,J,L,N)
            PRINT*, 'Inf STT ' // TRIM( LOCATION ) //
     &              '. SET STT TO BE LOWER LEVEL.'
            STT(I,J,L,N) =  STT(I,J,L-1,N)
!$OMP END CRITICAL

         ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

 100  FORMAT( ' STT(',i3,',',i3,',',i3,',',i3,') = ', f13.6 )

      ! Free pointer memory
      NULLIFY( STT )

      END SUBROUTINE CHECK_STT_025x03125
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_tracer
!
! !DESCRIPTION: Subroutine INIT\_TRACER initializes all module variables.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_TRACER( am_I_Root, Input_Opt, RC )
!
! !USES:
!
#if   defined( APM ) 
      USE APM_INIT_MOD,       ONLY : APM_NTRACERS
      USE APM_INIT_MOD,       ONLY : LAPM
#endif
      USE CMN_SIZE_MOD
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure
! 
! !REVISION HISTORY:
!  15 Jul 2004 - R. Yantosca - Initial version
!  25 Oct 2005 - R. Yantosca - Now allocate XNUMOL
!  16 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  05 Mar 2012 - M. Payer    - Added ProTeX headers
!  26 Feb 2013 - R. Yantosca - Now pass am_I_Root, Input_Opt, RC as arguments
!  07 Mar 2013 - R. Yantosca - Now assign RC=GIGC_SUCCESS outside of APM block
!  04 Apr 2013 - R. Yantosca - Removed STT (now in State_Chm)
!  25 Jun 2014 - R. Yantosca - Now is just used to set N_TRACERS & N_APMTRA
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: AS, N_TOT_TRC

#if   defined( APM ) 
      INTEGER :: N_APMTRAC1

      !=================================================================
      ! We are using APM microphysics!
      !=================================================================

      ! Force LAPM=T if we compile with APM=yes makefile option
      ! NOTE: LAPM is somewhat superfluous, since we now can use the
      ! C-preprocessor switch to activate the code for APM simulations
      ! (bmy, 2/17/11)
      LAPM      = .TRUE.

      ! Keep a shadow copy of Input_Opt%N_TRACERS
      N_TRACERS = Input_Opt%N_TRACERS

      ! Add additional tracers for APM
      CALL APM_NTRACERS( Input_Opt%N_TRACERS, N_APMTRAC1 )
      N_APMTRA  = N_APMTRAC1  
 
#else

      !=================================================================
      ! We are not using APM microphysics
      !=================================================================

      ! Keep a shadow copy of Input_Opt%N_TRACERS
      N_TRACERS = Input_Opt%N_TRACERS

#endif

      ! Assume success
      RC        = GIGC_SUCCESS

      END SUBROUTINE INIT_TRACER
!EOC
      END MODULE TRACER_MOD
