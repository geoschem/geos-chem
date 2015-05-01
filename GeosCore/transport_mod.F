!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: transport_mod
!
! !DESCRIPTION: Module TRANSPORT\_MOD is used to call the proper version of 
!  the TPCORE advection scheme for GCAP, GEOS-4, GEOS-5, or GEOS-5.7
!  nested-grid or global simulations.
!\\
!\\
! !INTERFACE: 
!
      MODULE TRANSPORT_MOD
! 
! !USES:
!
      USE PRECISION_MOD    ! For GEOS-Chem Precision (fp)

      IMPLICIT NONE
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_TRANSPORT
      PUBLIC  :: DO_TRANSPORT
      PUBLIC  :: INIT_TRANSPORT
      PUBLIC  :: INIT_GEOS5_WINDOW_TRANSPORT
      PUBLIC  :: INIT_GEOSFP_WINDOW_TRANSPORT
      PUBLIC  :: SET_TRANSPORT
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: GEOS4_GEOS5_GLOBAL_ADV      
      PRIVATE :: GCAP_GLOBAL_ADV             
      PRIVATE :: DO_GEOS5_WINDOW_TRANSPORT
      PRIVATE :: DO_GEOSFP_WINDOW_TRANSPORT      
      PRIVATE :: GET_AIR_MASS                
!
! !REVISION HISTORY:
!  10 Mar 2003 - Y. Wang, R. Yantosca - Initial version
!  (1 ) Now can select transport scheme for GEOS-3 winds.  Added code for PJC 
!        pressure fixer. (bdf, bmy, 5/8/03)
!  (2 ) Now delete DSIG array, it's obsolete.  Also added new PRIVATE function 
!        GET_AIR_MASS to compute air masses from the input/output pressures
!        from the new GEOS-4/fvDAS TPCORE. (bmy, 6/24/03)
!  (3 ) Now references DEBUG_MSG from "error_mod.f". (bmy, 8/7/03)
!  (4 ) Bug fix in DO_GLOBAL_TRANSPORT (bmy, 10/21/03)
!  (5 ) IORD, JORD, KORD are now module variables.  Now references 
!        "logical_mod.f" and "tracer_mod.f" (bmy, 7/20/04)
!  (6 ) Add mass-flux diagnostics to TPCORE_FVDAS (bdf, bmy, 9/28/04)
!  (7 ) Now references "diag_mod.f" (bmy, 9/28/04)
!  (8 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (9 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (10) Now flip arrays in call to TPCORE_FVDAS (bmy, 6/16/06)
!  (11) Added modifications for SUN compiler (bmy, 7/12/06)
!  (12) Bug fixes in DO_GLOBAL_TRANSPORT (bmy, 11/29/06)
!  (13) Split off GCAP, GEOS-3, GEOS-4/GEOS-5 specific calling sequences
!        into separate subroutines.  Also removed some obsolete module
!        variables. (bmy, 10/30/07)
!  (14) Modifications for GEOS-5 nested grid (yxw, dan, bmy, 11/6/08)
!  (15) Bug fix in mass balance in GCAP_GLOBAL_ADV and GEOS4_GEOS5_GLOBAL_ADV.
!        (ccc, 2/17/09)
!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTex Headers
!  08 Mar 2010 - C. Carouge  - Modify call to tpcore_fvdas. We do not re-order 
!                              mass fluxes diagnostics anymore.
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now reference new grid_mod.F90
!  21 Jun 2012 - R. Yantosca - Now use pointers to flip arrays in vertical
!  21 Jun 2012 - R. Yantosca - Comment out GEOS-3 window subroutine
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed "GEOS57" to "GEOSFP" in routine names
!  24 Nov 2014 - M. Yannetti - Added PRECISION_MOD
!  31 Mar 2015 - E. Lundgren - Modifications to change moist air mass to dry
!  15 Apr 2015 - R. Yantosca - Add TARGET spec to A_M2 array
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! (1 ) Ap     (REAL(fp) ) : Vertical coordinate array for TPCORE
      ! (2 ) A_M2   (REAL(fp) ) : Grid box surface areas [m2]
      ! (3 ) Bp     (REAL(fp) ) : Vertical coordinate array for TPCORE
      ! (4 ) IORD   (REAL(fp) ) : TPCORE E/W option flag
      ! (5 ) JORD   (REAL(fp) ) : TPCORE N/S option flag
      ! (6 ) KORD   (REAL(fp) ) : TPCORE vertical option flag
      ! (7 ) JLAST  (INTEGER) : For fvDAS TPCORE
      ! (8 ) MG     (INTEGER) : For fvDAS TPCORE
      ! (9 ) NG     (INTEGER) : For fvDAS TPCORE
      ! (10) N_ADJ  (INTEGER) : For fvDAS TPCORE
      !=================================================================
      INTEGER                       :: IORD,  JORD, KORD, JFIRST 
      INTEGER                       :: JLAST, NG,   MG,   N_ADJ
      REAL(fp), ALLOCATABLE         :: Ap(:)
      REAL(fp), ALLOCATABLE         :: Bp(:)
      REAL(fp), ALLOCATABLE, TARGET :: A_M2(:)

      ! This seems to be obsolete ???
      !REAL(fp), ALLOCATABLE :: STT_BC2(:,:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_transport
!
! !DESCRIPTION: Subroutine DO\_TRANSPORT is the driver routine for the proper 
!  TPCORE program for GEOS-3, GEOS-4/GEOS-5, or window simulations.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_TRANSPORT( am_I_Root, Input_Opt,
     &                         State_Met, State_Chm, RC )
!
! !USES:
!
      USE CMN_SIZE_MOD
      USE GRID_MOD,           ONLY : ITS_A_NESTED_GRID
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE TPCORE_BC_MOD,      ONLY : INIT_TPCORE_BC
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  10 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Removed IORD, JORD, KORD from the arg list.  Also now removed
!        reference to CMN, it's not needed. (bmy, 7/20/04)
!  (2 ) Now call separate routines for different met fields. (bmy, 10/30/07)
!  (3 ) Now references subroutine INIT_TPCORE_BC from tpcore_bc_mod.f and
!        DO_GEOS5_FVDAS_WINDOW_TRANSPORT from 
!        "tpcore_geos5_fvdas_window_mod.f90". (yxw, dan, bmy, 11/6/08)
!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  06 Oct 2010 - R. Yantosca - Treat MERRA in the same way as GEOS-5.
!  28 Feb 2012 - R. Yantosca - Treat GEOS-5.7 in the same way as MERRA
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  31 Mar 2015 - E. Lundgren - Make State_Met INTENT(INOUT) for AIRQNT calls
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      LOGICAL, SAVE     :: FIRST = .TRUE.

      !=================================================================
      ! DO_TRANSPORT begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN 

#if   defined( GRID05x0666 )

      ! GEOS-5.2.0 nested simulations
      CALL INIT_GEOS5_WINDOW_TRANSPORT
     &   ( am_I_Root, Input_Opt, RC )

#elif defined( GRID025x03125 )

      ! GEOS-5.7.2 (aka "FP") nested grids
      CALL INIT_GEOSFP_WINDOW_TRANSPORT
     &     ( am_I_Root, Input_Opt, RC )    
#else

      ! All global simulations
      CALL INIT_TRANSPORT                   
     &     ( am_I_Root, Input_Opt, RC )    

#endif
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Choose the proper version of TPCORE for the nested-grid window 
      ! region (usually 1x1 grids) or for the entire globe
      !=================================================================
      IF ( ITS_A_NESTED_GRID() ) THEN

#if defined(  GRID025x03125 )

         ! Nested-grid simulation with GEOS 5.7-x met
         CALL DO_GEOSFP_WINDOW_TRANSPORT
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

#elif defined( GRID05x0666 )

         ! Nested-grid simulation with GEOS-5.2.0 met
         CALL DO_GEOS5_WINDOW_TRANSPORT
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

#endif

      !=================================================================
      ! Choose the proper version of TPCORE for global simulations
      !=================================================================
      ELSE

#if defined( GCAP )

         ! Call TPCORE w/ proper settings for GCAP met
         CALL GCAP_GLOBAL_ADV
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

#else

         ! Call TPCORE w/ proper settings for GEOS-4/GEOS-5/MERRA/GEOS-5.7
         CALL GEOS4_GEOS5_GLOBAL_ADV
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )

#endif
      ENDIF

      ! Return to calling program
      END SUBROUTINE DO_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: geos4_geos5_global_adv
!
! !DESCRIPTION: Subroutine GEOS4\_GEOS5\_GLOBAL\_ADV is the driver routine 
!  for TPCORE with the GMAO GEOS-4 or GEOS-5 met fields.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GEOS4_GEOS5_GLOBAL_ADV( am_I_Root, Input_Opt,
     &                                   State_Met, State_Chm, RC )
!
! !USES:
!
      USE DIAG_MOD,           ONLY : MASSFLEW,  MASSFLNS,  MASSFLUP
      USE ERROR_MOD,          ONLY : IT_IS_NAN, DEBUG_MSG, SAFE_DIV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PJC_PFIX_MOD,       ONLY : DO_PJC_PFIX
      USE PRESSURE_MOD,       ONLY : SET_FLOATING_PRESSURE
      USE DAO_MOD,            ONLY : AIRQNT
      USE TIME_MOD,           ONLY : GET_TS_DYN
      USE TPCORE_BC_MOD,      ONLY : SAVE_GLOBAL_TPCORE_BC
      USE TPCORE_FVDAS_MOD,   ONLY : TPCORE_FVDAS
      USE UNITCONV_MOD
#if   defined( APM )
      USE TRACER_MOD,         ONLY : N_APMTRA
#endif

      USE CMN_SIZE_MOD             ! Size parameters        
      USE CMN_DIAG_MOD             ! NDxx flags
      USE CMN_GCTM_MOD             ! Physical constants
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  30 Oct 2007 - R. Yantosca - Initial version
!  (1 ) Split off the GEOS-4 & GEOS-5 relevant parts from the previous 
!        routine DO_GLOBAL_TRANSPORT (bmy, 10/30/07)
!  (2 ) Activate the call to SAVE_GLOBAL_TPCORE_BC (yxw, dan, bmy, 11/6/08)
!  (3 ) Bug fix in mass balance: only account for cells of STT with non-zero
!        concentrations when doing the computation (ccc, bmy, 2/17/09)
!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  16 Feb 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!  21 Jun 2012 - R. Yantosca - Now use pointers to flip indices in vertical
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  04 Jun 2013 - R. Yantosca - Replace pointer STT with State_Chm%Tracers
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE.
!  26 Mar 2015 - E. Lundgren - Pass State_Met to get_air_mass for dry mass calc
!  31 Mar 2015 - E. Lundgren - Move AIRQNT call from main to within transport
!  16 Apr 2015 - E. Lundgren - Convert mixing ratio to moist before advection
!                              and back to dry afterwards
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, L, L2, N, N_DYN, N_TOT_TRC
      REAL(fp)           :: A_DIFF, D_DYN, TR_DIFF
      REAL(fp)           :: AD_A(IIPAR,JJPAR,LLPAR)
      REAL(fp)           :: AD_B(IIPAR,JJPAR,LLPAR)
      REAL(fp)           :: P_TP1(IIPAR,JJPAR)
      REAL(fp)           :: P_TP2(IIPAR,JJPAR)
      REAL(fp)           :: P_TEMP(IIPAR,JJPAR)
      REAL(fp)           :: TR_A(IIPAR,JJPAR,LLPAR)
#if   defined( APM )
      REAL(fp)           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                         Input_Opt%N_TRACERS+N_APMTRA)  ! With APM
#else
      REAL(fp)           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                         Input_Opt%N_TRACERS)           ! W/o APM
#endif
      REAL(fp),  TARGET  :: XMASS(IIPAR,JJPAR,LLPAR) 
      REAL(fp),  TARGET  :: YMASS(IIPAR,JJPAR,LLPAR) 

      ! Pointers
      REAL(fp),  POINTER :: p_UWND   (:,:,:  )
      REAL(fp),  POINTER :: p_VWND   (:,:,:  )
      REAL(fp),  POINTER :: p_XMASS  (:,:,:  )
      REAL(fp),  POINTER :: p_YMASS  (:,:,:  )
      REAL(fp),  POINTER :: p_MOISTMW(:,:,:  )
      REAL(fp),  POINTER :: p_STT    (:,:,:,:)

      ! Variable to ensure mass conservation (ccc, 2/17/09)
      REAL(fp)           :: SUMADA

      ! For fields from Input_Opt
      LOGICAL            :: LFILL
      LOGICAL            :: LPRT
      LOGICAL            :: LWINDO
      INTEGER            :: N_TRACERS
      REAL(fp)           :: TCVV(Input_Opt%N_TRACERS)

      !=================================================================
      ! GEOS4_GEOS5_GLOBAL_ADV begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LFILL     = Input_Opt%LFILL
      LPRT      = Input_Opt%LPRT
      LWINDO    = Input_Opt%LWINDO
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV(1:N_TRACERS)

      ! Pick the proper # of tracers
#if   defined( APM )
      N_TOT_TRC = N_TRACERS + N_APMTRA   ! G-C advected tracers + APM tracers
#else
      N_TOT_TRC = N_TRACERS              ! G-C advected tracers only
#endif

      ! Dynamic timestep [s]
      N_DYN = GET_TS_DYN() * 60
      D_DYN = DBLE( N_DYN )

      !=================================================================
      ! Prepare variables for calls to PJC pressure-fixer and TPCORE
      !
      ! For GEOS-4 and GEOS-5 (hybrid grids), the pressure at the 
      ! bottom edge of grid box (I,J,L) is given by:
      !
      !    P(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
      !
      ! where Psurface is the true surface pressure (i.e. not PS-PTOP).
      ! and Ap(L), Bp(L) define the vertical grid (see pressure_mod.f)
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! True surface pressure at midpoint of dynamic timestep [hPa]
         P_TP1(I,J) = State_Met%PEDGE(I,J,1)

         ! True surface pressure at end of dynamic timestep [hPa]
         P_TP2(I,J) = State_Met%PSC2(I,J)    
   
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      !=================================================================
      ! Get the air and tracer mass before advection
      !=================================================================

      ! Dry airmass [kg] before transport
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         AD_B(I,J,L) = GET_AIR_MASS( I, J, L, P_TP1(I,J), State_Met )
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Tracer mass [kg] before transport
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TOT_TRC
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         TR_B(I,J,L,N) = State_Chm%Tracers(I,J,L,N)
     &                 * AD_B(I,J,L)
     &                 / TCVV(N)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Call the PJC/LLNL pressure fixer to get the adjusted air 
      ! masses XMASS and YMASS.  XMASS and YMASS need to be passed to 
      ! TPCORE_FVDAS in order to ensure mass conservation.
      !=================================================================

      ! NOTE: P_TP1 and P_TP2 are the true surface pressures!
      CALL DO_PJC_PFIX( D_DYN,        
     &                  P_TP1,        P_TP2,
     &                  State_Met%U,  State_Met%V,
     &                  XMASS,        YMASS )

      !=================================================================
      ! Call TPCORE_FVDAS to perform the advection
      !=================================================================

      ! Flip array indices in the vertical using pointer storage
      p_UWND    => State_Met%U      (:,:,LLPAR:1:-1  )
      p_VWND    => State_Met%V      (:,:,LLPAR:1:-1  )
      p_STT     => State_Chm%Tracers(:,:,LLPAR:1:-1,:) ! v/v
      p_XMASS   => XMASS            (:,:,LLPAR:1:-1  )
      p_YMASS   => YMASS            (:,:,LLPAR:1:-1  )
      p_MOISTMW => State_Met%MOISTMW(:,:,LLPAR:1:-1  )

      ! Convert tracer concentration from mol/mol dry air to 
      ! mol/mol moist air before advection (ewl, 4/16/15)
      CALL Convert_DryVV_to_MoistVV( am_I_Root, N_TRACERS, State_Met, 
     &                               State_Chm, RC )

      ! Do the advection
      ! Note: the mass flux diagnostic arrays (MASSFLEW, MASSFLNS and MASSFLUP)
      ! are incremented upside-down (level 1 = top of the atmosphere).
      ! The levels order is reversed only when written out in diag3.f
      ! (ccc, 3/8/10)
      CALL TPCORE_FVDAS( D_DYN,  Re,        IIPAR,    JJPAR,
     &                   LLPAR,  JFIRST,    JLAST,    NG,
     &                   MG,     N_TOT_TRC, Ap,       Bp,
     &                   p_UWND, p_VWND,    P_TP1,    P_TP2,
     &                   P_TEMP, p_STT,     IORD,     JORD,      
     &                   KORD,   N_ADJ,     p_XMASS,  p_YMASS,         
     &                   LFILL,  MASSFLEW,  MASSFLNS, MASSFLUP,    
     &                   A_M2,   TCVV,      ND24,     ND25,     
     &                   ND26,   p_MOISTMW )

      ! Convert tracer concentration from mol/mol moist air
      ! mol/mol dry air after advection (ewl, 4/16/15)
      CALL Convert_MoistVV_to_DryVV( am_I_Root, N_TRACERS, State_Met, 
     &                               State_Chm, RC ) 

      ! Free pointer memory
      NULLIFY( p_UWND, p_VWND, p_STT, p_XMASS, p_YMASS, p_MOISTMW )

      !=================================================================
      ! Reset surface pressure and ensure mass conservation
      !=================================================================

      ! Reset the floating surface pressure with P_TP2, the "true"
      ! surface pressure at the end of the dynamic timestep.
      CALL SET_FLOATING_PRESSURE( P_TP2 )

      ! Update State_Met air quantities with new pressure prior to updating 
      ! tracer concentrations (ewl, 3/31/15)
      CALL AIRQNT( am_I_Root, State_Met, RC )

      ! Adjust tracer to correct residual non-conservation of mass
      ! This was changed to be applied only for cells with tracer 
      ! concentration. (ccc, 11/20/08) 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, SUMADA, AD_A, TR_A, TR_DIFF )
      DO N = 1, N_TOT_TRC

         ! Zero summing variable
         SUMADA = 0.e+0_fp

         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Dry air mass [kg] after transport
            AD_A(I,J,L) = GET_AIR_MASS( I, J, L, P_TP2(I,J), State_Met )
         
            ! Tracer mass [kg] after transport
            TR_A(I,J,L) = State_Chm%Tracers(I,J,L,N)
     &                  * AD_A(I,J,L)
     &                  / TCVV(N)

            IF ( State_Chm%Tracers(I,J,L,N) > 0.e+0_fp   .or.
     &           State_Chm%Tracers(I,J,L,N) < 0.e+0_fp ) THEN
               SUMADA = SUMADA + AD_A(I,J,L)
            ENDIF
         ENDDO
         ENDDO
         ENDDO

         ! Residual mass difference [kg]: before - after
         TR_DIFF = SUM( TR_B(:,:,:,N) ) - SUM( TR_A )

         ! Convert from [kg] to [v/v]
         TR_DIFF = SAFE_DIV( TR_DIFF, SUMADA, 0.e+0_fp ) * TCVV(N)

         ! Add mass difference [v/v] back to STT
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            IF ( State_Chm%Tracers(I,J,L,N) > 0.e+0_fp   .or.
     &           State_Chm%Tracers(I,J,L,N) < 0.e+0_fp ) THEN
               State_Chm%Tracers(I,J,L,N) =
     &         State_Chm%Tracers(I,J,L,N) + TR_DIFF
            ENDIF

            State_Chm%Tracers(I,J,L,N) =
     &           MAX( State_Chm%Tracers(I,J,L,N), 0e+0_fp )

         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### G4_G5_GLOB_ADV: a TPCORE' ) 

      END SUBROUTINE GEOS4_GEOS5_GLOBAL_ADV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: gcap_global_adv
!
! !DESCRIPTION: Subroutine GCAP\_GLOBAL\_ADV is the driver routine for TPCORE 
!  with the GCAP/GISS met fields.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE GCAP_GLOBAL_ADV( am_I_Root, Input_Opt,
     &                            State_Met, State_Chm, RC )
!
! !USES:
!
      USE DIAG_MOD,           ONLY : MASSFLEW,  MASSFLNS,  MASSFLUP
      USE ERROR_MOD,          ONLY : IT_IS_NAN, DEBUG_MSG, SAFE_DIV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PJC_PFIX_MOD,       ONLY : DO_PJC_PFIX
      USE PRESSURE_MOD,       ONLY : SET_FLOATING_PRESSURE
      USE DAO_MOD,            ONLY : AIRQNT
      USE TIME_MOD,           ONLY : GET_TS_DYN
      USE TPCORE_FVDAS_MOD,   ONLY : TPCORE_FVDAS
      USE UNITCONV_MOD
#if   defined( APM )
      USE TRACER_MOD,         ONLY : N_APMTRA
#endif

      USE CMN_SIZE_MOD             ! Size parameters
      USE CMN_DIAG_MOD             ! NDxx flags
      USE CMN_GCTM_MOD             ! Physical constants
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  30 Oct 2007 - R. Yantosca - Initial version
!  (1 ) Split off the GCAP relevant parts from the previous routine
!        DO_GLOBAL_TRANSPORT (bmy, 10/30/07)
!  (2 ) Bug fix in mass balance: only account for cells of STT with non-zero
!        concentrations when doing the computation (ccc, bmy, 2/17/09)
!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  16 Feb 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!  21 Jun 2012 - R. Yantosca - Now use pointers to flip arrays in the vertical
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  04 Jun 2013 - R. Yantosca - Replace pointer STT with State_Chm%Tracers
!  04 Sep 2014 - R. Yantosca - Avoid div-by-zero errors in division w/ SUMADA
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE.
!  26 Mar 2015 - E. Lundgren - Pass State_Met to get_air_mass for dry mass calc
!  31 Mar 2015 - E. Lundgren - Move AIRQNT call from main to within transport
!  16 Apr 2015 - E. Lundgren - Convert mixing ratio to moist before advection
!                              and back to dry afterwards
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I, J, L, L2, N, N_DYN, N_TOT_TRC
      REAL(fp)           :: A_DIFF, D_DYN, TR_DIFF
      REAL(fp)           :: AD_A(IIPAR,JJPAR,LLPAR)
      REAL(fp)           :: AD_B(IIPAR,JJPAR,LLPAR)
      REAL(fp)           :: P_TP1(IIPAR,JJPAR)
      REAL(fp)           :: P_TP2(IIPAR,JJPAR)
      REAL(fp)           :: P_TEMP(IIPAR,JJPAR)
      REAL(fp)           :: TR_A(IIPAR,JJPAR,LLPAR)
#if   defined( APM )
      REAL(fp)           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                         Input_Opt%N_TRACERS+N_APMTRA)
#else                  
      REAL(fp)           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                         Input_Opt%N_TRACERS)
#endif
      REAL(fp)           :: VTMP(IIPAR,JJPAR,LLPAR)
      REAL(fp),  TARGET  :: XMASS(IIPAR,JJPAR,LLPAR) 
      REAL(fp),  TARGET  :: YMASS(IIPAR,JJPAR,LLPAR) 

      ! Pointers
      REAL(fp),  POINTER :: p_UWND   (:,:,:  )
      REAL(fp),  POINTER :: p_VWND   (:,:,:  )
      REAL(fp),  POINTER :: p_XMASS  (:,:,:  )
      REAL(fp),  POINTER :: p_YMASS  (:,:,:  )
      REAL(fp),  POINTER :: p_MOISTMW(:,:,:  )
      REAL(fp),  POINTER :: p_STT    (:,:,:,:) ! v/v

      ! For fields from Input_Opt
      LOGICAL            :: LFILL
      LOGICAL            :: LPRT
      LOGICAL            :: LWINDO
      INTEGER            :: N_TRACERS
      REAL(fp)           :: TCVV(Input_Opt%N_TRACERS)

      ! Variable to ensure mass conservation (ccc, 20/11/08)
      REAL(fp)           :: SUMADA

      !=================================================================
      ! GCAP_GLOBAL_ADV begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LFILL     = Input_Opt%LFILL
      LPRT      = Input_Opt%LPRT
      LWINDO    = Input_Opt%LWINDO
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV(1:N_TRACERS)

      ! Pick the proper # of tracers
#if   defined( APM )
      N_TOT_TRC = N_TRACERS + N_APMTRA   ! G-C advected tracers + APM tracers
#else
      N_TOT_TRC = N_TRACERS              ! G-C advected tracers only
#endif

      ! Dynamic timestep [s]
      N_DYN = GET_TS_DYN() * 60
      D_DYN = DBLE( N_DYN )

      !=================================================================
      ! Prepare variables for calls to PJC presure-fixer and TPCORE
      !
      ! For GCAP (hybrid grid, but expressed as a pure-sigma grid), the 
      ! pressure at the bottom edge grid box (I,J,L) is given by:
      !
      !    P(I,J,L) = Ap(L) + [ Bp(L) * ( Psurface(I,J) - PTOP ) ]
      !
      ! where Psurface is the true surface pressure (i.e. not PS-PTOP).
      ! and Ap(L), Bp(L) define the vertical grid (see pressure_mod.f)
      !
      ! Therefore, we construct the 3-D pressure edge arrays PLE_TP1
      ! and PLE_TP2 according to the above equation.  Note that PLE_TP1 
      ! and PLE_TP2 are inverted (i.e. L=1 is atm top) for compatibility 
      ! with TPCORE_FVDAS.
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Psurface - PTOP at midpoint of dynamic timestep [hPa]
         P_TP1(I,J) = State_Met%PEDGE(I,J,1) - PTOP

         ! Psurface - PTOP at end of dynamic timestep [hPa]
         P_TP2(I,J) = State_Met%PSC2(I,J) - PTOP

      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      !==============================================================
      ! Get the air & tracer mass before advection
      !==============================================================

      ! Dry airmass [kg] before transport
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         AD_B(I,J,L) = GET_AIR_MASS( I, J, L, P_TP1(I,J), State_Met )
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Tracer mass [kg] before transport
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TOT_TRC
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         TR_B(I,J,L,N) = State_Chm%Tracers(I,J,L,N)
     &                 * AD_B(I,J,L)
     &                 / TCVV(N)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Call the PJC/LLNL pressure fixer to get the adjusted air 
      ! masses XMASS and YMASS.  XMASS and YMASS need to be passed to 
      ! TPCORE_FVDAS in order to ensure mass conservation.
      !=================================================================

      ! NOTE: P_TP1+PTOP & P_TP2+PTOP are the true surface pressures!
      CALL DO_PJC_PFIX( D_DYN,
     &                  P_TP1+PTOP,   P_TP2+PTOP, 
     &                  State_Met%U,  State_Met%V,
     &                  XMASS,        YMASS )

      !=================================================================
      ! Call TPCORE_FVDAS to perform the advection
      !=================================================================

      ! Flip array indices in the vertical using pointer storage
      p_UWND    => State_Met%U      (:,:,LLPAR:1:-1  )
      p_VWND    => State_Met%V      (:,:,LLPAR:1:-1  )
      p_STT     => State_Chm%Tracers(:,:,LLPAR:1:-1,:)
      p_XMASS   => XMASS            (:,:,LLPAR:1:-1  )
      p_YMASS   => YMASS            (:,:,LLPAR:1:-1  )
      p_MOISTMW => State_Met%MOISTMW(:,:,LLPAR:1:-1  )

      ! Convert tracer concentration from mol/mol dry air to 
      ! mol/mol moist air before advection (ewl, 4/16/15)
      CALL Convert_DryVV_to_MoistVV( am_I_Root, N_TRACERS, State_Met, 
     &                               State_Chm, RC )

      ! Note: the mass flux diagnostic arrays (MASSFLEW, MASSFLNS and MASSFLUP)
      ! are incremented upside-down (level 1 = top of the atmosphere).
      ! The levels order is reversed only when written out in diag3.f
      ! (ccc, 3/8/10)
      ! Do the advection
      CALL TPCORE_FVDAS( D_DYN,   Re,        IIPAR,    JJPAR,
     &                   LLPAR,   JFIRST,    JLAST,    NG,
     &                   MG,      N_TOT_TRC, Ap,       Bp,
     &                   p_UWND,  p_VWND,    P_TP1,    P_TP2,
     &                   P_TEMP,  p_STT,     IORD,     JORD,      
     &                   KORD,    N_ADJ,     p_XMASS,  p_YMASS,
     &                   LFILL,   MASSFLEW,  MASSFLNS, MASSFLUP, 
     &                   A_M2,    TCVV,      ND24,     ND25,     
     &                   ND26,    p_MOISTMW )

      ! Convert tracer concentration from mol/mol moist air to 
      ! mol/mol dry air after advection (ewl, 4/16/15)      
      CALL Convert_MoistVV_to_DryVV( am_I_Root, N_TRACERS, State_Met, 
     &                               State_Chm, RC ) 

      ! Free pointer memory
      NULLIFY( p_UWND, p_VWND, p_STT, p_XMASS, p_YMASS, p_MOISTMW )

      !=================================================================
      ! Reset surface pressure and ensure mass conservation
      !=================================================================

      ! Reset the floating surface pressure with P_TP2+PTOP, the "true"
      ! surface pressure at the end of the dynamic timestep.
      CALL SET_FLOATING_PRESSURE( P_TP2 + PTOP )

      ! Update State_Met air quantities with new pressure prior to updating 
      ! tracer concentrations (ewl, 3/31/15)
      CALL AIRQNT( am_I_Root, State_Met, RC )

      ! Adjust tracer to correct residual non-conservation of mass
      ! This was changed to be applied only for cells with tracer 
      ! concentration. (ccc, 11/20/08) 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, SUMADA, AD_A, TR_A, TR_DIFF )
      DO N = 1, N_TOT_TRC

         ! Zero summing variable
         SUMADA = 0.e+0_fp

         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Dry air mass [kg] after transport
            AD_A(I,J,L) = GET_AIR_MASS( I, J, L, P_TP2(I,J), State_Met )
         
            ! Tracer mass [kg] after transport
            TR_A(I,J,L) = State_Chm%Tracers(I,J,L,N)
     &                  * AD_A(I,J,L)
     &                  / TCVV(N)

            ! We apply mass conservation only on cells with 
            ! nonzero tracer. (ccc, 11/20/08)
            IF ( State_Chm%Tracers(I,J,L,N) > 0.e+0_fp   .or.
     &           State_Chm%Tracers(I,J,L,N) < 0.e+0_fp ) THEN
               SUMADA = SUMADA + AD_A(I,J,L)
            ENDIF
         ENDDO
         ENDDO
         ENDDO

         ! Residual mass difference [kg]: before - after
         TR_DIFF = SUM( TR_B(:,:,:,N) ) - SUM( TR_A )

         ! Convert from [kg] to [v/v]
         TR_DIFF = SAFE_DIV( TR_DIFF, SUMADA, 0.e+0_fp ) * TCVV(N)

         ! Add mass difference [v/v] back to STT
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            IF ( State_Chm%Tracers(I,J,L,N) > 0.e+0_fp   .or. 
     &           State_Chm%Tracers(I,J,L,N) < 0.e+0_fp ) THEN
               State_Chm%Tracers(I,J,L,N) =
     &         State_Chm%Tracers(I,J,L,N) + TR_DIFF
            ENDIF

            State_Chm%Tracers(I,J,L,N) =
     &           MAX( State_Chm%Tracers(I,J,L,N), 0e+0_fp )
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### GCAP_GLOB_ADV: a TPCORE' ) 
      
      END SUBROUTINE GCAP_GLOBAL_ADV
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_geos5_window_transport
!
! !DESCRIPTION: Subroutine DO\_GEOS5\_WINDOW\_TRANSPORT is the driver program 
!  for the proper TPCORE program for the GEOS-5 nested-grid simulations. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_GEOS5_WINDOW_TRANSPORT( am_I_Root, Input_Opt,
     &                                      State_Met, State_Chm, RC )
!
! !USES:
!
      ! References to F90 modules
      USE DIAG_MOD,                  ONLY : MASSFLEW, MASSFLNS, MASSFLUP
      USE ERROR_MOD,                 ONLY : IT_IS_NAN,     DEBUG_MSG
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod,        ONLY : OptInput
      USE GIGC_State_Chm_Mod,        ONLY : ChmState
      USE GIGC_State_Met_Mod,        ONLY : MetState
      USE GRID_MOD,                  ONLY : GET_XOFFSET,   GET_YOFFSET 
      USE PJC_PFIX_GEOS5_WINDOW_MOD, ONLY : DO_PJC_PFIX_GEOS5_WINDOW
      USE PRESSURE_MOD,              ONLY : SET_FLOATING_PRESSURE
      USE DAO_MOD,                   ONLY : AIRQNT
      USE TIME_MOD,                  ONLY : GET_TS_DYN
      USE TPCORE_BC_MOD,             ONLY : I0_W, J0_W, I1_W, J1_W
      USE TPCORE_BC_MOD,             ONLY : I2_W, J2_W, IM_W, JM_W, IGZD
      USE TPCORE_BC_MOD,             ONLY : DO_WINDOW_TPCORE_BC
      USE TPCORE_GEOS5_WINDOW_MOD,   ONLY : TPCORE_GEOS5_WINDOW
      USE UNITCONV_MOD
#if   defined( APM )
      USE TRACER_MOD,                ONLY : N_APMTRA
#endif

      USE CMN_SIZE_MOD                    ! Size parameters
      USE CMN_DIAG_MOD                    ! NDxx flags
      USE CMN_GCTM_MOD                    ! Physical constants
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  10 Mar 2003 - R. Yantosca - Initial version

!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  16 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  21 Jun 2012 - R. Yantosca - Now use pointers to flip arrays in the vertical
!  09 Nov 2012 - M. Payer    - Replaced all met field arrays with State_Met
!                              derived type object
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  04 Jun 2013 - R. Yantosca - Replace pointer STT with State_Chm%Tracers
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE.
!  26 Mar 2015 - E. Lundgren - Pass State_Met to get_air_mass for dry mass calc
!  31 Mar 2015 - E. Lundgren - Move AIRQNT call from main to within transport
!  16 Apr 2015 - E. Lundgren - Convert mixing ratio to moist before advection
!                              and back to dry afterwards
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I0,J0
      INTEGER            :: I, J, L, L2, N, N_DYN, N_TOT_TRC
      REAL(fp)           :: A_DIFF, D_DYN, TR_DIFF
      REAL(fp)           :: AD_A(IIPAR,JJPAR,LLPAR)
      REAL(fp)           :: AD_B(IIPAR,JJPAR,LLPAR)
      REAL(fp)           :: P_TP1(IIPAR,JJPAR)
      REAL(fp)           :: P_TP2(IIPAR,JJPAR)
      REAL(fp)           :: P_TEMP(IIPAR,JJPAR)
      REAL(fp)           :: TR_A(IIPAR,JJPAR,LLPAR)
#if   defined( APM )
      REAL(fp)           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                         Input_Opt%N_TRACERS+N_APMTRA) 
#else
      REAL(fp)           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                         Input_Opt%N_TRACERS)
#endif
      REAL(fp),  TARGET  :: XMASS(IIPAR,JJPAR,LLPAR)
      REAL(fp),  TARGET  :: YMASS(IIPAR,JJPAR,LLPAR)

      ! Pointers
      REAL(fp),  POINTER :: p_UWND   (:,:,:  )
      REAL(fp),  POINTER :: p_VWND   (:,:,:  )
      REAL(fp),  POINTER :: p_XMASS  (:,:,:  )
      REAL(fp),  POINTER :: p_YMASS  (:,:,:  )
      REAL(fp),  POINTER :: p_MOISTMW(:,:,:  )
      REAL(fp),  POINTER :: p_STT    (:,:,:,:)

      ! For fields from Input_Opt
      LOGICAL            :: LFILL
      LOGICAL            :: LPRT
      LOGICAL            :: LWINDO
      INTEGER            :: N_TRACERS
      REAL(fp)           :: TCVV(Input_Opt%N_TRACERS)

      ! Variable to ensure mass conservation (ccc, 20/11/08)
      REAL(fp)           :: SUMADA

      !=================================================================
      ! DO_GEOS5_WINDOW_TRANSPORT begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LFILL     = Input_Opt%LFILL
      LPRT      = Input_Opt%LPRT
      LWINDO    = Input_Opt%LWINDO
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV(1:N_TRACERS)

      ! Pick the proper # of tracers
#if   defined( APM )
      N_TOT_TRC = N_TRACERS + N_APMTRA   ! G-C advected tracers + APM tracers
#else
      N_TOT_TRC = N_TRACERS              ! G-C advected tracers only
#endif

      ! Get nested-grid lon/lat offsets [# boxes]
      I0    = GET_XOFFSET( GLOBAL=.TRUE. )
      J0    = GET_YOFFSET( GLOBAL=.TRUE. )

      ! Dynamic timestep [s]
      N_DYN = GET_TS_DYN() * 60
      D_DYN = DBLE( N_DYN )

      !=================================================================
      ! Prepare variables for calls to PJC pressure-fixer and TPCORE
      !
      ! For GEOS-4 and GEOS-5 (hybrid grids), the pressure at the
      ! bottom edge of grid box (I,J,L) is given by:
      !
      !    P(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
      !
      ! where Psurface is the true surface pressure (i.e. not PS-PTOP).
      ! and Ap(L), Bp(L) define the vertical grid (see pressure_mod.f)
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! True surface pressure at midpoint of dynamic timestep [hPa]
         P_TP1(I,J) = State_Met%PEDGE(I,J,1)

         ! True surface pressure at end of dynamic timestep [hPa]
         P_TP2(I,J) = State_Met%PSC2(I,J)
  
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Get the air and tracer mass before advection
      !=================================================================

      ! Dry airmass [kg] before transport
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         AD_B(I,J,L) = GET_AIR_MASS( I, J, L, P_TP1(I,J), State_Met )
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Tracer mass [kg] before transport
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TOT_TRC
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         TR_B(I,J,L,N) = State_Chm%Tracers(I,J,L,N)
     &                 * AD_B(I,J,L)
     &                 / TCVV(N)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Call the PJC/LLNL pressure fixer to get the adjusted air
      ! masses XMASS and YMASS.  XMASS and YMASS need to be passed to
      ! TPCORE_FVDAS in order to ensure mass conservation.
      !=================================================================
      XMASS = 0e+0_fp !(dan)
      YMASS = 0e+0_fp

      ! NOTE: P_TP1 and P_TP2 are the true surface pressures!
      CALL DO_PJC_PFIX_GEOS5_WINDOW( D_DYN,
     &                               P_TP1,        P_TP2,
     &                               State_Met%U,  State_Met%V,
     &                               XMASS,        YMASS )

      IF ( LPRT ) CALL DEBUG_MSG( '### GEOS5_NESTED: a GET_TRA_MASS' )

      ! Impose TPCORE boundary conditions @ edges of nested grid
      CALL DO_WINDOW_TPCORE_BC
     &   ( am_I_Root, Input_Opt, State_Chm, RC )

      IF ( LPRT ) CALL DEBUG_MSG( '### GEOS5_NESTED: a DO_WIN_TP_BC' )

      ! Flip array indices in the vertical using pointer storage
      p_UWND  => State_Met%U        (:,:,LLPAR:1:-1  )
      p_VWND  => State_Met%V        (:,:,LLPAR:1:-1  )
      p_STT   => State_Chm%Tracers  (:,:,LLPAR:1:-1,:)
      p_XMASS => XMASS              (:,:,LLPAR:1:-1  )
      p_YMASS => YMASS              (:,:,LLPAR:1:-1  )
      p_MOISTMW => State_Met%MOISTMW(:,:,LLPAR:1:-1  )

      IF ( LPRT ) CALL DEBUG_MSG( '### GEOS5_NESTED: a array flip' )

      ! Convert tracer concentration from mol/mol moist air to 
      ! mol/mol dry air after advection (ewl, 4/16/15)
      CALL Convert_DryVV_to_MoistVV( am_I_Root, N_TRACERS, State_Met, 
     &                               State_Chm, RC )

      ! Note: the mass flux diagnostic arrays (MASSFLEW, MASSFLNS and MASSFLUP)
      ! are incremented upside-down (level 1 = top of the atmosphere).
      ! The levels order is reversed only when written out in diag3.f
      ! (ccc, 3/8/10)
      ! Do the advection
      CALL TPCORE_GEOS5_WINDOW( D_DYN,    Re,        IIPAR,    JJPAR,
     &                          LLPAR,    JFIRST,    JLAST,    NG,
     &                          MG,       N_TOT_TRC, Ap,       Bp,
     &                          p_UWND,   p_VWND,    P_TP1,    P_TP2,
     &                          P_TEMP,   p_STT,     IORD,     JORD,
     &                          KORD,     N_ADJ,     p_XMASS,  p_YMASS,
     &                          MASSFLEW, MASSFLNS,  MASSFLUP, A_M2,
     &                          TCVV,     ND24,      ND25,     ND26,  
     &                          p_MOISTMW )

      ! Convert tracer concentration from mol/mol moist air to 
      ! mol/mol dry air after advection (ewl, 4/16/15)
      CALL Convert_MoistVV_to_DryVV( am_I_Root, N_TRACERS, State_Met, 
     &                               State_Chm, RC ) 

      ! Free pointer memory
      NULLIFY( p_UWND, p_VWND, p_STT, p_XMASS, p_YMASS, p_MOISTMW ) 

      !=================================================================
      ! Reset surface pressure and ensure mass conservation
      !=================================================================

      ! Reset the floating surface pressure with P_TP2, the "true"
      ! surface pressure at the end of the dynamic timestep.
      CALL SET_FLOATING_PRESSURE( P_TP2 )

      ! Update State_Met air quantities with new pressure (ewl, 3/31/15)
      CALL AIRQNT( am_I_Root, State_Met, RC )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### G5_NESTED_ADV: a TPCORE' )

      END SUBROUTINE DO_GEOS5_WINDOW_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: do_geosfp_window_transport
!
! !DESCRIPTION: Subroutine DO\_GEOSFP\_WINDOW\_TRANSPORT is the driver program 
!  for the proper TPCORE program for the GEOS-5 nested-grid simulations. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE DO_GEOSFP_WINDOW_TRANSPORT( am_I_Root, Input_Opt,
     &                                       State_Met, State_Chm, RC )
!
! !USES:
!
      ! References to F90 modules
      USE DIAG_MOD,                  ONLY : MASSFLEW, MASSFLNS, MASSFLUP
      USE ERROR_MOD,                 ONLY : IT_IS_NAN,     DEBUG_MSG
      USE GRID_MOD,                  ONLY : GET_XOFFSET,   GET_YOFFSET
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod,        ONLY : OptInput
      USE GIGC_State_Chm_Mod,        ONLY : ChmState
      USE GIGC_State_Met_Mod,        ONLY : MetState
      USE PJC_PFIX_GEOSFP_WINDOW_MOD,ONLY : DO_PJC_PFIX_GEOSFP_WINDOW
      USE PRESSURE_MOD,              ONLY : SET_FLOATING_PRESSURE
      USE DAO_MOD,                   ONLY : AIRQNT
      USE TIME_MOD,                  ONLY : GET_TS_DYN
      USE TPCORE_BC_MOD,             ONLY : I0_W, J0_W, I1_W, J1_W
      USE TPCORE_BC_MOD,             ONLY : I2_W, J2_W, IM_W, JM_W, IGZD
      USE TPCORE_BC_MOD,             ONLY : DO_WINDOW_TPCORE_BC 
      USE TPCORE_GEOSFP_WINDOW_MOD,  ONLY : TPCORE_GEOSFP_WINDOW
      USE UNITCONV_MOD
#if   defined( APM )
      USE TRACER_MOD,                ONLY : N_APMTRA
#endif

      USE CMN_SIZE_MOD                    ! Size parameters
      USE CMN_DIAG_MOD                    ! NDxx flags
      USE CMN_GCTM_MOD                    ! Physical constants
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)    :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)    :: Input_Opt   ! Input Options object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
      TYPE(MetState), INTENT(INOUT) :: State_Met   ! Meteorology State object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT)   :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  10 Mar 2003 - R. Yantosca - Initial version

!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  16 Feb 2011 - R. Yantosca - Add modifications for APM from G. Luo
!  21 Jun 2012 - R. Yantosca - Now use pointers to flip arrays in the vertical
!  04 Jun 2013 - R. Yantosca - Replace pointer STT with State_Chm%Tracers
!  26 Feb 2015 - E. Lundgren - Replace GET_PEDGE with State_Met%PEDGE.
!  31 Mar 2015 - E. Lundgren - Move AIRQNT call from main to within transport
!  01 Apr 2015 - L. Zhang    - Add updates to skip the window region
!  15 Apr 2015 - R. Yantosca - Add more pointers to avoid array temporaries
!  16 Apr 2015 - E. Lundgren - Convert mixing ratio to moist before advection
!                              and back to dry afterwards
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER            :: I0,    J0
      INTEGER            :: IM_W1, JM_W1, I0_W1, J0_W1, BUFF_SIZE
      INTEGER            :: I, J, L, L2, N, N_DYN, N_TOT_TRC
      REAL(fp)           :: A_DIFF, D_DYN, TR_DIFF
      REAL(fp)           :: AD_A(IIPAR,JJPAR,LLPAR)
      REAL(fp)           :: AD_B(IIPAR,JJPAR,LLPAR)
      REAL(fp),  TARGET  :: P_TP1(IIPAR,JJPAR)
      REAL(fp),  TARGET  :: P_TP2(IIPAR,JJPAR)
      REAL(fp),  TARGET  :: P_TEMP(IIPAR,JJPAR)
      REAL(fp)           :: TR_A(IIPAR,JJPAR,LLPAR)
#if   defined( APM )
      REAL(fp)           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                      Input_Opt%N_TRACERS+N_APMTRA)
#else
      REAL(fp)           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                           Input_Opt%N_TRACERS)
#endif
      REAL(fp),  TARGET  :: XMASS(IIPAR,JJPAR,LLPAR)
      REAL(fp),  TARGET  :: YMASS(IIPAR,JJPAR,LLPAR)

      ! Pointers
      REAL(fp),  POINTER :: p_A_M2  (  :    ) => NULL()
      REAL(fp),  POINTER :: p_P_TP1 (:,:    ) => NULL()
      REAL(fp),  POINTER :: p_P_TP2 (:,:    ) => NULL()
      REAL(fp),  POINTER :: p_P_TEMP(:,:    ) => NULL()
      REAL(fp),  POINTER :: p_UWND  (:,:,:  ) => NULL()
      REAL(fp),  POINTER :: p_VWND  (:,:,:  ) => NULL()
      REAL(fp),  POINTER :: p_XMASS (:,:,:  ) => NULL()
      REAL(fp),  POINTER :: p_YMASS (:,:,:  ) => NULL()
      REAL(fp),  POINTER :: p_MOISTMW(:,:,: ) => NULL()
      REAL(fp),  POINTER :: p_MFLEW (:,:,:,:) => NULL()
      REAL(fp),  POINTER :: p_MFLNS (:,:,:,:) => NULL()
      REAL(fp),  POINTER :: p_MFLUP (:,:,:,:) => NULL()
      REAL(fp),  POINTER :: p_STT   (:,:,:,:) => NULL()

      ! For fields from Input_Opt
      LOGICAL            :: LFILL
      LOGICAL            :: LPRT
      LOGICAL            :: LWINDO
      INTEGER            :: N_TRACERS
      REAL(fp)           :: TCVV(Input_Opt%N_TRACERS)

      !=================================================================
      ! DO_GEOSFP_FVDAS_WINDOW_TRANSPORT begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LFILL     = Input_Opt%LFILL
      LPRT      = Input_Opt%LPRT
      LWINDO    = Input_Opt%LWINDO
      N_TRACERS = Input_Opt%N_TRACERS
      TCVV      = Input_Opt%TCVV(1:N_TRACERS)

      ! Pick the proper # of tracers
#if   defined( APM )
      N_TOT_TRC = N_TRACERS + N_APMTRA   ! G-C advected tracers + APM tracers
#else
      N_TOT_TRC = N_TRACERS              ! G-C advected tracers only
#endif

      ! Get nested-grid lon/lat offsets [# boxes]
      I0    = GET_XOFFSET( GLOBAL=.TRUE. )
      J0    = GET_YOFFSET( GLOBAL=.TRUE. )

      ! Dynamic timestep [s]
      N_DYN = GET_TS_DYN() * 60
      D_DYN = DBLE( N_DYN )

      ! (lzh, 09/01/2014)
      BUFF_SIZE = 2
      IM_W1     = IM_W + 2 * BUFF_SIZE
      JM_W1     = JM_W + 2 * BUFF_SIZE
      I0_W1     = I0_W     - BUFF_SIZE
      J0_W1     = J0_W     - BUFF_SIZE

      !=================================================================
      ! Prepare variables for calls to PJC pressure-fixer and TPCORE
      !
      ! For GEOS-4 and GEOS-5 (hybrid grids), the pressure at the
      ! bottom edge of grid box (I,J,L) is given by:
      !
      !    P(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
      !
      ! where Psurface is the true surface pressure (i.e. not PS-PTOP).
      ! and Ap(L), Bp(L) define the vertical grid (see pressure_mod.f)
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! True surface pressure at midpoint of dynamic timestep [hPa]
         P_TP1(I,J) = State_Met%PEDGE(I,J,1)

         ! True surface pressure at end of dynamic timestep [hPa]
         P_TP2(I,J) = State_Met%PSC2(I,J)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Get the air and tracer mass before advection
      !=================================================================

      ! Dry airmass [kg] before transport
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         AD_B(I,J,L) = GET_AIR_MASS( I, J, L, P_TP1(I,J), State_Met )
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Tracer mass [kg] before transport
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
      DO N = 1, N_TOT_TRC
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         TR_B(I,J,L,N) = State_Chm%Tracers(I,J,L,N)
     &                 * AD_B(I,J,L)
     &                 / TCVV(N)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Call the PJC/LLNL pressure fixer to get the adjusted air
      ! masses XMASS and YMASS.  XMASS and YMASS need to be passed to
      ! TPCORE_FVDAS in order to ensure mass conservation.
      !=================================================================
      XMASS = 0e+0_fp !(dan)
      YMASS = 0e+0_fp
      ! NOTE: P_TP1 and P_TP2 are the true surface pressures!  
      CALL DO_PJC_PFIX_GEOSFP_WINDOW( D_DYN,
     &                                P_TP1,        P_TP2,
     &                                State_Met%U,  State_Met%V,
     &                                XMASS,        YMASS )

      IF ( LPRT ) CALL DEBUG_MSG( '### GEOS5_NESTED: a GET_TRA_MASS' )

      ! Impose TPCORE boundary conditions @ edges of nested grid
      CALL DO_WINDOW_TPCORE_BC( am_I_Root, Input_Opt, State_Chm, RC )

      ! Flip array indices in the vertical using pointer storage

      ! Exclude the buffer zone (lzh, 4/1/2015)
      p_UWND  => State_Met%U(       I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1,
     &                              LLPAR   : 1            : -1    )

      p_VWND  => State_Met%V(       I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1, 
     &                              LLPAR   : 1            : -1    )

      p_STT   => State_Chm%Tracers( I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1,
     &                              LLPAR   : 1            : -1, : )

      p_XMASS => XMASS(             I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1,
     &                              LLPAR   : 1            : -1    )

      p_YMASS => YMASS(             I0_W1+1 : I0_W1+IM_W1, 
     &                              J0_W1+1 : J0_W1+JM_W1,
     &                              LLPAR   : 1            : -1    )

      p_MOISTMW => State_Met%MOISTMW( I0_W1+1 : I0_W1+IM_W1, 
     &                                J0_W1+1 : J0_W1+JM_W1,
     &                                LLPAR   : 1            : -1  )


      p_P_TP1 => P_TP1(             I0_W1+1 : I0_W1+IM_W1, 
     &                              J0_W1+1 : J0_W1+JM_W1          )

      p_P_TP2 => P_TP2(             I0_W1+1 : I0_W1+IM_W1, 
     &                              J0_W1+1 : J0_W1+JM_W1          )

      p_P_TEMP => P_TEMP(           I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1          )

      p_MFLEW => MASSFLEW(          I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1, :,    : )

      p_MFLNS => MASSFLNS(          I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1, :,    : )

      p_MFLUP => MASSFLUP(          I0_W1+1 : I0_W1+IM_W1,
     &                              J0_W1+1 : J0_W1+JM_W1, :,    : )

      p_A_M2  => A_M2(              J0_W1+1 : J0_W1+JM_W1          )


      ! Convert tracer concentration from mol/mol moist air to 
      ! mol/mol dry air after advection (ewl, 4/16/15)
      CALL Convert_DryVV_to_MoistVV( am_I_Root, N_TRACERS, State_Met, 
     &                               State_Chm, RC )

      ! Note: the mass flux diagnostic arrays (MASSFLEW, MASSFLNS and MASSFLUP)
      ! are incremented upside-down (level 1 = top of the atmosphere).
      ! The levels order is reversed only when written out in diag3.f
      ! (ccc, 3/8/10)
      ! Do the advection
!------------------------------------------------------------------------------
! Prior to 4/15/15:
! Now use pointers to avoid array temporaries (bmy, 4/15/15)
!      CALL TPCORE_GEOSFP_WINDOW( D_DYN,    Re,        IIPAR,    JJPAR,
!     &                           LLPAR,    JFIRST,    JLAST,    NG,
!     &                           MG,       N_TOT_TRC, Ap,       Bp,
!     &                           p_UWND,   p_VWND,    P_TP1,    P_TP2,
!     &                           P_TEMP,   p_STT,     IORD,     JORD,      
!     &                           KORD,     N_ADJ,     p_XMASS,  p_YMASS,
!     &                           MASSFLEW, MASSFLNS,  MASSFLUP, A_M2,
!     &                           TCVV,     ND24,      ND25,     ND26 )
!------------------------------------------------------------------------------

      CALL TPCORE_GEOSFP_WINDOW( 
     &      D_DYN,   Re,       IM_W1,   JM_W1,   LLPAR, 
     &      JFIRST,  JLAST,    NG,      MG,      N_TOT_TRC, 
     &      Ap,      Bp,       p_UWND,  p_VWND,  p_P_TP1, 
     &      p_P_TP2, p_P_TEMP, p_STT,   IORD,    JORD, 
     &      KORD,    N_ADJ,    p_XMASS, p_YMASS, p_MFLEW, 
     &      p_MFLNS, p_MFLUP,  p_A_M2,  TCVV,    ND24,     
     &      ND25,    ND26,     p_MOISTMW ) 

      ! Convert tracer concentration from mol/mol moist air to 
      ! mol/mol dry air after advection (ewl, 4/16/15)
      CALL Convert_MoistVV_to_DryVV( am_I_Root, N_TRACERS, State_Met, 
     &                               State_Chm, RC ) 

      ! Free pointer memory
      NULLIFY( p_UWND,  p_VWND,  p_STT,   p_XMASS, p_MOISTMW ) 
      NULLIFY( p_YMASS, p_P_TP1, p_P_TP2, p_P_TEMP )
      NULLIFY( p_MFLEW, p_MFLNS, p_MFLUP, p_A_M2   )

      !=================================================================
      ! Reset surface pressure and ensure mass conservation
      !=================================================================
      ! Reset the floating surface pressure with P_TP2, the "true"
      ! surface pressure at the end of the dynamic timestep.
      CALL SET_FLOATING_PRESSURE( P_TP2 )

      ! Update State_Met air quantities with new pressure (ewl, 3/31/15)
      CALL AIRQNT( am_I_Root, State_Met, RC )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### G5_NESTED_ADV: a TPCORE' )

      END SUBROUTINE DO_GEOSFP_WINDOW_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_air_mass
!
! !DESCRIPTION: Function GET\_AIR\_MASS returns the dry air mass based on the 
!  pressures returned before and after the call to the GEOS-4/fvDAS TPCORE 
!  code. (bmy, 6/24/03)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_AIR_MASS( I, J, L, P_SURF, State_Met ) 
     &                       RESULT( AIR_MASS )
!
! !USES:
!
      USE CMN_SIZE_MOD               ! Size parameters
      USE CMN_GCTM_MOD               ! g0_100
      USE GIGC_State_Met_Mod, ONLY : MetState
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I, J, L   ! GEOS-Chem lon, lat, level indices
      REAL(fp),  INTENT(IN) :: P_SURF    ! Surface pressure [hPa] at (I,J,L=1)
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object

! 
! !REVISION HISTORY: 
!  24 Jun 2003 - R. Yantosca - Initial version
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  26 Mar 2015 - E. Lundgren - Calculate dry rather than moist air mass
!  30 Apr 2015 - E. Lundgren - Use 1-SPHU for moist to dry conversion
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: L2
      REAL(fp)  :: P_BOT, P_TOP, AIR_MASS
      
      !=================================================================
      ! GET_AIR_MASS begins here!
      !=================================================================
      
      ! Index for Ap, Bp from atmosphere top down to surface
      ! since the Ap's and Bp's have been flipped for TPCORE
      L2       = ( LLPAR + 1 ) - L + 1
               
      ! Hybrid-grid formulation for air mass. 
      ! SPHU is defined as mass water vapor / mass moist air, and 
      ! therefore (1-SPHU) is mass dry air / mass moist air.
      ! Note that State_Met%SPHU is in [g/kg].
      P_BOT    = Ap(L2)   + ( Bp(L2)   * P_SURF )
      P_TOP    = Ap(L2-1) + ( Bp(L2-1) * P_SURF )
      AIR_MASS = ( P_BOT - P_TOP ) * G0_100 * A_M2(J)
     &             * ( 1.0e+0_fp - State_Met%SPHU(I,J,L) * 1.e-3_fp )

      END FUNCTION GET_AIR_MASS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_transport
!
! !DESCRIPTION: Subroutine SET\_TRANSPORT passes IORD, JORD, KORD values 
!  from "input\_mod.f".
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE SET_TRANSPORT( I_ORD, J_ORD, K_ORD )
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I_ORD  ! IORD option for E/W advection
      INTEGER, INTENT(IN) :: J_ORD  ! JORD option for N/S advection
      INTEGER, INTENT(IN) :: K_ORD  ! KORD option for vertical diffusion
! 
! !REVISION HISTORY: 
!  20 Jul 2004 - R. Yantosca - Initial version
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
      IORD = I_ORD
      JORD = J_ORD
      KORD = K_ORD 

      END SUBROUTINE SET_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_transport
!
! !DESCRIPTION: Subroutine INIT\_TRANSPORT initializes all module variables 
!  and arrays.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_TRANSPORT( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ERROR_MOD,          ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GRID_MOD,           ONLY : GET_AREA_M2, GET_YMID_R
      USE PRESSURE_MOD,       ONLY : GET_AP,      GET_BP
      USE TIME_MOD,           ONLY : GET_TS_DYN
      USE TPCORE_FVDAS_MOD,   ONLY : INIT_TPCORE

      USE CMN_SIZE_MOD         ! Size parameters
      USE CMN_GCTM_MOD         ! Re
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  10 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Now references GET_TS_DYN from "time_mod.f", INIT_TPCORE_FVDAS from
!        "tpcore_fvdas_mod.f90", and GET_YMID_R from "grid_mod.f".  Now also
!        include "CMN_SETUP".  (bdf, bmy, 4/28/03)
!  (2 ) Remove reference to DSIG, it's obsolete. (bmy, 6/24/03)
!  (3 ) Now references LEMBED & LTPFV from "logical_mod.f".  Now references
!        N_TRACERS from "tracer_mod.f". (bmy, 7/20/04)
!  (4 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
!  (5 ) Removed reference to USE_GEOS_4_TRANSPORT, STT_I1, STT_I2, STT_J1,
!        STT_J2, variables (bmy, 10/30/07)
!  (6 ) Deleted reference to CMN, it's not needed anymore (bmy, 11/6/08)
!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  28 Feb 2012 - R. Yantosca - Removed support for GEOS-3
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_M2(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID_R_W(I,J,L) from grid_mod.F90
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  05 Dec 2014 - M. Yannetti - Added REAL_N_DYN
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: J, K, L, N_DYN
      REAL(fp)  :: YMID_R(JJPAR)
      REAL(fp)  :: REAL_N_DYN

      ! For fields from Input_Opt
      LOGICAL :: LTRAN
      INTEGER :: N_TRACERS

      !=================================================================
      ! Initialize
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LTRAN     = Input_Opt%LTRAN
      N_TRACERS = Input_Opt%N_TRACERS

      !=================================================================
      ! Allocate arrays for TPCORE vertical coordinates 
      !
      ! For TPCORE v7.1.m (for GEOS-3 met fields):
      ! 
      !    P(I,J,L) = ( Ap(L) * PTOP ) + ( Bp(L) * ( Psurf(I,J)-PTOP ) )
      !
      ! For fvDAS TPCORE with for GEOS-4 or GEOS-5 met fields:
      !
      !    P(I,J,L) = Ap(L) + ( Bp(L) * Psurf(I,J) )
      !
      ! Also here Ap, Bp will be flipped since both TPCORE versions
      ! index levels from the atm. top downwards (bdf, bmy, 10/30/07)
      !=================================================================
      ALLOCATE( Ap( LLPAR+1 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Ap' )
 
      ALLOCATE( Bp( LLPAR+1 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Bp' )

      ! Flip Ap and Bp for TPCORE
      DO L = 1, LLPAR+1 

         ! As L runs from the surface up, 
         ! K runs from the top down
         K = ( LLPAR + 1 ) - L + 1

         Ap(L) = GET_AP(K)          ! Ap(L) is in [hPa] 
         Bp(L) = GET_BP(K)
      ENDDO

      !=================================================================
      ! Allocate arrays for surface area and layer thickness
      !=================================================================
      ALLOCATE( A_M2( JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'A_M2' )

      ! Surface area [m2]
      DO J = 1, JJPAR
         A_M2(J) = GET_AREA_M2( 1, J, 1 )
      ENDDO

      !=================================================================
      ! Additional setup for the GEOS-4/fvDAS version of TPCORE
      !=================================================================

      ! Initialize
      N_DYN = GET_TS_DYN() * 60
      N_ADJ = 0
      NG    = 0
      MG    = 0

      ! YMID_R is latitude of grid box center [radians]
      DO J = 1,JJPAR
         YMID_R(J) = GET_YMID_R( 1, J, 1 )
      ENDDO

      REAL_N_DYN = N_DYN

      ! Call INIT routine from "tpcore_fvdas_mod.f" 
      CALL INIT_TPCORE( IIPAR,  JJPAR, LLPAR,  JFIRST, JLAST, 
     &                  NG, MG, REAL_N_DYN, Re,     YMID_R )

      END SUBROUTINE INIT_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_geos5_window_transport

!
! !DESCRIPTION: Subroutine INIT\_GEOS5\_WINDOW\_TRANSPORT initializes all 
!  module variables and arrays for the GEOS-5 nested grid simulation.  
!  This routine is only called if we are using the GEOS-5 nested grid 
!  simulation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GEOS5_WINDOW_TRANSPORT( am_I_Root, Input_Opt, RC )
!
! !USES:
!
      USE ERROR_MOD,               ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod,      ONLY : OptInput
      USE GRID_MOD,                ONLY : GET_AREA_M2
      USE GRID_MOD,                ONLY : GET_YMID_R_W
      USE PRESSURE_MOD,            ONLY : GET_AP, GET_BP
      USE TIME_MOD,                ONLY : GET_TS_DYN
      USE TPCORE_FVDAS_MOD,        ONLY : INIT_TPCORE
      USE TPCORE_BC_MOD,           ONLY : I0_W, J0_W, I1_W, J1_W
      USE TPCORE_BC_MOD,           ONLY : I2_W, J2_W, IM_W, JM_W
      USE TPCORE_BC_MOD,           ONLY : IGZD, INIT_TPCORE_BC
      USE TPCORE_GEOS5_WINDOW_MOD, ONLY : INIT_GEOS5_WINDOW

      USE CMN_SIZE_MOD                ! Size parameters
      USE CMN_GCTM_MOD                ! Re
!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
! 
! !REVISION HISTORY: 
!  06 Jun 2008 - D. Chen & R. Yantosca - Initial version
!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID_R_W(I,J,L) from grid_mod.F90
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, RC
!  05 Dec 2014 - M. Yannetti - Added REAL_N_DYN
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL  :: LTRAN
      INTEGER  :: N_TRACERS
      INTEGER  :: J, K, L, N_DYN
      REAL(fp) :: REAL_N_DYN

      ! Arrays
      REAL(fp) :: YMID_R_W(0:JJPAR+1)


      !=================================================================
      ! Initialize
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LTRAN     = Input_Opt%LTRAN
      N_TRACERS = Input_Opt%N_TRACERS

      !=================================================================
      ! Allocate arrays for TPCORE vertical coordinates
      ! GEOS-5 nested grid simulation only!!!
      !
      ! For fvDAS TPCORE with for GEOS-5 met fields:
      !
      !    P(I,J,L) = Ap(L) + ( Bp(L) * Psurf(I,J) )
      !
      ! Also here Ap, Bp will be flipped since both TPCORE versions
      ! index levels from the atm. top downwards (bdf, bmy, 10/30/07)
      !=================================================================
      ALLOCATE( Ap( LLPAR+1 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Ap' )

      ALLOCATE( Bp( LLPAR+1 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Bp' )

      ! Flip Ap and Bp for TPCORE
      DO L = 1, LLPAR+1

         ! As L runs from the surface up,
         ! K runs from the top down
         K = ( LLPAR + 1 ) - L + 1

         Ap(L) = GET_AP(K)
         Bp(L) = GET_BP(K)
      ENDDO

      !=================================================================
      ! Allocate arrays for surface area and layer thickness
      !=================================================================
      ALLOCATE( A_M2( JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'A_M2' )

      ! Surface area [m2]
      DO J = 1, JJPAR
         A_M2(J) = GET_AREA_M2( 1, J, 1 )
      ENDDO

      !=================================================================
      ! Additional setup for the GEOS-4/fvDAS version of TPCORE
      !====================f=============================================

      ! Initialize
      N_DYN = GET_TS_DYN() * 60
      N_ADJ = 0
      NG    = 0
      MG    = 0

      ! YMID_R is latitude of grid box center [radians]
      DO J =0, JJPAR+1
         YMID_R_W(J) = GET_YMID_R_W( 1, J, 1 )
      ENDDO

      REAL_N_DYN = N_DYN

      ! Call INIT routine from "tpcore_geos5_fvdas_window_mod.f"
      CALL INIT_GEOS5_WINDOW( IIPAR, JJPAR,    LLPAR, JFIRST, 
     &                        JLAST, NG,       MG,    REAL_N_DYN , 
     &                        Re,    YMID_R_W )

      END SUBROUTINE INIT_GEOS5_WINDOW_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_geosfp_window_transport

!
! !DESCRIPTION: Subroutine INIT\_GEOSFP\_WINDOW\_TRANSPORT initializes all 
!  module variables and arrays for the GEOS-57 nested grid simulation.  
!  This routine is only called if we are using the GEOS-57 nested grid 
!  simulation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GEOSFP_WINDOW_TRANSPORT( am_I_Root,Input_Opt,RC )
!
! !USES:    
!
      USE ERROR_MOD,               ONLY : ALLOC_ERR
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod,      ONLY : OptInput
      USE GRID_MOD,                ONLY : GET_AREA_M2
      USE GRID_MOD,                ONLY : GET_YMID_R_W
      USE PRESSURE_MOD,            ONLY : GET_AP, GET_BP
      USE TIME_MOD,                ONLY : GET_TS_DYN
      USE TPCORE_FVDAS_MOD,        ONLY : INIT_TPCORE
      USE TPCORE_BC_MOD,           ONLY : I0_W, J0_W, I1_W, J1_W
      USE TPCORE_BC_MOD,           ONLY : I2_W, J2_W, IM_W, JM_W
      USE TPCORE_BC_MOD,           ONLY : IGZD, INIT_TPCORE_BC
      USE TPCORE_GEOSFP_WINDOW_MOD,ONLY : INIT_GEOSFP_WINDOW

      USE CMN_SIZE_MOD                ! Size parameters
      USE CMN_GCTM_MOD                ! Re

!
! !INPUT PARAMETERS:
!
      LOGICAL,        INTENT(IN)  :: am_I_Root   ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt   ! Input Options object
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC          ! Success or failure?
!
! !REVISION HISTORY: 
!  06 Jun 2008 - D. Chen & R. Yantosca - Initial version
!  26 Feb 2010 - R. Yantosca - Removed references to obsolete LEMBED switch
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!  01 Mar 2012 - R. Yantosca - Now use GET_AREA_CM2(I,J,L) from grid_mod.F90
!  01 Mar 2012 - R. Yantosca - Now use GET_YMID_R_W(I,J,L) from grid_mod.F90
!  25 Mar 2013 - R. Yantosca - Now accept am_I_Root, Input_Opt, State_Chm, RC
!  26 Sep 2013 - R. Yantosca - Renamed to INIT_GEOSFP_WINDOW_TRANSPORT
!  05 Dec 2014 - M. Yannetti - Added REAL_N_DYN
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      ! Scalars
      LOGICAL  :: LTRAN
      INTEGER  :: N_TRACERS
      INTEGER  :: BUFF_SIZE
      INTEGER  :: J,     K,     L,     N_DYN
      INTEGER  :: IM_W1, JM_W1, I0_W1, J0_W1
      REAL(fp) :: REAL_N_DYN

      ! Arrays
      REAL(fp) :: YMID_R_W(0:JJPAR+1)

      !=================================================================
      ! Initialize
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LTRAN     = Input_Opt%LTRAN
      N_TRACERS = Input_Opt%N_TRACERS

      ! Cast N_DYN to flexible precision
#if defined( USE_REAL8 )
      REAL_N_DYN = DBLE( N_DYN )
#else
      REAL_N_DYN = FLOAT( N_DYN )
#endif

      !=================================================================
      ! Allocate arrays for TPCORE vertical coordinates
      ! GEOS-5 nested grid simulation only!!!
      !
      ! For fvDAS TPCORE with for GEOS-5 met fields:
      !
      !    P(I,J,L) = Ap(L) + ( Bp(L) * Psurf(I,J) )
      !
      ! Also here Ap, Bp will be flipped since both TPCORE versions
      ! index levels from the atm. top downwards (bdf, bmy, 10/30/07)
      !=================================================================
      ALLOCATE( Ap( LLPAR+1 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Ap' )

      ALLOCATE( Bp( LLPAR+1 ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'Bp' )

      ! Flip Ap and Bp for TPCORE
      DO L = 1, LLPAR+1

         ! As L runs from the surface up,
         ! K runs from the top down
         K = ( LLPAR + 1 ) - L + 1

         Ap(L) = GET_AP(K)
         Bp(L) = GET_BP(K)
      ENDDO

      !=================================================================
      ! Allocate arrays for surface area and layer thickness
      !=================================================================
      ALLOCATE( A_M2( JJPAR ), STAT=RC )
      IF ( RC /= 0 ) CALL ALLOC_ERR( 'A_M2' )

      ! Surface area [m2]
      DO J = 1, JJPAR
         A_M2(J) = GET_AREA_M2( 1, J, 1 )
      ENDDO

      !=================================================================
      ! Additional setup for the GEOS-4/fvDAS version of TPCORE
      !=================================================================

      ! Initialize
      N_DYN = GET_TS_DYN() * 60
      N_ADJ = 0
      NG    = 0
      MG    = 0

      ! (lzh, 4/1/2015)
      BUFF_SIZE = 2
      IM_W1     = IM_W + 2 * BUFF_SIZE
      JM_W1     = JM_W + 2 * BUFF_SIZE
      I0_W1     = I0_W - BUFF_SIZE
      J0_W1     = J0_W - BUFF_SIZE

      ! YMID_R is latitude of grid box center [radians]
      DO J =0, JJPAR+1
         YMID_R_W(J) = GET_YMID_R_W( 1, J, 1 )
      ENDDO

      REAL_N_DYN = N_DYN

      ! Call INIT routine from "tpcore_geosfp_fvdas_window_mod.f"
      ! subroutine "init_geosfp_windows" is in "tpcore_geosfp_window_mod.F90"
!------------------------------------------------------------------------------
! Prior to 4/1/15:
!      CALL INIT_GEOSFP_WINDOW( IIPAR, JJPAR,    LLPAR, JFIRST, 
!     &                         JLAST, NG,       MG,    REAL_N_DYN ,
!     &                         Re,    YMID_R_W )
!------------------------------------------------------------------------------
      CALL INIT_GEOSFP_WINDOW( IM_W1, 
     &                         JM_W1,    
     &                         LLPAR, 
     &                         JFIRST, 
     &                         JLAST, 
     &                         NG,       
     &                         MG,    
!------------------------------------------------------------------------------
! Prior to 4/23/15:
! Now use REAL_N_DYN instead of DBLE( N_DYN ), which causes problems
! when compiling with PRECISION=4 (bmy, 4/23/15)
!     &                         DBLE( N_DYN ), 
!------------------------------------------------------------------------------
     &                         REAL_N_DYN,
     &                         Re,  
     &                         YMID_R_W( J0_W1:(J0_W1+JM_W1+1) ) )
      END SUBROUTINE INIT_GEOSFP_WINDOW_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: cleanup_transport
!
! !DESCRIPTION: Subroutine CLEANUP\_TRANSPORT deallocates all module arrays. 
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP_TRANSPORT
! 
! !REVISION HISTORY: 
!  10 Mar 2003 - R. Yantosca - Initial version
!  (1 ) Remove reference to DSIG, it's obsolete. (bmy, 6/24/03)
!  (2 ) Remove obsolete embedded chemistry arrays (bmy, 10/30/07)
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
      IF ( ALLOCATED( Ap     ) ) DEALLOCATE( Ap     )
      IF ( ALLOCATED( A_M2   ) ) DEALLOCATE( A_M2   )
      IF ( ALLOCATED( Bp     ) ) DEALLOCATE( Bp     )

      END SUBROUTINE CLEANUP_TRANSPORT
!EOC
      END MODULE TRANSPORT_MOD
