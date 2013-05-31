!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      IMPLICIT NONE
#     include "define.h"
      PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
      PUBLIC  :: CLEANUP_TRANSPORT
      PUBLIC  :: DO_TRANSPORT
      PUBLIC  :: INIT_TRANSPORT
      PUBLIC  :: INIT_GEOS5_WINDOW_TRANSPORT
      PUBLIC  :: INIT_GEOS57_WINDOW_TRANSPORT
      PUBLIC  :: SET_TRANSPORT
!
! !PRIVATE MEMBER FUNCTIONS:
!
      PRIVATE :: GEOS4_GEOS5_GLOBAL_ADV      
      PRIVATE :: GCAP_GLOBAL_ADV             
      PRIVATE :: DO_GEOS5_WINDOW_TRANSPORT
      PRIVATE :: DO_GEOS57_WINDOW_TRANSPORT      
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
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! MODULE VARIABLES:
      !
      ! (1 ) Ap     (REAL*8 ) : Vertical coordinate array for TPCORE
      ! (2 ) A_M2   (REAL*8 ) : Grid box surface areas [m2]
      ! (3 ) Bp     (REAL*8 ) : Vertical coordinate array for TPCORE
      ! (4 ) IORD   (REAL*8 ) : TPCORE E/W option flag
      ! (5 ) JORD   (REAL*8 ) : TPCORE N/S option flag
      ! (6 ) KORD   (REAL*8 ) : TPCORE vertical option flag
      ! (7 ) JLAST  (INTEGER) : For fvDAS TPCORE
      ! (8 ) MG     (INTEGER) : For fvDAS TPCORE
      ! (9 ) NG     (INTEGER) : For fvDAS TPCORE
      ! (10) N_ADJ  (INTEGER) : For fvDAS TPCORE
      !=================================================================
      INTEGER             :: IORD,  JORD, KORD, JFIRST 
      INTEGER             :: JLAST, NG,   MG,   N_ADJ
      REAL*8, ALLOCATABLE :: Ap(:)
      REAL*8, ALLOCATABLE :: A_M2(:)
      REAL*8, ALLOCATABLE :: Bp(:)

      ! This seems to be obsolete ???
      !REAL*8, ALLOCATABLE :: STT_BC2(:,:,:)

      CONTAINS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
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

#elif defined( GRID025x03125)

      ! GEOS-5.7.2 (aka "FP") nested grids
      CALL INIT_GEOS57_WINDOW_TRANSPORT
     &     ( am_I_Root, Input_Opt, RC )    
#else

      ! All global simulations
      CALL INIT_TRANSPORT                   
     &     ( am_I_Root, Input_Opt, RC )    

#endif
         FIRST = .FALSE.
      ENDIF

      ! Choose the proper version of TPCORE for the nested-grid window 
      ! region (usually 1x1 grids) or for the entire globe
      IF ( ITS_A_NESTED_GRID() ) THEN
         CALL DO_GEOS5_WINDOW_TRANSPORT
     &      ( am_I_Root, Input_Opt, State_Met, State_Chm, RC )
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      USE DIAG_MOD,           ONLY : MASSFLEW, MASSFLNS, MASSFLUP
      USE ERROR_MOD,          ONLY : IT_IS_NAN, DEBUG_MSG, SAFE_DIV
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PJC_PFIX_MOD,       ONLY : DO_PJC_PFIX
      USE PRESSURE_MOD,       ONLY : GET_PEDGE, SET_FLOATING_PRESSURE
      USE TIME_MOD,           ONLY : GET_TS_DYN
      USE TPCORE_BC_MOD,      ONLY : SAVE_GLOBAL_TPCORE_BC
      USE TPCORE_FVDAS_MOD,   ONLY : TPCORE_FVDAS
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
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER          :: I, J, L, L2, N, N_DYN, N_TOT_TRC
      REAL*8           :: A_DIFF, D_DYN, TR_DIFF
      REAL*8           :: AD_A(IIPAR,JJPAR,LLPAR)
      REAL*8           :: AD_B(IIPAR,JJPAR,LLPAR)
      REAL*8           :: P_TP1(IIPAR,JJPAR)
      REAL*8           :: P_TP2(IIPAR,JJPAR)
      REAL*8           :: P_TEMP(IIPAR,JJPAR)
      REAL*8           :: TR_A(IIPAR,JJPAR,LLPAR)
#if   defined( APM )
      REAL*8           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                         Input_Opt%N_TRACERS+N_APMTRA)  ! With APM
#else
      REAL*8           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                         Input_Opt%N_TRACERS)           ! W/o APM
#endif
      REAL*8,  TARGET  :: XMASS(IIPAR,JJPAR,LLPAR) 
      REAL*8,  TARGET  :: YMASS(IIPAR,JJPAR,LLPAR) 

      ! Pointers
      REAL*8,  POINTER :: p_UWND (:,:,:  )
      REAL*8,  POINTER :: p_VWND (:,:,:  )
      REAL*8,  POINTER :: p_XMASS(:,:,:  )
      REAL*8,  POINTER :: p_YMASS(:,:,:  )
      REAL*8,  POINTER :: p_STT  (:,:,:,:)
      REAL*8,  POINTER :: STT    (:,:,:,:)

      ! Variable to ensure mass conservation (ccc, 2/17/09)
      REAL*8           :: SUMADA

      ! For fields from Input_Opt
      LOGICAL          :: LFILL
      LOGICAL          :: LMFCT
      LOGICAL          :: LPRT
      LOGICAL          :: LWINDO
      INTEGER          :: N_TRACERS
      REAL*8           :: TCVV(Input_Opt%N_TRACERS)

      !=================================================================
      ! GEOS4_GEOS5_GLOBAL_ADV begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LFILL     = Input_Opt%LFILL
      LMFCT     = Input_Opt%LMFCT
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

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

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
         P_TP1(I,J) = GET_PEDGE(I,J,1)

         ! True surface pressure at end of dynamic timestep [hPa]
         P_TP2(I,J) = State_Met%PSC2(I,J)    
   
      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      !=================================================================
      ! Get the air and tracer mass before advection
      !=================================================================

      ! Airmass [kg] before transport
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         AD_B(I,J,L) = GET_AIR_MASS( I, J, L, P_TP1(I,J) )
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
         TR_B(I,J,L,N) = STT(I,J,L,N) * AD_B(I,J,L) / TCVV(N)
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
      p_UWND  => State_Met%U(:,:,LLPAR:1:-1  )
      p_VWND  => State_Met%V(:,:,LLPAR:1:-1  )
      p_STT   => STT  (:,:,LLPAR:1:-1,:)
      p_XMASS => XMASS(:,:,LLPAR:1:-1  )
      p_YMASS => YMASS(:,:,LLPAR:1:-1  )

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
     &                   ND26 )

      ! Free pointer memory
      NULLIFY( p_UWND, p_VWND, p_STT, p_XMASS, p_YMASS )

      !=================================================================
      ! Reset surface pressure and ensure mass conservation
      !=================================================================

      ! Reset the floating surface pressure with P_TP2, the "true"
      ! surface pressure at the end of the dynamic timestep.
      CALL SET_FLOATING_PRESSURE( P_TP2 )

      ! Adjust tracer to correct residual non-conservation of mass
      ! This was changed to be applied only for cells with tracer 
      ! concentration. (ccc, 11/20/08) 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, SUMADA, AD_A, TR_A, TR_DIFF )
      DO N = 1, N_TOT_TRC

         ! Zero summing variable
         SUMADA = 0.d0

         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Air mass [kg] after transport
            AD_A(I,J,L) = GET_AIR_MASS( I, J, L, P_TP2(I,J) )
         
            ! Tracer mass [kg] after transport
            TR_A(I,J,L) = STT(I,J,L,N) * AD_A(I,J,L) / TCVV(N)

            ! We apply mass conservation only on cells w/
            ! nonzero tracer (ccc, 2/17/09)
            IF ( STT(I,J,L,N) > 0.d0 .or. STT(I,J,L,N) < 0.d0 ) THEN
               SUMADA = SUMADA + AD_A(I,J,L)
            ENDIF
         ENDDO
         ENDDO
         ENDDO

         ! Residual mass difference [kg]: before - after
         TR_DIFF = SUM( TR_B(:,:,:,N) ) - SUM( TR_A )

         ! Convert from [kg] to [v/v]
         TR_DIFF = SAFE_DIV(TR_DIFF, SUMADA, 0.d0) * TCVV(N)

         ! Add mass difference [v/v] back to STT
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            IF ( STT(I,J,L,N) > 0.d0 .or. STT(I,J,L,N) < 0.d0 ) THEN
               STT(I,J,L,N) = STT(I,J,L,N) + TR_DIFF
            ENDIF

            STT(I,J,L,N) = MAX( STT(I,J,L,N), 0d0 )
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### G4_G5_GLOB_ADV: a TPCORE' ) 

      END SUBROUTINE GEOS4_GEOS5_GLOBAL_ADV
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      USE DIAG_MOD,           ONLY : MASSFLEW, MASSFLNS, MASSFLUP
      USE ERROR_MOD,          ONLY : IT_IS_NAN, DEBUG_MSG
      USE GIGC_ErrCode_Mod
      USE GIGC_Input_Opt_Mod, ONLY : OptInput
      USE GIGC_State_Chm_Mod, ONLY : ChmState
      USE GIGC_State_Met_Mod, ONLY : MetState
      USE PJC_PFIX_MOD,       ONLY : DO_PJC_PFIX
      USE PRESSURE_MOD,       ONLY : GET_PEDGE, SET_FLOATING_PRESSURE
      USE TIME_MOD,           ONLY : GET_TS_DYN
      USE TPCORE_FVDAS_MOD,   ONLY : TPCORE_FVDAS
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
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
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

!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER          :: I, J, L, L2, N, N_DYN, N_TOT_TRC
      REAL*8           :: A_DIFF, D_DYN, TR_DIFF
      REAL*8           :: AD_A(IIPAR,JJPAR,LLPAR)
      REAL*8           :: AD_B(IIPAR,JJPAR,LLPAR)
      REAL*8           :: P_TP1(IIPAR,JJPAR)
      REAL*8           :: P_TP2(IIPAR,JJPAR)
      REAL*8           :: P_TEMP(IIPAR,JJPAR)
      REAL*8           :: TR_A(IIPAR,JJPAR,LLPAR)
#if   defined( APM )
      REAL*8           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                         Input_Opt%N_TRACERS+N_APMTRA)
#else                  
      REAL*8           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                         Input_Opt%N_TRACERS)
#endif
      REAL*8           :: VTMP(IIPAR,JJPAR,LLPAR)
      REAL*8,  TARGET  :: XMASS(IIPAR,JJPAR,LLPAR) 
      REAL*8,  TARGET  :: YMASS(IIPAR,JJPAR,LLPAR) 

      ! Pointers
      REAL*8,  POINTER :: p_UWND (:,:,:  )
      REAL*8,  POINTER :: p_VWND (:,:,:  )
      REAL*8,  POINTER :: p_XMASS(:,:,:  )
      REAL*8,  POINTER :: p_YMASS(:,:,:  )
      REAL*8,  POINTER :: p_STT  (:,:,:,:)
      REAL*8,  POINTER :: STT    (:,:,:,:)

      ! For fields from Input_Opt
      LOGICAL          :: LFILL
      LOGICAL          :: LMFCT
      LOGICAL          :: LPRT
      LOGICAL          :: LWINDO
      INTEGER          :: N_TRACERS
      REAL*8           :: TCVV(Input_Opt%N_TRACERS)

      ! Variable to ensure mass conservation (ccc, 20/11/08)
      REAL*8           :: SUMADA

      !=================================================================
      ! GCAP_GLOBAL_ADV begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LFILL     = Input_Opt%LFILL
      LMFCT     = Input_Opt%LMFCT
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

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

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
         P_TP1(I,J) = GET_PEDGE(I,J,1) - PTOP

         ! Psurface - PTOP at end of dynamic timestep [hPa]
         P_TP2(I,J) = State_Met%PSC2(I,J) - PTOP

      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      !==============================================================
      ! Get the air & tracer mass before advection
      !==============================================================

      ! Airmass [kg] before transport
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         AD_B(I,J,L) = GET_AIR_MASS( I, J, L, P_TP1(I,J) )
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
         TR_B(I,J,L,N) = STT(I,J,L,N) * AD_B(I,J,L) / TCVV(N)
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
      p_UWND  => State_Met%U(:,:,LLPAR:1:-1  )
      p_VWND  => State_Met%V(:,:,LLPAR:1:-1  )
      p_STT   => STT  (:,:,LLPAR:1:-1,:)
      p_XMASS => XMASS(:,:,LLPAR:1:-1  )
      p_YMASS => YMASS(:,:,LLPAR:1:-1  )

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
     &                   ND26 )

      ! Free pointer memory
      NULLIFY( p_UWND, p_VWND, p_STT, p_XMASS, p_YMASS )

      !=================================================================
      ! Reset surface pressure and ensure mass conservation
      !=================================================================

      ! Reset the floating surface pressure with P_TP2+PTOP, the "true"
      ! surface pressure at the end of the dynamic timestep.
      CALL SET_FLOATING_PRESSURE( P_TP2 + PTOP )

      ! Adjust tracer to correct residual non-conservation of mass
      ! This was changed to be applied only for cells with tracer 
      ! concentration. (ccc, 11/20/08) 
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N, SUMADA, AD_A, TR_A, TR_DIFF )
      DO N = 1, N_TOT_TRC

         ! Zero summing variable
         SUMADA = 0.d0

         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Air mass [kg] after transport
            AD_A(I,J,L) = GET_AIR_MASS( I, J, L, P_TP2(I,J) )
         
            ! Tracer mass [kg] after transport
            TR_A(I,J,L) = STT(I,J,L,N) * AD_A(I,J,L) / TCVV(N)

            ! We apply mass conservation only on cells with 
            ! nonzero tracer. (ccc, 11/20/08)
            IF ( STT(I,J,L,N) > 0.d0 .or. STT(I,J,L,N) < 0.d0 ) THEN
               SUMADA = SUMADA + AD_A(I,J,L)
            ENDIF
         ENDDO
         ENDDO
         ENDDO

         ! Residual mass difference [kg]: before - after
         TR_DIFF = SUM( TR_B(:,:,:,N) ) - SUM( TR_A )

         ! Convert from [kg] to [v/v]
         TR_DIFF = TR_DIFF / SUMADA * TCVV(N)

         ! Add mass difference [v/v] back to STT
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            IF ( STT(I,J,L,N) > 0.d0 .or. STT(I,J,L,N) < 0.d0 ) THEN
               STT(I,J,L,N) = STT(I,J,L,N) + TR_DIFF
            ENDIF

            STT(I,J,L,N) = MAX( STT(I,J,L,N), 0d0 )
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      ! Free pointer
      NULLIFY( STT )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### GCAP_GLOB_ADV: a TPCORE' ) 
      
      END SUBROUTINE GCAP_GLOBAL_ADV
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      USE PRESSURE_MOD,              ONLY : GET_PEDGE
      USE PRESSURE_MOD,              ONLY : SET_FLOATING_PRESSURE
      USE TIME_MOD,                  ONLY : GET_TS_DYN
      USE TPCORE_BC_MOD,             ONLY : I0_W, J0_W, I1_W, J1_W
      USE TPCORE_BC_MOD,             ONLY : I2_W, J2_W, IM_W, JM_W, IGZD
      USE TPCORE_BC_MOD,             ONLY : DO_WINDOW_TPCORE_BC
      USE TPCORE_GEOS5_WINDOW_MOD,   ONLY : TPCORE_GEOS5_WINDOW
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
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER          :: I0,J0
      INTEGER          :: I, J, L, L2, N, N_DYN, N_TOT_TRC
      REAL*8           :: A_DIFF, D_DYN, TR_DIFF
      REAL*8           :: AD_A(IIPAR,JJPAR,LLPAR)
      REAL*8           :: AD_B(IIPAR,JJPAR,LLPAR)
      REAL*8           :: P_TP1(IIPAR,JJPAR)
      REAL*8           :: P_TP2(IIPAR,JJPAR)
      REAL*8           :: P_TEMP(IIPAR,JJPAR)
      REAL*8           :: TR_A(IIPAR,JJPAR,LLPAR)
#if   defined( APM )
      REAL*8           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                         Input_Opt%N_TRACERS+N_APMTRA) 
#else
      REAL*8           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                         Input_Opt%N_TRACERS)
#endif
      REAL*8,  TARGET  :: XMASS(IIPAR,JJPAR,LLPAR)
      REAL*8,  TARGET  :: YMASS(IIPAR,JJPAR,LLPAR)

      ! Pointers
      REAL*8,  POINTER :: p_UWND (:,:,:  )
      REAL*8,  POINTER :: p_VWND (:,:,:  )
      REAL*8,  POINTER :: p_XMASS(:,:,:  )
      REAL*8,  POINTER :: p_YMASS(:,:,:  )
      REAL*8,  POINTER :: p_STT  (:,:,:,:)
      REAL*8,  POINTER :: STT    (:,:,:,:)

      ! For fields from Input_Opt
      LOGICAL          :: LFILL
      LOGICAL          :: LMFCT
      LOGICAL          :: LPRT
      LOGICAL          :: LWINDO
      INTEGER          :: N_TRACERS
      REAL*8           :: TCVV(Input_Opt%N_TRACERS)

      ! Variable to ensure mass conservation (ccc, 20/11/08)
      REAL*8           :: SUMADA

      !=================================================================
      ! DO_GEOS5_FVDAS_WINDOW_TRANSPORT begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LFILL     = Input_Opt%LFILL
      LMFCT     = Input_Opt%LMFCT
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

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

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
         P_TP1(I,J) = GET_PEDGE(I,J,1)

         ! True surface pressure at end of dynamic timestep [hPa]
         P_TP2(I,J) = State_Met%PSC2(I,J)
  
      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Get the air and tracer mass before advection
      !=================================================================

      ! Airmass [kg] before transport
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         AD_B(I,J,L) = GET_AIR_MASS( I, J, L, P_TP1(I,J) )
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
         TR_B(I,J,L,N) = STT(I,J,L,N) * AD_B(I,J,L) / TCVV(N)
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
      XMASS = 0d0 !(dan)
      YMASS = 0d0

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
      p_UWND  => State_Met%U(:,:,LLPAR:1:-1  )
      p_VWND  => State_Met%V(:,:,LLPAR:1:-1  )
      p_STT   => STT  (:,:,LLPAR:1:-1,:)
      p_XMASS => XMASS(:,:,LLPAR:1:-1  )
      p_YMASS => YMASS(:,:,LLPAR:1:-1  )

      IF ( LPRT ) CALL DEBUG_MSG( '### GEOS5_NESTED: a array flip' )

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
     &                          TCVV,     ND24,      ND25,     ND26 )

      ! Free pointer memory
      NULLIFY( p_UWND, p_VWND, p_STT, p_XMASS, p_YMASS, STT ) 

      !=================================================================
      ! Reset surface pressure and ensure mass conservation
      !=================================================================

      ! Reset the floating surface pressure with P_TP2, the "true"
      ! surface pressure at the end of the dynamic timestep.
      CALL SET_FLOATING_PRESSURE( P_TP2 )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### G5_NESTED_ADV: a TPCORE' )

      END SUBROUTINE DO_GEOS5_WINDOW_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
      SUBROUTINE DO_GEOS57_WINDOW_TRANSPORT( am_I_Root, Input_Opt,
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
      USE PJC_PFIX_GEOS57_WINDOW_MOD,ONLY : DO_PJC_PFIX_GEOS57_WINDOW
      USE PRESSURE_MOD,              ONLY : GET_PEDGE
      USE PRESSURE_MOD,              ONLY : SET_FLOATING_PRESSURE
      USE TIME_MOD,                  ONLY : GET_TS_DYN
      USE TPCORE_BC_MOD,             ONLY : I0_W, J0_W, I1_W, J1_W
      USE TPCORE_BC_MOD,             ONLY : I2_W, J2_W, IM_W, JM_W, IGZD
      USE TPCORE_BC_MOD,             ONLY : DO_WINDOW_TPCORE_BC 
      USE TPCORE_GEOS57_WINDOW_MOD,  ONLY : TPCORE_GEOS57_WINDOW
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
      TYPE(MetState), INTENT(IN)    :: State_Met   ! Meteorology State object
!
! !INPUT/OUTPUT PARAMETERS:
!
      TYPE(ChmState), INTENT(INOUT) :: State_Chm   ! Chemistry State object
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER          :: I0,J0
      INTEGER          :: I, J, L, L2, N, N_DYN, N_TOT_TRC
      REAL*8           :: A_DIFF, D_DYN, TR_DIFF
      REAL*8           :: AD_A(IIPAR,JJPAR,LLPAR)
      REAL*8           :: AD_B(IIPAR,JJPAR,LLPAR)
      REAL*8           :: P_TP1(IIPAR,JJPAR)
      REAL*8           :: P_TP2(IIPAR,JJPAR)
      REAL*8           :: P_TEMP(IIPAR,JJPAR)
      REAL*8           :: TR_A(IIPAR,JJPAR,LLPAR)
#if   defined( APM )
      REAL*8           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                         Input_Opt%N_TRACERS+N_APMTRA)
#else
      REAL*8           :: TR_B(IIPAR,JJPAR,LLPAR,
     &                         Input_Opt%N_TRACERS)
#endif
      REAL*8,  TARGET  :: XMASS(IIPAR,JJPAR,LLPAR)
      REAL*8,  TARGET  :: YMASS(IIPAR,JJPAR,LLPAR)

      ! Pointers
      REAL*8,  POINTER :: p_UWND (:,:,:  )
      REAL*8,  POINTER :: p_VWND (:,:,:  )
      REAL*8,  POINTER :: p_XMASS(:,:,:  )
      REAL*8,  POINTER :: p_YMASS(:,:,:  )
      REAL*8,  POINTER :: p_STT  (:,:,:,:)
      REAL*8,  POINTER :: STT(:,:,:,:)

      ! For fields from Input_Opt
      LOGICAL          :: LFILL
      LOGICAL          :: LMFCT
      LOGICAL          :: LPRT
      LOGICAL          :: LWINDO
      INTEGER          :: N_TRACERS
      REAL*8           :: TCVV(Input_Opt%N_TRACERS)

      !=================================================================
      ! DO_GEOS57_FVDAS_WINDOW_TRANSPORT begins here!
      !=================================================================

      ! Assume success
      RC        =  GIGC_SUCCESS

      ! Copy values from Input_Opt
      LFILL     = Input_Opt%LFILL
      LMFCT     = Input_Opt%LMFCT
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

      ! Initialize GEOS-Chem tracer array [kg] from Chemistry State object
      ! (mpayer, 12/6/12)
      STT => State_Chm%Tracers

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
         P_TP1(I,J) = GET_PEDGE(I,J,1)

         ! True surface pressure at end of dynamic timestep [hPa]
         P_TP2(I,J) = State_Met%PSC2(I,J)

      ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !=================================================================
      ! Get the air and tracer mass before advection
      !=================================================================

      ! Airmass [kg] before transport
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
      DO L = 1, LLPAR
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         AD_B(I,J,L) = GET_AIR_MASS( I, J, L, P_TP1(I,J) )
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
         TR_B(I,J,L,N) = STT(I,J,L,N) * AD_B(I,J,L) / TCVV(N)
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
      XMASS = 0d0 !(dan)
      YMASS = 0d0
      ! NOTE: P_TP1 and P_TP2 are the true surface pressures!  
      CALL DO_PJC_PFIX_GEOS57_WINDOW( D_DYN,
     &                                P_TP1,        P_TP2,
     &                                State_Met%U,  State_Met%V,
     &                                XMASS,        YMASS )

      IF ( LPRT ) CALL DEBUG_MSG( '### GEOS5_NESTED: a GET_TRA_MASS' )

      ! Impose TPCORE boundary conditions @ edges of nested grid
      CALL DO_WINDOW_TPCORE_BC( am_I_Root, Input_Opt, State_Chm, RC )

      ! Flip array indices in the vertical using pointer storage
      p_UWND  => State_Met%U(:,:,LLPAR:1:-1  )
      p_VWND  => State_Met%V(:,:,LLPAR:1:-1  )
      p_STT   => STT  (:,:,LLPAR:1:-1,:)
      p_XMASS => XMASS(:,:,LLPAR:1:-1  )
      p_YMASS => YMASS(:,:,LLPAR:1:-1  )

      ! Note: the mass flux diagnostic arrays (MASSFLEW, MASSFLNS and MASSFLUP)
      ! are incremented upside-down (level 1 = top of the atmosphere).
      ! The levels order is reversed only when written out in diag3.f
      ! (ccc, 3/8/10)
      ! Do the advection
      CALL TPCORE_GEOS57_WINDOW( D_DYN,    Re,        IIPAR,    JJPAR,
     &                           LLPAR,    JFIRST,    JLAST,    NG,
     &                           MG,       N_TOT_TRC, Ap,       Bp,
     &                           p_UWND,   p_VWND,    P_TP1,    P_TP2,
     &                           P_TEMP,   p_STT,     IORD,     JORD,      
     &                           KORD,     N_ADJ,     p_XMASS,  p_YMASS,
     &                           MASSFLEW, MASSFLNS,  MASSFLUP, A_M2,
     &                           TCVV,     ND24,      ND25,     ND26 )

      ! Free pointer memory
      NULLIFY( p_UWND, p_VWND, p_STT, p_XMASS, p_YMASS, STT )

      !=================================================================
      ! Reset surface pressure and ensure mass conservation
      !=================================================================
      ! Reset the floating surface pressure with P_TP2, the "true"
      ! surface pressure at the end of the dynamic timestep.
      CALL SET_FLOATING_PRESSURE( P_TP2 )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### G5_NESTED_ADV: a TPCORE' )

      END SUBROUTINE DO_GEOS57_WINDOW_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_air_mass
!
! !DESCRIPTION: Function GET\_AIR\_MASS returns the air mass based on the 
!  pressures returned before and after the call to the GEOS-4/fvDAS TPCORE 
!  code. (bmy, 6/24/03)
!\\
!\\
! !INTERFACE:
!
      FUNCTION GET_AIR_MASS( I, J, L, P_SURF ) RESULT( AIR_MASS )
!
! !USES:
!
      USE CMN_SIZE_MOD               ! Size parameters
      USE CMN_GCTM_MOD               ! g0_100
!
! !INPUT PARAMETERS:
!
      INTEGER, INTENT(IN) :: I, J, L   ! GEOS-Chem lon, lat, level indices
      REAL*8,  INTENT(IN) :: P_SURF    ! Surface pressure [hPa] at (I,J,L=1)

! 
! !REVISION HISTORY: 
!  24 Jun 2003 - R. Yantosca - Initial version
!  26 Feb 2010 - R. Yantosca - Added ProTeX headers
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: L2
      REAL*8  :: P_BOT, P_TOP, AIR_MASS
      
      !=================================================================
      ! GET_AIR_MASS begins here!
      !=================================================================
      
      ! Index for Ap, Bp from atmosphere top down to surface
      ! since the Ap's and Bp's have been flipped for TPCORE
      L2       = ( LLPAR + 1 ) - L + 1
               
      ! Hybrid-grid formulation for air mass
      P_BOT    = Ap(L2)   + ( Bp(L2)   * P_SURF )
      P_TOP    = Ap(L2-1) + ( Bp(L2-1) * P_SURF )
      AIR_MASS = ( P_BOT - P_TOP ) * G0_100 * A_M2(J)

      END FUNCTION GET_AIR_MASS
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: J, K, L, N_DYN
      REAL*8  :: YMID_R(JJPAR)

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

      ! Call INIT routine from "tpcore_fvdas_mod.f" 
      CALL INIT_TPCORE( IIPAR,  JJPAR, LLPAR,  JFIRST, JLAST, 
     &                  NG, MG, DBLE( N_DYN ), Re,     YMID_R )


      END SUBROUTINE INIT_TRANSPORT

!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: J, K, L, N_DYN
      REAL*8  :: YMID_R_W(0:JJPAR+1)

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

      ! YMID_R is latitude of grid box center [radians]
      DO J =0, JJPAR+1
         YMID_R_W(J) = GET_YMID_R_W( 1, J, 1 )
      ENDDO

      ! Call INIT routine from "tpcore_geos5_fvdas_window_mod.f"
      CALL INIT_GEOS5_WINDOW( IIPAR, JJPAR,    LLPAR, JFIRST, 
     &                        JLAST, NG,       MG,    DBLE( N_DYN ), 
     &                        Re,    YMID_R_W )

      END SUBROUTINE INIT_GEOS5_WINDOW_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_geos57_window_transport

!
! !DESCRIPTION: Subroutine INIT\_GEOS57\_WINDOW\_TRANSPORT initializes all 
!  module variables and arrays for the GEOS-57 nested grid simulation.  
!  This routine is only called if we are using the GEOS-57 nested grid 
!  simulation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE INIT_GEOS57_WINDOW_TRANSPORT( am_I_Root,Input_Opt,RC )
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
      USE TPCORE_GEOS57_WINDOW_MOD,ONLY : INIT_GEOS57_WINDOW

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
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      INTEGER :: J, K, L, N_DYN
      REAL*8  :: YMID_R_W(0:JJPAR+1)

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

      ! YMID_R is latitude of grid box center [radians]
      DO J =0, JJPAR+1
         YMID_R_W(J) = GET_YMID_R_W( 1, J, 1 )
      ENDDO

      ! Call INIT routine from "tpcore_geos57_fvdas_window_mod.f"
      ! subroutine "init_geos57_windows" is in "tpcore_geos57_window_mod.F90"
      CALL INIT_GEOS57_WINDOW( IIPAR, JJPAR,    LLPAR, JFIRST, 
     &                         JLAST, NG,       MG,    DBLE( N_DYN ),
     &                         Re,    YMID_R_W )

      END SUBROUTINE INIT_GEOS57_WINDOW_TRANSPORT
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
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
