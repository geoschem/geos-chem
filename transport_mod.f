! $Id: transport_mod.f,v 1.23 2009/05/29 14:31:06 ccarouge Exp $
      MODULE TRANSPORT_MOD
!
!******************************************************************************
!  Module TRANSPORT_MOD is used to call the proper version of TPCORE for
!  GEOS-3, GEOS-4, GEOS-5, or GCAP nested-grid or global simulations.
!  (yxw, bmy, 3/10/03, 2/17/09)
! 
!  Module Variables:
!  ============================================================================
!  (1 ) Ap     (REAL*8 )       : Vertical coordinate array for TPCORE
!  (2 ) A_M2   (REAL*8 )       : Grid box surface areas [m2]
!  (3 ) Bp     (REAL*8 )       : Vertical coordinate array for TPCORE
!  (4 ) IORD   (REAL*8 )       : TPCORE E/W option flag
!  (5 ) JORD   (REAL*8 )       : TPCORE N/S option flag
!  (6 ) KORD   (REAL*8 )       : TPCORE vertical option flag
!  (7 ) JLAST  (INTEGER)       : For fvDAS TPCORE
!  (8 ) MG     (INTEGER)       : For fvDAS TPCORE
!  (9 ) NG     (INTEGER)       : For fvDAS TPCORE
!  (10) N_ADJ  (INTEGER)       : For fvDAS TPCORE
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_TRANSPORT                : Driver: calls global or window TPCORE
!  (2 ) GEOS4_GEOS5_GLOBAL_ADV      : TPCORE driver for GEOS-4/GEOS-5 met
!  (3 ) GEOS3_GLOBAL_ADV            : TPCORE driver routine for GEOS-3 met
!  (4 ) GCAP_GLOBAL_ADV             : TPCORE driver routine for GCAP met
!  (5 ) DO_WINDOW_TRANSPORT         : Calls nested-grid version of TPCORE
!  (6 ) GET_AIR_MASS                : Computes air mass from TPCORE in/out P's
!  (7 ) SET_TRANSPORT               : Gets IORD, JORD, KORD from "input_mod.f"
!  (8 ) INIT_TRANSPORT              : Initializes module arrays
!  (9 ) INIT_GEOS5_WINDOW_TRANSPORT : Initializes module arrays
!  (10) CLEANUP_TRANSPORT           : Deallocates module arrays
!
!  GEOS-Chem modules referenced by transport_mod.f
!  ============================================================================
!  (1 ) dao_mod.f              : Module w/ arrays for DAO met fields
!  (2 ) diag_mod.f             : Module w/ GEOS-CHEM diagnostic arrays
!  (3 ) error_mod.f            : Module w/ I/O error/NaN check routines
!  (4 ) grid_mod.f             : Module w/ horizontal grid information
!  (5 ) logical_mod.f          : Module w/ GEOS-CHEM logical switches
!  (6 ) pjc_pfix_mod.f         : Module w/ Phil Cameron-Smith P-fixer
!  (7 ) pressure_mod.f         : Module w/ routines to compute P(I,J,L)
!  (8 ) time_mod.f             : Module w/ routines to compute date/time
!  (9 ) tpcore_mod.f           : Module w/ TPCORE for GEOS1,GEOSS,GEOS3
!  (10) tpcore_bc_mod.f        : Module w/ TPCORE boundary cond. routines
!  (11) tpcore_window_mod.f    : Module w/ TPCORE for nested-grid windows
!  (12) tpcore_fvdas_mod.f90   : Module w/ TPCORE for GEOS-4/fvDAS
!  (12) tpcore_geos5_window_mod.f90: Module for GEOS-5 nested grid
!  (13) tracer_mod.f           : Module w/ GEOS-CHEM tracer array STT etc
!
!  NOTES:
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
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "transport_mod.f"
      !=================================================================

      ! Make everything PRIVATE...
      PRIVATE

      ! ... except these routines
      PUBLIC :: CLEANUP_TRANSPORT
      PUBLIC :: DO_TRANSPORT
      PUBLIC :: INIT_TRANSPORT
      PUBLIC :: INIT_GEOS5_WINDOW_TRANSPORT
      PUBLIC :: SET_TRANSPORT

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      INTEGER             :: IORD,  JORD, KORD, JFIRST 
      INTEGER             :: JLAST, NG,   MG,   N_ADJ
      REAL*8, ALLOCATABLE :: Ap(:)
      REAL*8, ALLOCATABLE :: A_M2(:)
      REAL*8, ALLOCATABLE :: Bp(:)
      REAL*8, ALLOCATABLE :: STT_BC2(:,:,:)

      !=================================================================
      ! MODULE ROUTINES -- follow below the "CONTAINS" statement
      !=================================================================
      CONTAINS

!------------------------------------------------------------------------------
      
      SUBROUTINE DO_TRANSPORT
!
!******************************************************************************
!  Subroutine DO_TRANSPORT is the driver routine for the proper TPCORE program
!  for GEOS-3, GEOS-4/GEOS-5, or window simulations. (bmy, 3/10/03, 10/30/07)
! 
!  NOTES:
!  (1 ) Removed IORD, JORD, KORD from the arg list.  Also now removed
!        reference to CMN, it's not needed. (bmy, 7/20/04)
!  (2 ) Now call separate routines for different met fields. (bmy, 10/30/07)
!  (3 ) Now references subroutine INIT_TPCORE_BC from tpcore_bc_mod.f and
!        DO_GEOS5_FVDAS_WINDOW_TRANSPORT from 
!        "tpcore_geos5_fvdas_window_mod.f90". (yxw, dan, bmy, 11/6/08)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : ITS_A_NESTED_GRID
      USE TPCORE_BC_MOD, ONLY : INIT_TPCORE_BC

#     include "CMN_SIZE" ! Size parameters

      ! Local variables
      LOGICAL, SAVE     :: FIRST = .TRUE.

      !=================================================================
      ! DO_TRANSPORT begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN 
#if   defined( GRID05x0666 )
      CALL INIT_GEOS5_WINDOW_TRANSPORT
#else
      CALL INIT_TRANSPORT
#endif
         FIRST = .FALSE.
      ENDIF

      ! Choose the proper version of TPCORE for the nested-grid window 
      ! region (usually 1x1 grids) or for the entire globe
      IF ( ITS_A_NESTED_GRID() ) THEN
#if   defined( GRID1x1 )
         CALL DO_WINDOW_TRANSPORT
#elif defined( GRID05x0666 )
         CALL DO_GEOS5_WINDOW_TRANSPORT
#endif
      ELSE

#if defined( GEOS_4 ) || defined( GEOS_5 )

         ! Call TPCORE w/ proper settings for GEOS-4/GEOS-5 met
         CALL GEOS4_GEOS5_GLOBAL_ADV

#elif defined( GEOS_3 )

         ! Call TPCORE w/ proper settings for GEOS-3 met
         CALL GEOS3_GLOBAL_ADV

#elif defined( GCAP )

         ! Call TPCORE w/ proper settings for GCAP met
         CALL GCAP_GLOBAL_ADV

#endif
      ENDIF

      ! Return to calling program
      END SUBROUTINE DO_TRANSPORT

!------------------------------------------------------------------------------
      
      SUBROUTINE GEOS4_GEOS5_GLOBAL_ADV
!     
!******************************************************************************
!  Subroutine GEOS4_GEOS_5_GLOBAL_ADV is the driver routine for TPCORE with 
!  the GMAO GEOS-4 or GEOS-5 met fields. (bmy, 10/30/07, 2/17/09)
!     
!  NOTES:
!  (1 ) Split off the GEOS-4 & GEOS-5 relevant parts from the previous 
!        routine DO_GLOBAL_TRANSPORT (bmy, 10/30/07)
!  (2 ) Activate the call to SAVE_GLOBAL_TPCORE_BC (yxw, dan, bmy, 11/6/08)
!  (3 ) Bug fix in mass balance: only account for cells of STT with non-zero
!        concentrations when doing the computation (ccc, bmy, 2/17/09)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,          ONLY : PSC2, UWND, VWND
      USE DIAG_MOD,         ONLY : MASSFLEW, MASSFLNS, MASSFLUP
      USE ERROR_MOD,        ONLY : IT_IS_NAN, DEBUG_MSG, SAFE_DIV
      USE LOGICAL_MOD,      ONLY : LEMBED, LFILL, LMFCT, LPRT, LWINDO
      USE PJC_PFIX_MOD,     ONLY : DO_PJC_PFIX
      USE PRESSURE_MOD,     ONLY : GET_PEDGE, SET_FLOATING_PRESSURE
      USE TIME_MOD,         ONLY : GET_TS_DYN
      USE TPCORE_BC_MOD,    ONLY : SAVE_GLOBAL_TPCORE_BC
      USE TPCORE_FVDAS_MOD, ONLY : TPCORE_FVDAS
      USE TRACER_MOD,       ONLY : N_TRACERS, STT, TCVV

#     include "CMN_SIZE"         ! Size parameters        
#     include "CMN_GCTM"         ! Physical constants
#     include "CMN_DIAG"         ! NDxx flags
      
      ! Local variables
      INTEGER                   :: I, J, L, L2, N, N_DYN
      REAL*8                    :: A_DIFF, D_DYN, TR_DIFF
      REAL*8                    :: AD_A(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: AD_B(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: P_TP1(IIPAR,JJPAR)
      REAL*8                    :: P_TP2(IIPAR,JJPAR)
      REAL*8                    :: P_TEMP(IIPAR,JJPAR)
      REAL*8                    :: TR_A(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: TR_B(IIPAR,JJPAR,LLPAR,N_TRACERS)
      REAL*8                    :: UTMP(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: VTMP(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: XMASS(IIPAR,JJPAR,LLPAR) 
      REAL*8                    :: YMASS(IIPAR,JJPAR,LLPAR) 

      ! Variable to ensure mass conservation (ccc, 2/17/09)
      REAL*8                    :: SUMADA

      !=================================================================
      ! GEOS4_GEOS5_GLOBAL_ADV begins here!
      !=================================================================

      ! Save boundary conditions (global grid) for future nested run
      IF ( LWINDO ) CALL SAVE_GLOBAL_TPCORE_BC
       
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
         P_TP2(I,J) = PSC2(I,J)    
   
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
      DO N = 1, N_TRACERS
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
      CALL DO_PJC_PFIX( D_DYN, P_TP1, P_TP2, 
     &                  UWND,  VWND,  XMASS, YMASS )

      !=================================================================
      ! Call TPCORE_FVDAS to perform the advection
      !=================================================================

      ! Store winds in UTMP, VTMP to preserve UWND, VWND 
      ! for diagnostics.  Flip in the vertial for TPCORE.
      UTMP(:,:,1:LLPAR) = UWND(:,:,LLPAR:1:-1)
      VTMP(:,:,1:LLPAR) = VWND(:,:,LLPAR:1:-1)

      ! Do the advection
      CALL TPCORE_FVDAS( D_DYN,    Re,        IIPAR,    JJPAR,
     &                   LLPAR,    JFIRST,    JLAST,    NG,
     &                   MG,       N_TRACERS, Ap,       Bp,
     &                   UTMP,     VTMP,      P_TP1,    P_TP2,
     &                   P_TEMP,   STT(:,:,LLPAR:1:-1,:),       
     &                   IORD,     JORD,      KORD,     N_ADJ,     
     &                   XMASS(:,:,LLPAR:1:-1),    
     &                   YMASS(:,:,LLPAR:1:-1),
     &                   MASSFLEW(:,:,LLPAR:1:-1,:), 
     &                   MASSFLNS(:,:,LLPAR:1:-1,:),  
     &                   MASSFLUP(:,:,LLPAR:1:-1,:),    A_M2,
     &                   TCVV,     ND24,      ND25,     ND26 )

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
      DO N = 1, N_TRACERS

         ! Zero summing variable
         SUMADA = 0.d0

         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Air mass [kg] after transport
! prior to 4/15/09 (ccc)
!            IF ( N == 1 ) THEN
               AD_A(I,J,L) = GET_AIR_MASS( I, J, L, P_TP2(I,J) )
!            ENDIF
         
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
         !-------------------------------------------------------------------
         ! Prior to 2/17/09:
         ! This was changed to take into account only cells 
         ! w/ nonzero tracer. (ccc, 2/17/09)
         ! TR_DIFF = TR_DIFF / SUM( AD_A ) * TCVV(N)
         !-------------------------------------------------------------------
         TR_DIFF = SAFE_DIV(TR_DIFF, SUMADA, 0.d0) * TCVV(N)

         ! Add mass difference [v/v] back to STT
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            !------------------------------------------------------------------
            ! Prior to 2/17/09:
            ! Only apply change to cells w/ nonzero tracer (ccc, 2/17/09)
            ! STT(I,J,L,N) = STT(I,J,L,N) + TR_DIFF
            !------------------------------------------------------------------
            IF ( STT(I,J,L,N) > 0.d0 .or. STT(I,J,L,N) < 0.d0 ) THEN
               STT(I,J,L,N) = STT(I,J,L,N) + TR_DIFF
            ENDIF

            STT(I,J,L,N) = MAX( STT(I,J,L,N), 0d0 )
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### G4_G5_GLOB_ADV: a TPCORE' ) 

      ! Return to calling program
      END SUBROUTINE GEOS4_GEOS5_GLOBAL_ADV

!------------------------------------------------------------------------------

      SUBROUTINE GEOS3_GLOBAL_ADV
!
!******************************************************************************
!  Subroutine GEOS3_GLOBAL_ADV is the driver routine for TPCORE with the 
!  GMAO GEOS-3 met fields. (bmy, 10/30/07)
!     
!  NOTES:
!  (1 ) Split off the GEOS-3 relevant parts from the previous routine
!        DO_GLOBAL_TRANSPORT (bmy, 10/30/07)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,          ONLY : PSC2, UWND, VWND
      USE DIAG_MOD,         ONLY : MASSFLEW, MASSFLNS, MASSFLUP
      USE ERROR_MOD,        ONLY : IT_IS_NAN, DEBUG_MSG
      USE LOGICAL_MOD,      ONLY : LEMBED, LFILL, LMFCT, LPRT, LWINDO
      USE PRESSURE_MOD,     ONLY : GET_PEDGE, SET_FLOATING_PRESSURE
      USE TIME_MOD,         ONLY : GET_TS_DYN
      USE TPCORE_BC_MOD,    ONLY : SAVE_GLOBAL_TPCORE_BC
      USE TPCORE_MOD,       ONLY : TPCORE
      USE TRACER_MOD,       ONLY : N_TRACERS, STT, TCVV

#     include "CMN_SIZE"         ! Size parameters
#     include "CMN_GCTM"         ! Physical constants
#     include "CMN_DIAG"         ! NDxx flags

      ! Local variables
      INTEGER                   :: I, J, L, N_DYN
      REAL*8                    :: A_DIFF, D_DYN, TR_DIFF
      REAL*8                    :: P_TP1(IIPAR,JJPAR)
      REAL*8                    :: P_TP2(IIPAR,JJPAR)
      REAL*8                    :: UTMP(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: VTMP(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: WW(IIPAR,JJPAR,LLPAR) 

      ! Parameters
      INTEGER, PARAMETER        :: IGD=0, J1=3
      REAL*8,  PARAMETER        :: Umax=200d0

      !=================================================================
      ! GEOS3_GLOBAL_ADV begins here!
      !=================================================================

      ! Save boundary conditions (global grid) for future nested run
      IF ( LWINDO ) CALL SAVE_GLOBAL_TPCORE_BC

      ! Dynamic timestep [s]
      N_DYN = GET_TS_DYN() * 60
      D_DYN = DBLE( N_DYN )

      !=================================================================
      ! Prepare variables for calls to TPCORE
      !
      ! For GEOS-3 (pure sigma grid), the pressure at the bottom edge 
      ! of grid box (I,J,L) is given by:
      !
      !    P(I,J,L) = Ap(L) + [ Bp(L) * ( Psurface(I,J) - PTOP ) ]
      !
      ! where Psurface is the true surface pressure (i.e. not PS-PTOP).
      ! and Ap(L), Bp(L) define the vertical grid (see pressure_mod.f)
      !
      ! Therefore, we construct the 2-D surface pressure arrays that
      ! are required for TPCORE:
      !
      !    P_TP1(I,J) = ( Psurf(I,J) - PTOP ) at midpt of dyn timestep
      !    P_TP2(I,J) = ( Psurf(I,J) - PTOP ) at end   of dyn timestep
      !
      ! Note that P_TP1 and P_TP2 need to be ( Psurface - PTOP ) in 
      ! order to be consistent with the definition of Ap and Bp that
      ! defines the GEOS-3 pure-sigma grid.
      !
      ! Also note that TPCORE will call the pressure fixer internally, 
      ! which will ensure mass conservation.  However, we must reset 
      ! the floating pressure with the output pressure after advection.
      !=================================================================

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J )
      DO J = 1, JJPAR
      DO I = 1, IIPAR

         ! Psurface - PTOP at midpt of dynamic timestep [hPa]
         P_TP1(I,J) = GET_PEDGE(I,J,1) - PTOP

         ! Psurface - PTOP at end   of dynamic timestep [hPa]
         P_TP2(I,J) = PSC2(I,J)        - PTOP

      ENDDO
      ENDDO
!$OMP END PARALLEL DO
      
      !==============================================================
      ! Call TPCORE to perform the advection
      !==============================================================

      ! Store winds in UTMP, VTMP to preserve UWND, VWND 
      ! for diagnostics.  Flip in the vertical for TPCORE.
      UTMP(:,:,1:LLPAR) = UWND(:,:,LLPAR:1:-1)
      VTMP(:,:,1:LLPAR) = VWND(:,:,LLPAR:1:-1)

      ! TPCORE v7.1.m transport scheme.  Pass P_TP1 and P_TP2, as are
      ! computed above.  P_TP2 will be overwritten with the new value
      ! of (Psurface-PTOP) after transport.  NOTE: TPCORE assumes that 
      ! L=1 is atm top, so flip in the vertical. (bmy, 10/30/07)
      CALL TPCORE( IGD,   STT(:,:,LLPAR:1:-1,:),
     &             P_TP1, P_TP2, UTMP, VTMP,  WW,    
     &             N_DYN, IORD,  JORD, KORD,  N_TRACERS, 
     &             IIPAR, JJPAR, J1,   LLPAR, Ap,   
     &             Bp,    PTOP,  Re,   LFILL, LMFCT, Umax )

      !=================================================================
      ! Reset surface pressure in order to ensure mass conservation
      !=================================================================

      ! Reset floating pressure w/ pressure adjusted by TPCORE
      ! P_TP2 is PS-PTOP, so reset the pressure with P_TP2+PTOP.
      CALL SET_FLOATING_PRESSURE( P_TP2 + PTOP )
           
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### GEOS3_GLOB_ADV: a TPCORE' ) 

      ! Return to calling program
      END SUBROUTINE GEOS3_GLOBAL_ADV

!------------------------------------------------------------------------------

      SUBROUTINE GCAP_GLOBAL_ADV
!
!******************************************************************************
!  Subroutine GCAP_GLOBAL_ADV is the driver routine for TPCORE with the 
!  GCAP / GISS met fields. (bmy, 10/30/07, 2/17/09)
!     
!  NOTES:
!  (1 ) Split off the GCAP relevant parts from the previous routine
!        DO_GLOBAL_TRANSPORT (bmy, 10/30/07)
!  (2 ) Bug fix in mass balance: only account for cells of STT with non-zero
!        concentrations when doing the computation (ccc, bmy, 2/17/09)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,          ONLY : PSC2, UWND, VWND
      USE DIAG_MOD,         ONLY : MASSFLEW, MASSFLNS, MASSFLUP
      USE ERROR_MOD,        ONLY : IT_IS_NAN, DEBUG_MSG
      USE LOGICAL_MOD,      ONLY : LEMBED, LFILL, LMFCT, LPRT, LWINDO
      USE PJC_PFIX_MOD,     ONLY : DO_PJC_PFIX
      USE PRESSURE_MOD,     ONLY : GET_PEDGE, SET_FLOATING_PRESSURE
      USE TIME_MOD,         ONLY : GET_TS_DYN
      USE TPCORE_FVDAS_MOD, ONLY : TPCORE_FVDAS
      USE TRACER_MOD,       ONLY : N_TRACERS, STT, TCVV

#     include "CMN_SIZE"         ! Size parameters
#     include "CMN_GCTM"         ! Physical constants
#     include "CMN_DIAG"         ! NDxx flags

      ! Local variables
      INTEGER                   :: I, J, L, L2, N, N_DYN
      REAL*8                    :: A_DIFF, D_DYN, TR_DIFF
      REAL*8                    :: AD_A(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: AD_B(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: P_TP1(IIPAR,JJPAR)
      REAL*8                    :: P_TP2(IIPAR,JJPAR)
      REAL*8                    :: P_TEMP(IIPAR,JJPAR)
      REAL*8                    :: TR_A(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: TR_B(IIPAR,JJPAR,LLPAR,N_TRACERS)
      REAL*8                    :: UTMP(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: VTMP(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: XMASS(IIPAR,JJPAR,LLPAR) 
      REAL*8                    :: YMASS(IIPAR,JJPAR,LLPAR) 

      ! Variable to ensure mass conservation (ccc, 20/11/08)
      REAL*8                    :: SUMADA

      !=================================================================
      ! GCAP_GLOBAL_ADV begins here!
      !=================================================================

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
         P_TP2(I,J) = PSC2(I,J)        - PTOP

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
      DO N = 1, N_TRACERS
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
      CALL DO_PJC_PFIX( D_DYN, P_TP1+PTOP, P_TP2+PTOP, 
     &                  UWND,  VWND,       XMASS,     YMASS )

      !=================================================================
      ! Call TPCORE_FVDAS to perform the advection
      !=================================================================

      ! Store winds in UTMP, VTMP to preserve UWND, VWND 
      ! for diagnostics.  Flip in the vertical for TPCORE.
      UTMP(:,:,1:LLPAR) = UWND(:,:,LLPAR:1:-1)
      VTMP(:,:,1:LLPAR) = VWND(:,:,LLPAR:1:-1)

      ! Do the advection
      CALL TPCORE_FVDAS( D_DYN,    Re,        IIPAR,    JJPAR,
     &                   LLPAR,    JFIRST,    JLAST,    NG,
     &                   MG,       N_TRACERS, Ap,       Bp,
     &                   UTMP,     VTMP,      P_TP1,    P_TP2,
     &                   P_TEMP,   STT(:,:,LLPAR:1:-1,:),       
     &                   IORD,     JORD,      KORD,     N_ADJ,     
     &                   XMASS(:,:,LLPAR:1:-1),    
     &                   YMASS(:,:,LLPAR:1:-1),
     &                   MASSFLEW(:,:,LLPAR:1:-1,:), 
     &                   MASSFLNS(:,:,LLPAR:1:-1,:),  
     &                   MASSFLUP(:,:,LLPAR:1:-1,:),    A_M2,
     &                   TCVV,     ND24,      ND25,     ND26 )

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
      DO N = 1, N_TRACERS

         ! Zero summing variable
         SUMADA = 0.d0

         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Air mass [kg] after transport
! prior to 4/15/09 (ccc)
!            IF ( N == 1 ) THEN
               AD_A(I,J,L) = GET_AIR_MASS( I, J, L, P_TP2(I,J) )
!            ENDIF
         
            ! Tracer mass [kg] after transport
            TR_A(I,J,L) = STT(I,J,L,N) * AD_A(I,J,L) / TCVV(N)

            ! We apply mass conservation only on cells with 
            ! nonzero tracer. (ccc, 11/20/08)
            IF ( STT(I,J,L,N) > 0.d0 ) THEN
               SUMADA = SUMADA + AD_A(I,J,L)
            ENDIF
         ENDDO
         ENDDO
         ENDDO

         ! Residual mass difference [kg]: before - after
         TR_DIFF = SUM( TR_B(:,:,:,N) ) - SUM( TR_A )

         ! Convert from [kg] to [v/v]
         !-------------------------------------------------------------------
         ! Prior to 2/17/09:
         ! This was changed to take into account only cells w/
         ! nonzero tracer. (ccc, 2/17/09)
         ! TR_DIFF = TR_DIFF / SUM( AD_A ) * TCVV(N)
         !-------------------------------------------------------------------
         TR_DIFF = TR_DIFF / SUMADA * TCVV(N)

         ! Add mass difference [v/v] back to STT
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            !----------------------------------------------------------------
            ! Prior to 2/17/09:
            ! This was changed to take into account only cells 
            ! w/ nonzero tracer. (ccc, 2/17/09)
            ! STT(I,J,L,N) = STT(I,J,L,N) + TR_DIFF
            !----------------------------------------------------------------
            IF ( STT(I,J,L,N) > 0.d0 ) THEN
               STT(I,J,L,N) = STT(I,J,L,N) + TR_DIFF
            ENDIF

            STT(I,J,L,N) = MAX( STT(I,J,L,N), 0d0 )
         ENDDO
         ENDDO
         ENDDO
      ENDDO
!$OMP END PARALLEL DO

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### GCAP_GLOB_ADV: a TPCORE' ) 
      
      ! Return to calling program
      END SUBROUTINE GCAP_GLOBAL_ADV

!------------------------------------------------------------------------------

      SUBROUTINE DO_GEOS5_WINDOW_TRANSPORT
!
!******************************************************************************
!  Subroutine GEOS4_GEOS_5_GLOBAL_ADV is the driver routine for TPCORE with
!  the GMAO GEOS-4 or GEOS-5 met fields. (bmy, 10/30/07)
!
!  NOTES:
!  (1 ) 
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,                   ONLY : PSC2,     UWND,     VWND
      USE DIAG_MOD,                  ONLY : MASSFLEW, MASSFLNS, MASSFLUP
      USE ERROR_MOD,                 ONLY : IT_IS_NAN,     DEBUG_MSG
      USE GRID_MOD,                  ONLY : GET_XOFFSET,   GET_YOFFSET 
      USE LOGICAL_MOD,               ONLY : LEMBED, LFILL, LMFCT
      USE LOGICAL_MOD,               ONLY : LPRT,   LWINDO
      USE PJC_PFIX_GEOS5_WINDOW_MOD, ONLY : DO_PJC_PFIX_GEOS5_WINDOW
      USE PRESSURE_MOD,              ONLY : GET_PEDGE
      USE PRESSURE_MOD,              ONLY : SET_FLOATING_PRESSURE
      USE TIME_MOD,                  ONLY : GET_TS_DYN
      USE TPCORE_BC_MOD,             ONLY : I0_W, J0_W, I1_W, J1_W
      USE TPCORE_BC_MOD,             ONLY : I2_W, J2_W, IM_W, JM_W, IGZD
      USE TPCORE_BC_MOD,             ONLY : DO_WINDOW_TPCORE_BC
      USE TPCORE_WINDOW_MOD,         ONLY : TPCORE_WINDOW
      USE TPCORE_GEOS5_WINDOW_MOD,   ONLY : TPCORE_GEOS5_WINDOW
      USE TRACER_MOD,                ONLY : N_TRACERS, STT, TCVV

#     include "CMN_SIZE"         ! Size parameters
#     include "CMN_GCTM"         ! Physical constants
#     include "CMN_DIAG"         ! NDxx flags

      ! Local variables
      INTEGER                   :: I0,J0
      INTEGER                   :: I, J, L, L2, N, N_DYN
      REAL*8                    :: A_DIFF, D_DYN, TR_DIFF
      REAL*8                    :: AD_A(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: AD_B(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: P_TP1(IIPAR,JJPAR)
      REAL*8                    :: P_TP2(IIPAR,JJPAR)
      REAL*8                    :: P_TEMP(IIPAR,JJPAR)
      REAL*8                    :: TR_A(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: TR_B(IIPAR,JJPAR,LLPAR,N_TRACERS)
      REAL*8                    :: UTMP(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: VTMP(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: XMASS(IIPAR,JJPAR,LLPAR)
      REAL*8                    :: YMASS(IIPAR,JJPAR,LLPAR)

      !=================================================================
      ! DO_GEOS5_FVDAS_WINDOW_TRANSPORT begins here!
      !=================================================================

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
         P_TP1(I,J) = GET_PEDGE(I,J,1)

         ! True surface pressure at end of dynamic timestep [hPa]
         P_TP2(I,J) = PSC2(I,J)
  
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
      DO N = 1, N_TRACERS
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
      CALL DO_PJC_PFIX_GEOS5_WINDOW( D_DYN, P_TP1, P_TP2,
     &                               UWND,  VWND,  XMASS, YMASS )

      IF ( LPRT ) CALL DEBUG_MSG( '### GEOS5_NESTED: a GET_TRA_MASS' )

      ! Impose TPCORE boundary conditions @ edges of 1x1 grid
      CALL DO_WINDOW_TPCORE_BC

      ! for diagnostics.  Flip in the vertial for TPCORE.
      UTMP(:,:,1:LLPAR) = UWND(:,:,LLPAR:1:-1)
      VTMP(:,:,1:LLPAR) = VWND(:,:,LLPAR:1:-1)

      ! Do the advection
      CALL TPCORE_GEOS5_WINDOW( D_DYN, Re, IIPAR, JJPAR,
     &                          LLPAR,    JFIRST,    JLAST,    NG,
     &                          MG,       N_TRACERS, Ap,       Bp,
     &                          UTMP,     VTMP,      P_TP1,    P_TP2,
     &                          P_TEMP,   STT(:,:,LLPAR:1:-1,:),
     &                          IORD,     JORD,      KORD,     N_ADJ,
     &                          XMASS(:,:,LLPAR:1:-1),
     &                          YMASS(:,:,LLPAR:1:-1),
     &                          MASSFLEW(:,:,LLPAR:1:-1,:),
     &                          MASSFLNS(:,:,LLPAR:1:-1,:),
     &                          MASSFLUP(:,:,LLPAR:1:-1,:),    A_M2,
     &                          TCVV,     ND24,      ND25,     ND26 )


      !=================================================================
      ! Reset surface pressure and ensure mass conservation
      !=================================================================

      ! Reset the floating surface pressure with P_TP2, the "true"
      ! surface pressure at the end of the dynamic timestep.
      CALL SET_FLOATING_PRESSURE( P_TP2 )

      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### G5_NESTED_ADV: a TPCORE' )

      ! Return to calling program
      END SUBROUTINE DO_GEOS5_WINDOW_TRANSPORT

!------------------------------------------------------------------------------

      SUBROUTINE DO_WINDOW_TRANSPORT
!
!******************************************************************************
!  Subroutine DO_WINDOW_TRANSPORT is the driver program for the proper TPCORE
!  program for nested-grid window simulations. (yxw, bmy, 8/7/03, 10/3/05)
!
!  NOTES:
!  (1 ) Now references DEBUG_MSG from "error_mod.f" (bmy, 8/7/03)
!  (2 ) Removed IORD, JORD, KORD from the arg list, since these are now
!        module variables.  Now reference LFILL, LMFCT, LPRT from 
!        "logical_mod.f".  Now reference STT, N_TRACERS from "tracer_mod.f".
!  (3 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,           ONLY : PSC2, UWND, VWND
      USE ERROR_MOD,         ONLY : DEBUG_MSG
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE LOGICAL_MOD,       ONLY : LFILL, LMFCT, LPRT
      USE PRESSURE_MOD,      ONLY : GET_PEDGE, SET_FLOATING_PRESSURE
      USE TIME_MOD,          ONLY : GET_TS_DYN
      USE TPCORE_BC_MOD,     ONLY : I0_W, J0_W, I1_W, J1_W
      USE TPCORE_BC_MOD,     ONLY : I2_W, J2_W, IM_W, JM_W, IGZD 
      USE TPCORE_BC_MOD,     ONLY : DO_WINDOW_TPCORE_BC
      USE TPCORE_WINDOW_MOD, ONLY : TPCORE_WINDOW
      USE TRACER_MOD,        ONLY : STT, N_TRACERS

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"  ! Re

      ! Local variables
      INTEGER             :: I, I0, J, J0, N_DYN
      REAL*8              :: P_TP1(IIPAR,JJPAR)
      REAL*8              :: P_TP2(IIPAR,JJPAR)
      REAL*8              :: UTMP(IIPAR,JJPAR,LLPAR)
      REAL*8              :: VTMP(IIPAR,JJPAR,LLPAR)
      REAL*8              :: WW(IIPAR,JJPAR,LLPAR) 

      ! Parameters
      INTEGER, PARAMETER  :: IGD=0, J1=3
      REAL*8,  PARAMETER  :: Umax=150d0

      !=================================================================
      ! DO_WINDOW_TRANSPORT begins here!
      !=================================================================

      ! Get nested-grid lon/lat offsets [# boxes]
      I0    = GET_XOFFSET( GLOBAL=.TRUE. )
      J0    = GET_YOFFSET( GLOBAL=.TRUE. )

      ! Dynamic timestep [s]
      N_DYN = GET_TS_DYN() * 60
       
      ! Impose TPCORE boundary conditions @ edges of 1x1 grid
      CALL DO_WINDOW_TPCORE_BC 

      ! Flip UWND, VWND, STT in vertical dimension
      ! Now store into temp arrays to preserve UWND, VWND for diags
      UTMP(:,:,1:LLPAR)          = UWND(:,:,LLPAR:1:-1)
      VTMP(:,:,1:LLPAR)          = VWND(:,:,LLPAR:1:-1)
      STT (:,:,1:LLPAR,:)        = STT (:,:,LLPAR:1:-1,:)

      ! Set temp arrays for passing pressures to TPCORE
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         P_TP1(I,J) = GET_PEDGE(I,J,1) - PTOP
         P_TP2(I,J) = PSC2(I,J) - PTOP
      ENDDO
      ENDDO

      ! Call the nested-grid window version of TPCORE (v.7.1)
      ! Use the pressures at the middle and the end of the 
      ! dynamic timestep (P = PS-PTOP; P_TEMP = PSC2-PTOP).
      CALL TPCORE_WINDOW( IGD,   STT,       P_TP1, P_TP2,  UTMP,  
     &                    VTMP,  WW,        N_DYN, IORD,   JORD,   
     &                    KORD,  N_TRACERS, IIPAR, JJPAR,  J1,    
     &                    I0,    J0,        I0_W,  J0_W,   I1_W,  
     &                    J1_W,  I2_W,      J2_W,  IM_W,   JM_W,  
     &                    IGZD,  LLPAR,     AP,    BP,     PTOP,  
     &                    Re,    LFILL,     LMFCT, Umax )

      ! Reset floating pressure w/ output of TPCORE
      CALL SET_FLOATING_PRESSURE( P_TP2 + PTOP )

      ! Re-Flip STT in the vertical dimension
      STT(:,:,1:LLPAR,:) = STT(:,:,LLPAR:1:-1,:)
              
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### MAIN: a TPCORE_WINDOW' )

      ! Return to calling program
      END SUBROUTINE DO_WINDOW_TRANSPORT

!------------------------------------------------------------------------------

      FUNCTION GET_AIR_MASS( I, J, L, P_SURF ) RESULT( AIR_MASS )
!
!******************************************************************************
!  Function GET_AIR_MASS returns the air mass based on the pressures returned
!  before and after the call to the GEOS-4/fvDAS TPCORE code. (bmy, 6/24/03)
! 
!  Arguments as Input:
!  ============================================================================
!  (1-3) I, J, L (INTEGER) : Grid box longitude, latitude, altitude indices 
!  (4  ) PSURF   (REAL*8 ) : Surface pressure at grid box (I,J,L=1)
!
!  NOTES:
!******************************************************************************
!     
#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"  ! g0_100

      ! Arguments
      INTEGER, INTENT(IN) :: I, J, L
      REAL*8,  INTENT(IN) :: P_SURF

      ! Local variables
      INTEGER             :: L2
      REAL*8              :: P_BOT, P_TOP, AIR_MASS
      
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

      ! Return to calling program
      END FUNCTION GET_AIR_MASS

!------------------------------------------------------------------------------

      SUBROUTINE SET_TRANSPORT( I_ORD, J_ORD, K_ORD )
!
!******************************************************************************
!  Subroutine SET_TRANSPORT passes IORD, JORD, KORD values from "input_mod.f"
!  to "transport_mod.f". (bmy, 7/20/04)
!
!  Arguments as Input:
!  ============================================================================
!  (1 ) I_ORD (INTEGER) : TPCORE E/W     transport option flag
!  (2 ) J_ORD (INTEGER) : TPCORE N/S     transport option flag
!  (3 ) K_ORD (INTEGER) : TPCORE up/down transport option flag
!
!  NOTES:
!******************************************************************************
!
      ! Arguments
      INTEGER, INTENT(IN) :: I_ORD, J_ORD, K_ORD

      !=================================================================
      ! SET_TRANSPORT begins here!
      !=================================================================

      ! Assign module variables
      IORD = I_ORD
      JORD = J_ORD
      KORD = K_ORD 

      ! Return to calling program
      END SUBROUTINE SET_TRANSPORT

!------------------------------------------------------------------------------

      SUBROUTINE INIT_TRANSPORT
!
!******************************************************************************
!  Subroutine INIT_TRANSPORT initializes all module variables and arrays. 
!  (bmy, 3/10/03, 10/30/07)
!
!  NOTES:
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
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,        ONLY : ALLOC_ERR
      USE GRID_MOD,         ONLY : GET_AREA_M2, GET_YMID_R
      USE LOGICAL_MOD,      ONLY : LEMBED,      LTPFV,     LTRAN
      USE PRESSURE_MOD,     ONLY : GET_AP,      GET_BP
      USE TIME_MOD,         ONLY : GET_TS_DYN
      USE TPCORE_FVDAS_MOD, ONLY : INIT_TPCORE
      USE TRACER_MOD,       ONLY : N_TRACERS

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"  ! Re

      ! Local variables
      INTEGER             :: AS, J, K, L, N_DYN
      REAL*8              :: YMID_R(JJPAR)

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
      ALLOCATE( Ap( LLPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Ap' )
 
      ALLOCATE( Bp( LLPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Bp' )

      ! Flip Ap and Bp for TPCORE
      DO L = 1, LLPAR+1 

         ! As L runs from the surface up, 
         ! K runs from the top down
         K = ( LLPAR + 1 ) - L + 1

#if   defined( GEOS_3 )
         Ap(L) = GET_AP(K) / PTOP   ! Ap(L) = 1 for all levels L
#else
         Ap(L) = GET_AP(K)          ! Ap(L) is in [hPa] 
#endif

         Bp(L) = GET_BP(K)
      ENDDO

      !=================================================================
      ! Allocate arrays for surface area and layer thickness
      !=================================================================
      ALLOCATE( A_M2( JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'A_M2' )

      ! Surface area [m2]
      DO J = 1, JJPAR
         A_M2(J) = GET_AREA_M2( J )
      ENDDO

      !=================================================================
      ! Additional setup for the GEOS-4/fvDAS version of TPCORE
      !=================================================================
#if   !defined( GEOS_3 )

      ! Initialize
      N_DYN = GET_TS_DYN() * 60
      N_ADJ = 0
      NG    = 0
      MG    = 0

      ! YMID_R is latitude of grid box center [radians]
      DO J = 1,JJPAR
         YMID_R(J) = GET_YMID_R(J)
      ENDDO

      ! Call INIT routine from "tpcore_fvdas_mod.f" 
      CALL INIT_TPCORE( IIPAR,  JJPAR, LLPAR,  JFIRST, JLAST, 
     &                  NG, MG, DBLE( N_DYN ), Re,     YMID_R )

#endif

      ! Return to calling program
      END SUBROUTINE INIT_TRANSPORT

!------------------------------------------------------------------------------

      SUBROUTINE INIT_GEOS5_WINDOW_TRANSPORT
!
!******************************************************************************
!  Subroutine INIT_GEOS5_WINDOW_TRANSPORT initializes all module variables and 
!  arrays for the GEOS-5 nested grid simulation.  This routine is only called
!  if we are using the GEOS-5 nested grid simulation. (dan, bmy, 11/6/08)
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,               ONLY : ALLOC_ERR
      USE GRID_MOD,                ONLY : GET_AREA_M2
      USE GRID_MOD,                ONLY : GET_YMID_R_W
      USE LOGICAL_MOD,             ONLY : LEMBED, LTPFV, LTRAN
      USE PRESSURE_MOD,            ONLY : GET_AP, GET_BP
      USE TIME_MOD,                ONLY : GET_TS_DYN
      USE TPCORE_FVDAS_MOD,        ONLY : INIT_TPCORE
      USE TPCORE_BC_MOD,           ONLY : I0_W, J0_W, I1_W, J1_W
      USE TPCORE_BC_MOD,           ONLY : I2_W, J2_W, IM_W, JM_W
      USE TPCORE_BC_MOD,           ONLY : IGZD, INIT_TPCORE_BC
      USE TPCORE_GEOS5_WINDOW_MOD, ONLY : INIT_GEOS5_WINDOW
      USE TRACER_MOD,              ONLY : N_TRACERS

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN_GCTM"  ! Re

      ! Local variables
      INTEGER             :: AS, J, K, L, N_DYN
      REAL*8              :: YMID_R_W(0:JJPAR+1)

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
      ALLOCATE( Ap( LLPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Ap' )

      ALLOCATE( Bp( LLPAR+1 ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'Bp' )

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
      ALLOCATE( A_M2( JJPAR ), STAT=AS )
      IF ( AS /= 0 ) CALL ALLOC_ERR( 'A_M2' )

      ! Surface area [m2]
      DO J = 1, JJPAR
         A_M2(J) = GET_AREA_M2( J )
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
         YMID_R_W(J) = GET_YMID_R_W(J)
      ENDDO

      ! Call INIT routine from "tpcore_geos5_fvdas_window_mod.f"
      CALL INIT_GEOS5_WINDOW( IIPAR, JJPAR,    LLPAR, JFIRST, 
     &                        JLAST, NG,       MG,    DBLE( N_DYN ), 
     &                        Re,    YMID_R_W )

      ! Return to calling program
      END SUBROUTINE INIT_GEOS5_WINDOW_TRANSPORT

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_TRANSPORT
!
!******************************************************************************
!  Subroutine CLEANUP_TRANSPORT deallocates all module arrays. 
!  (bmy, 3/10/03, 10/30/07)
!
!  NOTES:
!  (1 ) Remove reference to DSIG, it's obsolete. (bmy, 6/24/03)
!  (2 ) Remove obsolete embedded chemistry arrays (bmy, 10/30/07)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_TRANSPORT begins here!
      !=================================================================
      IF ( ALLOCATED( Ap     ) ) DEALLOCATE( Ap     )
      IF ( ALLOCATED( A_M2   ) ) DEALLOCATE( A_M2   )
      IF ( ALLOCATED( Bp     ) ) DEALLOCATE( Bp     )

      ! Return to calling program
      END SUBROUTINE CLEANUP_TRANSPORT

!------------------------------------------------------------------------------

      END MODULE TRANSPORT_MOD
