! $Id: transport_mod.f,v 1.13 2006/08/14 17:58:20 bmy Exp $
      MODULE TRANSPORT_MOD
!
!******************************************************************************
!  Module TRANSPORT_MOD is used to call the proper version of TPCORE for
!  GEOS-1, GEOS-STRAT, GEOS-3 or GEOS-4 nested-grid or global simulations.
!  (yxw, bmy, 3/10/03, 7/12/06)
! 
!  Module Variables:
!  ============================================================================
!  (1 ) Ap     (REAL*8 )     : Vertical coordinate array for TPCORE
!  (2 ) A_M2   (REAL*8 )     : Grid box surface areas [m2]
!  (3 ) Bp     (REAL*8 )     : Vertical coordinate array for TPCORE
!  (4 ) IORD   (REAL*8 )     : TPCORE E/W option flag
!  (5 ) JORD   (REAL*8 )     : TPCORE N/S option flag
!  (6 ) KORD   (REAL*8 )     : TPCORE vertical option flag
!  (7 ) STT_I1 (REAL*8 )     : Array for NH embedded chem boundary condition 
!  (8 ) STT_I2 (REAL*8 )     : Array for NH embedded chem boundary condition
!  (9 ) STT_J1 (REAL*8 )     : Array for NH embedded chem boundary condition
!  (10) STT_J2 (REAL*8 )     : Array for NH embedded chem boundary condition
!  (11) JLAST  (INTEGER)     : For GEOS-4 TPCORE
!  (12) MG     (INTEGER)     : For GEOS-4 TPCORE
!  (13) NG     (INTEGER)     : For GEOS-4 TPCORE
!  (14) N_ADJ  (INTEGER)     : For GEOS-4 TPCORE
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_TRANSPORT         : Driver which calls global or window TPCORE
!  (2 ) DO_GLOBAL_TRANSPORT  : Calls either GEOS-1/S/3 or fvDAS TPCORE (global)
!  (3 ) DO_WINDOW_TRANSPORT  : Calls nested-grid window version of TPCORE
!  (4 ) GET_AIR_MASS         : Computes air mass from TPCORE in/out pressures
!  (5 ) SET_TRANSPORT        : Gets IORD, JORD, KORD values from "input_mod.f"
!  (6 ) INIT_TRANSPORT       : Initializes module arrays
!  (7 ) CLEANUP_TRANSPORT    : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by transport_mod.f
!  ============================================================================
!  (1 ) dao_mod.f            : Module containing arrays for DAO met fields
!  (2 ) diag_mod.f           : Module containing GEOS-CHEM diagnostic arrays
!  (3 ) error_mod.f          : Module containing I/O error/NaN check routines
!  (4 ) grid_mod.f           : Module containing horizontal grid information
!  (5 ) logical_mod.f        : Module containing GEOS-CHEM logical switches
!  (6 ) pjc_pfix_mod.f       : Module containing Phil Cameron-Smith P-fixer
!  (7 ) pressure_mod.f       : Module containing routines to compute P(I,J,L)
!  (8 ) time_mod.f           : Module containing routines to compute date/time
!  (9 ) tpcore_mod.f         : Module containing TPCORE for GEOS1,GEOSS,GEOS3
!  (10) tpcore_bc_mod.f      : Module containing TPCORE boundary cond. routines
!  (11) tpcore_window_mod.f  : Module containing TPCORE for nested-grid windows
!  (12) tpcore_fvdas_mod.f90 : Module containing TPCORE for GEOS-4/fvDAS
!  (13) tracer_mod.f         : Module containing GEOS-CHEM tracer array STT etc
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
      PUBLIC :: SET_TRANSPORT

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      LOGICAL             :: USE_GEOS_4_TRANSPORT
      INTEGER             :: IORD,  JORD, KORD, JFIRST 
      INTEGER             :: JLAST, NG,   MG,   N_ADJ
      REAL*8, ALLOCATABLE :: Ap(:)
      REAL*8, ALLOCATABLE :: A_M2(:)
      REAL*8, ALLOCATABLE :: Bp(:)
      REAL*8, ALLOCATABLE :: DSIG(:)
      REAL*8, ALLOCATABLE :: STT_I1(:,:,:) 
      REAL*8, ALLOCATABLE :: STT_I2(:,:,:) 
      REAL*8, ALLOCATABLE :: STT_J1(:,:,:) 
      REAL*8, ALLOCATABLE :: STT_J2(:,:,:) 
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
!  for GEOS-3, GEOS-4, or window simulations. (bmy, 3/10/03, 7/20/04)
! 
!  NOTES:
!  (1 ) Removed IORD, JORD, KORD from the arg list.  Also now removed
!        reference to CMN, it's not needed. (bmy, 7/20/04)
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : ITS_A_NESTED_GRID

      ! Local variables
      LOGICAL, SAVE :: FIRST = .TRUE.

      !=================================================================
      ! DO_TRANSPORT begins here!
      !=================================================================

      ! First-time initialization
      IF ( FIRST ) THEN 
         CALL INIT_TRANSPORT
         FIRST = .FALSE.
      ENDIF

      ! Choose the proper version of TPCORE for the nested-grid window 
      ! region (usually 1x1 grids) or for the entire globe
      IF ( ITS_A_NESTED_GRID() ) THEN
         CALL DO_WINDOW_TRANSPORT
      ELSE
         CALL DO_GLOBAL_TRANSPORT
      ENDIF

      ! Return to calling program
      END SUBROUTINE DO_TRANSPORT

!------------------------------------------------------------------------------
      
      SUBROUTINE DO_GLOBAL_TRANSPORT
!
!******************************************************************************
!  Subroutine DO_TRANSPORT is the driver routine for the proper TPCORE 
!  program for GEOS-1, GEOS-STRAT, GEOS-3 or GEOS-4 global simulations. 
!  (bdf, bmy, 3/10/03, 7/12/06)
! 
!  NOTES:
!  (1 ) Now references routine TPCORE_FVDAS from "tpcore_fvdas_mod.f90".
!        Also now use logical flag USE_GEOS_4_TRANSPORT to decide which
!        version of TPCORE is used.  Now call routine DO_PJC_PFIX from
!        "pjc_pfix_mod.f" which calls the Phil Cameron-Smith pressure fixer
!        for the GEOS-4/fvDAS transport scheme. (bdf, bmy, 5/8/03)
!  (2 ) Now call GET_AIR_MASS to compute air masses based on the input/output
!        pressures of the GEOS-4 version of TPCORE (bmy, 6/24/03)
!  (3 ) Now references DEBUG_MSG from "error_mod.f". (bmy, 8/7/03)
!  (4 ) Bug fix: rewrote first parallel DO-loop to avoid NaN's.  Now also make
!        sure to pass surface pressures which are consistent with the Ap and
!        Bp coordinates which define the vertical grid to both TPCORE and
!        DO_PJC_PFIX. (bmy, 10/27/03)
!  (5 ) Removed IORD, JORD, KORD from the arg list, since these are now
!        module variables.  Now get LFILL, LMFCT, LPRT, LEMBED, LWINDO from 
!        "logical_mod.f".  Now references STT, N_TRACERS, TCVV from 
!        "tracer_mod.f".  Now parallelized embedded chemistry BC's.  
!        (bmy, 7/20/04) 
!  (6 ) Now references MASSFLEW, MASSFLNS, MASSFLUP from "diag_mod.f".
!        Also references ND24, ND25, ND26 from "CMN_DIAG". (bmy, 9/28/04)
!  (7 ) For GEOS-3 transport, we don't have to flip the STT array before & 
!        after the call to TPCORE because we now call TPCORE with the array 
!        mask statement STT(:,:,LLPAR:1:-1,:).  Also modified for GEOS-5 
!        and GCAP met fields. (swu, bmy, 5/25/05)
!  (8 ) Now make sure all USE statements are USE, ONLY (bmy, 10/3/05)
!  (9 ) Now do flipping of arrays in call to TPCORE_FVDAS (bmy, 6/16/06)
!  (10) Rewrote some parallel loops for the SUN compiler (bmy, 7/14/06)
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
      USE TPCORE_BC_MOD,    ONLY : SAVE_GLOBAL_TPCORE_BC
      USE TPCORE_MOD,       ONLY : TPCORE
      USE TPCORE_FVDAS_MOD, ONLY : TPCORE_FVDAS
      USE TRACER_MOD,       ONLY : N_TRACERS, STT, TCVV

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! IEBD1, IEBD2, JEBD1, JEBD2
#     include "CMN_GCTM"  ! Re, G0_100
#     include "CMN_DIAG"  ! ND24, ND25, ND26
      
      ! Local variables
      INTEGER             :: I, J, L, L2, N, N_DYN
      REAL*8              :: A_DIFF, D_DYN, TR_DIFF
      REAL*8              :: AD_A(IIPAR,JJPAR,LLPAR)
      REAL*8              :: AD_B(IIPAR,JJPAR,LLPAR)
      REAL*8              :: P_TP1(IIPAR,JJPAR)
      REAL*8              :: P_TP2(IIPAR,JJPAR)
      REAL*8              :: P_TEMP(IIPAR,JJPAR)
      REAL*8              :: TR_A(IIPAR,JJPAR,LLPAR)
      REAL*8              :: TR_B(IIPAR,JJPAR,LLPAR,N_TRACERS)
      REAL*8              :: UTMP(IIPAR,JJPAR,LLPAR)
      REAL*8              :: VTMP(IIPAR,JJPAR,LLPAR)
      REAL*8              :: WW(IIPAR,JJPAR,LLPAR) 
      REAL*8              :: XMASS(IIPAR,JJPAR,LLPAR) 
      REAL*8              :: YMASS(IIPAR,JJPAR,LLPAR) 

      ! Parameters
      INTEGER, PARAMETER  :: IGD=0, J1=3
      REAL*8,  PARAMETER  :: Umax=200d0

      !=================================================================
      ! DO_GLOBAL_TRANSPORT begins here!
      !=================================================================

      ! Save boundary conditions (global grid) for future nested run
      IF ( LWINDO ) CALL SAVE_GLOBAL_TPCORE_BC

      !=================================================================
      ! Set boundary conditions for embedded chemistry.  For latitude 
      ! boxes adjacent to the embedded box, chemistry is turned off, 
      ! but we still need reasonable tracer concentrations, because 
      ! transport is still turned on for the entire globe.  We keep 
      ! the same initial conditions for these adjacent boxes for the 
      ! duration of the run.
      !=================================================================
      IF ( LEMBED ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
         DO L = 1, LLPAR

            !-------------------------
            ! E/W boundary conditions
            !-------------------------
            DO J = 1, JJPAR
               IF ( IEBD1 > 1 ) THEN
                  STT_I1(J,L,N) = STT(IEBD1-1,J,L,N)
               ENDIF
         
               IF ( IEBD2 < IIPAR ) THEN
                  STT_I2(J,L,N) = STT(IEBD2+1,J,L,N)
               ENDIF
            ENDDO

            !-------------------------
            ! N/S boundary conditions
            !-------------------------
            DO I = 1, IIPAR              
               IF ( JEBD1 > 1  ) THEN
                  STT_J1(I,L,N) = STT(I,JEBD1-1,L,N)
               ENDIF
         
               IF ( JEBD2 < JJPAR ) THEN
                  STT_J2(I,L,N) = STT(I,JEBD2+1,L,N)
               ENDIF
            ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF

      !=================================================================
      ! Prepare variables for calls to PJC P-fixer and TPCORE
      !
      ! For GEOS-4 (hybrid grid), the pressure at the bottom edge 
      ! grid box (I,J,L) is given by:
      !
      !    P(I,J,L) = Ap(L) + Bp(L) * Psurf(I,J)
      !
      ! where:
      !
      !    Ap(L) and Bp(L) are defined in "pressure_mod.f"
      !    Psurf(I,J) = "true" surface pressure at surface box (I,J)
      !
      ! However, for GEOS-1, GEOS-STRAT, and GEOS-3, these are pure
      ! sigma grids, and the pressure at the bottom edge of level L
      ! is given by:
      !
      !    P(I,J,L) = Ap(L) + Bp(L) * ( Psurf(I,J) - PTOP )
      !
      ! where:
      !
      !    Ap(L)      = PTOP for all L
      !    Bp(L)      = bottom sigma edge of level L
      !    Psurf(I,J) = "true" surface pressure at surface box (I,J)
      !    PTOP       = model top pressure
      !
      ! When passing pressures to TPCORE, we must make sure that they
      ! are consistent with the definition of the corresponding Ap and
      ! Bp vertical coordinates.  This means:
      !
      !    GEOS-4    : pass Psurf(I,J)            to TPCORE
      !    GEOS-3    : pass ( Psurf(I,J) - PTOP ) to TPCORE
      !    GEOS-STRAT: pass ( Psurf(I,J) - PTOP ) to TPCORE
      !    GEOS-1    : pass ( Psurf(I,J) - PTOP ) to TPCORE
      !
      ! where Psurf(I,J) is the true surface pressure at box (I,J) 
      ! and PTOP the model top pressure.  
      !
      ! Also, the PJC P-fixer driver routine, DO_PJC_PFIX, now accepts 
      ! the true surface pressure instead of Psurf-PTOP.  This means:
      !
      !    GEOS-4    : pass P_TP1,      P_TP2      to DO_PJC_PFIX
      !    GEOS-3    : pass P_TP1+PTOP, P_TP2+PTOP to DO_PJC_PFIX
      !=================================================================
       
      ! Dynamic timestep [s]
      N_DYN = GET_TS_DYN() * 60
      D_DYN = DBLE( N_DYN )

      ! P_TP1 = PS - PTOP at the middle of the dynamic timestep
      ! P_TP2 = PS - PTOP at the end    of the dynamic timestep      
      DO J = 1, JJPAR
      DO I = 1, IIPAR

#if defined( GEOS_4 ) || defined( GEOS_5 )

         ! *** For GEOS-4, GEOS-5 winds ***
         ! We need to have P_TP1 and P_TP2 as the true surface pressure, 
         ! in order to be consistent with the Ap and Bp coordinates which 
         ! define the GEOS-4 hybrid grid.
         P_TP1(I,J) = GET_PEDGE(I,J,1)
         P_TP2(I,J) = PSC2(I,J)       

#else 

         ! *** For GCAP, GEOS-3, GEOS-STRAT, GEOS-1 winds *** 
         ! We need to have P_TP1 and P_TP2 to be ( true sfc pressure - PTOP )
         ! in order to be consistent with the Ap and Bp coordinates which 
         ! define the pure-sigma grid.  
         P_TP1(I,J) = GET_PEDGE(I,J,1) - PTOP
         P_TP2(I,J) = PSC2(I,J)        - PTOP

#endif
      ENDDO
      ENDDO
      
      ! Select proper version of TPCORE
      IF ( USE_GEOS_4_TRANSPORT ) THEN

         !==============================================================
         ! Use GEOS-4/fvDAS version of TPCORE
         ! (compatible with GEOS-3, GEOS-4, or GCAP winds)
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

         !----------------------------
         ! Apply pressure fixer 
         !----------------------------

#if   defined( GEOS_4 ) || defined( GEOS_5 )

         ! *** For GEOS-4 and GEOS-5 winds ***
         ! Call PJC P-fixer to get adjusted air masses (XMASS, YMASS)
         ! Pass "true" surface pressures P_TP1 and P_TP2 to DO_PJC_PFIX.
         CALL DO_PJC_PFIX( D_DYN, P_TP1, P_TP2, 
     &                     UWND,  VWND,  XMASS, YMASS )

#else
         ! *** For GCAP and GEOS-3 winds ***
         ! Call PJC P-fixer to get adjusted air masses (XMASS, YMASS)
         ! P_TP1 and P_TP2 are ( Psurface - PTOP ), so we must call 
         ! DO_PJC_PFIX with ( P_TP1 + PTOP ) and ( P_TP2 + PTOP ). 
         CALL DO_PJC_PFIX( D_DYN, P_TP1+PTOP, P_TP2+PTOP, 
     &                     UWND,  VWND,       XMASS, YMASS )

#endif

         !----------------------------
         ! Call transport code
         !----------------------------

         ! Flip arrays in vertical dimension for TPCORE
         ! Store winds in UTMP, VTMP to preserve UWND, VWND for diagnostics
         UTMP(:,:,1:LLPAR) = UWND(:,:,LLPAR:1:-1)
         VTMP(:,:,1:LLPAR) = VWND(:,:,LLPAR:1:-1)

         ! GEOS-4/fvDAS transport (the output pressure is P_TEMP)
         ! NOTE: P_TP1 and P_TP2 must be consistent with the definition
         ! of Ap and Bp.  For GEOS-4, P_TP1 and P_TP2 must be the "true"
         ! surface pressure, but for GEOS-3, they must be ( Psurface -PTOP ).  
         CALL TPCORE_FVDAS( D_DYN,    Re,        IIPAR,    JJPAR,
     &                      LLPAR,    JFIRST,    JLAST,    NG,
     &                      MG,       N_TRACERS, Ap,       Bp,
     &                      UTMP,     VTMP,      P_TP1,    P_TP2,
     &                      P_TEMP,   STT(:,:,LLPAR:1:-1,:),       
     &                      IORD,     JORD,      KORD,     N_ADJ,     
     &                      XMASS(:,:,LLPAR:1:-1),    
     &                      YMASS(:,:,LLPAR:1:-1),
     &                      MASSFLEW(:,:,LLPAR:1:-1,:), 
     &                      MASSFLNS(:,:,LLPAR:1:-1,:),  
     &                      MASSFLUP(:,:,LLPAR:1:-1,:),    A_M2,
     &                      TCVV,     ND24,      ND25,     ND26 )

         !----------------------------
         ! Reset surface pressure
         !----------------------------

#if   defined( GEOS_4 ) || defined( GEOS_5 )

         ! *** For GEOS-4 or GEOS-5 winds ***
         ! P_TP2 is the "true" surface pressure at the end of the 
         ! dynamic timestep.  Reset the pressure with P_TP2.  This will 
         ! be the pressure at the start of the next dynamic timestep.
         CALL SET_FLOATING_PRESSURE( P_TP2 )
#else

         ! *** For GCAP and GEOS-3 winds ***
         ! P_TP2 is the "true" surface pressure at the end of the dynamic 
         ! timestep - PTOP.  Reset the pressure with P_TP2+PTOP.  This 
         ! will be the pressure at the start of the next dynamic timestep.
         CALL SET_FLOATING_PRESSURE( P_TP2 + PTOP )

#endif

         ! Adjust tracer to correct residual non-conservation of mass
         DO N = 1, N_TRACERS

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Air mass [kg] after transport
               IF ( N == 1 ) THEN
                  AD_A(I,J,L) = GET_AIR_MASS( I, J, L, P_TP2(I,J) )
               ENDIF
         
               ! Tracer mass [kg] after transport
               TR_A(I,J,L) = STT(I,J,L,N) * AD_A(I,J,L) / TCVV(N)
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO

            ! Residual mass difference [kg]: before - after
            TR_DIFF = SUM( TR_B(:,:,:,N) ) - SUM( TR_A )

            ! Convert from [kg] to [v/v]
            TR_DIFF = TR_DIFF / SUM( AD_A ) * TCVV(N)

            ! Add mass difference [v/v] back to STT
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR
               STT(I,J,L,N) = STT(I,J,L,N) + TR_DIFF               
               STT(I,J,L,N) = MAX( STT(I,J,L,N), 0d0 )
            ENDDO
            ENDDO
            ENDDO
!$OMP END PARALLEL DO
         ENDDO

      ELSE
      
         !==============================================================
         ! Use TPCORE version 7.1.m
         ! (compatible with GEOS-1, GEOS-STRAT, or GEOS-3)
         !==============================================================

         ! Flip arrays in vertical dimension
         ! Store winds in UTMP, VTMP to preserve UWND, VWND for diagnostics
         UTMP(:,:,1:LLPAR  ) = UWND(:,:,LLPAR:1:-1  )
         VTMP(:,:,1:LLPAR  ) = VWND(:,:,LLPAR:1:-1  )

         ! TPCORE v7.1.m transport scheme (output pressure is P_TP2)
         ! The pressures P_TP1 and P_TP2 are PS-PTOP, in order to
         ! be consistent with the definition of Ap and Bp for GEOS-3
         ! GEOS-STRAT, and GEOS-1 winds. (bmy, 10/27/03)
         CALL TPCORE( IGD,   STT(:,:,LLPAR:1:-1,:),
     &                P_TP1, P_TP2, UTMP, VTMP,  WW,    
     &                N_DYN, IORD,  JORD, KORD,  N_TRACERS, 
     &                IIPAR, JJPAR, J1,   LLPAR, Ap,   
     &                Bp,    PTOP,  Re,   LFILL, LMFCT, Umax )

         ! Reset floating pressure w/ pressure adjusted by TPCORE.  Here
         ! P_TP2 is PS-PTOP, so reset the pressure with P_TP2+PTOP. This
         ! will be the pressure at the start of the next dynamic timestep.
         CALL SET_FLOATING_PRESSURE( P_TP2 + PTOP )
           
      ENDIF

      !=================================================================
      ! Reset the boundary conditions at adjacent boxes for embedded 
      ! chemistry, so transport in the non-chemistry southern hemisphere 
      ! does not transport weird things into the northern hemisphere.
      !=================================================================
      IF ( LEMBED ) THEN

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
         DO N = 1, N_TRACERS
         DO L = 1, LLPAR

            !-------------------------
            ! E/W boundary conditions
            !-------------------------
            DO J = 1, JJPAR
               IF ( IEBD1 > 1  ) THEN
                  STT(IEBD1-1,J,L,N) = STT_I1(J,L,N)
               ENDIF
         
               IF ( IEBD2 < IIPAR ) THEN
                  STT(IEBD2+1,J,L,N) = STT_I2(J,L,N)
               ENDIF
            ENDDO

            !-------------------------
            ! E/W boundary conditions
            !-------------------------
            DO I = 1, IIPAR
               IF ( JEBD1 > 1  ) THEN
                  STT(I,JEBD1-1,L,N) = STT_J1(I,L,N)
               ENDIF
               
               IF ( JEBD2 < JJPAR ) THEN
                  STT(I,JEBD2+1,L,N) = STT_J2(I,L,N)
               ENDIF
            ENDDO

         ENDDO
         ENDDO
!$OMP END PARALLEL DO

      ENDIF               
                    
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_GLOBAL_TRANSPORT: a TPCORE') 

      ! Return to calling program
      END SUBROUTINE DO_GLOBAL_TRANSPORT

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
!  (bmy, 3/10/03, 5/25/05)
!
!  NOTES:
!  (1 ) Now references GET_TS_DYN from "time_mod.f", INIT_TPCORE_FVDAS from
!        "tpcore_fvdas_mod.f90", and GET_YMID_R from "grid_mod.f".  Now also
!        include "CMN_SETUP".  (bdf, bmy, 4/28/03)
!  (2 ) Remove reference to DSIG, it's obsolete. (bmy, 6/24/03)
!  (3 ) Now references LEMBED & LTPFV from "logical_mod.f".  Now references
!        N_TRACERS from "tracer_mod.f". (bmy, 7/20/04)
!  (4 ) Now modified for GEOS-5 and GCAP met fields (swu, bmy, 5/25/05)
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
#     include "CMN"       ! IEBD1, JEBD1, IEBD2, JEBD2
#     include "CMN_GCTM"  ! Re

      ! Local variables
      INTEGER             :: AS, J, K, L, N_DYN
      REAL*8              :: YMID_R(JJPAR)

#if   defined( GEOS_4 ) || defined( GEOS_5 ) || defined( GCAP )

      ! For GEOS-4, GEOS-5, or GCAP winds, use the fvDAS transport routines
      USE_GEOS_4_TRANSPORT = .TRUE.

#elif defined( GEOS_3 ) 

      ! For GEOS-3 winds, select either the GEOS-4/fvDAS transport
      ! (if LTPFV=T) or the existing TPCORE 7.1m (if LTPFV=F)
      USE_GEOS_4_TRANSPORT = LTPFV

#else

      ! We can't use the GEOS-4/fvDAS transport for GEOS-1/GEOS-STRAT
      USE_GEOS_4_TRANSPORT = .FALSE.
      
#endif
      
      !=================================================================
      ! Allocate arrays for embedded chemistry boundary conditions
      !=================================================================
      IF ( LEMBED ) THEN 
         ALLOCATE( STT_I1( IIPAR, LLPAR, N_TRACERS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT_I1' )
         STT_I1 = 0d0

         ALLOCATE( STT_I2( IIPAR, LLPAR, N_TRACERS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT_I2' ) 
         STT_I2 = 0d0

         ALLOCATE( STT_J1( JJPAR, LLPAR, N_TRACERS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT_J1' ) 
         STT_J1 = 0d0

         ALLOCATE( STT_J2( JJPAR, LLPAR, N_TRACERS ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT_J2' ) 
         STT_J2 = 0d0
      ENDIF

      !=================================================================
      ! Allocate arrays for TPCORE vertical coordinates 
      !
      ! For TPCORE v7.1.m (for GEOS-1, GEOS-STRAT, GEOS-3 met fields):
      ! 
      !    P(I,J,L) = ( Ap(L) * PTOP ) + ( Bp(L) * ( Psurf(I,J)-PTOP ) )
      !
      ! For GEOS-4/fvDAS TPCORE (for GEOS-3 or GEOS-4 only):
      !
      !    P(I,J,L) = Ap(L) + ( Bp(L) * Psurf(I,J) )
      !
      ! Also here Ap, Bp will be flipped since both TPCORE versions
      ! index levels from the atm. top downwards (bdf, bmy, 7/20/04)
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

         IF ( USE_GEOS_4_TRANSPORT ) THEN
            Ap(L) = GET_AP(K)        ! Ap(L) is in [hPa]
         ELSE
            Ap(L) = GET_AP(K) / PTOP ! Ap(L) = 1 for all levels L
         ENDIF

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
      IF ( USE_GEOS_4_TRANSPORT ) THEN

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
     &                     NG, MG, DBLE( N_DYN ), Re,     YMID_R )

      ENDIF

      ! Return to calling program
      END SUBROUTINE INIT_TRANSPORT

!------------------------------------------------------------------------------

      SUBROUTINE CLEANUP_TRANSPORT
!
!******************************************************************************
!  Subroutine CLEANUP_TRANSPORT deallocates all module arrays. 
!  (bmy, 3/10/03, 6/24/03)
!
!  NOTES:
!  (1 ) Remove reference to DSIG, it's obsolete. (bmy, 6/24/03)
!******************************************************************************
!
      !=================================================================
      ! CLEANUP_TRANSPORT begins here!
      !=================================================================
      IF ( ALLOCATED( Ap     ) ) DEALLOCATE( Ap     )
      IF ( ALLOCATED( A_M2   ) ) DEALLOCATE( A_M2   )
      IF ( ALLOCATED( Bp     ) ) DEALLOCATE( Bp     )
      IF ( ALLOCATED( STT_I1 ) ) DEALLOCATE( STT_I1 )
      IF ( ALLOCATED( STT_I2 ) ) DEALLOCATE( STT_I2 )
      IF ( ALLOCATED( STT_J1 ) ) DEALLOCATE( STT_J1 )
      IF ( ALLOCATED( STT_J2 ) ) DEALLOCATE( STT_J2 )

      ! Return to calling program
      END SUBROUTINE CLEANUP_TRANSPORT

!------------------------------------------------------------------------------

      END MODULE TRANSPORT_MOD
