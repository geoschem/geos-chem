! $Id: transport_mod.f,v 1.1 2003/06/30 20:26:02 bmy Exp $
      MODULE TRANSPORT_MOD
!
!******************************************************************************
!  Module TRANSPORT_MOD is used to call the proper version of TPCORE for
!  GEOS-1, GEOS-STRAT, GEOS-3 or GEOS-4 nested-grid or global simulations.
!  (yxw, bmy, 3/10/03, 6/24/03)
! 
!  Module Variables:
!  ============================================================================
!  (1 ) Ap     (REAL*8)      : Vertical coordinate array for TPCORE
!  (2 ) A_M2   (REAL*8)      : Grid box surface areas [m2]
!  (3 ) Bp     (REAL*8)      : Vertical coordinate array for TPCORE
!  (4 ) STT_I1 (REAL*8)      : Array for NH embedded chem boundary condition 
!  (5 ) STT_I2 (REAL*8)      : Array for NH embedded chem boundary condition
!  (6 ) STT_J1 (REAL*8)      : Array for NH embedded chem boundary condition
!  (7 ) STT_J2 (REAL*8)      : Array for NH embedded chem boundary condition
!
!  Module Routines:
!  ============================================================================
!  (1 ) DO_TRANSPORT         : Driver which calls global or window TPCORE
!  (2 ) DO_GLOBAL_TRANSPORT  : Calls either GEOS-1/S/3 or fvDAS TPCORE (global)
!  (3 ) DO_WINDOW_TRANSPORT  : Calls nested-grid window version of TPCORE
!  (4 ) GET_AIR_MASS         : Computes air mass from TPCORE in/out pressures
!  (5 ) INIT_TPCORE_CALL     : Initializes module arrays
!  (6 ) CLEANUP_TPCORE_CALL  : Deallocates module arrays
!
!  GEOS-CHEM modules referenced by transport_mod.f
!  ============================================================================
!  (1 ) dao_mod.f            : Module containing arrays for DAO met fields
!  (2 ) error_mod.f          : Module containing I/O error/NaN check routines
!  (3 ) grid_mod.f           : Module containing horizontal grid information
!  (4 ) pjc_pfix_mod.f       : Module containing Phil Cameron-Smith P-fixer
!  (5 ) pressure_mod.f       : Module containing routines to compute P(I,J,L)
!  (6 ) time_mod.f           : Module containing routines to compute date/time
!  (7 ) tpcore_mod.f         : Module containing TPCORE for GEOS1,GEOSS,GEOS3
!  (8 ) tpcore_bc_mod.f      : Module containing TPCORE boundary cond. routines
!  (9 ) tpcore_window_mod.f  : Module containing TPCORE for nested-grid windows
!  (10) tpcore_fvdas_mod.f90 : Module containing TPCORE for GEOS-4/fvDAS
!
!  NOTES:
!  (1 ) Now can select transport scheme for GEOS-3 winds.  Added code for PJC 
!        pressure fixer. (bdf, bmy, 5/8/03)
!  (2 ) Now delete DSIG array, it's obsolete.  Also added new PRIVATE function 
!        GET_AIR_MASS to compute air masses from the input/output pressures
!        from the new GEOS-4/fvDAS TPCORE. (bmy, 6/24/03)
!******************************************************************************
!
      IMPLICIT NONE

      !=================================================================
      ! MODULE PRIVATE DECLARATIONS -- keep certain internal variables 
      ! and routines from being seen outside "transport_mod.f"
      !=================================================================

      ! PRIVATE module variables
      PRIVATE :: Ap, Bp, STT_I1, STT_I2, STT_J1, STT_J2
      PRIVATE :: JFIRST, JLAST,  NG,     MG,     N_ADJ
      !-----------------------------------------------------------------
      ! Prior to 6/24/03:
      !PRIVATE :: A_M2,   DSIG,   USE_GEOS_4_TRANSPORT
      !-----------------------------------------------------------------
      PRIVATE :: A_M2,   USE_GEOS_4_TRANSPORT

      ! PRIVATE module routines
      PRIVATE :: INIT_TRANSPORT
      PRIVATE :: DO_GLOBAL_TRANSPORT
      PRIVATE :: DO_WINDOW_TRANSPORT
      PRIVATE :: GET_AIR_MASS

      !=================================================================
      ! MODULE VARIABLES
      !=================================================================
      LOGICAL             :: USE_GEOS_4_TRANSPORT
      INTEGER             :: JFIRST, JLAST, NG, MG, N_ADJ
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
      
      SUBROUTINE DO_TRANSPORT( IORD, JORD, KORD )
!
!******************************************************************************
!  Subroutine DO_TRANSPORT is the driver routine for the proper TPCORE program
!  for GEOS-3, GEOS-4, or window simulations. (bmy, 3/10/03)
! 
!  Arguments as Input:
!  ===========================================================================
!  (1 ) IORD (INTEGER) : TPCORE E/W      transport option flag 
!  (2 ) JORD (INTEGER) : TPCORE N/S      transport option flag
!  (3 ) KORD (INTEGER) : TPCORE vertical transport option flag
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE GRID_MOD, ONLY : ITS_A_NESTED_GRID

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! LWINDO ??? -- replace later

      ! Arguments
      INTEGER, INTENT(IN) :: IORD, JORD, KORD

      ! Local variables
      LOGICAL, SAVE       :: FIRST = .TRUE. 

      !=================================================================
      ! DO_TRANSPORT begins here!
      !=================================================================

      ! First-time initialization 
      IF ( FIRST ) THEN
         CALL INIT_TRANSPORT
         FIRST = .FALSE.
      ENDIF

      !=================================================================
      ! Choose the proper version of TPCORE for the nested-grid window 
      ! region (usually 1x1 grids) or for the entire globe
      !=================================================================
      IF ( ITS_A_NESTED_GRID() ) THEN
         CALL DO_WINDOW_TRANSPORT( IORD, JORD, KORD )
      ELSE
         CALL DO_GLOBAL_TRANSPORT( IORD, JORD, KORD )
      ENDIF

      ! Return to calling program
      END SUBROUTINE DO_TRANSPORT

!------------------------------------------------------------------------------
      
      SUBROUTINE DO_GLOBAL_TRANSPORT( IORD, JORD, KORD )
!
!******************************************************************************
!  Subroutine DO_TRANSPORT is the driver routine for the proper TPCORE |
!  program for GEOS-1, GEOS-STRAT, GEOS-3 or GEOS-4 global simulations. 
!  (bdf, bmy, 3/10/03, 6/24/03)
! 
!  Arguments as Input:
!  ===========================================================================
!  (1 ) IORD (INTEGER) : TPCORE E/W      transport option flag 
!  (2 ) JORD (INTEGER) : TPCORE N/S      transport option flag
!  (3 ) KORD (INTEGER) : TPCORE vertical transport option flag
!
!  NOTES:
!  (1 ) Now references routine TPCORE_FVDAS from "tpcore_fvdas_mod.f90".
!        Also now use logical flag USE_GEOS_4_TRANSPORT to decide which
!        version of TPCORE is used.  Now call routine DO_PJC_PFIX from
!        "pjc_pfix_mod.f" which calls the Phil Cameron-Smith pressure fixer
!        for the GEOS-4/fvDAS transport scheme. (bdf, bmy, 5/8/03)
!  (2 ) Now call GET_AIR_MASS to compute air masses based on the input/output
!        pressures of the GEOS-4 version of TPCORE (bmy, 6/24/03)
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,          ONLY : PSC2, UWND, VWND
      USE PJC_PFIX_MOD,     ONLY : DO_PJC_PFIX
      USE PRESSURE_MOD,     ONLY : GET_PEDGE, SET_FLOATING_PRESSURE
      USE TIME_MOD,         ONLY : GET_TS_DYN
      USE TPCORE_BC_MOD,    ONLY : SAVE_GLOBAL_TPCORE_BC
      USE TPCORE_MOD,       ONLY : TPCORE
      USE TPCORE_FVDAS_MOD, ONLY : TPCORE_FVDAS
      USE ERROR_MOD,        ONLY : IT_IS_NAN

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! STT, NTRACE, LPRT, LWINDO, TCVV
#     include "CMN_GCTM"  ! Re, G0_100
#     include "CMN_SETUP" ! LFILL, LMFCT

      ! Arguments
      INTEGER, INTENT(IN) :: IORD, JORD, KORD

      ! Local variables
      INTEGER             :: I, J, L, N, N_DYN
      REAL*8              :: A_DIFF, D_DYN, TR_DIFF
      REAL*8              :: AD_A(IIPAR,JJPAR,LLPAR)
      REAL*8              :: AD_B(IIPAR,JJPAR,LLPAR)
      REAL*8              :: P_TP1(IIPAR,JJPAR)
      REAL*8              :: P_TP2(IIPAR,JJPAR)
      REAL*8              :: P_TEMP(IIPAR,JJPAR)
      REAL*8              :: TR_A(IIPAR,JJPAR,LLPAR)
      REAL*8              :: TR_B(IIPAR,JJPAR,LLPAR,NNPAR)
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
         IF ( IEBD1 > 1 ) THEN
            STT_I1(:,:,1:NTRACE) = STT((IEBD1-1),:,:,1:NTRACE)
         ENDIF
   
         IF ( IEBD2 < IIPAR ) THEN
            STT_I2(:,:,1:NTRACE) = STT((IEBD2+1),:,:,1:NTRACE)
         ENDIF

         IF ( JEBD1 > 1  ) THEN
            STT_J1(:,:,1:NTRACE) = STT(:,(JEBD1-1),:,1:NTRACE)
          ENDIF
    
         IF ( JEBD2 < JJPAR ) THEN
            STT_J2(:,:,1:NTRACE) = STT(:,(JEBD2+1),:,1:NTRACE)
         ENDIF
      ENDIF

      !=================================================================
      ! Prepare some variables for call to TPCORE
      !=================================================================
       
      ! Dynamic timestep [s]
      N_DYN = GET_TS_DYN() * 60
      D_DYN = DBLE( N_DYN )

      ! P_TP1 = PS - PTOP at the middle of the dynamic timestep
      ! P_TP2 = PS - PTOP at the end    of the dynamic timestep
      DO J = 1, JJPAR
      DO I = 1, IIPAR
         P_TP1(I,J) = GET_PEDGE(I,J,1) - PTOP
         P_TP2(I,J) = PSC2(I,J)        - PTOP
      ENDDO
      ENDDO

      ! Select proper version of TPCORE
      IF ( USE_GEOS_4_TRANSPORT ) THEN

         !==============================================================
         ! GEOS-3 or GEOS-4: Use GEOS-4/fvDAS version of TPCORE
         !==============================================================
!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L, N )
!$OMP+SCHEDULE( DYNAMIC )
         DO N = 1, NTRACE
         DO L = 1, LLPAR
         DO J = 1, JJPAR
         DO I = 1, IIPAR

            ! Airmass [kg] before transport
            IF ( N == 1 ) THEN
               !-------------------------------------------------------------
               ! Prior to 6/24/03:
               ! Call GET_AIR_MASS to get air mass based on P_TP1
               !AD_B(I,J,L) = P_TP1(I,J) * DSIG(L) * G0_100 * A_M2(J)
               !-------------------------------------------------------------
               AD_B(I,J,L) = GET_AIR_MASS( I, J, L, P_TP1(I,J) )
            ENDIF

            ! Tracer mass [kg] before transport
            TR_B(I,J,L,N) = STT(I,J,L,N) * AD_B(I,J,L) / TCVV(N)
         ENDDO
         ENDDO
         ENDDO
         ENDDO
!$OMP END PARALLEL DO

         ! Call PJC P-fixer -- get adjusted air masses (XMASS, YMASS)
         CALL DO_PJC_PFIX( D_DYN, P_TP1, P_TP2, 
     &                     UWND,  VWND,  XMASS, YMASS )

         ! Flip arrays in vertical dimension for TPCORE
         ! Store winds in UTMP, VTMP to preserve UWND, VWND for diagnostics
         UTMP (:,:,1:LLPAR)          = UWND (:,:,LLPAR:1:-1)
         VTMP (:,:,1:LLPAR)          = VWND (:,:,LLPAR:1:-1)
         XMASS(:,:,1:LLPAR)          = XMASS(:,:,LLPAR:1:-1)
         YMASS(:,:,1:LLPAR)          = YMASS(:,:,LLPAR:1:-1)
         STT  (:,:,1:LLPAR,1:NTRACE) = STT  (:,:,LLPAR:1:-1,1:NTRACE)

         ! GEOS-4/fvDAS transport (the output pressure is P_TEMP)
         CALL TPCORE_FVDAS( D_DYN,  Re,     IIPAR,   JJPAR,
     &                      LLPAR,  JFIRST, JLAST,   NG,
     &                      MG,     NTRACE, Ap,      Bp,
     &                      UTMP,   VTMP,   P_TP1,   P_TP2,
     &                      P_TEMP, STT(:,:,:,1:NTRACE),
     &                      IORD,   JORD,   KORD,    N_ADJ,
     &                      XMASS,  YMASS )

         ! Reset the floating pressure w/ the met field pressure,
         ! this is correct for the new GEOS-4/fvDAS transport
         CALL SET_FLOATING_PRESSURE( P_TP2 + PTOP )

         ! Re-Flip STT in the vertical dimension
         STT(:,:,1:LLPAR,1:NTRACE) = STT(:,:,LLPAR:1:-1,1:NTRACE)

         ! Adjust tracer to correct residual non-conservation of mass
         DO N = 1, NTRACE

!$OMP PARALLEL DO
!$OMP+DEFAULT( SHARED )
!$OMP+PRIVATE( I, J, L )
            DO L = 1, LLPAR
            DO J = 1, JJPAR
            DO I = 1, IIPAR

               ! Air mass [kg] after transport
               IF ( N == 1 ) THEN
                  !----------------------------------------------------------
                  ! Prior to 6/24/03:
                  ! Call GET_AIR_MASS to get air mass from P_TP2
                  !AD_A(I,J,L) = P_TP2(I,J) * DSIG(L) * G0_100 * A_M2(J)
                  !----------------------------------------------------------
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
         ! GEOS-1, GEOS-STRAT, or GEOS-3: Use TPCORE version 7.1.m
         !==============================================================

         ! Flip arrays in vertical dimension
         ! Store winds in UTMP, VTMP to preserve UWND, VWND for diagnostics
         UTMP (:,:,1:LLPAR)          = UWND (:,:,LLPAR:1:-1)
         VTMP (:,:,1:LLPAR)          = VWND (:,:,LLPAR:1:-1)
         STT  (:,:,1:LLPAR,1:NTRACE) = STT  (:,:,LLPAR:1:-1,1:NTRACE)

         ! TPCORE v7.1.m transport scheme (output pressure is P_TP2)
         CALL TPCORE( IGD,   STT(:,:,:,1:NTRACE), P_TP1,   
     &                P_TP2, UTMP, VTMP,  WW,     N_DYN, 
     &                IORD,  JORD, KORD,  NTRACE, IIPAR, 
     &                JJPAR, J1,   LLPAR, Ap,     Bp,     
     &                PTOP,  Re,   LFILL, LMFCT,  Umax )

         ! Reset floating pressure w/ pressure adjusted by TPCORE
         CALL SET_FLOATING_PRESSURE( P_TP2 + PTOP )
           
         ! Re-Flip STT in the vertical dimension
         STT(:,:,1:LLPAR,1:NTRACE) = STT(:,:,LLPAR:1:-1,1:NTRACE)

      ENDIF

      !=================================================================
      ! Reset the boundary conditions at adjacent boxes for embedded 
      ! chemistry, so transport in the non-chemistry southern hemisphere 
      ! does not transport weird things into the northern hemisphere.
      !=================================================================
      IF ( LEMBED ) THEN
         IF ( IEBD1 > 1  ) THEN
            STT((IEBD1-1),:,:,1:NTRACE) = STT_I1(:,:,1:NTRACE)
         ENDIF

         IF ( IEBD2 < IIPAR ) THEN
             STT((IEBD2+1),:,:,1:NTRACE) = STT_I2(:,:,1:NTRACE)
          ENDIF
 
         IF ( JEBD1 > 1  ) THEN
            STT(:,(JEBD1-1),:,1:NTRACE) = STT_J1(:,:,1:NTRACE)
         ENDIF

         IF ( JEBD2 < JJPAR ) THEN
            STT(:,(JEBD2+1),:,1:NTRACE) = STT_J2(:,:,1:NTRACE)
         ENDIF
      ENDIF               
                    
      !### Debug
      IF ( LPRT ) CALL DEBUG_MSG( '### DO_GLOBAL_TRANSPORT: a TPCORE') 

      ! Return to calling program
      END SUBROUTINE DO_GLOBAL_TRANSPORT

!------------------------------------------------------------------------------

      SUBROUTINE DO_WINDOW_TRANSPORT( IORD, JORD, KORD )
!
!******************************************************************************
!  Subroutine DO_WINDOW_TRANSPORT is the driver program for the proper TPCORE
!  program for nested-grid window simulations. (yxw, bmy, 3/10/03)
!
!  Arguments as Input:
!  ===========================================================================
!  (1 ) IORD (INTEGER) : TPCORE E/W      transport option flag 
!  (2 ) JORD (INTEGER) : TPCORE N/S      transport option flag
!  (3 ) KORD (INTEGER) : TPCORE vertical transport option flag
!
!  NOTES:
!******************************************************************************
!
      ! References to F90 modules
      USE DAO_MOD,           ONLY : PSC2, UWND, VWND
      USE GRID_MOD,          ONLY : GET_XOFFSET, GET_YOFFSET
      USE PRESSURE_MOD,      ONLY : GET_PEDGE,   SET_FLOATING_PRESSURE
      USE TIME_MOD,          ONLY : GET_TS_DYN
      USE TPCORE_BC_MOD
      USE TPCORE_WINDOW_MOD, ONLY : TPCORE_WINDOW

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! STT, NTRACE, LPRT, LWINDO
#     include "CMN_GCTM"  ! Re
#     include "CMN_SETUP" ! LFILL, LMFCT

      ! Arguments
      INTEGER, INTENT(IN) :: IORD, JORD, KORD

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
      STT (:,:,1:LLPAR,1:NTRACE) = STT (:,:,LLPAR:1:-1,1:NTRACE)

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
      CALL TPCORE_WINDOW( IGD,   STT(:,:,:,1:NTRACE),  P_TP1,  
     &                    P_TP2, UTMP,  VTMP,  WW,     N_DYN,  
     &                    IORD,  JORD,  KORD,  NTRACE, IIPAR, 
     &                    JJPAR, J1,    I0,    J0,     I0_W, 
     &                    J0_W,  I1_W,  J1_W,  I2_W,   J2_W, 
     &                    IM_W,  JM_W,  IGZD,  LLPAR,  AP,     
     &                    BP,    PTOP,  Re,    LFILL,  LMFCT,  
     &                    Umax )

      ! Reset floating pressure w/ output of TPCORE
      CALL SET_FLOATING_PRESSURE( P_TP2 + PTOP )

      ! Re-Flip STT in the vertical dimension
      STT(:,:,1:LM,1:NTRACE) = STT(:,:,LM:1:-1,1:NTRACE)
              
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

      SUBROUTINE INIT_TRANSPORT
!
!******************************************************************************
!  Subroutine INIT_TPCORE_CALL initializes all module variables and arrays. 
!  (bmy, 3/10/03, 6/24/03)
!
!  NOTES:
!  (1 ) Now references GET_TS_DYN from "time_mod.f", INIT_TPCORE_FVDAS from
!        "tpcore_fvdas_mod.f90", and GET_YMID_R from "grid_mod.f".  Now also
!        include "CMN_SETUP".  (bdf, bmy, 4/28/03)
!  (2 ) Remove reference to DSIG, it's obsolete. (bmy, 6/24/03)
!******************************************************************************
!
      ! References to F90 modules
      USE ERROR_MOD,        ONLY : ALLOC_ERR
      USE GRID_MOD,         ONLY : GET_AREA_M2, GET_YMID_R
      USE PRESSURE_MOD,     ONLY : GET_AP,      GET_BP
      USE TIME_MOD,         ONLY : GET_TS_DYN
      USE TPCORE_FVDAS_MOD, ONLY : INIT_TPCORE

#     include "CMN_SIZE"  ! Size parameters
#     include "CMN"       ! LEMBED, NTRACE
#     include "CMN_GCTM"  ! Re
#     include "CMN_SETUP" ! LTPFV

      ! Local variables
      INTEGER :: AS, J, K, L, N_DYN
      REAL*8  :: YMID_R(JJPAR)

      !=================================================================
      ! INIT_TRANSPORT begins here
      !=================================================================
 
#if   defined( GEOS_4 ) 

      ! For GEOS-4/fvDAS fields, use the fvDAS transport routines
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
         ALLOCATE( STT_I1( IIPAR, LLPAR, NTRACE ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT_I1' )
         STT_I1 = 0d0

         ALLOCATE( STT_I2( IIPAR, LLPAR, NTRACE ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT_I2' ) 
         STT_I2 = 0d0

         ALLOCATE( STT_J1( JJPAR, LLPAR, NTRACE ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT_J1' ) 
         STT_J1 = 0d0

         ALLOCATE( STT_J2( JJPAR, LLPAR, NNPAR ), STAT=AS )
         IF ( AS /= 0 ) CALL ALLOC_ERR( 'STT_J2' ) 
         STT_J2 = 0d0
      ENDIF

      !=================================================================
      ! Allocate arrays for TPCORE vertical coordinates 
      !
      ! For TPCORE v7.1.m (for GEOS-1, GEOS-STRAT, GEOS-3 met fields):
      ! 
      !    P(I,J,L) = ( Ap(L) * PTOP ) + ( Bp(L) * Psurface(I,J) )
      !
      ! For GEOS-4/fvDAS TPCORE (for GEOS-3 or GEOS-4 only):
      !
      !    P(I,J,L) = Ap(L) + ( Bp(L) * Psurface(I,J) )
      !
      ! Also here Ap, Bp will be flipped since both TPCORE versions
      ! index levels from the atm. top downwards (bdf, bmy, 5/8/03)
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

      !-----------------------------------------------------------------
      ! Prior to 6/24/03:
      !ALLOCATE( DSIG( LLPAR ), STAT=AS )
      !IF ( AS /= 0 ) CALL ALLOC_ERR( 'DSIG' )
      !-----------------------------------------------------------------

      ! Surface area [m2]
      DO J = 1, JJPAR
         A_M2(J) = GET_AREA_M2( J )
      ENDDO

      !-----------------------------------------------------------------
      ! Prior to 6/24/03:
      !! Layer thickness [unitless]
      !DO L = 1, LLPAR
      !   DSIG(L) = GET_BP( L ) - GET_BP( L+1 )
      !ENDDO
      !-----------------------------------------------------------------

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
      !-----------------------------------------------------------------
      ! Prior to 6/24/03:
      !IF ( ALLOCATED( DSIG   ) ) DEALLOCATE( DSIG   )
      !-----------------------------------------------------------------
      IF ( ALLOCATED( STT_I1 ) ) DEALLOCATE( STT_I1 )
      IF ( ALLOCATED( STT_I2 ) ) DEALLOCATE( STT_I2 )
      IF ( ALLOCATED( STT_J1 ) ) DEALLOCATE( STT_J1 )
      IF ( ALLOCATED( STT_J2 ) ) DEALLOCATE( STT_J2 )

      ! Return to calling program
      END SUBROUTINE CLEANUP_TRANSPORT

!------------------------------------------------------------------------------

      END MODULE TRANSPORT_MOD
