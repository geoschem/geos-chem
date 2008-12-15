! $Id: tpcore_fvdas_mod.f90,v 1.9 2008/12/15 15:55:14 bmy Exp $
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: Tpcore_FvDas_Mod
!
! !DESCRIPTION: \subsection*{Overview}
!  Module Tpcore\_Fvdas\_Mod contains routines for the TPCORE
!  transport scheme, as implemented in the GMI model (cf. John Tannahill),
!  based on Lin \ Rood 1995.  The Harvard Atmospheric Chemistry Modeling Group 
!  has added modifications to implement the Philip-Cameron Smith pressure 
!  fixer for mass conservation.  Mass flux diagnostics have also been added.
!
!\subsection*{References}
!  \begin{enumerate}
!  \item Lin, S.-J., and R. B. Rood, 1996: \emph{Multidimensional flux 
!         form semi-Lagrangian transport schemes}, 
!         \underline{ Mon. Wea. Rev.}, \textbf{124}, 2046-2070.
!  \item Lin, S.-J., W. C. Chao, Y. C. Sud, and G. K. Walker, 1994: 
!         \emph{A class of the van Leer-type transport schemes and its 
!         applications to the moisture transport in a General Circulation 
!         Model}, \underline{Mon. Wea. Rev.}, \textbf{122}, 1575-1593.
!  \end{enumerate}
!
!\subsection*{Selecting E/W, N/S and vertical advection options}
!
!  The flags IORD, JORD, KORD select which transport schemes are used in the 
!  E/W, N/S, and vertical directions, respectively.  Here is a list of the
!  possible values that IORD, JORD, KORD may be set to (original notes from
!  S-J Lin):
!   
!  \begin{enumerate}
!  \item 1st order upstream scheme (too diffusive, not a real option; 
!         it can be used for debugging purposes; this is THE only known 
!         "linear" monotonic advection scheme.).
!  \item 2nd order van Leer (full monotonicity constraint; 
!         see Lin et al 1994, MWR)
!  \item monotonic PPM* (Collela \& Woodward 1984)
!  \item semi-monotonic PPM (same as 3, but overshoots are allowed)
!  \item positive-definite PPM (constraint on the subgrid distribution is
!         only strong enough to prevent generation of negative values;
!         both overshoots \& undershoots are possible).
!  \item un-constrained PPM (nearly diffusion free; faster but
!         positivity of the subgrid distribution is not quaranteed. Use
!         this option only when the fields and winds are very smooth.
!  \item Huynh/Van Leer/Lin full monotonicity constraint.  Only KORD can be 
!         set to 7 to enable the use of Huynh's 2nd monotonicity constraint 
!         for piece-wise parabolic distribution.
!  \end {enumerate}
!
!  Recommended values:
!
!  \begin{itemize}
!  \item IORD=JORD=3 for high horizontal resolution.
!  \item KORD=3 or 7  
!  \end{itemize}  
!
!  The implicit numerical diffusion decreases as \_ORD increases.
!  DO NOT use option 4 or 5 for non-positive definite scalars
!  (such as Ertel Potential Vorticity).
!\\
!\\
! In GEOS-Chem we have been using IORD=3, JORD=3, KORD=7.  We have tested
! the OpenMP parallelization with these options.  GEOS-Chem users who wish to
! use different (I,J,K)ORD options should consider doing single-procsessor
! vs. multi-processor tests to test the implementation of the parallelization.
!
!\subsection*{GEOS-4 and GEOS-5 Hybrid Grid Definition}
!  
!  For GEOS-4 and GEOS-5 met fields, the pressure at the bottom edge of 
!  grid box (I,J,L) is defined as follows:
!  
!     $$P_{edge}(I,J,L) = A_{k}(L) + [ B_{k}(L) * P_{surface}(I,J) ]$$
!  
!  where
!
!  \begin{itemize}
!  \item $P_{surface}$(I,J) is the "true" surface pressure at lon,lat (I,J)
!  \item $A_{k}$(L) has the same units as surface pressure [hPa]
!  \item $B_{k}$(L) is a unitless constant given at level edges
!  \end{itemize}
!
!  $A_{k}(L)$ and $B_{k}(L)$ are supplied to us by GMAO.
!\\
!\\
! !REMARKS:
!  Ak(L) and Bk(L) are defined at layer edges.
!
!
!                  /////////////////////////////////
!              / \ ------ Model top P=ak(1) --------- ak(1), bk(1)
!               |
!    delp(1)    |  ........... q(i,j,1) ............
!               |
!              \ / ---------------------------------  ak(2), bk(2)
!
!
!
!              / \ ---------------------------------  ak(k), bk(k)
!               |
!    delp(k)    |  ........... q(i,j,k) ............
!               |
!              \ / ---------------------------------  ak(k+1), bk(k+1)
!
!
!
!              / \ ---------------------------------  ak(km), bk(km)
!               |
!    delp(km)   |  ........... q(i,j,km) .........
!               |
!              \ / -----Earth's surface P=Psfc ------ ak(km+1), bk(km+1)
!                 //////////////////////////////////
!
! Note: surface pressure can be of any unit (e.g., pascal or mb) as 
! long as it is consistent with the definition of (ak, bk) defined above.
!
! Winds (u,v), ps, and q are assumed to be defined at the same points.
! 
! The latitudes are given to the initialization routine: init_tpcore.
!
! !INTERFACE: 
!
MODULE Tpcore_FvDas_Mod
!
! !USES:
! 
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  ::  Init_Tpcore
  PUBLIC  ::  Exit_Tpcore
  PUBLIC  ::  Tpcore_FvDas
!
! !PRIVATE MEMBER FUNCTIONS:
! 
  PRIVATE ::  Average_Const_Poles
  PRIVATE ::  Set_Cross_Terms
  PRIVATE ::  Calc_Vert_Mass_Flux
  PRIVATE ::  Set_Jn_Js
  PRIVATE ::  Calc_Advec_Cross_Terms
  PRIVATE ::  Qckxyz
  PRIVATE ::  Set_Lmts
  PRIVATE ::  Set_Press_Terms
  PRIVATE ::  Calc_Courant
  PRIVATE ::  Calc_Divergence
  PRIVATE ::  Do_Divergence_Pole_Sum
  PRIVATE ::  Do_Cross_Terms_Pole_I2d2
  PRIVATE ::  Xadv_Dao2
  PRIVATE ::  Yadv_Dao2
  PRIVATE ::  Do_Yadv_Pole_I2d2
  PRIVATE ::  Do_Yadv_Pole_Sum
  PRIVATE ::  Xtp
  PRIVATE ::  Xmist
  PRIVATE ::  Fxppm
  PRIVATE ::  Lmtppm
  PRIVATE ::  Ytp
  PRIVATE ::  Ymist
  PRIVATE ::  Do_Ymist_Pole1_I2d2
  PRIVATE ::  Do_Ymist_Pole2_I2d2
  PRIVATE ::  Fyppm
  PRIVATE ::  Do_Fyppm_Pole_I2d2
  PRIVATE ::  Do_Ytp_Pole_Sum
  PRIVATE ::  Fzppm
  PRIVATE ::  Average_Press_Poles
!
! !PRIVATE DATA MEMBERS:
!  
  REAL*8, ALLOCATABLE, SAVE :: dtdx5(:)
  REAL*8, ALLOCATABLE, SAVE :: dtdy5(:)
  REAL*8, ALLOCATABLE, SAVE :: cosp(:)
  REAL*8, ALLOCATABLE, SAVE :: cose(:)
  REAL*8, ALLOCATABLE, SAVE :: gw(:)
  REAL*8, ALLOCATABLE, SAVE :: rgw(:)
  REAL*8, ALLOCATABLE, SAVE :: DLAT(:)
  REAL*8, ALLOCATABLE, SAVE :: RGW_25(:) 
!
! !AUTHOR:
! Original code from Shian-Jiann Lin, GMAO 
! Modified for GMI model by John Tannahill, LLNL (jrt@llnl.gov)
! Implemented into GEOS-Chem by Claire Carouge (ccarouge@seas.harvard.edu)
! ProTeX documentation added by Bob Yantosca (yantosca@seas.harvard.edu)
! OpenMP parallelization added by Bob Yantosca (yantosca@seas.harvard.edu)
!
! !REVISION HISTORY:
! 05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                             Yeh with the TPCORE routines from the GMI model.
!                             This eliminates the polar overshoot in the
!                             stratosphere.
! 05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                             Declare all REAL variables as REAL*8.  Added
!                             OpenMP parallel loops in various routines (and
!                             made some modifications to facilitate OpenMP).
!EOP
!------------------------------------------------------------------------------

CONTAINS

!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Init_Tpcore
!
! !DESCRIPTION: Subroutine Init\_Tpcore allocates and initializes all module
!  variables,
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Init_Tpcore( IM, JM, KM, JFIRST, JLAST, NG, MG, dt, ae, clat )
!
! !USES:
!
#   include "CMN_GCTM"                ! Physical constants etc.
!
! !INPUT PARAMETERS: 
!
    INTEGER, INTENT(IN)  :: IM        ! Global E-W dimension
    INTEGER, INTENT(IN)  :: JM        ! Global N-S dimension
    INTEGER, INTENT(IN)  :: KM        ! Vertical dimension
    INTEGER, INTENT(IN)  :: NG        ! large ghost width
    INTEGER, INTENT(IN)  :: MG        ! small ghost width
    REAL*8,  INTENT(IN)  :: dt        ! Time step in seconds
    REAL*8,  INTENT(IN)  :: ae        ! Earth's radius (m)
    REAL*8,  INTENT(IN)  :: clat(JM)  ! latitude in radian
!
! !OUTPUT PARAMETERS:
!
    INTEGER, INTENT(OUT) :: JFIRST    ! Local first index for N-S direction
    INTEGER, INTENT(OUT) :: JLAST     ! Local last  index for N-S direction
!
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    REAL*8  :: elat(jm+1)      ! cell edge latitude in radian
    REAL*8  :: sine(jm+1)  
    REAL*8  :: SINE_25(JM+1)   ! 
    REAL*8  :: dlon
    !----------------------------------------  
    ! Prior to 12/12/08:
    ! Use PI from CMN_GCTM (bmy, 12/12/08)
    !REAL*8  :: pi
    !----------------------------------------
    INTEGER :: I, J

    ! NOTE: since we are not using MPI parallelization, we can set JFIRST 
    ! and JLAST to the global grid limits in latitude. (bmy, 12/3/08)
    jfirst = 1
    jlast  = jm
    
    if ( jlast - jfirst < 2 ) then
       write(*,*) 'Minimum size of subdomain is 3'
    endif
    
    !----------------
    ! Allocate arrays
    !----------------
    
    ALLOCATE( cosp  ( JM ) ) 
    ALLOCATE( cose  ( JM ) ) 
    ALLOCATE( gw    ( JM ) ) 
    ALLOCATE( rgw   ( JM ) ) 
    ALLOCATE( dtdx5 ( JM ) ) 
    ALLOCATE( dtdy5 ( JM ) ) 
    ALLOCATE( DLAT  ( JM ) )    ! For PJC pressure-fixer 
    ALLOCATE( RGW_25( JM ) )    ! Latitude factor, for PJC pressure-fixer
    
    !----------------------------------------
    ! Prior to 12/12/08:
    ! Use PI from CMN_GCTM (bmy, 12/12/08)
    !PI   = 4.0d0 * ATAN(1.0d0)
    !----------------------------------------
    
    !----------------------------------------
    ! Prior to 12/12/08:
    ! Use double precision (bmy, 12/12/08)
    !dlon = 2.d0 * PI / float(im)
    !----------------------------------------
    dlon = 2.d0 * PI / DBLE( IM )
    
    ! S. Pole
    elat(1)    = -0.5d0*PI
    sine(1)    = -1.0d0
    SINE_25(1) = -1.0d0
    cose(1)    =  0.0d0
    
    do j=2,jm
       elat(j)    = 0.5d0*(clat(j-1) + clat(j))
       sine(j)    = SIN( elat(j) )
       SINE_25(J) = SIN( CLAT(J) )
       cose(j)    = COS( elat(j) ) 
    enddo
    
    ! N. Pole
    elat(jm+1)    = 0.5d0*PI       
    sine(jm+1)    = 1.0d0
    SINE_25(JM+1) = 1.0d0
    
    ! Polar cap (S. Pole)
    dlat(1) = 2.d0*(elat(2) - elat(1))  
    do j=2,jm-1
       dlat(j) = elat(j+1) - elat(j)
    enddo

    ! Polar cap (N. Pole)
    dlat(jm) = 2.0d0*(elat(jm+1) - elat(jm))    
    
    do j=1,jm
       gw(j)     = sine(j+1) - sine(j)
       cosp(j)   = gw(j) / dlat(j)
       rgw(j)    = 1.0d0 / gw(j)
       RGW_25(J) = 1.0d0 / ( SINE_25(J+1) - SINE_25(J) )

       dtdx5(j)  = 0.5d0 * dt / (dlon*ae*cosp(j))
       dtdy5(j)  = 0.5d0 * dt / (ae*dlat(j))
    enddo
      
    ! Echo info to stdout
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    WRITE( 6, '(a)' ) &
 'TPCORE_FVDAS (based on GMI) Tracer Transport Module successfully initialized'
    WRITE( 6, '(a)' ) REPEAT( '=', 79 )
    
  END SUBROUTINE Init_Tpcore
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Exit_Tpcore
!
! !DESCRIPTION: Subroutine Exit\_Tpcore deallocates all module variables.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Exit_Tpcore
!
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.
!EOP
!------------------------------------------------------------------------------
!BOC
    
    ! Deallocate arrays only if they are allocated
    IF ( ALLOCATED( COSP   ) ) DEALLOCATE( COSP   ) 
    IF ( ALLOCATED( COSE   ) ) DEALLOCATE( COSE   ) 
    IF ( ALLOCATED( GW     ) ) DEALLOCATE( GW     ) 
    IF ( ALLOCATED( RGW    ) ) DEALLOCATE( RGW    ) 
    IF ( ALLOCATED( RGW_25 ) ) DEALLOCATE( RGW_25 )  ! (bdf, bmy, 10/29/04)
    IF ( ALLOCATED( DTDX5  ) ) DEALLOCATE( DTDX5  ) 
    IF ( ALLOCATED( DTDY5  ) ) DEALLOCATE( DTDY5  ) 
    IF ( ALLOCATED( DLAT   ) ) DEALLOCATE( DLAT   )

  END SUBROUTINE Exit_Tpcore
!EOC  
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Tpcore_FvDas
!
! !DESCRIPTION: Subroutine Tpcore\_FvDas takes horizontal winds on sigma 
!  (or hybrid sigma-p) surfaces and calculates mass fluxes, and then updates 
!   the 3D mixing ratio fields one time step (tdt).  The basic scheme is a 
!   Multi-Dimensional Flux Form Semi-Lagrangian (FFSL) based on the van Leer 
!   or PPM (see Lin and Rood, 1995).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Tpcore_FvDas( dt,       ae,       IM,      JM,     KM,       &
                           JFIRST,   JLAST,    ng,      mg,     nq,       &
                           ak,       bk,       u,       v,      ps1,      &
                           ps2,      ps,       q,       iord,   jord,     &
                           kord,     n_adj,    XMASS,   YMASS,  MASSFLEW, &
                           MASSFLNS, MASSFLUP, AREA_M2, TCVV,   ND24,     &
                           ND25,     ND26 )
!
! !USES:
!
    ! Include file w/ physical constants
#   include "CMN_GCTM"
!
! !INPUT PARAMETERS: 
!
    ! Transport time step [s]
    REAL*8,  INTENT(IN)    :: dt                    

    ! Earth's radius [m]
    REAL*8,  INTENT(IN)    :: ae                    

    ! Global E-W, N-S, and vertical dimensions
    INTEGER, INTENT(IN)    :: IM         
    INTEGER, INTENT(IN)    :: JM         
    INTEGER, INTENT(IN)    :: KM         

    ! Latitude indices for local first box and local last box
    ! (NOTE: for global grids these are 1 and JM, respectively)
    INTEGER, INTENT(IN)    :: JFIRST     
    INTEGER, INTENT(IN)    :: JLAST 

    ! Primary ghost region
    ! (NOTE: only required for MPI parallelization; use 0 otherwise)
    INTEGER, INTENT(IN)    :: ng
         
    ! Secondary ghost region
    ! (NOTE: only required for MPI parallelization; use 0 otherwise)
    INTEGER, INTENT(IN)    :: mg         

    ! Ghosted latitudes (3 required by PPM)
    ! (NOTE: only required for MPI parallelization; use 0 otherwise)
    INTEGER, INTENT(IN)    :: nq         

    ! Flags to denote E-W, N-S, and vertical transport schemes
    INTEGER, INTENT(IN)    :: iord       
    INTEGER, INTENT(IN)    :: jord       
    INTEGER, INTENT(IN)    :: kord       

    ! Number of adjustments to air_mass_flux (0 = no adjustment)
    INTEGER, INTENT(IN)    :: n_adj      
    
    ! Ak and Bk coordinates to specify the hybrid grid
    ! (see the REMARKS section below)
    REAL*8,  INTENT(IN)    :: ak(KM+1)              
    REAL*8,  INTENT(IN)    :: bk(KM+1)              

    ! u-wind (m/s) at mid-time-level (t=t+dt/2)
    REAL*8,  INTENT(IN)    :: u(IM,JFIRST:JLAST,KM) 

    ! E/W and N/S mass fluxes [kg/s]
    ! (These are computed by the pressure fixer, and passed into TPCORE)
    REAL*8,  INTENT(IN)    :: XMASS(IM,JM,KM)
    REAL*8,  INTENT(IN)    :: YMASS(IM,JM,KM)

    ! Grid box surface area for mass flux diag [m2]
    REAL*8,  INTENT(IN)    :: AREA_M2(JM)        

    ! Tracer masses for flux diag
    REAL*8,  INTENT(IN)    :: TCVV(NQ)              

    ! Diagnostic flags
    INTEGER, INTENT(IN)    :: ND24    ! Turns on E/W     flux diagnostic
    INTEGER, INTENT(IN)    :: ND25    ! Turns on N/S     flux diagnostic
    INTEGER, INTENT(IN)    :: ND26    ! Turns on up/down flux diagnostic 
!
! !INPUT/OUTPUT PARAMETERS: 
!
    ! V-wind (m/s) at mid-time-level (t=t+dt/2)
    REAL*8,  INTENT(INOUT) :: v(IM, JFIRST-MG:JLAST+MG, KM) 

    ! surface pressure at current time
    REAL*8,  INTENT(INOUT) :: ps1(IM, JFIRST:JLAST)  

    ! surface pressure at future time=t+dt
    REAL*8,  INTENT(INOUT) :: ps2(IM, JFIRST:JLAST)  

    ! Tracer "mixing ratios" [v/v]
    REAL*8,  INTENT(INOUT) :: q(IM, JFIRST-NG:JLAST+NG, KM, NQ)  

    ! E/W, N/S, and up/down diagnostic mass fluxes
    REAL*8,  INTENT(INOUT) :: MASSFLEW(IM,JM,KM,NQ)  ! for ND24 diagnostic
    REAL*8,  INTENT(INOUT) :: MASSFLNS(IM,JM,KM,NQ)  ! for ND25 diagnostic
    REAL*8,  INTENT(INOUT) :: MASSFLUP(IM,JM,KM,NQ)  ! for ND26 diagnostic 
!
! !OUTPUT PARAMETERS:
!
    ! "Predicted" surface pressure [hPa]
    REAL*8,  INTENT(OUT)   :: ps(IM,JFIRST:JLAST)  
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO) 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    LOGICAL, PARAMETER :: FILL = .true.                 ! Fill negatives ?
    INTEGER, PARAMETER :: ADVEC_CONSRV_OPT = 2          ! 2=floating pressure
    LOGICAL, PARAMETER :: CROSS = .true.
!
! !LOCAL VARIABLES:
!
    INTEGER            :: rj2m1
    INTEGER            :: j1p, j2p
    INTEGER            :: jn (km)
    INTEGER            :: js (km)
    INTEGER            :: il, ij, ik, iq, k, j, i
    INTEGER            :: num, k2m1
                       
    REAL*8             :: dap   (km)
    REAL*8             :: dbk   (km)
    REAL*8             :: cx(im,jfirst-ng:jlast+ng,km)  ! E-W CFL # on C-grid
    REAL*8             :: cy(im,jfirst:jlast+mg,km)     ! N-S CFL # on C-grid
    REAL*8             :: delp1(im, jm, km)  
    REAL*8             :: delp2(im, jm, km)  
    REAL*8             :: delpm(im, jm, km)
    REAL*8             :: pu   (im, jm, km)
    REAL*8             :: dpi(im, jm, km)
    REAL*8             :: geofac  (jm)
    REAL*8             :: geofac_pc
    REAL*8             :: dp
    REAL*8             :: dps_ctm(im,jm)
    REAL*8             :: ua (im, jm, km)
    REAL*8             :: va (im, jm, km)
    REAL*8             :: wz(im, jm, km)
    REAL*8             :: dq1(im,jfirst-ng:jlast+ng,km)  
    REAL*8             :: qqu(im, jm, km)
    REAL*8             :: qqv(im, jm, km)
    REAL*8             :: adx(im, jm, km)
    REAL*8             :: ady(im, jm, km)
    REAL*8             :: fx (im, jm, km)
    REAL*8             :: fy (im, jm+1, km)           ! one more for edges
    REAL*8             :: fz  (im, jm, km)
    REAL*8             :: qtemp (im, jm, km)
    REAL*8             :: DTC(IM,JM,KM)               ! up/down flux temp array
    REAL*8             :: TRACE_DIFF                  ! up/down flux variable
                       
    LOGICAL, SAVE      :: first = .true.
    
    !     ----------------------------------------------------
    !     ilmt : controls various options in E-W     advection
    !     jlmt : controls various options in N-S     advection
    !     klmt : controls various options in vertcal advection
    !     ----------------------------------------------------
    
    INTEGER, SAVE      :: ilmt, jlmt, klmt
    INTEGER            :: js2g0, jn2g0
    
    !     ----------------
    !     Begin execution.
    !     ----------------
    
    ! Add definition of j1p and j2p for enlarge polar cap. (ccc, 11/20/08)
    j1p = 3
    j2p = jm - j1p + 1
    
    ! Average surf. pressures in the polar cap. (ccc, 11/20/08)
    CALL Average_Press_Poles( area_m2, ps1, 1, im, 1, jm, 1, im, 1, jm )
    CALL Average_Press_Poles( area_m2, ps2, 1, im, 1, jm, 1, im, 1, jm )
    
    
    ! Calculation of some geographic factors. (ccc, 11/20/08)
    rj2m1 = jm - 1
    dp    = PI / rj2m1
    
    do ij = 1, jm
       geofac(ij) = dp / (2.0d0 * area_m2(ij)/(sum(area_m2) * im) * im)
    end do
    
    geofac_pc =  &
         dp / (2.0d0 * (Sum (area_m2(1:2))/(sum(area_m2) * im)) * im)
    
    
    if (first) then
       
       first = .false.
       
     ! =============
       call Set_Lmts  &
     ! =============
            (ilmt, jlmt, klmt, im, jm, iord, jord, kord)
       
    end if
    
    ! Pressure calculations. (ccc, 11/20/08)
    do ik=1,km
       dap(ik) = ak(ik+1) - ak(ik)
       dbk(ik) = bk(ik+1) - bk(ik)
    enddo
      
  ! ====================
    call Set_Press_Terms  &
  ! ====================
         (dap, dbk, ps1, ps2, delp1, delpm, pu, &
         1, jm, 1, im, 1, jm, &
         j1p, j2p, 1, im, 1, jm, 1, km)
    !
    !...intent(in)  dap - difference in ai across layer (mb)
    !...intent(in)  dbk - difference in bi across layer (mb)
    !...intent(in)  pres1 - surface pressure at t1 (mb)
    !...intent(in)  pres2 - surface pressure at t1+tdt (mb)
    !...intent(out) delp1 - pressure thickness at t1 (mb)
    !...intent(out) delpm - pressure thickness at t1+tdt/2 (mb)
    !...intent(out) pu - pressure at edges of box for "u" (mb)
    !
    
    if (j1p /= 1+1) then
       
       do iq = 1, nq
          
       !  ========================
          call Average_Const_Poles  &
       !  ========================
               (dap, dbk, area_m2, ps1, q(:,:,:,iq), &
               1, jm, im, &
               1, im, 1, jm, 1, km, 1, im, 1, jm)
          
       end do
       
    end if
    
    
  ! =================
    call Calc_Courant  &
  ! =================
         (cose, delpm, pu, xmass, ymass, cx, cy, &
         j1p, j2p, &
         1, jm, 1, im, 1, jm, 1, im, 1, jm, 1, km)
    
  ! ====================
    call Calc_Divergence  &
  ! ====================
         (.true., geofac_pc, geofac, dpi, xmass, ymass, &
         j1p, j2p, 1, im, &
         1, jm, 1, im, 1, jm, 1, im, 1, jm, 1, km)
    
    dps_ctm(:,:) = Sum (dpi(:,:,:), dim=3)
    
  ! ====================
    call Set_Cross_Terms  &
  ! ====================
         (cx, cy, ua, va, &
         j1p, j2p, 1, im, 1, jm, &
         1, im, 1, jm, 1, im, 1, jm, 1, km, CROSS)
    
  ! ========================
    call Calc_Vert_Mass_Flux  &
  ! ========================
         (dbk, dps_ctm, dpi, wz, &
         1, im, 1, jm, 1, km)
    
    !.sds2.. have all mass flux here: east-west(xmass), 
    !        north-south(ymass), vertical(wz)
    !.sds2.. save omega (vertical flux) as diagnostic
      
    
  ! ==============
    call Set_Jn_Js  &
  ! ==============
           (jn, js, cx, &
           1, im, 1, jm, 1, jm, j1p, j2p, &
           1, im, 1, jm, 1, km)
    
    
    !.sds.. convert to "mass"
    do iq = 1, nq
       dq1(:,:,:) = q(:,:,:,iq) * delp1(:,:,:)
       
     ! ===========================
       call Calc_Advec_Cross_Terms  &
     ! ===========================
            (jn, js, q(:,:,:,iq), qqu, qqv, ua, va, &
            j1p, j2p, im, 1, jm, 1, im, 1, jm, &
            1, im, 1, jm, 1, km, CROSS)

       !.sds.. notes on arrays
       !  q (in)    - species mixing ratio
       !  qqu (out) - concentration contribution from E-W 
       !             advection cross terms(mixing ratio)
       !  qqv (out) - concentration contribution from N-S 
       !             advection cross terms(mixing ratio)
       !  ua  (in)  - average of Courant numbers from il and il+1
       !  va  (in)  - average of Courant numbers from ij and ij+1
              
     ! ----------------------------------------------------
     !  Add advective form E-W operator for E-W cross terms.
     ! ----------------------------------------------------
         
     ! ==============
       call Xadv_Dao2  &
     ! ==============
            (2, jn, js, adx, qqv, ua, &
            1, im, 1, jm, 1, jm, j1p, j2p, &
            1, im, 1, jm, 1, km)
       !.sds notes on output arrays
       !  adx (out)- cross term due to E-W advection (mixing ratio)
       !  qqv (in) - concentration contribution from N-S 
       !             advection (mixing ratio)
       !  ua  (in) - average of Courant numbers from il and il+1
       !.sds
         
     ! ----------------------------------------------------
     ! Add advective form N-S operator for N-S cross terms.
     ! ----------------------------------------------------
         
     ! ==============
       call Yadv_Dao2  &
     ! ==============
            (2, ady, qqu, va, &
            1, im, 1, jm, &
            j1p, j2p, 1, im, 1, jm, 1, im, 1, jm, 1, km)

       !.sds notes on output arrays
       !  ady (out)- cross term due to N-S advection (mixing ratio)
       !  qqu (in) - concentration contribution from N-S advection 
       !             (mixing ratio)
       !  va  (in) - average of Courant numbers from il and il+1
       !.sds
       !
       !.bmy notes: use a polar cap of 2 boxes (i.e. the "2" as
       ! the first argument to YADV_DAO2.  The older TPCORE only had
       ! a polar cap of 1 box (just the Pole itself).  Claire figured
       ! this out.  (bmy, 12/11/08)
       
       !... update constituent array qq1 by adding in cross terms 
       !           - use in fzppm
       q(:,:,:,iq) = q(:,:,:,iq) + ady + adx
       
       
     ! ========
       call Xtp  &
     ! ========
            (ilmt, jn, js, pu, cx, dq1, qqv, xmass, fx, &
            j1p, j2p, im, 1, jm, 1, im, 1, jm, &
            1, im, 1, jm, 1, km, IORD)
       
       !.sds notes on output arrays
       !  pu  (in)    - pressure at edges in "u" (mb)
       !  crx (in)    - Courant number in E-W direction
       !  dq1 (inout) - species density (mb) - updated with the E-W flux 
       !                fx in Xtp)
       !  qqv (inout) - concentration contribution from N-S advection 
       !                (mixing ratio)
       !  xmass(in)   - horizontal mass flux in E-W direction (mb)
       !  fx  (out)   - species E-W mass flux
       !.sds
         
     ! ========
       call Ytp  &
     ! ========
            (jlmt, geofac_pc, geofac, cy, dq1, qqu, qqv, ymass, fy, &
            j1p, j2p, 1, im, 1, jm, im, &
            1, im, 1, jm, 1, im, 1, jm, 1, km, jord)
         
       !.sds notes on output arrays
       !  cy (in)     - Courant number in N-S direction
       !  dq1 (inout) - species density (mb) - updated with the N-S flux 
       !                (fy in Ytp)
       !  qqu (in)    - concentration contribution from E-W advection 
       !                (mixing ratio)
       !  qqv (inout) - concentration contribution from N-S advection 
       !                (mixing ratio)  
       !  ymass(in)   - horizontal mass flux in E-W direction (mb)
       !  fy  (out)   - species N-S mass flux (need to mult by geofac)
       !.sds
         
       qtemp(:,:,:) = q(:,:,:,iq)
     ! ==========
       call Fzppm  &
     ! ==========
            (klmt, delp1, dpi, wz, dq1, q(:,:,:,iq), fz, &
            j1p, 1, jm, 1, im, 1, jm, &
            im, km, 1, im, 1, jm, 1, km)
         
       !.sds notes on output arrays
       !   wz  (in) : vertical mass flux
       !   dq1 (inout) : species density (mb)
       !   q (in) : species concentration (mixing ratio)
       !.sds
         
         
       if (FILL) then
        ! ===========
          call Qckxyz  &
        ! ===========
               (dq1, &
               j1p, j2p, 1, jm, &
               1, im, 1, jm, 1, im, 1, jm, 1, km)
       end if
       
       if (advec_consrv_opt == 0) then
          
          !----------------------------------------------------------------
          ! Prior to 12/5/08:
          ! Replace these with explicit DO loops to facilitate
          ! OpenMP parallelization (bmy, 12/5/08)
          !do ik = 1, km
          !   
          !   delp2(:,:,ik) =  &
          !        dap(ik) +  &
          !        (dbk(ik) * (ps1(:,:) +  &
          !        dps_ctm(:,:)))
          !     
          !end do
          !----------------------------------------------------------------
          
          !$OMP PARALLEL DO           &
          !$OMP DEFAULT( SHARED     ) &
          !$OMP PRIVATE( IK, IJ, IL )
          do ik = 1, km
          do ij = 1, jm
          do il = 1, im
             delp2(il,ij,ik) =  &
                  dap(ik) +  &
                  (dbk(ik) * (ps1(il,ij) +  &
                  dps_ctm(il,ij)))
               
          end do
          end do
          end do
          !$OMP END PARALLEL DO

       else if ((advec_consrv_opt == 1) .or.  &
            (advec_consrv_opt == 2)) then
          
          !----------------------------------------------------------------
          ! Prior to 12/5/08:
          ! Replace these with explicit DO loops to facilitate
          ! OpenMP parallelization (bmy, 12/5/08)
          !do il = 1, im
          !   
          !   delp2(:,:,ik) =  &
          !        dap(ik) +  &
          !        (dbk(ik) * ps2(:,:))
          !   
          !end do
          !----------------------------------------------------------------

          !$OMP PARALLEL DO           &
          !$OMP DEFAULT( SHARED     ) &
          !$OMP PRIVATE( IK, IJ, IL )
          do ik = 1, km
          do ij = 1, jm
          do il = 1, im
             
             delp2(il,ij,ik) =  &
                  dap(ik) +  &
                  (dbk(ik) * ps2(il,ij))
             
          end do
          end do
          end do
          !$OMP END PARALLEL DO
            
       end if
       
       ! Calculate surf. pressure at t+dt. (ccc, 11/20/08)
       ps = ak(1)+sum(delp2,dim=3)
         
         
       q(:,:,:,iq) =  &
            dq1 / delp2
         
         
       if (j1p /= 2) then
            
          q(:,2,:,iq) = q(:,1,:,iq)
          q(:,jm-1,:,iq)  = q(:,jm,:,iq)
            
       end if
       
       ! Calculate fluxes for diag. (ccc, 11/20/08)
       !--------------------------------------------------------------
       ! Prior to 12/11/08:
       ! Set with J1P and J2P for extended polar cap (bmy, 12/11/08)
       !js2g0  = max(2,jfirst)          !  No ghosting
       !jn2g0  = min(jm-1,jlast)        !  No ghosting
       !--------------------------------------------------------------
       JS2G0  = MAX( J1P, JFIRST )     !  No ghosting
       JN2G0  = MIN( J2P, JLAST  )     !  No ghosting

       !======================================================================
       ! MODIFICATION by Harvard Atmospheric Chemistry Modeling Group
       !  
       ! Implement ND24 diag: E/W flux of tracer [kg/s]  (ccarouge 12/2/08)  
       !
       !  The unit conversion is:
       !
       !  Mass    P diff     100      1       area of     kg tracer     1
       ! ------ = in grid *  ---  *  ---   *  grid box * ----------- * ---
       !  time    box         1       g       AREA_M2      kg air       s
       !
       !   kg       hPa       Pa     s^2        m^2          1          1 
       !  ----  =  -----  * ----- * -----  *   -----   *   ------  * --------
       !   s         1       hPa      m          1          TCVV      DeltaT
       !======================================================================
       IF ( ND24 > 0 ) THEN

          ! Zero temp array
          DTC = 0d0

          !$OMP PARALLEL DO        &
          !$OMP DEFAULT( SHARED  ) &
          !$OMP PRIVATE( I, J, K ) 
          DO K = 1,     KM
          DO J = JS2G0, JN2G0
          DO I = 1,     IM

             ! Compute mass flux
             DTC(I,J,K)         = ( FX(I,J,K)  * AREA_M2(J)  * 100.d0 ) / &
                                  ( TCVV(IQ)   * DT          * 9.8d0  )

             ! Save into MASSFLEW diagnostic array
             MASSFLEW(I,J,K,IQ) = MASSFLEW(I,J,K,IQ) + DTC(I,J,K)

          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ENDIF

       !======================================================================
       ! MODIFICATION by Harvard Atmospheric Chemistry Modeling Group
       !  
       ! Implement ND25 diag: N/S flux of tracer [kg/s] 
       ! (bdf, bmy, 9/28/04, ccarouge 12/12/08)
       !
       ! NOTE, the unit conversion is the same as desciribed above for the
       ! ND24 E-W diagnostics.  However, we also need to multiply by the 
       ! RGW_25 area factor.
       !======================================================================
       IF ( ND25 > 0 ) THEN

          ! Zero temp array
          DTC = 0d0

          !$OMP PARALLEL DO        &
          !$OMP DEFAULT( SHARED  ) &
          !$OMP PRIVATE( I, J, K ) 
          DO K = 1, KM
          DO J = 1, JM 
          DO I = 1, IM 

             ! Compute mass flux
             DTC(I,J,K)    = ( FY(I,J,K) * AREA_M2(J) * RGW_25(J) * 1d2 ) / & 
                             ( TCVV(IQ)  * DT         * 9.8d0           ) 

             ! Save into MASSFLNS diagnostic array
             MASSFLNS(I,J,K,IQ) = MASSFLNS(I,J,K,IQ) + DTC(I,J,K) 

          ENDDO
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO

       ENDIF

       !======================================================================
       ! MODIFICATION by Harvard Atmospheric Chemistry Modeling Group
       !  
       ! Implement ND26 diag: Up/down flux of tracer [kg/s] 
       ! (bmy, bdf, 9/28/04, ccarouge 12/2/08)
       !
       ! The vertical transport done in qmap.  We need to find the difference 
       ! in order to to interpret transport.
       !
       ! Break up diagnostic into up & down fluxes using the surface boundary 
       ! conditions.  Start from top down (really surface up for flipped 
       ! TPCORE)
       !
       ! By construction, MASSFLUP is flux into the bottom of the box. The 
       ! flux at the bottom of KM (the surface box) is not zero by design. 
       ! (phs, 3/4/08)
       !======================================================================
       IF ( ND26 > 0 ) THEN
          
          ! Zero temp array
          DTC = 0d0

          !-----------------
          ! start with top
          !-----------------
          K = 1
          
          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( I, J )  
          DO J  = 1, JM
          DO I  = 1, IM

             ! Compute mass flux
             DTC(I,J,K)         = ( Q(I,J,K,IQ)  * DELP1(I,J,K)   -   &
                                    QTEMP(I,J,K) * DELP2(I,J,K)  ) *  &
                                  (100d0) * AREA_M2(J) / ( 9.8d0 * TCVV(IQ) )
                
             ! top layer should have no residual.  the small residual is 
             ! from a non-pressure fixed flux diag.  The z direction may 
             ! be off by a few percent.
             !
             ! Uncomment now, since this is upflow into the box from its
             ! bottom (phs, 3/4/08)
             MASSFLUP(I,J,K,IQ) = MASSFLUP(I,J,K,IQ) + DTC(I,J,K)/DT
          ENDDO
          ENDDO
          !$OMP END PARALLEL DO
          
          !----------------------------------------------------
          ! Get the other fluxes using a mass balance equation
          !----------------------------------------------------
          DO K  = 2, KM

             !$OMP PARALLEL DO                 &
             !$OMP DEFAULT( SHARED )           &
             !$OMP PRIVATE( I, J, TRACE_DIFF )
             DO J  = 1, JM
             DO I  = 1, IM

                ! Compute tracer difference
                TRACE_DIFF         = ( Q(I,J,K,IQ)     * DELP1(I,J,K)  -  &
                                       QTEMP(I,J,K) * DELP2(I,J,K) ) *  &
                                       (100D0) * AREA_M2(J) /           &
                                       ( 9.8D0* TCVV(IQ) )
                
                ! Compute mass flux
                DTC(I,J,K)         = DTC(I,J,K-1) + TRACE_DIFF
                  
                ! Save to the MASSFLUP diagnostic array 
                MASSFLUP(I,J,K,IQ) = MASSFLUP(I,J,K,IQ) + DTC(I,J,K)/DT
                
             ENDDO
             ENDDO
             !$OMP END PARALLEL DO

          ENDDO
       ENDIF
    ENDDO
      
  END SUBROUTINE Tpcore_FvDas
!EOC  
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Average_Const_Poles
!
! !DESCRIPTION: Subroutine Average\_Const\_Poles averages the species 
!  concentrations at the Poles when the Polar cap is enlarged.  It makes the 
!  last two latitudes equal.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Average_Const_Poles( dap ,   dbk,   rel_area, pctm1, const1, &
                                  JU1_GL, J2_GL, I2_GL,    I1,    I2,     &
                                  JU1,    J2,    K1,       K2,    ILO,    &
                                  IHI,    JULO,  JHI )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices of the South Pole and North Pole
    INTEGER, INTENT(IN)   :: JU1_GL, J2_GL

    ! Global max longitude index
    INTEGER, INTENT(IN)   :: I2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)   :: I1,  I2
    INTEGER, INTENT(IN)   :: JU1, J2
    INTEGER, INTENT(IN)   :: K1,  K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)   :: ILO,  IHI
    INTEGER, INTENT(IN)   :: JULO, JHI

    ! Pressure difference across layer from (ai * pt) term [hPa]
    REAL*8,  INTENT(IN)   :: dap(K1:K2)

    ! Difference in bi across layer - the dSigma term
    REAL*8,  INTENT(IN)   :: dbk(K1:K2)

    ! Relative surface area of grid box [fraction]
    REAL*8,  INTENT(IN)   :: rel_area(JU1:J2)

    ! CTM surface pressure at t1 [hPa]
    REAL*8,  INTENT(IN)   :: pctm1( ILO:IHI, JULO:JHI )
!
! !INPUT/OUTPUT PARAMETERS: 
!
    ! Species concentration, known at zone center [mixing ratio]
    REAL*8, INTENT(INOUT) :: const1( I1:I2, JU1:J2, K1:K2 )
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO) 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.
!EOP
!------------------------------------------------------------------------------
!BOC  
!
! !LOCAL VARIABLES:
! 
    INTEGER :: ik, il
    
    REAL*8  :: meanq
    REAL*8  :: sum1, sum2

!   -----------------------------------------------------------------
!   delp1n : pressure thickness at North Pole, the psudo-density in a
!            hydrostatic system at t1 (mb)
!   delp1s : pressure thickness at South Pole, the psudo-density in a
!            hydrostatic system at t1 (mb)
!   -----------------------------------------------------------------

    REAL*8  :: delp1n(i1:i2, j2-1:j2,    k1:k2)
    REAL*8  :: delp1s(i1:i2,  ju1:ju1+1, k1:k2)


!   ----------------
!   Begin execution.
!   ----------------

!   =================
    if (ju1 == ju1_gl) then
!   =================

       !$OMP PARALLEL DO                         &
       !$OMP DEFAULT( SHARED )                   &
       !$OMP PRIVATE( IL, IK, SUM1, SUM2, MEANQ )
       do ik = k1, k2
          
          delp1s(i1:i2,ju1:ju1+1,ik) =                &
               dap(ik) +                             &
               (dbk(ik) * pctm1(i1:i2,ju1:ju1+1))
          
          sum1=0.0d0
          sum2=0.0d0
          do il = i1, i2
             sum1 = sum1 +                            &
                    Sum (const1  (il,ju1:ju1+1,ik) *  &
                         delp1s  (il,ju1:ju1+1,ik) *  &
                         rel_area(ju1:ju1+1))         &
                  / (sum(rel_area) * i2_gl)

             sum2 = sum2 +                           &
                     Sum (delp1s  (il,ju1:ju1+1,ik) * & 
                       rel_area(ju1:ju1+1))        &
                  / (sum(rel_area) * i2_gl) 
          enddo

          meanq = sum1 / sum2

          const1(:,ju1:ju1+1,ik) = meanq

       end do
       !$OMP END PARALLEL DO

    end if


!   ================
    if (j2 == j2_gl) then
!   ================

       !$OMP PARALLEL DO                         &
       !$OMP DEFAULT( SHARED )                   &
       !$OMP PRIVATE( IL, IK, SUM1, SUM2, MEANQ )
       do ik = k1, k2

          delp1n(i1:i2,j2-1:j2,ik) =               &
               dap(ik) +                           &
               (dbk(ik) * pctm1(i1:i2,j2-1:j2))
            
          sum1=0.0d0
          sum2=0.0d0
          do il = i1, i2
             sum1 = sum1 +                         & 
                    Sum (const1  (il,j2-1:j2,ik) * & 
                         delp1n  (il,j2-1:j2,ik) * & 
                         rel_area(j2-1:j2))        & 
                  / (sum(rel_area) * i2_gl)
            
             sum2 = sum2 +                         &
                    Sum (delp1n  (il,j2-1:j2,ik) * & 
                         rel_area(j2-1:j2))        &
                  / (sum(rel_area) * i2_gl)      
          enddo


          meanq = sum1 / sum2

          const1(:,j2-1:j2,ik) = meanq

       end do
       !$OMP END PARALLEL DO
         
    end if

  END SUBROUTINE Average_Const_Poles
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Cross_Terms
!
! !DESCRIPTION: Subroutine Set\_Cross\_Terms sets the cross terms for 
!  E-W horizontal advection.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Cross_Terms( crx,   cry,   ua, va, J1P,   J2P,   &
                              I1_GL, I2_GL, JU1_GL, J2_GL, ILO,   &
                              IHI,   JULO,  JHI,    I1,    I2,    &
                              JU1,   J2,    K1,     K2,    CROSS )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)   :: J1P,    J2P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)   :: I1_GL,  I2_GL
    INTEGER, INTENT(IN)   :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)   :: I1,     I2
    INTEGER, INTENT(IN)   :: JU1,    J2
    INTEGER, INTENT(IN)   :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)   :: ILO,    IHI
    INTEGER, INTENT(IN)   :: JULO,   JHI

    ! Courant number in E-W direction
    REAL*8,  INTENT(IN) :: crx(ILO:IHI, JULO:JHI, K1:K2)

    ! Courant number in N-S direction
    REAL*8,  INTENT(IN) :: cry(ILO:IHI, JULO:JHI, K1:K2)

    ! Logical switch.  If CROSS=T then cross-terms will be computed.
    LOGICAL, INTENT(IN) :: CROSS
!
! !OUTPUT PARAMETERS:
!
    ! Average of Courant numbers from il and il+1
    REAL*8, INTENT(OUT) :: ua(ILO:IHI, JULO:JHI, K1:K2)

    ! Average of Courant numbers from ij and ij+1
    REAL*8, INTENT(OUT) :: va(ILO:IHI, JULO:JHI, K1:K2)

! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO) 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!     

    ! Grid box indices for lon & lat
    INTEGER :: il, ij, ik
    
    !     ----------------
    !     Begin execution.
    !     ----------------
    
    
    if (.not. CROSS) then
       
       ua(:,:,:) = 0.0d0
       va(:,:,:) = 0.0d0
       
    else
       
       !---------------------------------------------------------------------
       ! Prior to 12/4/08:
       ! Rewrite DO loop to parallelize over vertical levels (bmy, 12/4/08)
       !do ij = j1p, j2p
       !   do il = i1, i2-1
       !      
       !      ua(il,ij,:) = 0.5d0 * (crx(il,ij,:) + crx(il+1,ij,:))
       !      
       !   end do
       !   ua(i2,ij,:) = 0.5d0 * (crx(i2,ij,:) + crx(1,ij,:))
       !
       !end do
       !
       !
       !do ij = j1p, j2p
       !   do il = i1, i2
       !      
       !      va(il,ij,:) = 0.5d0 * (cry(il,ij,:) + cry(il,ij+1,:))
       !      
       !   end do
       !end do
       !---------------------------------------------------------------------
       
       !$OMP PARALLEL DO            &
       !$OMP DEFAULT( SHARED )      &
       !$OMP PRIVATE( IL, IJ, IK )
       do ik = k1, k2

          do ij = j1p, j2p
             do il = i1, i2-1
             
                ua(il,ij,ik) = 0.5d0 * (crx(il,ij,ik) + crx(il+1,ij,ik))
             
             end do
             ua(i2,ij,ik) = 0.5d0 * (crx(i2,ij,ik) + crx(1,ij,ik))

          end do
       
       
          do ij = j1p, j2p
             do il = i1, i2
             
                va(il,ij,ik) = 0.5d0 * (cry(il,ij,ik) + cry(il,ij+1,ik))
             
             end do
          end do

       enddo
       !$OMP END PARALLEL DO
       
!      =============================
       call Do_Cross_Terms_Pole_I2d2  &
!      =============================
            (cry, va, &
             i1_gl, i2_gl, ju1_gl, j2_gl, j1p,  &
             ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)
       
       
    end if
    
  END SUBROUTINE Set_Cross_Terms
!EOC  
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Vert_Mass_Flux
!
! !DESCRIPTION: Subroutine Calc\_Vert\_Mass\_Flux calculates the vertical 
!  mass flux.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_Vert_Mass_Flux( dbk, dps_ctm, dpi, wz, I1,  &
                                  I2,  JU1,     J2,  K1, K2 )
!
! !INPUT PARAMETERS: 
!
    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)   :: I1,  I2
    INTEGER, INTENT(IN)   :: JU1, J2
    INTEGER, INTENT(IN)   :: K1,  K2

    ! Difference in bi across layer - the dSigma term
    REAL*8,  INTENT(IN)  :: dbk(K1:K2)

    ! CTM surface pressure tendency; sum over vertical of dpi
    ! calculated from original mass fluxes [hPa]
    REAL*8,  INTENT(IN)  :: dps_ctm(I1:I2, JU1:J2)

    ! Divergence at a grid point; used to calculate vertical motion [mb]
    REAL*8,  INTENT(IN)  :: dpi(I1:I2, JU1:J2, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Large scale mass flux (per time step tdt) in the vertical
    ! direction as diagnosed from the hydrostatic relationship [hPa]
    REAL*8, INTENT(OUT) :: wz(I1:I2, JU1:J2, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO) 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops
!EOP
!------------------------------------------------------------------------------
!BOC  
!
! !LOCAL VARIABLES:
!    
    INTEGER :: ik, ij, il
       
!   ----------------
!   Begin execution.
!   ----------------
     
!   --------------------------------------------------
!   Compute vertical mass flux from mass conservation.
!   --------------------------------------------------
    
    !---------------------------------------------------------------------
    ! Prior to 12/5/08:
    ! Need to add explicit IJ and IL loops for OpenMP parallelization
    ! (bmy, 12/5/08)
    ! 
    !wz(:,:,k1) =  &
    !     dpi(:,:,k1) -  &
    !     (dbk(k1) * dps_ctm(i1:i2,ju1:j2))
    !
    !wz(:,:,k2) = 0.0d0
    !
    !
    !do ik = k1 + 1, k2 - 1
    !   
    !   wz(:,:,ik) =  &
    !        wz (:,:,ik-1) +  &
    !        dpi(:,:,ik)   -  &
    !        (dbk(ik) * dps_ctm(i1:i2,ju1:j2))
    !
    !end do
    !---------------------------------------------------------------------

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( IJ, IL )
    do ij = ju1, j2
    do il = i1,  i2
       wz(il,ij,k1) =  &
            dpi(il,ij,k1) -  &
            (dbk(k1) * dps_ctm(il,ij))
    
       wz(il,ij,k2) = 0.0d0
    end do
    end do
    !$OMP END PARALLEL DO
    
    do ik = k1 + 1, k2 - 1
       
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IJ, IL )
       do ij = ju1, j2
       do il = i1,  i2

          wz(il,ij,ik) =  &
               wz (il,ij,ik-1) +  &
               dpi(il,ij,ik)   -  &
               (dbk(ik) * dps_ctm(il,ij))
       end do
       end do
       !$OMP END PARALLEL DO

    end do

    
  END SUBROUTINE Calc_Vert_Mass_Flux
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Jn_Js
!
! !DESCRIPTION: Subroutine Set\_Jn\_Js determines Jn and Js, by looking 
!  where Courant number is > 1.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Jn_Js( jn,  js,     crx,   ILO, IHI, JULO, &
                        JHI, JU1_GL, J2_GL, J1P, J2P, I1,   &
                        I2,  JU1,    J2,    K1,  K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)   :: J1P,    J2P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)   :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)   :: I1,     I2
    INTEGER, INTENT(IN)   :: JU1,    J2
    INTEGER, INTENT(IN)   :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)   :: ILO,    IHI
    INTEGER, INTENT(IN)   :: JULO,   JHI

    ! Courant number in E-W direction
    REAL*8,  INTENT(IN)  :: crx(ILO:IHI, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Northward of latitude index = jn; Courant numbers could be > 1,
    ! so use the flux-form semi-Lagrangian scheme
    INTEGER, INTENT(OUT) :: jn(K1:K2)

    ! Southward of latitude index = js; Courant numbers could be > 1,
    ! so use the flux-form semi-Lagrangian scheme
    INTEGER, INTENT(OUT) :: js(K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO) 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REMARKS:
!   We cannot parallelize this subroutine because there is a CYCLE statement
!   within the outer loop.
!
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  
!EOP
!------------------------------------------------------------------------------
!BOC  
!
! !LOCAL VARIABLES:
    
    INTEGER :: il, ij, ik
    INTEGER :: jn0, js0
    INTEGER :: jst, jend
    
    
    !     ----------------
    !     Begin execution.
    !     ----------------
    
    js0  = (j2_gl + 1 ) / 2
    jn0  = j2_gl - js0 + 1
    
    jst  = Max (ju1, j1p)
    jend = Min (j2,  js0)
    
    ikloop1: do ik = k1, k2
       
       js(ik) = j1p
 
       do ij = jend, jst, -1
          do il = i1, i2
             
             if (Abs (crx(il,ij,ik)) > 1.0d0) then
                
                js(ik) = ij
                
!               =============
                cycle ikloop1
!               =============
                
             end if
             
          end do
       end do
       
    end do ikloop1
    
    
    jst  = Max (ju1, jn0)
    jend = Min (j2,  j2p)
    
    ikloop2: do ik = k1, k2
       
       jn(ik) = j2p
       
       do ij = jst, jend
          do il = i1, i2
             
             if (Abs (crx(il,ij,ik)) > 1.0d0) then
                
                jn(ik) = ij
                
!               =============
                cycle ikloop2
!               =============

             end if
             
          end do
       end do

    end do ikloop2

  END SUBROUTINE Set_Jn_Js
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Advec_Cross_Terms
!
! !DESCRIPTION: Subroutine Calc\_Advec\_Cross\_Terms calculates the advective 
!  cross terms.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_Advec_Cross_Terms( jn,     js,    qq1,  qqu,  qqv,   &
                                     ua,     va,    J1P,  J2P,  I2_GL, &
                                     JU1_GL, J2_GL, ILO,  IHI,  JULO,  &
                                     JHI,    I1,    I2,   JU1,  J2,    &
                                     K1,     K2,    CROSS )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: J1P,    J2P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  ::         I2_GL
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI

    ! Northward of latitude index = jn, Courant numbers could be > 1,
    ! so use the flux-form semi-Lagrangian scheme
    INTEGER, INTENT(IN)  :: Jn(K1:K2)

    ! Southward of latitude index = js, Courant numbers could be > 1,
    ! so use the flux-form semi-Lagrangian scheme
    INTEGER, INTENT(IN)  :: Js(K1:K2)

    ! Species concentration (mixing ratio)
    REAL*8,  INTENT(IN)  :: qq1(ILO:IHI, JULO:JHI, K1:K2)

    ! Average of Courant numbers from il and il+1
    REAL*8,  INTENT(IN)  :: ua (ILO:IHI, JULO:JHI, K1:K2)

    ! Average of Courant numbers from ij and ij+1
    REAL*8,  INTENT(IN)  :: va (ILO:IHI, JULO:JHI, K1:K2)

    ! Logical switch: If CROSS=T then cross-terms are being computed
    LOGICAL, INTENT(IN)  :: CROSS
!
! !OUTPUT PARAMETERS:
!
    ! Concentration contribution from E-W advection [mixing ratio]
    REAL*8,  INTENT(OUT) :: qqu(ILO:IHI, JULO:JHI, K1:K2)

    ! concentration contribution from N-S advection [mixing ratio]
    REAL*8,  INTENT(OUT) :: qqv(ILO:IHI, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO) 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel do loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!    
    INTEGER :: i, imp, il, ij, ik, iu
    INTEGER :: jv, iuw, iue
    REAL*8  :: ril, rij, riu
    REAL*8  :: ru
    REAL*8  :: qtmp(-i2/3:i2+i2/3, julo:jhi, k1:k2)
    
    !     ----------------
    !     Begin execution.
    !     ----------------
    
    !----------------------------------------------------------------
    ! Prior to 12/5/08
    ! Now add explicit IJ and IK loops for OpenMP parallelization
    ! (bmy, 12/5/08)
    !do i = 1, i2
    !   qtmp(i,:,:) = qq1(i,:,:)
    !enddo
    !
    !do il = -i2/3, 0
    !   qtmp(il,:,:) = qq1(i2+il,:,:)
    !enddo
    !
    !do il = i2+1,i2+i2/3
    !   qtmp(il,:,:) = qq1(il-i2,:,:)
    !enddo
    !----------------------------------------------------------------

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( IK, IJ )
    do ik = k1,   k2
    do ij = julo, jhi
       do i = 1, i2
          qtmp(i,ij,ik) = qq1(i,ij,ik)
       enddo
    
       do il = -i2/3, 0
          qtmp(il,ij,ik) = qq1(i2+il,ij,ik)
       enddo
    
       do il = i2+1,i2+i2/3
          qtmp(il,ij,ik) = qq1(il-i2,ij,ik)
       enddo
    enddo
    enddo
    !$OMP END PARALLEL DO

!   ================
    if (.not. CROSS) then
!   ================
       
       qqu(:,:,:) = qq1(:,:,:)
       qqv(:,:,:) = qq1(:,:,:)
       
       
!   ====
    else
!   ====

       qqu(:,:,:) = 0.0d0
       qqv(:,:,:) = 0.0d0
       
       !$OMP PARALLEL DO                            &
       !$OMP DEFAULT( SHARED )                      &
       !$OMP PRIVATE( IL, IJ, IK, IU, RIU, RU, RIL )
       do ik = k1, k2
          do ij = j1p, j2p
             
             if ((ij <= js(ik)) .or. (ij >= jn(ik))) then

!             ----------------------------------------------------------
!             In Polar area, so need to deal with large courant numbers.
!             ----------------------------------------------------------

                do il = i1, i2
                   
!c?
                   iu  = ua(il,ij,ik)
                   riu = iu
                   ru  = ua(il,ij,ik) - riu
                   iu  = il - iu
                   
                   if (ua(il,ij,ik) >= 0.0d0) then
                      
                      qqu(il,ij,ik) =  &
                           qtmp(iu,ij,ik) +  &
                           ru * (qtmp(iu-1,ij,ik) - qtmp(iu,ij,ik))
                      
                   else
                      
                      qqu(il,ij,ik) =  &
                           qtmp(iu,ij,ik) +  &
                           ru * (qtmp(iu,ij,ik) - qtmp(iu+1,ij,ik))
                      
                   end if
                   
                   qqu(il,ij,ik) = qqu(il,ij,ik) - qtmp(il,ij,ik)
                   
                end do
                
             else  ! js < ij < jn
                
                !             ---------------------------
                !             Do interior area (use PPM).
                !             ---------------------------
                
                do il = i1, i2
                   
                   ril = il
                   iu  = ril - ua(il,ij,ik)
                   
                   qqu(il,ij,ik) =  &
                        ua(il,ij,ik) *  &
                        (qtmp(iu,ij,ik) - qtmp(iu+1,ij,ik))
                   
                end do
                
             end if

          end do
          
       end do
       !$OMP END PARALLEL DO

       !$OMP PARALLEL DO                   &
       !$OMP DEFAULT( SHARED )             &
       !$OMP PRIVATE( IL, IJ, IK, RIJ, JV )
       do ik = k1, k2
          do ij = j1p, j2p
             do il = i1, i2
                
!c?
                rij = ij
                jv  = rij - va(il,ij,ik)
                
                qqv(il,ij,ik) =  &
                     va(il,ij,ik) *  &
                     (qtmp(il,jv,ik) - qtmp(il,jv+1,ik))

             end do
          end do
       end do
       !$OMP END PARALLEL DO

       !----------------------------------------------------------------
       ! Prior to 12/5/08
       ! Now add explicit IJ and IK loops for OpenMP parallelization
       ! (bmy, 12/5/08)
       !qqu(i1:i2,ju1:j2,:) =  &
       !     qtmp(i1:i2,ju1:j2,:) + (0.5d0 * qqu(i1:i2,ju1:j2,:))
       !
       !qqv(i1:i2,ju1:j2,:) =  &
       !     qtmp(i1:i2,ju1:j2,:) + (0.5d0 * qqv(i1:i2,ju1:j2,:))
       !----------------------------------------------------------------

       !$OMP PARALLEL DO           &
       !$OMP DEFAULT( SHARED )     &
       !$OMP PRIVATE( IK, IJ, IL )
       do ik = k1,  k2
       do ij = ju1, j2
       do il = i1,  i2
          qqu(il,ij,ik) =  &
               qtmp(il,ij,ik) + (0.5d0 * qqu(il,ij,ik))
       
          qqv(il,ij,ik) =  &
               qtmp(il,ij,ik) + (0.5d0 * qqv(il,ij,ik))
       enddo
       enddo
       enddo
       !$OMP END PARALLEL DO
       
!   ======
    end if
!   ======
    
    
  END SUBROUTINE Calc_Advec_Cross_Terms
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Qckxyz
!
! !DESCRIPTION: Subroutine Qckxyz routine checks for "filling".
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Qckxyz( dq1, J1P, J2P,  JU1_GL, J2_GL, &
                     ILO, IHI, JULO, JHI,    I1,    &
                     I2,  JU1, J2,   K1,     K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: J1P,    J2P

    ! Global min & max latitude (J) indices
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI
!
! !INPUT/OUTPUT PARAMETERS: 
!
    ! Species density [hPa]
    REAL*8,  INTENT(INOUT) :: dq1(ILO:IHI, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO) 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
    LOGICAL, PARAMETER :: FILL_DIAG = .false.
!
! LOCAL VARIABLES:
!     
    INTEGER :: il, ij, ik
    INTEGER :: ip
    INTEGER :: k1p1, k2m1
    REAL*8  :: dup, qup
    REAL*8  :: qly
    REAL*8  :: sum


!     ----------------
!     Begin execution.
!     ----------------

    ip = 0


!     ----------
!     Top layer.
!     ----------

    k1p1 = k1 + 1

    !$OMP PARALLEL DO          &
    !$OMP DEFAULT( SHARED )    &
    !$OMP PRIVATE( IJ, IL, IP )
    do ij = j1p, j2p
       do il = i1, i2

          if (dq1(il,ij,k1) < 0.0d0) then

             ip = ip + 1

             dq1(il,ij,k1p1) = dq1(il,ij,k1p1) + dq1(il,ij,k1)
             dq1(il,ij,k1)   = 0.0d0
            
          end if

       end do
    end do
    !$OMP END PARALLEL DO

    
    do ik = k1 + 1, k2 - 1

       !$OMP PARALLEL DO                         &
       !$OMP DEFAULT( SHARED )                   &
       !$OMP PRIVATE( IJ, IL, IP, QUP, QLY, DUP )
       do ij = j1p, j2p
          do il = i1, i2
             
             if (dq1(il,ij,ik) < 0.0d0) then
                
                ip = ip + 1

!             -----------
!             From above.
!             -----------

                qup =  dq1(il,ij,ik-1)
                qly = -dq1(il,ij,ik)
                dup =  Min (qly, qup)

                dq1(il,ij,ik-1) = qup - dup
                dq1(il,ij,ik)   = dup - qly

!             -----------
!             From below.
!             -----------

                dq1(il,ij,ik+1) = dq1(il,ij,ik+1) + dq1(il,ij,ik)
                dq1(il,ij,ik)   = 0.0d0
                
             end if

          end do
       end do
       !$OMP END PARALLEL DO

    end do


!     -------------
!     Bottom layer.
!     -------------

    sum  = 0.0d0

    k2m1 = k2 - 1

    ! NOTE: Sum seems to be not used in the loop below!
    !$OMP PARALLEL DO                          &
    !$OMP DEFAULT( SHARED )                    &
    !$OMP PRIVATE( IJ, IL, IP, QUP, QLY, DUP ) &
    !$OMP REDUCTION( +:SUM ) 
    do ij = j1p, j2p
       do il = i1, i2

          if (dq1(il,ij,k2) < 0.0d0) then

             ip = ip + 1

!           -----------
!           From above.
!           -----------

             qup =  dq1(il,ij,k2m1)
             qly = -dq1(il,ij,k2)
             dup = Min (qly, qup)
             
             dq1(il,ij,k2m1) = qup - dup

!           -------------------------
!           From "below" the surface.
!           -------------------------

             sum = sum + qly - dup

             dq1(il,ij,k2) = 0.0d0

          end if

       end do
    end do
    !$OMP END PARALLEL DO

! We don't want to replace zero values by 1e-30. (ccc, 11/20/08)
!!     =======================================
!    where ((dq1(i1:i2,j1p:j2p,:) < 1.0d-30))  &
!         dq1(i1:i2,j1p:j2p,:) = 1.0d-30
!!     =======================================

  END SUBROUTINE Qckxyz
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Lmts
!
! !DESCRIPTION: Subroutine Set\_Lmts sets ILMT, JLMT, KLMT.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Lmts( ilmt, jlmt, klmt, I2_GL, J2_GL, iord, jord, kord )
!
! !INPUT PARAMETERS: 
!
    ! Global maximum longitude (I) and longitude (J) indices
    INTEGER, INTENT(IN)  :: I2_GL, J2_GL

    ! Flags to denote E-W, N-S, and vertical transport schemes
    ! (See REMARKS section of routine Tpcore_FvDas for more info)
    INTEGER, INTENT(IN)  :: iord, jord, kord
!
! !OUTPUT PARAMETERS:
!
    ! Controls various options in E-W advection
    INTEGER, INTENT(OUT) :: ilmt

    ! Controls various options in N-S advection
    INTEGER, INTENT(OUT) :: jlmt

    ! Controls various options in vertical advection
    INTEGER, INTENT(OUT) :: klmt
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO) 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

    INTEGER :: j2_glm1

!     ----------------
!     Begin execution.
!     ----------------

    j2_glm1 = j2_gl - 1

!c?
    if (IORD <= 0) then
       if (i2_gl >= 144) then
          ilmt = 0
       else if (i2_gl >= 72) then
          ilmt = 1
       else
          ilmt = 2
       end if
    else
       ilmt = IORD - 3
    end if
    
    
!c?
    if (JORD <= 0) then
       if (j2_glm1 >= 90) then
          jlmt = 0
       else if (j2_glm1 >= 45) then
          jlmt = 1
       else
          jlmt = 2
       end if
    else
       jlmt = JORD - 3
    end if

    klmt = Max ((KORD-3), 0)
    
  END SUBROUTINE Set_Lmts
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Set_Press_Terms
!
! !DESCRIPTION: Subroutine Set\_Press\_Terms sets the pressure terms: 
!  DELP1, DELPM, PU.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Set_Press_Terms( dap,   dbk,  pres1,  pres2, delp1,   &
                              delpm, pu,   JU1_GL, J2_GL, ILO,     &
                              IHI,   JULO, JHI,    J1P,   J2P,     &
                              I1,    I2,   JU1,    J2,    K1,  K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: J1P,    J2P

    ! Global min & max latitude (J) indices
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI

    ! Pressure difference across layer from (ai * pt) term [hPa]
    REAL*8,  INTENT(IN)  :: dap(K1:K2)

    ! Difference in bi across layer - the dSigma term
    REAL*8,  INTENT(IN)  :: dbk(K1:K2)

    ! Surface pressure at t1 [hPa]
    REAL*8,  INTENT(IN)  :: pres1(ILO:IHI, JULO:JHI)
    
    ! Surface pressure at t1+tdt [hPa]
    REAL*8,  INTENT(IN)  :: pres2(ILO:IHI, JULO:JHI)
!
! !OUTPUT PARAMETERS:
!
    ! Pressure thickness, the pseudo-density in a 
    ! hydrostatic system at t1 [hPa]
    REAL*8, INTENT(OUT) :: delp1(ILO:IHI, JULO:JHI, K1:K2)

    ! Pressure thickness, the pseudo-density in a 
    ! hydrostatic system at t1+tdt/2 (approximate) [hPa]
    REAL*8, INTENT(OUT) :: delpm(ILO:IHI, JULO:JHI, K1:K2)

    ! Pressure at edges in "u" [hPa]
    REAL*8, INTENT(OUT) :: pu(ILO:IHI, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO) 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: il, ij, ik

!   ----------------
!   Begin execution.
!   ----------------

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( IK     ) 
    do ik = k1, k2
    
       delp1(:,:,ik) = dap(ik) + (dbk(ik) * pres1(:,:))
    
       delpm(:,:,ik) =  &
                dap(ik) +  &
                (dbk(ik) * 0.5d0 * (pres1(:,:) + pres2(:,:)))
    
    end do
    !$OMP END PARALLEL DO
    
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( IL, IJ ) 
    do ij = j1p, j2p
       pu(1,ij,:) = 0.5d0 * (delpm(1,ij,:) + delpm(i2,ij,:))
       do il = i1+1, i2
          
          pu(il,ij,:) = 0.5d0 * (delpm(il,ij,:) + delpm(il-1,ij,:))
    
       end do
    end do
    !$OMP END PARALLEL DO


  END SUBROUTINE Set_Press_Terms
!EOC
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Courant
!
! !DESCRIPTION: Subroutine Calc\_Courant calculates courant numbers from 
!  the horizontal mass fluxes.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_Courant( cose, delpm, pu,     xmass, ymass, crx, cry,  & 
                           J1P,  J2P,   JU1_GL, J2_GL, ILO,   IHI, JULO, &
                           JHI,  I1,    I2,     JU1,   J2,    K1,  K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: J1P,    J2P

    ! Global min & max latitude (J) indices
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI

    ! Cosine of grid box edges
    REAL*8,  INTENT(IN)  :: cose (JU1_GL:J2_GL)

    ! Pressure thickness, the pseudo-density in a hydrostatic system
    ! at t1+tdt/2 (approximate) (mb)
    REAL*8,  INTENT(IN)  :: delpm(ILO:IHI, JULO:JHI, K1:K2)

    ! pressure at edges in "u"  (mb)
    REAL*8,  INTENT(IN)  :: pu   (iLO:IHI, JULO:JHI, K1:K2)

    ! horizontal mass flux in E-W and N-S directions [hPa]
    REAL*8,  INTENT(IN)  :: xmass(ILO:IHI, JULO:JHI, K1:K2)
    REAL*8,  INTENT(IN)  :: ymass(ILO:IHI, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Courant numbers in E-W and N-S directions
    REAL*8,  INTENT(OUT) :: crx(ILO:IHI, JULO:JHI, K1:K2)
    REAL*8,  INTENT(OUT) :: cry(ILO:IHI, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO) 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    INTEGER :: ij, ik

!   ----------------
!   Begin execution.
!   ----------------

    crx(:,:,:) = 0.0d0
    cry(:,:,:) = 0.0d0

!-----------------------------------------------------------------------------
! Prior to 12/4/08:
! We need to add an outer IK loop for OpenMP parallelization.  
! Preserve original code here! (bmy, 12/4/08)
!!   -----------------------------------
!!   Calculate E-W horizontal mass flux.
!!   -----------------------------------
!
!    do ij = j1p, j2p
!    
!       crx(:,ij,:) =  &
!            xmass(:,ij,:) / pu(:,ij,:)
!    
!    end do
!
!
!!   -----------------------------------
!!   Calculate N-S horizontal mass flux.
!!   -----------------------------------
!
!    do ij = j1p, j2p+1
!       
!       cry(:,ij,:) =  &
!            ymass(:,ij,:) /  &
!            ((0.5d0 * cose(ij)) *  &
!            (delpm(:,ij,:) + delpm(:,ij-1,:)))
!       
!    end do
!-----------------------------------------------------------------------------


    ! Rewrite DO loops above so that we can parallelize them.  Add the
    ! outer loop over IK and then parallelize over that. (bmy, 12/4/08)
    !
    ! NOTE: in the N-S mass flux loop, note that we have an index of
    ! an array that is IJ-1.  This is not allowed if IJ is the outermost
    ! loop and the loop we are parallelizing over.  However, if we 
    ! parallelize the higher-level IK loop then this is allowed. 
    ! (bmy, 12/4/08)

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( IJ, IK )
    do ik = k1, k2

!      -----------------------------------
!      Calculate E-W horizontal mass flux.
!      -----------------------------------

       do ij = j1p, j2p
    
          crx(:,ij,ik) =  &
               xmass(:,ij,ik) / pu(:,ij,ik)
    
       end do


!      -----------------------------------
!      Calculate N-S horizontal mass flux.
!      -----------------------------------

       do ij = j1p, j2p+1
       
          cry(:,ij,ik) =  &
               ymass(:,ij,ik) /  &
               ((0.5d0 * cose(ij)) *  &
               (delpm(:,ij,ik) + delpm(:,ij-1,ik)))
       
       end do
       
    enddo
    !$OMP END PARALLEL DO


  END SUBROUTINE Calc_Courant
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Calc_Divergence

!
! !DESCRIPTION: Subroutine Calc\_Divergence calculates the divergence.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Calc_Divergence( do_reduction, geofac_pc, geofac, dpi,   &
                              xmass,        ymass,     J1P,    J2P,   &
                              I1_GL,        I2_GL,     JU1_GL, J2_GL, &
                              ILO,          IHI,       JULO,   JHI,   &
                              I1,           I2,        JU1,    J2,    &
                              K1,           K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: J1P,    J2P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: I1_GL,  I2_GL
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI

    ! Set to F if called on Master or T if called by Slaves
    ! (NOTE: This is only for MPI parallelization, for OPENMP it should be F)
    LOGICAL, INTENT(IN)  :: do_reduction

    ! Special geometrical factor (geofac) for Polar cap
    REAL*8 , INTENT(IN)  :: geofac_pc

    ! Geometrical factor for meridional advection; geofac uses correct 
    ! spherical geometry, and replaces acosp as the meridional geometrical 
    ! factor in TPCORE
    REAL*8 , INTENT(IN)  :: geofac(JU1_GL:J2_GL)

    ! Horizontal mass flux in E/W and N/S directions [hPa]
    REAL*8 , INTENT(IN)  :: xmass(ILO:IHI, JULO:JHI, K1:K2)
    REAL*8 , INTENT(IN)  :: ymass(ILO:IHI, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Divergence at a grid point; used to calculate vertical motion [hPa]
    REAL*8,  INTENT(OUT) :: dpi(I1:I2, JU1:J2, K1:K2)
!      
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    INTEGER :: il, ij, ik

!   ----------------
!   Begin execution.
!   ----------------

!------------------------------------------------------------------------------
! Prior to 12/4/08:
! We need to add an outer IK loop for OpenMP parallelization.  
! Preserve original code here! (bmy, 12/4/08)
!   -------------------------
!   Calculate N-S divergence.
!!   -------------------------
!
!    do ij = j1p, j2p
!
!       dpi(:,ij,:) =  &
!            (ymass(:,ij,:) - ymass(:,ij+1,:)) *  &
!            geofac(ij)
!
!    end do
!
!
!!     -------------------------
!!     Calculate E-W divergence.
!!     -------------------------
!
!    do ij = j1p, j2p
!       do il = i1, i2-1
!          
!          dpi(il,ij,:) =  &
!               dpi(il,ij,:) +  &
!               xmass(il,ij,:) - xmass(il+1,ij,:)
!
!       end do
!       dpi(i2,ij,:) =  &
!            dpi(i2,ij,:) +  &
!            xmass(i2,ij,:) - xmass(1,ij,:)
!    end do
!------------------------------------------------------------------------------

    ! Rewrite DO loops above so that we can parallelize them.  Add the
    ! outer loop over IK and then parallelize over that. (bmy, 12/4/08)
    !
    ! NOTE: in the N-S mass flux loop, note that we have an index of
    ! an array that is IJ+1.  This is not allowed if IJ is the outermost
    ! loop and the loop we are parallelizing over.  However, if we 
    ! parallelize the higher-level IK loop then this is allowed. 
    ! (bmy, 12/4/08)

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( IL, IJ, IK )
    do ik = k1, k2

!      -------------------------
!      Calculate N-S divergence.
!      -------------------------
       
       do ij = j1p, j2p

          dpi(:,ij,ik) =  &
               (ymass(:,ij,ik) - ymass(:,ij+1,ik)) *  &
               geofac(ij)

       end do

!      -------------------------
!      Calculate E-W divergence.
!      -------------------------

       do ij = j1p, j2p
          do il = i1, i2-1
          
             dpi(il,ij,ik) =  &
                  dpi(il,ij,ik) +  &
                  xmass(il,ij,ik) - xmass(il+1,ij,ik)

          end do

          dpi(i2,ij,ik) =  &
               dpi(i2,ij,ik) +  &
               xmass(i2,ij,ik) - xmass(1,ij,ik)
       end do
    enddo
    !$OMP END PARALLEL DO
    

!   ===========================
    call Do_Divergence_Pole_Sum  &
!   ===========================
         (do_reduction, geofac_pc, dpi, ymass, &
         i1_gl, i2_gl, j1p, j2p, &
         ju1_gl, j2_gl, ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)


    if (j1p /= ju1_gl+1) then

!       --------------------------------------------
!       Polar cap enlarged:  copy dpi to polar ring.
!       --------------------------------------------
       
       !--------------------------------------------------------------
       ! Prior to 12/4/08: 
       ! We need to add an outer IK loop for OpenMP parallelization
       ! Preserve original code here! (bmy, 12/4/08)
       !dpi(:,ju1+1,:) = dpi(:,ju1,:)
       !dpi(:,j2-1,:)  = dpi(:,j2,:)
       !--------------------------------------------------------------

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK     )
       do ik = k1, k2
          dpi(:,ju1+1,ik) = dpi(:,ju1,ik)
          dpi(:,j2-1,ik)  = dpi(:,j2,ik)
       enddo
       !$OMP END PARALLEL DO

    end if


  END SUBROUTINE Calc_Divergence
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Divergence_Pole_Sum
!
! !DESCRIPTION: Subroutine Do\_Divergence\_Pole\_Sum sets the divergence 
!  at the Poles.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Divergence_Pole_Sum( do_reduction, geofac_pc, dpi, ymass, &
                                     I1_GL,        I2_GL,     J1P, J2P,   &
                                     JU1_GL,       J2_GL,     ILO, IHI,   &
                                     JULO,         JHI,       I1,  I2,    &
                                     JU1,          J2,        K1,  K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: J1P,    J2P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: I1_GL,  I2_GL
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI

    ! Set to T if called on Master or F if called by slaves
    ! NOTE: This seems not to be used here....)
    LOGICAL, INTENT(IN)   :: do_reduction

    ! Special geometrical factor (geofac) for Polar cap
    REAL*8,  INTENT(in)   :: geofac_pc

    ! Horizontal mass flux in N-S direction [hPa]
    REAL*8,  INTENT(IN)   :: ymass(ILO:IHI, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Divergence at a grid point; used to calculate vertical motion [hPa]
    REAL*8,  INTENT(OUT)  :: dpi(I1:I2, JU1:J2, K1:K2)

! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent. Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    INTEGER :: il, ik
    REAL*8  :: ri2
    REAL*8  :: mean_np(k1:k2)
    REAL*8  :: mean_sp(k1:k2)
    REAL*8  :: sumnp  (k1:k2)
    REAL*8  :: sumsp  (k1:k2)
   

!   ----------------
!   Begin execution.
!   ----------------

    ri2 = i2_gl


!   ==================
    if (ju1 == ju1_gl) then
!   ==================

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IL ) 
       do ik = k1, k2

          sumsp(ik) = 0.0d0

          do il = i1, i2

             sumsp(ik) = sumsp(ik) + ymass(il,j1p,ik)

          end do

       end do
       !$OMP END PARALLEL DO
       
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IL ) 
       do ik = k1, k2
          
          mean_sp(ik) = -sumsp(ik) / ri2 * geofac_pc

          do il = i1, i2

             dpi(il,ju1,ik) = mean_sp(ik)
             
          end do
          
       end do
       !$OMP END PARALLEL DO

!   ======
    end if
!   ======


!   ================
    if (j2 == j2_gl) then
!   ================

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IL ) 
       do ik = k1, k2

          sumnp(ik) = 0.0d0

          do il = i1, i2

             sumnp(ik) = sumnp(ik) + ymass(il,j2p+1,ik)
            
          end do

       end do
       !$OMP END PARALLEL DO

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IL ) 
       do ik = k1, k2

          mean_np(ik) = sumnp(ik) / ri2 * geofac_pc

          do il = i1, i2

             dpi(il,j2,ik) = mean_np(ik)

          end do

       end do
       !$OMP END PARALLEL DO

!   ======
    end if
!   ======

  END SUBROUTINE Do_Divergence_Pole_Sum
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Cross_Terms_Pole_I2d2
!
! !DESCRIPTION: Subroutine Do\_Cross\_Terms\_Pole\_I2d2 sets "va" at the Poles.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Cross_Terms_Pole_I2d2( cry,   va,  I1_GL, I2_GL, JU1_GL, &
                                       J2_GL, J1P, ILO,   IHI,   JULO,   &
                                       JHI,   I1,  I2,    JU1,   J2,     &
                                       K1,    K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edge of the South polar cap
    ! J1P=JU1_GL+1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: J1P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: I1_GL,  I2_GL
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI

    ! Courant number in N-S direction
    REAL*8,  INTENT(IN)  :: cry(ILO:IHI, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Average of Courant numbers from ij and ij+1
    REAL*8,  INTENT(OUT) :: va(ILO:IHI, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    INTEGER :: i2d2
    INTEGER :: il, ik

!   ----------------
!   Begin execution.
!   ----------------

    i2d2 = i2_gl / 2


!   ====================
    if (j1p == ju1_gl+1) then
!   ====================

!------------------------------------------------------------------------------
! Prior to 12/4/08:
! We need to add outer IK loops OpenMP parallelization.  
! Preserve original code here! (bmy, 12/4/08)
!!       ---------------------------------------------
!!       Polar Cap NOT Enlarged:
!!       Get cross terms for N-S horizontal advection.
!!       ---------------------------------------------
!
!!      ==================
!       if (ju1 == ju1_gl) then
!!      ==================
!
!          do il = i1, i2d2
!
!             va(il,ju1,:) =  &
!                  0.5d0 * (cry(il,ju1+1,:) - cry(il+i2d2,ju1+1,:))
!
!             va(il+i2d2,ju1,:) = -va(il,ju1,:)
!
!          end do
!
!!      ======
!       end if
!!      ======
!
!
!!      ================
!       if (j2 == j2_gl) then
!!      ================
!
!          do il = i1, i2d2
!
!             va(il,j2,:) =  &
!                  0.5d0 * (cry(il,j2,:) - cry(il+i2d2,j2-1,:))
!
!             va(il+i2d2,j2,:) = -va(il,j2,:)
!
!          end do
!
!!      ======
!       end if
!!      ======
!
!!   ======
!    end if
!!   ======
!------------------------------------------------------------------------------

       ! Add outermost IK loops so that we can parallelize
       ! (bmy, 12/4/08) 

!      ---------------------------------------------
!      Polar Cap NOT Enlarged:
!      Get cross terms for N-S horizontal advection.
!      ---------------------------------------------

!      ==================
       if (ju1 == ju1_gl) then
!      ==================

          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( IK, IL )
          do ik = k1, k2
          do il = i1, i2d2

             va(il,ju1,ik) =  &
                  0.5d0 * (cry(il,ju1+1,ik) - cry(il+i2d2,ju1+1,ik))

             va(il+i2d2,ju1,ik) = -va(il,ju1,ik)

          end do
          end do
          !$OMP END PARALLEL DO

!      ======
       end if
!      ======


!      ================
       if (j2 == j2_gl) then
!      ================

          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( IK, IL )
          do ik = k1, k2
          do il = i1, i2d2

             va(il,j2,ik) =  &
                  0.5d0 * (cry(il,j2,ik) - cry(il+i2d2,j2-1,ik))

             va(il+i2d2,j2,ik) = -va(il,j2,ik)

          end do
          end do
          !$OMP END PARALLEL DO

!      ======
       end if
!      ======

!   ======
    end if
!   ======


  END SUBROUTINE Do_Cross_Terms_Pole_I2d2
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Xadv_Dao2
!
! !DESCRIPTION: Subroutine Xadv\_Dao2 is the advective form E-W operator for 
!  computing the adx (E-W) cross term.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Xadv_Dao2( iad,    jn,    js,  adx,  qqv, &
                        ua,     ILO,   IHI, JULO, JHI, &
                        JU1_GL, J2_GL, J1P, J2P,  I1,  &
                        I2,     JU1,   J2,  K1,   K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: J1P,    J2P

    ! Global min & max latitude (J) indices
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI

    ! if iad = 1, use 1st order accurate scheme;
    ! if iad = 2, use 2nd order accurate scheme
    INTEGER, INTENT(IN)  :: iad

    ! Northward of latitude index = jn, Courant numbers could be > 1,
    ! so use the flux-form semi-Lagrangian scheme
    INTEGER, INTENT(IN)  :: jn(K1:K2)

    ! southward of latitude index = js, Courant numbers could be > 1,
    ! so use the flux-form semi-Lagrangian scheme
    INTEGER, INTENT(IN)  :: js(K1:K2)
 
    ! Concentration contribution from N-S advection [mixing ratio]
    REAL*8,  INTENT(IN)  :: qqv(ILO:IHI, JULO:JHI, K1:K2)

    ! Average of Courant numbers from il and il+1
    REAL*8,  INTENT(IN)  :: ua(ILO:IHI, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Cross term due to E-W advection [mixing ratio]
    REAL*8,  INTENT(OUT) :: adx(ILO:IHI, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: il, ij, ik, iu
    INTEGER :: imp, iue, iuw
    REAL*8  :: a1, b1, c1
    REAL*8  :: rdiff
    REAL*8  :: ril, riu
    real*8  :: ru
    
    ! Arrays
    REAL*8  :: qtmp(-i2/3:i2+i2/3, julo:jhi, k1:k2)
      
    !     ----------------
    !     Begin execution.
    !     ----------------
    
    ! Zero output array
    adx = 0d0

    !-----------------------------------------------------------------------
    ! Prior to 12/5/08:
    ! We need to add outer IJ and IK loops for OpenMP parallelization.
    ! Preserve original code here. (bmy, 12/5/08)
    !do il=1,i2
    !   qtmp(il,:,:) = qqv(il,:,:)
    !enddo
    !  
    !do il=-i2/3,0
    !   qtmp(il,:,:) = qqv(i2+il,:,:)
    !enddo
    !
    !do il=i2+1,i2+i2/3
    !   qtmp(il,:,:) = qqv(il-i2,:,:)
    !enddo
    !-----------------------------------------------------------------------

    !$OMP PARALLEL DO          &
    !$OMP DEFAULT( SHARED )    &
    !$OMP PRIVATE( IK, IJ, IL )
    do ik = k1, k2
       do ij = julo, jhi
  
          do il=1,i2
             qtmp(il,ij,ik) = qqv(il,ij,ik)
          enddo
          
          do il=-i2/3,0
             qtmp(il,ij,ik) = qqv(i2+il,ij,ik)
          enddo
          
          do il=i2+1,i2+i2/3
             qtmp(il,ij,ik) = qqv(il-i2,ij,ik)
          enddo
       enddo
    enddo
    !$OMP END PARALLEL DO

!   =============
    if (iad == 1) then
!   =============
         
       !       ----------
       !       1st order.
       !       ----------

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IJ, IL, IU, RIU, RU, RDIFF, RIL )
       do ik = k1, k2
          do ij = j1p, j2p
               

             if ((ij <= js(ik)) .or. (ij >= jn(ik))) then
                  
                !             --------------
                !             In Polar area.
                !             --------------
                  
                do il = i1, i2
                     
                   iu  = ua(il,ij,ik)
                   riu = iu
                   ru  = ua(il,ij,ik) - riu
                   iu  = il - iu
                     
                   if (ua(il,ij,ik) >= 0.0d0) then
                      rdiff = qtmp(iu-1,ij,ik) - qtmp(iu,ij,ik)
                   else
                      rdiff = qtmp(iu,ij,ik)   - qtmp(iu+1,ij,ik)
                   end if
                     
                   adx(il,ij,ik) = (qtmp(iu,ij,ik) - qtmp(il,ij,ik)) +  &
                                   (ru * rdiff)
                     
                end do
                  
             else  ! js < ij < jn
                  
                !             ----------------
                !             Eulerian upwind.
                !             ----------------
                  
                do il = i1, i2
                     
                   ril = il
                   iu  = ril - ua(il,ij,ik)
                     
                   adx(il,ij,ik) = ua(il,ij,ik) *  &
                        (qtmp(iu,ij,ik) - qtmp(iu+1,ij,ik))
                     
                end do
                  
             end if
               
          end do
       end do
       !$OMP END PARALLEL DO
         
         
!   ==================
    else if (iad == 2) then
!   ==================
         

       !$OMP PARALLEL DO                                   &
       !$OMP DEFAULT( SHARED )                             &
       !$OMP PRIVATE( IK, IJ, IL, IU, RIU, RU, A1, B1, C1 )
       do ik = k1, k2
          do ij = j1p, j2p
               

             if ((ij <= js(ik)) .or. (ij >= jn(ik))) then
                  
                !             --------------
                !             In Polar area.
                !             --------------
                  
                do il = i1, i2
                     
                   iu  = Nint (ua(il,ij,ik))
                   riu = iu
                   ru  = riu - ua(il,ij,ik)
                   iu  = il - iu
                     
                   a1 = 0.5d0 * (qtmp(iu+1,ij,ik) + qtmp(iu-1,ij,ik)) -  &
                                 qtmp(iu,ij,ik)
                     
                   b1 = 0.5d0 * (qtmp(iu+1,ij,ik) - qtmp(iu-1,ij,ik))
                     
                   c1 = qtmp(iu,ij,ik) - qtmp(il,ij,ik)
                     
                   adx(il,ij,ik) = (ru * ((a1 * ru) + b1)) + c1
                     
                end do
                  
             else  ! js < ij < jn
                  
                !             ----------------
                !             Eulerian upwind.
                !             ----------------
                  
                do il = i1, i2
                   
                   iu  = Nint (ua(il,ij,ik))
                   riu = iu
                   ru  = riu - ua(il,ij,ik)
                   iu  = il - iu
                     
                   a1 = 0.5d0 * (qtmp(iu+1,ij,ik) + qtmp(iu-1,ij,ik)) -  &
                                 qtmp(iu,ij,ik)
                   
                   b1 = 0.5d0 * (qtmp(iu+1,ij,ik) - qtmp(iu-1,ij,ik))
                     
                   c1 = qtmp(iu,ij,ik) - qtmp(il,ij,ik)
                     
                   adx(il,ij,ik) = (ru * ((a1 * ru) + b1)) + c1
                     
                end do
                  
             end if
               
          end do
       end do
       !$OMP END PARALLEL DO
         
!   ======
    end if
!   ======
      
      
    if (ju1 == ju1_gl) then
       
       !---------------------------------------------------------------
       ! Prior to 12/4/08:
       ! We need to rewrite the DO loop below for OpenMP.  
       ! Preserve original code here! (bmy, 12/4/08)
       !adx(i1:i2,ju1,:) = 0.0d0
       !  
       !if (j1p /= ju1_gl+1) then
       !     
       !   adx(i1:i2,ju1+1,:) = 0.0d0
       !     
       !end if
       !---------------------------------------------------------------

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK     )
       do ik = k1, k2

          adx(i1:i2,ju1,ik) = 0.0d0
         
          if (j1p /= ju1_gl+1) then
            
             adx(i1:i2,ju1+1,ik) = 0.0d0
            
          end if
          
       enddo
       !$OMP END PARALLEL DO
  
    end if
      
      
    if (j2 == j2_gl) then

       !---------------------------------------------------------------
       ! Prior to 12/4/08:
       ! We need to rewrite the DO loop below for OpenMP.  
       ! Preserve original code here! (bmy, 12/4/08)
       !adx(i1:i2,j2,:) = 0.0d0
       !  
       !if (j1p /= ju1_gl+1) then
       !     
       !   adx(i1:i2,j2-1,:) = 0.0d0
       !     
       !end if
       !---------------------------------------------------------------

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK     )
       do ik = k1, k2

          adx(i1:i2,j2,ik) = 0.0d0
         
          if (j1p /= ju1_gl+1) then
            
             adx(i1:i2,j2-1,ik) = 0.0d0
            
          end if

       enddo
       !$OMP END PARALLEL DO

    end if
      
  
  END SUBROUTINE Xadv_Dao2
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Yadv_Dao2
!
! !DESCRIPTION: Subroutine Yadv\_Dao2 is the advective form N-S operator 
!  for computing the ady (N-S) cross term.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Yadv_Dao2( iad,   ady,    qqu,   va,  I1_GL, &
                        I2_GL, JU1_GL, J2_GL, J1P, J2P,   &
                        ILO,   IHI,    JULO,  JHI, I1,    &
                        I2,    JU1,    J2,    K1,  K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: J1P,    J2P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: I1_GL,  I2_GL
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI

    ! If iad = 1, use 1st order accurate scheme;
    ! If iad = 2, use 2nd order accurate scheme
    INTEGER, INTENT(IN)  :: iad

    ! Concentration contribution from E-W advection [mixing ratio]
    REAL*8,  INTENT(IN)  :: qqu(ILO:IHI, JULO:JHI, K1:K2)

    ! Average of Courant numbers from ij and ij+1
    REAL*8,  INTENT(IN)  :: va(ILO:IHI, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Cross term due to N-S advection (mixing ratio)
    REAL*8,  INTENT(OUT) :: ady(ILO:IHI, JULO:JHI, K1:K2)

! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: il, ij, ik
    INTEGER :: jv
    REAL*8  :: a1, b1, c1
    REAL*8  :: rij, rjv
    REAL*8  :: rv

    ! Arrays 
    ! We may need a small ghost zone depending 
    ! on the polar cap used
    REAL*8  :: qquwk(ilo:ihi, julo-2:jhi+2, k1:k2)   
                                                       
!     ----------------
!     Begin execution.
!     ----------------

    ! Zero output array
    ady = 0d0

    ! Make work array
    do ij = julo, jhi
       qquwk(:,ij,:) = qqu(:,ij,:)
    end do


    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! This routine creates a ghost zone in latitude in case of 
    ! not enlarged polar cap 
    ! (ccc, 11/20/08)
!   ======================
    call Do_Yadv_Pole_I2d2 &
!   ======================
         (qqu, qquwk, &
          i1_gl, i2_gl, ju1_gl, j2_gl, j1p, &
          ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)

      
!   =============
    if (iad == 1) then
!   =============
         
       !       ----------
       !       1st order.
       !       ----------
         
       !$OMP PARALLEL DO                   &
       !$OMP DEFAULT( SHARED )             &
       !$OMP PRIVATE( IK, IJ, IL, RIJ, JV ) 
       do ik = k1, k2
          do ij = j1p-1, j2p+1
             do il = i1, i2
!c?
                rij = ij
                jv  = rij - va(il,ij,ik)
                  
                ady(il,ij,ik) = va(il,ij,ik) *  &
                     (qquwk(il,jv,ik) - qquwk(il,jv+1,ik))
                  
             end do
          end do
       end do
       !$OMP END PARALLEL DO
         
         
!   ==================
    else if (iad == 2) then
!   ==================
         
       !$OMP PARALLEL DO                                   &
       !$OMP DEFAULT( SHARED )                             &
       !$OMP PRIVATE( IK, IJ, IL, JV, RJV, RV, A1, B1, C1 )
       do ik = k1, k2
          do ij = j1p-1, j2p+1
             do il = i1, i2
!c?
                jv  = Nint (va(il,ij,ik))
                rjv = jv
                rv  = rjv - va(il,ij,ik)
                jv  = ij - jv
                  
                a1 = 0.5d0 * (qquwk(il,jv+1,ik) + qquwk(il,jv-1,ik)) -  &
                              qquwk(il,jv,ik)
                  
                b1 = 0.5d0 * (qquwk(il,jv+1,ik) - qquwk(il,jv-1,ik))
                  
                c1 = qquwk(il,jv,ik) - qquwk(il,ij,ik)
                  
                ady(il,ij,ik) = (rv * ((a1 * rv) + b1)) + c1
                  
             end do
          end do
       end do
       !$OMP END PARALLEL DO
         
    end if
      
      
!   =====================
    call Do_Yadv_Pole_Sum &
!   =====================
         ( ady, &
           i1_gl, i2_gl, ju1_gl, j2_gl, j1p, &
           ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)
      

  END SUBROUTINE Yadv_Dao2
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Yadv_Pole_I2d2
!
! !DESCRIPTION: Subroutine Do\_Yadv\_Pole\_I2d2 sets "qquwk" at the Poles.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Yadv_Pole_I2d2 ( qqu, qquwk, I1_GL, I2_GL, JU1_GL, J2_GL, &
                                 J1P, ILO,   IHI,   JULO,  JHI,    I1,    &
                                 I2,  JU1,   J2,    K1,    K2  )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the South polar cap
    ! J1P=JU1_GL+1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: J1P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: I1_GL,  I2_GL
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI
      
    ! concentration contribution from E-W advection [mixing ratio]
    REAL*8,  INTENT(IN)  :: qqu(ILO:IHI, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! qqu working array [mixing ratio]
    REAL*8,  INTENT(OUT) :: qquwk(ILO:IHI, JULO-2:JHI+2, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: i2d2
    INTEGER :: il, ij, ik
    INTEGER :: inb
      
      
!   ----------------
!   Begin execution.
!   ----------------

    i2d2 = i2_gl / 2


!   ====================
    if (j1p == ju1_gl+1) then
!   ====================

!      -----------------------
!      Polar Cap NOT Enlarged.
!      -----------------------

!      ==================
       if (ju1 == ju1_gl) then
!      ==================

          !-----------------------------------------------------------------
          ! Prior to 12/4/08:
          ! We need to add an outer IK loop for OpenMP parallelization
          ! Preserve original code here (bmy, 12/4/08)
          !do il = i1, i2d2
          !   do inb = 1, 2
          !       
          !      qquwk(il,     ju1-inb,:) = qqu(il+i2d2,ju1+inb,:)
          !      qquwk(il+i2d2,ju1-inb,:) = qqu(il,     ju1+inb,:)
          !       
          !   end do
          !end do
          !-----------------------------------------------------------------

          !$OMP PARALLEL DO            &
          !$OMP DEFAULT( SHARED      ) &    
          !$OMP PRIVATE( IK, IL, INB ) 
          do ik = k1, k2
          do il = i1, i2d2
             do inb = 1, 2
                 
                qquwk(il,     ju1-inb,ik) = qqu(il+i2d2,ju1+inb,ik)
                qquwk(il+i2d2,ju1-inb,ik) = qqu(il,     ju1+inb,ik)
                 
             end do
          end do
          end do
          !$OMP END PARALLEL DO


!      ======
       end if
!      ======


!      ================
       if (j2 == j2_gl) then
!      ================

          !-----------------------------------------------------------------
          ! Prior to 12/4/08:
          ! We need to add an outer IK loop for OpenMP parallelization
          ! Preserve original code here (bmy, 12/4/08)
          !do il = i1, i2d2
          !   do inb = 1, 2
          !       
          !      qquwk(il,     j2+inb,:) = qqu(il+i2d2,j2-inb,:)
          !      qquwk(il+i2d2,j2+inb,:) = qqu(il,     j2-inb,:)
          !       
          !   end do
          !end do
          !-----------------------------------------------------------------

          !$OMP PARALLEL DO            &
          !$OMP DEFAULT( SHARED      ) &    
          !$OMP PRIVATE( IK, IL, INB ) 
          do ik = k1, k2
          do il = i1, i2d2
             do inb = 1, 2
                 
                qquwk(il,     j2+inb,ik) = qqu(il+i2d2,j2-inb,ik)
                qquwk(il+i2d2,j2+inb,ik) = qqu(il,     j2-inb,ik)
                 
             end do
          end do
          end do
          !$OMP END PARALLEL DO


!      ======
       end if
!      ======

!   ======
    end if
!   ======


  END SUBROUTINE Do_Yadv_Pole_I2d2
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Yadv_Pole_Sum
!
! !DESCRIPTION: Subroutine Do\_Yadv\_Pole\_Sum sets the cross term due to 
!  N-S advection at the Poles.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Yadv_Pole_Sum( ady, I1_GL, I2_GL, JU1_GL, J2_GL, J1P, &
                               ILO, IHI,   JULO,  JHI,    I1,    I2,  &
                               JU1, J2,    K1,    K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude index at the edge of the South polar cap
    ! J1P=JU1_GL+1; for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)    :: J1P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)    :: I1_GL,  I2_GL
    INTEGER, INTENT(IN)    :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)    :: I1,     I2
    INTEGER, INTENT(IN)    :: JU1,    J2
    INTEGER, INTENT(IN)    :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)    :: ILO,    IHI
    INTEGER, INTENT(IN)    :: JULO,   JHI
!
! !OUTPUT PARAMETERS:
!
    ! Cross term due to N-S advection (mixing ratio)
    REAL*8,  INTENT(INOUT) :: ady(ILO:IHI, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.  Also make a logical
!                               to test if we are using an extended polar cap.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 

    ! Scalars
    INTEGER :: il, ik
      
    ! Arrays
    REAL*8  :: sumnp(k1:k2)
    REAL*8  :: sumsp(k1:k2)
      
    ! Add 
    LOGICAL :: IS_EXT_POLAR_CAP
    

    ! ----------------
    ! Begin execution.
    ! ----------------

    ! Test if we are using extended polar caps (i.e. the S pole and next N 
    ! latitude and N. Pole and next S latitude).  Do this outside the loops.
    ! (bmy, 12/11/08)
    IS_EXT_POLAR_CAP = ( J1P == JU1_GL + 2 )

    ! ------------
    !  South Pole
    ! ------------

!   ==================
    if (ju1 == ju1_gl) then
!   ==================

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IL )
       do ik = k1, k2
            
          sumsp(ik) = 0.0d0
            
          !------------------------------
          ! Prior to 12/11/08:
          !if (j1p /= ju1_gl+1) then
          !------------------------------
          if ( IS_EXT_POLAR_CAP ) then

             ! For a 2-latitude polar cap (S. Pole + next Northward latitude)
             do il = i1, i2
                  
                sumsp(ik) = sumsp(ik) + ady(il,ju1+1,ik)
                  
             end do
               
          else
               
             ! For a 1-latitude polar cap (S. Pole only)
             do il = i1, i2
                  
                sumsp(ik) = sumsp(ik) + ady(il,ju1,ik)
                  
             end do
               
          end if
            
       end do
       !$OMP END PARALLEL DO
         
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IL )
       do ik = k1, k2

          sumsp(ik) = sumsp(ik) / i2_gl
          
          !------------------------------
          ! Prior to 12/11/08:
          !if (j1p /= ju1_gl+1) then
          !------------------------------
          if ( IS_EXT_POLAR_CAP ) then
             
             ! For a 2-latitude polar cap (S. Pole + next Northward latitude)
             do il = i1, i2
                
                ady(il,ju1+1,ik) = sumsp(ik)
                ady(il,ju1,ik)   = sumsp(ik)
                  
             end do
               
          else
               
             ! For a 1-latitude polar cap (S. Pole only)
             do il = i1, i2
                  
                ady(il,ju1,ik) = sumsp(ik)
                  
             end do
               
          end if
            
       end do
       !$OMP END PARALLEL DO
  
!   ======
    end if
!   ======

    ! ------------
    !  North Pole
    ! ------------

!   ================
    if (j2 == j2_gl) then
!   ================

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IL )
       do ik = k1, k2
            
          sumnp(ik) = 0.0d0

          !-------------------------------
          ! Prior to 12/11/08:
          !if (j1p /= ju1_gl+1) then
          !-------------------------------
          if ( IS_EXT_POLAR_CAP ) then
               
             ! For a 2-latitude polar cap (N. Pole + next Southward latitude)
             do il = i1, i2
                  
                sumnp(ik) = sumnp(ik) + ady(il,j2-1,ik)
                  
             end do
               
          else
               
             ! For a 1-latitude polar cap (N. Pole only)
             do il = i1, i2
                  
                sumnp(ik) = sumnp(ik) + ady(il,j2,ik)
                  
             end do
               
          end if
            
       end do
       !$OMP END PARALLEL DO
         
       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IL )
       do ik = k1, k2
            
          sumnp(ik) = sumnp(ik) / i2_gl
            
          !-------------------------------
          ! Prior to 12/11/08:
          !if (j1p /= ju1_gl+1) then
          !-------------------------------
          if ( IS_EXT_POLAR_CAP ) then

             ! For a 2-latitude polar cap (N. Pole + next Southward latitude)
             do il = i1, i2
                  
                ady(il,j2-1,ik) = sumnp(ik)
                ady(il,j2,ik)   = sumnp(ik)
                  
             end do
               
          else
               
             ! For a 1-latitude polar cap (N. Pole only)
             do il = i1, i2
                  
                ady(il,j2,ik) = sumnp(ik)
                  
             end do
               
          end if
            
       end do
       !$OMP END PARALLEL DO
         
!   ======
    end if
!   ======


  END SUBROUTINE Do_Yadv_Pole_Sum
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Xtp
!
! !DESCRIPTION: Subroutine Xtp does horizontal advection in the E-W direction.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Xtp( ilmt, jn,  js,    pu,     crx,   dq1, qqv, xmass, fx,  &
                  J1P,  J2P, I2_GL, JU1_GL, J2_GL, ILO, IHI, JULO,  JHI, &
                  I1,   I2,  JU1,   J2,     K1,    K2,  iord )

!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)    :: J1P,    J2P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)    ::         I2_GL
    INTEGER, INTENT(IN)    :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)    :: I1,     I2
    INTEGER, INTENT(IN)    :: JU1,    J2
    INTEGER, INTENT(IN)    :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)    :: ILO,    IHI
    INTEGER, INTENT(IN)    :: JULO,   JHI

    ! Controls various options in E-W advection
    INTEGER, INTENT(IN)    :: ilmt

    ! Northward of latitude index = jn, Courant numbers could be > 1,
    ! so use the flux-form semi-Lagrangian scheme
    INTEGER, INTENT(IN)    :: jn(K1:K2)

    ! Southward of latitude index = js, Courant numbers could be > 1,
    ! so use the flux-form semi-Lagrangian scheme
    INTEGER, INTENT(IN)    :: js(K1:K2)

    ! Option for E-W transport scheme.  See module header for more info.
    INTEGER, INTENT(IN)    :: iord

    ! pressure at edges in "u" [hPa]
    REAL*8,  INTENT(IN)    :: pu(ILO:IHI, JULO:JHI, K1:K2)

    ! Courant number in E-W direction
    REAL*8,  INTENT(IN)    :: crx(ILO:IHI, JULO:JHI, K1:K2)

    ! Horizontal mass flux in E-W direction [hPa]
    REAL*8,  INTENT(IN)    :: xmass(ILO:IHI, JULO:JHI, K1:K2)
!
! !INPUT/OUTPUT PARAMETERS: 
!
    ! Species density [hPa]
    REAL*8,  INTENT(INOUT) :: dq1(ILO:IHI, JULO:JHI, K1:K2)

    ! Concentration contribution from N-S advection [mixing ratio]
    REAL*8,  INTENT(INOUT) :: qqv(ILO:IHI, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! E-W flux [mixing ratio]
    REAL*8,  INTENT(OUT)   :: fx(ILO:IHI, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: il, ij, ik, ic
    INTEGER :: iu, ix, iuw, iue, imp
    INTEGER :: jvan
    REAL*8  :: rc
    REAL*8  :: ric, ril

    ! Arrays
    INTEGER :: isav(i1:i2)
    REAL*8  :: dcx(-i2/3:i2+i2/3, julo:jhi, k1:k2)
    REAL*8  :: qtmp(-i2/3:i2+i2/3, julo:jhi, k1:k2)


    !     ----------------
    !     Begin execution.
    !     ----------------

    dcx(:,:,:) = 0.0d0
    fx(:,:,:) = 0.0d0
    
    imp = i2+1
    
    ! NOTE: these loops do not parallelize well (bmy, 12/5/08)

    !!ccc Populate qtmp
    do il=i1,i2
       qtmp(il,:,:) = qqv(il,:,:)
    enddo
    
    do il = -i2/3,0
       qtmp(il,:,:) = qqv(i2+il,:,:)
    enddo
    
    do il = i2+1,i2+i2/3
       qtmp(il,:,:) = qqv(il-i2,:,:)
    enddo

    if (IORD /= 1) then
       qtmp(i1-1,:,:) = qqv(i2,:,:)
       qtmp(i1-2,:,:) = qqv(i2-1,:,:)
       qtmp(i2+1,:,:) = qqv(i1,:,:)
       qtmp(i2+2,:,:) = qqv(i1+1,:,:)

!      ==========
       call Xmist &
!      ==========
            (dcx, qtmp, &
             j1p, j2p, i2_gl, ju1_gl, j2_gl, ilo, ihi, julo, jhi, &
             i1, i2, ju1, j2, k1, k2)
    end if


    jvan = Max (1, j2_gl / 18)


!   ==============
    !$OMP PARALLEL DO                                          &
    !$OMP DEFAULT( SHARED )                                    &
    !$OMP PRIVATE( IK, IJ, IL, RIL, IU, IC, ISAV, RIC, RC, IX )
    do ik = k1, k2
       do ij = j1p, j2p
!      =================

!         ======================================
          if ((ij > js(ik)) .and. (ij < jn(ik))) then
!         ======================================

!           ------------------------------------------------------
!           Do horizontal Eulerian advection in the E-W direction.
!           ------------------------------------------------------

             if ((IORD == 1) .or.  &
                  (ij == j1p) .or. (ij == j2p)) then
                  
                do il = i1, i2
                   ril = il
                   iu  = ril - crx(il,ij,ik)
                   
                   fx(il,ij,ik) = qtmp(iu,ij,ik)
                end do
                  
             else
                  
                if ((IORD == 2) .or.  &
                     (ij <= (j1p+jvan)) .or. (ij >= (j2p-jvan))) then
                     
                   do il = i1, i2
                      ril = il
                      iu  = ril - crx(il,ij,ik)
                      
                      fx(il,ij,ik) =  &
                           qtmp(iu,ij,ik) +  &
                           (dcx(iu,ij,ik) *  &
                           (Sign (1.0d0, crx(il,ij,ik)) - crx(il,ij,ik)))
                   end do
                     
                else

!                  ==========
                   call Fxppm  &
                        (ij, ik, ilmt, crx, dcx, fx, qtmp, &
                        -i2/3, i2+i2/3, julo, jhi, i1, i2, k1, k2)
!  qtmp (inout) - can be updated
!                  ==========

                end if

             end if

             !---------------------------------------------------------------
             ! Prior to 12/5/08:
             ! We need to write this as an explicit loop over IL
             ! to facilitate OpenMP parallelization.  Preserve original
             ! code here. (bmy, 12/5/08)
             !fx(i1:i2,ij,ik) = fx(i1:i2,ij,ik) * xmass(i1:i2,ij,ik)
             !---------------------------------------------------------------
             do il = i1, i2  
                fx(il,ij,ik) = fx(il,ij,ik) * xmass(il,ij,ik)
             enddo

!         ====
          else
!         ====

!           ------------------------------------------------------------
!           Do horizontal Conservative (flux-form) Semi-Lagrangian
!           advection in the E-W direction (van Leer at high latitudes).
!           ------------------------------------------------------------

             if ((IORD == 1) .or.  &
                  (ij == j1p) .or. (ij == j2p)) then

                do il = i1, i2
                   ic       = crx(il,ij,ik)
                   isav(il) = il - ic
                   ril      = il
                   iu       = ril - crx(il,ij,ik)
                   ric      = ic
                   rc       = crx(il,ij,ik) - ric
                   
                   fx(il,ij,ik) = rc * qtmp(iu,ij,ik)
                end do
                
             else
                  
                do il = i1, i2
                   ic       = crx(il,ij,ik)
                   isav(il) = il - ic
                   ril      = il
                   iu       = ril - crx(il,ij,ik)
                   ric      = ic
                   rc       = crx(il,ij,ik) - ric
                     
                   fx(il,ij,ik) =  &
                        rc *  &
                        (qtmp(iu,ij,ik) +  &
                        (dcx(iu,ij,ik) * (Sign (1.0d0, rc) - rc)))
                end do

             end if
               
             do il = i1, i2
                  
                if (crx(il,ij,ik) > 1.0d0) then

                   do ix = isav(il), il - 1
                      fx(il,ij,ik) = fx(il,ij,ik) + qtmp(ix,ij,ik)
                   end do

                else if (crx(il,ij,ik) < -1.0d0) then
                                        
                   do ix = il, isav(il) - 1
                      fx(il,ij,ik) = fx(il,ij,ik) - qtmp(ix,ij,ik)
                   end do
                     
                end if

             end do

             !---------------------------------------------------------------
             ! Prior to 12/5/08:
             ! We need to write this as an explicit loop over IL
             ! to facilitate OpenMP parallelization.  Preserve original
             ! code here. (bmy, 12/5/08)
             !fx(i1:i2,ij,ik) = pu(i1:i2,ij,ik) * fx(i1:i2,ij,ik)
             !---------------------------------------------------------------
             do il = i1, i2
                fx(il,ij,ik) = pu(il,ij,ik) * fx(il,ij,ik)
             enddo

!         ======
          end if
!         ======

!      ======
       end do
    end do
    !$OMP END PARALLEL DO
!   ======

    ! NOTE: This loop does not parallelize well (bmy, 12/5/08)
    do ij = j1p, j2p
       do il = i1, i2-1
            
          dq1(il,ij,:) = dq1(il,ij,:) +  &
                         (fx(il,ij,:) - fx(il+1,ij,:))
    
       end do
       dq1(i2,ij,:) = dq1(i2,ij,:) +  &
                      (fx(i2,ij,:) - fx(i1,ij,:))
    end do
    
  END SUBROUTINE Xtp
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Xmist
!
! !DESCRIPTION: Subroutine Xmist computes the linear tracer slope in the 
!  E-W direction. It uses the Lin et. al. 1994 algorithm.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Xmist( dcx,  qqv, J1P, J2P, I2_GL, JU1_GL, J2_GL, ILO, IHI, &
                    JULO, JHI, I1,  I2,  JU1,   J2,     K1,    K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: J1P,    J2P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  ::         I2_GL
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI

    ! Concentration contribution from N-S advection [mixing ratio]
    REAL*8,  INTENT(IN)  :: qqv(-I2/3:I2+I2/3, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Slope of concentration distribution in E-W direction [mixing ratio]
    REAL*8,  INTENT(OUT) :: dcx(-I2/3:I2+I2/3, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: il, ij, ik
    REAL*8  :: pmax, pmin
    REAL*8  :: r24
    REAL*8  :: tmp
      

!   ----------------
!   Begin execution.
!   ----------------

    r24 = 1.0d0 / 24.0d0

    !$OMP PARALLEL DO                            &
    !$OMP DEFAULT( SHARED )                      &
    !$OMP PRIVATE( IK, IJ, IL, TMP, PMAX, PMIN )
    do ik = k1, k2
       do ij = j1p+1, j2p-1
          do il = i1, i2

             tmp =  &
                  ((8.0d0 * (qqv(il+1,ij,ik) - qqv(il-1,ij,ik))) +  &
                  qqv(il-2,ij,ik) - qqv(il+2,ij,ik)) *  &
                  r24

             pmax =  &
                  Max (qqv(il-1,ij,ik), qqv(il,ij,ik), qqv(il+1,ij,ik)) -  &
                  qqv(il,ij,ik)

             pmin =  &
                  qqv(il,ij,ik) -  &
                  Min (qqv(il-1,ij,ik), qqv(il,ij,ik), qqv(il+1,ij,ik))

             dcx(il,ij,ik) = Sign (Min (Abs (tmp), pmax, pmin), tmp)

          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !--------------------------------------------------------------------
    ! Prior to 12/4/08:
    ! We need to add outer IK and IJ loops for OpenMP parallelization.
    ! Preserve original code here (bmy, 12/4/08)
    !! Populate ghost zones of dcx (ccc, 11/20/08)
    !do il = -i2/3, 0
    !   dcx(il,:,:) = dcx(i2+il,:,:)
    !enddo
    !  
    !do il = i2+1, i2+i2/3
    !   dcx(il,:,:) = dcx(il-i2,:,:)
    !enddo
    !--------------------------------------------------------------------

    ! Populate ghost zones of dcx (ccc, 11/20/08)

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( IK, IJ, IL ) 
    do ik = k1,   k2
    do ij = julo, jhi

       do il = -i2/3, 0
          dcx(il,ij,ik) = dcx(i2+il,ij,ik)
       enddo
      
       do il = i2+1, i2+i2/3
          dcx(il,ij,ik) = dcx(il-i2,ij,ik)

       enddo
    enddo
    enddo
    !$OMP END PARALLEL DO

  END SUBROUTINE Xmist
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fxppm
!
! !DESCRIPTION: Subroutine Fxppm is the 1D "outer" flux form operator based 
!  on the Piecewise Parabolic Method (PPM; see also Lin and Rood 1996) for 
!  computing the fluxes in the E-W direction.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Fxppm( ij,  ik,  ilmt, crx, dcx, fx, qqv,     &
                    ILO, IHI, JULO, JHI, I1,  I2, K1,  K2 )
!
! !INPUT PARAMETERS: 
!
    ! Local min & max longitude (I) and altitude (K) indices
    INTEGER, INTENT(IN)    :: I1,     I2
    INTEGER, INTENT(IN)    :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)    :: ILO,    IHI
    INTEGER, INTENT(IN)    :: JULO,   JHI

    ! Latitude (IJ) and altitude (IK) indices
    INTEGER, INTENT(IN)    :: ij,     ik

    ! Controls various options in E-W advection
    INTEGER, INTENT(IN)    :: ilmt

    ! Courant number in E-W direction
    REAL*8,  INTENT(IN)    :: crx(I1:I2, JULO:JHI, K1:K2)
!
! !INPUT/OUTPUT PARAMETERS: 
!
      ! Concentration contribution from N-S advection [mixing ratio]
    REAL*8,  INTENT(INOUT) :: qqv(ILO:IHI, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Slope of concentration distribution in E-W direction (mixing ratio)
    REAL*8,  INTENT(OUT)   :: dcx(ILO:IHI, JULO:JHI, K1:K2)

    ! E-W flux [mixing ratio]
    REAL*8,  INTENT(OUT)   :: fx(I1:I2, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
!
! !REMARKS:
!   This routine is called from w/in a OpenMP parallel loop fro

! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.
!                               Also remove the allocatable arrays, which
!                               interfere w/ OpenMP parallelization.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    !-------------------------------------------------------------------------
    ! Prior to 12/5/08:
    ! Remove this (explanation below).
    !LOGICAL,             SAVE :: first = .true.
    !-------------------------------------------------------------------------

    INTEGER                   :: il
    INTEGER                   :: ilm1
    INTEGER                   :: lenx
    REAL*8                    :: r13, r23
    REAL*8                    :: rval
 
    !------------------------------------------------------------------------
    ! Prior to 12/5/08:
    ! NOTE: It is a bad idea to make these arrays allocatable.  The way this
    ! was implemented, it tried to create these arrays once for each thread.
    ! This led to a segmentation fault.  Better to just define these arrays
    ! with the appropriate dimensions.  Also note, we don't really need to
    ! use SAVE since these arrays are being reset to zero on each call
    ! to Fxppm. (bmy, 12/5/08)
    !
    !! Arrays
    !REAL*8, ALLOCATABLE, SAVE :: a6(:)
    !REAL*8, ALLOCATABLE, SAVE :: al(:)
    !REAL*8, ALLOCATABLE, SAVE :: ar(:)
    !REAL*8, ALLOCATABLE, SAVE :: a61(:)
    !REAL*8, ALLOCATABLE, SAVE :: al1(:)
    !REAL*8, ALLOCATABLE, SAVE :: ar1(:)
    !REAL*8, ALLOCATABLE, SAVE :: dcxi1(:)
    !REAL*8, ALLOCATABLE, SAVE :: qqvi1(:)
    !------------------------------------------------------------------------

    ! Arrays
    REAL*8                    :: a6( ILO:IHI )
    REAL*8                    :: al( ILO:IHI )
    REAL*8                    :: ar( ILO:IHI )
    REAL*8                    :: a61(   (IHI-1) - (ILO+1) + 1 )
    REAL*8                    :: al1(   (IHI-1) - (ILO+1) + 1 )
    REAL*8                    :: ar1(   (IHI-1) - (ILO+1) + 1 )
    REAL*8                    :: dcxi1( (IHI-1) - (ILO+1) + 1 )
    REAL*8                    :: qqvi1( (IHI-1) - (ILO+1) + 1 )

    !     ----------------
    !     Begin execution.
    !     ----------------

!------------------------------------------------------------------------------
! Prior to 12/5/08:
! Remove the ALLOCATE command, since we are now declaring these as regular
! subroutine arrays and not making them allocatable. (bmy, 12/5/08)
!!   ==========
!    if (first) then
!!   ==========
!         
!       first = .false.
!         
!       Allocate (a6(ilo:ihi))
!       Allocate (al(ilo:ihi))
!       Allocate (ar(ilo:ihi))
!       a6 = 0.0d0; al = 0.0d0; ar = 0.0d0
!         
!       Allocate (a61((ihi-1)-(ilo+1)+1))
!       Allocate (al1((ihi-1)-(ilo+1)+1))
!       Allocate (ar1((ihi-1)-(ilo+1)+1))
!       a61 = 0.0d0; al1 = 0.0d0; ar1 = 0.0d0
!         
!       Allocate (dcxi1((ihi-1)-(ilo+1)+1))
!       Allocate (qqvi1((ihi-1)-(ilo+1)+1))
!       dcxi1 = 0.0d0; qqvi1 = 0.0d0
!       
!    end if 
!------------------------------------------------------------------------------
     
    ! Zero arrays (bmy, 12/5/08)
    a6    = 0.0d0
    al    = 0.0d0
    ar    = 0.0d0
    a61   = 0.0d0
    al1   = 0.0d0
    ar1   = 0.0d0
    dcxi1 = 0.0d0
    qqvi1 = 0.0d0

    r13 = 1.0d0 / 3.0d0
    r23 = 2.0d0 / 3.0d0


    do il = ilo + 1, ihi
       
       rval = 0.5d0 * (qqv(il-1,ij,ik) + qqv(il,ij,ik)) +  &
                      (dcx(il-1,ij,ik) - dcx(il,ij,ik)) * r13

       al(il)   = rval
       ar(il-1) = rval

    end do
      

    do il = ilo + 1, ihi - 1
       a6(il) = 3.0d0 *  &
                (qqv(il,ij,ik) + qqv(il,ij,ik) - (al(il) + ar(il)))
    end do


!   ==============
    if (ilmt <= 2) then
!   ==============

       a61(:) = 0.0d0
       al1(:) = 0.0d0
       ar1(:) = 0.0d0
         
       dcxi1(:) = 0.0d0
       qqvi1(:) = 0.0d0
         
       lenx = 0
         
       do il = ilo + 1, ihi - 1
            
          lenx = lenx + 1
            
          a61(lenx)   = a6(il)
          al1(lenx)   = al(il)
          ar1(lenx)   = ar(il)
            
          dcxi1(lenx) = dcx(il,ij,ik)
          qqvi1(lenx) = qqv(il,ij,ik)
            
       end do
         
!      ===========
       call Lmtppm  &
            (lenx, ilmt, a61, al1, ar1, dcxi1, qqvi1)
!      ===========
         
       lenx = 0
         
       do il = ilo + 1, ihi - 1
            
          lenx = lenx + 1
            
          a6(il)   = a61(lenx)
          al(il)   = al1(lenx)
          ar(il)   = ar1(lenx)
            
          dcx(il,ij,ik) = dcxi1(lenx)
          qqv(il,ij,ik) = qqvi1(lenx)
            
       end do

       ! Populate ghost zones of qqv and dcx with new values (ccc, 11/20/08)
       do il = -i2/3,0
          dcx(il,ij,ik) = dcx(i2+il,ij,ik)
          qqv(il,ij,ik) = qqv(i2+il,ij,ik)
       enddo

       do il = i2+1, i2+i2/3
          dcx(il,ij,ik) = dcx(il-i2,ij,ik)
          qqv(il,ij,ik) = qqv(il-i2,ij,ik)            
       enddo
         
    end if
      
      
    do il = i1+1, i2
         
       if (crx(il,ij,ik) > 0.0d0) then
          
          ilm1 = il - 1
            
          fx(il,ij,ik) =  &
               ar(ilm1) +  &
               0.5d0 * crx(il,ij,ik) *  &
               (al(ilm1) - ar(ilm1) +  &
               (a6(ilm1) * (1.0d0 - (r23 * crx(il,ij,ik)))))
            
       else
            
          fx(il,ij,ik) =  &
               al(il) -  &
               0.5d0 * crx(il,ij,ik) *  &
               (ar(il) - al(il) +  &
               (a6(il) * (1.0d0 + (r23 * crx(il,ij,ik)))))
            
       end if
         
    end do

    ! First box case (ccc, 11/20/08)
    if (crx(i1,ij,ik) > 0.0d0) then
         
       ilm1 = i2
       
       fx(i1,ij,ik) =  &
            ar(ilm1) +  &
            0.5d0 * crx(i1,ij,ik) *  &
            (al(ilm1) - ar(ilm1) +  &
            (a6(ilm1) * (1.0d0 - (r23 * crx(i1,ij,ik)))))
         
    else
         
       fx(i1,ij,ik) =  &
            al(i1) -  &
            0.5d0 * crx(i1,ij,ik) *  &
            (ar(i1) - al(i1) +  &
            (a6(i1) * (1.0d0 + (r23 * crx(i1,ij,ik)))))
         
    end if
      
      
      
  END SUBROUTINE Fxppm
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Lmtppm
!
! !DESCRIPTION: Subroutine Lmtppm enforces the full monotonic, semi-monotonic, 
!  or the positive-definite constraint to the sub-grid parabolic distribution
!  of the Piecewise Parabolic Method (PPM).
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Lmtppm( lenx, lmt, a6, al, ar, dc, qa )
!
! !INPUT PARAMETERS: 

    ! If 0 => full monotonicity;
    ! If 1 => semi-monotonic constraint (no undershoots);
    ! If 2 => positive-definite constraint
    INTEGER, INTENT(IN)    :: lmt

    ! Vector length
    INTEGER, INTENT(IN)    :: lenx
!
! !INPUT/OUTPUT PARAMETERS: 
!
    ! Curvature of the test parabola
    REAL*8,  INTENT(INOUT) :: a6(lenx)

    ! Left edge value of the test parabola
    REAL*8,  INTENT(INOUT) :: al(lenx)

    ! Right edge value of the test parabola
    REAL*8,  INTENT(INOUT) :: ar(lenx)

    ! 0.5 * mismatch
    REAL*8,  INTENT(INOUT) :: dc(lenx)

    ! Cell-averaged value
    REAL*8,  INTENT(INOUT) :: qa(lenx)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: il
    REAL*8  :: a6da
    REAL*8  :: da1, da2
    REAL*8  :: fmin, ftmp
    REAL*8  :: r12


!   ----------------
!   Begin execution.
!   ----------------

    r12 = 1.0d0 / 12.0d0


!   =============
    if (lmt == 0) then
!   =============

!      ----------------
!      Full constraint.
!      ----------------

       do il = 1, lenx
          
          if (dc(il) == 0.0d0) then
               
             a6(il) = 0.0d0
             al(il) = qa(il)
             ar(il) = qa(il)
               
          else

             da1  = ar(il) - al(il)
             da2  = da1    * da1
             a6da = a6(il) * da1
               
             if (a6da < -da2) then
                
                a6(il) = 3.0d0 * (al(il) - qa(il))
                ar(il) = al(il) - a6(il)
                  
             else if (a6da > da2) then
                  
                a6(il) = 3.0d0 * (ar(il) - qa(il))
                al(il) = ar(il) - a6(il)
                  
             end if
               
          end if
            
       end do
         
         
!   ==================
    else if (lmt == 1) then
!   ==================

!      --------------------------
!      Semi-monotonic constraint.
!      --------------------------

       do il = 1, lenx
            
          if (Abs (ar(il) - al(il)) < -a6(il)) then
               
             if ((qa(il) < ar(il)) .and. (qa(il) < al(il))) then
                
                a6(il) = 0.0d0
                al(il) = qa(il)
                ar(il) = qa(il)
                  
             else if (ar(il) > al(il)) then
                
                a6(il) = 3.0d0 * (al(il) - qa(il))
                ar(il) = al(il) - a6(il)
                  
             else
                  
                a6(il) = 3.0d0 * (ar(il) - qa(il))
                al(il) = ar(il) - a6(il)
                  
             end if
               
          end if
            
       end do


!   ==================
    else if (lmt == 2) then
!   ==================

       do il = 1, lenx
            
          if (Abs (ar(il) - al(il)) < -a6(il)) then
               
             ftmp = ar(il) - al(il)
               
             fmin = qa(il) +  &
                    0.25d0 * (ftmp * ftmp) / a6(il) +  &
                    a6(il) * r12
               
             if (fmin < 0.0d0) then
                  
                if ((qa(il) < ar(il)) .and. (qa(il) < al(il))) then
                     
                   a6(il) = 0.0d0
                   al(il) = qa(il)
                   ar(il) = qa(il)
                     
                else if (ar(il) > al(il)) then
                     
                   a6(il) = 3.0d0 * (al(il) - qa(il))
                   ar(il) = al(il) - a6(il)
                     
                else
                     
                   a6(il) = 3.0d0 * (ar(il) - qa(il))
                   al(il) = ar(il) - a6(il)
                     
                end if
                  
             end if

          end if

       end do

    end if


  END SUBROUTINE Lmtppm
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Ytp
!
! !DESCRIPTION: Subroutine Ytp does horizontal advection in the N-S direction.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Ytp( jlmt,  geofac_pc, geofac, cry,  dq1,   qqu,    qqv,    &
                  ymass, fy,        J1P,    J2P,  I1_GL, I2_GL,  JU1_GL, &
                  J2_GL, ilong,     ILO,    IHI,  JULO,  JHI,    I1,     &
                  I2,    JU1,       J2,     K1,   K2,    jord )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)    :: J1P,    J2P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)    :: I1_GL,  I2_GL
    INTEGER, INTENT(IN)    :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)    :: I1,     I2
    INTEGER, INTENT(IN)    :: JU1,    J2
    INTEGER, INTENT(IN)    :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)    :: ILO,    IHI
    INTEGER, INTENT(IN)    :: JULO,   JHI

    ! ???
    INTEGER, INTENT(IN)    :: ilong

    ! Controls various options in N-S advection
    INTEGER, INTENT(IN)    :: jlmt

    ! N-S transport scheme (see module header for more info)
    INTEGER, INTENT(IN)    :: jord       

    ! special geometrical factor (geofac) for Polar cap
    REAL*8,  INTENT(IN)    :: geofac_pc

    ! geometrical factor for meridional advection; geofac uses correct 
    ! spherical geometry, and replaces acosp as the  meridional geometrical 
    ! factor in tpcore
    REAL*8,  INTENT(IN)    :: geofac(JU1_GL:J2_GL)

    ! Courant number in N-S direction
    REAL*8,  INTENT(IN)    :: cry(ILO:IHI, JULO:JHI, K1:K2)
    
    ! Concentration contribution from E-W advection [mixing ratio]
    REAL*8,  INTENT(IN)    :: qqu(ILO:IHI, JULO:JHI, K1:K2)

    ! Horizontal mass flux in N-S direction [hPa]
    REAL*8,  INTENT(IN)    :: ymass(ILO:IHI, JULO:JHI, K1:K2)
!
! !INPUT/OUTPUT PARAMETERS: 
!
    ! Species density [hPa]
    REAL*8,  INTENT(INOUT) :: dq1(ILO:IHI, JULO:JHI, K1:K2)

    ! Concentration contribution from N-S advection [mixing ratio]
    REAL*8,  INTENT(INOUT) :: qqv(ILO:IHI, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! N-S flux [mixing ratio]
    REAL*8,  INTENT(OUT)   :: fy(ILO:IHI, JULO:JHI+1, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: il, ij, ik
    INTEGER :: jv
    REAL*8  :: rj1p

    ! Arrays
    REAL*8  :: dcy(ilo:ihi, julo:jhi, k1:k2)


!   ----------------
!   Begin execution.
!   ----------------

    dcy(:,:,:) = 0.0d0
    fy(:,:,:) = 0.0d0

    rj1p = j1p
      

!   ==============
    if (JORD == 1) then
!   ==============

       !$OMP PARALLEL DO              &
       !$OMP DEFAULT( SHARED )        &
       !$OMP PRIVATE( IK, IJ, IL, JV )
       do ik = k1, k2
          do ij = j1p, j2p+1
             do il = i1, i2
!c?
                jv = rj1p - cry(il,ij,ik)
                  
                qqv(il,ij,ik) = qqu(il,jv,ik)
                  
             end do
          end do
       end do
       !$OMP END PARALLEL DO

!   ====
    else
!   ====

!      ==========
       call Ymist  &
!     ==========
            (4, dcy, qqu, &
             i1_gl, i2_gl, ju1_gl, j2_gl, j1p, &
             ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)


       if ((JORD <= 0) .or. (JORD >= 3)) then
            
!         ==========
          call Fyppm  &
!         ==========
               (jlmt, cry, dcy, qqu, qqv, &
                j1p, j2p, i1_gl, i2_gl, ju1_gl, j2_gl, ilong, &
                ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)
            
       else

          !$OMP PARALLEL DO               &
          !$OMP DEFAULT( SHARED )         &
          !$OMP PRIVATE( IK, IJ, IL, JV )
          do ik = k1, k2
             do ij = j1p, j2p+1
                do il = i1, i2
!c?
                   jv = rj1p - cry(il,ij,ik)
                     
                   qqv(il,ij,ik) =  &
                        qqu(il,jv,ik) +  &
                        ((Sign (1.0d0, cry(il,ij,ik)) - cry(il,ij,ik)) *  &
                        dcy(il,jv,ik))

                end do
             end do
          end do
          !$OMP END PARALLEL DO
          
       end if
         
    end if

    !-----------------------------------------------------------------------
    ! Prior to 12/5/08:
    ! We need to add an outer IK loop for OpenMP parallelization.  
    ! Preserve original code here (bmy, 12/5/08)
    !do ij = j1p, j2p+1
    !   qqv(i1:i2,ij,:) = qqv(i1:i2,ij,:) * ymass(i1:i2,ij,:)
    !end do
    !-----------------------------------------------------------------------

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( IK, IJ )
    do ik = k1,  k2
    do ij = j1p, j2p+1
       qqv(i1:i2,ij,ik) = qqv(i1:i2,ij,ik) * ymass(i1:i2,ij,ik)
    end do
    end do
    !$OMP END PARALLEL DO

    !.sds.. save N-S species flux as diagnostic
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( IK, IJ )
    do ik = k1,k2
       do ij = i1,i2
          fy(ij,j1p:j2p+1,ik) = qqv(ij,j1p:j2p+1,ik) * geofac(j1p:j2p+1)
       enddo
    enddo
    !$OMP END PARALLEL DO

    !--------------------------------------------------------------------
    ! Prior to 12/5/08:
    ! We need to add an outer IK loop for OpenMP parallelization.  
    ! Preserve original code here (bmy, 12/5/08)
    !!... meridional flux update
    !do ij = j1p, j2p
    !     
    !   dq1(i1:i2,ij,:) =  &
    !        dq1(i1:i2,ij,:) +  &
    !        (qqv(i1:i2,ij,:) - qqv(i1:i2,ij+1,:)) * geofac(ij)
    !     
    !end do
    !--------------------------------------------------------------------

    !... meridional flux update
    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( IK, IJ ) 
    do ik = k1,  k2
    do ij = j1p, j2p
         
       dq1(i1:i2,ij,ik) =  &
            dq1(i1:i2,ij,ik) +  &
            (qqv(i1:i2,ij,ik) - qqv(i1:i2,ij+1,ik)) * geofac(ij)
         
    end do
    end do
    !$OMP END PARALLEL DO

!   ====================
    call Do_Ytp_Pole_Sum  &
!   ====================
         (geofac_pc, dq1, qqv, fy, &
          i1_gl, i2_gl, ju1_gl, j2_gl, j1p, j2p, &
          ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)

  END SUBROUTINE Ytp
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ymist
!
! !DESCRIPTION: Subroutine Ymist computes the linear tracer slope in the N-S 
!  direction.  It uses the Lin et. al. 1994 algorithm.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Ymist( id,    dcy, qqu, I1_GL, I2_GL, JU1_GL, &
                    J2_GL, J1P, ILO, IHI,   JULO,  JHI,    &
                    I1,    I2,  JU1, J2,    K1,    K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude index at the edge of the South polar cap
    ! J1P=JU1_GL+1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: J1P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: I1_GL,  I2_GL
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI

    ! The "order" of the accuracy in the computed linear "slope"
    ! (or mismatch, Lin et al. 1994); it is either 2 or 4.
    INTEGER, INTENT(IN)  :: id

    ! Concentration contribution from E-W advection (mixing ratio)
    REAL*8,  INTENT(IN)  :: qqu(ILO:IHI, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Slope of concentration distribution in N-S direction [mixing ratio]
    REAL*8,  INTENT(OUT) :: dcy(ILO:IHI, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: il, ij, ik
    REAL*8  :: pmax, pmin
    REAL*8  :: r24
    REAL*8  :: tmp

    ! Arrays
    !ccc Test for indexes going further than poles : 
    ! I suppose the values for these indexes are 0. 
    ! It should work as the pole values are re-calculated in the 
    ! pole functions.
    REAL*8 :: qtmp(ilo:ihi, julo-2:jhi+2,k1:k2)

!   ----------------
!   Begin execution.
!   ----------------

    r24  = 1.0d0 / 24.0d0

    !ccc Populate qtmp
    qtmp = 0.
    do ij=ju1,j2
       qtmp(:,ij,:) = qqu(:,ij,:)
    enddo

!   ============
    if (id == 2) then
!   ============

       !$OMP PARALLEL DO                           &
       !$OMP DEFAULT( SHARED )                     &
       !$OMP PRIVATE( IK, IJ, IL, TMP, PMAX, PMIN )
       do ik = k1, k2
          do ij = ju1 - 1, j2 - 1
             do il = i1, i2

                tmp  = 0.25d0 * (qtmp(il,ij+2,ik) - qtmp(il,ij,ik))
                  
                pmax =  &
                  Max (qtmp(il,ij,ik), qtmp(il,ij+1,ik), qtmp(il,ij+2,ik)) -  &
                  qtmp(il,ij+1,ik)
                
                pmin =  &
                     qtmp(il,ij+1,ik) -  &
                     Min (qtmp(il,ij,ik), qtmp(il,ij+1,ik), qtmp(il,ij+2,ik))

                dcy(il,ij+1,ik) = Sign (Min (Abs (tmp), pmin, pmax), tmp)

             end do
          end do
       end do
       !$OMP END PARALLEL DO

!   ====
    else
!   ====

!      ========================
       call Do_Ymist_Pole1_I2d2 &
!      ========================
            (dcy, qtmp, &
             i1_gl, i2_gl, ju1_gl, j2_gl, &
             ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)

       !$OMP PARALLEL DO                           &
       !$OMP DEFAULT( SHARED )                     &
       !$OMP PRIVATE( IK, IJ, IL, TMP, PMAX, PMIN )
       do ik = k1, k2
          do ij = ju1 - 2, j2 - 2
             do il = i1, i2
                  
                tmp  = ((8.0d0 * (qtmp(il,ij+3,ik) - qtmp(il,ij+1,ik))) +  &
                         qtmp(il,ij,ik) - qtmp(il,ij+4,ik)) *  &
                         r24

                pmax =  &
                  Max (qtmp(il,ij+1,ik), qtmp(il,ij+2,ik), qtmp(il,ij+3,ik))  &
                       - qtmp(il,ij+2,ik)

                pmin =  &
                     qtmp(il,ij+2,ik) -  &
                     Min (qtmp(il,ij+1,ik), qtmp(il,ij+2,ik), qtmp(il,ij+3,ik))

                dcy(il,ij+2,ik) = Sign (Min (Abs (tmp), pmin, pmax), tmp)
                  
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    end if


!   ========================
    call Do_Ymist_Pole2_I2d2 &
!   ========================
         (dcy, qtmp, &
          i1_gl, i2_gl, ju1_gl, j2_gl, j1p, &
          ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)


  END SUBROUTINE Ymist
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Ymist_Pole1_I2d2
!
! !DESCRIPTION: Subroutine Do\_Ymist\_Pole1\_I2d2 sets "dcy" at the Poles.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Ymist_Pole1_I2d2( dcy,   qqu, I1_GL, I2_GL, JU1_GL,   &
                                  J2_GL, ILO, IHI,   JULO,  JHI,      &
                                  I1,    I2,  JU1,   J2,    K1,   K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global min & max longitude (I) and latitude (J) indices
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: I1_GL,  I2_GL
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI

    ! Concentration contribution from E-W advection [mixing ratio]
    REAL*8,  INTENT(IN)  :: qqu(ILO:IHI, JULO-2:JHI+2, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Slope of concentration distribution in N-S direction [mixing ratio]
    REAL*8, INTENT(OUT)  :: dcy(ILO:IHI, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    INTEGER :: i2d2
    INTEGER :: il, ik
    REAL*8  :: pmax, pmin
    REAL*8  :: r24
    REAL*8  :: tmp
      

!   ----------------
!   Begin execution.
!   ----------------

    i2d2 = i2_gl / 2
      
    r24  = 1.0d0 / 24.0d0


!   ==================
    if (ju1 == ju1_gl) then
!   ==================

       !$OMP PARALLEL DO                       &
       !$OMP DEFAULT( SHARED )                 &
       !$OMP PRIVATE( IK, IL, TMP, PMAX, PMIN ) 
       do ik = k1, k2
          do il = i1, i2d2
               
             tmp  =  &
                  ((8.0d0 * (qqu(il,ju1+2,ik) - qqu(il,ju1,ik))) +  &
                  qqu(il+i2d2,ju1+1,ik) - qqu(il,ju1+3,ik)) *  &
                  r24
               
             pmax = Max (qqu(il,ju1,ik), qqu(il,ju1+1,ik),  &
                  qqu(il,ju1+2,ik)) -  &
                  qqu(il,ju1+1,ik)
             
             pmin = qqu(il,ju1+1,ik) -  &
                  Min (qqu(il,ju1,ik), qqu(il,ju1+1,ik),  &
                  qqu(il,ju1+2,ik))
               
             dcy(il,ju1+1,ik) =  &
                  Sign (Min (Abs (tmp), pmin, pmax), tmp)
               
          end do
       end do
       !$OMP END PARALLEL DO
        
       !$OMP PARALLEL DO                       &
       !$OMP DEFAULT( SHARED )                 &
       !$OMP PRIVATE( IK, IL, TMP, PMAX, PMIN ) 
        do ik = k1, k2
          do il = i1 + i2d2, i2
             
             tmp  =  &
                  ((8.0d0 * (qqu(il,ju1+2,ik) - qqu(il,ju1,ik))) +  &
                  qqu(il-i2d2,ju1+1,ik) - qqu(il,ju1+3,ik)) *  &
                  r24
             
             pmax = Max (qqu(il,ju1,ik), qqu(il,ju1+1,ik),  &
                  qqu(il,ju1+2,ik)) -  &
                  qqu(il,ju1+1,ik)
             
             pmin = qqu(il,ju1+1,ik) -  &
                  Min (qqu(il,ju1,ik), qqu(il,ju1+1,ik),  &
                  qqu(il,ju1+2,ik))
             
             dcy(il,ju1+1,ik) =  &
                  Sign (Min (Abs (tmp), pmin, pmax), tmp)
               
          end do
       end do
       !$OMP END PARALLEL DO

!   ======
    end if
!   ======


!   ================
    if (j2 == j2_gl) then
!   ================

       !$OMP PARALLEL DO                       &
       !$OMP DEFAULT( SHARED )                 &
       !$OMP PRIVATE( IK, IL, TMP, PMAX, PMIN ) 
       do ik = k1, k2
          do il = i1, i2d2
               
             tmp  =  &
                  ((8.0d0 * (qqu(il,j2,ik) - qqu(il,j2-2,ik))) +  &
                  qqu(il,j2-3,ik) - qqu(il+i2d2,j2-1,ik)) *  &
                  r24
               
             pmax = Max (qqu(il,j2-2,ik), qqu(il,j2-1,ik),  &
                  qqu(il,j2,ik)) -  &
                  qqu(il,j2-1,ik)
               
             pmin = qqu(il,j2-1,ik) -  &
                  Min (qqu(il,j2-2,ik), qqu(il,j2-1,ik),  &
                  qqu(il,j2,ik))
               
             dcy(il,j2-1,ik) =  &
                  Sign (Min (Abs (tmp), pmin, pmax), tmp)
               
          end do
       end do
       !$OMP END PARALLEL DO

       !$OMP PARALLEL DO                       &
       !$OMP DEFAULT( SHARED )                 &
       !$OMP PRIVATE( IK, IL, TMP, PMAX, PMIN ) 
       do ik = k1, k2
          do il = i1 + i2d2, i2
               
             tmp  =  &
                  ((8.0d0 * (qqu(il,j2,ik) - qqu(il,j2-2,ik))) +  &
                  qqu(il,j2-3,ik) - qqu(il-i2d2,j2-1,ik)) *  &
                  r24
             
             pmax = Max (qqu(il,j2-2,ik), qqu(il,j2-1,ik),  &
                  qqu(il,j2,ik)) -  &
                  qqu(il,j2-1,ik)
               
             pmin = qqu(il,j2-1,ik) -  &
                  Min (qqu(il,j2-2,ik), qqu(il,j2-1,ik),  &
                  qqu(il,j2,ik))
               
             dcy(il,j2-1,ik) =  &
                  Sign (Min (Abs (tmp), pmin, pmax), tmp)
             
          end do
       end do
       !$OMP END PARALLEL DO

!   ======
    end if
!   ======


  END SUBROUTINE Do_Ymist_Pole1_I2d2
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Ymist_Pole2_I2d2
!
! !DESCRIPTION: Subroutine Do\_Ymist\_Pole2\_I2d2 sets "dcy" at the Poles.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Ymist_Pole2_I2d2( dcy,   qqu, I1_GL, I2_GL, JU1_GL, &
                                  J2_GL, J1P, ILO,   IHI,   JULO,   &
                                  JHI,   I1,  I2,    JU1,   J2,     &
                                  K1,    K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude index at the edge of the South polar cap
    ! J1P=JU1_GL+1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: J1P 

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: I1_GL,  I2_GL
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI

    ! Concentration contribution from E-W advection [mixing ratio]
    REAL*8,  INTENT(IN)  :: qqu(ILO:IHI, JULO-2:JHI+2, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Slope of concentration distribution in N-S direction [mixing ratio]
    REAL*8,  INTENT(OUT) :: dcy(ILO:IHI, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: i2d2
    INTEGER :: il, ik
    REAL*8  :: pmax, pmin
    REAL*8  :: tmp

!   ----------------
!   Begin execution.
!   ----------------

    i2d2 = i2_gl / 2


!   ==================
    if (ju1 == ju1_gl) then
!   ==================
         
       if (j1p /= ju1_gl+1) then
            
          dcy(i1:i2,ju1,:) = 0.0d0
            
       else

!         -----------------------------------------------
!         Determine slope in South Polar cap for scalars.
!         -----------------------------------------------

          !$OMP PARALLEL DO                        &
          !$OMP DEFAULT( SHARED )                  &
          !$OMP PRIVATE( IK, IL, TMP, PMAX, PMIN )
          do ik = k1, k2
             do il = i1, i2d2
                  
                tmp  =  &
                     0.25d0 *  &
                     (qqu(il,ju1+1,ik) - qqu(il+i2d2,ju1+1,ik))

                pmax =  &
                     Max (qqu(il,ju1+1,ik), qqu(il,ju1,ik),  &
                     qqu(il+i2d2,ju1+1,ik)) -  &
                     qqu(il,ju1,ik)
                
                pmin =  &
                     qqu(il,ju1,ik) -  &
                     Min (qqu(il,ju1+1,ik), qqu(il,ju1,ik),  &
                     qqu(il+i2d2,ju1+1,ik))

                dcy(il,ju1,ik) =  &
                     Sign (Min (Abs (tmp), pmax, pmin), tmp)
                
             end do
          end do
          !$OMP END PARALLEL DO

          !----------------------------------------------------------------
          ! Prior to 12/5/08:
          ! We need to add and outer IK loop for OpenMP parallelization.
          ! Preserve original code here (bmy, 12/5/08)
          !do il = i1 + i2d2, i2
          !   dcy(il,ju1,:) = -dcy(il-i2d2,ju1,:)
          !end do
          !----------------------------------------------------------------
          
          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( IK, IL )
          do ik = k1,        k2
          do il = i1 + i2d2, i2
             dcy(il,ju1,ik) = -dcy(il-i2d2,ju1,ik)
          end do
          end do
          !$OMP END PARALLEL DO

       end if

!   ======
    end if
!   ======


!   ================
    if (j2 == j2_gl) then
!   ================

       if (j1p /= ju1_gl+1) then

          dcy(i1:i2,j2,:) = 0.0d0
            
       else

!         -----------------------------------------------
!         Determine slope in North Polar cap for scalars.
!         -----------------------------------------------

          !$OMP PARALLEL DO                        &
          !$OMP DEFAULT( SHARED )                  &
          !$OMP PRIVATE( IK, IL, TMP, PMAX, PMIN )
          do ik = k1, k2
             do il = i1, i2d2
                 
                tmp  =  &
                     0.25d0 *  &
                     (qqu(il+i2d2,j2-1,ik) - qqu(il,j2-1,ik))
                
                pmax =  &
                     Max (qqu(il+i2d2,j2-1,ik), qqu(il,j2,ik),  &
                     qqu(il,j2-1,ik)) -  &
                     qqu(il,j2,ik)

                pmin =  &
                     qqu(il,j2,ik) -  &
                     Min (qqu(il+i2d2,j2-1,ik), qqu(il,j2,ik),  &
                     qqu(il,j2-1,ik))

                dcy(il,j2,ik) =  &
                     Sign (Min (Abs (tmp), pmax, pmin), tmp)

             end do
          end do
          !$OMP END PARALLEL DO

          !----------------------------------------------------------------
          ! Prior to 12/5/08:
          ! We need to add and outer IK loop for OpenMP parallelization.
          ! Preserve original code here (bmy, 12/5/08)
          !do il = i1 + i2d2, i2
          !   dcy(il,j2,:) = -dcy(il-i2d2,j2,:)
          !end do
          !----------------------------------------------------------------

          !$OMP PARALLEL DO       &
          !$OMP DEFAULT( SHARED ) &
          !$OMP PRIVATE( IK, IL )!
          do ik = k1,        k2
          do il = i1 + i2d2, i2
             dcy(il,j2,ik) = -dcy(il-i2d2,j2,ik)
          end do
          end do
          !$OMP END PARALLEL DO

       end if

!   ======
    end if
!   ======


  END SUBROUTINE Do_Ymist_Pole2_I2d2
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fyppm
!
! !DESCRIPTION: Subroutine Fyppm is the 1D "outer" flux form operator based 
!  on the Piecewise Parabolic Method (PPM; see also Lin and Rood 1996) for 
!  computing the fluxes in the N-S direction.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Fyppm( jlmt,  cry,   dcy,    qqu,   qqv,   j1p, j2p,    &
                    i1_gl, i2_gl, ju1_gl, j2_gl, ilong, ilo, ihi,    &
                    julo,  jhi,   i1,     i2,    ju1,   j2,  k1, k2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)  :: J1P,    J2P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: I1_GL,  I2_GL
    INTEGER, INTENT(IN)  :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)  :: I1,     I2
    INTEGER, INTENT(IN)  :: JU1,    J2
    INTEGER, INTENT(IN)  :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)  :: ILO,    IHI
    INTEGER, INTENT(IN)  :: JULO,   JHI
    
    ! ILONG ??
    INTEGER, INTENT(IN)  :: ilong

    ! Controls various options in N-S advection
    INTEGER, INTENT(IN)  :: jlmt

    ! Courant number in N-S direction
    REAL*8,  INTENT(IN)  :: cry(ILO:IHI, JULO:JHI, K1:K2)

    ! Slope of concentration distribution in N-S direction [mixing ratio]
    REAL*8,  INTENT(IN)  :: dcy(ILO:IHI, JULO:JHI, K1:K2)

    ! Concentration contribution from E-W advection [mixing ratio]
    REAL*8,  INTENT(IN)  :: qqu(ILO:IHI, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Concentration contribution from N-S advection [mixing ratio]
    REAL*8,  INTENT(OUT) :: qqv(ILO:IHI, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: ijm1
    INTEGER :: il, ij, ik
    INTEGER :: lenx
    REAL*8  :: r13, r23

    ! Arrays
    REAL*8  :: a61 (ilong*((JHI-1)-(JULO+1)+1))
    REAL*8  :: al1 (ilong*((JHI-1)-(JULO+1)+1))
    REAL*8  :: ar1 (ilong*((JHI-1)-(JULO+1)+1))
    REAL*8  :: dcy1(ilong*((JHI-1)-(JULO+1)+1))
    REAL*8  :: qqu1(ilong*((JHI-1)-(JULO+1)+1))
    REAL*8  :: a6(ILO:IHI, JULO:JHI, K1:K2)
    REAL*8  :: al(ILO:IHI, JULO:JHI, K1:K2)
    REAL*8  :: ar(ILO:IHI, JULO:JHI, K1:K2)
   
    ! NOTE: The code was writtein with I1:I2 as the first dimension of AL,
    ! AR, A6, AL1, A61, AR1.  However, the limits should really should be 
    ! ILO:IHI.  In practice, however, for a global grid (and OpenMP 
    ! parallelization) ILO=I1 and IHI=I2.  Nevertheless, we will change the 
    ! limits to ILO:IHI to be consistent and to avoid future problems. 
    ! (bmy, 12/5/08)

!   ----------------
!   Begin execution.
!   ----------------

    a6(:,:,:) = 0.0d0; al(:,:,:) = 0.0d0; ar(:,:,:) = 0.0d0


    r13 = 1.0d0 / 3.0d0
    r23 = 2.0d0 / 3.0d0

    !-----------------------------------------------------------------------
    ! Prior to 12/5/08:
    ! We need to add IK and IL loops for OpenMP parallelization.
    ! Preserve original code here (bmy, 12/5/08)
    !do ij = julo + 1, jhi
    !   al(i1:i2,ij,:) =  &
    !        0.5d0 * (qqu(i1:i2,ij-1,:) + qqu(i1:i2,ij,:)) +  &
    !        (dcy(i1:i2,ij-1,:) - dcy(i1:i2,ij,:)) * r13
    !end do
    !-----------------------------------------------------------------------

    !$OMP PARALLEL DO       &
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( IK, IJ )
    do ik = k1,     k2
    do ij = julo+1, jhi
    do il = ilo,    ihi
       al(il,ij,ik) =  &
            0.5d0 * (qqu(il,ij-1,ik) + qqu(il,ij,ik)) +  &
            (dcy(il,ij-1,ik) - dcy(il,ij,ik)) * r13
    end do
    end do
    end do
    !$OMP END PARALLEL DO

    !-------------------------------------------------------------------------
    ! Prior to 12/5/08:
    ! We need to add IK and IL loops for OpenMP parallelization.
    ! Preserve original code here (bmy, 12/5/08)
    ! NOTE: This DO loop doesn't parallelize, so leave it alone (bmy, 12/5/08)
    !do ij = julo, jhi - 1
    !   ar(i1:i2,ij,:) = al(i1:i2,ij+1,:)
    !end do
    !-------------------------------------------------------------------------

    !$OMP PARALLEL DO           &
    !$OMP DEFAULT( SHARED     ) &
    !$OMP PRIVATE( IK, IJ, IL )
    do ik = k1,   k2
    do ij = julo, jhi-1
    do il = ilo,  ihi
       ar(il,ij,ik) = al(il,ij+1,ik)
    end do
    end do
    end do
    !$OMP END PARALLEL DO

!   =======================
    call Do_Fyppm_Pole_I2d2 &
!   =======================
         (al, ar, &
          i1_gl, i2_gl, ju1_gl, j2_gl, &
          ilo, ihi, julo, jhi, i1, i2, ju1, j2, k1, k2)

    !-----------------------------------------------------------------------
    ! Prior to 12/5/08:
    ! We need to add IK and IL loops for OpenMP parallelization.
    ! Preserve original code here (bmy, 12/5/08)
    !do ij = julo + 1, jhi - 1
    !
    !   a6(i1:i2,ij,:) =  &
    !        3.0d0 *  &
    !        (qqu(i1:i2,ij,:) + qqu(i1:i2,ij,:) -  &
    !        (al(i1:i2,ij,:) + ar(i1:i2,ij,:)))
    !   
    !end do
    !-----------------------------------------------------------------------

    !$OMP PARALLEL DO       &  
    !$OMP DEFAULT( SHARED ) &
    !$OMP PRIVATE( IK, IJ )
    do ik = k1,     k2
    do ij = julo+1, jhi-1
    do il = ilo,    ihi

       a6(il,ij,ik) =  &
            3.0d0 *  &
            (qqu(il,ij,ik) + qqu(il,ij,ik) -  &
            (al(il,ij,ik) + ar(il,ij,ik)))
       
    end do
    end do
    end do
    !$OMP END PARALLEL DO

!   ==============
    if (jlmt <= 2) then
!   ==============
       
       !$OMP PARALLEL DO                                           &
       !$OMP DEFAULT( SHARED )                                     &
       !$OMP PRIVATE( IK, LENX, IJ, IL, A61, AL1, AR1, DCY1, QQU1 )
       do ik = k1, k2

          lenx = 0

          do ij = julo + 1, jhi - 1
             !=== Prior to 12/5/08
             !do il = i1, i2
             do il = ilo, ihi

                lenx = lenx + 1

                a61 (lenx) = a6 (il,ij,ik)
                al1 (lenx) = al (il,ij,ik)
                ar1 (lenx) = ar (il,ij,ik)
                dcy1(lenx) = dcy(il,ij,ik)
                qqu1(lenx) = qqu(il,ij,ik)

             end do
          end do

!         ===========
          call Lmtppm  &
               (lenx, jlmt, a61, al1, ar1, dcy1, qqu1)
!         ===========

          lenx = 0

          do ij = julo + 1, jhi - 1
             !=== Prior to 12/5/08
             !do il = i1, i2
             do il = ilo, ihi

                lenx = lenx + 1
                  
                a6(il,ij,ik) = a61(lenx)
                al(il,ij,ik) = al1(lenx)
                ar(il,ij,ik) = ar1(lenx)

             end do
          end do

       end do
       !$OMP END PARALLEL DO

    end if


    !$OMP PARALLEL DO                &
    !$OMP DEFAULT( SHARED )          &
    !$OMP PRIVATE( IK, IJ, IJM1, IL )
    do ik = k1, k2
       do ij = j1p, j2p+1

          ijm1 = ij - 1

          !=== Prior to 12/5/08
          !do il = i1, i2
          do il = ilo, ihi

             if (cry(il,ij,ik) > 0.0d0) then

                qqv(il,ij,ik) =  &
                     ar(il,ijm1,ik) +  &
                     0.5d0 * cry(il,ij,ik) *  &
                     (al(il,ijm1,ik) - ar(il,ijm1,ik) +  &
                     (a6(il,ijm1,ik) * (1.0d0 - (r23 * cry(il,ij,ik)))))

             else

                qqv(il,ij,ik) =  &
                     al(il,ij,ik) -  &
                     0.5d0 * cry(il,ij,ik) *  &
                     (ar(il,ij,ik) - al(il,ij,ik) +  &
                     (a6(il,ij,ik) * (1.0d0 + (r23 * cry(il,ij,ik)))))

             end if

          end do

       end do
    end do
    !$OMP END PARALLEL DO


  END SUBROUTINE Fyppm
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Fyppm_Pole_I2d2
!
! !DESCRIPTION: Subroutine Do\_Fyppm\_Pole\_I2d2  sets "al" \& "ar" at 
!  the Poles.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Fyppm_Pole_I2d2( al,  ar,  I1_GL, I2_GL, JU1_GL, J2_GL, &
                                 ILO, IHI, JULO,  JHI,   I1,     I2,    &
                                 JU1, J2,  K1,    K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)    :: I1_GL,  I2_GL
    INTEGER, INTENT(IN)    :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)    :: I1,     I2
    INTEGER, INTENT(IN)    :: JU1,    J2
    INTEGER, INTENT(IN)    :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)    :: ILO,    IHI
    INTEGER, INTENT(IN)    :: JULO,   JHI
!
! !OUTPUT PARAMETERS:
!
    ! Left (al) and right (ar) edge values of the test parabola
    REAL*8,  INTENT(INOUT) :: al(ILO:IHI, JULO:JHI, K1:K2)
    REAL*8,  INTENT(INOUT) :: ar(ILO:IHI, JULO:JHI, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: i2d2
    INTEGER :: il, ik


!   ----------------
!   Begin execution.
!   ----------------

    i2d2 = i2_gl / 2


!   ==================
    if (ju1 == ju1_gl) then
!   ==================

       !---------------------------------------------------------
       ! Prior to 12/5/08:
       ! We need to add an IK loop for OpenMP parallelization.
       ! Preserve original code here. (bmy, 12/5/08)
       !do il = i1, i2d2
       !   al(il,     ju1,:) = al(il+i2d2,ju1+1,:)
       !   al(il+i2d2,ju1,:) = al(il,     ju1+1,:)
       !end do
       !---------------------------------------------------

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IL ) 
       do ik = k1, k2
       do il = i1, i2d2
          al(il,     ju1,ik) = al(il+i2d2,ju1+1,ik)
          al(il+i2d2,ju1,ik) = al(il,     ju1+1,ik)
       end do
       end do
       !$OMP END PARALLEL DO

!   ======
    end if
!   ======


!   ================
    if (j2 == j2_gl) then
!   ================

       !---------------------------------------------------------
       ! Prior to 12/5/08:
       ! We need to add an IK loop for OpenMP parallelization.
       ! Preserve original code here. (bmy, 12/5/08)
       !do il = i1, i2d2
       !   ar(il,     j2,:)  = ar(il+i2d2,j2-1,:)
       !   ar(il+i2d2,j2,:)  = ar(il,     j2-1,:)
       !end do
       !---------------------------------------------------------

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IL ) 
       do ik = k1, k2
       do il = i1, i2d2
          ar(il,     j2,ik)  = ar(il+i2d2,j2-1,ik)
          ar(il+i2d2,j2,ik)  = ar(il,     j2-1,ik)
       end do
       end do
      !$OMP END PARALLEL DO

!   ======
    end if
!   ======


  END SUBROUTINE Do_Fyppm_Pole_I2d2
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Do_Ytp_Pole_Sum
!
! !DESCRIPTION: Subroutine Do\_Ytp\_Pole\_Sum sets "dq1" at the Poles.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Do_Ytp_Pole_Sum( geofac_pc, dq1,    qqv,   fy,  I1_GL,  &
                              I2_GL,     JU1_GL, J2_GL, J1P, J2P,    &
                              ILO,       IHI,    JULO,  JHI, I1,     &
                              I2,        JU1,    J2,    K1,  K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude indices at the edges of the S/N polar caps
    ! J1P=JU1_GL+1; J2P=J2_GL-1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2; J2P=J2_GL-2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)    :: J1P,    J2P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)    :: I1_GL,  I2_GL
    INTEGER, INTENT(IN)    :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)    :: I1,     I2
    INTEGER, INTENT(IN)    :: JU1,    J2
    INTEGER, INTENT(IN)    :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)    :: ILO,    IHI
    INTEGER, INTENT(IN)    :: JULO,   JHI

    ! Special geometrical factor (geofac) for Polar cap
    REAL*8,  INTENT(IN)    :: geofac_pc

    ! Concentration contribution from N-S advection [mixing ratio]
    REAL*8,  INTENT(IN)    :: qqv(ILO:IHI, JULO:JHI, K1:K2)
!
! !INPUT/OUTPUT PARAMETERS: 
!
    ! Species density [hPa]
    REAL*8,  INTENT(INOUT) :: dq1(ILO:IHI, JULO:JHI, K1:K2)

    ! N-S mass flux [mixing ratio]
    REAL*8,  INTENT(INOUT) :: fy (ILO:IHI, JULO:JHI+1, K1:K2)
!
! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.  Added
!                               OpenMP parallel DO loops.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: il, ik
    REAL*8  :: ri2

    ! Arrays
    REAL*8  :: dq_np(K1:K2)
    REAL*8  :: dq_sp(K1:K2)
    REAL*8  :: dqik (K1:K2)
    REAL*8  :: sumnp(K1:K2)
    REAL*8  :: sumsp(K1:K2)


!   ----------------
!   Begin execution.
!   ----------------

    ri2 = i2_gl

    dqik(:) = 0.0d0

!... Do south pole
!   ==================
    if (ju1 == ju1_gl) then
!   ==================

!... Integrate N-S flux around polar cap lat circle for each level

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IL )
       do ik = k1, k2
          sumsp(ik) = 0.0d0
            
          do il = i1, i2
             sumsp(ik) = sumsp(ik) + qqv(il,j1p,ik)
          enddo
            
       enddo
       !$OMP END PARALLEL DO

!... wrap in E-W direction
       if (i1 == i1_gl) then
          dqik(:) = dq1(i1,ju1,:)
       endif

!... normalize and set inside polar cap

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IL )
       do ik = k1, k2

          dq_sp(ik) = dqik(ik) - (sumsp(ik) / ri2 * geofac_pc)

          do il = i1, i2
             dq1(il,ju1,ik) = dq_sp(ik)
!... save polar flux
             fy(il,ju1,ik) = - (sumsp(ik) / ri2 * geofac_pc)
          enddo

          if (j1p /= ju1_gl+1) then
             do il = i1, i2
                dq1(il,ju1+1,ik) = dq_sp(ik)
!... save polar flux
                fy(il,ju1+1,ik) = - (sumsp(ik) / ri2 * geofac_pc)
             enddo

          endif

       enddo
       !$OMP END PARALLEL DO

!   ======
    endif
!   ======


!... Do north pole
!   ================
    if (j2 == j2_gl) then
!   ================

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IL )
       do ik = k1, k2
          sumnp(ik) = 0.0d0

          do il = i1, i2
             sumnp(ik) = sumnp(ik) + qqv(il,j2p+1,ik)
          enddo

       enddo
       !$OMP END PARALLEL DO

       if (i1 == i1_gl) then
          dqik(:) = dq1(i1,j2,:)
       endif

       !$OMP PARALLEL DO       &
       !$OMP DEFAULT( SHARED ) &
       !$OMP PRIVATE( IK, IL )
       do ik = k1, k2

          dq_np(ik) = dqik(ik) + (sumnp(ik) / ri2 * geofac_pc)

          do il = i1, i2
             dq1(il,j2,ik) = dq_np(ik)
!... save polar flux
             fy(il,j2+1,ik) = (sumnp(ik) / ri2* geofac_pc)
          enddo

          if (j1p /= ju1_gl+1) then

             do il = i1, i2
                dq1(il,j2-1,ik) = dq_np(ik)
!... save polar flux
                fy(il,j2,ik) = (sumnp(ik) / ri2* geofac_pc)
             enddo

          endif

       enddo
       !$OMP END PARALLEL DO

!   ======
    endif
!   ======


  END SUBROUTINE Do_Ytp_Pole_Sum
!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Fzppm
!
! !DESCRIPTION: Subroutine Fzppm is the 1D "outer" flux form operator based 
!  on the Piecewise Parabolic Method (PPM; see also Lin and Rood 1996) for 
!  computing the fluxes in the vertical direction.
!\\
!\\
!  Fzppm was modified by S.-J. Lin, 12/14/98, to allow the use of the KORD=7
!  (klmt=4) option.  KORD=7 enforces the 2nd monotonicity constraint of
!  Huynh (1996).  Note that in Huynh's original scheme, two constraints are
!  necessary for the preservation of monotonicity.  To use Huynh's
!  algorithm, it was modified as follows.  The original PPM is still used to
!  obtain the first guesses for the cell edges, and as such Huynh's 1st
!  constraint is no longer needed.  Huynh's median function is also replaced
!  by a simpler yet functionally equivalent in-line algorithm.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Fzppm( klmt,  delp1,  dpi,   wz,  dq1, qq1,  fz,      &
                    J1P,   JU1_GL, J2_GL, ILO, IHI, JULO, JHI,     &
                    ILONG, IVERT,  I1,    I2,  JU1, J2,   K1,  K2 )
!
! !INPUT PARAMETERS: 
!
    ! Global latitude index at the edges of the South polar cap
    ! J1P=JU1_GL+1 for a polar cap of 1 latitude band
    ! J1P=JU1_GL+2 for a polar cap of 2 latitude bands
    INTEGER, INTENT(IN)    :: J1P

    ! Global min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)    :: JU1_GL, J2_GL

    ! Local min & max longitude (I), latitude (J), altitude (K) indices
    INTEGER, INTENT(IN)    :: I1,     I2
    INTEGER, INTENT(IN)    :: JU1,    J2
    INTEGER, INTENT(IN)    :: K1,     K2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)    :: ILO,    IHI
    INTEGER, INTENT(IN)    :: JULO,   JHI

    ! Dimensions in longitude & altitude ???
    INTEGER, INTENT(IN)    :: ilong,  ivert

    ! Controls various options in vertical advection
    INTEGER, INTENT(IN)    :: klmt

    ! Pressure thickness, the pseudo-density in a 
    ! hydrostatic system at t1 [hPa]
    REAL*8,  INTENT(IN)    :: delp1(ILO:IHI, JULO:JHI, K1:K2)

    ! Large scale mass flux (per time step tdt) in the vertical
    ! direction as diagnosed from the hydrostatic relationship [hPa]
    REAL*8,  INTENT(IN)    :: wz(I1:I2, JU1:J2, K1:K2)

    ! Species concentration [mixing ratio]
    REAL*8,  INTENT(IN)    :: qq1(ILO:IHI, JULO:JHI, K1:K2)
!
! !INPUT/OUTPUT PARAMETERS: 
!
    ! Species density [hPa]
    REAL*8,  INTENT(INOUT) :: dq1(ILO:IHI, JULO:JHI, K1:K2)
!
! !OUTPUT PARAMETERS:
!
    ! Work array 
    REAL*8,  INTENT(OUT)   :: dpi(I1:I2, JU1:J2, K1:K2)

    ! Vertical flux [mixing ratio]
    REAL*8,  INTENT(OUT)   :: fz(ILO:IHI, JULO:JHI,  K1:K2)

! !AUTHOR:
!   Original code from Shian-Jiann Lin, DAO 
!   John Tannahill, LLNL (jrt@llnl.gov)
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: il, ij, ik
    INTEGER :: k1p1, k1p2
    INTEGER :: k2m1, k2m2
    INTEGER :: lenx
    REAL*8  :: a1, a2
    REAL*8  :: aa, bb
    REAL*8  :: c0, c1, c2
    REAL*8  :: cm, cp
    REAL*8  :: fac1, fac2
    REAL*8  :: lac
    REAL*8  :: qmax, qmin
    REAL*8  :: qmp
    REAL*8  :: r13, r23
    REAL*8  :: tmp

    ! Arrays
    REAL*8  :: a61  (ilong*(ivert-4))
    REAL*8  :: al1  (ilong*(ivert-4))
    REAL*8  :: ar1  (ilong*(ivert-4))
    REAL*8  :: dca1 (ilong*(ivert-4))
    REAL*8  :: qq1a1(ilong*(ivert-4))
    REAL*8  :: a6   (i1:i2, k1:k2)
    REAL*8  :: al   (i1:i2, k1:k2)
    REAL*8  :: ar   (i1:i2, k1:k2)
    REAL*8  :: dca  (i1:i2, k1:k2)
    REAL*8  :: dlp1a(i1:i2, k1:k2)
    REAL*8  :: qq1a (i1:i2, k1:k2)
    REAL*8  :: wza  (i1:i2, k1:k2)
    REAL*8  :: dc   (i1:i2, ju1:j2, k1:k2)


!     ----------------
!     Begin execution.
!     ----------------

    a6(:,:) = 0.0d0
    al(:,:) = 0.0d0
    ar(:,:) = 0.0d0
    dc(:,:,:) = 0.0d0
!.sds... diagnostic vertical flux for species - set top to 0.0
    fz(:,:,:) = 0.0


    k1p1 = k1 + 1
    k1p2 = k1 + 2
      
    k2m1 = k2 - 1
    k2m2 = k2 - 2

    r13  = 1.0d0 / 3.0d0
    r23  = 2.0d0 / 3.0d0


!   -------------------
!   Compute dc for PPM.
!   -------------------

    do ik = k1, k2m1
       dpi(:,:,ik) = qq1(i1:i2,ju1:j2,ik+1) - qq1(i1:i2,ju1:j2,ik)
    end do


    do ik = k1p1, k2m1

       ! NOTE: We can't parallelize over IK because there are IK-1 and IK+1 
       ! terms in the loop below.  Parallelize over IJ instead. (bmy, 12/5/08)

       !$OMP PARALLEL DO                                   &
       !$OMP DEFAULT( SHARED )                             &
       !$OMP PRIVATE( IJ, IL, C0, C1, C2, TMP, QMAX, QMIN )
       do ij = ju1, j2
          do il = i1, i2

             c0 = delp1(il,ij,ik) /  &
                  (delp1(il,ij,ik-1) + delp1(il,ij,ik) +  &
                  delp1(il,ij,ik+1))

             c1 = (delp1(il,ij,ik-1) + (0.5d0 * delp1(il,ij,ik))) /  &
                  (delp1(il,ij,ik+1) + delp1(il,ij,ik))
               
             c2 = (delp1(il,ij,ik+1) + (0.5d0 * delp1(il,ij,ik))) /  &
                  (delp1(il,ij,ik-1) + delp1(il,ij,ik))
               
             tmp = c0 *  &
                   ((c1 * dpi(il,ij,ik)) +  &
                   (c2 * dpi(il,ij,ik-1)))
               
             qmax =  &
                  Max (qq1(il,ij,ik-1),  &
                  qq1(il,ij,ik),  &
                  qq1(il,ij,ik+1)) -  &
                  qq1(il,ij,ik)
             
             qmin =  &
                  qq1(il,ij,ik) -  &
                  Min (qq1(il,ij,ik-1),  &
                       qq1(il,ij,ik),  &
                       qq1(il,ij,ik+1))

             dc(il,ij,ik) = Sign (Min (Abs (tmp), qmax, qmin), tmp)

          end do
       end do
       !$OMP END PARALLEL DO

    end do


!c?
!   -------------------------------------
!   Loop over latitudes (to save memory).
!   -------------------------------------

!   =======================
    !$OMP PARALLEL DO                                                 &
    !$OMP DEFAULT( SHARED )                                           &
    !$OMP PRIVATE( IL,   IJ,   IK,    DCA,  WZA,  DLP1A, QQ1A, FAC1 ) &
    !$OMP PRIVATE( FAC2, AA,   BB,    AL,   AR,   C1,    C2,   A1   ) &
    !$OMP PRIVATE( A2,   LENX, A6,    QMP,  LAC,  QMIN,  QMAX, AL1  ) &
    !$OMP PRIVATE( AR1,  DCA1, QQ1A1, A61,  CM,   CP                ) 
    ijloop: do ij = ju1, j2
!   =======================

       if (((ij == ju1_gl+1) .or. (ij == j2_gl-1)) .and.  &
            (j1p /= ju1_gl+1)) then
!         ============
          cycle ijloop
!         ============
       end if


       !----------------------------------------------------------------
       ! Prior to 12/5/08:
       ! Replace these with explicit loops below to facilitate
       ! OpenMP parallelization.  Preserve original code here.
       ! (bmy, 12/5/08)
       !
       !dca(:,:) = dc(:,ij,:)  ! the monotone slope
       !wza(:,:) = wz(:,ij,:)
       !
       !dlp1a(:,:) = delp1(i1:i2,ij,:)
       !qq1a (:,:) = qq1  (i1:i2,ij,:)
       !----------------------------------------------------------------

       do ik = k1, k2
       do il = i1, i2
          dca  (il,ik) = dc   (il,ij,ik)  ! the monotone slope
          wza  (il,ik) = wz   (il,ij,ik)
          dlp1a(il,ik) = delp1(il,ij,ik)
          qq1a (il,ik) = qq1  (il,ij,ik)
       enddo
       enddo

!     ----------------------------------------------------------------
!     Compute first guesses at cell interfaces.  First guesses are
!     required to be continuous.  Three-cell parabolic subgrid
!     distribution at model top; two-cell parabolic with zero gradient
!     subgrid distribution at the surface.
!     ----------------------------------------------------------------

!     ---------------------------
!     First guess top edge value.
!     ---------------------------

       do il = i1, i2

!         ------------------------------------------------------------
!         Three-cell PPM; compute a, b, & c of q = aP^2 + bP + c using
!         cell averages and dlp1a.
!         ------------------------------------------------------------

          fac1 = dpi(il,ij,k1p1) -  &
                 dpi(il,ij,k1) * (dlp1a(il,k1p1) + dlp1a(il,k1p2)) /  &
                 (dlp1a(il,k1) + dlp1a(il,k1p1))

          fac2 = (dlp1a(il,k1p1) + dlp1a(il,k1p2)) *  &
                 (dlp1a(il,k1) + dlp1a(il,k1p1) + dlp1a(il,k1p2))
            
          aa = 3.0d0 * fac1 / fac2

          bb =  &
               2.0d0 * dpi(il,ij,k1) / (dlp1a(il,k1) + dlp1a(il,k1p1)) -  &
               r23 * aa * (2.0d0 * dlp1a(il,k1) + dlp1a(il,k1p1))
            
          al(il,k1) = qq1a(il,k1) -  &
               dlp1a(il,k1) *  &
               (r13 * aa * dlp1a(il,k1) +  &
               0.5d0 * bb)
            
          al(il,k1p1) = dlp1a(il,k1) * (aa * dlp1a(il,k1) + bb) +  &
               al(il,k1)

!         ---------------------
!         Check if change sign.
!         ---------------------

          if ((qq1a(il,k1) * al(il,k1)) <= 0.0d0) then
             
             al (il,k1) = 0.0d0
             dca(il,k1) = 0.0d0
               
          else
               
             dca(il,k1) = qq1a(il,k1) - al(il,k1)
               
          end if

       end do

!      -------
!      Bottom.
!      -------

       do il = i1, i2

!         ---------------------------------------------------
!         2-cell PPM with zero gradient right at the surface.
!         ---------------------------------------------------

          fac1 = dpi(il,ij,k2m1) * (dlp1a(il,k2) * dlp1a(il,k2)) /  &
                 ((dlp1a(il,k2) + dlp1a(il,k2m1)) *  &
                 (2.0d0 * dlp1a(il,k2) + dlp1a(il,k2m1)))
            
          ar(il,k2) = qq1a(il,k2) + fac1
          al(il,k2) = qq1a(il,k2) - (fac1 + fac1)
            
          if ((qq1a(il,k2) * ar(il,k2)) <= 0.0d0) then
             ar(il,k2) = 0.0d0
          end if
            
          dca(il,k2) = ar(il,k2) - qq1a(il,k2)
            
       end do


!     ----------------------------------------
!     4th order interpolation in the interior.
!     ----------------------------------------

       do ik = k1p2, k2m1
          do il = i1, i2
               
             c1 = (dpi(il,ij,ik-1) * dlp1a(il,ik-1)) /  &
                  (dlp1a(il,ik-1) + dlp1a(il,ik))
               
             c2 = 2.0d0 /  &
                  (dlp1a(il,ik-2) + dlp1a(il,ik-1) +  &
                  dlp1a(il,ik)   + dlp1a(il,ik+1))
               
             a1 = (dlp1a(il,ik-2) + dlp1a(il,ik-1)) /  &
                  (2.0d0 * dlp1a(il,ik-1) + dlp1a(il,ik))
               
             a2 = (dlp1a(il,ik) + dlp1a(il,ik+1)) /  &
                  (2.0d0 * dlp1a(il,ik) + dlp1a(il,ik-1))
               
             al(il,ik) =  &
                  qq1a(il,ik-1) + c1 +  &
                  c2 *  &
                  (dlp1a(il,ik) * (c1 * (a1 - a2) + a2 * dca(il,ik-1)) -  &
                  dlp1a(il,ik-1) * a1 * dca(il,ik))
               
          end do
       end do

       !-----------------------------------------------------------------
       ! Prior to 12/5/08:
       !  Replace these with explicit loops below to facilitate
       ! OpenMP parallelization.  Preserve original code here.
       ! (bmy, 12/5/08)
       !
       !do ik = k1, k2m1
       !   ar(:,ik) = al(:,ik+1)
       !end do         
       !-----------------------------------------------------------------

       do ik = k1, k2m1
       do il = i1, i2
          ar(il,ik) = al(il,ik+1)
       end do
       end do

!      ---------------------------------------
!      Top & Bottom 2 layers always monotonic.
!      ---------------------------------------

       lenx = i2 - i1 + 1
         
       do ik = k1, k1p1
            
          do il = i1, i2
             
             a6(il,ik) =  &
                  3.0d0 * (qq1a(il,ik) + qq1a(il,ik) -  &
                  (al(il,ik)  + ar(il,ik)))
          end do
            
!         ===========
          call Lmtppm  &
               (lenx, 0, a6(i1,ik), al(i1,ik), ar(i1,ik),  &
               dca(i1,ik), qq1a(i1,ik))
!         ===========

       end do


       do ik = k2m1, k2
            
          do il = i1, i2
             
             a6(il,ik) =  &
                  3.0d0 * (qq1a(il,ik) + qq1a(il,ik) -  &
                  (al(il,ik)  + ar(il,ik)))
          end do
            
!         ===========
          call Lmtppm  &
               (lenx, 0, a6(i1,ik), al(i1,ik), ar(i1,ik),  &
                dca(i1,ik), qq1a(i1,ik))
!         ===========

       end do


!      ---------------------------
!      Interior depending on klmt.
!      ---------------------------

!      ==============
       if (klmt == 4) then
!       ==============

!         -------------------------------
!         KORD=7, Huynh's 2nd constraint.
!         -------------------------------

          !-----------------------------------------------------------------
          ! Prior to 12/5/08:
          ! Replace these with explicit loops below to facilitate
          ! OpenMP parallelization.  Preserve original code here.
          ! (bmy, 12/5/08)
          !
          !do ik = k1p1, k2m1
          !   dca(:,ik) = dpi(:,ij,ik) - dpi(:,ij,ik-1)
          !end do
          !-----------------------------------------------------------------

          do ik = k1p1, k2m1
          do il = i1,   i2
             dca(il,ik) = dpi(il,ij,ik) - dpi(il,ij,ik-1)
          end do
          end do


          do ik = k1p2, k2m2
             do il = i1, i2

!             ------------
!             Right edges.
!             ------------

                qmp   = qq1a(il,ik) + (2.0d0 * dpi(il,ij,ik-1))
                lac   = qq1a(il,ik) +  &
                     (1.5d0 * dca(il,ik-1)) + (0.5d0 * dpi(il,ij,ik-1))
                qmin  = Min (qq1a(il,ik), qmp, lac)
                qmax  = Max (qq1a(il,ik), qmp, lac)
                
                ar(il,ik) = Min (Max (ar(il,ik), qmin), qmax)

!             -----------
!             Left edges.
!             -----------
                
                qmp   = qq1a(il,ik) - (2.0d0 * dpi(il,ij,ik))
                lac   = qq1a(il,ik) +  &
                        (1.5d0 * dca(il,ik+1)) - (0.5d0 * dpi(il,ij,ik))
                qmin  = Min (qq1a(il,ik), qmp, lac)
                qmax  = Max (qq1a(il,ik), qmp, lac)
                
                al(il,ik) = Min (Max (al(il,ik), qmin), qmax)

!             -------------
!             Recompute a6.
!             -------------

                a6(il,ik) =  &
                     3.0d0 * (qq1a(il,ik) + qq1a(il,ik) -  &
                     (ar(il,ik)  + al(il,ik)))
             end do
          end do


!      ===================
       else if (klmt <= 2) then
!      ===================

          lenx = 0
            
          do ik = k1p2, k2m2
             do il = i1, i2
                  
                lenx = lenx + 1
                
                al1  (lenx) = al  (il,ik)
                ar1  (lenx) = ar  (il,ik)
                dca1 (lenx) = dca (il,ik)
                qq1a1(lenx) = qq1a(il,ik)
                  
                a61  (lenx) = 3.0d0 * (qq1a1(lenx) + qq1a1(lenx) -  &
                             (al1(lenx)  + ar1(lenx)))
             end do
          end do

!         ===========
          call Lmtppm  &
               (lenx, klmt, a61, al1, ar1, dca1, qq1a1)
!         ===========

          lenx = 0
            
          do ik = k1p2, k2m2
             do il = i1, i2
                  
                lenx = lenx + 1
                
                a6  (il,ik) = a61  (lenx)
                al  (il,ik) = al1  (lenx)
                ar  (il,ik) = ar1  (lenx)
                dca (il,ik) = dca1 (lenx)
                qq1a(il,ik) = qq1a1(lenx)
                
             end do
          end do


       end if

       
       do ik = k1, k2m1
          do il = i1, i2
               
             if (wza(il,ik) > 0.0d0) then
                
                cm = wza(il,ik) / dlp1a(il,ik)
                
                dca(il,ik+1) =  &
                     ar(il,ik) +  &
                     0.5d0 * cm *  &
                     (al(il,ik) - ar(il,ik) +  &
                     a6(il,ik) * (1.0d0 - r23 * cm))

             else

                cp = wza(il,ik) / dlp1a(il,ik+1)

                dca(il,ik+1) =  &
                     al(il,ik+1) +  &
                     0.5d0 * cp *  &
                     (al(il,ik+1) - ar(il,ik+1) -  &
                     a6(il,ik+1) * (1.0d0 + r23 * cp))
                  
             end if
               
          end do
       end do
         
         
       !-----------------------------------------------------------------
       ! Prior to 12/5/08:
       ! Replace these with explicit loops below to facilitate
       ! OpenMP parallelization.  Preserve original code here.
       ! (bmy, 12/5/08)
       !
       !do ik = k1, k2m1
       !   dca(:,ik+1) = wza(:,ik) * dca(:,ik+1)
       !   !.sds.. save vertical flux for species as diagnostic
       !   fz(i1:i2,ij,ik+1) = dca(:,ik+1)
       !end do
       !
       !dq1(i1:i2,ij,k1) = dq1(i1:i2,ij,k1) - dca(:,k1p1)
       !dq1(i1:i2,ij,k2) = dq1(i1:i2,ij,k2) + dca(:,k2)
       !
       !do ik = k1p1, k2m1
       !     
       !   dq1(i1:i2,ij,ik) =  &
       !        dq1(i1:i2,ij,ik) + dca(:,ik) - dca(:,ik+1)
       !     
       !end do
       !-----------------------------------------------------------------

       do ik = k1, k2m1
       do il = i1, i2
          dca(il,ik+1) = wza(il,ik) * dca(il,ik+1)
          !.sds.. save vertical flux for species as diagnostic
          fz(il,ij,ik+1) = dca(il,ik+1)
       end do
       end do

       do il = i1, i2 
          dq1(il,ij,k1) = dq1(il,ij,k1) - dca(il,k1p1)
          dq1(il,ij,k2) = dq1(il,ij,k2) + dca(il,k2)
       enddo

       do ik = k1p1, k2m1
       do il = i1,   i2

          dq1(il,ij,ik) =  &
               dq1(il,ij,ik) + dca(il,ik) - dca(il,ik+1)
            
       end do
       end do

!   =============
    end do ijloop
    !$OMP END PARALLEL DO
!   =============


  END SUBROUTINE Fzppm

!EOC 
!------------------------------------------------------------------------------
!          Harvard University Atmospheric Chemistry Modeling Group            !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Average_Press_Poles
!
! !DESCRIPTION: Subroutine Average\_Press\_Poles averages pressure at the 
!  Poles when the Polar cap is enlarged.  It makes the last two latitudes 
!  equal.
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE Average_Press_Poles( area_1D, press, I1,  I2,   JU1,  &
                                  J2,      ILO,   IHI, JULO, JHI )
!
! !INPUT PARAMETERS: 
!
    ! Local min & max longitude (I), latitude (J)    
    INTEGER, INTENT(IN)   :: I1,     I2
    INTEGER, INTENT(IN)   :: JU1,    J2

    ! Local min & max longitude (I) and latitude (J) indices
    INTEGER, INTENT(IN)   :: ILO,    IHI
    INTEGER, INTENT(IN)   :: JULO,   JHI

    ! Surface area of grid box
    REAL*8,  INTENT(IN)   :: AREA_1D(JU1:J2)
!
! !INPUT/OUTPUT PARAMETERS: 
!
    ! Surface pressure [hPa]
    REAL*8, INTENT(INOUT) :: press(ILO:IHI, JULO:JHI)
!
! !AUTHOR:
!   Philip Cameron-Smith and John Tannahill, GMI project @ LLNL (2003) 
!   Implemented into GEOS-Chem by Claire Carouge (ccarouge@seas.harvard.edu)
!
! !REMARKS:
!  Subroutine from pjc_pfix.  Call this one once everything is working fine.
! 
! !REVISION HISTORY: 
!   05 Dec 2008 - C. Carouge  - Replaced TPCORE routines by S-J Lin and Kevin
!                               Yeh with the TPCORE routines from GMI model.
!                               This eliminates the polar overshoot in the
!                               stratosphere.
!   05 Dec 2008 - R. Yantosca - Updated documentation and added ProTeX headers.
!                               Declare all REAL variables as REAL*8.  Also 
!                               make sure all numerical constants are declared
!                               with the "D" double-precision exponent.
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
! 
    ! Scalars
    INTEGER :: I, J
    REAL*8  :: meanp
    REAL*8  :: REL_AREA(JU1:J2)
    REAL*8  :: SUM_AREA

    !----------------
    !Begin execution.
    !----------------
  
!    ! Calculate rel_area for each lat. (ccc, 11/20/08)
!    DO J = JU1, J2
!       REL_AREA(J) = AREA_1D(J) / ( SUM( AREA_1D ) * I2 )
!    ENDDO
!
!    ! Average pressure in the S. Polar zone
!    meanp = 0.0d0
!    DO I = i1,i2
!       meanp = meanp + SUM( rel_area(JU1:JU1+1)  * press(I,JU1:JU1+1) )
!    ENDDO
!    meanp = meanp / ( SUM( rel_area(JU1:JU1+1) ) * I2 )
!      
!    press(i1:i2,ju1:ju1+1) = meanp
!
!    ! Average pressures in the N. Polar zone
!    meanp = 0.0d0
!    DO I = i1,i2
!       meanp = meanp + SUM( rel_area(J2-1:J2)  * press(I,J2-1:J2) )
!    ENDDO
!    meanp = meanp / ( SUM( rel_area(J2-1:J2) ) * I2 )
!
!    press(I1:I2,J2-1:J2) = meanp


    ! Compute the sum of surface area 
    SUM_AREA = SUM( AREA_1D ) * DBLE( I2 )

    ! Calculate rel_area for each lat. (ccc, 11/20/08)
    DO J = JU1, J2
       REL_AREA(J) = AREA_1D(J) / SUM_AREA
    ENDDO

    !--------------
    ! South Pole 
    !--------------

    ! Surface area of the S. Polar cap
    SUM_AREA = SUM( rel_area( JU1:JU1+1 ) ) * DBLE( I2 )

    ! Zero
    meanp = 0.d0

    ! Sum pressure * surface area over the S. Polar cap
    DO J = JU1, JU1+1
    DO I = I1,  I2
       meanp = meanp + ( rel_area(J)  * press(I,J) )
    ENDDO
    ENDDO

    ! Normalize pressure in all grid boxes w/in the S. Polar cap
    press( :, JU1:JU1+1 ) = meanp / SUM_AREA

    !--------------
    ! North Pole 
    !--------------

    ! Surface area of the N. Polar cap
    SUM_AREA = SUM( rel_area( J2-1:J2 ) ) * DBLE( I2 )

    ! Zero
    meanp = 0.d0

    ! Sum pressure * surface area over the N. Polar cap
    DO J = J2-1, J2
    DO I = I1,   I2
       meanp = meanp + ( rel_area(J)  * press(I,J) )
    ENDDO
    ENDDO

    ! Normalize pressure in all grid boxes w/in the N. Polar cap
    press( :, J2-1:J2 ) = meanp / SUM_AREA

  END SUBROUTINE Average_Press_Poles

END MODULE Tpcore_FvDas_mod
!EOC
